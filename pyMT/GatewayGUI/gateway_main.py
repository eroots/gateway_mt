"""Project-manager / launcher window for the pyMT GUI suite.

This module provides :class:`NewProject`, a small PyQt5 window that sits
"above" the other GUIs.  It shows a table of named *datasets* -- each a
bundle of list / data / response / model / resolution file paths -- lets
that table be saved to and loaded from a pyMT project file (``.pymt``,
the same pystart format read by
:meth:`pyMT.GUI_common.classes.FileInputParser.read_pystart`), and offers
three buttons that launch the individual GUIs preloaded with the selected
dataset(s):

* **Data Plot** -> :class:`pyMT.DataGUI.data_plot.DataMain`
* **Model Viewer** -> :class:`pyMT.ModelGUI.model_viewer.ModelWindow`
* **Mesh Designer** -> :class:`pyMT.ModelGUI.mesh_designer.model_viewer_2d`

The window itself is built from ``new_project.ui`` via
:func:`PyQt5.uic.loadUiType`; widgets are wired by ``objectName``.  The
console entry point ``gateway_mt`` (see ``setup.py``) calls :func:`main`.

Notes
-----
``gateway_main.ui`` and ``new_project_treeFormat.ui`` in this directory
are unused alternative layouts kept for reference; only
``new_project.ui`` is loaded.
"""

import os
os.environ['QT_SCREEN_SCALE_FACTOR'] = '1'
os.environ['QT_ENABLE_HIGHDPI_SCALING'] = '1'
from PyQt5.uic import loadUiType
from PyQt5 import QtWidgets, QtCore, QtGui
from pyMT.GUI_common.classes import FileInputParser
from pyMT.DataGUI.data_plot import DataMain
from pyMT.ModelGUI.model_viewer import ModelWindow
from pyMT.ModelGUI.mesh_designer import model_viewer_2d as MeshDesigner
import sys

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources
from pyMT import resources
from pyMT.WSExceptions import WSFileError


path = os.path.dirname(os.path.realpath(__file__))
with pkg_resources.path(resources, 'model_viewer.jpg') as p:
    model_viewer_jpg = str(p)
with pkg_resources.path(resources, 'data_plot.jpg') as p:
    data_plot_jpg = str(p)
with pkg_resources.path(resources, 'mesh_designer.jpg') as p:
    mesh_designer_jpg = str(p)
Ui_NewProject, QNewProject = loadUiType(os.path.join(path, 'new_project.ui'))


class NewProject(QNewProject, Ui_NewProject):
    """Project-manager window: edit, save and load a table of datasets.

    Each row of ``datasetTable`` is one dataset -- a name plus up to five
    file paths (list, data, response, model, resolution).  The table can
    be written to / read from a ``.pymt`` project file, and the currently
    selected row(s) can be handed to the Data Plot, Model Viewer or Mesh
    Designer GUIs via the three launch buttons.

    Built from ``new_project.ui``; multiple-row selection is enabled so
    Data Plot can be launched with several datasets at once.

    Attributes
    ----------
    datasetTable : PyQt5.QtWidgets.QTableWidget
        The dataset table; column 0 is the name, columns 1-5 the files.
    file_to_index, index_to_file : dict
        Two-way map between the file-role name (``'list'``, ``'data'``,
        ``'response'``, ``'model'``, ``'resolution'``) and its column.
    project_path, project_name : str
        Directory and full path of the currently loaded project file;
        used to resolve relative paths before launching a GUI.
    dp_windows, model_windows, md_windows : list
        Live references to spawned Data Plot / Model Viewer / Mesh
        Designer windows, kept so they are not garbage-collected.
    table_modified : bool
        ``True`` once the user has edited a cell; gates the
        "save first?" prompt in :meth:`open_project`.
    """

    def __init__(self):
        super(NewProject, self).__init__()
        self.setupUi(self)
        self.setWindowFlag(QtCore.Qt.WindowMinimizeButtonHint, True)
        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.any_datasets = False
        self.datasets = {}
        self.current_dir = ''
        self.project_path = ''
        self.project_name = ''
        self.dp_windows = []
        self.model_windows = []
        self.md_windows = []
        self.table_modified = False
        self.file_to_index = {'dataset': 0, 'list': 1, 'data': 2, 'response': 3, 'model': 4, 'resolution': 5}
        self.index_to_file = {0: 'dataset', 1: 'list', 2: 'data', 3: 'response', 4: 'model', 5: 'resolution'}
        self.setup_widgets()
        ret = self.open_project()
        if not ret:
            self.init_dataset_table()
        self.datasetTable.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

    def dataset_index(self):
        """Row indices of the currently selected datasets.

        Returns
        -------
        list of int
            Sorted, de-duplicated row numbers of every selected cell
            (empty if nothing is selected).
        """
        idx = (index.row() for index in self.datasetTable.selectedIndexes())
        idx = set(idx)
        idx = sorted(idx)
        return idx
        # return sorted(set(index.row() for index in
                      # self.datasetTable.selectedIndexes()))

    @property
    def empty_table(self):
        """bool: ``True`` when every cell of ``datasetTable`` is blank."""
        for row in range(self.datasetTable.rowCount()):
            for col in range(self.datasetTable.columnCount()):
                item = self.datasetTable.item(row, col)
                if item:
                    if item.text():
                        return False
        return True

    def dataset_dict(self):
        """Every dataset in the table as a nested dict.

        Returns
        -------
        dict
            ``{dataset_name: {file_role: path, ...}, ...}`` over all
            rows, where ``file_role`` is one of ``'list'``, ``'data'``,
            ``'response'``, ``'model'``, ``'resolution'``.  Only
            populated cells contribute a key.

        See Also
        --------
        active_datasets : the same shape, restricted to selected rows.
        """
        ds_dict = {}
        for row in range(self.datasetTable.rowCount()):
            ds_name = self.datasetTable.item(row, 0).text()
            ds_dict.update({ds_name: {}})
            for col in range(1, self.datasetTable.columnCount()):
                item = self.datasetTable.item(row, col)
                if item:
                    ds_dict[ds_name].update({self.index_to_file[col]: item.text()})
        return ds_dict

    def active_datasets(self):
        """Selected datasets as a nested dict.

        Returns
        -------
        dict
            Same shape as :meth:`dataset_dict` but only for the rows in
            :meth:`dataset_index`.  This is what the launch handlers pass
            to the individual GUIs.
        """
        ds_dict = {}
        active_rows = self.dataset_index()
        for row in active_rows:
            ds_name = self.datasetTable.item(row, 0).text()
            ds_dict.update({ds_name: {}})
            for col in range(1, self.datasetTable.columnCount()):
                item = self.datasetTable.item(row, col)
                if item:
                    ds_dict[ds_name].update({self.index_to_file[col]: item.text()})
        return ds_dict

    def setup_widgets(self):
        """Connect the buttons / table signals and set the launch icons.

        Called once from :meth:`__init__`.
        """
        self.addDataset.clicked.connect(self.add_dataset)
        self.deleteDataset.clicked.connect(self.delete_dataset)
        self.datasetTable.cellChanged.connect(self.modify_table_cell)
        self.saveProject.clicked.connect(self.save_project)
        self.openProject.clicked.connect(self.open_project)
        self.browseButton.clicked.connect(self.browse_files)
        self.projectName.setText('')
        self.launchDataPlot.clicked.connect(self.launch_data_plot)
        self.launchDataPlot.setIcon(QtGui.QIcon(data_plot_jpg))
        self.launchDataPlot.setIconSize(QtCore.QSize(128,128))
        self.launchModelViewer.clicked.connect(self.launch_model_viewer)
        self.launchModelViewer.setIcon(QtGui.QIcon(model_viewer_jpg))
        self.launchModelViewer.setIconSize(QtCore.QSize(128,128))
        self.launchMeshDesigner.clicked.connect(self.launch_mesh_designer)
        self.launchMeshDesigner.setIcon(QtGui.QIcon(mesh_designer_jpg))
        self.launchMeshDesigner.setIconSize(QtCore.QSize(128,128))

    def sort_files(self, files):
        """Classify a list of paths into dataset file roles by extension.

        Parameters
        ----------
        files : list of str
            Arbitrary file paths, e.g. the result of a file-browser
            dialog.

        Returns
        -------
        dict
            ``{role: path}`` where ``role`` is one of ``'model'``,
            ``'data'``, ``'response'``, ``'list'``, ``'resolution'``.
            If several paths map to the same role the last one wins.

        Notes
        -----
        ``.rho`` / ``.model`` / ``.mod`` / ``.zani`` / ``.ani`` -> model,
        ``.lst`` / ``.list`` -> list, ``.reso`` / ``.resolution`` ->
        resolution, ``.resp`` -> response.  For ``.dat`` / ``.data`` /
        ``.adat`` / ``.gdat`` the first line decides data vs response: a
        ``.gdat`` starting with ``#`` (or any file whose first line does
        not contain ``"response"``) is treated as data.  pyMT-written
        ModEM response files do not contain the word ``"response"``, so
        they are classified as data unless given a ``.resp`` extension.
        """
        ret_dict = {}

        for file in files:
            name, ext = os.path.splitext(file)
            if ext in ('.rho', '.model', '.mod', '.zani', '.ani'):
                ret_dict.update({'model': file})
            elif ext in ('.dat', '.data', '.adat', '.gdat'):
                with open(file, 'r') as f:
                    line = f.readline()
                if ext == '.gdat':
                    if line.startswith('#'):
                        ret_dict.update({'data': file})
                    else:
                        ret_dict.update({'response': file})
                elif 'response' in line:
                    ret_dict.update({'response': file})
                else:
                    ret_dict.update({'data': file})
            elif ext in ('.lst', '.list'):
                ret_dict.update({'list': file})
            elif ext in ('.reso', '.resolution'):
                ret_dict.update({'resolution': file})
            elif ext in ('.resp',):
                ret_dict.update({'response': file})
        return ret_dict

    def init_dataset_table(self, files_dict=None):
        """(Re)build ``datasetTable`` from a project dict, or with one blank row.

        Parameters
        ----------
        files_dict : dict, optional
            ``{dataset_name: {role: path, ...}, ...}`` as returned by
            :meth:`pyMT.GUI_common.classes.FileInputParser.read_pystart`.
            If ``None`` (or empty) the table is reset to a single row
            named ``'ds1'`` with blank file cells.

        Notes
        -----
        The caller is responsible for disconnecting
        ``datasetTable.cellChanged`` around this call if the
        :meth:`modify_table_cell` slot should not fire for the
        programmatic edits (see :meth:`open_project`).
        """
        header = ('Dataset', 'List', 'Data', 'Response', 'Model', 'Resolution') #, 'Path')
        self.datasetTable.setColumnCount(len(header))
        for ii, label in enumerate(header):
            self.datasetTable.setHorizontalHeaderItem(ii, QtWidgets.QTableWidgetItem(label))
        if files_dict:
            row = 0
            self.datasetTable.setRowCount(len(files_dict))
            for dataset_name, startup in files_dict.items():
                try:
                    self.datasetTable.setItem(row, 0, QtWidgets.QTableWidgetItem(dataset_name))
                    self.datasetTable.setItem(row, 1, QtWidgets.QTableWidgetItem(startup.get('list', '')))
                    self.datasetTable.setItem(row, 2, QtWidgets.QTableWidgetItem(startup.get('data', '')))
                    self.datasetTable.setItem(row, 3, QtWidgets.QTableWidgetItem(startup.get('response', '')))
                    self.datasetTable.setItem(row, 4, QtWidgets.QTableWidgetItem(startup.get('model', '')))
                    self.datasetTable.setItem(row, 5, QtWidgets.QTableWidgetItem(startup.get('resolution', '')))
                except KeyError:
                    pass
                row += 1
        else:
            self.datasetTable.setRowCount(1)
            self.datasetTable.setItem(0, 0, QtWidgets.QTableWidgetItem('ds1'))
            for col in range(1, 6):
                self.datasetTable.setItem(0, col, QtWidgets.QTableWidgetItem(''))

    def update_dataset_table(self, files):
        """Drop browsed files into the selected row, one per role column.

        Parameters
        ----------
        files : list of str
            Paths to classify with :meth:`sort_files` and write into the
            first selected row.  Only reached from :meth:`browse_files`,
            which guarantees exactly one selected row.
        """
        sorted_files = self.sort_files(files=files)
        for key, value in sorted_files.items():
            self.datasetTable.setItem(self.dataset_index()[0], self.file_to_index[key], QtWidgets.QTableWidgetItem(value))

    def add_dataset(self):
        """Append a new empty row named ``ds<n>`` (``addDataset`` button)."""
        self.datasetTable.setRowCount(self.datasetTable.rowCount() + 1)
        self.datasetTable.setItem(self.datasetTable.rowCount() - 1, 0,
                                  QtWidgets.QTableWidgetItem('ds{}'.format(self.datasetTable.rowCount())))

    def delete_dataset(self):
        """Remove every selected row (``deleteDataset`` button).

        Rows are removed high index first so earlier removals do not
        shift the indices still to be deleted.
        """
        for idx in reversed(sorted(self.dataset_index())):
            self.datasetTable.removeRow(idx)

    def modify_table_cell(self, row, col):
        """``datasetTable.cellChanged`` slot: mark dirty, enforce unique names.

        Parameters
        ----------
        row, col : int
            The cell that changed.

        Notes
        -----
        When a name cell (column 0) is edited to duplicate another row's
        name, the user is prompted (via ``QInputDialog.getText``) for a
        replacement.  Known issue: the prompt reopens until a non-empty
        name is entered, so it cannot be cancelled.
        """
        self.table_modified = True
        item = self.datasetTable.item(row, col)
        if item:
            if item.column() == 0:
                for row in range(self.datasetTable.rowCount()):
                    if row != item.row():
                        if item.text() == self.datasetTable.item(row, 0).text():
                            text = ''
                            while True:
                                if text:
                                    item.setText(text)
                                    break
                                else:
                                    text, ok = QtWidgets.QInputDialog.getText(self, '', 'Data sets must have unique names. Enter a new name: ')

    def save_project(self):
        """Collect the whole table and hand it to :meth:`write_project`.

        Every non-empty row must have a name; a nameless row raises a
        warning dialog and aborts the save.  (``saveProject`` button.)
        """
        all_datasets = {}
        for row in range(self.datasetTable.rowCount()):
            dataset = []
            for col in range(self.datasetTable.columnCount()):
                item = self.datasetTable.item(row, col)
                if item:
                    dataset.append(self.datasetTable.item(row, col).text())
                else:
                    dataset.append('')
            if not dataset[0]:
                QtWidgets.QMessageBox.warning(self, '...', 'Dataset in row {} needs a name!'.format(row + 1))
                return
            if any(dataset):
                all_datasets.update({dataset[0]: {'list': dataset[1],
                                                  'data': dataset[2],
                                                  'resp': dataset[3],
                                                  'model': dataset[4],
                                                  'resolution': dataset[5]}})
        self.write_project(all_datasets)

    def open_project(self):
        """Prompt for a project file and load it into the table.

        Called from :meth:`__init__` (so the window opens straight onto a
        chosen project) and from the ``openProject`` button.  If the
        current table has unsaved edits the user is first asked whether
        to save; choosing *Cancel* aborts and returns ``False``.

        Returns
        -------
        bool or None
            ``True`` if a project file was chosen and loaded, ``False``
            if the user cancelled, ``None`` if no file was picked.

        Notes
        -----
        A *save* dialog (with overwrite confirmation disabled) is used to
        pick the file to open, so that typing a new name here starts a
        fresh project rather than erroring -- see :meth:`write_project`.
        """
        files_dict = None

        if self.empty_table or not self.table_modified:
            pass
        else:
            retval = QtWidgets.QMessageBox.question(self, '', 'Do you want to save the current project first?',
                                                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel)
            if retval == QtWidgets.QMessageBox.Yes:
                self.save_project()
            elif retval == QtWidgets.QMessageBox.No:
                pass
            elif retval == QtWidgets.QMessageBox.Cancel:
                return False
        project_file = QtWidgets.QFileDialog.getSaveFileName(self, 'Open Project', self.project_path,
                                                             'pyMT Project Files (*.pymt);; All Files (*)',
                                                             options=QtWidgets.QFileDialog.DontConfirmOverwrite)[0]
        if project_file:
            if os.path.exists(project_file):
                files_dict = FileInputParser.read_pystart(project_file)
            self.current_dir = os.path.abspath(project_file)
            self.projectName.setText(os.path.basename(project_file))
            self.datasetTable.cellChanged.disconnect()
            self.init_dataset_table(files_dict)
            self.datasetTable.cellChanged.connect(self.modify_table_cell)
            self.project_path = os.path.dirname(project_file)
            self.project_name = project_file
            return True

    def write_project(self, all_datasets):
        """Write ``all_datasets`` to a user-chosen ``.pymt`` file.

        Parameters
        ----------
        all_datasets : dict
            ``{name: {'list': ..., 'data': ..., 'resp': ..., 'model':
            ..., 'resolution': ...}}`` as built by :meth:`save_project`.

        Notes
        -----
        Output is pystart format -- ``% <name>`` then one ``<role>
        <path>`` line per non-empty file.  Blank paths and missing roles
        are skipped.  On success ``table_modified`` is cleared.
        """
        project_file = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Project', self.project_path + self.project_name, 'pyMT Project Files (*.pymt);; All Files (*)')[0]
        if project_file:
            self.current_dir = os.path.abspath(project_file)
            with open(project_file, 'w') as f:
                for name, dataset in all_datasets.items():
                    f.write('% {}\n'.format(name))
                    for file_type in ('list', 'data', 'resp', 'model', 'resolution'):
                        try:
                            file = dataset[file_type]
                            if file:
                                f.write('{} {}\n'.format(file_type, file))
                        except KeyError:
                            pass
                self.projectName.setText(os.path.basename(project_file))
            self.table_modified = False

    def browse_files(self):
        """Open a file picker and load the chosen files into the selected row.

        Requires exactly one selected row (otherwise a warning dialog is
        shown).  Files are routed to columns by :meth:`sort_files` via
        :meth:`update_dataset_table`.  (``browseButton``.)
        """
        if len(self.dataset_index()) == 1:
            fnames = QtWidgets.QFileDialog.getOpenFileNames(self, 'Browse Files',
                                                           self.current_dir,
                                                           'All Files (*.rho *.model *.zani *.dat *.gdat *.data *.resp *.lst *.reso);;' +\
                                                           'ModEM Files (*.rho *.model *.dat *.lst *.reso *.data *.resp);;' +\
                                                           'EM3DANI Files (*.mod *.dat *.adat *.resp *.lst);;' +\
                                                           'MT3DANI Files (*.ani *.zani *.dat);; All Files (*)')[0]
            if fnames:
                self.current_dir = os.path.abspath(fnames[0])
                self.update_dataset_table(fnames)
        else:
            QtWidgets.QMessageBox.warning(self, '', 'You must select (highlight) a single dataset to add to!')

    def launch_data_plot(self):
        """Open a Data Plot window on the selected dataset(s).

        Relative paths are resolved against ``project_path`` first.  Any
        number of datasets may be selected.  The new
        :class:`~pyMT.DataGUI.data_plot.DataMain` is kept in
        ``dp_windows``; a :class:`~pyMT.WSExceptions.WSFileError` during
        load is shown as a warning dialog.  (``launchDataPlot``.)
        """
        if len(self.dataset_index()) == 0:
            QtWidgets.QMessageBox.warning(self, '', 'Select (highlight) the desired dataset(s) first.')
            return
        active_ds = self.active_datasets()
        try:
            for ds in active_ds.keys():
                for key in active_ds[ds].keys():
                    if active_ds[ds][key] and not os.path.isabs(active_ds[ds][key]):
                        active_ds[ds].update({key: self.project_path + '/' + active_ds[ds][key]})
            if (not any(list(active_ds.values()))) or (not active_ds):
                QtWidgets.QMessageBox.warning(self, '', 'Selected data set is empty.')
                return
            dp_main = DataMain(dataset_dict=active_ds, edi_locs_from=self.get_locs_from())
            dp_main.setWindowIcon(QtGui.QIcon(data_plot_jpg))
            self.dp_windows.append(dp_main)
            self.dp_windows[-1].show()
        except WSFileError as e:
            QtWidgets.QMessageBox.warning(self, 'File Not Found', e.message)

    def launch_model_viewer(self):
        """Open a Model Viewer window on the single selected dataset.

        Requires exactly one selected row and a ``model`` file in it.
        Relative paths are resolved against ``project_path``.  The new
        :class:`~pyMT.ModelGUI.model_viewer.ModelWindow` is kept in
        ``model_windows``.  (``launchModelViewer``.)
        """
        if len(self.dataset_index()) == 0:
            QtWidgets.QMessageBox.warning(self, '', 'Select (highlight) the desired dataset(s) first.')
            return
        elif len(self.dataset_index()) > 1:
            QtWidgets.QMessageBox.warning(self, '', 'Model Viewer can only load one data set at a time (for now)...')
            return
        active_ds = self.active_datasets()
        ds = list(active_ds.values())[0]
        model = ds.get('model', None)
        if model:
            for key in ds.keys():
                if ds[key] and not os.path.isabs(ds[key]):
                    ds[key] = self.project_path + '/' + ds[key]
            try:
                mv_main = ModelWindow(files=ds)
                mv_main.setWindowIcon(QtGui.QIcon(model_viewer_jpg))
                mv_main.show()
                mv_main.closed.connect(self.close_mv)
                self.model_windows.append(mv_main)
                self.model_windows[-1].show()
            except WSFileError as e:
                QtWidgets.QMessageBox.warning(self, 'File Not Found', e.message)
        else:
            QtWidgets.QMessageBox.warning(self, '', 'A model must be available in the selected data set to use Model Viewer!')

    def close_mv(self, mv):
        """Slot for ``ModelWindow.closed``.

        Currently a no-op, so closed Model Viewer windows stay in
        ``model_windows`` (known issue).
        """
        pass

    def launch_mesh_designer(self):
        """Open a Mesh Designer window on the single selected dataset.

        Requires exactly one selected row with a ``model`` and/or
        ``data`` file.  The new
        :class:`~pyMT.ModelGUI.mesh_designer.model_viewer_2d` is kept in
        ``md_windows``.  (``launchMeshDesigner``.)
        """
        if len(self.dataset_index()) == 0:
            QtWidgets.QMessageBox.warning(self, '', 'Select (highlight) the desired dataset(s) first.')
            return
        elif len(self.dataset_index()) > 1:
            QtWidgets.QMessageBox.warning(self, '', 'Mesh Designer can only load one data set at a time (for now)...')
            return
        active_ds = self.active_datasets()
        ds = list(active_ds.values())[0]
        model = ds.get('model', None)
        data = ds.get('data', None)
        if model or data:
            for key in ds.keys():
                if ds[key] and not os.path.isabs(ds[key]):
                    ds[key] = self.project_path + '/' + ds[key]
            try:
                md_main = MeshDesigner(model=model, dataset=data, path=self.project_path)
                md_main.setWindowIcon(QtGui.QIcon(model_viewer_jpg))
                self.md_windows.append(md_main)
                self.md_windows[-1].show()
            except WSFileError as e:
                QtWidgets.QMessageBox.warning(self, 'File Not Found', e.message)
        else:
            QtWidgets.QMessageBox.warning(self, '', 'Either a model or data file (or both) must be available to use Mesh Designer!')

    def get_locs_from(self):
        """str: where Data Plot should read EDI site locations from.

        The current text of the ``ediLocs`` combo box, lower-cased
        (e.g. ``'definemeas'``); passed as ``edi_locs_from`` to
        :class:`~pyMT.DataGUI.data_plot.DataMain`.
        """
        return self.ediLocs.currentText().lower()


def main():
    """Console entry point (``gateway_mt``): run the project manager."""
    app = QtWidgets.QApplication(sys.argv)  # Starts GUI event loop
    mainGUI = NewProject()  # Instantiate a GUI window
    mainGUI.show()  # Show it.
    ret = app.exec_()
    sys.exit(ret)  # Properly close the loop when the window is closed.


if __name__ == '__main__':
    main()
