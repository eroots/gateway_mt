"""3-D / 2-D viewer for ModEM resistivity models.

``model_viewer`` is a PyQt window for inspecting a ModEM model: a 3-D
PyVista (`pyvistaqt.QtInteractor`) render view with orthogonal slice
planes, volume clipping, isosurfaces, an optional mesh overlay, and an
interactively-picked transect, alongside a matplotlib tab that draws the
plan / X / Y slices (through :class:`pyMT.gplot.MapView`) plus a separate
transect figure.  It can overlay station locations (from a data file)
and shade the model by a resolution / sensitivity volume.

ModEM-only, and largely superseded by ParaView for general 3-D
inspection; kept for the tightly-integrated 2-D slice + transect views.

Launching
---------
Installed as the ``model_viewer`` console script (:func:`main`).  Takes
positional file paths in any order; :func:`pyMT.utils.sort_files`
classifies them by extension into ``{'model', 'data'/'response',
'list', 'resolution', ...}``.  Only a model file is required.

Relation to the rest of pyMT
----------------------------
Reads models / data through :mod:`pyMT.data_structures`; reuses
:class:`pyMT.gplot.MapView` for the 2-D matplotlib panels; colour maps
from :mod:`pyMT.e_colours.colourmaps`; menus from
:mod:`pyMT.GUI_common.classes`.

Notes
-----
Importing this module sets global PyVista theme state
(``pv.set_plot_theme('dark')`` and a few ``pv.global_theme`` fields).
"""
from os import path as ospath
import sys
from copy import deepcopy
import numpy as np
import pyvista as pv
pv.set_plot_theme('dark')
pv.global_theme.allow_empty_mesh = True
pv.global_theme.font.size = 10
pv.global_theme.font.family = 'arial'
from pyvistaqt import QtInteractor
import pyMT.data_structures as DS
from pyMT.utils import sort_files, check_file
import pyMT.e_colours.colourmaps as cm
from pyMT.GUI_common import classes
from pyMT import gplot
from pyMT import resources
from PyQt5.uic import loadUiType
from pyMT.IO import debug_print
from scipy.interpolate import RegularGridInterpolator as RGI
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from PyQt5 import QtWidgets, QtGui, QtCore
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources


path = ospath.dirname(ospath.realpath(__file__))

with pkg_resources.path(resources, 'model_viewer.jpg') as p:
    model_viewer_jpg = str(p)
UI_ModelWindow, QModelWindow = loadUiType(ospath.join(path, 'model_viewer.ui'))


def model_to_rectgrid(model, resolution=None, rho_axis='rho_x'):
    """Convert a :class:`pyMT.data_structures.Model` to a PyVista grid.

    Parameters
    ----------
    model : pyMT.data_structures.Model
        The model to convert.  Its cell-edge coordinates ``dx`` / ``dy``
        / ``dz`` become the grid axes (``dz`` is negated so depth points
        down).
    resolution : pyMT.data_structures.Model, optional
        A resolution / sensitivity volume on the same mesh.  When given,
        two extra cell-data arrays are added: ``'RawResolution'``
        (:math:`\\log_{10}` of the raw values) and ``'Resolution'`` (the
        values divided by cell volume, clipped at the median, then
        standardised into roughly ``[0, 1]`` for use as a PyVista opacity
        field).  When omitted, ``'Resolution'`` is all ones.
    rho_axis : str, optional
        Model attribute to take resistivities from -- ``'rho_x'``
        (default), ``'rho_y'``, ``'rho_z'``, or ``'a / b'`` for the ratio
        of two attributes (anisotropic models only).

    Returns
    -------
    pyvista.RectilinearGrid
        Cell data ``'Resistivity'`` holds :math:`\\log_{10}\\rho` with the
        model's ``(nx, ny, nz)`` array transposed / flipped into PyVista's
        ``(y, x, -z)`` Fortran order.

    Raises
    ------
    ValueError
        If ``rho_axis`` names a ratio and one of the attributes is empty
        (e.g. ``'rho_x/rho_y'`` on an isotropic model).
    """
    if '/' in rho_axis:
        rho_axis = rho_axis.split('/')
        v1 = getattr(model, rho_axis[0].strip())
        v2 = getattr(model, rho_axis[1].strip())
        vals = v1 / v2
    else:
        vals = getattr(model, rho_axis, 'vals')
    grid = pv.RectilinearGrid(np.array(model.dy),
                              np.array(model.dx),
                              -np.array(model.dz[::-1]))
    vals = np.log10(np.swapaxes(np.flip(vals, 2), 0, 1)).flatten(order='F')
    grid.cell_data['Resistivity'] = vals
    if resolution:
        X, Y, Z = np.meshgrid(resolution.yCS, resolution.xCS, resolution.zCS)
        volume = np.swapaxes(np.flip(X * Y * Z, 2), 0, 1).flatten(order='F')
        resvals = np.swapaxes(np.flip(resolution.vals, 2), 0, 1).flatten(order='F')
        grid.cell_data['RawResolution'] = np.log10(resvals)
        resvals = resvals / (volume.ravel())# ** (2/3))
        resvals[resvals > np.median(resvals.ravel())] = np.median(resvals.ravel())
        resvals = resvals - np.mean(resvals)
        resvals = resvals / np.std(resvals)
        resvals = 0.5 + resvals * np.sqrt(0.5)
        resvals[resvals > 1] = 1
        grid.cell_data['Resolution'] = resvals
    else:
        grid.cell_data['Resolution'] = np.ones(vals.shape)
        grid.cell_data['Resolution'][0] = .99
    return grid


class ModelWindow(QModelWindow, UI_ModelWindow):
    """Main viewer window (``model_viewer.ui``).

    Owns the model / dataset, a :class:`pyvistaqt.QtInteractor` 3-D view
    (``self.vtk_widget``), and a :class:`pyMT.gplot.MapView` driving the
    2-D matplotlib panels (``self.map``).  The menu actions, sliders and
    text boxes read GUI widgets, update ``self.model`` / ``self.clip_model``
    / ``self.rect_grid`` / colour state, and call one of the redraw
    entry points -- :meth:`render_3D` (3-D), :meth:`update_2D_X` /
    :meth:`update_2D_Y` / :meth:`update_plan_view` (2-D), or
    :meth:`update_all` (both).

    Parameters
    ----------
    files : dict
        ``{'model': ..., 'data'/'response': ..., 'list': ...,
        'resolution': ...}`` from :func:`pyMT.utils.sort_files`.  Only
        ``'model'`` is required.
    parent : QtWidgets.QWidget, optional

    Attributes
    ----------
    dataset : pyMT.data_structures.Dataset
    model, clip_model : pyMT.data_structures.Model
        The full model and the clipped copy actually rendered.
    rect_grid : pyvista.RectilinearGrid
        ``clip_model`` as a PyVista grid (see :func:`model_to_rectgrid`).
    map : pyMT.gplot.MapView
        Backs the 2-D plan / X / Y panels.
    x_slice, y_slice, z_slice : int
        Current slice indices (properties; the setter also updates the
        slice labels).
    cax, rho_cax, resolution_cax, aniso_cax : list of float
        Colour-axis ranges; ``cax`` is the active one.
    use_scalar : {'Resistivity', 'RawResolution'}
    rho_axis : str
        Passed to :func:`model_to_rectgrid`.

    Notes
    -----
    The slice widgets follow one shape and are not documented
    individually: ``<axis>_text_change`` / ``<axis>_slider_change`` read
    the edit box / slider, clamp the value, set the ``<axis>_slice``
    property and call ``<axis>_slice_change``, which re-adds that slice
    actor in 3-D and refreshes the matching 2-D panel.  ``clip_x`` /
    ``clip_y`` / ``clip_z`` likewise just validate the two clip counts
    for that axis (the actual re-clip happens in :meth:`clip_volume`).
    """

    closed = QtCore.pyqtSignal(QModelWindow)

    def __init__(self, files, parent=None):
        """Load the model / data, build the 3-D and 2-D views, wire the
        widgets and draw the initial slices."""
        super(ModelWindow, self).__init__()
        self.setupUi(self)
        # This is needed to ensure the window gets deleted, and doesn't trigger a doom loop with the VTK window trying to update itself after closure
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
        # Add the model
        self.setWindowTitle('Model Viewer - {}'.format(ospath.abspath(files['model'])))
        if files.get('data', None):
            self.dataset = DS.Dataset(modelfile=files['model'],
                                        datafile=files['data'])
        elif files.get('response', None):
            self.dataset = DS.Dataset(modelfile=files['model'],
                                        datafile=files['response'])
        else:
            self.dataset = DS.Dataset(modelfile=files['model'])
            print('No data file given. Site locations will not be available.')
        # If there is a list, get it but don't load it until needed
        if files.get('list', None):
            self.list_file = files['list']
        else:
            self.list_file = None
        self.raw_data_initalized = False
        if self.dataset.model.dimensionality == '2d':
            self.dataset.data.locations[:, 1] += self.dataset.model.dy[0]
        if files.get('resolution', None):
            try:
                self.resolution = DS.Model(files['resolution'])
                self.use_resolution = 'Resolution'
            except ValueError:
                with open(files['resolution'], 'r') as f:
                    r = [float(x.strip()) for x in f.readlines()]
                    self.resolution = deepcopy(self.dataset.model)
                    self.resolution.vals = np.transpose(np.reshape(r, [self.resolution.nz,
                                                                       self.resolution.ny,
                                                                       self.resolution.nx]), [2, 1, 0])
                    self.use_resolution = 'Resolution'
            except KeyError:
                self.resolution = None
        else:
            self.resolution = None
        # NOTE: this unconditionally overrides the 'Resolution' value set
        # above, so a resolution file has no effect until the user ticks
        # the resolution checkbox.  Left as-is (changing it would alter
        # first-render output).
        self.use_resolution = 1
        self.use_scalar = 'Resistivity'
        self.rho_axis = 'rho_x'
        self.model = self.dataset.model
        self.clip_model = deepcopy(self.model)
        self.clip_resolution = deepcopy(self.resolution)
        # Make sure the frame fills the tab
        vlayout3D = QtWidgets.QVBoxLayout()
        self.vtk_widget = QtInteractor(self.frame3D)
        self.vtk_widget.background_color = 'paraview'
        vlayout3D.addWidget(self.vtk_widget)
        self.frame3D.setLayout(vlayout3D)
        self.pv_default_background = 'paraview'
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        self.colourMenu = classes.ColourMenu(self)
        mainMenu.addMenu(self.colourMenu)
        self.coordinateMenu = classes.CoordinateMenu(self)
        mainMenu.addMenu(self.coordinateMenu)
        self.spatial_units = self.coordinateMenu.default_units
        self.coordinate_system = self.coordinateMenu.default_system
        self.dataset.spatial_units = self.spatial_units
        self.model.spatial_units = self.spatial_units
        self.clip_model.spatial_units = self.spatial_units
        self.exitButton = QtWidgets.QAction('Exit', self)
        self.exitButton.setShortcut('Ctrl+Q')
        # Trigger close instead of closeEvent.
        # Triggering close calls closeEvent but with the proper input (QCloseEvent)
        self.exitButton.triggered.connect(self.close)
        fileMenu.addAction(self.exitButton)
        self.mesh_group = QtWidgets.QActionGroup(self)
        self.mesh_group.addAction(self.meshOn)
        self.mesh_group.addAction(self.meshOff)
        self.mesh_group.setExclusive(True)

        self.plot_data = {'stations': []}
        has_data = files.get('data', None)
        has_response = files.get('response', None)
        if has_data or has_response:
            self.plot_locations = True
            self.locs_3D = np.zeros((len(self.dataset.data.locations), 3))
            self.plot_data['stations'] = pv.PolyData(self.locs_3D)
            self.locs_3D[:, 0] = self.dataset.data.locations[:, 1]
            self.locs_3D[:, 1] = self.dataset.data.locations[:, 0]
        else:
            self.locs_3D = np.zeros((1, 3))
            self.plot_locations = False
        self.rect_grid = model_to_rectgrid(self.clip_model, self.clip_resolution, rho_axis=self.rho_axis)
        self.lut = 32
        # self.colourmap = 'turbo_r'
        self.colourmap = self.colourMenu.action_group.checkedAction().text()
        self.cmap = cm.get_cmap(self.colourmap, self.lut, self.colourMenu.map.invert_cmap.isChecked())
        self.cax, self.rho_cax, self.resolution_cax, self.aniso_cax = [1, 5], [1, 5], [-6, -2], [-2, 2]
        self.actors = {'X': [], 'Y': [], 'Z': [], 'transect': [], 'isosurface': []}
        self.contours = []
        self.slices = {'X': [], 'Y': [], 'Z': [], 'transect': []}
        self.transect_plot = {'easting': [], 'northing': []}
        self.transect_picking = False
        self.n_interp = 200
        self.show_transect_markers = True
        self.show_transect_line = True
        self.transect_artists = {'line': [], 'markers': []}
        self.faces = []
        self.init_widgets()
        # self.clip_volume()
        self.connect_widgets()
        self.render_3D()
        self.view_xy()
        # self.vtk_widget.view_isometric()
        # self.vtk_widget.view_xy()

        # Set up 2-D views
        self.canvas = {'2D': [], 'transect': []}
        self.toolbar = {'2D': [], 'transect': []}
        self.fig_2D = Figure()
        self.spec = {'X': [], 'Y': [], 'Z': []}
        self.spec['Y'] = GridSpec(12, 12).new_subplotspec((0, 8), colspan=4, rowspan=8)
        self.spec['X'] = GridSpec(12, 12).new_subplotspec((8, 0), colspan=8, rowspan=4)
        self.spec['Z'] = GridSpec(12, 12).new_subplotspec((0, 0), colspan=8, rowspan=8)
        self.spec['colorbar'] = GridSpec(12, 12).new_subplotspec((11, 9), colspan=3, rowspan=1)
        self.fig_2D.add_subplot(self.spec['Z'])
        self.fig_2D.add_subplot(self.spec['Y'])
        self.fig_2D.add_subplot(self.spec['X'])
        self.fig_2D.add_subplot(self.spec['colorbar'])
        self.fig_2D.subplots_adjust(top=0.95, bottom=0.1,
                                    left=0.06, right=0.94,
                                    hspace=0.08, wspace=0.02)
        self.map = gplot.MapView(figure=self.fig_2D)
        # self.map.window['colorbar'] = self.
        self.init_map(self.dataset)
        self.add_mpl(self.map.window['figure'], '2D')
        self.fig_transect = Figure()
        self.add_mpl(self.fig_transect, 'transect')
        self.fig_transect.add_subplot(111)
        # Click event hooks
        self.cid = {'transect_2D': []}

        self.update_plan_view()
        self.update_2D_X()
        self.update_2D_Y()
        self.update_2D_colorbar()

    def closeEvent(self, event):
        """Finalize the VTK widget on close so it doesn't keep trying to
        redraw itself (works for File->Exit, Ctrl+Q and the window's X,
        and for closing / reopening multiple windows)."""
        super().closeEvent(event)
        self.vtk_widget.Finalize()
        self.close()

    @property
    def x_slice(self):
        """int: Current X (northing) slice index.  Setting it also
        refreshes the ``xSliceLabel`` text.  ``y_slice`` / ``z_slice``
        behave the same for the other two axes."""
        return self._x_slice

    @x_slice.setter
    def x_slice(self, value):
        self._x_slice = value
        self.xSliceLabel.setText('{:<4.4g} {}'.format(self.model.dx[self.x_slice],
                                                      self.dataset.spatial_units))

    @property
    def y_slice(self):
        return self._y_slice

    @y_slice.setter
    def y_slice(self, value):
        self._y_slice = value
        self.ySliceLabel.setText('{:<4.4g} {}'.format(self.model.dy[self.y_slice],
                                                      self.dataset.spatial_units))

    @property
    def z_slice(self):
        return self._z_slice

    @z_slice.setter
    def z_slice(self, value):
        self._z_slice = value
        self.zSliceLabel.setText('{:<4.4g} {}'.format(self.model.dz[self.z_slice],
                                                      self.dataset.spatial_units))

    def add_mpl(self, fig, tab='2D'):
        """Embed ``fig`` (canvas + navigation toolbar) in the ``'2D'`` or
        ``'transect'`` tab layout."""
        self.canvas[tab] = FigureCanvas(fig)
        if tab == '2D':
            frame = self.widget2D
        else:
            frame = self.widgetTransect
        frame.addWidget(self.canvas[tab])
        self.toolbar[tab] = NavigationToolbar(canvas=self.canvas[tab],
                                              parent=self.canvas[tab], coordinates=True)
        self.canvas[tab].draw()
        frame.addWidget(self.toolbar[tab])

    def init_map(self, dataset):
        """Load ``dataset`` into ``self.map`` and compute the initial site
        locations for the 2-D plan panel."""
        self.map.dataset = dataset
        self.map.data = dataset.data
        self.map.raw_data = dataset.raw_data
        self.map.response = dataset.response
        self.map.site_names = dataset.data.site_names
        self.map.model = self.model
        self.map.colourmap = self.colourmap
        # self.map.model_cax = self.cax
        self.map._active_sites = []
        self.map.site_locations['generic'] = self.map.get_locations(
            sites=self.map.generic_sites)
        self.map.site_locations['active'] = self.map.get_locations(
            sites=self.map.active_sites)
        self.map.site_locations['all'] = self.map.get_locations(self.map.site_names)

    def update_2D_X(self):
        """Redraw the bottom (X / northing) 2-D panel for the current
        ``x_slice`` via :meth:`pyMT.gplot.MapView.plot_x_slice`.
        :meth:`update_2D_Y` is the equivalent for the right-hand panel."""
        self.map.window['axes'][2].clear()
        self.map.plot_x_slice(x_slice=self.x_slice,
                              ax=self.map.window['axes'][2],
                              rho_axis=self.rho_axis)
        bounds = np.array([self.clip_model.dy[0], self.clip_model.dy[-1],
                           self.clip_model.dz[0], self.clip_model.dz[-1]])
        self.map.set_axis_limits(ax=self.map.window['axes'][2], bounds=bounds)
        if not self.map.window['axes'][2].yaxis_inverted():
            self.map.window['axes'][2].invert_yaxis()
        self.canvas['2D'].draw()

    def update_2D_Y(self):
        self.map.window['axes'][1].clear()
        self.map.plot_y_slice(y_slice=self.y_slice,
                              ax=self.map.window['axes'][1],
                              orientation='zx',
                              rho_axis=self.rho_axis)
        bounds = np.array([self.clip_model.dz[0], self.clip_model.dz[-1],
                           self.clip_model.dx[0], self.clip_model.dx[-1]])
        self.map.set_axis_limits(ax=self.map.window['axes'][1], bounds=bounds)
        self.map.window['axes'][1].yaxis.tick_right()
        self.map.window['axes'][1].yaxis.set_label_position('right')
        self.canvas['2D'].draw()

    def update_plan_view(self):
        """Redraw the plan (Z / depth) 2-D panel for the current
        ``z_slice``: the model layer via
        :meth:`pyMT.gplot.MapView.plot_plan_view`, the X / Y slice
        indicator lines, station markers, the transect line / markers,
        and the depth title.  Full redraw every call."""
        self.map.window['axes'][0].clear()
        if self.map.window['colorbar']:
            self.map.window['colorbar'].remove()
            self.map.window['colorbar'] = None
        # self.toolbar.push_current()
        bounds = np.array([self.clip_model.dy[0], self.clip_model.dy[-1],
                           self.clip_model.dx[0], self.clip_model.dx[-1]])
        self.map.plot_plan_view(z_slice=self.z_slice, ax=self.map.window['axes'][0], rho_axis=self.rho_axis)
        self.update_slice_indicators('xy')
        if self.plot_locations:
            self.map.plot_locations()
            # print(self.map.data.locations)
        self.map.set_axis_limits(bounds=bounds)
        self.map.window['axes'][0].get_xaxis().set_visible(False)
        self.plot_transect_markers(redraw=True)
        self.plot_transect_line(redraw=True)
        self.update_plan_title()
        # self.map.window['axes'][0].set_aspect('equal')
        self.canvas['2D'].draw()

    def update_2D_colorbar(self):
        """Redraw the shared horizontal colour bar (4th axis of the 2-D
        figure) for the current ``cax`` / ``cmap`` / ``lut``."""
        self.map.window['axes'][3].clear()
        fake_vals = np.linspace(self.cax[0], self.cax[1], self.lut)
        fake_im = self.map.window['axes'][3].scatter(range(self.lut), range(self.lut), s=10, marker='.',
                                                     c=fake_vals, cmap=self.cmap)
        fake_im.set_visible(False)
        self.map.window['figure'].colorbar(mappable=fake_im, cax=self.map.window['axes'][3],
                                           orientation='horizontal')
        self.map.window['axes'][3].set_xlabel(r'$Log_{10}$ Resistivity ($\Omega$m)', fontsize=self.map.label_fontsize)
        self.canvas['2D'].draw()

    def update_plan_title(self):
        """Set the plan-panel title to the layer number and its centre
        depth (``depthLabelCenter``) or depth range (``depthLabelRange``)."""
        if self.depthLabelCenter.isChecked():
            self.map.window['axes'][0].set_title('Layer ' + str(self.z_slice) +
                                                 ', Depth = '+'{:.3g}'.format((self.model.dz[self.z_slice]+self.model.dz[self.z_slice-1])/2)+' {}'.format(self.spatial_units), fontsize=self.map.label_fontsize)
        elif self.depthLabelRange.isChecked():
            self.map.window['axes'][0].set_title('Layer ' + str(self.z_slice) +
                                                 ', Depth = '+'{:.3g}'.format(self.model.dz[self.z_slice-1])+u'\u2013'+'{:.3g} {}'.format(self.model.dz[self.z_slice], self.spatial_units),
                                                  fontsize=self.map.label_fontsize)
        self.canvas['2D'].draw()

    def update_slice_indicators(self, direction='xy', redraw=True):
        """Draw or move the white dashed lines on the plan panel that mark
        where the X and Y slices cut.  ``direction`` selects which
        line(s); ``redraw=False`` just moves existing lines instead of
        re-creating them."""
        x_slice_loc = self.model.dx[self.x_slice]
        y_slice_loc = self.model.dy[self.y_slice]
        x_slice_limits = (self.clip_model.dy[0], self.clip_model.dy[-1])
        y_slice_limits = (self.clip_model.dx[0], self.clip_model.dx[-1])
        if 'y' in direction.lower() and redraw:
            self.plot_data['Y_indicator'], = self.map.window['axes'][0].plot((y_slice_loc, y_slice_loc),
                                                                             y_slice_limits,
                                                                             'w--')
        if 'x' in direction.lower() and redraw:
            self.plot_data['X_indicator'], = self.map.window['axes'][0].plot(x_slice_limits,
                                                                             (x_slice_loc, x_slice_loc),
                                                                             'w--')
        if 'x' in direction.lower() and not redraw:
            self.plot_data['X_indicator'].set_ydata((x_slice_loc, x_slice_loc))
        if 'y' in direction.lower() and not redraw:
            self.plot_data['Y_indicator'].set_xdata((y_slice_loc, y_slice_loc))

    def init_widgets(self):
        """Set the initial state of every slice / clip / interpolation
        widget from the model dimensions, and enable the resolution /
        elevation checkboxes only when that data is available."""
        self.x_slice = 1
        self.y_slice = 1
        self.z_slice = 1
        self.xSliceEdit.setText('1')
        self.ySliceEdit.setText('1')
        self.zSliceEdit.setText('1')
        self.xSliceLabel.setText('1')
        self.ySliceLabel.setText('1')
        self.zSliceLabel.setText('1')
        self.xSliceSlider.setValue(1)
        self.ySliceSlider.setValue(1)
        self.zSliceSlider.setValue(1)
        self.xSliceSlider.setSingleStep(1)
        self.ySliceSlider.setSingleStep(1)
        self.zSliceSlider.setSingleStep(1)
        self.xSliceSlider.setPageStep(5)
        self.ySliceSlider.setPageStep(5)
        self.zSliceSlider.setPageStep(5)
        self.xSliceSlider.setMinimum(1)
        self.ySliceSlider.setMinimum(1)
        self.zSliceSlider.setMinimum(1)
        self.xSliceSlider.setMaximum(self.model.nx - 1)
        self.ySliceSlider.setMaximum(self.model.ny - 1)
        self.zSliceSlider.setMaximum(self.model.nz - 1)
        self.xSliceCheckbox.setCheckState(2)
        self.ySliceCheckbox.setCheckState(2)
        self.zSliceCheckbox.setCheckState(2)
        # Init clip boxes
        self.x_clip = [0, 0]
        self.y_clip = [0, 0]
        self.z_clip = [0, 0]
        self.x0ClipEdit.setText('0')
        self.y0ClipEdit.setText('0')
        self.z0ClipEdit.setText('0')
        self.x1ClipEdit.setText('0')
        self.y1ClipEdit.setText('0')
        self.z1ClipEdit.setText('0')
        self.nInterp.setText(str(self.n_interp))
        if not self.resolution:
            self.resolutionCheckBox.setEnabled(False)
        if self.plot_locations:
            if np.any([site.locations['elev'] for site in self.dataset.data.sites.values()]):
                self.elevationCheckBox.setEnabled(True)


    def connect_widgets(self):
        """Wire every slider / edit box / checkbox / menu action to its
        handler.  Coordinate-system options beyond ``local`` are disabled
        unless a list file was supplied; the rho-axis group is enabled
        only for anisotropic models."""
        self.xSliceEdit.editingFinished.connect(self.x_text_change)
        self.ySliceEdit.editingFinished.connect(self.y_text_change)
        self.zSliceEdit.editingFinished.connect(self.z_text_change)
        self.xSliceSlider.valueChanged.connect(self.x_slider_change)
        self.ySliceSlider.valueChanged.connect(self.y_slider_change)
        self.zSliceSlider.valueChanged.connect(self.z_slider_change)
        # Clipping widgets
        self.x0ClipEdit.editingFinished.connect(self.clip_x)
        self.x1ClipEdit.editingFinished.connect(self.clip_x)
        self.y0ClipEdit.editingFinished.connect(self.clip_y)
        self.y1ClipEdit.editingFinished.connect(self.clip_y)
        self.z0ClipEdit.editingFinished.connect(self.clip_z)
        self.z1ClipEdit.editingFinished.connect(self.clip_z)
        self.clipVolume.pressed.connect(self.clip_volume)
        # Slice Checkboxes
        self.xSliceCheckbox.stateChanged.connect(self.x_slice_change)
        self.ySliceCheckbox.stateChanged.connect(self.y_slice_change)
        self.zSliceCheckbox.stateChanged.connect(self.z_slice_change)
        # Colour options
        self.colourMenu.action_group.triggered.connect(self.change_cmap)
        self.colourMenu.limits.triggered.connect(self.set_clim)
        self.colourMenu.lut.triggered.connect(self.set_lut)
        self.colourMenu.map.invert_cmap.triggered.connect(self.change_cmap)
        # self.colourMenu.all_maps[self.colourmap].setChecked(True)
        # Coordinate options
        self.coordinateMenu.systemGroup.triggered.connect(self.load_raw_data)
        self.coordinateMenu.unitGroup.triggered.connect(self.set_spatial_units)
        self.coordinateMenu.all_systems['LatLong'].setEnabled(False)
        if not self.list_file:
            self.coordinateMenu.all_systems['Lambert'].setEnabled(False)
            self.coordinateMenu.all_systems['UTM'].setEnabled(False)
        # 2-D Plot options
        self.mesh_group.triggered.connect(self.show_mesh)
        self.depthTitles = QtWidgets.QActionGroup(self)
        self.depthTitles.addAction(self.depthLabelRange)
        self.depthTitles.addAction(self.depthLabelCenter)
        self.depthTitles.triggered.connect(self.update_plan_title)
        # View buttons
        self.viewXY.triggered.connect(self.view_xy)
        self.viewXZ.triggered.connect(self.view_xz)
        self.viewYZ.triggered.connect(self.view_yz)
        # Interpolation widgets
        self.selectPoints.clicked.connect(self.transect_pick)
        self.nInterp.editingFinished.connect(self.set_nInterp)
        self.interpCheckbox.clicked.connect(self.toggle_transect3D)
        # Isosurface widgets
        self.isoPlot.stateChanged.connect(self.update_isosurface)
        self.isoRecalculate.clicked.connect(self.generate_isosurface)
        self.isoOpacity.valueChanged.connect(self.update_isosurface)
        # Resolution toggles
        self.resolutionCheckBox.clicked.connect(self.set_resolution_plot)
        # Station location toggles
        self.elevationCheckBox.clicked.connect(self.set_station_elevation)
        # PyVista backgrounds
        self.pv_background_group = QtWidgets.QActionGroup(self)
        self.pv_background_group.addAction(self.actionDefault)
        self.pv_background_group.addAction(self.actionGrey)
        self.pv_background_group.addAction(self.actionBlack)
        self.pv_background_group.addAction(self.actionBlue)
        self.pv_background_group.addAction(self.actionWhite)
        self.pv_background_group.setExclusive(True)
        self.pv_background_group.triggered.connect(self.change_background)
        # Rho axis selection actions
        self.rho_axis_group = QtWidgets.QActionGroup(self)
        self.rho_axis_group.addAction(self.action_plot_rho_x)
        self.rho_axis_group.addAction(self.action_plot_rho_y)
        self.rho_axis_group.addAction(self.action_plot_rho_z)
        self.rho_axis_group.addAction(self.action_plot_rho_xy)
        self.rho_axis_group.addAction(self.action_plot_rho_xz)
        self.rho_axis_group.addAction(self.action_plot_rho_yz)
        self.rho_axis_group.setExclusive(True)

        if self.dataset.model.isotropic:
            self.rho_axis_group.setEnabled(False)
        else:
            self.rho_axis_group.triggered.connect(self.change_rho_axis)

    def set_spatial_units(self):
        """Switch the model / dataset / grid between metres and
        kilometres, rescale the stored camera position to match, and
        redraw everything."""
        if self.spatial_units == self.coordinateMenu.unitGroup.checkedAction().text().lower():
            return
        self.spatial_units = self.coordinateMenu.unitGroup.checkedAction().text().lower()
        self.dataset.spatial_units = self.spatial_units
        self.model.spatial_units = self.spatial_units
        self.clip_model.spatial_units = self.spatial_units
        self.rect_grid = model_to_rectgrid(self.clip_model, self.clip_resolution, rho_axis=self.rho_axis)
        self.map.set_locations()
        self.locs_3D[:, 0] = self.map.site_locations['all'][:, 1]
        self.locs_3D[:, 1] = self.map.site_locations['all'][:, 0]
        camera = self.vtk_widget.camera_position
        if self.spatial_units == 'km':
            camera = [[x / 1000 for x in p1] for p1 in camera]        
        else:
            camera = [[x * 1000 for x in p1] for p1 in camera]        
        self.update_all(reset_camera=False)
        self.vtk_widget.camera_position = camera


    def load_raw_data(self):
        """Coordinate-system menu handler: lazily load the raw EDI data
        (needed for UTM / lat-long / Lambert) via a progress popup the
        first time, then :meth:`set_coordinate_system`."""
        if not self.raw_data_initalized:
            if self.list_file:
                try:
                    self.popup = classes.PopUpProgress(file_dict={'list': self.list_file})
                    self.popup.finished.connect(self.collect_from_popup)
                    self.popup.start_progress()

                except FileNotFoundError:
                    QtWidgets.QMessageBox.warning(self, '...', 'List file not found - Coordinate system cannot be changed.')
        else:
            self.set_coordinate_system()

    def set_coordinate_system(self):
        """Reproject the model and station locations into the selected
        coordinate system (local / UTM / lat-long / Lambert), re-clip and
        redraw.  Warns and reverts to local if the target can't be
        verified."""
        can_do = 1
        coord_system = self.coordinateMenu.systemGroup.checkedAction().text().lower()
        if coord_system == 'local':
            if self.map.verify_coordinate_system('local'):
                self.map.coordinate_system = 'local'
        elif coord_system == 'utm':
            if self.map.verify_coordinate_system('utm'):
                self.map.coordinate_system = 'utm'
            else:
                can_do = 0
        elif coord_system == 'latlong':
            if self.map.verify_coordinate_system('latlong'):
                self.map.coordinate_system = 'latlong'
            else:
                can_do = 0
        elif coord_system == 'lambert':
            if self.map.verify_coordinate_system('lambert'):
                self.map.coordinate_system = 'lambert'
            else:
                can_do = 0
        if can_do:
            self.model = self.map.model
            self.map.set_locations()
            self.locs_3D[:, 0] = self.map.site_locations['all'][:, 1]
            self.locs_3D[:, 1] = self.map.site_locations['all'][:, 0]
            self.clip_volume()
            self.vtk_widget.reset_camera()
        else:
            QtWidgets.QMessageBox.warning(self, '...', 'Unable to change coordinate system.')
            self.coordinateMenu.all_systems['Local'].setChecked(True)

    def collect_from_popup(self):
        """Progress-popup callback: store the loaded ``raw_data`` and
        continue into :meth:`set_coordinate_system`."""
        self.dataset.raw_data = self.popup.ret_val['raw_data']
        self.raw_data_initalized = True
        self.map.site_data['raw_data'] = self.dataset.raw_data
        self.set_coordinate_system()

    def change_rho_axis(self):
        """Switch which model attribute (or ratio) is displayed
        (anisotropic models), pick the matching colour axis, rebuild the
        grid and redraw."""
        self.rho_axis = self.rho_axis_group.checkedAction().text().lower()
        if '/' in self.rho_axis:
            self.cax = self.aniso_cax
        else:
            self.cax = self.rho_cax
        self.rect_grid = model_to_rectgrid(self.clip_model, self.clip_resolution, rho_axis=self.rho_axis)
        self.update_all(reset_camera=False)

    def change_background(self):
        """Set the PyVista render background from the menu selection."""
        color = self.pv_background_group.checkedAction().text().lower()
        if color != 'default':
            self.vtk_widget.background_color = color
        else:
            self.vtk_widget.background_color = deepcopy(self.pv_default_background)

    def set_station_elevation(self):
        """Toggle whether station markers sit at their true elevation
        (from the data) or on the z = 0 plane, and redraw the 3-D view."""
        if self.elevationCheckBox.checkState():
            self.locs_3D[:, 2] = np.array([-self.dataset.data.sites[site].locations['elev'] / 1000 for site in self.dataset.data.site_names])
        else:
            self.locs_3D[:, 2] = 0
        self.render_3D(reset_camera=False)

    def set_resolution_plot(self):
        """Tri-state resolution checkbox: 0 = plain resistivity,
        1 = resistivity with resolution as opacity, 2 = raw resolution as
        the scalar.  Sets ``use_scalar`` / ``use_resolution`` / ``cax``
        accordingly and redraws the 3-D view."""
        if self.resolutionCheckBox.checkState() == 0:
            self.use_scalar = 'Resistivity'
            self.use_resolution = 1
            self.cax = self.rho_cax
        elif self.resolutionCheckBox.checkState() == 1:
            self.use_scalar = 'Resistivity'
            self.use_resolution = 'Resolution'
            self.cax = self.rho_cax
        elif self.resolutionCheckBox.checkState() == 2:
            self.use_scalar = 'RawResolution'
            self.use_resolution = 1
            self.cax = self.resolution_cax
        self.render_3D(reset_camera=False)

    def generate_isosurface(self):
        """Compute resistivity isosurfaces between the low / high spin-box
        values (``isoNum`` of them) and show them (:meth:`update_isosurface`)."""
        low_val = self.isoLow.value()
        high_val = self.isoHigh.value()
        if not (low_val <= high_val):
            low_val = high_val
            self.isoLow.setValue(high_val)
        if low_val == high_val:
            self.isoNum.setValue(1)
        self.contours = pv.DataSetFilters.cell_data_to_point_data(self.rect_grid).contour(isosurfaces=self.isoNum.value(),
                                                                  rng=(np.log10(low_val), np.log10(high_val)))
        # self.actors['isosurface'] = self.vtk_widget.add_points(self.contours,
        #                                                      cmap=self.cmap,
        #                                                      clim=self.cax)
        self.isoPlot.setCheckState(2)
        self.update_isosurface()

    def update_isosurface(self):
        """Add / update the isosurface actor from ``self.contours`` at the
        current opacity (or unchecks ``isoPlot`` if there are no
        contours)."""
        if self.isoPlot.checkState():
            if self.contours:
                opacity = self.isoOpacity.value()
                try:
                    self.actors['isosurface'] = self.vtk_widget.add_mesh(self.contours,
                                                                         cmap=self.cmap,
                                                                         clim=self.cax,
                                                                         opacity=opacity,
                                                                         name='isosurface')
                except RuntimeError:
                    self.isoPlot.setCheckState(0)
                    QtWidgets.QMessageBox.warning(self, '...', 'No model values fit the requested contour.')
                    return
            else:
                self.isoPlot.setCheckState(0)
        self.show_bounds()

    def show_mesh(self):
        """Toggle the cell-edge mesh overlay on the 2-D panels and
        redraw."""
        if self.meshOn.isChecked():
            self.map.mesh = True
        else:
            self.map.mesh = False
        self.update_all()

    def set_nInterp(self):
        """Apply the transect interpolation-point count from ``nInterp``
        and regenerate the transect if one is picked."""
        try:
            N = int(self.nInterp.text())
            self.n_interp = N
            if len(self.transect_plot['easting']) > 1:
                self.generate_transect2D()
                self.generate_transect3D()
                # self.plot_transect()
        except ValueError:
            self.nInterp.setText(str(self.n_interp))
            self.debugLabel.setText('# of interp points must be an integer')

    def transect_pick(self):
        """"Select points" button: on the first press, connect the plan
        panel click handler (:meth:`click_2D`) so the user can lay down
        transect vertices; on the second press, disconnect it and build
        the 2-D / 3-D transect from the picked points."""
        if not self.cid['transect_2D']:
            self.cid['transect_2D'] = self.canvas['2D'].mpl_connect('button_release_event', self.click_2D)
            self.tabWidget_2.setCurrentIndex(1)
            self.selectPoints.setStyleSheet('QPushButton {background-color: lightgreen; color: black;}')
            self.transect_plot['easting'] = []
            self.transect_plot['northing'] = []
            self.debugLabel.setText('Picking turned on.')
        else:
            self.canvas['2D'].mpl_disconnect(self.cid['transect_2D'])
            self.selectPoints.setStyleSheet('QPushButton {background-color: None; color: black;}')
            if len(self.transect_plot['easting']) > 1:
                self.generate_transect2D()
                self.generate_transect3D()
            self.cid['transect_2D'] = 0
            self.interpCheckbox.setCheckState(2)
            self.toggle_transect3D()
            self.tabWidget_2.setCurrentIndex(0)

    def save_picks(self, picks):
        """Callback for PyVista cell-picking (currently unused -- the
        commented-out ``enable_cell_picking`` path in
        :meth:`transect_pick`)."""
        self.picks.append(picks)
        print('Saving picks')

    def change_cmap(self):
        """Rebuild the colour map from the Colour menu (name + invert +
        LUT) and redraw both views."""
        self.colourmap = self.colourMenu.action_group.checkedAction().text()
        self.cmap = cm.get_cmap(self.colourMenu.action_group.checkedAction().text(),
                                self.lut,
                                self.colourMenu.map.invert_cmap.isChecked())
        self.map.colourmap = self.colourmap
        self.map.lut = self.lut
        self.update_all(reset_camera=False)

    def set_lut(self):
        """Prompt for the number of colour levels and apply it via
        :meth:`change_cmap`."""
        inputs, ret = self.colourMenu.set_lut(initial=self.lut)
        if ret:
            self.lut = inputs
            self.change_cmap()

    def set_clim(self):
        """Prompt for the ``[lower, upper]`` colour-axis limits and apply
        them to both views."""
        inputs, ret = self.colourMenu.set_clim(initial_1=str(self.cax[0]),
                                               initial_2=str(self.cax[1]))
        lower, upper = inputs
        if ret and [lower, upper] != self.cax:
            try:
                self.cax[0] = float(lower)
                self.cax[1] = float(upper)
                self.map.model_cax = self.cax
                self.update_all()
            except ValueError:
                print('Bad inputs to clim')

    def clip_x(self):
        """Validate the two X clip-count edit boxes (0/1 = layers trimmed
        from each end); reject if they exceed ``nx``.  ``clip_y`` /
        ``clip_z`` do the same for the other axes.  The trim itself
        happens in :meth:`clip_volume`."""
        try:
            x0 = int(self.x0ClipEdit.text())
            x1 = int(self.x1ClipEdit.text())
        except ValueError:
            self.x0ClipEdit.setText(str(self.x_clip[0]))
            self.x1ClipEdit.setText(str(self.x_clip[1]))
            x0 = int(self.x0ClipEdit.text())
            x1 = int(self.x1ClipEdit.text())
        if x0 + x1 > self.model.nx:
            self.x0ClipEdit.setText(str(self.x_clip[0]))
            self.x1ClipEdit.setText(str(self.x_clip[1]))
        else:
            self.x_clip[0] = x0
            self.x_clip[1] = x1

    def clip_y(self):
        try:
            y0 = int(self.y0ClipEdit.text())
            y1 = int(self.y1ClipEdit.text())
        except ValueError:
            self.y0ClipEdit.setText(str(self.y_clip[0]))
            self.y1ClipEdit.setText(str(self.y_clip[1]))
            y0 = int(self.y0ClipEdit.text())
            y1 = int(self.y1ClipEdit.text())
        if y0 + y1 > self.model.ny:
            self.y0ClipEdit.setText(str(self.y_clip[0]))
            self.y1ClipEdit.setText(str(self.y_clip[1]))
        else:
            self.y_clip[0] = y0
            self.y_clip[1] = y1

    def clip_z(self):
        try:
            z0 = int(self.z0ClipEdit.text())
            z1 = int(self.z1ClipEdit.text())
        except ValueError:
            self.z0ClipEdit.setText(str(self.z_clip[0]))
            self.z1ClipEdit.setText(str(self.z_clip[1]))
            z0 = int(self.z0ClipEdit.text())
            z1 = int(self.z1ClipEdit.text())
        if z0 + z1 > self.model.nz:
            self.z0ClipEdit.setText(str(self.z_clip[0]))
            self.z1ClipEdit.setText(str(self.z_clip[1]))
        else:
            self.z_clip[0] = z0
            self.z_clip[1] = z1

    def x_text_change(self):
        try:
            x = int(self.xSliceEdit.text())
        except ValueError:
            self.xSliceEdit.setText(str(self.x_slice))
            return
        if x < 1 or x > self.model.nx:
            self.xSliceEdit.setText(str(self.x_slice))
        else:
            self.x_slice = x
            self.xSliceSlider.setValue(x)
            self.x_slice_change()

    def y_text_change(self):
        try:
            y = int(self.ySliceEdit.text())
        except ValueError:
            self.ySliceEdit.setText(str(self.y_slice))
            return
        if y < 1 or y > self.model.ny:
            self.ySliceEdit.setText(str(self.y_slice))
        else:
            self.y_slice = y
            self.ySliceSlider.setValue(y)
            self.y_slice_change()

    def z_text_change(self):
        try:
            z = int(self.zSliceEdit.text())
        except ValueError:
            self.zSliceEdit.setText(str(self.z_slice))
            return
        if z < 1 or z > self.model.nz:
            self.zSliceEdit.setText(str(self.z_slice))
        else:
            self.z_slice = z
            self.zSliceSlider.setValue(z)
            self.z_slice_change()

    def x_slider_change(self, value):
        self.xSliceEdit.setText(str(value))
        self.x_slice = value
        self.x_slice_change()

    def y_slider_change(self, value):
        self.ySliceEdit.setText(str(value))
        self.y_slice = value
        self.y_slice_change()

    def z_slider_change(self, value):
        self.zSliceEdit.setText(str(value))
        self.z_slice = value
        self.z_slice_change()

    def x_slice_change(self):
        """Re-take the X slice: add / remove its 3-D actor per the
        ``xSliceCheckbox`` and refresh the 2-D X panel + the plan-view
        indicator line.  ``y_slice_change`` / ``z_slice_change`` are the
        equivalents for the other axes (``z`` refreshes the plan panel)."""
        if self.xSliceCheckbox.checkState():
            self.slices['X'] = self.generate_slice('x')
            self.actors['X'] = self.vtk_widget.add_mesh(self.slices['X'],
                                                        cmap=self.cmap,
                                                        clim=self.cax,
                                                        scalars=self.use_scalar,
                                                        opacity=self.use_resolution,
                                                        use_transparency=False,
                                                        name='x_slice')
        else:
            self.vtk_widget.remove_actor(self.actors['X'])
        self.update_2D_X()
        self.update_slice_indicators(direction='x', redraw=False)
        self.fig_2D.canvas.draw()
        # self.vtk_widget.update()
        # self.show_bounds()

    def y_slice_change(self):
        # self.vtk_widget.remove_actor(self.actors['Y'])
        if self.ySliceCheckbox.checkState():
            # y_slice = self.clip_model.dy[self.y_slice] 
            # self.actors['Y'] = self.rect_grid.slice(normal='x',
            #                                         origin=(y_slice, 0, 0))
            self.slices['Y'] = self.generate_slice('y')
            self.actors['Y'] = self.vtk_widget.add_mesh(self.slices['Y'],
                                                        cmap=self.cmap,
                                                        clim=self.cax,
                                                        scalars=self.use_scalar,
                                                        opacity=self.use_resolution,
                                                        use_transparency=False,
                                                        name='y_slice')
        else:
            self.vtk_widget.remove_actor(self.actors['Y'])
        self.update_2D_Y()
        self.update_slice_indicators(direction='y', redraw=False)
        self.fig_2D.canvas.draw()
        self.vtk_widget.update()
        # self.show_bounds()

    def z_slice_change(self):
        # self.vtk_widget.remove_actor(self.actors['Z'])
        if self.zSliceCheckbox.checkState():
            # z_slice = self.clip_model.dz[self.z_slice] 
            # self.actors['Z'] = self.rect_grid.slice(normal='z', origin=(0, 0, z_slice))
            self.slices['Z'] = self.generate_slice('z')
            self.actors['Z'] = self.vtk_widget.add_mesh(self.slices['Z'],
                                                        cmap=self.cmap,
                                                        clim=self.cax,
                                                        scalars=self.use_scalar,
                                                        opacity=self.use_resolution,
                                                        use_transparency=False,
                                                        name='z_slice')
        else:
            self.vtk_widget.remove_actor(self.actors['Z'])
        self.update_plan_view()
        # self.vtk_widget.update()
        # self.show_bounds()

    def generate_slice(self, normal='X'):
        """Return the ``rect_grid`` slice normal to ``'x'`` / ``'y'`` /
        ``'z'`` at the current slice index, its origin offset to the
        clipped-volume corner."""
        ox = max(0, self.model.dx[self.x_clip[0]])
        oy = max(0, self.model.dy[self.y_clip[0]])
        if normal == 'x':
            slice_loc = self.model.dx[self.x_slice]
            gen_slice = self.rect_grid.slice(normal='y',
                                             origin=(oy, slice_loc, -self.model.dz[self.z_clip[0] + 1]))
        elif normal == 'y':
            slice_loc = self.model.dy[self.y_slice]
            gen_slice = self.rect_grid.slice(normal='x',
                                             origin=(slice_loc, ox, -self.model.dz[self.z_clip[0] + 1]))
        elif normal == 'z':
            slice_loc = -self.model.dz[self.z_slice]
            gen_slice = self.rect_grid.slice(normal=normal,
                                             origin=(oy, ox, slice_loc))
        return gen_slice

    def clip_volume(self):
        """Rebuild ``clip_model`` / ``clip_resolution`` / ``rect_grid``
        by deleting the requested number of layers from each end of each
        axis, then re-range the slice sliders and redraw."""
        self.clip_model = deepcopy(self.model)
        self.clip_resolution = deepcopy(self.resolution)
        for ix in range(self.x_clip[0]):
            self.clip_model.dx_delete(0)
            if self.resolution:
                self.clip_resolution.dx_delete(0)
        for ix in range(self.x_clip[1]):
            self.clip_model.dx_delete(self.clip_model.nx)
            if self.resolution:
                self.clip_resolution.dx_delete(self.clip_model.nx)
        for iy in range(self.y_clip[0]):
            self.clip_model.dy_delete(0)
            if self.resolution:
                self.clip_resolution.dy_delete(0)
        for iy in range(self.y_clip[1]):
            self.clip_model.dy_delete(self.clip_model.ny)
            if self.resolution:
                self.clip_resolution.dy_delete(self.clip_model.ny)
        for iz in range(self.z_clip[0]):
            self.clip_model.dz_delete(0)
            if self.resolution:
                self.clip_resolution.dz_delete(0)
        for iz in range(self.z_clip[1]):
            self.clip_model.dz_delete(self.clip_model.nz)
            if self.resolution:
                self.clip_resolution.dz_delete(self.clip_model.nz)
        # self.map.model = self.clip_model
        self.rect_grid = model_to_rectgrid(self.clip_model, self.clip_resolution, rho_axis=self.rho_axis)
        self.validate_slice_locations()
        self.xSliceSlider.setMinimum(self.x_clip[0] + 1)
        self.ySliceSlider.setMinimum(self.y_clip[0] + 1)
        self.zSliceSlider.setMinimum(self.z_clip[0] + 1)
        self.xSliceSlider.setMaximum(self.model.nx - self.x_clip[1] - 1)
        self.ySliceSlider.setMaximum(self.model.ny - self.y_clip[1] - 1)
        self.zSliceSlider.setMaximum(self.model.nz - self.z_clip[1] - 1)
        self.update_all()

    def validate_slice_locations(self):
        """Clamp each slice index into the clipped-volume range,
        temporarily disconnecting the slider / edit-box signals so the
        corrections don't retrigger a redraw."""
        self.xSliceSlider.valueChanged.disconnect(self.x_slider_change)
        self.ySliceSlider.valueChanged.disconnect(self.y_slider_change)
        self.zSliceSlider.valueChanged.disconnect(self.z_slider_change)
        self.xSliceEdit.editingFinished.disconnect(self.x_text_change)
        self.ySliceEdit.editingFinished.disconnect(self.y_text_change)
        self.zSliceEdit.editingFinished.disconnect(self.z_text_change)
        if self.x_slice <= self.x_clip[0]:
            self.x_slice = self.x_clip[0] + 1
            self.xSliceSlider.setValue(self.x_slice)
            self.xSliceEdit.setText(str(self.x_slice))
        if self.y_slice <= self.y_clip[0]:
            self.y_slice = self.y_clip[0] + 1
            self.ySliceSlider.setValue(self.y_slice)
            self.ySliceEdit.setText(str(self.y_slice))
        if self.z_slice <= self.z_clip[0]:
            self.z_slice = self.z_clip[0] + 1
            self.zSliceSlider.setValue(self.z_slice)
            self.zSliceEdit.setText(str(self.z_slice))
        if self.x_slice > self.model.nx - self.x_clip[1]:
            self.x_slice = self.model.nx - self.x_clip[1] - 1
            self.xSliceSlider.setValue(self.x_slice)
            self.xSliceEdit.setText(str(self.x_slice))
        if self.y_slice > self.model.ny - self.y_clip[1]:
            self.y_slice = self.model.ny - self.y_clip[1] - 1
            self.ySliceSlider.setValue(self.y_slice)
            self.ySliceEdit.setText(str(self.y_slice))
        if self.z_slice > self.model.nz - self.z_clip[1]:
            self.z_slice = self.model.nz - self.z_clip[1] - 1
            self.zSliceSlider.setValue(self.z_slice)
            self.zSliceEdit.setText(str(self.z_slice))
        self.xSliceSlider.valueChanged.connect(self.x_slider_change)
        self.ySliceSlider.valueChanged.connect(self.y_slider_change)
        self.zSliceSlider.valueChanged.connect(self.z_slider_change)
        self.xSliceEdit.editingFinished.connect(self.x_text_change)
        self.ySliceEdit.editingFinished.connect(self.y_text_change)
        self.zSliceEdit.editingFinished.connect(self.z_text_change)

    def render_3D(self, reset_camera=True):
        """Full 3-D redraw: clear the render window and re-add the enabled
        slice actors, the transect, and the station points, then show the
        axis bounds (and reset the camera unless ``reset_camera`` is
        ``False``).  Slices are regenerated even when their checkbox is
        off so they stay trimmed to the clip box."""
        self.vtk_widget.clear()
        if self.xSliceCheckbox.checkState():
            self.slices['X'] = self.generate_slice('x')
            self.actors['X'] = self.vtk_widget.add_mesh(self.slices['X'],
                                                        scalars=self.use_scalar,
                                                        cmap=self.cmap,
                                                        clim=self.cax,
                                                        opacity=self.use_resolution,
                                                        use_transparency=False,
                                                        name='x_slice')
        if self.ySliceCheckbox.checkState():
            self.slices['Y'] = self.generate_slice('y')
            self.actors['Y'] = self.vtk_widget.add_mesh(self.slices['Y'],
                                                        scalars=self.use_scalar,
                                                        cmap=self.cmap,
                                                        clim=self.cax,
                                                        opacity=self.use_resolution,
                                                        use_transparency=False,
                                                        name='y_slice')
        if self.zSliceCheckbox.checkState():
            self.slices['Z'] = self.generate_slice('z')
            self.actors['Z'] = self.vtk_widget.add_mesh(self.slices['Z'],
                                                        scalars=self.use_scalar,
                                                        cmap=self.cmap,
                                                        clim=self.cax,
                                                        opacity=self.use_resolution,
                                                        use_transparency=False,
                                                        name='z_slice')
        if len(self.transect_plot['easting']) > 1:
            self.generate_transect3D(redraw=False)
        # Checkbox for sites?
        if self.plot_locations:
            self.actors['stations'] = self.vtk_widget.add_mesh(self.plot_data['stations'],
                                                               style='points',
                                                               color='white',
                                                               name='stations')
        # self.vtk_widget.update()
        self.show_bounds()
        if reset_camera:
            self.vtk_widget.reset_camera()

    def show_bounds(self):
        """(Re)draw the labelled 3-D axes / grid box around the clipped
        model extent."""
        self.vtk_widget.show_bounds(bounds=[self.clip_model.dy[0], self.clip_model.dy[-1],
                                            self.clip_model.dx[0], self.clip_model.dx[-1],
                                           -self.clip_model.dz[-1], self.clip_model.dz[0]],
                                    axes_ranges=[round(self.clip_model.dy[0], -1), round(self.clip_model.dy[-1],-1),
                                                 round(self.clip_model.dx[0], -1), round(self.clip_model.dx[-1],-1),
                                                -round(self.clip_model.dz[-1], -1), round(self.clip_model.dz[0],-1)],
                                    xtitle='Easting ({})'.format(self.spatial_units),
                                    ytitle='Northing ({})'.format(self.spatial_units),
                                    ztitle='Depth ({})'.format(self.spatial_units),
                                    font_family='arial',
                                    grid=True,
                                    location='outer',
                                    use_3d_text=False)

    def view_xy(self):
        """Snap the 3-D camera to the top-down (XY) view.  ``view_xz`` /
        ``view_yz`` are the side views."""
        self.vtk_widget.view_xy()

    def view_xz(self):
        self.vtk_widget.view_xz()

    def view_yz(self):
        self.vtk_widget.view_yz()

    def update_all(self, reset_camera=True):
        """Redraw everything: all three 2-D panels, the colour bar and
        the 3-D view."""
        self.update_2D_X()
        self.update_2D_Y()
        self.update_plan_view()
        self.update_2D_colorbar()
        self.render_3D(reset_camera)

    def click_2D(self, event):
        """Plan-panel click handler while transect picking is on: append
        the clicked (easting, northing) to the transect vertex list."""
        if self.map.window['axes'][0] == event.inaxes:
            self.transect_plot['easting'].append(event.xdata)
            self.transect_plot['northing'].append(event.ydata)
            self.debugLabel.setText('{}'.format(self.transect_plot))
        else:
            print('Outside valid axis')
            return

    def interpolate_between_points(self, E, N):
        """Densify a polyline of ``(E, N)`` vertices to ~``n_interp``
        points, distributing them between segments by length.

        Returns
        -------
        qx, qy : numpy.ndarray
            The interpolated easting / northing coordinates.
        """
        qx = []
        qy = []
        distance = []
        for ii in range(len(E) - 1):
            distance.append(np.sqrt((E[ii] - E[ii + 1]) ** 2 +
                                    (N[ii] - N[ii + 1]) ** 2))
        distance = np.array(distance)
        total_distance = np.sum(distance)
        N_per_segment = np.int64(np.ceil(self.n_interp * distance / total_distance))
        for ii in range(len(E) - 1):
            if E[ii] > E[ii + 1]:
                E0, E1 = E[ii + 1], E[ii]
                N0, N1 = N[ii + 1], N[ii]
            else:
                E0, E1 = E[ii], E[ii + 1]
                N0, N1 = N[ii], N[ii + 1]
            X = np.linspace(E0, E1, N_per_segment[ii])
            Y = np.interp(X, (E0, E1), (N0, N1))
            if E[ii] > E[ii + 1]:
                X = np.flip(X, axis=0)
                Y = np.flip(Y, axis=0)
            qx.append(X)
            qy.append(Y)
        qx = np.concatenate(qx).ravel()
        qy = np.concatenate(qy).ravel()
        return qx, qy

    def generate_transect2D(self):
        """Build the 2-D transect section: densify the picked polyline,
        sample :math:`\\log_{10}\\rho` along it at every model depth with a
        :class:`~scipy.interpolate.RegularGridInterpolator`, store the
        result and hand it to :meth:`plot_transect2D`."""
        if not (self.transect_plot['easting'] and self.transect_plot['northing']):
            return
        qx, qy = self.interpolate_between_points(self.transect_plot['easting'],
                                                 self.transect_plot['northing'])
        qz = np.array(self.clip_model.dz)
        query_points = np.zeros((len(qx) * len(qz), 3))
        faces = np.zeros(((len(qx) - 1) * (len(qz) - 1), 5))
        cc = 0
        cc2 = 0
        for ix in range(len(qx)):
            for ii, iz in enumerate(qz):
                query_points[cc, :] = np.array((qx[ix], qy[ix], iz))
                if ii < len(qz) - 1 and ix < len(qx) - 1:
                    faces[cc2, 0] = 4
                    faces[cc2, 1] = cc
                    faces[cc2, 2] = cc + 1
                    faces[cc2, 3] = len(qz) + cc + 1
                    faces[cc2, 4] = len(qz) + cc
                    cc2 += 1
                cc += 1

        cell_N, cell_E, cell_z = self.clip_model.cell_centers()
        rho = getattr(self.clip_model, self.rho_axis)
        vals = np.transpose(np.log10(rho), [1, 0, 2])
        self.interpolator = RGI((cell_E, cell_N, cell_z),
                                vals, bounds_error=False, fill_value=None)
        interp_vals = self.interpolator(query_points)
        self.interp_vals = np.reshape(interp_vals, [len(qx), len(qz)])
        self.qx = qx
        self.qy = qy
        self.qz = qz
        self.query_points = query_points
        self.faces = faces
        self.plot_transect2D(redraw=True)

    def generate_transect3D(self, redraw=True):
        """Build the 3-D transect surface by slicing ``rect_grid`` along
        each polyline segment and clipping to the model depth range, then
        :meth:`plot_transect3D`."""
        query_points = deepcopy(self.query_points)
        query_points[:, 2] *= -1
        qp = np.vstack([self.transect_plot['easting'],
                        self.transect_plot['northing'],
                        np.zeros((len(self.transect_plot['easting'])))]).T
        for ii in range(len(qp) - 1):
            line = pv.Line(pointa=qp[ii, :], pointb=qp[ii+1, :], resolution=10)
            if ii == 0:
                self.slices['transect'] = self.rect_grid.slice_along_line(line).clip_box(invert=False,
                                                                                         bounds=line.bounds[:4] +
                                                                                         (-self.clip_model.dz[-1],
                                                                                           self.clip_model.dz[0]))
            else:
                self.slices['transect'].merge(self.rect_grid.slice_along_line(line).clip_box(invert=False,
                                                                                             bounds=line.bounds[:4] +
                                                                                           (-self.clip_model.dz[-1],
                                                                                             self.clip_model.dz[0])),
                                              merge_points=False,
                                              inplace=True)
        self.plot_transect3D(redraw=redraw)

    def plot_transect2D(self, redraw=True):
        """Draw the interpolated transect section (linear distance vs
        depth) on the transect figure, then overlay the pick line /
        markers on the plan panel."""
        self.fig_transect.axes[0].clear()
        linear_x = np.zeros((len(self.qx)) + 1)
        linear_x[1:-1] = np.sqrt((self.qx[1:] - self.qx[:-1]) ** 2 + (self.qy[1:] - self.qy[:-1]) ** 2)
        linear_x = np.cumsum(linear_x)
        linear_x[-1] += (linear_x[-1] - linear_x[-1])  # legacy no-op, kept for parity
        self.fig_transect.axes[0].clear()
        self.fig_transect.axes[0].pcolormesh(linear_x, np.array(list(self.qz) + [self.qz[-1] + 1]),
                                             self.interp_vals.T,
                                             vmin=self.cax[0], vmax=self.cax[1],
                                             cmap=self.cmap)
        self.fig_transect.axes[0].invert_yaxis()
        self.fig_transect.axes[0].set_xlabel('Linear Distance ({})'.format(self.spatial_units))
        self.fig_transect.axes[0].set_ylabel('Depth ({})'.format(self.spatial_units))
        self.canvas['transect'].draw()
        self.plot_transect_line()
        self.plot_transect_markers()

    def plot_transect3D(self, redraw=False):
        """Add or remove the 3-D transect-surface actor per the
        ``interpCheckbox`` state."""
        if self.interpCheckbox.checkState():
            self.actors['transect'] = self.vtk_widget.add_mesh(self.slices['transect'],
                                                               style='surface',
                                                               cmap=self.cmap,
                                                               clim=self.cax,
                                                               name='transect')
        else:
            self.vtk_widget.remove_actor(self.actors['transect'])
        if redraw:
            self.show_bounds()

    def plot_transect_line(self, redraw=False):
        """Draw / move the dashed transect line on the plan panel.
        :meth:`plot_transect_markers` does the same for the vertex
        markers."""
        if self.show_transect_line and self.transect_plot['easting']:
            if not self.transect_artists['line'] or redraw:
                self.transect_artists['line'], = self.map.window['axes'][0].plot(self.qx, self.qy, 'k--')
            else:
                self.transect_artists['line'].set_xdata(self.qx)
                self.transect_artists['line'].set_ydata(self.qy)
                # self.debugLabel.setText('Redrawing lines')
            self.canvas['2D'].draw()

    def plot_transect_markers(self, redraw=False):
        if self.show_transect_markers and self.transect_plot['easting']:
            if not self.transect_artists['markers'] or redraw:
                self.transect_artists['markers'], = self.map.window['axes'][0].plot(self.transect_plot['easting'],
                                                                                    self.transect_plot['northing'],
                                                                                    'yv')
            else:
                self.transect_artists['markers'].set_xdata(self.transect_plot['easting'])
                self.transect_artists['markers'].set_ydata(self.transect_plot['northing'])
            self.canvas['2D'].draw()

    def toggle_transect3D(self):
        """``interpCheckbox`` handler: show / hide the 3-D transect
        surface (unchecks itself if no transect has been generated yet)."""
        if self.slices['transect'] != []:
            if self.interpCheckbox.checkState():
                self.plot_transect3D(redraw=True)
            else:
                self.vtk_widget.remove_actor(self.actors['transect'])
                self.show_bounds()
        else:
            self.interpCheckbox.setCheckState(False)


def main():
    """Console entry point: classify the positional file arguments with
    :func:`pyMT.utils.sort_files`, require a model file, and show the
    :class:`ModelWindow`."""
    files = sys.argv[1:]
    for file in files:
        if not check_file(file):
            print('File {} not found.'.format(file))
            return
    files = sort_files(files=files)
    try:
        model = DS.Model(files['model'])
    except KeyError:
        print('Model must be specified')
        return
    app = QtWidgets.QApplication(sys.argv)
    viewer = ModelWindow(files)
    viewer.setWindowIcon(QtGui.QIcon(model_viewer_jpg))
    viewer.show()
    ret = app.exec_()
    sys.exit(ret)


if __name__ == '__main__':
    main()
