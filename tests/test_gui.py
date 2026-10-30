"""The GUIs, run offscreen: they load the synthetic project and draw it.

These are smoke tests: they catch crashes and wiring mistakes, not visual
details. Model Viewer's window needs a real OpenGL display, so only its
model-to-VTK conversion is tested here.
"""
import os

import numpy as np
import pytest

from conftest import quiet, synthetic

pytestmark = pytest.mark.usefixtures('qapp')


@pytest.fixture
def no_dialogs(monkeypatch):
    """Answer file dialogs with 'cancel' so constructors that ask for a file don't block."""
    from PyQt5 import QtWidgets
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', staticmethod(lambda *a, **k: ('', '')))
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getOpenFileName', staticmethod(lambda *a, **k: ('', '')))
    monkeypatch.setattr(QtWidgets.QMessageBox, 'warning', staticmethod(lambda *a, **k: None))


@pytest.fixture
def canned_dialogs(monkeypatch, tmp_path):
    """Answer every modal dialog the GUI actions can raise; file prompts get ``answers['save']``."""
    from PyQt5 import QtWidgets
    from gateway_mt.GUI_common import classes as gc
    answers = {'save': str(tmp_path / 'out')}
    QMB = QtWidgets.QMessageBox
    patches = [
        (QtWidgets.QFileDialog, 'getSaveFileName', staticmethod(lambda *a, **k: (answers['save'], 'All Files (*)'))),
        (QtWidgets.QFileDialog, 'getOpenFileName', staticmethod(lambda *a, **k: ('', ''))),
        (QMB, 'question', staticmethod(lambda *a, **k: QMB.Yes)),
        (QMB, 'warning', staticmethod(lambda *a, **k: QMB.Ok)),
        (QMB, 'information', staticmethod(lambda *a, **k: QMB.Ok)),
        (QtWidgets.QInputDialog, 'getDouble', staticmethod(lambda *a, **k: (float(a[3]) if len(a) > 3 else 1.0, True))),
        (QtWidgets.QInputDialog, 'getText', staticmethod(lambda *a, **k: ('', True))),
        (gc.MyPopupDialog, 'get_file', staticmethod(lambda *a, **k: (answers['save'], True))),
        (gc.TwoInputDialog, 'get_inputs', staticmethod(lambda l1, l2, initial_1=None, initial_2=None, **k:
                                                       ((initial_1, initial_2), True))),
        (gc.FileDialog, 'write_file', staticmethod(lambda *a, **k: ('', False))),
        (QtWidgets.QDialog, 'exec_', lambda self, *a, **k: 1),
    ]
    for obj, name, value in patches:
        monkeypatch.setattr(obj, name, value)
    return answers


def _project_files():
    """The synthetic .pymt project, with its relative paths resolved."""
    from gateway_mt.GUI_common.classes import FileInputParser
    files = FileInputParser.read_pystart(synthetic('synthetic.pymt'))
    return {name: {k: synthetic(v) for k, v in paths.items()} for name, paths in files.items()}


# --------------------------------------------------------------------------- #
# Gateway (project window)
# --------------------------------------------------------------------------- #
@pytest.fixture
def gateway(no_dialogs):
    from gateway_mt.GatewayGUI import gateway_main
    window = gateway_main.NewProject()
    yield window
    window.close()


def test_gateway_starts_with_one_empty_dataset(gateway):
    assert gateway.datasetTable.rowCount() == 1
    assert gateway.datasetTable.item(0, 0).text() == 'ds1'


def test_gateway_sorts_files_into_columns(gateway):
    files = [synthetic('synthetic.model'), synthetic('edis', 'synthetic.lst'), synthetic('synthetic_modem.dat')]
    assert gateway.sort_files(files) == {'model': files[0], 'list': files[1], 'data': files[2]}


def test_gateway_project_save_and_reload(gateway, monkeypatch, tmp_path):
    from PyQt5 import QtWidgets
    from gateway_mt.GUI_common.classes import FileInputParser
    project = str(tmp_path / 'test.pymt')
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', staticmethod(lambda *a, **k: (project, '')))
    gateway.write_project({'run1': {'list': 'a.lst', 'data': 'a.dat', 'resp': 'a_resp.dat', 'model': 'a.model'}})
    assert FileInputParser.read_pystart(project) == {
        'run1': {'list': 'a.lst', 'data': 'a.dat', 'response': 'a_resp.dat', 'model': 'a.model'}}


def test_gateway_project_keeps_the_resolution_file(gateway, monkeypatch, tmp_path):
    from PyQt5 import QtWidgets
    from gateway_mt.GUI_common.classes import FileInputParser
    project = str(tmp_path / 'test.pymt')
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', staticmethod(lambda *a, **k: (project, '')))
    gateway.write_project({'run1': {'model': 'a.model', 'resolution': 'a.res'}})
    assert FileInputParser.read_pystart(project)['run1'].get('resolution') == 'a.res'


@pytest.mark.parametrize('answer, expected', [(('renamed', True), 'renamed'), (('', False), 'ds1_2')])
def test_gateway_duplicate_names_can_be_renamed_or_cancelled(gateway, monkeypatch, answer, expected):
    from PyQt5 import QtWidgets
    monkeypatch.setattr(QtWidgets.QInputDialog, 'getText', staticmethod(lambda *a, **k: answer))
    table = gateway.datasetTable
    table.setRowCount(2)
    table.setItem(1, 0, QtWidgets.QTableWidgetItem('ds1'))         # same name as row 0
    assert table.item(1, 0).text() == expected
    assert table.item(0, 0).text() == 'ds1'


def test_gateway_gives_mesh_designer_resolved_paths(gateway, monkeypatch):
    """Paths in a project are relative to it; Mesh Designer must get them resolved."""
    from PyQt5 import QtWidgets
    from gateway_mt.GatewayGUI import gateway_main
    opened = {}
    monkeypatch.setattr(gateway_main, 'MeshDesigner', lambda **k: opened.update(k) or QtWidgets.QMainWindow())
    gateway.project_path = synthetic()
    gateway.datasetTable.setItem(0, gateway.file_to_index['model'], QtWidgets.QTableWidgetItem('synthetic.model'))
    gateway.datasetTable.selectRow(0)
    gateway.launch_mesh_designer()
    assert os.path.exists(opened['model'])


def test_gateway_forgets_closed_model_viewers(gateway):
    viewer = object()
    gateway.model_windows.append(viewer)
    gateway.close_mv(viewer)
    assert viewer not in gateway.model_windows


# --------------------------------------------------------------------------- #
# Data Plot
# --------------------------------------------------------------------------- #
@pytest.fixture
def data_plot(no_dialogs):
    import gateway_mt.DataGUI.data_plot as dp
    with quiet():
        window = dp.DataMain(_project_files())
    yield window
    window.close()


@pytest.mark.parametrize('args', [['data_plot'], ['data_plot', '-b']])
def test_data_plot_browses_for_a_start_file(monkeypatch, args):
    """With no options, or -b, data_plot asks for the start file and loads it."""
    import gateway_mt.DataGUI.data_plot as dp
    files = _project_files()
    monkeypatch.setattr(dp.FileInputParser, 'get_files_dialog', staticmethod(lambda: (files, True)))
    assert dp.parse_commandline(args) == files


def test_data_plot_loads_the_project(data_plot):
    assert data_plot.dataset.data.NS == 30
    assert data_plot.site_names == ['SYN00{}'.format(i) for i in range(1, 7)]   # first page of sites
    assert len(data_plot.dpm.axes) == 6
    assert len(data_plot.dpm.sites['response']) == 6


def test_data_plot_pages_through_sites(data_plot):
    with quiet():
        data_plot.Forward()
    assert data_plot.site_names[0] == 'SYN007'
    with quiet():
        data_plot.Back()
    assert data_plot.site_names[0] == 'SYN001'


def test_data_plot_writes_the_data_it_shows(data_plot, canned_dialogs, tmp_path):
    """File > Write > ModEM writes a file that reads back as the loaded data."""
    import gateway_mt.data_structures as DS
    canned_dialogs['save'] = str(tmp_path / 'written')
    with quiet():
        data_plot.write_ModEM()
        written = [f for f in tmp_path.iterdir() if f.name.startswith('written')]
        assert written, 'no file written'
        back = DS.Data(datafile=str(written[0]))
    original = data_plot.dataset.data
    assert list(back.site_names) == list(original.site_names)
    for name in ('SYN001', 'SYN030'):
        for comp in original.components:
            used = original.sites[name].used_error[comp] != original.REMOVE_FLAG
            np.testing.assert_allclose(back.sites[name].data[comp][used],
                                       original.sites[name].data[comp][used], rtol=1e-6, atol=1e-12)


def test_data_plot_rotation_round_trip(data_plot, canned_dialogs):
    before = data_plot.dataset.data.sites['SYN001'].data['ZXYR'].copy()
    with quiet():
        data_plot.set_azimuth(25)
    assert data_plot.dataset.azimuth == 25
    with quiet():
        data_plot.set_azimuth(0)
    np.testing.assert_allclose(data_plot.dataset.data.sites['SYN001'].data['ZXYR'], before, atol=1e-12)


@pytest.fixture
def typed_limits(monkeypatch):
    """Make every TwoInputDialog accept with the text '2' and '10'."""
    from PyQt5 import QtWidgets
    from gateway_mt.GUI_common import classes as gc
    monkeypatch.setattr(QtWidgets.QDialog, 'exec_', lambda self, *a, **k: 1)
    monkeypatch.setattr(gc.TwoInputDialog, 'inputs', lambda self: ('2', '10'))


def test_data_plot_limits_compare_numbers_not_text(data_plot, typed_limits):
    """'10' sorts before '2' as text; the limits must still be accepted."""
    with quiet():
        data_plot.set_rho_limits()
        data_plot.show_map()
        clim = data_plot.map_view.set_cax(1, 5)
    assert data_plot.dpm.ax_lim_dict['rho'] == [2.0, 10.0]
    assert clim == (2.0, 10.0)


def test_1d_window_layer_table_and_csv(data_plot, monkeypatch, tmp_path):
    """Build a layer in the 1-D window, mistype a value, and save the synthetic model."""
    from PyQt5 import QtWidgets
    with quiet():
        data_plot.launch_1D_window()
    window = data_plot.modeling_window
    table = window.synthetic_window.layerTable
    table.setItem(1, 0, QtWidgets.QTableWidgetItem('5'))           # 5 km layer, default resistivities
    table.setItem(1, 1, QtWidgets.QTableWidgetItem('abc'))         # not a number: put back
    assert float(table.item(1, 1).text()) == 1000
    out = str(tmp_path / 'synthetic.csv')
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getSaveFileName', staticmethod(lambda *a, **k: (out, '')))
    window.write_synthetic_model()
    rows = [line.split(',') for line in open(out).read().splitlines()]
    assert all(len(row) == 3 for row in rows)                     # depth, rho_x, rho_y


def test_1d_window_needs_data(no_dialogs):
    import gateway_mt.data_structures as DS
    from gateway_mt.GUI_common.windows import Modeling1D
    with pytest.raises(ValueError):
        Modeling1D(dataset=DS.Dataset())


def test_data_plot_map_view(data_plot):
    with quiet():
        data_plot.show_map()
    map_view = data_plot.map_view.map
    assert len(map_view.site_names) == 30
    assert map_view.window['axes'][0].collections          # stations drawn


# --------------------------------------------------------------------------- #
# Mesh Designer
# --------------------------------------------------------------------------- #
@pytest.fixture
def mesh_designer_window(no_dialogs):
    import gateway_mt.data_structures as DS
    from gateway_mt.ModelGUI import mesh_designer
    with quiet():
        model = DS.Model(synthetic('synthetic.model'))
        data = DS.Dataset(datafile=synthetic('synthetic_modem.dat'))
        window = mesh_designer.model_viewer_2d(model=model, dataset=data)
    yield window
    window.close()


def test_mesh_designer_opens_a_model_with_data(mesh_designer_window):
    assert mesh_designer_window.model.nx == 32
    assert np.asarray(mesh_designer_window.mesh_plot.get_array()).size == 32 * 32


def test_mesh_designer_adds_and_removes_padding(mesh_designer_window, canned_dialogs):
    window = mesh_designer_window
    nx, ny = window.model.nx, window.model.ny
    for box in (window.padLeft, window.padRight, window.padTop, window.padBottom):
        box.setChecked(True)
    with quiet():
        window.add_pads()
    assert (window.model.nx, window.model.ny) == (nx + 2, ny + 2)
    with quiet():
        window.remove_pads()
    assert (window.model.nx, window.model.ny) == (nx, ny)


def test_mesh_designer_writes_a_readable_model(mesh_designer_window, canned_dialogs, tmp_path):
    import gateway_mt.data_structures as DS
    canned_dialogs['save'] = str(tmp_path / 'designed.model')
    with quiet():
        mesh_designer_window.write_modem()
        back = DS.Model(str(tmp_path / 'designed.model'))
    assert (back.nx, back.ny, back.nz) == (32, 32, 35)
    np.testing.assert_allclose(back.vals, mesh_designer_window.model.vals, rtol=1e-5)


# --------------------------------------------------------------------------- #
# Model Viewer (conversion only; the window needs OpenGL)
# --------------------------------------------------------------------------- #
def test_model_viewer_converts_a_model_to_a_vtk_grid():
    pytest.importorskip('pyvista')
    import gateway_mt.data_structures as DS
    from gateway_mt.ModelGUI.model_viewer import model_to_rectgrid
    with quiet():
        model = DS.Model(synthetic('synthetic.model'))
    grid = model_to_rectgrid(model)
    assert grid.n_cells == model.nx * model.ny * model.nz
    resistivity = grid.cell_data['Resistivity']               # log10(ohm-m)
    assert resistivity.min() == pytest.approx(np.log10(model.vals.min()))
    assert resistivity.max() == pytest.approx(np.log10(model.vals.max()))
