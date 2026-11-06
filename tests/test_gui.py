"""The GUIs, run offscreen: they load the synthetic project and draw it.

These are smoke tests: they catch crashes and wiring mistakes, not visual
details. Model Viewer's window needs a real OpenGL display, so it isn't
built here; its helpers are tested in test_model_tools.
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


def test_launcher_starts_what_is_installed(monkeypatch):
    """Each program runs as its own process; one that isn't installed is greyed out."""
    from gateway_mt import launcher
    started = []
    monkeypatch.setattr(launcher, 'find_command', {'gateway_mt': '/bin/gateway_mt', 'mtplot': '/bin/mtplot'}.get)
    monkeypatch.setattr(launcher.subprocess, 'Popen', lambda command: started.append(command))
    window = launcher.Launcher()
    assert list(window.buttons) == ['Gateway', 'mtinv', 'mtplot']
    assert not window.buttons['mtinv'].isEnabled() and 'installed' in window.buttons['mtinv'].toolTip()
    window.buttons['mtinv'].click()
    assert started == []
    window.buttons['mtplot'].click()
    window.buttons['Gateway'].click()
    assert started == [['/bin/mtplot'], ['/bin/gateway_mt']]


def test_launcher_warns_if_a_program_will_not_start(monkeypatch):
    from PyQt5 import QtWidgets
    from gateway_mt import launcher

    def fail(command):
        raise OSError('broken')
    warnings = []
    monkeypatch.setattr(launcher, 'find_command', lambda command: '/bin/' + command)
    monkeypatch.setattr(launcher.subprocess, 'Popen', fail)
    monkeypatch.setattr(QtWidgets.QMessageBox, 'warning', staticmethod(lambda *a: warnings.append(a[2])))
    launcher.Launcher().buttons['mtinv'].click()
    assert warnings == ['Could not start mtinv: broken']


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


def test_start_file_with_raw_data_needs_a_list_file(tmp_path):
    from gateway_mt.GUI_common.classes import FileInputParser
    start = tmp_path / 'pystart'
    start.write_text('% with_list\nlist a.lst\nraw edis\n% without_list\ndata a.dat\nraw edis\n')
    with quiet():
        assert FileInputParser.read_pystart(str(start)) is False
    start.write_text('% with_list\nlist a.lst\nraw edis\n')
    assert FileInputParser.read_pystart(str(start)) == {'with_list': {'list': 'a.lst', 'raw_path': 'edis'}}


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
    with open(out) as f:
        rows = [line.split(',') for line in f.read().splitlines()]
    assert all(len(row) == 3 for row in rows)                     # depth, rho_x, rho_y
    # The response is the 1-D response of the table (5 km of 100 ohm-m over the 1000 ohm-m half-space)
    from helpers import one_d_impedance
    table.setItem(1, 1, QtWidgets.QTableWidgetItem('100'))
    table.setItem(1, 2, QtWidgets.QTableWidgetItem('100'))
    site = window.synthetic_window.site
    Z = one_d_impedance(site.periods, [100, 1000], [5000])
    np.testing.assert_allclose(site.data['ZXYR'], Z.real, rtol=1e-10)
    np.testing.assert_allclose(site.data['ZXYI'], -Z.imag, rtol=1e-10)


def test_1d_inversion_on_the_3d_model_layers(data_plot):
    """The 3-D model's depths become layer thicknesses, and the result is drawn against depth."""
    with quiet():
        data_plot.launch_1D_window()
        settings = data_plot.modeling_window.inversion_settings
        settings.numIteration.setValue(2)
        settings.check3DLayers.setChecked(True)
        settings.start_inversion()
    model = settings.results['Model']
    model_3d = data_plot.dataset.model
    edges = np.array(model_3d.dz) * (1000 if model_3d.spatial_units == 'km' else 1)   # in m
    np.testing.assert_allclose(model['dz'][:-1], np.diff(edges))          # + the appended half-space
    np.testing.assert_allclose(model['depth'][:-1], edges[:-1])
    inverted = [l for l in data_plot.modeling_window.synthetic_window.axis.get_lines() if l.get_label() == 'Inv.'][0]
    np.testing.assert_allclose(inverted.get_ydata()[:-1], model['depth'] / 1000)


def test_1d_window_needs_data(no_dialogs):
    import gateway_mt.data_structures as DS
    from gateway_mt.GUI_common.windows import Modeling1D
    with pytest.raises(ValueError):
        Modeling1D(dataset=DS.Dataset())


def test_data_plot_period_tolerances(data_plot):
    raw = data_plot.dataset.raw_data
    assert (data_plot.periodToleranceFlag.value(), data_plot.periodToleranceRemove.value()) == (10, 20)
    data_plot.periodToleranceFlag.setValue(5)
    data_plot.periodToleranceRemove.setValue(15)
    assert (raw.tol, raw.remove_tol) == pytest.approx((0.05, 0.15))


def test_data_plot_map_view(data_plot):
    with quiet():
        data_plot.show_map()
    map_view = data_plot.map_view.map
    assert len(map_view.site_names) == 30
    assert map_view.window['axes'][0].collections          # stations drawn


def test_map_viewer_reads_its_controls(data_plot):
    with quiet():
        data_plot.show_map()
        view = data_plot.map_view
        view.toggle_dataPhaseTensor.click()
    view.toggle_dataPseudo.setChecked(True)
    view.toggle_dataInduction.setChecked(True)
    view.arrowType.setCurrentIndex(view.arrowType.findText('Both'))
    settings = view.map_settings()
    assert settings.phase_tensors == ['data'] and settings.induction_arrows == ['data']
    assert settings.pseudosection == ['data'] and settings.pseudosection_fill == 'RhoXY'
    assert settings.arrow_type == ['R', 'I'] and settings.plan_slice is None
    view.toggle_planView.setChecked(True)
    settings = view.map_settings()
    assert settings.plan_slice == 0 and settings.pseudosection == []
    assert not view.toggle_dataPseudo.isChecked()


def _select(site_list, names):
    for row in range(site_list.count()):
        site_list.item(row).setSelected(site_list.item(row).text() in names)


@pytest.mark.parametrize('pages, removed, expected', [
    (0, ['SYN003'], ['SYN001', 'SYN002', 'SYN007', 'SYN004', 'SYN005', 'SYN006']),
    (5, ['SYN030'], ['SYN025', 'SYN026', 'SYN027', 'SYN028', 'SYN029', 'SYN024']),   # at the end: fill from before
])
def test_data_plot_removes_and_restores_sites(data_plot, pages, removed, expected):
    """Removed sites on screen are replaced by the next ones; re-adding them restores the data."""
    original = data_plot.dataset.data.sites[removed[0]]
    with quiet():
        for _ in range(pages):
            data_plot.Forward()
        _select(data_plot.siteList, removed)
        data_plot.remove_sites()
    assert data_plot.dpm.site_names == expected
    assert data_plot.dataset.data.NS == 29 and data_plot.dataset.raw_data.NS == 29
    assert [data_plot.rmSitesList.item(0).text()] == removed
    with quiet():
        _select(data_plot.rmSitesList, removed)
        data_plot.add_sites()
    assert data_plot.dataset.data.NS == 30 and data_plot.dataset.raw_data.NS == 30
    assert data_plot.dataset.data.sites[removed[0]] is original
    assert data_plot.rmSitesList.count() == 0


def test_data_plot_drops_subplots_when_sites_run_out(data_plot):
    """Removing all but three sites leaves three subplots, none repeated."""
    removed = ['SYN{:03d}'.format(ii) for ii in range(1, 28)]
    with quiet():
        _select(data_plot.siteList, removed)
        data_plot.remove_sites()
    assert data_plot.dpm.site_names == ['SYN028', 'SYN029', 'SYN030']
    assert data_plot.site_names == data_plot.dpm.site_names
    assert len(data_plot.dpm.axes) == 3 and data_plot.numSubplots.value() == 3
    with quiet():
        _select(data_plot.rmSitesList, removed)
        data_plot.add_sites()
        data_plot.numSubplots.setValue(6)
    assert data_plot.dpm.site_names == ['SYN028', 'SYN029', 'SYN030', 'SYN001', 'SYN002', 'SYN003']


class _Pick:
    """Just enough of a matplotlib pick event for DataMain.onpick."""
    def __init__(self, ax, ind, button, raw):
        self.mouseevent = type('MouseEvent', (), {'inaxes': ax, 'button': button})()
        self.ind = [ind]
        properties = {'markeredgewidth': 0 if raw else 1, 'linestyle': 'None'}
        self.artist = type('Artist', (), {'properties': lambda self: properties})()


def test_data_plot_adds_and_removes_periods_by_clicking(data_plot):
    """Left-click a raw-data point to add its period; right-click a data point to remove it."""
    data, raw_site = data_plot.dataset.data, data_plot.dataset.raw_data.sites['SYN001']
    ax = data_plot.dpm.axes[0]
    ind = next(i for i, p in enumerate(raw_site.periods)
               if np.min(np.abs(data.periods - p) / p) > 0.2)              # a raw period not in the data
    period = raw_site.periods[ind]
    with quiet():
        data_plot.onpick(_Pick(ax, ind, button=1, raw=True))
    assert data.NP == 19 and np.any(np.isclose(data.periods, period))
    assert all(np.any(np.isclose(site.periods, period)) for site in data.sites.values())
    with quiet():
        data_plot.onpick(_Pick(ax, ind, button=1, raw=True))                 # already there: no change
    assert data.NP == 19
    with quiet():
        data_plot.onpick(_Pick(ax, list(data.sites['SYN001'].periods).index(period), button=3, raw=False))
    assert data.NP == 18 and not np.any(np.isclose(data.periods, period))


@pytest.mark.parametrize('keys, text, n_sites, n_periods, n_comps', [
    ('', '3', 1, 1, 1),
    ('Shift', '3', 1, 18, 1),
    ('Control', '3', 1, 1, 12),
    ('Alt', '3', 30, 1, 1),
    ('', '-2', 1, 1, 1),
])
def test_data_plot_error_tree_edits(data_plot, keys, text, n_sites, n_periods, n_comps):
    """An edit sets the error multiplier at one point, widened to all periods / components / sites by
    Shift / Ctrl / Alt. A negative value divides the error instead and resets the multiplier."""
    from PyQt5 import QtCore
    data = data_plot.dataset.data
    before = {name: {comp: (site.errmap[comp].copy(), site.used_error[comp].copy())
                     for comp in data.components} for name, site in data.sites.items()}
    item = data_plot.tree_dict['SYN002'].child(4)                       # SYN002, 5th period
    column = 2 + 2                                                     # 3rd component
    modifiers = {'': QtCore.Qt.NoModifier, 'Shift': QtCore.Qt.ShiftModifier,
                 'Control': QtCore.Qt.ControlModifier, 'Alt': QtCore.Qt.AltModifier}
    data_plot.stored_key_presses = modifiers[keys]
    data_plot.current_tree_item, data_plot.current_tree_col = item, column
    data_plot.old_val = item.text(column)
    with quiet():
        item.setText(column, text)                                     # emits itemChanged -> post_edit_error
    changed = [(name, comp, ii) for name, site in data.sites.items() for comp in data.components
               for ii in np.flatnonzero(site.used_error[comp] != before[name][comp][1])]
    assert len({c[0] for c in changed}) == n_sites
    assert len({c[2] for c in changed}) == n_periods
    assert len({c[1] for c in changed}) == n_comps
    name, comp, ii = changed[0]
    old_map, old_error = before[name][comp][0][ii], before[name][comp][1][ii]
    if text == '3':
        assert data.sites[name].errmap[comp][ii] == 3
        assert data.sites[name].used_error[comp][ii] == pytest.approx(old_error / old_map * 3)
    else:
        assert data.sites[name].errmap[comp][ii] == 1
        assert data.sites[name].used_error[comp][ii] == pytest.approx(old_error / 2)


def test_data_plot_error_tree_rejects_zero(data_plot):
    data = data_plot.dataset.data
    before = data.sites['SYN002'].used_error[data.components[0]].copy()
    item = data_plot.tree_dict['SYN002'].child(0)
    data_plot.current_tree_item, data_plot.current_tree_col = item, 2
    data_plot.old_val = item.text(2)
    with quiet():
        item.setText(2, '0')
    assert item.text(2) == data_plot.old_val
    np.testing.assert_array_equal(data.sites['SYN002'].used_error[data.components[0]], before)


def test_replacement_sites():
    from gateway_mt.DataGUI.data_plot import replacement_sites
    names = ['S{}'.format(i) for i in range(10)]
    assert replacement_sites(names, ['S2', 'S3', 'gone'], 1) == ['S4']           # next after the plotted ones
    assert replacement_sites(names, ['S7', 'S8', 'S9', 'gone'], 2) == ['S6', 'S5']  # at the end: before them
    assert replacement_sites(names, ['S6', 'S7', 'S8'], 2) == ['S9', 'S5']         # after, then before
    assert replacement_sites(['S0', 'S1'], ['S0', 'gone'], 3) == ['S1']                # not enough: fewer


def test_error_tree_values():
    from gateway_mt.DataGUI.data_plot import error_multiplier, period_label, period_from_label
    assert error_multiplier('4') == (4, 4)
    assert error_multiplier('-2') == (1, -2)
    for bad in ('0', '1.5', 'x'):
        with pytest.raises(ValueError):
            error_multiplier(bad)
    assert period_label(0.01) == -100 and period_label(10.0) == 10
    assert period_from_label(str(period_label(0.01))) == pytest.approx(0.01)
    assert period_from_label('10.0') == 10


def test_data_plot_keeps_removed_sites_with_their_dataset(no_dialogs):
    """Each dataset has its own removed sites; switching shows that dataset's list."""
    import gateway_mt.DataGUI.data_plot as dp
    files = _project_files()
    (name, paths), = files.items()
    with quiet():
        window = dp.DataMain({'a': dict(paths), 'b': dict(paths)})
    try:
        a, b = window.stored_datasets['a'], window.stored_datasets['b']
        with quiet():
            _select(window.siteList, ['SYN010'])
            window.remove_sites()
            window.currentDataset.setCurrentIndex(1)
        assert window.dataset is b and window.rmSitesList.count() == 0
        assert b.data.NS == 30 and b.data.site_names.count('SYN010') == 1
        with quiet():
            window.currentDataset.setCurrentIndex(0)
        assert [window.rmSitesList.item(i).text() for i in range(window.rmSitesList.count())] == ['SYN010']
        with quiet():
            _select(window.rmSitesList, ['SYN010'])
            window.add_sites()
        assert a.data.NS == a.raw_data.NS == 30 and not a.removed_sites
    finally:
        window.close()


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


def test_mesh_designer_regenerate_mesh_sizes(mesh_designer_window, monkeypatch):
    """Cell sizes start from the station spread; a huge mesh is only built if confirmed."""
    import gateway_mt.utils as utils
    from PyQt5 import QtWidgets
    window = mesh_designer_window
    locs = window.site_locations
    assert float(window.minX.text()) == pytest.approx(utils.default_cell_size(locs[:, 0]), abs=0.5)
    asked = []
    monkeypatch.setattr(QtWidgets.QMessageBox, 'question',
                        staticmethod(lambda *a, **k: asked.append(a[2]) or QtWidgets.QMessageBox.No))
    nx = window.model.nx
    window.minX.setText('10')
    window.minY.setText('10')
    with quiet():
        window.regenerate_mesh()
    assert asked and window.model.nx == nx                         # asked, declined, nothing changed
    window.minX.setText('2000')
    window.minY.setText('2000')
    with quiet():
        window.regenerate_mesh()
    assert len(asked) == 1
    assert window.model.nx == utils.regular_mesh_size(locs[:, 0], 2000) - 1


def test_mesh_designer_writes_a_readable_model(mesh_designer_window, canned_dialogs, tmp_path):
    import gateway_mt.data_structures as DS
    canned_dialogs['save'] = str(tmp_path / 'designed.model')
    with quiet():
        mesh_designer_window.write_modem()
        back = DS.Model(str(tmp_path / 'designed.model'))
    assert (back.nx, back.ny, back.nz) == (32, 32, 35)
    np.testing.assert_allclose(back.vals, mesh_designer_window.model.vals, rtol=1e-5)


# --------------------------------------------------------------------------- #
# Model Viewer (the window needs OpenGL; its helpers are tested in test_model_tools)
# --------------------------------------------------------------------------- #
def test_model_viewer_helpers_still_import_from_it():
    from gateway_mt import model_tools
    from gateway_mt.ModelGUI import model_viewer
    for name in ('model_to_rectgrid', 'sample_path', 'nearest_point', 'project_onto_path',
                 'export_to_paraview', 'write_section_csv'):
        assert getattr(model_viewer, name) is getattr(model_tools, name)
