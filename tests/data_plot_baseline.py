"""Before/after verification harness for ``pyMT.DataGUI.data_plot``.

A GUI has no natural numeric output, so this harness runs it head-less
(``QT_QPA_PLATFORM=offscreen``): it builds a real ``DataMain`` from the
SWZ ModEM fixtures, drives a fixed script of menu / toggle handlers with
the modal dialogs monkey-patched to canned answers, and freezes

* the ``DataPlotManager`` / ``MapView`` configuration state and the
  plotted artist data after the script,
* the ``Dataset`` state (periods / sites / components / a few sites'
  errors / RMS),
* the text of every file the GUI writes (list + ModEM / WSINV3DMT /
  GoFEM data), timestamps blanked,
* one entry per driven step recording whether it raised.

Usage
-----
    python tests/data_plot_baseline.py save      # snapshot
    python tests/data_plot_baseline.py compare   # re-run and diff

Outputs land in ``tests/baseline_data_plot/``.
"""

import os
import re
import sys
import pickle
import zlib
import traceback
import warnings

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import numpy as np
import matplotlib
matplotlib.use('Agg')

from PyQt5 import QtWidgets, QtCore

warnings.filterwarnings('ignore')

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(THIS_DIR, os.pardir))
CTD = os.path.join(REPO, 'test_data', 'claude_test_data')
PYSTART = os.path.join(CTD, 'pytest.pymt')

OUT_DIR = os.path.join(THIS_DIR, 'baseline_data_plot')
SCRATCH = os.path.join(OUT_DIR, 'scratch')
PICKLE_PATH = os.path.join(OUT_DIR, 'data_plot_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'data_plot_baseline_summary.txt')

_DATE_RE = re.compile(r'\d{4}[-/]\d{2}[-/]\d{2}[ T]\d{2}:\d{2}:\d{2}(\.\d+)?')
_DATE_RE2 = re.compile(r'\d{4}-\d{2}-\d{2}')


# --------------------------------------------------------------------------- #
# snapshot machinery (shared shape with the other baseline harnesses)
# --------------------------------------------------------------------------- #
def _store(results, label, func, *args, **kwargs):
    try:
        value = func(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001
        value = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
        traceback.print_exc()
    results[label] = value
    return value


def _describe(value):
    if isinstance(value, str):
        if len(value) > 200:
            return 'text len={} crc32={:08x}'.format(len(value), zlib.crc32(value.encode()))
        return repr(value)
    try:
        arr = np.asarray(value)
        if arr.dtype == object:
            return 'object array shape={}'.format(arr.shape)
        flat = arr.ravel()
        if np.issubdtype(arr.dtype, np.number):
            head = ', '.join('{:.6g}'.format(x) for x in flat[:4].astype(float))
        else:
            head = str(flat[:4])
        return 'shape={} dtype={} head=[{}]'.format(arr.shape, arr.dtype, head)
    except Exception:  # noqa: BLE001
        return repr(value)[:120]


def _diff(o, n, rtol, atol):
    if isinstance(o, str) or isinstance(n, str):
        if o == n:
            return ''
        if isinstance(o, str) and isinstance(n, str):
            ol, nl = o.splitlines(), n.splitlines()
            for i, (a, b) in enumerate(zip(ol, nl)):
                if a != b:
                    return 'text differs at line {}: {!r} -> {!r}'.format(i + 1, a[:80], b[:80])
            return 'text differs (length {} -> {})'.format(len(ol), len(nl))
        return 'type changed {} -> {}'.format(type(o).__name__, type(n).__name__)
    try:
        o_arr = np.asarray(o, dtype=float)
        n_arr = np.asarray(n, dtype=float)
    except (TypeError, ValueError):
        o_seq, n_seq = np.asarray(o, dtype=object), np.asarray(n, dtype=object)
        if o_seq.shape != n_seq.shape:
            return 'object shape {} -> {}'.format(o_seq.shape, n_seq.shape)
        for i in range(o_seq.size):
            sub = _diff(o_seq.ravel()[i], n_seq.ravel()[i], rtol, atol)
            if sub:
                return 'element {}: {}'.format(i, sub)
        return ''
    if o_arr.shape != n_arr.shape:
        return 'shape {} -> {}'.format(o_arr.shape, n_arr.shape)
    if o_arr.size == 0:
        return ''
    if not np.allclose(o_arr, n_arr, rtol=rtol, atol=atol, equal_nan=True):
        return 'max abs diff = {:.3e}'.format(np.nanmax(np.abs(o_arr - n_arr)))
    return ''


def _tol_for(label):
    loose = ('/axes/', 'written/', 'locations', 'rms')
    return (1e-6, 1e-6) if any(k in label for k in loose) else (1e-9, 1e-12)


def _norm_text(txt):
    txt = _DATE_RE.sub('<DATE>', txt)
    txt = _DATE_RE2.sub('<DATE>', txt)
    txt = txt.replace(SCRATCH.replace('\\', '/'), '<SCRATCH>').replace(SCRATCH, '<SCRATCH>')
    return txt


# --------------------------------------------------------------------------- #
# monkey-patched dialogs
# --------------------------------------------------------------------------- #
class _Patches:
    """Install canned answers for every modal the driven handlers can hit."""

    def __init__(self, save_path_provider):
        self._save_path = save_path_provider
        self._orig = []

    def _set(self, obj, name, fn):
        self._orig.append((obj, name, getattr(obj, name)))
        setattr(obj, name, fn)

    def __enter__(self):
        from pyMT.GUI_common import classes as gc
        QID = QtWidgets.QInputDialog
        QMB = QtWidgets.QMessageBox
        QFD = QtWidgets.QFileDialog
        QD = QtWidgets.QDialog
        self._set(QID, 'getDouble', staticmethod(lambda *a, **k: (float(a[3]) if len(a) > 3 else 1.0, True)))
        self._set(QID, 'getInt', staticmethod(lambda *a, **k: (int(a[3]) if len(a) > 3 else 1, True)))
        self._set(QID, 'getItem', staticmethod(lambda *a, **k: (a[3][0] if len(a) > 3 and a[3] else '', True)))
        self._set(QID, 'getText', staticmethod(lambda *a, **k: ('', True)))
        self._set(QMB, 'question', staticmethod(lambda *a, **k: QMB.Yes))
        self._set(QMB, 'warning', staticmethod(lambda *a, **k: QMB.Ok))
        self._set(QMB, 'information', staticmethod(lambda *a, **k: QMB.Ok))
        self._set(QFD, 'getSaveFileName', staticmethod(lambda *a, **k: (self._save_path(), True)))
        self._set(QFD, 'getOpenFileName', staticmethod(lambda *a, **k: ('', False)))
        # custom QDialog subclasses in GUI_common -- all block on .exec_()
        self._set(gc.MyPopupDialog, 'get_file',
                  staticmethod(lambda *a, **k: (self._save_path(), True)))
        self._set(gc.TwoInputDialog, 'get_inputs',
                  staticmethod(lambda label_1, label_2, initial_1=None, initial_2=None, **k:
                               ((initial_1, initial_2), True)))
        self._set(gc.FileDialog, 'write_file',
                  staticmethod(lambda *a, **k: ('', False)))
        self._set(QD, 'exec_', lambda self_, *a, **k: 1)   # blanket fallback
        return self

    def __exit__(self, *exc):
        for obj, name, val in reversed(self._orig):
            setattr(obj, name, val)
        self._orig = []


# --------------------------------------------------------------------------- #
# state snapshot
# --------------------------------------------------------------------------- #
DPM_ATTRS = ('components', 'scale', 'errors', 'which_errors', 'wrap_phase',
             'show_outliers', 'plot_flagged_data', 'outlier_thresh',
             'min_ylim', 'max_ylim', 'link_axes_bounds', 'pt_units')
MAP_ATTRS = ('rho_cax', 'phase_cax', 'phase_split_cax', 'delta_cax', 'diff_cax',
             'tipper_cax', 'skew_cax', 'model_cax', 'aniso_cax', 'rms_cax',
             'colourmap', 'invert_cmap', 'lut', 'coordinate_system', 'annotate_sites',
             'markersize', 'site_marker', 'site_fill', 'plot_rms', 'rms_plot_style',
             'induction_scale', 'induction_cutoff', 'induction_error_tol',
             'pt_scale', 'min_pt_ratio', 'pt_difference_definition', 'spatial_units',
             'interpolant', 'include_outliers', 'allowed_std')
SAMPLE_SITES = ('e91003', 'SWZ001M', 'SWZ050M', 'SWZ131L')


def _snapshot(results, prefix, gui):
    dpm, mv, ds = gui.dpm, gui.map_view.map, gui.dataset

    for a in DPM_ATTRS:
        results['{}/dpm/{}'.format(prefix, a)] = np.asarray(getattr(dpm, a), dtype=object)
    results['{}/dpm/site_names'.format(prefix)] = np.asarray(dpm.site_names, dtype=object)
    results['{}/dpm/site_counts'.format(prefix)] = np.asarray(
        [len(dpm.sites[k]) for k in ('raw_data', 'data', 'response', 'smoothed_data', '1d')])
    results['{}/dpm/toggles'.format(prefix)] = np.asarray(
        [int(bool(dpm.toggles[k])) for k in sorted(dpm.toggles)])
    results['{}/dpm/n_axes'.format(prefix)] = np.asarray(len(dpm.axes))
    results['{}/dpm/tiling'.format(prefix)] = np.asarray(dpm.tiling)

    for i, ax in enumerate(dpm.axes):
        lines = [np.asarray(ln.get_xydata(), dtype=float) for ln in ax.get_lines()]
        if lines:
            results['{}/dpm/axes/{}/xy'.format(prefix, i)] = np.concatenate(lines) if len(lines) > 1 else lines[0]
        results['{}/dpm/axes/{}/ylim'.format(prefix, i)] = np.asarray(ax.get_ylim(), dtype=float)
        results['{}/dpm/axes/{}/title'.format(prefix, i)] = np.asarray(str(ax.get_title()))

    for a in MAP_ATTRS:
        results['{}/map/{}'.format(prefix, a)] = np.asarray(getattr(mv, a), dtype=object)
    for key in ('all', 'active', 'generic'):
        loc = mv.site_locations.get(key)
        if loc is not None and len(loc):
            arr = np.asarray(loc, dtype=float)
            # 'generic' row order follows MapView.generic_sites, a set difference;
            # sort rows so set-iteration order isn't a spurious diff.
            results['{}/map/locations/{}'.format(prefix, key)] = arr[np.lexsort(arr.T[::-1])]

    results['{}/ds/NP_NS_NR'.format(prefix)] = np.asarray(
        [ds.data.NP, ds.data.NS, ds.data.NR])
    results['{}/ds/inv_type'.format(prefix)] = np.asarray(ds.data.inv_type or 0)
    results['{}/ds/azimuth'.format(prefix)] = np.asarray(float(ds.azimuth))
    results['{}/ds/components'.format(prefix)] = np.asarray(sorted(ds.data.components), dtype=object)
    for name in SAMPLE_SITES:
        if name in ds.data.sites:
            s = ds.data.sites[name]
            results['{}/ds/{}/periods'.format(prefix, name)] = np.asarray(s.periods, dtype=float)
            for comp in sorted(s.components):
                results['{}/ds/{}/data/{}'.format(prefix, name, comp)] = np.asarray(s.data[comp], dtype=float)
                results['{}/ds/{}/used_error/{}'.format(prefix, name, comp)] = np.asarray(s.used_error[comp], dtype=float)
    if ds.rms:
        results['{}/ds/rms_total'.format(prefix)] = np.asarray(
            [ds.rms['Station'][s]['Total'] for s in sorted(ds.data.site_names)], dtype=float)


# --------------------------------------------------------------------------- #
# the driven script
# --------------------------------------------------------------------------- #
def _run(results):
    from pyMT.GUI_common.classes import FileInputParser
    import pyMT.DataGUI.data_plot as dp

    os.makedirs(SCRATCH, exist_ok=True)
    for f in os.listdir(SCRATCH):
        os.remove(os.path.join(SCRATCH, f))

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    files = FileInputParser.read_pystart(PYSTART)

    holder = {'path': ''}

    with _Patches(lambda: holder['path']):
        try:
            gui = dp.DataMain(files)
            results['construct'] = np.asarray('ok')
        except Exception as exc:  # noqa: BLE001
            results['construct'] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
            traceback.print_exc()
            return
        _snapshot(results, 'initial', gui)

        steps = []

        def step(name, fn):
            steps.append(name)
            _store(results, 'step/{}'.format(name), lambda: (fn(), 'ok')[1])

        mm = gui.map_view

        # --- data-panel toggles / options ---
        step('toggle_response_off', lambda: (gui.toggleResponse.setChecked(False), gui.toggle_response(0)))
        step('toggle_response_on', lambda: (gui.toggleResponse.setChecked(True), gui.toggle_response(2)))
        step('scaling_log', lambda: (gui.actionSqrtScaling.setChecked(False),
                                     gui.actionLogScaling.setChecked(True), gui.change_scaling(0)))
        step('scaling_sqrt', lambda: (gui.actionLogScaling.setChecked(False),
                                      gui.actionSqrtScaling.setChecked(True), gui.change_scaling(0)))
        step('components_rho', lambda: (setattr(gui.dpm, 'components', ['RhoXY', 'RhoYX']),
                                        gui.update_dpm(updated_sites=gui.dpm.site_names,
                                                       updated_comps=['RhoXY', 'RhoYX'])))
        step('components_phase', lambda: (setattr(gui.dpm, 'components', ['PhaXY']),
                                          gui.update_dpm(updated_sites=gui.dpm.site_names,
                                                         updated_comps=['PhaXY'])))
        step('phase_wrap_off', lambda: (gui.actionPhase_Wrap.setChecked(False), gui.set_phase_wrap()))
        step('components_z', lambda: (setattr(gui.dpm, 'components', ['ZXYR']),
                                      gui.update_dpm(updated_sites=gui.dpm.site_names,
                                                     updated_comps=['ZXYR'])))
        step('forward', gui.Forward)
        step('back', gui.Back)
        step('sort_south_north', lambda: (gui.sortSites.setCurrentIndex(2), gui.sort_sites(2)))
        step('sort_default', lambda: (gui.sortSites.setCurrentIndex(0), gui.sort_sites(0)))
        step('toggle_outliers_off', lambda: gui.toggle_outliers(0))
        step('toggle_outliers_on', lambda: gui.toggle_outliers(2))
        step('flag_a_period', lambda: _flag_one_period(gui))
        step('unflag_a_period', lambda: (gui.dataset.data.reset_errors(), gui.update_dpm()))
        step('azimuth_10', lambda: gui.set_azimuth(10))
        step('azimuth_0', lambda: gui.set_azimuth(0))

        # --- map window ---
        step('map_pt_data', lambda: (mm.toggle_dataPhaseTensor.setChecked(True), mm.data_phase_tensor()))
        step('map_pt_none', lambda: (mm.toggle_nonePhaseTensor.setChecked(True), mm.none_phase_tensor()))
        step('map_induction_data', lambda: (mm.toggle_dataInduction.setChecked(True), mm.update_map()))
        step('map_pseudo_data', lambda: (mm.toggle_dataInduction.setChecked(False),
                                         mm.toggle_dataPseudo.setChecked(True), mm.update_map()))
        step('map_pseudo_off', lambda: (mm.toggle_dataPseudo.setChecked(False), mm.update_map()))
        step('map_coord_utm', lambda: (mm.action_coordUTM.setChecked(True), mm.coord_system()))
        step('map_coord_local', lambda: (mm.action_coordLocal.setChecked(True), mm.coord_system()))
        step('map_planview', lambda: (mm.toggle_planView.setChecked(True), mm.update_map()))
        step('map_planview_off', lambda: (mm.toggle_planView.setChecked(False), mm.update_map()))

        # --- file writes (the real regression signal) ---
        def save_path_for(tag, ext):
            return os.path.join(SCRATCH, tag + ext)

        holder['path'] = save_path_for('out_list', '.lst')
        step('write_list', gui.write_list)
        holder['path'] = save_path_for('out_modem', '.dat')
        step('write_modem', gui.write_ModEM)
        holder['path'] = save_path_for('out_wsinv', '.data')
        step('write_wsinv', gui.write_wsinv3dmt)
        holder['path'] = save_path_for('out_gofem', '.gdat')
        step('write_gofem', gui.write_gofem)

        results['steps_run'] = np.asarray(steps, dtype=object)
        _snapshot(results, 'final', gui)
        gui.close()

    for fn in sorted(os.listdir(SCRATCH)):
        p = os.path.join(SCRATCH, fn)
        try:
            with open(p, 'r', errors='replace') as fh:
                txt = _norm_text(fh.read())
        except Exception as exc:  # noqa: BLE001
            results['written/{}'.format(fn)] = 'EXCEPTION: {}'.format(exc)
            continue
        if len(txt) <= 100_000:
            results['written/{}'.format(fn)] = txt
        else:
            lines = txt.splitlines()
            results['written/{}'.format(fn)] = 'crc32={:08x} nlines={}\nHEAD:\n{}\nTAIL:\n{}'.format(
                zlib.crc32(txt.encode()), len(lines),
                '\n'.join(lines[:40]), '\n'.join(lines[-40:]))


def _flag_one_period(gui):
    """Flag period index 5 on the first plotted site, mimicking onpick_flagged."""
    ds = gui.dataset.data
    site = gui.dpm.site_names[0]
    for comp in ds.IMPEDANCE_COMPONENTS:
        if comp in ds.sites[site].used_error:
            ds.sites[site].used_error[comp][5] = ds.REMOVE_FLAG
    gui.update_dpm()


# --------------------------------------------------------------------------- #
def collect():
    results = {}
    _run(results)
    return results


def save():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = collect()
    with open(PICKLE_PATH, 'wb') as f:
        pickle.dump(results, f)
    with open(SUMMARY_PATH, 'w', encoding='utf-8') as f:
        for label in sorted(results):
            f.write('{:72s} {}\n'.format(label, _describe(results[label])))
    n_exc = sum(1 for v in results.values() if isinstance(v, str) and v.startswith('EXCEPTION'))
    print('\nSaved {} entries to {}'.format(len(results), PICKLE_PATH))
    if n_exc:
        print('WARNING: {} entries recorded an exception (see summary).'.format(n_exc))


def compare():
    if not os.path.exists(PICKLE_PATH):
        sys.exit('No baseline at {} -- run "save" first.'.format(PICKLE_PATH))
    with open(PICKLE_PATH, 'rb') as f:
        old = pickle.load(f)
    new = collect()
    only_old = sorted(set(old) - set(new))
    only_new = sorted(set(new) - set(old))
    mismatches, exc_changes = [], []
    for label in sorted(set(old) & set(new)):
        o, n = old[label], new[label]
        o_exc = isinstance(o, str) and o.startswith('EXCEPTION')
        n_exc = isinstance(n, str) and n.startswith('EXCEPTION')
        if o_exc or n_exc:
            if o != n:
                exc_changes.append((label, o, n))
            continue
        msg = _diff(o, n, *_tol_for(label))
        if msg:
            mismatches.append((label, msg))
    print('\n=== data_plot baseline comparison ===')
    print('entries: {} old / {} new'.format(len(old), len(new)))
    for k in only_old:
        print('  - missing now: {}'.format(k))
    for k in only_new:
        print('  + new: {}'.format(k))
    for label, o, n in exc_changes:
        print('  ! {}\n      old: {}\n      new: {}'.format(label, o, n))
    for label, msg in mismatches:
        print('  * {}: {}'.format(label, msg))
    if not (only_old or only_new or exc_changes or mismatches):
        print('\nALL MATCH -- no differences within tolerance.')
    else:
        sys.exit(1)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'save'
    if mode == 'save':
        save()
    elif mode == 'compare':
        compare()
    else:
        sys.exit('usage: python tests/data_plot_baseline.py [save|compare]')
