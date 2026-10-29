"""Before/after verification baseline for ``pyMT.ModelGUI.mesh_designer``.

``mesh_designer`` is a matplotlib (not VTK) PyQt GUI, so it *can* run
head-less (``QT_QPA_PLATFORM=offscreen`` + Agg).  This harness builds a
real ``model_viewer_2d`` from the sample ModEM model + data, drives a
fixed script of mesh / depth / colour / background handlers with the
modal dialogs monkey-patched, and freezes:

* the plan-view ``pcolormesh`` array / clim and the overlaid line data,
* the model cell-edge arrays (``dx`` / ``dy`` / ``dz``) and
  ``(nx, ny, nz)``, plus ``vals`` statistics,
* the text of a written ModEM model + covariance,
* one entry per driven step recording whether it raised.

Usage
-----
    python tests/mesh_designer_baseline.py save      # snapshot
    python tests/mesh_designer_baseline.py compare   # re-run and diff

Outputs land in ``tests/baseline_mesh_designer/``.
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

from PyQt5 import QtWidgets

warnings.filterwarnings('ignore')

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(THIS_DIR, os.pardir))
SYN = os.path.join(THIS_DIR, 'data', 'synthetic')
MODEL_FILE = os.path.join(SYN, 'synthetic.model')
DATA_FILE = os.path.join(SYN, 'synthetic_modem.dat')

OUT_DIR = os.path.join(THIS_DIR, 'baseline_mesh_designer')
SCRATCH = os.path.join(OUT_DIR, 'scratch')
PICKLE_PATH = os.path.join(OUT_DIR, 'mesh_designer_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'mesh_designer_baseline_summary.txt')

_DATE_RE = re.compile(r'\d{4}[-/]\d{2}[-/]\d{2}([ T]\d{2}:\d{2}:\d{2})?')


# --------------------------------------------------------------------------- #
def _store(results, label, func):
    try:
        results[label] = func()
    except Exception as exc:  # noqa: BLE001
        results[label] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
        traceback.print_exc()


def _describe(value):
    if isinstance(value, str):
        if len(value) > 200:
            return 'text len={} crc32={:08x}'.format(len(value), zlib.crc32(value.encode()))
        return repr(value)
    try:
        arr = np.asarray(value)
        if np.issubdtype(arr.dtype, np.number):
            flat = arr.ravel().astype(float)
            head = ', '.join('{:.6g}'.format(x) for x in flat[:4])
            return 'shape={} dtype={} head=[{}]'.format(arr.shape, arr.dtype, head)
        return 'shape={} dtype={} {}'.format(arr.shape, arr.dtype, arr.ravel()[:3])
    except Exception:  # noqa: BLE001
        return repr(value)[:120]


def _diff(o, n, rtol=1e-9, atol=1e-12):
    if isinstance(o, str) or isinstance(n, str):
        if o == n:
            return ''
        if isinstance(o, str) and isinstance(n, str):
            ol, nl = o.splitlines(), n.splitlines()
            for i, (a, b) in enumerate(zip(ol, nl)):
                if a != b:
                    return 'text differs at line {}: {!r} -> {!r}'.format(i + 1, a[:80], b[:80])
            return 'text differs (len {} -> {})'.format(len(ol), len(nl))
        return 'type {} -> {}'.format(type(o).__name__, type(n).__name__)
    try:
        o_arr, n_arr = np.asarray(o, dtype=float), np.asarray(n, dtype=float)
    except (TypeError, ValueError):
        o_arr, n_arr = np.asarray(o), np.asarray(n)
        return '' if o_arr.shape == n_arr.shape and np.array_equal(o_arr, n_arr) \
            else '{!r} -> {!r}'.format(o_arr.ravel()[:4], n_arr.ravel()[:4])
    if o_arr.shape != n_arr.shape:
        return 'shape {} -> {}'.format(o_arr.shape, n_arr.shape)
    if o_arr.size and not np.allclose(o_arr, n_arr, rtol=rtol, atol=atol, equal_nan=True):
        return 'max abs diff = {:.3e}'.format(np.nanmax(np.abs(o_arr - n_arr)))
    return ''


def _norm_text(txt):
    txt = _DATE_RE.sub('<DATE>', txt)
    return txt.replace(SCRATCH.replace('\\', '/'), '<SCRATCH>').replace(SCRATCH, '<SCRATCH>')


# --------------------------------------------------------------------------- #
class _Patches:
    def __init__(self, save_path_provider):
        self._save_path = save_path_provider
        self._orig = []

    def _set(self, obj, name, fn):
        self._orig.append((obj, name, getattr(obj, name)))
        setattr(obj, name, fn)

    def __enter__(self):
        QID, QMB, QFD = (QtWidgets.QInputDialog, QtWidgets.QMessageBox, QtWidgets.QFileDialog)
        self._set(QID, 'getDouble', staticmethod(lambda *a, **k: (float(a[3]) if len(a) > 3 else 1.0, True)))
        self._set(QID, 'getInt', staticmethod(lambda *a, **k: (int(a[3]) if len(a) > 3 else 1, True)))
        self._set(QID, 'getItem', staticmethod(lambda *a, **k: (a[3][0] if len(a) > 3 and a[3] else '', True)))
        self._set(QMB, 'question', staticmethod(lambda *a, **k: QMB.Yes))
        self._set(QMB, 'warning', staticmethod(lambda *a, **k: QMB.Ok))
        self._set(QFD, 'getSaveFileName', staticmethod(lambda *a, **k: (self._save_path(), 'x')))
        return self

    def __exit__(self, *exc):
        for obj, name, val in reversed(self._orig):
            setattr(obj, name, val)
        self._orig = []


# --------------------------------------------------------------------------- #
def _snapshot(results, prefix, gui):
    m = gui.model
    results['{}/model/nxnynz'.format(prefix)] = np.asarray([m.nx, m.ny, m.nz])
    results['{}/model/dx'.format(prefix)] = np.asarray(m.dx, dtype=float)
    results['{}/model/dy'.format(prefix)] = np.asarray(m.dy, dtype=float)
    results['{}/model/dz'.format(prefix)] = np.asarray(m.dz, dtype=float)
    v = np.asarray(m.vals, dtype=float)
    results['{}/model/vals_stats'.format(prefix)] = np.asarray(
        [v.min(), v.max(), v.mean(), v.std(), v.size])

    ax = gui.plan_view
    results['{}/plan/xlim'.format(prefix)] = np.asarray(ax.get_xlim(), dtype=float)
    results['{}/plan/ylim'.format(prefix)] = np.asarray(ax.get_ylim(), dtype=float)
    results['{}/plan/aspect'.format(prefix)] = np.asarray(str(ax.get_aspect()))
    qm = getattr(gui, 'mesh_plot', None)
    if qm is not None:
        arr = np.asarray(qm.get_array(), dtype=float)
        results['{}/plan/mesh_array'.format(prefix)] = arr
        results['{}/plan/mesh_clim'.format(prefix)] = np.asarray(qm.get_clim(), dtype=float)
    for i, ln in enumerate(ax.get_lines()):
        results['{}/plan/line{}'.format(prefix, i)] = np.asarray(ln.get_xydata(), dtype=float)
    results['{}/labels'.format(prefix)] = np.asarray(
        [str(gui.nxLabel.text()), str(gui.nyLabel.text()), str(gui.nzLabel.text())], dtype=object)


# --------------------------------------------------------------------------- #
def _run(results):
    import pyMT.data_structures as DS
    import pyMT.ModelGUI.mesh_designer as md

    os.makedirs(SCRATCH, exist_ok=True)
    for f in os.listdir(SCRATCH):
        os.remove(os.path.join(SCRATCH, f))

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    model = DS.Model(MODEL_FILE)
    dataset = DS.Dataset(datafile=DATA_FILE)

    holder = {'path': ''}
    with _Patches(lambda: holder['path']):
        try:
            gui = md.model_viewer_2d(model=model, dataset=dataset)
            results['construct'] = np.asarray('ok')
        except Exception as exc:  # noqa: BLE001
            results['construct'] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
            traceback.print_exc()
            return
        _snapshot(results, 'initial', gui)

        steps = []

        def step(name, fn):
            steps.append(name)
            if os.environ.get('MD_TRACE'):
                print('  step: {}'.format(name), file=sys.stderr, flush=True)
            _store(results, 'step/{}'.format(name), lambda: (fn(), 'ok')[1])

        step('z_slice_2', lambda: (gui.zSlice.setValue(2), gui.validate_z_slice(2)))
        step('z_slice_0', lambda: (gui.zSlice.setValue(0), gui.validate_z_slice(0)))
        step('colourmap_jet', lambda: (gui.actionJet.setChecked(True), gui.set_colourmap()))
        step('colourmap_bwr', lambda: (gui.actionBwr.setChecked(True), gui.set_colourmap()))
        step('mesh_colour_white', lambda: (gui.action_meshWhite.setChecked(True), gui.set_mesh_colour()))
        step('marker_size', gui.set_marker_size)   # getDouble stubbed to (default, True)
        step('bg_rho', lambda: (gui.bgRho.setValue(1000.0),
                                setattr(gui, 'background_resistivity', 1000.0),
                                gui.set_background_rho()))
        step('smooth', lambda: (gui.sigmaX.setValue(1), gui.sigmaY.setValue(1),
                                gui.sigmaZ.setValue(1), gui.smooth_model()))
        step('add_pad_all', lambda: ([w.setChecked(True) for w in
                                      (gui.padLeft, gui.padRight, gui.padTop, gui.padBottom)],
                                     gui.add_pads()))
        step('remove_pad_all', lambda: gui.remove_pads())
        step('insert_y_line', lambda: _click(gui, 1, gui.model.dy[3] + 10.0, gui.model.dx[3] + 10.0))
        step('insert_x_line', lambda: _click(gui, 3, gui.model.dy[3] + 10.0, gui.model.dx[3] + 10.0))
        step('delete_y_line', lambda: _click_ctrl(gui, 1, gui.model.dy[4], gui.model.dx[4]))
        step('min_depth_edit', lambda: (gui.minDepth.setText('100'), gui.min_depth()))
        step('max_depth_edit', lambda: (gui.maxDepth.setText('500000'), gui.max_depth()))
        step('lock_aspect', lambda: (gui.actionLock_Aspect_Ratio.setChecked(True),
                                     gui.redraw_workaround()))
        # NOTE: regenerate_mesh (Model.generate_mesh) and generate_depths_decades
        # (Model.generate_zmesh) are heavy data_structures operations -- one hangs
        # on this fixture, the other needs a populated decade list -- and are out
        # of scope for the mesh_designer cleanup, so they are not driven here.

        holder['path'] = os.path.join(SCRATCH, 'out.model')
        step('write_modem', gui.write_modem)

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
                zlib.crc32(txt.encode()), len(lines), '\n'.join(lines[:40]), '\n'.join(lines[-40:]))


def _click(gui, button, xdata, ydata):
    """Synthesise the mouse-release event mesh_designer.click expects."""
    class _Ev:
        pass
    ev = _Ev()
    ev.inaxes = gui.plan_view
    ev.xdata, ev.ydata, ev.button = xdata, ydata, button
    gui.click(ev)


def _click_ctrl(gui, button, xdata, ydata):
    """A click with the Control modifier held (mesh-line delete)."""
    from PyQt5 import QtCore
    orig = QtWidgets.QApplication.keyboardModifiers
    QtWidgets.QApplication.keyboardModifiers = staticmethod(lambda: QtCore.Qt.ControlModifier)
    try:
        _click(gui, button, xdata, ydata)
    finally:
        QtWidgets.QApplication.keyboardModifiers = orig


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
            f.write('{:56s} {}\n'.format(label, _describe(results[label])))
    n_exc = sum(1 for v in results.values() if isinstance(v, str) and v.startswith('EXCEPTION'))
    print('Saved {} entries to {}'.format(len(results), PICKLE_PATH))
    if n_exc:
        print('WARNING: {} entries recorded an exception.'.format(n_exc))


def compare():
    if not os.path.exists(PICKLE_PATH):
        sys.exit('No baseline at {} -- run "save" first.'.format(PICKLE_PATH))
    with open(PICKLE_PATH, 'rb') as f:
        old = pickle.load(f)
    new = collect()
    only_old = sorted(set(old) - set(new))
    only_new = sorted(set(new) - set(old))
    mism, exc_changes = [], []
    for label in sorted(set(old) & set(new)):
        o, n = old[label], new[label]
        o_exc = isinstance(o, str) and o.startswith('EXCEPTION')
        n_exc = isinstance(n, str) and n.startswith('EXCEPTION')
        if o_exc or n_exc:
            if o != n:
                exc_changes.append((label, o, n))
            continue
        msg = _diff(o, n)
        if msg:
            mism.append((label, msg))
    print('\n=== mesh_designer baseline comparison ===')
    print('entries: {} old / {} new'.format(len(old), len(new)))
    for k in only_old:
        print('  - missing now: {}'.format(k))
    for k in only_new:
        print('  + new: {}'.format(k))
    for label, o, n in exc_changes:
        print('  ! {}\n      old: {}\n      new: {}'.format(label, o, n))
    for label, msg in mism:
        print('  * {}: {}'.format(label, msg))
    if not (only_old or only_new or exc_changes or mism):
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
        sys.exit('usage: python tests/mesh_designer_baseline.py [save|compare]')
