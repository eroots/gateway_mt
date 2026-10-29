"""Before/after verification baseline for ``gateway_mt.GatewayGUI.gateway_main``.

``gateway_main`` is a thin PyQt launcher: a table of named datasets, a
save / load of that table as a ``.pymt`` project file, and three buttons
that spawn the other GUIs.  The spawning is covered by the individual GUI
baselines and by a manual smoke-test checklist (``ModelWindow`` segfaults
head-less on this machine).  This harness freezes the module's *pure*
logic, which runs fine under ``QT_QPA_PLATFORM=offscreen``:

* ``NewProject`` construction (with the open-project dialog stubbed out),
* ``sort_files`` extension classification for a spread of inputs,
* the dataset-table round trip: ``init_dataset_table`` ->
  ``dataset_dict`` / ``active_datasets`` / ``dataset_index`` /
  ``empty_table`` / ``get_locs_from``,
* the text ``write_project`` emits for a known table, and the
  ``FileInputParser.read_pystart`` parse of that same file.

Usage
-----
    python tests/gateway_gui_baseline.py save      # snapshot
    python tests/gateway_gui_baseline.py compare   # re-run and diff

Outputs land in ``tests/baseline_gateway_gui/``.
"""

import os
import sys
import pickle
import zlib
import traceback
import warnings

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import numpy as np

from PyQt5 import QtWidgets

warnings.filterwarnings('ignore')

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(THIS_DIR, os.pardir))
SYN = os.path.join(THIS_DIR, 'data', 'synthetic')

DATA_FILE = os.path.join(SYN, 'synthetic_modem.dat')
RESP_FILE = os.path.join(SYN, 'synthetic_modem_resp.dat')
LIST_FILE = os.path.join(SYN, 'edis', 'synthetic.lst')
MODEL_FILE = os.path.join(SYN, 'synthetic.model')

OUT_DIR = os.path.join(THIS_DIR, 'baseline_gateway_gui')
SCRATCH = os.path.join(OUT_DIR, 'scratch')
PICKLE_PATH = os.path.join(OUT_DIR, 'gateway_gui_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'gateway_gui_baseline_summary.txt')

_APP = None  # live QApplication reference (see _run)


# --------------------------------------------------------------------------- #
def _norm(txt):
    """Make written paths machine-independent."""
    for raw in (REPO, REPO.replace('\\', '/')):
        txt = txt.replace(raw, '<REPO>')
    return txt.replace('\\', '/')


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
            head = ', '.join('{:.6g}'.format(x) for x in flat[:6])
            return 'shape={} dtype={} head=[{}]'.format(arr.shape, arr.dtype, head)
        return 'shape={} dtype={} {}'.format(arr.shape, arr.dtype, arr.ravel()[:4])
    except Exception:  # noqa: BLE001
        return repr(value)[:160]


def _diff(o, n):
    if isinstance(o, str) or isinstance(n, str):
        if o == n:
            return ''
        if isinstance(o, str) and isinstance(n, str):
            ol, nl = o.splitlines(), n.splitlines()
            for i, (a, b) in enumerate(zip(ol, nl)):
                if a != b:
                    return 'text differs at line {}: {!r} -> {!r}'.format(i + 1, a[:90], b[:90])
            return 'text differs (len {} -> {})'.format(len(ol), len(nl))
        return 'type {} -> {}'.format(type(o).__name__, type(n).__name__)
    try:
        o_arr, n_arr = np.asarray(o), np.asarray(n)
    except Exception:  # noqa: BLE001
        return '' if o == n else '{!r} -> {!r}'.format(o, n)
    if o_arr.shape != n_arr.shape:
        return 'shape {} -> {}'.format(o_arr.shape, n_arr.shape)
    if not np.array_equal(o_arr, n_arr):
        return '{!r} -> {!r}'.format(o_arr.ravel()[:6], n_arr.ravel()[:6])
    return ''


# --------------------------------------------------------------------------- #
def _new_project(gm):
    """Construct ``NewProject`` head-less with the open-project dialog stubbed."""
    orig = QtWidgets.QFileDialog.getSaveFileName
    QtWidgets.QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: ('', ''))
    try:
        return gm.NewProject()
    finally:
        QtWidgets.QFileDialog.getSaveFileName = orig


def _fill_table(w, files_dict):
    """Repopulate ``datasetTable`` from a name -> {type: path} dict."""
    w.datasetTable.cellChanged.disconnect()
    w.init_dataset_table(files_dict)
    w.datasetTable.cellChanged.connect(w.modify_table_cell)


# --------------------------------------------------------------------------- #
def _run(results):
    import gateway_mt.GatewayGUI.gateway_main as gm
    from gateway_mt.GUI_common.classes import FileInputParser

    os.makedirs(SCRATCH, exist_ok=True)
    for f in os.listdir(SCRATCH):
        os.remove(os.path.join(SCRATCH, f))

    # NB: keep a live reference -- a garbage-collected QApplication segfaults
    # the first QWidget built after it.
    global _APP
    _APP = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)

    try:
        w = _new_project(gm)
        results['construct'] = np.asarray('ok')
    except Exception as exc:  # noqa: BLE001
        results['construct'] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
        traceback.print_exc()
        return

    # -- construction state -------------------------------------------------- #
    header = [str(w.datasetTable.horizontalHeaderItem(c).text())
              for c in range(w.datasetTable.columnCount())]
    results['init/rows'] = np.asarray(w.datasetTable.rowCount())
    results['init/cols'] = np.asarray(w.datasetTable.columnCount())
    results['init/headers'] = np.asarray(header, dtype=object)
    results['init/empty_table'] = np.asarray(bool(w.empty_table))
    results['init/get_locs_from'] = np.asarray(str(w.get_locs_from()), dtype=object)
    results['init/file_to_index'] = np.asarray(repr(w.file_to_index), dtype=object)
    results['init/index_to_file'] = np.asarray(repr(w.index_to_file), dtype=object)

    # -- sort_files -------------------------------------------------------- #
    def _sf(paths):
        return _norm(repr(w.sort_files(paths)))

    _store(results, 'sort_files/real_mixed', lambda: _sf(
        [DATA_FILE, RESP_FILE, LIST_FILE, MODEL_FILE]))
    _store(results, 'sort_files/data_only', lambda: _sf([DATA_FILE]))
    _store(results, 'sort_files/resp_by_line', lambda: _sf([RESP_FILE]))
    _store(results, 'sort_files/ext_only', lambda: _sf(
        ['/x/a.rho', '/x/b.model', '/x/c.zani', '/x/d.lst', '/x/e.list',
         '/x/f.reso', '/x/g.resolution', '/x/h.resp']))
    _store(results, 'sort_files/unknown_ext', lambda: _sf(['/x/a.txt', '/x/b.foo']))

    # -- dataset-table round trip ----------------------------------------- #
    files_dict = {
        'ds_full': {'list': LIST_FILE, 'data': DATA_FILE,
                    'response': RESP_FILE, 'model': MODEL_FILE, 'resolution': ''},
        'ds_data_only': {'data': DATA_FILE},
        'ds_empty': {},
    }
    _fill_table(w, files_dict)
    results['table/rows_after_fill'] = np.asarray(w.datasetTable.rowCount())
    results['table/empty_table'] = np.asarray(bool(w.empty_table))
    _store(results, 'table/dataset_dict', lambda: _norm(repr(w.dataset_dict())))

    w.datasetTable.clearSelection()
    results['table/dataset_index_none'] = np.asarray(list(w.dataset_index()), dtype=object)

    w.datasetTable.selectRow(0)
    results['table/dataset_index_row0'] = np.asarray(list(w.dataset_index()), dtype=object)
    _store(results, 'table/active_datasets_row0', lambda: _norm(repr(w.active_datasets())))

    w.datasetTable.selectRow(1)
    _store(results, 'table/active_datasets_row1', lambda: _norm(repr(w.active_datasets())))

    # -- write_project / read_pystart round trip ------------------------- #
    project_path = os.path.join(SCRATCH, 'proj.pymt')
    all_datasets = {
        'first': {'list': 'a.lst', 'data': 'a.dat', 'resp': 'a_resp.dat',
                  'model': 'a.model', 'resolution': ''},
        'second': {'list': '', 'data': 'b.data', 'resp': '',
                   'model': '', 'resolution': 'b.reso'},
    }
    orig = QtWidgets.QFileDialog.getSaveFileName
    QtWidgets.QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (project_path, 'x'))
    try:
        _store(results, 'write_project/call', lambda: (w.write_project(all_datasets), 'ok')[1])
    finally:
        QtWidgets.QFileDialog.getSaveFileName = orig

    if os.path.exists(project_path):
        with open(project_path, 'r') as fh:
            results['write_project/text'] = _norm(fh.read())
        _store(results, 'read_pystart/roundtrip',
               lambda: _norm(repr(FileInputParser.read_pystart(project_path))))
    else:
        results['write_project/text'] = 'EXCEPTION: project file not written'

    w.close()


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
            f.write('{:44s} {}\n'.format(label, _describe(results[label])))
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
    print('\n=== gateway_gui baseline comparison ===')
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
        sys.exit('usage: python tests/gateway_gui_baseline.py [save|compare]')
