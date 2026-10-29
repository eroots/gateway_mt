"""Before/after verification baseline for ``pyMT.ModelGUI.model_viewer``.

The GUI itself (``ModelWindow``) needs a live VTK / OpenGL context and
segfaults when run head-less on this machine, so -- as CLAUDE.md
prescribes for GUIs -- it is covered by a manual smoke-test checklist,
not an automated harness.

What *is* automatable is the one reusable, render-free entry point:
``model_to_rectgrid``.  This harness freezes its output
(``pyvista.RectilinearGrid`` cell-data arrays, dimensions, bounds) for
the sample ModEM model, with and without a resolution volume, across
every ``rho_axis`` option.

Usage
-----
    python tests/model_viewer_baseline.py save      # snapshot
    python tests/model_viewer_baseline.py compare   # re-run and diff

Outputs land in ``tests/baseline_model_viewer/``.
"""

import os
import sys
import pickle
import traceback
import warnings
from copy import deepcopy

import numpy as np

import pyMT.data_structures as DS
from pyMT.ModelGUI.model_viewer import model_to_rectgrid

warnings.filterwarnings('ignore')

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(THIS_DIR, os.pardir))
MODEL_FILE = os.path.join(THIS_DIR, 'data', 'synthetic', 'synthetic.model')

OUT_DIR = os.path.join(THIS_DIR, 'baseline_model_viewer')
PICKLE_PATH = os.path.join(OUT_DIR, 'model_viewer_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'model_viewer_baseline_summary.txt')

# The sample model is isotropic, so only ``rho_x`` (== ``vals``) is
# meaningful; ``rho_y`` / ``rho_z`` and the anisotropic ratios need an
# anisotropic model, which we don't have a fixture for.
RHO_AXES = ('rho_x',)


# --------------------------------------------------------------------------- #
def _store(results, label, func):
    try:
        results[label] = func()
    except Exception as exc:  # noqa: BLE001
        results[label] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
        traceback.print_exc()


def _freeze_grid(results, prefix, grid):
    results['{}/dimensions'.format(prefix)] = np.asarray(grid.dimensions)
    results['{}/n_cells'.format(prefix)] = np.asarray(grid.n_cells)
    results['{}/bounds'.format(prefix)] = np.asarray(grid.bounds, dtype=float)
    for key in sorted(grid.cell_data.keys()):
        results['{}/cell_data/{}'.format(prefix, key)] = np.asarray(grid.cell_data[key], dtype=float)


def _make_resolution(model):
    """A stand-in resolution volume: the model's own cell sizes as a
    positive scalar field, so the ``resolution`` branch of
    ``model_to_rectgrid`` is exercised deterministically."""
    res = deepcopy(model)
    xcs = np.asarray(res.xCS)[:, None, None]
    ycs = np.asarray(res.yCS)[None, :, None]
    zcs = np.asarray(res.zCS)[None, None, :]
    res.vals = (xcs * ycs * zcs) / float(np.mean(res.xCS) * np.mean(res.yCS) * np.mean(res.zCS))
    return res


def collect():
    results = {}
    model = DS.Model(MODEL_FILE)
    results['model/dimensionality'] = np.asarray(str(model.dimensionality))
    results['model/nxnynz'] = np.asarray([model.nx, model.ny, model.nz])
    results['model/isotropic'] = np.asarray(bool(model.isotropic))

    for rho_axis in RHO_AXES:
        tag = rho_axis.replace('/', '_over_')
        _store(results, 'no_res/{}/_call'.format(tag),
               lambda ra=rho_axis: _freeze_grid(results, 'no_res/{}'.format(
                   ra.replace('/', '_over_')), model_to_rectgrid(model, rho_axis=ra)) or 'ok')

    res = _make_resolution(model)
    _store(results, 'with_res/rho_x/_call',
           lambda: _freeze_grid(results, 'with_res/rho_x',
                                model_to_rectgrid(model, resolution=res, rho_axis='rho_x')) or 'ok')

    return results


# --------------------------------------------------------------------------- #
def _describe(value):
    if isinstance(value, str):
        return repr(value)
    arr = np.asarray(value)
    if np.issubdtype(arr.dtype, np.number):
        flat = arr.ravel().astype(float)
        head = ', '.join('{:.6g}'.format(x) for x in flat[:4])
        return 'shape={} dtype={} head=[{}]'.format(arr.shape, arr.dtype, head)
    return 'shape={} dtype={} {}'.format(arr.shape, arr.dtype, arr.ravel()[:2])


def _diff(o, n, rtol=1e-9, atol=1e-12):
    if isinstance(o, str) or isinstance(n, str):
        return '' if o == n else '{!r} -> {!r}'.format(o, n)
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


def save():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = collect()
    with open(PICKLE_PATH, 'wb') as f:
        pickle.dump(results, f)
    with open(SUMMARY_PATH, 'w', encoding='utf-8') as f:
        for label in sorted(results):
            f.write('{:52s} {}\n'.format(label, _describe(results[label])))
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
    mism = []
    for label in sorted(set(old) & set(new)):
        msg = _diff(old[label], new[label])
        if msg:
            mism.append((label, msg))
    print('\n=== model_viewer baseline comparison ===')
    print('entries: {} old / {} new'.format(len(old), len(new)))
    for k in only_old:
        print('  - missing now: {}'.format(k))
    for k in only_new:
        print('  + new: {}'.format(k))
    for label, msg in mism:
        print('  * {}: {}'.format(label, msg))
    if not (only_old or only_new or mism):
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
        sys.exit('usage: python tests/model_viewer_baseline.py [save|compare]')
