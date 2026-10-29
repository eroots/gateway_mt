"""Before/after verification baseline for ``pyMT.utils``.

This script is part of the ``cleanup/utils`` refactor workflow.  It exercises the
public functions of :mod:`pyMT.utils` against real sample data (EDI files, a
ModEM data file, and -- if it loads -- a GoFEM data file) plus a set of frozen
synthetic inputs, and stores every result.

Usage
-----
Snapshot the current behaviour (run this on the pre-refactor code)::

    python tests/utils_baseline.py save

Re-run after refactoring and diff against the snapshot::

    python tests/utils_baseline.py compare

Outputs land in ``tests/baseline/``:

* ``utils_baseline.pkl``      -- pickled ``{label: value}`` mapping (the snapshot)
* ``utils_baseline_summary.txt`` -- human-readable shapes / sample values
* ``geotools_filter.png`` / ``compute_bost1D.png`` -- curve snapshots for eyeballing

The comparison uses ``numpy.allclose`` with a loose tolerance for anything that
passes through ``pyproj`` and a tight tolerance for the pure-numeric helpers.
"""

import os
import sys
import pickle
import traceback
from collections import OrderedDict

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pyMT.utils as utils
from pyMT.data_structures import RawData, Data


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(THIS_DIR, os.pardir, 'test_data', 'claude_test_data')
EDI_LIST = os.path.join(DATA_DIR, 'edis', 'allall.lst')
MODEM_FILE = os.path.join(DATA_DIR, 'swz_test_modem.dat')
GOFEM_FILE = os.path.join(DATA_DIR, 'swz_test_gofem.dat.gdat')
FREQ_FILE = os.path.join(DATA_DIR, 'frequencies')
RECEIVER_FILE = os.path.join(DATA_DIR, 'receivers_projected.csv')

OUT_DIR = os.path.join(THIS_DIR, 'baseline')
PICKLE_PATH = os.path.join(OUT_DIR, 'utils_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'utils_baseline_summary.txt')

# Components / invariants understood by compute_rho / compute_phase.
RHO_PHA_COMPS = ('xy', 'yx', 'xx', 'yy', 'det', 'gav', 'aav', 'ssq')
BOST_METHODS = ('phase', 'bostick', 'phase-bostick')

# Sites picked out for the per-site numeric snapshot (kept small + deterministic).
SAMPLE_SITES = ('SWZ001M', 'SWZ050M', 'e91003', 'g86001')


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _store(results, label, func, *args, **kwargs):
    """Call ``func`` and record the result (or the exception text) under ``label``."""
    try:
        value = func(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001 - we want every failure captured, not raised
        value = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
        traceback.print_exc()
    results[label] = value
    return value


def _frozen_synthetic_inputs():
    """Deterministic synthetic arrays for the pure-numeric helpers."""
    rng = np.random.RandomState(12345)
    periods = np.logspace(-3, 3, 40)
    smooth = 100 * (1 + 0.5 * np.sin(np.log10(periods)))
    spiky = smooth.copy()
    spiky[::7] *= 3.0
    return {
        'periods': periods,
        'edges': np.array([0.0, 1.0, 3.0, 7.0, 15.0, 31.0]),
        'smooth_curve': smooth,
        'spiky_curve': spiky,
        'noisy': smooth + rng.normal(scale=5.0, size=smooth.shape),
        'to_round': np.array([0.0, 1.23456789e-4, 9.87654321, -1234.56789, 3.14159265]),
        'locs': np.random.RandomState(7).uniform(-5e4, 5e4, size=(25, 2)),
        'vec_row': np.random.RandomState(3).uniform(size=(2, 9)),
        'dms_strings': ['48:14:26.50', '-114:03:11.1', '0:30:0', '145:0:0.0'],
        'dd_values': np.array([48.24069, -114.05308, 0.5, 145.0]),
        'wanted_p': np.array([1e-2, 1.0, 55.0, 900.0]),
    }


# --------------------------------------------------------------------------- #
# data loading
# --------------------------------------------------------------------------- #
def load_datasets():
    datasets = {}
    print('Loading EDIs from {} ...'.format(EDI_LIST))
    datasets['raw'] = RawData(listfile=EDI_LIST)
    print('  -> {} sites'.format(datasets['raw'].NS))

    print('Loading ModEM data from {} ...'.format(MODEM_FILE))
    datasets['modem'] = Data(datafile=MODEM_FILE, file_format='modem')
    print('  -> {} sites'.format(len(datasets['modem'].site_names)))

    try:
        print('Loading GoFEM data from {} ...'.format(GOFEM_FILE))
        datasets['gofem'] = Data(datafile=GOFEM_FILE, file_format='gofem')
        print('  -> {} sites'.format(len(datasets['gofem'].site_names)))
    except Exception as exc:  # noqa: BLE001
        print('  -> GoFEM load failed ({}); skipping that source.'.format(exc))
        datasets['gofem'] = None
    return datasets


# --------------------------------------------------------------------------- #
# the actual exercises
# --------------------------------------------------------------------------- #
def exercise_pure_numeric(results):
    fx = _frozen_synthetic_inputs()
    for key, arr in fx.items():
        results['input/{}'.format(key)] = arr

    _store(results, 'edge2center/edges', utils.edge2center, fx['edges'])
    _store(results, 'edge2center/periods', utils.edge2center, fx['periods'])
    _store(results, 'truncate/array', utils.truncate, fx['to_round'])
    _store(results, 'truncate/list', utils.truncate, list(fx['to_round']))
    _store(results, 'truncate/scalar', utils.truncate, 1234.56789)
    _store(results, 'normalize/default', utils.normalize, fx['noisy'])
    _store(results, 'normalize/explicit', utils.normalize, fx['noisy'],
           lower=10.0, upper=200.0, explicit_bounds=True)
    _store(results, 'normalize_range', utils.normalize_range, fx['noisy'],
           lower_range=0.0, upper_range=300.0, lower_norm=0.0, upper_norm=1.0)
    _store(results, 'percdiff', utils.percdiff, fx['smooth_curve'], fx['noisy'])
    _store(results, 'remove_outliers', utils.remove_outliers, fx['spiky_curve'])
    _store(results, 'remove_outliers/size3_t0.5', utils.remove_outliers,
           fx['spiky_curve'], size=3, threshold=0.5)
    _store(results, 'geotools_filter/loglog', utils.geotools_filter,
           fx['periods'], fx['spiky_curve'])
    _store(results, 'geotools_filter/linear', utils.geotools_filter,
           fx['periods'], fx['spiky_curve'], use_log=False)
    _store(results, 'geotools_filter/fwidth2', utils.geotools_filter,
           fx['periods'], fx['spiky_curve'], fwidth=2)
    _store(results, 'closest_periods', utils.closest_periods,
           fx['periods'], fx['wanted_p'])
    _store(results, 'linear_distance', utils.linear_distance,
           fx['locs'][:, 0], fx['locs'][:, 1])
    _store(results, 'flatten_list', utils.flatten_list, [[1, 2], [3, 4, 5], [6]])
    _store(results, 'skin_depth', utils.skin_depth, 100.0, fx['periods'])
    _store(results, 'is_all_empties/true', utils.is_all_empties, {'a': [], 'b': []})
    _store(results, 'is_all_empties/false', utils.is_all_empties, {'a': [], 'b': [1]})

    # geometry on locations
    _store(results, 'center_locs/locs', lambda a: utils.center_locs(a)[0], fx['locs'])
    _store(results, 'center_locs/center', lambda a: np.array(utils.center_locs(a)[1]), fx['locs'])
    _store(results, 'rotate_locs/30', utils.rotate_locs, fx['locs'], 30.0)
    _store(results, 'rotate_locs/-42.5', utils.rotate_locs, fx['locs'], -42.5)
    _store(results, 'project_to_line/20', utils.project_to_line, fx['locs'], 20.0)

    # DMS <-> DD
    _store(results, 'dms2dd', lambda s: np.array([utils.dms2dd(x) for x in s]),
           fx['dms_strings'])
    _store(results, 'dms_to_dd', lambda s: np.array([utils.dms_to_dd(x) for x in s]),
           fx['dms_strings'])
    _store(results, 'dd2dms', lambda v: np.array([utils.dd2dms(x) for x in v]),
           fx['dd_values'])
    _store(results, 'dd_to_dms', lambda v: np.array([utils.dd_to_dms(x) for x in v]),
           fx['dd_values'])
    _store(results, 'parse_dms', lambda s: np.array(utils.parse_dms(s)), '48:14:26.5')

    # projection helpers on a synthetic lon/lat cluster
    lon = np.array([-114.05, -114.00, -113.95, -113.90])
    lat = np.array([48.20, 48.24, 48.28, 48.32])
    results['input/lon'] = lon
    results['input/lat'] = lat
    _store(results, 'geo2utm', lambda: np.array(utils.geo2utm(lat[0], lon[0])))
    _store(results, 'compute_longorigin', lambda: np.array(
        [utils.compute_longorigin(L) for L in (0.0, 15.0, -21.0, 123.0)]))
    _store(results, 'zone_convert', lambda: np.array(
        [utils.zone_convert((x, y)) for x, y in zip(lon, lat)]))
    _store(results, 'project', lambda: np.array(
        [utils.project((x, y))[2:] for x, y in zip(lon, lat)]))
    _store(results, 'unproject', lambda: np.array(
        utils.unproject(11, 'T', 500000.0, 5343000.0)))
    _store(results, 'generic_projection/latlong_to_lambert', utils.generic_projection,
           lat, lon, 'latlong', 'lambert')
    _store(results, 'to_lambert', utils.to_lambert, lat, lon)

    # decorators / validation
    @utils.enforce_input(x=np.ndarray, n=int)
    def _enf(x=None, n=0):
        return np.asarray(x) * n

    _store(results, 'enforce_input', _enf, x=[1.0, 2.0, 3.0], n='4')
    _store(results, 'validate_input/ok', utils.validate_input, '3.5', float)
    _store(results, 'validate_input/bad', utils.validate_input, 'abc', float)
    _store(results, 'to_list/scalar', lambda v: np.array(utils.to_list(v)), 5)
    _store(results, 'to_list/iterable', lambda v: np.array(utils.to_list(v)), (1, 2, 3))
    _store(results, 'list2np', utils.list2np, [1.0, 2.0, 3.0])
    _store(results, 'np2list', lambda a: np.array(utils.np2list(a)), fx['periods'][:5])
    _store(results, 'check_extention/missing', utils.check_extention, 'foo', 'dat')
    _store(results, 'check_extention/present', utils.check_extention, 'foo.dat', 'dat')
    _store(results, 'fileparts', lambda p: np.array(utils.fileparts(p)),
           os.path.join('some', 'dir', 'file.edi'))

    # convert2impedance round-trip style call
    _store(results, 'convert2impedance/xy', lambda: np.array(utils.convert2impedance(
        fx['smooth_curve'], 45.0 + 0 * fx['smooth_curve'], fx['periods'], 'xy')))


def exercise_site_functions(results, tag, site, periods):
    """Run the MT computations that operate on a single Site."""
    for comp in RHO_PHA_COMPS:
        _store(results, '{}/compute_rho/{}'.format(tag, comp),
               lambda c=comp: np.array(utils.compute_rho(site, calc_comp=c)))
        _store(results, '{}/compute_rho_err/{}'.format(tag, comp),
               lambda c=comp: np.array(utils.compute_rho(site, calc_comp=c, errtype='used_error')))
        _store(results, '{}/compute_phase/{}'.format(tag, comp),
               lambda c=comp: np.array(utils.compute_phase(site, calc_comp=c)))
        _store(results, '{}/compute_phase_err/{}'.format(tag, comp),
               lambda c=comp: np.array(utils.compute_phase(site, calc_comp=c, errtype='used_error')))
    for method in BOST_METHODS:
        _store(results, '{}/compute_bost1D/{}'.format(tag, method),
               lambda m=method: np.array(utils.compute_bost1D(site, method=m, comp='det')))
    _store(results, '{}/compute_MT_determinant'.format(tag),
           lambda: np.array(utils.compute_MT_determinant(site)))
    _store(results, '{}/compute_gav'.format(tag), lambda: np.array(utils.compute_gav(site)))
    _store(results, '{}/compute_aav'.format(tag), lambda: np.array(utils.compute_aav(site)))
    _store(results, '{}/compute_ssq'.format(tag), lambda: np.array(utils.compute_ssq(site)))
    _store(results, '{}/skin_depth'.format(tag),
           lambda: utils.skin_depth(utils.compute_rho(site, calc_comp='det')[0], periods))


def exercise_misfit(results, raw):
    """calculate_misfit needs two Site objects; build a perturbed 'response'."""
    from copy import deepcopy
    name = SAMPLE_SITES[0]
    if name not in raw.sites:
        name = raw.site_names[0]
    data_site = raw.sites[name]
    resp_site = deepcopy(data_site)
    rng = np.random.RandomState(99)
    for comp in resp_site.components:
        resp_site.data[comp] = data_site.data[comp] * (1 + 0.05 * rng.standard_normal(
            data_site.data[comp].shape))
    try:
        misfit, comp_misfit = utils.calculate_misfit(data_site, resp_site)
        for comp, arr in misfit.items():
            results['calculate_misfit/misfit/{}'.format(comp)] = np.asarray(arr)
        for comp, arr in comp_misfit.items():
            results['calculate_misfit/comp_misfit/{}'.format(comp)] = np.asarray(arr)
    except Exception as exc:  # noqa: BLE001
        results['calculate_misfit'] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)


def exercise_projection_crosscheck(results, raw, modem):
    """Run real EDI lat/longs through project()/geo2utm and freeze the result.

    Also records the RawData projected X/Y and the ModEM local X/Y for the same
    sites so the geographic -> local pipeline can be eyeballed in the summary.
    """
    common = [s for s in raw.site_names if s in set(modem.site_names)][:12]
    results['crosscheck/common_sites'] = np.array(common)
    if not common:
        return
    lat = np.array([raw.sites[s].locations['Lat'] for s in common])
    lon = np.array([raw.sites[s].locations['Long'] for s in common])
    results['crosscheck/edi_lat'] = lat
    results['crosscheck/edi_lon'] = lon
    _store(results, 'crosscheck/project_xy',
           lambda: np.array([utils.project((lo, la))[2:] for lo, la in zip(lon, lat)]))
    _store(results, 'crosscheck/geo2utm_xy',
           lambda: np.array([utils.geo2utm(la, lo)[:2] for la, lo in zip(lat, lon)]))
    results['crosscheck/rawdata_xy'] = np.array(
        [[raw.sites[s].locations['Y'], raw.sites[s].locations['X']] for s in common])
    modem_locs = np.array(modem.get_locs())
    idx = [modem.site_names.index(s) for s in common]
    results['crosscheck/modem_xy'] = modem_locs[idx]


def exercise_regrid(results):
    model_file = os.path.join(THIS_DIR, os.pardir, 'test_data', 'testfile_ModEM.model')
    if not os.path.exists(model_file):
        results['regrid_model'] = 'SKIPPED: no model file'
        return
    from pyMT.data_structures import Model
    mod = Model(modelfile=model_file)
    new_x = utils.edge2center(np.asarray(mod.dx))[::2]
    new_y = utils.edge2center(np.asarray(mod.dy))[::2]
    new_z = utils.edge2center(np.asarray(mod.dz))[::2]
    # edge2center already halved the count; feed edge-like arrays back in
    _store(results, 'regrid_model',
           utils.regrid_model, mod,
           np.asarray(mod.dx)[::2], np.asarray(mod.dy)[::2], np.asarray(mod.dz)[::2])


def make_plots(results):
    fx_p = results.get('input/periods')
    if fx_p is not None and not isinstance(results.get('geotools_filter/loglog'), str):
        plt.figure(figsize=(7, 4))
        plt.loglog(fx_p, results['input/spiky_curve'], 'o-', label='input', alpha=0.5)
        plt.loglog(fx_p, results['geotools_filter/loglog'], 'r-', label='geotools_filter')
        plt.xlabel('period'); plt.ylabel('value'); plt.legend(); plt.title('geotools_filter')
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, 'geotools_filter.png'), dpi=110)
        plt.close()

    key = 'raw:{}/compute_bost1D/phase'.format(SAMPLE_SITES[0])
    bost = results.get(key)
    if bost is not None and not isinstance(bost, str):
        bostick, depth, rhofit, phase = bost
        plt.figure(figsize=(5, 6))
        plt.plot(bostick, depth, 'b-')
        plt.gca().invert_yaxis()
        plt.xscale('log')
        plt.xlabel('Bostick resistivity'); plt.ylabel('depth (km)')
        plt.title(key)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, 'compute_bost1D.png'), dpi=110)
        plt.close()


def collect():
    results = OrderedDict()
    datasets = load_datasets()
    raw = datasets['raw']
    modem = datasets['modem']

    exercise_pure_numeric(results)

    for name in SAMPLE_SITES:
        if name in raw.sites:
            exercise_site_functions(results, 'raw:{}'.format(name),
                                    raw.sites[name], raw.sites[name].periods)
    for name in SAMPLE_SITES:
        if name in getattr(modem, 'sites', {}):
            exercise_site_functions(results, 'modem:{}'.format(name),
                                    modem.sites[name], modem.sites[name].periods)
    if datasets.get('gofem') is not None:
        gofem = datasets['gofem']
        for name in list(gofem.sites)[:2]:
            exercise_site_functions(results, 'gofem:{}'.format(name),
                                    gofem.sites[name], gofem.sites[name].periods)

    exercise_misfit(results, raw)
    exercise_projection_crosscheck(results, raw, modem)
    exercise_regrid(results)

    make_plots(results)
    return results


# --------------------------------------------------------------------------- #
# save / compare
# --------------------------------------------------------------------------- #
def _describe(value):
    if isinstance(value, str):
        return value
    arr = np.asarray(value, dtype=object) if value is None else value
    try:
        arr = np.asarray(value)
        if arr.dtype == object:
            return 'object array shape={}'.format(arr.shape)
        flat = arr.ravel()
        if np.issubdtype(arr.dtype, np.complexfloating):
            head = ', '.join('{:.6g}{:+.6g}j'.format(x.real, x.imag) for x in flat[:4])
        elif np.issubdtype(arr.dtype, np.number):
            head = ', '.join('{:.6g}'.format(x) for x in flat[:4].astype(float))
        else:
            head = str(flat[:4])
        return 'shape={} dtype={} head=[{}]'.format(arr.shape, arr.dtype, head)
    except Exception:  # noqa: BLE001
        return repr(value)[:120]


def save():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = collect()
    with open(PICKLE_PATH, 'wb') as f:
        pickle.dump(results, f)
    with open(SUMMARY_PATH, 'w') as f:
        for label, value in results.items():
            f.write('{:60s} {}\n'.format(label, _describe(value)))
    n_exc = sum(1 for v in results.values() if isinstance(v, str) and v.startswith('EXCEPTION'))
    print('\nSaved {} entries to {}'.format(len(results), PICKLE_PATH))
    print('Summary -> {}'.format(SUMMARY_PATH))
    if n_exc:
        print('WARNING: {} entries recorded an exception (see summary).'.format(n_exc))


def _diff(o, n, rtol, atol):
    """Return a short mismatch description, or '' if o and n agree within tolerance.

    Handles numeric arrays, ragged/object arrays, and nested sequences.
    """
    # Base cases: strings and non-numeric scalars compare exactly.
    if isinstance(o, (str, bytes, np.str_)) or isinstance(n, (str, bytes, np.str_)):
        return '' if o == n else '{!r} -> {!r}'.format(o, n)
    if np.isscalar(o) and np.isscalar(n) and not isinstance(o, (int, float, complex,
                                                                np.number)):
        return '' if o == n else '{!r} -> {!r}'.format(o, n)
    try:
        o_arr = np.asarray(o, dtype=float)
        n_arr = np.asarray(n, dtype=float)
    except (TypeError, ValueError):
        o_seq, n_seq = np.asarray(o, dtype=object), np.asarray(n, dtype=object)
        if o_seq.shape != n_seq.shape:
            return 'object shape {} -> {}'.format(o_seq.shape, n_seq.shape)
        if o_seq.ndim == 0:
            return '' if o_seq.item() == n_seq.item() else '{!r} -> {!r}'.format(
                o_seq.item(), n_seq.item())
        flat_o, flat_n = o_seq.ravel(), n_seq.ravel()
        for i in range(flat_o.size):
            sub = _diff(flat_o[i], flat_n[i], rtol, atol)
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
    loose = ('project', 'unproject', 'geo2utm', 'zone_convert', 'lambert',
             'generic_projection', 'crosscheck')
    return (1e-6, 1e-6) if any(k in label for k in loose) else (1e-9, 1e-12)


def compare():
    if not os.path.exists(PICKLE_PATH):
        sys.exit('No baseline at {} -- run "save" first.'.format(PICKLE_PATH))
    with open(PICKLE_PATH, 'rb') as f:
        old = pickle.load(f)
    new = collect()

    old_keys, new_keys = set(old), set(new)
    only_old = sorted(old_keys - new_keys)
    only_new = sorted(new_keys - old_keys)
    mismatches = []
    exc_changes = []

    for label in sorted(old_keys & new_keys):
        o, n = old[label], new[label]
        o_exc = isinstance(o, str) and o.startswith('EXCEPTION')
        n_exc = isinstance(n, str) and n.startswith('EXCEPTION')
        if o_exc or n_exc:
            if o != n:
                exc_changes.append((label, o, n))
            continue
        if isinstance(o, str) or isinstance(n, str):
            if repr(o) != repr(n):
                mismatches.append((label, 'string/other value changed'))
            continue
        rtol, atol = _tol_for(label)
        msg = _diff(o, n, rtol, atol)
        if msg:
            mismatches.append((label, msg))

    print('\n=== utils baseline comparison ===')
    print('entries: {} old / {} new'.format(len(old), len(new)))
    if only_old:
        print('\nMissing now (were in baseline):')
        for k in only_old:
            print('  - {}'.format(k))
    if only_new:
        print('\nNew (not in baseline):')
        for k in only_new:
            print('  + {}'.format(k))
    if exc_changes:
        print('\nException status changed:')
        for label, o, n in exc_changes:
            print('  ! {}\n      old: {}\n      new: {}'.format(label, o, n))
    if mismatches:
        print('\nNumeric/value mismatches:')
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
        sys.exit('usage: python tests/utils_baseline.py [save|compare]')
