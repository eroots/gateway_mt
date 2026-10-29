"""Before/after verification baseline for ``gateway_mt.data_structures``.

Part of the ``cleanup/data_structures`` refactor workflow.  Loads the sample
datasets through every ``data_structures`` read path, exercises the public API
of ``RawData`` / ``Data`` / ``Model`` / ``Dataset`` / ``Site`` /
``PhaseTensor`` / ``CART``, and stores every result for ``numpy.allclose``
comparison after refactoring.

Usage
-----
    python tests/data_structures_baseline.py save      # snapshot current behaviour
    python tests/data_structures_baseline.py compare    # re-run and diff

Outputs land in ``tests/baseline_ds/``:

* ``ds_baseline.pkl``          -- pickled ``{label: value}`` snapshot
* ``ds_baseline_summary.txt``  -- human-readable shapes / sample values
"""

import os
import sys
import pickle
import traceback
import warnings
from collections import OrderedDict
from copy import deepcopy

import numpy as np

import gateway_mt.data_structures as ds
from gateway_mt.data_structures import RawData, Data, Model, Dataset

warnings.filterwarnings('ignore')

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.join(THIS_DIR, 'data', 'synthetic')

EDI_LIST = os.path.join(SYN, 'edis', 'synthetic.lst')
MODEM_FILE = os.path.join(SYN, 'synthetic_modem.dat')
GOFEM_FILE = os.path.join(SYN, 'synthetic_gofem.gdat')
MODEL_FILE = os.path.join(SYN, 'synthetic.model')
WS_FILE = os.path.join(SYN, 'synthetic_ws.data')
MARE_FILE = os.path.join(SYN, 'synthetic_mare2dem.emdata')

OUT_DIR = os.path.join(THIS_DIR, 'baseline_ds')
PICKLE_PATH = os.path.join(OUT_DIR, 'ds_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'ds_baseline_summary.txt')

SAMPLE_SITES = ('SYN001', 'SYN015', 'SYN029', 'SYN030')
PT_DIFF_DEFS = ('booker', 'heise', 'hill', 'castro', 'roots_1', 'roots_2')


# --------------------------------------------------------------------------- #
# snapshot machinery (shared shape with tests/utils_baseline.py)
# --------------------------------------------------------------------------- #
def _store(results, label, func, *args, **kwargs):
    try:
        value = func(*args, **kwargs)
    except Exception as exc:  # noqa: BLE001
        value = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
        traceback.print_exc()
    results[label] = value
    return value


def _site_dicts(results, prefix, site):
    """Freeze the per-site array dicts and scalar props of a Site."""
    for attr in ('data', 'errors', 'used_error', 'errmap', 'flags'):
        d = getattr(site, attr, None)
        if not isinstance(d, dict):
            continue
        for comp in sorted(d):
            results['{}/{}/{}'.format(prefix, attr, comp)] = np.asarray(d[comp])
    results['{}/periods'.format(prefix)] = np.asarray(site.periods)
    results['{}/components'.format(prefix)] = np.array(sorted(site.components))
    for prop in ('NP', 'NR', 'swift_skew', 'det_rho', 'ssq_rho'):
        def _get(p=prop):
            v = getattr(site, p)
            return np.asarray(v() if callable(v) else v)
        _store(results, '{}/{}'.format(prefix, prop), _get)


def _pt_params(results, prefix, pt):
    for name in ('phi', 'phi_error', 'det_phi', 'skew_phi', 'phi_1', 'phi_2', 'phi_3',
                 'phi_max', 'phi_min', 'alpha', 'beta', 'Lambda', 'azimuth'):
        if hasattr(pt, name):
            _store(results, '{}/{}'.format(prefix, name), lambda n=name: np.asarray(getattr(pt, n)))


# --------------------------------------------------------------------------- #
# exercises
# --------------------------------------------------------------------------- #
def exercise_rawdata(results):
    rd = RawData(listfile=EDI_LIST)
    results['raw/NS'] = np.asarray(rd.NS)
    results['raw/site_names'] = np.array(rd.site_names)
    results['raw/locations'] = np.asarray(rd.locations)
    results['raw/origin'] = np.asarray(rd.origin)
    results['raw/azimuth'] = np.asarray(rd.azimuth)
    results['raw/master_periods'] = np.asarray(sorted(rd.master_periods))
    results['raw/narrow_periods'] = np.asarray(sorted(rd.narrow_periods))
    for azi in (0, 30, -45):
        _store(results, 'raw/get_locs/azi{}'.format(azi), lambda a=azi: np.asarray(rd.get_locs(azi=a)))
    _store(results, 'raw/average_rho', lambda: np.asarray(rd.average_rho(comp='det')[0]))

    for name in SAMPLE_SITES:
        if name in rd.sites:
            _site_dicts(results, 'raw:{}'.format(name), rd.sites[name])
            site = rd.sites[name]
            if site.phase_tensors:
                for idx in (0, len(site.phase_tensors) // 2, len(site.phase_tensors) - 1):
                    _pt_params(results, 'raw:{}/PT[{}]'.format(name, idx), site.phase_tensors[idx])

    # RawData.get_data -> the {name: Site} dict a Data object is built from
    _store(results, 'raw/get_data/keys', lambda: np.array(sorted(rd.get_data(sites=list(SAMPLE_SITES[:2])))))
    got = None
    try:
        got = rd.get_data(sites=list(SAMPLE_SITES[:2]))
    except Exception:  # noqa: BLE001
        pass
    if got:
        for name, site in got.items():
            _site_dicts(results, 'raw_get_data:{}'.format(name), site)
    return rd


def exercise_data(results, tag, datafile, file_format):
    d = Data(datafile=datafile, file_format=file_format)
    results['{}/NP'.format(tag)] = np.asarray(d.NP)
    results['{}/NS'.format(tag)] = np.asarray(d.NS)
    results['{}/NR'.format(tag)] = np.asarray(d.NR)
    results['{}/periods'.format(tag)] = np.asarray(d.periods)
    results['{}/components'.format(tag)] = np.array(sorted(d.components))
    results['{}/site_names'.format(tag)] = np.array(d.site_names)
    _store(results, '{}/inv_type'.format(tag), lambda: np.asarray(d.inv_type))
    _store(results, '{}/used_components'.format(tag), lambda: np.array(sorted(d.used_components)))
    _store(results, '{}/has_flagged_data'.format(tag), lambda: np.asarray(d.has_flagged_data))
    _store(results, '{}/locations'.format(tag), lambda: np.asarray(d.get_locs()))
    _store(results, '{}/average_rho'.format(tag), lambda: np.asarray(d.average_rho()[0]))
    for name in SAMPLE_SITES:
        if name in getattr(d, 'sites', {}):
            _site_dicts(results, '{}:{}'.format(tag, name), d.sites[name])
    return d


def exercise_data_ops(results, d, tag):
    """Rotation / error-floor / outlier ops on copies, so the base object is untouched."""
    rot = deepcopy(d)
    _store(results, '{}/rotate_sites(30)/locations'.format(tag),
           lambda: (rot.rotate_sites(30), np.asarray(rot.get_locs()))[1])
    name = next((s for s in SAMPLE_SITES if s in rot.sites), rot.site_names[0])
    for comp in sorted(rot.sites[name].components):
        results['{}/rotate30:{}/{}'.format(tag, name, comp)] = np.asarray(rot.sites[name].data[comp])

    ef = deepcopy(d)
    _store(results, '{}/apply_error_floor'.format(tag), lambda: ef.apply_error_floor())
    for comp in sorted(ef.sites[name].components):
        results['{}/errfloor:{}/{}'.format(tag, name, comp)] = np.asarray(ef.sites[name].used_error[comp])

    site = deepcopy(d.sites[name])
    _store(results, '{}/site.calculate_error_floor'.format(tag),
           lambda: np.asarray([np.asarray(v) for v in site.calculate_error_floor().values()], dtype=object))
    _store(results, '{}/site.detect_outliers'.format(tag),
           lambda: (site.detect_outliers(10), np.asarray([np.asarray(v) for v in site.errmap.values()],
                                                         dtype=object))[1])
    _store(results, '{}/site.get_data'.format(tag),
           lambda: np.asarray([np.asarray(v) for v in
                               d.sites[name].get_data().data.values()], dtype=object))


def exercise_pt_diffs(results, rd):
    """All PhaseTensor.__sub__ difference definitions on a real tensor pair."""
    a_name, b_name = SAMPLE_SITES[0], SAMPLE_SITES[1]
    if a_name not in rd.sites or b_name not in rd.sites:
        return
    sa, sbt = rd.sites[a_name], rd.sites[b_name]
    if not sa.phase_tensors or not sbt.phase_tensors:
        return
    n = min(len(sa.phase_tensors), len(sbt.phase_tensors))
    idx = n // 2
    for ddef in PT_DIFF_DEFS:
        def _do(ddef=ddef):
            pa = deepcopy(sa.phase_tensors[idx])
            pb = deepcopy(sbt.phase_tensors[idx])
            pa.difference_definition = ddef
            diff = pa - pb
            return np.asarray([np.asarray(diff.phi)] +
                              [np.asarray(getattr(diff, k)) for k in sorted(vars(diff))
                               if k.startswith('diff_') and np.isscalar(getattr(diff, k))],
                              dtype=object)
        _store(results, 'PT_sub/{}'.format(ddef), _do)


def exercise_model(results):
    m = Model(modelfile=MODEL_FILE)
    for a in ('nx', 'ny', 'nz'):
        results['model/{}'.format(a)] = np.asarray(getattr(m, a))
    for a in ('dx', 'dy', 'dz', 'xCS', 'yCS', 'zCS'):
        results['model/{}'.format(a)] = np.asarray(getattr(m, a))
    results['model/center'] = np.asarray(m.center)
    results['model/vals_stats'] = np.asarray([np.min(m.vals), np.max(m.vals),
                                              np.mean(m.vals), np.std(m.vals),
                                              len(np.unique(m.vals))])
    _store(results, 'model/cell_centers', lambda: np.asarray(m.cell_centers(), dtype=object))
    _store(results, 'model/is_half_space', lambda: np.asarray(m.is_half_space()))
    _store(results, 'model/conductance(0,10)', lambda: np.asarray(m.conductance(0, 10)))
    _store(results, 'model/depth_weighted_rho(0,10)', lambda: np.asarray(m.depth_weighted_rho(0, 10)))
    _store(results, 'model/elevation', lambda: np.asarray(m.elevation))

    ins = deepcopy(m)
    _store(results, 'model/dz_insert', lambda: (ins.dz_insert(1, (ins.dz[0] + ins.dz[1]) / 2),
                                                np.asarray(ins.dz))[1])
    dele = deepcopy(m)
    _store(results, 'model/dx_delete', lambda: (dele.dx_delete(1), np.asarray(dele.dx))[1])

    loc = deepcopy(m)
    _store(results, 'model/to_local/dx', lambda: (loc.to_local(), np.asarray(loc.dx))[1])
    return m


def exercise_dataset_rms(results):
    dset = Dataset(datafile=MODEM_FILE, modelfile=MODEL_FILE)
    results['dataset/azimuth'] = np.asarray(dset.azimuth)
    results['dataset/site_names'] = np.array(list(dset.site_names or []))
    # Build a deterministic pseudo-response so calculate_RMS has something to chew on.
    dset.response = deepcopy(dset.data)
    rng = np.random.RandomState(1234)
    for s in dset.response.site_names:
        for comp in dset.response.sites[s].components:
            arr = dset.response.sites[s].data[comp]
            dset.response.sites[s].data[comp] = arr * (1 + 0.05 * rng.standard_normal(arr.shape))
    try:
        rms = dset.calculate_RMS()
    except Exception as exc:  # noqa: BLE001
        rms = None
        results['dataset/calculate_RMS'] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
    if isinstance(rms, dict):
        results['dataset/RMS/Total'] = np.asarray(rms['Total'])
        results['dataset/RMS/Component'] = np.asarray(
            [rms['Component'][c] for c in sorted(rms['Component'])])
        results['dataset/RMS/Period_Total'] = np.asarray(rms['Period']['Total'])
        first_site = sorted(rms['Station'])[0]
        results['dataset/RMS/Station[0]/Total'] = np.asarray(rms['Station'][first_site]['Total'])


def exercise_roundtrip(results, tag, d, file_format, ext):
    out = os.path.join(OUT_DIR, 'roundtrip_{}{}'.format(tag, ext))
    try:
        d.write(out, file_format=file_format)
        d2 = Data(datafile=out, file_format=file_format)
    except Exception as exc:  # noqa: BLE001
        results['roundtrip/{}'.format(tag)] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)
        return
    results['roundtrip/{}/periods_close'.format(tag)] = np.asarray(
        np.allclose(np.sort(d.periods), np.sort(d2.periods)))
    name = next((s for s in SAMPLE_SITES if s in d2.sites), d2.site_names[0])
    for comp in sorted(d2.sites[name].components):
        results['roundtrip/{}/{}/{}'.format(tag, name, comp)] = np.asarray(d2.sites[name].data[comp])
        results['roundtrip/{}/{}/{}_err'.format(tag, name, comp)] = np.asarray(
            d2.sites[name].used_error[comp])


def collect():
    results = OrderedDict()
    rd = exercise_rawdata(results)
    exercise_pt_diffs(results, rd)

    modem = exercise_data(results, 'modem', MODEM_FILE, 'modem')
    exercise_data_ops(results, modem, 'modem')
    exercise_roundtrip(results, 'modem', modem, 'modem', '.dat')

    try:
        gofem = exercise_data(results, 'gofem', GOFEM_FILE, 'gofem')
        exercise_data_ops(results, gofem, 'gofem')
    except Exception as exc:  # noqa: BLE001
        results['gofem/load'] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)

    for tag, path, fmt in (('ws', WS_FILE, 'wsinv3dmt'), ('mare', MARE_FILE, 'mare2dem')):
        if os.path.exists(path):
            try:
                exercise_data(results, tag, path, fmt)
            except Exception as exc:  # noqa: BLE001
                results['{}/load'.format(tag)] = 'EXCEPTION: {}: {}'.format(type(exc).__name__, exc)

    exercise_model(results)
    exercise_dataset_rms(results)
    return results


# --------------------------------------------------------------------------- #
# save / compare  (same logic as tests/utils_baseline.py)
# --------------------------------------------------------------------------- #
def _describe(value):
    if isinstance(value, str):
        return value
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


def _diff(o, n, rtol, atol):
    if isinstance(o, (str, bytes, np.str_)) or isinstance(n, (str, bytes, np.str_)):
        return '' if o == n else '{!r} -> {!r}'.format(o, n)
    if np.isscalar(o) and np.isscalar(n) and not isinstance(o, (int, float, complex, np.number)):
        return '' if o == n else '{!r} -> {!r}'.format(o, n)
    try:
        o_arr = np.asarray(o, dtype=float)
        n_arr = np.asarray(n, dtype=float)
    except (TypeError, ValueError):
        o_seq, n_seq = np.asarray(o, dtype=object), np.asarray(n, dtype=object)
        if o_seq.shape != n_seq.shape:
            return 'object shape {} -> {}'.format(o_seq.shape, n_seq.shape)
        if o_seq.ndim == 0:
            try:
                same = bool(o_seq.item() == n_seq.item())
            except (ValueError, TypeError):
                same = repr(o_seq.item()) == repr(n_seq.item())
            return '' if same else 'scalar value changed'
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
    loose = ('locations', 'get_locs', 'origin', 'to_local', 'to_UTM', 'rotate',
             'project', 'center', 'PT', 'phase_tensor', 'RMS', 'conductance', 'average_rho')
    return (1e-6, 1e-6) if any(k in label for k in loose) else (1e-9, 1e-12)


def save():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = collect()
    with open(PICKLE_PATH, 'wb') as f:
        pickle.dump(results, f)
    with open(SUMMARY_PATH, 'w') as f:
        for label, value in results.items():
            f.write('{:64s} {}\n'.format(label, _describe(value)))
    n_exc = sum(1 for v in results.values() if isinstance(v, str) and v.startswith('EXCEPTION'))
    print('\nSaved {} entries to {}'.format(len(results), PICKLE_PATH))
    print('Summary -> {}'.format(SUMMARY_PATH))
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

    print('\n=== data_structures baseline comparison ===')
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
        sys.exit('usage: python tests/data_structures_baseline.py [save|compare]')
