"""Before/after verification baseline for ``pyMT.IO``.

Part of the ``cleanup/IO`` refactor workflow.  Exercises the reader and
writer functions of :mod:`pyMT.IO` against the sample datasets (and a few
generated fixtures), and stores every result for ``numpy.allclose`` /
text comparison after refactoring.

Usage
-----
    python tests/io_baseline.py save      # snapshot current behaviour
    python tests/io_baseline.py compare    # re-run and diff

Outputs land in ``tests/baseline_io/``:

* ``io_baseline.pkl``         -- pickled ``{label: value}`` snapshot
* ``io_baseline_summary.txt`` -- human-readable shapes / sample values
* ``scratch/``                -- files written by the round-trip checks
"""

import os
import re
import sys
import pickle
import shutil
import zlib
import traceback
import warnings
from collections import OrderedDict

import numpy as np

import pyMT.IO as IO
from pyMT.data_structures import Data, Model

warnings.filterwarnings('ignore')

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.join(THIS_DIR, 'data', 'synthetic')

EDI_LIST = os.path.join(SYN, 'edis', 'synthetic.lst')
EDI_DIR = os.path.join(SYN, 'edis')
MODEM_FILE = os.path.join(SYN, 'synthetic_modem.dat')
GOFEM_FILE = os.path.join(SYN, 'synthetic_gofem.gdat')
WS_FILE = os.path.join(SYN, 'synthetic_ws.data')
MARE_FILE = os.path.join(SYN, 'synthetic_mare2dem.emdata')
MODEL_FILE = os.path.join(SYN, 'synthetic.model')

OUT_DIR = os.path.join(THIS_DIR, 'baseline_io')
SCRATCH = os.path.join(OUT_DIR, 'scratch')
PICKLE_PATH = os.path.join(OUT_DIR, 'io_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'io_baseline_summary.txt')

SAMPLE_SITES = ('SYN001', 'SYN015', 'SYN029', 'SYN030')

# Volatile substrings to blank out of written text before snapshotting.
_DATE_RE = re.compile(r'\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(\.\d+)?')
_SCRATCH_RE = None  # set in main once SCRATCH is known


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


def _norm_text(path):
    """Read a written file and blank out volatile bits for stable comparison."""
    with open(path, 'r', errors='replace') as f:
        txt = f.read()
    txt = _DATE_RE.sub('<DATE>', txt)
    txt = txt.replace(SCRATCH.replace('\\', '/'), '<SCRATCH>')
    txt = txt.replace(SCRATCH, '<SCRATCH>')
    return txt


def _freeze_sites(results, prefix, sites):
    """Freeze the {name: {data/errors/periods/locations/...}} dict from read_data / read_raw_data."""
    results['{}/site_names'.format(prefix)] = np.array(sorted(sites))
    for name in SAMPLE_SITES:
        if name not in sites:
            continue
        s = sites[name]
        for key in ('data', 'errors', 'errmap'):
            d = s.get(key) if isinstance(s, dict) else None
            if not isinstance(d, dict):
                continue
            for comp in sorted(d):
                results['{}:{}/{}/{}'.format(prefix, name, key, comp)] = np.asarray(d[comp])
        results['{}:{}/periods'.format(prefix, name)] = np.asarray(s['periods'])
        results['{}:{}/azimuth'.format(prefix, name)] = np.asarray(s.get('azimuth', np.nan))
        locs = s.get('locations', {})
        for lk in sorted(locs):
            try:
                results['{}:{}/loc/{}'.format(prefix, name, lk)] = np.asarray(locs[lk], dtype=float)
            except (TypeError, ValueError):
                results['{}:{}/loc/{}'.format(prefix, name, lk)] = np.asarray(str(locs[lk]))


# --------------------------------------------------------------------------- #
# exercises
# --------------------------------------------------------------------------- #
def exercise_read_data(results):
    for tag, path, fmt in (('modem', MODEM_FILE, 'modem'),
                           ('gofem', GOFEM_FILE, 'gofem'),
                           ('ws', WS_FILE, 'wsinv3dmt'),
                           ('mare', MARE_FILE, 'mare2dem')):
        def _do(path=path, fmt=fmt):
            sites, other = IO.read_data(datafile=path, file_format=fmt)
            return sites, other
        got = _store(results, 'read_data/{}/_call'.format(tag), _do)
        if isinstance(got, tuple):
            sites, other = got
            _freeze_sites(results, 'read_data/{}'.format(tag), sites)
            for k in sorted(other):
                try:
                    results['read_data/{}/other/{}'.format(tag, k)] = np.asarray(other[k], dtype=float)
                except (TypeError, ValueError):
                    results['read_data/{}/other/{}'.format(tag, k)] = np.asarray(str(other[k]))
        # drop the raw tuple (unpicklable-friendly / noisy)
        if not isinstance(results.get('read_data/{}/_call'.format(tag)), str):  # keep EXCEPTION entries
            results.pop('read_data/{}/_call'.format(tag), None)


def exercise_read_raw_data(results):
    names = IO.read_sites(EDI_LIST)
    _store(results, 'read_raw_data/_names', lambda: np.array(sorted(names)))
    def _do():
        return IO.read_raw_data(names, datpath=EDI_DIR)
    got = _store(results, 'read_raw_data/_call', _do)
    if isinstance(got, dict):
        _freeze_sites(results, 'read_raw_data', got)
    if not isinstance(results.get('read_raw_data/_call'), str):  # keep EXCEPTION entries
        results.pop('read_raw_data/_call', None)


def exercise_read_model(results):
    def _do():
        return IO.read_model(modelfile=MODEL_FILE, file_format='modem3d')
    got = _store(results, 'read_model/_call', _do)
    if isinstance(got, tuple):
        mod, dim = got
        results['read_model/dimensionality'] = np.asarray(str(dim))
        for k in ('xCS', 'yCS', 'zCS'):
            results['read_model/{}'.format(k)] = np.asarray(mod[k], dtype=float)
        v = np.asarray(mod['vals'], dtype=float)
        results['read_model/vals_stats'] = np.asarray(
            [v.min(), v.max(), v.mean(), v.std(), len(np.unique(v)), v.shape[0], v.shape[1], v.shape[2]])
    if not isinstance(results.get('read_model/_call'), str):  # keep EXCEPTION entries
        results.pop('read_model/_call', None)


def exercise_misc_readers(results):
    _store(results, 'read_sites/list', lambda: np.array(sorted(IO.read_sites(EDI_LIST))))
    _store(results, 'read_sites/single_edi',
           lambda: np.array(IO.read_sites(os.path.join(EDI_DIR, 'SYN001.edi'))))
    _store(results, 'read_data_header/modem',
           lambda: np.asarray(IO.read_data_header(MODEM_FILE, file_format='modem')))
    _store(results, 'read_data_header/ws',
           lambda: np.asarray(IO.read_data_header(WS_FILE, file_format='wsinv3dmt')))
    _store(results, 'read_data_header/mare',
           lambda: np.asarray(IO.read_data_header(MARE_FILE, file_format='mare2dem')))
    for it in (1, 2, 3, 4, 5):
        _store(results, 'get_components/invType={}'.format(it),
               lambda it=it: np.array(IO.get_components(invType=it)[0]))
    for nr in (4, 8, 12):
        _store(results, 'get_components/NR={}'.format(nr),
               lambda nr=nr: np.array(IO.get_components(NR=nr)[0]))
    # read_freqset on a generated freqset file
    fdir = os.path.join(SCRATCH, 'fs')
    os.makedirs(fdir, exist_ok=True)
    with open(os.path.join(fdir, 'freqset'), 'w') as f:
        f.write('5\n')
        for p in (0.01, 0.1, 1.0, 10.0, 100.0):
            f.write('{}\n'.format(p))
    _store(results, 'read_freqset', lambda: np.asarray(IO.read_freqset(fdir + os.sep)))


def exercise_write_roundtrips(results):
    specs = (('modem', MODEM_FILE, 'modem', '.dat'),
             ('ws', WS_FILE, 'wsinv3dmt', '.data'),
             ('gofem', GOFEM_FILE, 'gofem', '.gdat'),
             ('mare', MARE_FILE, 'mare2dem', '.emdata'))
    # GoFEM re-reading needs a receiver file next to the .gdat
    _RECV = os.path.join(SYN, 'receivers_projected.csv')
    if os.path.exists(_RECV):
        shutil.copy(_RECV, os.path.join(SCRATCH, 'receivers_projected.csv'))
    for tag, path, fmt, ext in specs:
        out = os.path.join(SCRATCH, 'rt_{}{}'.format(tag, ext))
        def _do(path=path, fmt=fmt, out=out):
            d = Data(datafile=path, file_format=fmt)
            IO.write_data(data=d, outfile=out, file_format=fmt)
            sites, other = IO.read_data(datafile=out, file_format=fmt)
            return d, sites, other
        got = _store(results, 'roundtrip/{}/_call'.format(tag), _do)
        if isinstance(got, tuple):
            d, sites, other = got
            name = next((s for s in SAMPLE_SITES if s in sites), sorted(sites)[0])
            for comp in sorted(sites[name].get('data', {})):
                results['roundtrip/{}/{}/{}'.format(tag, name, comp)] = np.asarray(sites[name]['data'][comp])
                results['roundtrip/{}/{}/{}_err'.format(tag, name, comp)] = np.asarray(sites[name]['errors'][comp])
            results['roundtrip/{}/periods'.format(tag)] = np.asarray(sites[name]['periods'])
            if os.path.exists(out):
                results['roundtrip/{}/written_text'.format(tag)] = _norm_text(out)
        if not isinstance(results.get('roundtrip/{}/_call'.format(tag)), str):  # keep EXCEPTION entries
            results.pop('roundtrip/{}/_call'.format(tag), None)


def exercise_write_model_roundtrip(results):
    out = os.path.join(SCRATCH, 'rt_model.model')
    def _do():
        m = Model(modelfile=MODEL_FILE)
        IO.write_model(m, out, file_format='modem')
        mod, dim = IO.read_model(modelfile=out, file_format='modem3d')
        return mod, dim
    got = _store(results, 'roundtrip/model/_call', _do)
    if isinstance(got, tuple):
        mod, dim = got
        for k in ('xCS', 'yCS', 'zCS'):
            results['roundtrip/model/{}'.format(k)] = np.asarray(mod[k], dtype=float)
        v = np.asarray(mod['vals'], dtype=float)
        results['roundtrip/model/vals_stats'] = np.asarray([v.min(), v.max(), v.mean(), v.std()])
        if os.path.exists(out):
            results['roundtrip/model/written_text_head'] = _norm_text(out)[:4000]
    if not isinstance(results.get('roundtrip/model/_call'), str):  # keep EXCEPTION entries
        results.pop('roundtrip/model/_call', None)


def exercise_covariance_roundtrip(results):
    out = os.path.join(SCRATCH, 'rt.cov')
    nx, ny, nz = 6, 7, 8
    def _do():
        IO.write_covariance(out, NX=nx, NY=ny, NZ=nz, exceptions=None,
                            sigma_x=0.3, sigma_y=0.4, sigma_z=0.5, num_smooth=2)
        return IO.read_covariance(out)
    got = _store(results, 'covariance/_call', _do)
    if isinstance(got, tuple):
        rnx, rny, rnz, sx, sy, sz, nsm, cov_exc = got
        results['covariance/dims'] = np.asarray([rnx, rny, rnz, nsm])
        results['covariance/sigmas'] = np.asarray([np.ravel(sx)[0], np.ravel(sy)[0],
                                                   np.ravel(sz)[0]])
        if os.path.exists(out):
            results['covariance/written_text'] = _norm_text(out)
    if not isinstance(results.get('covariance/_call'), str):  # keep EXCEPTION entries
        results.pop('covariance/_call', None)


def exercise_vtk_and_csv(results):
    d = Data(datafile=MODEM_FILE, file_format='modem')
    m = Model(modelfile=MODEL_FILE)
    origin = (500000.0, 5340000.0)

    mv = os.path.join(SCRATCH, 'model')          # model_to_vtk appends '_model.vtk'
    _store(results, 'model_to_vtk/_call',
           lambda: IO.model_to_vtk(m, outfile=mv, origin=origin, UTM='16N'))
    mvf = mv + '_model.vtk'
    if os.path.exists(mvf):
        results['model_to_vtk/text_head'] = _norm_text(mvf)[:6000]
    if not isinstance(results.get('model_to_vtk/_call'), str):  # keep EXCEPTION entries
        results.pop('model_to_vtk/_call', None)

    sv = os.path.join(SCRATCH, 'sites')
    _store(results, 'sites_to_vtk/_call',
           lambda: IO.sites_to_vtk(d, origin=origin, outfile=sv, UTM='16N'))
    svf = sv + '_sites.vtk'
    if os.path.exists(svf):
        results['sites_to_vtk/text'] = _norm_text(svf)
    if not isinstance(results.get('sites_to_vtk/_call'), str):  # keep EXCEPTION entries
        results.pop('sites_to_vtk/_call', None)

    pt = os.path.join(SCRATCH, 'pt')
    _store(results, 'write_phase_tensors/_call',
           lambda: IO.write_phase_tensors(d, pt, verbose=True))
    if os.path.exists(pt + '.csv'):
        results['write_phase_tensors/text'] = _norm_text(pt + '.csv')
    if not isinstance(results.get('write_phase_tensors/_call'), str):  # keep EXCEPTION entries
        results.pop('write_phase_tensors/_call', None)

    ia = os.path.join(SCRATCH, 'ia')
    _store(results, 'write_induction_arrows/_call',
           lambda: IO.write_induction_arrows(d, ia, verbose=True))
    if os.path.exists(ia + '.csv'):
        results['write_induction_arrows/text'] = _norm_text(ia + '.csv')
    if not isinstance(results.get('write_induction_arrows/_call'), str):  # keep EXCEPTION entries
        results.pop('write_induction_arrows/_call', None)

    loc = os.path.join(SCRATCH, 'locs.csv')
    _store(results, 'write_locations/_call',
           lambda: IO.write_locations(d, out_file=loc, file_format='csv', verbose=1))
    if os.path.exists(loc):
        results['write_locations/text'] = _norm_text(loc)
    if not isinstance(results.get('write_locations/_call'), str):  # keep EXCEPTION entries
        results.pop('write_locations/_call', None)

    # write_edi for one site (from raw data, which carries Lat/Long)
    from pyMT.data_structures import RawData
    rd = RawData(listfile=EDI_LIST)
    sname = next((s for s in SAMPLE_SITES if s in rd.sites), rd.site_names[0])
    edi = os.path.join(SCRATCH, sname + '_out.edi')
    _store(results, 'write_edi/_call', lambda: IO.write_edi(rd.sites[sname], edi))
    if os.path.exists(edi):
        results['write_edi/text'] = _norm_text(edi)
    if not isinstance(results.get('write_edi/_call'), str):  # keep EXCEPTION entries
        results.pop('write_edi/_call', None)

    # write_response (legacy WS format) + write_list
    resp = os.path.join(SCRATCH, 'resp.dat')
    _store(results, 'write_response/_call', lambda: IO.write_response(d, resp))
    if os.path.exists(resp):
        results['write_response/text'] = _norm_text(resp)
    if not isinstance(results.get('write_response/_call'), str):  # keep EXCEPTION entries
        results.pop('write_response/_call', None)

    lst = os.path.join(SCRATCH, 'out.lst')
    _store(results, 'write_list/_call', lambda: IO.write_list(d, lst))
    if os.path.exists(lst):
        results['write_list/text'] = _norm_text(lst)
    if not isinstance(results.get('write_list/_call'), str):  # keep EXCEPTION entries
        results.pop('write_list/_call', None)


def collect():
    results = OrderedDict()
    if os.path.isdir(SCRATCH):
        shutil.rmtree(SCRATCH)
    os.makedirs(SCRATCH, exist_ok=True)
    exercise_read_data(results)
    exercise_read_raw_data(results)
    exercise_read_model(results)
    exercise_misc_readers(results)
    exercise_write_roundtrips(results)
    exercise_write_model_roundtrip(results)
    exercise_covariance_roundtrip(results)
    exercise_vtk_and_csv(results)
    return results


# --------------------------------------------------------------------------- #
# save / compare  (same logic as the other baseline harnesses)
# --------------------------------------------------------------------------- #
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
    if isinstance(o, str) or isinstance(n, str):
        if o == n:
            return ''
        if isinstance(o, str) and isinstance(n, str):
            # locate first differing line for readability
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
    loose = ('loc/', 'other/origin', 'vtk', 'model_to_vtk', 'sites_to_vtk',
             'phase_tensors', 'induction_arrows', 'write_locations', 'write_edi')
    return (1e-6, 1e-6) if any(k in label for k in loose) else (1e-9, 1e-12)


def save():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = collect()
    with open(PICKLE_PATH, 'wb') as f:
        pickle.dump(results, f)
    with open(SUMMARY_PATH, 'w') as f:
        for label, value in results.items():
            f.write('{:66s} {}\n'.format(label, _describe(value)))
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

    print('\n=== IO baseline comparison ===')
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
        print('\nNumeric/text mismatches:')
        for label, msg in mismatches:
            print('  * {}: {}'.format(label, msg))
    if not (only_old or only_new or exc_changes or mismatches):
        print('\nALL MATCH -- no differences within tolerance.')
    else:
        sys.exit(1)


if __name__ == '__main__':
    _SCRATCH_RE = SCRATCH
    mode = sys.argv[1] if len(sys.argv) > 1 else 'save'
    if mode == 'save':
        save()
    elif mode == 'compare':
        compare()
    else:
        sys.exit('usage: python tests/io_baseline.py [save|compare]')
