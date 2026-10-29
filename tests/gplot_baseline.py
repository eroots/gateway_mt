"""Before/after verification baseline for ``pyMT.gplot``.

Part of the ``cleanup/gplot`` refactor workflow.  ``gplot`` produces
matplotlib artists rather than numbers, so this harness drives every
public ``DataPlotManager`` / ``MapView`` method against the sample data
and freezes the *artist state* it produces (line xy-data, collection
offsets, quiver vectors, patch vertices, image arrays, axis
labels/limits, ``format_coord`` strings) plus the handful of pure-numeric
helpers.  Everything is stored for ``numpy.allclose`` comparison after
refactoring.

Usage
-----
    python tests/gplot_baseline.py save      # snapshot current behaviour
    python tests/gplot_baseline.py compare   # re-run and diff

Outputs land in ``tests/baseline_gplot/``:

* ``gplot_baseline.pkl``         -- pickled ``{label: value}`` snapshot
* ``gplot_baseline_summary.txt`` -- human-readable shapes / sample values
"""

import os
import sys
import pickle
import zlib
import traceback
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.colors
from matplotlib.figure import Figure

import pyMT.gplot as gplot
from pyMT.data_structures import Dataset, Data, RawData, Response, Model

warnings.filterwarnings('ignore')

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.join(THIS_DIR, 'data', 'synthetic')

EDI_LIST = os.path.join(SYN, 'edis', 'synthetic.lst')
MODEM_FILE = os.path.join(SYN, 'synthetic_modem.dat')
RESP_FILE = os.path.join(SYN, 'synthetic_modem_resp.dat')
MODEL_FILE = os.path.join(SYN, 'synthetic.model')

OUT_DIR = os.path.join(THIS_DIR, 'baseline_gplot')
PICKLE_PATH = os.path.join(OUT_DIR, 'gplot_baseline.pkl')
SUMMARY_PATH = os.path.join(OUT_DIR, 'gplot_baseline_summary.txt')

# A fixed handful of sites present in the ModEM data file.
SAMPLE_SITES = ['SYN001', 'SYN015', 'SYN029', 'SYN030']
PERIOD_IDX = 10


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
        if o_seq.ndim == 0:  # scalar non-numeric object (None, dict, ...): no recursion
            return '' if o_seq.item() == n_seq.item() else '{!r} -> {!r}'.format(o_seq.item(), n_seq.item())
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
    loose = ('/loc', 'get_locations', 'induction', 'phase_tensor', 'phase_bar',
             'pseudosection', 'interpolate', 'plan_view', '_slice', 'format_coord')
    return (1e-6, 1e-6) if any(k in label for k in loose) else (1e-9, 1e-12)


# --------------------------------------------------------------------------- #
# artist introspection
# --------------------------------------------------------------------------- #
def _snapshot_axes(results, prefix, ax):
    """Freeze every drawable on ``ax`` into comparable arrays."""
    results['{}/title'.format(prefix)] = np.asarray(str(ax.get_title()))
    results['{}/xlabel'.format(prefix)] = np.asarray(str(ax.get_xlabel()))
    results['{}/ylabel'.format(prefix)] = np.asarray(str(ax.get_ylabel()))
    results['{}/xlim'.format(prefix)] = np.asarray(ax.get_xlim(), dtype=float)
    results['{}/ylim'.format(prefix)] = np.asarray(ax.get_ylim(), dtype=float)
    results['{}/yscale'.format(prefix)] = np.asarray(str(ax.get_yscale()))
    results['{}/aspect'.format(prefix)] = np.asarray(str(ax.get_aspect()))

    for i, ln in enumerate(ax.get_lines()):
        xy = np.asarray(ln.get_xydata(), dtype=float)
        results['{}/line{}/xy'.format(prefix, i)] = xy
        results['{}/line{}/color'.format(prefix, i)] = np.asarray(
            matplotlib.colors.to_rgba(ln.get_color()), dtype=float)

    for i, coll in enumerate(ax.collections):
        # Build one combined [offsets | array | UV] block and sort its rows so
        # that set-ordering nondeterminism in site lists (MapView.generic_sites
        # is a set difference) does not show up as a spurious diff.
        cols, tags = [], []
        try:
            off = np.asarray(coll.get_offsets(), dtype=float)
            if off.ndim == 2 and off.shape[0] > 1:
                cols.append(off)
                tags.append(('offsets', 2))
        except Exception:  # noqa: BLE001
            off = None
        arr = coll.get_array()
        if arr is not None:
            arr = np.asarray(arr, dtype=float).reshape(-1)
        if hasattr(coll, 'get_UVC'):
            u, v = coll.get_UVC()
            uv = np.column_stack([np.asarray(u, dtype=float), np.asarray(v, dtype=float)])
        else:
            uv = None
        if cols and arr is not None and arr.shape[0] == cols[0].shape[0]:
            cols.append(arr.reshape(-1, 1)); tags.append(('array', 1))
        if cols and uv is not None and uv.shape[0] == cols[0].shape[0]:
            cols.append(uv); tags.append(('UV', 2))
        if cols:
            block = np.hstack(cols)
            order = np.lexsort(block.T[::-1])
            block = block[order]
            c = 0
            for name, w in tags:
                results['{}/coll{}/{}'.format(prefix, i, name)] = block[:, c:c + w]
                c += w
        else:
            if off is not None:
                results['{}/coll{}/offsets'.format(prefix, i)] = np.atleast_2d(off)
            if arr is not None:
                results['{}/coll{}/array'.format(prefix, i)] = arr
            if uv is not None:
                results['{}/coll{}/UV'.format(prefix, i)] = uv
        try:
            results['{}/coll{}/clim'.format(prefix, i)] = np.asarray(coll.get_clim(), dtype=float)
        except Exception:  # noqa: BLE001
            pass
        if hasattr(coll, 'get_segments'):
            segs = coll.get_segments()
            if segs:
                seg = np.asarray(np.concatenate(segs), dtype=float)
                results['{}/coll{}/segments'.format(prefix, i)] = seg[np.lexsort(seg.T[::-1])]

    for i, patch in enumerate(ax.patches):
        try:
            xy = patch.get_xy()
            results['{}/patch{}/xy'.format(prefix, i)] = np.atleast_2d(np.asarray(xy, dtype=float))
        except Exception:  # noqa: BLE001
            pass
        results['{}/patch{}/facecolor'.format(prefix, i)] = np.asarray(
            patch.get_facecolor(), dtype=float)

    for i, im in enumerate(ax.images):
        results['{}/image{}/array'.format(prefix, i)] = np.asarray(im.get_array(), dtype=float)
        results['{}/image{}/clim'.format(prefix, i)] = np.asarray(im.get_clim(), dtype=float)
        results['{}/image{}/extent'.format(prefix, i)] = np.asarray(im.get_extent(), dtype=float)

    texts = [(str(t.get_text()), float(t.get_position()[0]), float(t.get_position()[1]))
             for t in ax.texts]
    if texts:
        results['{}/texts'.format(prefix)] = np.asarray(sorted(texts), dtype=object)


# --------------------------------------------------------------------------- #
# fixtures
# --------------------------------------------------------------------------- #
def _load():
    np.random.seed(0)
    dataset = Dataset(listfile=EDI_LIST, datafile=MODEM_FILE, responsefile=RESP_FILE)
    data = Data(datafile=MODEM_FILE, listfile=EDI_LIST)
    raw = RawData(listfile=EDI_LIST)
    resp = Response(datafile=RESP_FILE, listfile=EDI_LIST)
    raw.locations = raw.get_locs(mode='utm')
    data.locations = raw.locations
    resp.locations = raw.locations
    model, _ = Model(modelfile=MODEL_FILE), None
    return dataset, data, raw, resp, model


# --------------------------------------------------------------------------- #
# DataPlotManager exercises
# --------------------------------------------------------------------------- #
DPM_COMPONENT_SETS = {
    'offdiag': ['ZXYR', 'ZXYI', 'ZYXR', 'ZYXI'],
    'diag': ['ZXXR', 'ZXXI', 'ZYYR', 'ZYYI'],
    'rho': ['RHOXY', 'RHOYX'],
    'pha': ['PHAXY', 'PHAYX'],
    'tip': ['TZXR', 'TZYR'],
    'pt': ['PTXY', 'PTYX'],
    'bost': ['BOSTXY'],
    'beta': ['BETA'],
}


def exercise_dataplotmanager(results):
    dataset, *_ = _load()
    sites = dataset.get_sites(site_names=SAMPLE_SITES, dTypes='all')

    # property / helper coverage that doesn't need a plot
    dpm0 = gplot.DataPlotManager(fig=Figure())
    dpm0.sites = sites
    dpm0.components = ['ZXYR']
    results['dpm/site_names'] = np.asarray(dpm0.site_names, dtype=object)
    results['dpm/num_sites'] = np.asarray(dpm0.num_sites)
    for n in (1, 2, 3, 5, 6, 9, 12):
        dpm0.sites = {k: v[:n] for k, v in sites.items()}
        results['dpm/gen_tiling/{}'.format(n)] = np.asarray(dpm0.gen_tiling())
    dpm0.sites = sites
    for comp in ('ZXYR', 'RHOXY', 'PHAXY', 'TZXR', 'PTXY', 'BOSTXY', 'BETA'):
        dpm0.components = [comp]
        results['dpm/units/{}'.format(comp)] = np.asarray(str(dpm0.units))
        results['dpm/log10_yscale/{}'.format(comp)] = np.asarray(bool(dpm0.log10_yscale))
    for i, comp in enumerate(['ZXYR', 'ZXYI', 'RHOXY', 'PHAYX', 'tzx', 'ssq', 'det']):
        results['dpm/designated_colour/{}'.format(comp)] = np.asarray(
            matplotlib.colors.to_rgba(dpm0.get_designated_colour(comp, i)), dtype=float)

    for cname, comps in DPM_COMPONENT_SETS.items():
        # the scale knob only bites on raw impedance; sweep it fully there and
        # just spot-check one alt scale for the derived quantities.
        scales = (('sqrt(periods)', 'periods', 'log10', 'none')
                  if cname in ('offdiag', 'diag') else ('sqrt(periods)', 'log10'))
        for scale in scales:
            for errors in ('mapped', 'none'):
                tag = '{}/{}/{}'.format(cname, scale.replace('(', '').replace(')', ''), errors)

                def _do(comps=comps, scale=scale, errors=errors):
                    fig = Figure()
                    dpm = gplot.DataPlotManager(fig=fig)
                    dpm.components = list(comps)
                    dpm.scale = scale
                    dpm.errors = errors
                    dpm.show_outliers = True
                    dpm.plot_flagged_data = False
                    dpm.sites = sites
                    dpm.plot_data()
                    return dpm

                dpm = _store(results, 'dpm/plot/{}/_call'.format(tag), _do)
                if isinstance(dpm, str):
                    continue
                for ai, ax in enumerate(dpm.axes):
                    _snapshot_axes(results, 'dpm/plot/{}/ax{}'.format(tag, ai), ax)
                results.pop('dpm/plot/{}/_call'.format(tag), None)

    # plot_site direct return values (ax, max, min, artist)
    for cname, comps in DPM_COMPONENT_SETS.items():
        def _do(comps=comps):
            fig = Figure()
            ax = fig.add_subplot(111)
            dpm = gplot.DataPlotManager(fig=fig)
            dpm.axes = [ax]          # plot_site -> set_axis_scales touches self.axes
            dpm.components = list(comps)
            dpm.errors = 'mapped'
            # one site only: set_axis_scales assumes len(axes) == len(site_names)
            dpm.sites = {k: v[:1] for k, v in sites.items()}
            site = dpm.sites['data'][0]
            _, mx, mn, _ = dpm.plot_site(site, dType='data', ax=ax)
            return np.asarray([mx, mn], dtype=float)
        _store(results, 'dpm/plot_site/{}/maxmin'.format(cname), _do)

    # redraw_single_axis.  Note: it indexes self.sites['1d'], a key that
    # Dataset.get_sites() (and DataPlotManager.__init__) do not create -> KeyError
    # today.  Freeze that, and also run a working variant with '1d' injected.
    def _do_redraw_bare():
        fig = Figure()
        dpm = gplot.DataPlotManager(fig=fig)
        dpm.components = ['ZXYR', 'ZYXR']
        dpm.sites = sites
        dpm.plot_data()
        dpm.redraw_single_axis(axnum=2)
        return np.asarray('ok')
    _store(results, 'dpm/redraw_single/no_1d_key', _do_redraw_bare)

    def _do_redraw():
        fig = Figure()
        dpm = gplot.DataPlotManager(fig=fig)
        dpm.components = ['ZXYR', 'ZYXR']
        dpm.sites = dict(sites, **{'1d': []})
        dpm.plot_data()
        dpm.redraw_single_axis(axnum=2)
        return dpm
    dpm = _store(results, 'dpm/redraw_single/_call', _do_redraw)
    if not isinstance(dpm, str):
        _snapshot_axes(results, 'dpm/redraw_single/ax2', dpm.axes[2])
        results.pop('dpm/redraw_single/_call', None)

    def _do_link():
        fig = Figure()
        dpm = gplot.DataPlotManager(fig=fig)
        dpm.components = ['ZXYR']
        dpm.sites = sites
        dpm.plot_data()
        dpm.link_axes(y_bounds=[-1.0, 1.0])
        return np.asarray([ax.get_ylim() for ax in dpm.axes], dtype=float)
    _store(results, 'dpm/link_axes/ylims', _do_link)


# --------------------------------------------------------------------------- #
# MapView exercises
# --------------------------------------------------------------------------- #
def _fresh_mapview(data, raw, resp, colourmap='turbo', **attrs):
    fig = Figure()
    ax = fig.add_subplot(111)
    mv = gplot.MapView(fig=fig)
    mv.window['figure'] = fig
    mv.window['axes'] = [ax]
    mv.site_data['data'] = data
    mv.site_data['raw_data'] = raw
    mv.site_data['response'] = resp
    mv.site_names = list(data.site_names)
    mv.dataset = None
    mv.colourmap = colourmap
    mv.interpolant = 'linear'
    # ``plot_plan_view`` reads ``self.use_colorbar`` (missing 'u'); shim it so the
    # model plots run today.  Post-fix the harness keeps colourbars off, so the
    # snapshot is unaffected either way.
    mv.use_colorbar = False
    mv.use_colourbar = False
    mv.site_locations['all'] = data.locations
    for k, v in attrs.items():
        setattr(mv, k, v)
    return mv, ax


def exercise_mapview_helpers(results):
    dataset, data, raw, resp, model = _load()
    mv, ax = _fresh_mapview(data, raw, resp)

    # pure-numeric geometry helpers
    phi = np.array([[1.2, 0.3], [-0.2, 0.8]])
    results['mv/generate_ellipse'] = np.asarray(mv.generate_ellipse(phi), dtype=float)
    results['mv/generate_rectangle'] = np.asarray(
        mv.generate_rectangle(0.4, 0.2, 0.5), dtype=float)
    ax_x, ax_y, ret = mv.generate_arc(phi, 0.1, 1.2)
    results['mv/generate_arc'] = np.asarray([ax_x, ax_y], dtype=float)
    results['mv/generate_arc/ret'] = np.asarray(bool(ret))
    results['mv/resize_ellipse_aspect'] = np.asarray(
        mv.resize_ellipse_aspect(0.6, 1.2), dtype=float)
    results['mv/rescale_ellipse'] = np.asarray(
        mv.rescale_ellipse(np.arange(5.0), np.arange(5.0), 2.5), dtype=float)

    for fp in ('beta', 'absbeta', 'alpha', 'azimuth', 'delta', 'diff_delta',
               'Lambda', 'phi_split_pt', 'phi_split', 'dimensionality',
               'diff_beta', 'phi_max', 'phi_min', 'det_phi'):
        for pt in ('phi', 'phi_a', 'ua', 'va'):
            try:
                lo, hi = mv.pt_fill_limits(fp, pt)
                results['mv/pt_fill_limits/{}/{}'.format(fp, pt)] = np.asarray([lo, hi], dtype=float)
            except Exception as exc:  # noqa: BLE001
                results['mv/pt_fill_limits/{}/{}'.format(fp, pt)] = 'EXCEPTION: {}: {}'.format(
                    type(exc).__name__, exc)

    for p in ('alpha', 'beta', 'phi_1', 'phi_2', 'phi_3', 'det_phi',
              'azimuth', 'phi_max', 'phi_min', 'phi_split_pt', 'nonsense'):
        results['mv/get_label/{}'.format(p)] = np.asarray(str(mv.get_label(p)))

    # get_tensor_params for a real site
    site = data.sites[SAMPLE_SITES[0]]
    for fp in ('beta', 'phi_max', 'azimuth'):
        def _do_phi(fp=fp):
            _, phi_, _, _, _, _ = mv.get_tensor_params('phi', fp, site=site, period_idx=PERIOD_IDX)
            return np.asarray(phi_, dtype=float).ravel()

        def _do_scalars(fp=fp):
            _, _, mx, mn, azi, fill = mv.get_tensor_params('phi', fp, site=site, period_idx=PERIOD_IDX)
            return np.asarray([mx, mn, azi, fill], dtype=float)
        _store(results, 'mv/get_tensor_params/{}/phi'.format(fp), _do_phi)
        _store(results, 'mv/get_tensor_params/{}/scalars'.format(fp), _do_scalars)

    # get_locations (local coordinate system)
    def _do_locs():
        return np.asarray(mv.get_locations(sites=mv.site_names), dtype=float)
    _store(results, 'mv/get_locations/local', _do_locs)

    # interpolate
    def _do_interp():
        pts = np.column_stack([data.locations[:, 1], data.locations[:, 0],
                               np.zeros(data.locations.shape[0])])
        vals = np.linspace(1.0, 5.0, pts.shape[0])
        gv, gx, gy, xf, yf = mv.interpolate(pts, vals, 25)
        return np.asarray(gv, dtype=float)
    _store(results, 'mv/interpolate/linear', _do_interp)


def exercise_mapview_plots(results):
    dataset, data, raw, resp, model = _load()

    # plot_locations (no RMS)
    def _do_loc():
        mv, ax = _fresh_mapview(data, raw, resp)
        mv.plot_rms = False
        mv.set_locations()
        mv.plot_locations()
        return mv, ax
    got = _store(results, 'mv/plot_locations/_call', _do_loc)
    if not isinstance(got, str):
        _snapshot_axes(results, 'mv/plot_locations', got[1])
        results.pop('mv/plot_locations/_call', None)

    # plot_annotate
    def _do_ann():
        mv, ax = _fresh_mapview(data, raw, resp)
        mv.active_sites = SAMPLE_SITES
        mv.annotate_sites = 'active'
        mv.plot_annotate()
        return mv, ax
    got = _store(results, 'mv/plot_annotate/_call', _do_ann)
    if not isinstance(got, str):
        _snapshot_axes(results, 'mv/plot_annotate', got[1])
        results.pop('mv/plot_annotate/_call', None)

    # induction arrows
    for dt in ('data', 'response'):
        for conv in ('parkinson', 'weise'):
            def _do(dt=dt, conv=conv):
                mv, ax = _fresh_mapview(data, raw, resp)
                mv.induction_scale = 5
                mv.induction_cutoff = 2
                mv.plot_induction_arrows(data_type=[dt], normalize=False,
                                         period_idx=PERIOD_IDX, arrow_type=['R', 'I'],
                                         arrow_convention=conv)
                return mv, ax
            got = _store(results, 'mv/induction/{}/{}/_call'.format(dt, conv), _do)
            if not isinstance(got, str):
                _snapshot_axes(results, 'mv/induction/{}/{}'.format(dt, conv), got[1])
                results.pop('mv/induction/{}/{}/_call'.format(dt, conv), None)

    # phase tensor ellipses
    for fp in ('beta', 'phi_max', 'phi_min', 'azimuth', 'phi_split_pt'):
        def _do(fp=fp):
            mv, ax = _fresh_mapview(data, raw, resp)
            mv.min_pt_ratio = 0.2
            mv.pt_scale = 1.0
            mv.plot_phase_tensor(data_type=['data'], normalize=True,
                                 fill_param=fp, period_idx=PERIOD_IDX, cb_slot=0)
            return mv, ax
        got = _store(results, 'mv/phase_tensor/{}/_call'.format(fp), _do)
        if not isinstance(got, str):
            _snapshot_axes(results, 'mv/phase_tensor/{}'.format(fp), got[1])
            results.pop('mv/phase_tensor/{}/_call'.format(fp), None)

    # phase tensor difference (data - response)
    def _do_ptdiff():
        mv, ax = _fresh_mapview(data, raw, resp)
        mv.pt_difference_definition = 'Booker'
        mv.plot_phase_tensor(data_type=['data', 'response'], normalize=True,
                             fill_param='beta', period_idx=PERIOD_IDX, cb_slot=0)
        return mv, ax
    got = _store(results, 'mv/phase_tensor_diff/_call', _do_ptdiff)
    if not isinstance(got, str):
        _snapshot_axes(results, 'mv/phase_tensor_diff', got[1])
        results.pop('mv/phase_tensor_diff/_call', None)

    # phase bar
    def _do_bar():
        mv, ax = _fresh_mapview(data, raw, resp)
        mv.plot_phase_bar(data_type=['data'], normalize=True,
                          fill_param='beta', period_idx=PERIOD_IDX, cb_slot=0)
        return mv, ax
    got = _store(results, 'mv/phase_bar/_call', _do_bar)
    if not isinstance(got, str):
        _snapshot_axes(results, 'mv/phase_bar', got[1])
        results.pop('mv/phase_bar/_call', None)

    # plan pseudosection
    for fp in ('rhoxy', 'phaxy', 'tip', 'beta', 'azimuth', 'bostxy', 'phi_min'):
        def _do(fp=fp):
            mv, ax = _fresh_mapview(data, raw, resp)
            mv.plan_pseudosection(data_type=['data'], fill_param=fp,
                                  n_interp=30, period_idx=PERIOD_IDX, cb_slot=0)
            return mv, ax
        got = _store(results, 'mv/plan_pseudo/{}/_call'.format(fp), _do)
        if not isinstance(got, str):
            _snapshot_axes(results, 'mv/plan_pseudo/{}'.format(fp), got[1])
            results.pop('mv/plan_pseudo/{}/_call'.format(fp), None)

    # plan pseudosection difference
    def _do_pdiff():
        mv, ax = _fresh_mapview(data, raw, resp)
        mv.plan_pseudosection(data_type=['data', 'response'], fill_param='rhoxy',
                              n_interp=30, period_idx=PERIOD_IDX, cb_slot=0)
        return mv, ax
    got = _store(results, 'mv/plan_pseudo_diff/_call', _do_pdiff)
    if not isinstance(got, str):
        _snapshot_axes(results, 'mv/plan_pseudo_diff', got[1])
        results.pop('mv/plan_pseudo_diff/_call', None)

    # tensor pseudosections along a profile
    prof = SAMPLE_SITES
    for meth in ('tensor_ellipse_pseudosection', 'tensor_bar_pseudosection'):
        def _do(meth=meth):
            mv, ax = _fresh_mapview(data, raw, resp)
            periods = (data.periods[2], data.periods[-3], 0)
            getattr(mv, meth)(sites=list(prof), data_type='data',
                              fill_param='beta', periods=periods,
                              x_axis='linear')
            return mv, ax
        got = _store(results, 'mv/{}/_call'.format(meth), _do)
        if not isinstance(got, str):
            _snapshot_axes(results, 'mv/{}'.format(meth), got[1])
            results.pop('mv/{}/_call'.format(meth), None)


def exercise_mapview_model(results):
    dataset, data, raw, resp, model = _load()

    for meth, kw in (('plot_plan_view', dict(z_slice=4)),
                     ('plot_x_slice', dict(x_slice=10)),
                     ('plot_y_slice', dict(y_slice=10))):
        def _do(meth=meth, kw=kw):
            mv, ax = _fresh_mapview(data, raw, resp)
            mv.model = model
            getattr(mv, meth)(ax=ax, rho_axis='rho_x', **kw)
            return mv, ax
        got = _store(results, 'mv/model/{}/_call'.format(meth), _do)
        if not isinstance(got, str):
            _snapshot_axes(results, 'mv/model/{}'.format(meth), got[1])
            # format_coord readout at a fixed cursor point
            fc = got[0].window['axes'][0].format_coord
            results['mv/model/{}/format_coord'.format(meth)] = np.asarray(
                str(fc(float(np.mean(got[1].get_xlim())),
                       float(np.mean(got[1].get_ylim())))))
            results.pop('mv/model/{}/_call'.format(meth), None)


# --------------------------------------------------------------------------- #
# collect / save / compare
# --------------------------------------------------------------------------- #
def collect():
    results = {}
    exercise_dataplotmanager(results)
    exercise_mapview_helpers(results)
    exercise_mapview_plots(results)
    exercise_mapview_model(results)
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

    print('\n=== gplot baseline comparison ===')
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
    mode = sys.argv[1] if len(sys.argv) > 1 else 'save'
    if mode == 'save':
        save()
    elif mode == 'compare':
        compare()
    else:
        sys.exit('usage: python tests/gplot_baseline.py [save|compare]')
