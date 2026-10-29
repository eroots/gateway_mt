"""Plot raw data, inversion data and / or a response to a multi-page PDF.

Usage::

    plot_dataset [-l LIST] [-d DATA] [-r RESPONSE] [-o OUT.pdf]
                 [--plot rho-phase | impedance] [--sites SITE ...]

Any combination works: raw data alone (the EDI / j-format files in a list),
an inversion data file before inverting, data against a response, or all
three. Four sites per page; each site shows apparent resistivity and phase
(default) or the impedance tensor (``--plot impedance``), with the tipper
underneath when present. Inversion data is drawn as points with error bars,
the response as lines, and raw data as small faint points (or as the main
points when it is the only input). Any data format pyMT reads is accepted;
the format is taken from the file. Installed as the ``plot_dataset`` script.

The page layout is shared with ``autogen_report``.
"""
import argparse
import os
import sys

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

import pyMT.data_structures as WSDS
import pyMT.utils as utils

SITES_PER_PAGE = 4
COLOURS = {'xy': 'b', 'yx': 'r', 'xx': 'b', 'yy': 'r', 'zx': 'b', 'zy': 'r'}
ALPHA = {'xy': 0.9, 'yx': 0.9, 'xx': 0.4, 'yy': 0.4, 'zx': 0.9, 'zy': 0.9}
# Which error attribute each kind of site carries (responses have none)
ERRORS = {'raw_data': 'errors', 'data': 'used_error', 'response': None}


# --------------------------------------------------------------------------- #
# Drawing helpers
# --------------------------------------------------------------------------- #
def _flagged(site, kind, key):
    """True where ``key`` is flagged as removed (REMOVE_FLAG error)."""
    if ERRORS[kind] is None or key not in getattr(site, ERRORS[kind]):
        return np.zeros(len(site.periods), dtype=bool)
    return np.asarray(getattr(site, ERRORS[kind])[key]) == site.REMOVE_FLAG


def _draw(ax, kind, primary, periods, values, errors, flagged, colour, alpha, filled=True, dashed=False,
          label=None):
    """Points with error bars (primary points), faint points (secondary raw data) or a line (response)."""
    values = np.where(flagged, np.nan, np.real(values))
    if kind == 'response':
        ax.plot(periods, values, linestyle='--' if dashed else '-', color=colour, alpha=alpha, label=label)
    elif kind == primary:
        errors = None if errors is None else np.where(flagged, np.nan, errors)
        ax.errorbar(periods, values, yerr=errors, linestyle='', marker='o', markersize=4,
                    color=colour, mec='k', mew=0.5, mfc=colour if filled else 'w', alpha=alpha, label=label)
    else:
        ax.plot(periods, values, linestyle='', marker='o', markersize=2.5, color=colour,
                mfc=colour if filled else 'none', alpha=0.35 * alpha)


def _rho_phase_panels(axes, sites, primary, has_diag):
    rho_ax, pha_ax = axes
    off_diagonal = []
    for kind, site in sites.items():
        err = ERRORS[kind] or 'none'
        for comp in ('xy', 'yx') + (('xx', 'yy') if has_diag else ()):
            key = 'Z{}R'.format(comp.upper())
            if key not in site.components:
                continue
            flagged = _flagged(site, kind, key)
            rho, rho_err, _ = utils.compute_rho(site, calc_comp=comp, errtype=err)
            _draw(rho_ax, kind, primary, site.periods, rho, rho_err, flagged, COLOURS[comp], ALPHA[comp],
                  label=comp.upper() if kind == primary else None)
            if comp in ('xy', 'yx'):
                off_diagonal += list(np.real(rho)[~flagged])
                pha, pha_err = utils.compute_phase(site, calc_comp=comp, errtype=ERRORS[kind], wrap=1)
                _draw(pha_ax, kind, primary, site.periods, pha, pha_err, flagged, COLOURS[comp], ALPHA[comp])
    rho_ax.set_yscale('log')
    # Limits from the off-diagonal curves; small diagonal values are clipped
    off_diagonal = np.array([x for x in off_diagonal if np.isfinite(x) and x > 0])
    if off_diagonal.size:
        rho_ax.set_ylim(10 ** (np.floor(np.log10(off_diagonal.min())) - 1),
                        10 ** (np.ceil(np.log10(off_diagonal.max())) + 0.5))
    rho_ax.set_ylabel(r'$\rho_a$ ($\Omega$m)')
    pha_ax.set_ylim(0, 120)
    pha_ax.set_ylabel('Phase (deg)')


def _impedance_panels(axes, sites, primary, has_diag):
    groups = [('xy', 'yx')] + ([('xx', 'yy')] if has_diag else [])
    for ax, comps in zip(axes, groups):
        for kind, site in sites.items():
            root_t = np.sqrt(site.periods)
            for comp in comps:
                for part, filled, dashed in (('R', True, False), ('I', False, True)):
                    key = 'Z{}{}'.format(comp.upper(), part)
                    if key not in site.components:
                        continue
                    errors = getattr(site, ERRORS[kind])[key] * root_t if ERRORS[kind] else None
                    _draw(ax, kind, primary, site.periods, site.data[key] * root_t, errors,
                          _flagged(site, kind, key), COLOURS[comp], ALPHA[comp], filled=filled, dashed=dashed,
                          label='{} {}'.format(comp.upper(), 'Re' if part == 'R' else 'Im')
                          if kind == primary else None)
        ax.set_ylabel(r'Z$\cdot\sqrt{T}$' + ('\n(off-diagonal)' if comps[0] == 'xy' else '\n(diagonal)'))


def _tipper_panel(ax, sites, primary):
    for kind, site in sites.items():
        for comp in ('zx', 'zy'):
            for part, filled, dashed in (('R', True, False), ('I', False, True)):
                key = 'T{}{}'.format(comp.upper(), part)
                if key not in site.components:
                    continue
                errors = getattr(site, ERRORS[kind])[key] if ERRORS[kind] else None
                _draw(ax, kind, primary, site.periods, site.data[key], errors, _flagged(site, kind, key),
                      COLOURS[comp], ALPHA[comp], filled=filled, dashed=dashed,
                      label='{} {}'.format(comp.upper(), 'Re' if part == 'R' else 'Im')
                      if kind == primary else None)
    ax.axhline(0, color='0.6', linewidth=0.5)
    ax.set_ylabel('Tipper')


def _has(site, keys, kind='data'):
    return site is not None and all(k in site.components and not _flagged(site, kind, k).all() for k in keys)


def plot_page(fig, site_groups, plot_type='rho-phase', site_rms=None):
    """Draw up to four sites (2 x 2) on ``fig``.

    ``site_groups`` is a list of ``{kind: Site}`` dicts, one per site, with
    ``kind`` in ``'raw_data'``, ``'data'``, ``'response'``. The first kind
    present out of data / raw data / response is drawn as the main points.
    """
    outer = GridSpec(2, 2, figure=fig, hspace=0.25, wspace=0.25)
    kinds_seen = set()
    for kk, sites in enumerate(site_groups):
        primary = next(k for k in ('data', 'raw_data', 'response') if k in sites)
        ref = sites[primary]
        kinds_seen.update(sites)
        has_diag = _has(ref, ('ZXXR', 'ZYYR'), primary)
        has_tipper = _has(ref, ('TZXR', 'TZYR'), primary)
        if not _has(ref, ('ZXYR', 'ZYXR'), primary):
            main_heights = []   # tipper-only data
        elif plot_type == 'rho-phase':
            main_heights = [3, 2]
        else:
            main_heights = [2, 2] if has_diag else [3]
        heights = main_heights + ([2] if has_tipper else [])
        if not heights:
            continue
        inner = outer[kk].subgridspec(len(heights), 1, height_ratios=heights, hspace=0.08)
        axes = [fig.add_subplot(inner[0])]
        axes += [fig.add_subplot(inner[ii], sharex=axes[0]) for ii in range(1, len(heights))]
        if main_heights and plot_type == 'rho-phase':
            _rho_phase_panels(axes[:2], sites, primary, has_diag)
        elif main_heights:
            _impedance_panels(axes[:len(main_heights)], sites, primary, has_diag)
        if has_tipper:
            _tipper_panel(axes[-1], sites, primary)
        for ax in axes:
            ax.set_xscale('log')
            ax.grid(True, which='major', alpha=0.3)
            if ax is not axes[-1]:
                plt.setp(ax.get_xticklabels(), visible=False)
        axes[-1].set_xlabel('Period (s)')
        title = ref.name
        if site_rms and ref.name in site_rms:
            title += '  (RMS {:.2f})'.format(site_rms[ref.name])
        axes[0].set_title(title)
        if kk == 0:
            axes[0].legend(loc='best', fontsize=7, ncol=2)
            if has_tipper and len(axes) > 1:
                axes[-1].legend(loc='best', fontsize=7, ncol=2)
    handles = []
    if 'data' in kinds_seen:
        handles.append(Line2D([], [], marker='o', linestyle='', color='0.4', mec='k', label='data'))
    if 'raw_data' in kinds_seen:
        primary_raw = 'data' not in kinds_seen
        handles.append(Line2D([], [], marker='o', linestyle='', color='0.4', mec='k' if primary_raw else '0.4',
                              markersize=6 if primary_raw else 3, alpha=1 if primary_raw else 0.4,
                              label='raw data'))
    if 'response' in kinds_seen:
        handles.append(Line2D([], [], color='0.4', label='response'))
    fig.legend(handles=handles, loc='upper center', ncol=len(handles), frameon=False)


def write_pdf(dataset, outfile, kinds, plot_type='rho-phase', site_names=None, site_rms=None):
    """Write ``kinds`` (any of 'raw_data', 'data', 'response') of ``dataset`` to a PDF.

    Returns the number of pages written.
    """
    if site_names is None:
        site_names = next(getattr(dataset, k).site_names for k in ('data', 'response', 'raw_data') if k in kinds)
    pages = 0
    with PdfPages(outfile) as pdf:
        for ii in range(0, len(site_names), SITES_PER_PAGE):
            fig = plt.figure(figsize=(12, 16))
            groups = []
            for name in site_names[ii:ii + SITES_PER_PAGE]:
                groups.append({k: getattr(dataset, k).sites[name] for k in ('raw_data', 'data', 'response')
                               if k in kinds and name in getattr(dataset, k).sites})
            plot_page(fig, [g for g in groups if g], plot_type, site_rms)
            pdf.savefig(fig)
            plt.close(fig)
            pages += 1
    return pages


# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='plot_dataset',
        description='Plot raw data, inversion data and / or a response to a PDF, four sites per page. '
                    'Give any combination of -l, -d and -r; data formats are detected from the files.')
    parser.add_argument('-l', '--list', help='site list -- plots the raw EDI / j-format data')
    parser.add_argument('-d', '--data', help='inversion data file')
    parser.add_argument('-r', '--response', help='inversion response (predicted data) file')
    parser.add_argument('-o', '--output', help='output PDF (default: <first input name>_plots.pdf)')
    parser.add_argument('--plot', choices=('rho-phase', 'impedance'), default='rho-phase',
                        help='plot apparent resistivity and phase (default) or the impedance tensor; '
                             'the tipper is added in both cases when present')
    parser.add_argument('--sites', nargs='+', metavar='SITE', help='only plot these sites')
    parser.add_argument('-f', '--force', action='store_true', help='overwrite the output file')
    args = parser.parse_args(argv)
    if not (args.list or args.data or args.response):
        parser.error('give at least one of -l (raw data), -d (data) or -r (response)')
    return args


def main(argv=None):
    args = parse_args(argv)
    for label, path in (('list', args.list), ('data', args.data), ('response', args.response)):
        if path and not os.path.isfile(path):
            print('plot_dataset: {} file {} not found'.format(label, path))
            return 1
    first = args.data or args.response or args.list
    outfile = args.output or '{}_plots.pdf'.format(os.path.splitext(os.path.basename(first))[0])
    if os.path.exists(outfile) and not args.force:
        print('plot_dataset: {} already exists (use --force to overwrite)'.format(outfile))
        return 1
    dataset = WSDS.Dataset(listfile=args.list or '', datafile=args.data or '', responsefile=args.response or '')
    # (A list-only Dataset also fills .data from the raw data; plot only what was asked for.)
    kinds = [k for k, given in (('raw_data', args.list), ('data', args.data), ('response', args.response)) if given]
    site_rms = None
    if args.data and args.response and dataset.rms:
        site_rms = {site: info['Total'] for site, info in dataset.rms['Station'].items()}
    site_names = None
    if args.sites:
        available = set().union(*[getattr(dataset, k).site_names for k in kinds])
        missing = [s for s in args.sites if s not in available]
        if missing:
            print('plot_dataset: sites not found: {}'.format(', '.join(missing)))
            return 1
        site_names = args.sites
    pages = write_pdf(dataset, outfile, kinds, args.plot, site_names, site_rms)
    print('Wrote {} ({} page{})'.format(outfile, pages, '' if pages == 1 else 's'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
