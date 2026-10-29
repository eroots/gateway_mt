"""Write a misfit report for an inversion: a text summary plus a PDF of misfit curves.

Usage::

    autogen_report -d DATA -r RESPONSE [-l LIST] [-m MODEL] [-o REPORT_DIR]
                   [--plot rho-phase | impedance]

Writes ``REPORT_DIR/summary.txt`` (input files, inversion type, sites and
periods, overall / per-component RMS, worst-fitting sites, the period table,
and blank fields to fill in by hand) and ``REPORT_DIR/misfit_curves.pdf``:
data (points with error bars) against response (lines), four sites per page.
Each site shows apparent resistivity and phase (default) or the impedance
tensor (``--plot impedance``), with the tipper underneath when present.
Works with any data format pyMT reads (ModEM, WSINV3DMT, GoFEM, ...); the
format is taken from the file. Installed as the ``autogen_report`` script.
"""
import argparse
import datetime
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
COLOURS = {'xy': 'b', 'yx': 'r', 'xx': 'b', 'yy': 'r', 'tzx': 'b', 'tzy': 'r'}
ALPHA = {'xy': 0.9, 'yx': 0.9, 'xx': 0.4, 'yy': 0.4, 'tzx': 0.9, 'tzy': 0.9}


def describe_components(components):
    """Short description of the data types present, e.g. 'full impedance, tipper'."""
    comps = set(components)
    parts = []
    if {'ZXXR', 'ZYYR', 'ZXYR', 'ZYXR'} <= comps:
        parts.append('full impedance')
    elif {'ZXYR', 'ZYXR'} <= comps:
        parts.append('off-diagonal impedance')
    if {'TZXR', 'TZYR'} <= comps:
        parts.append('tipper')
    if {'PTXX', 'PTYY'} <= comps:
        parts.append('phase tensor')
    return ', '.join(parts) or 'unknown'


def write_summary(dataset, outfile, files):
    """Write the text summary. ``files`` maps 'data' / 'response' / 'list' / 'model' to paths."""
    rms = dataset.rms
    with open(outfile, 'w') as f:
        f.write('Project Name: \n')
        f.write('Directory: {}\n'.format(os.getcwd()))
        f.write('Report generated: {}\n'.format(datetime.date.today().isoformat()))
        f.write('Site list: {} \n'.format(files.get('list') or '--'))
        f.write('Data file: {} \n'.format(files['data']))
        f.write('Response file: {} \n'.format(files['response']))
        if files.get('model'):
            f.write('Model file: {} ({} x {} x {} cells)\n'.format(files['model'], dataset.model.nx,
                                                                  dataset.model.ny, dataset.model.nz))
        else:
            f.write('Model file: --\n')
        f.write('Inversion type: {} ({})\n'.format(dataset.data.inv_type,
                                                   describe_components(dataset.data.components)))
        f.write('Azimuth: {}\n'.format(dataset.data.azimuth))
        f.write('Number of sites: {}\n'.format(dataset.data.NS))
        f.write('Error floors: \n')
        f.write('Starting half-space rho: \n')
        f.write('\n')
        f.write('Final RMS: {:.5f} \n'.format(rms['Total']))
        f.write('RMS by component:\n')
        for comp in dataset.data.components:
            if comp in rms['Component']:
                f.write('{:>10}{:12.3f}\n'.format(comp, rms['Component'][comp]))
        worst = sorted(rms['Station'].items(), key=lambda item: item[1]['Total'], reverse=True)[:5]
        f.write('Worst-fitting sites:\n')
        for site, info in worst:
            f.write('{:>10}{:12.3f}\n'.format(site, info['Total']))
        f.write('\n')
        f.write('{} Inverted Periods: \n'.format(dataset.data.NP))
        f.write('{:>20}{:>20}{:>20}\n'.format('Period (s)', 'Frequency (Hz)', 'Log10(Period)'))
        for p in dataset.data.periods:
            f.write('{:20.5f}{:20.5f}{:20.5f}\n'.format(p, 1 / p, np.log10(p)))
        f.write('Notes:')


# --------------------------------------------------------------------------- #
# Misfit curves
# --------------------------------------------------------------------------- #
def _flagged(site, comp):
    """True where ``comp`` was excluded from the inversion (REMOVE_FLAG error)."""
    return np.asarray(site.used_error[comp]) == site.REMOVE_FLAG


def _plot_pair(ax, periods, data, errors, resp, flagged, colour, alpha, filled=True, dashed=False, label=None):
    """Data as points with error bars, response as a line; flagged points are left out."""
    data = np.where(flagged, np.nan, data)
    errors = np.where(flagged, np.nan, errors)
    ax.errorbar(periods, data, yerr=errors, linestyle='', marker='o', markersize=4,
                color=colour, mec='k', mew=0.5, mfc=colour if filled else 'w', alpha=alpha, label=label)
    ax.plot(periods, resp, linestyle='--' if dashed else '-', color=colour, alpha=alpha)


def _rho_phase_panels(axes, d, r, has_diag):
    rho_ax, pha_ax = axes
    off_diagonal = []
    for comp in ('xy', 'yx') + (('xx', 'yy') if has_diag else ()):
        flagged = _flagged(d, 'Z{}R'.format(comp.upper()))
        rho, rho_err, _ = utils.compute_rho(d, calc_comp=comp, errtype='used_error')
        rho_resp, _, _ = utils.compute_rho(r, calc_comp=comp)
        _plot_pair(rho_ax, d.periods, np.real(rho), rho_err, np.real(rho_resp), flagged,
                   COLOURS[comp], ALPHA[comp], label=comp.upper())
        if comp in ('xy', 'yx'):
            off_diagonal += list(np.real(rho)[~flagged]) + list(np.real(rho_resp))
            pha, pha_err = utils.compute_phase(d, calc_comp=comp, errtype='used_error', wrap=1)
            pha_resp, _ = utils.compute_phase(r, calc_comp=comp, wrap=1)
            _plot_pair(pha_ax, d.periods, pha, pha_err, pha_resp, flagged, COLOURS[comp], ALPHA[comp])
    rho_ax.set_yscale('log')
    # Limits from the off-diagonal curves; small diagonal values are clipped
    off_diagonal = np.array([x for x in off_diagonal if np.isfinite(x) and x > 0])
    if off_diagonal.size:
        rho_ax.set_ylim(10 ** (np.floor(np.log10(off_diagonal.min())) - 1),
                        10 ** (np.ceil(np.log10(off_diagonal.max())) + 0.5))
    rho_ax.set_ylabel(r'$\rho_a$ ($\Omega$m)')
    pha_ax.set_ylim(0, 120)
    pha_ax.set_ylabel('Phase (deg)')


def _impedance_panels(axes, d, r, has_diag):
    groups = [('xy', 'yx')] + ([('xx', 'yy')] if has_diag else [])
    root_t = np.sqrt(d.periods)
    for ax, comps in zip(axes, groups):
        for comp in comps:
            for part, filled, dashed in (('R', True, False), ('I', False, True)):
                key = 'Z{}{}'.format(comp.upper(), part)
                _plot_pair(ax, d.periods, d.data[key] * root_t, d.used_error[key] * root_t,
                           r.data[key] * root_t, _flagged(d, key), COLOURS[comp], ALPHA[comp],
                           filled=filled, dashed=dashed,
                           label='{} {}'.format(comp.upper(), 'Re' if part == 'R' else 'Im'))
        ax.set_ylabel(r'Z$\cdot\sqrt{T}$' + ('\n(off-diagonal)' if comps[0] == 'xy' else '\n(diagonal)'))


def _tipper_panel(ax, d, r):
    for comp in ('tzx', 'tzy'):
        for part, filled, dashed in (('R', True, False), ('I', False, True)):
            key = '{}{}'.format(comp.upper(), part)
            _plot_pair(ax, d.periods, d.data[key], d.used_error[key], r.data[key], _flagged(d, key),
                       COLOURS[comp], ALPHA[comp], filled=filled, dashed=dashed,
                       label='{} {}'.format(comp[1:].upper(), 'Re' if part == 'R' else 'Im'))
    ax.axhline(0, color='0.6', linewidth=0.5)
    ax.set_ylabel('Tipper')


def plot_misfit_page(fig, data_sites, resp_sites, site_rms, plot_type='rho-phase'):
    """Draw up to four sites (2 x 2) on ``fig``, each as a column of panels."""
    outer = GridSpec(2, 2, figure=fig, hspace=0.25, wspace=0.25)
    for kk, (d, r) in enumerate(zip(data_sites, resp_sites)):
        has_diag = {'ZXXR', 'ZYYR'} <= set(d.components)
        has_tipper = {'TZXR', 'TZYR'} <= set(d.components) and not _flagged(d, 'TZXR').all()
        if not {'ZXYR', 'ZYXR'} <= set(d.components):
            main_heights = []   # tipper-only data
        elif plot_type == 'rho-phase':
            main_heights = [3, 2]
        else:
            main_heights = [2, 2] if has_diag else [3]
        heights = main_heights + ([2] if has_tipper else [])
        inner = outer[kk].subgridspec(len(heights), 1, height_ratios=heights, hspace=0.08)
        axes = [fig.add_subplot(inner[0])]
        axes += [fig.add_subplot(inner[ii], sharex=axes[0]) for ii in range(1, len(heights))]
        if main_heights and plot_type == 'rho-phase':
            _rho_phase_panels(axes[:2], d, r, has_diag)
        elif main_heights:
            _impedance_panels(axes[:len(main_heights)], d, r, has_diag)
        if has_tipper:
            _tipper_panel(axes[-1], d, r)
        for ax in axes:
            ax.set_xscale('log')
            ax.grid(True, which='major', alpha=0.3)
            if ax is not axes[-1]:
                plt.setp(ax.get_xticklabels(), visible=False)
        axes[-1].set_xlabel('Period (s)')
        axes[0].set_title('{}  (RMS {:.2f})'.format(d.name, site_rms.get(d.name, np.nan)))
        if kk == 0:
            axes[0].legend(loc='best', fontsize=7, ncol=2)
            if has_tipper:
                axes[-1].legend(loc='best', fontsize=7, ncol=2)
    fig.legend(handles=[Line2D([], [], marker='o', linestyle='', color='0.4', label='data'),
                        Line2D([], [], color='0.4', label='response')],
               loc='upper center', ncol=2, frameon=False)


def generate_misfit_plots(dataset, outfile, plot_type='rho-phase'):
    """Write all sites to one multi-page PDF, four per page. Returns the number of pages."""
    site_rms = {site: info['Total'] for site, info in dataset.rms['Station'].items()}
    names = dataset.data.site_names
    pages = 0
    with PdfPages(outfile) as pdf:
        for ii in range(0, len(names), SITES_PER_PAGE):
            fig = plt.figure(figsize=(12, 16))
            chunk = names[ii:ii + SITES_PER_PAGE]
            plot_misfit_page(fig, [dataset.data.sites[s] for s in chunk],
                             [dataset.response.sites[s] for s in chunk], site_rms, plot_type)
            pdf.savefig(fig)
            plt.close(fig)
            pages += 1
    return pages


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='autogen_report',
        description='Write a misfit report (text summary + a PDF of misfit curves) for an inversion. '
                    'Any data format pyMT reads is accepted; the format is taken from the file.')
    parser.add_argument('-d', '--data', required=True, help='inversion data file')
    parser.add_argument('-r', '--response', required=True, help='inversion response (predicted data) file')
    parser.add_argument('-l', '--list', help='site list; names the sites for formats that do not store '
                                             'names (e.g. WSINV3DMT)')
    parser.add_argument('-m', '--model', help='model file (listed in the summary)')
    parser.add_argument('-o', '--out-dir', default='Report', help='output folder (default: Report)')
    parser.add_argument('--plot', choices=('rho-phase', 'impedance'), default='rho-phase',
                        help='plot apparent resistivity and phase (default) or the impedance tensor; '
                             'the tipper is added in both cases when present')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwrite an existing report in the output folder')
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    for label, path in (('data', args.data), ('response', args.response),
                        ('list', args.list), ('model', args.model)):
        if path and not os.path.isfile(path):
            print('autogen_report: {} file {} not found'.format(label, path))
            return 1
    summary = os.path.join(args.out_dir, 'summary.txt')
    pdf_file = os.path.join(args.out_dir, 'misfit_curves.pdf')
    if (os.path.exists(summary) or os.path.exists(pdf_file)) and not args.force:
        print('autogen_report: a report already exists in {} (use --force to overwrite)'.format(args.out_dir))
        return 1
    dataset = WSDS.Dataset(datafile=args.data, responsefile=args.response,
                           listfile=args.list or '', modelfile=args.model or '')
    if not dataset.rms:
        print('autogen_report: could not compute an RMS -- do the data and response files match?')
        return 1
    comps = set(dataset.data.components)
    if not ({'ZXYR', 'ZYXR'} <= comps or {'TZXR', 'TZYR'} <= comps):
        print('autogen_report: no impedance or tipper data to plot (phase-tensor data is not supported yet)')
        return 1
    if args.plot == 'rho-phase' and not {'ZXYR', 'ZYXR'} <= comps:
        print('autogen_report: no impedance data, so plotting the tipper only')
    os.makedirs(args.out_dir, exist_ok=True)
    write_summary(dataset, summary, {'data': args.data, 'response': args.response,
                                     'list': args.list, 'model': args.model})
    pages = generate_misfit_plots(dataset, pdf_file, args.plot)
    print('Wrote {} and {} ({} page{})'.format(summary, pdf_file, pages, '' if pages == 1 else 's'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
