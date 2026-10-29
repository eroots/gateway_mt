"""Write a misfit report for an inversion: a text summary plus a PDF of misfit curves.

Usage::

    autogen_report -d DATA -r RESPONSE [-l LIST] [-m MODEL] [-o REPORT_DIR]
                   [--plot rho-phase | impedance]

Writes ``REPORT_DIR/summary.txt`` (input files, inversion type, sites and
periods, overall / per-component RMS, worst-fitting sites, the period table,
and blank fields to fill in by hand) and ``REPORT_DIR/misfit_curves.pdf``:
data (points with error bars) against response (lines), four sites per page,
in the same layout as ``plot_dataset``. Each site shows apparent resistivity
and phase (default) or the impedance tensor (``--plot impedance``), with the
tipper underneath when present; with a list file the raw data is added as
faint points.
Works with any data format gateway_mt reads (ModEM, WSINV3DMT, GoFEM, ...); the
format is taken from the file. Installed as the ``autogen_report`` script.
"""
import argparse
import datetime
import os
import sys

import numpy as np

import gateway_mt.data_structures as WSDS
from gateway_mt.scripts import plot_dataset


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


def generate_misfit_plots(dataset, outfile, plot_type='rho-phase', include_raw=False):
    """Write data vs. response for every site to ``outfile`` (see :mod:`gateway_mt.scripts.plot_dataset`).

    Returns the number of pages written.
    """
    site_rms = {site: info['Total'] for site, info in dataset.rms['Station'].items()}
    kinds = (['raw_data'] if include_raw else []) + ['data', 'response']
    return plot_dataset.write_pdf(dataset, outfile, kinds, plot_type, site_rms=site_rms)



def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='autogen_report',
        description='Write a misfit report (text summary + a PDF of misfit curves) for an inversion. '
                    'Any data format gateway_mt reads is accepted; the format is taken from the file.')
    parser.add_argument('-d', '--data', required=True, help='inversion data file')
    parser.add_argument('-r', '--response', required=True, help='inversion response (predicted data) file')
    parser.add_argument('-l', '--list', help='site list; adds the raw data to the plots and names the sites '
                                             'for formats that do not store names (e.g. WSINV3DMT)')
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
    pages = generate_misfit_plots(dataset, pdf_file, args.plot, include_raw=bool(args.list))
    print('Wrote {} and {} ({} page{})'.format(summary, pdf_file, pages, '' if pages == 1 else 's'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
