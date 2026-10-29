"""Write a misfit report for an inversion: a text summary plus misfit-curve PDFs.

Usage::

    autogen_report -d DATA -r RESPONSE [-l LIST] [-m MODEL] [-o REPORT_DIR]

Writes ``REPORT_DIR/summary.txt`` (input files, inversion type, sites and
periods, overall / per-component RMS, worst-fitting sites, the period table,
and blank fields to fill in by hand) and one multi-page PDF per component
group in ``REPORT_DIR/misfit_curves/`` (data vs. response, six sites per
page). Works with any data format pyMT reads (ModEM, WSINV3DMT, GoFEM, ...);
the format is taken from the file. Installed as the ``autogen_report``
script.
"""
import argparse
import datetime
import os
import sys

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from numpy import log10

import pyMT.data_structures as WSDS
import pyMT.gplot as gplot

# Component groups plotted by generate_misfit_plots, and the data components each one needs
COMPONENT_GROUPS = {'impedance_diagonal': ('ZXXR', 'ZXXI', 'ZYYR', 'ZYYI'),
                    'impedance_off-diagonal': ('ZXYR', 'ZXYI', 'ZYXR', 'ZYXI'),
                    'phase_diagonal': ('PhaXX', 'PhaYY'),
                    'rho_diagonal': ('RhoXX', 'RhoYY'),
                    'rho_off-diagonal': ('RhoXY', 'RhoYX'),
                    'phase_off-diagonal': ('PhaXY', 'PhaYX'),
                    'tipper': ('TZXR', 'TZXI', 'TZYR', 'TZYI')}
REQUIRED = {'impedance_diagonal': ('ZXXR', 'ZYYR'),
            'impedance_off-diagonal': ('ZXYR', 'ZYXR'),
            'phase_diagonal': ('ZXXR', 'ZYYR'),
            'rho_diagonal': ('ZXXR', 'ZYYR'),
            'rho_off-diagonal': ('ZXYR', 'ZYXR'),
            'phase_off-diagonal': ('ZXYR', 'ZYXR'),
            'tipper': ('TZXR', 'TZYR')}


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
            f.write('{:20.5f}{:20.5f}{:20.5f}\n'.format(p, 1 / p, log10(p)))
        f.write('Notes:')


def generate_misfit_plots(dataset, out_dir='Report'):
    """One multi-page PDF per component group in ``out_dir/misfit_curves``, six sites per page."""
    path = os.path.join(out_dir, 'misfit_curves')
    os.makedirs(path, exist_ok=True)
    fig = plt.figure(figsize=(26, 15))
    dpm = gplot.DataPlotManager(fig=fig)
    components = {name: comps for name, comps in COMPONENT_GROUPS.items()
                  if set(REQUIRED[name]) <= set(dataset.data.components)}
    dpm.show_outliers = True
    dpm.markersize = 6
    dpm.outlier_thresh = 5
    dpm.wrap = True
    written = []
    for comps in components.keys():
        if 'phase' in comps or 'rho' in comps:
            dpm.link_axes_bounds = True
        else:
            dpm.link_axes_bounds = False
        print('Generating {} plots...'.format(comps))
        dpm.components = components[comps]
        outfile = os.path.join(path, comps + '_misfit.pdf')
        with PdfPages(outfile) as pp:
            for ii in range(0, len(dataset.data.site_names), 6):
                if comps == 'tipper':
                    dpm.scale = 'none'
                else:
                    dpm.scale = 'sqrt(periods)'
                sites = dataset.data.site_names[ii:ii + 6]
                dpm.sites = dataset.get_sites(site_names=sites, dTypes='all')
                dpm.which_errors = ['data']
                dpm.plot_data()
                dpm.fig.savefig(pp, format='pdf', dpi=1000, bbox_inches='tight')
        written.append(outfile)
    plt.close(fig)
    return written


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='autogen_report',
        description='Write a misfit report (text summary + misfit-curve PDFs) for an inversion. '
                    'Any data format pyMT reads is accepted; the format is taken from the file.')
    parser.add_argument('-d', '--data', required=True, help='inversion data file')
    parser.add_argument('-r', '--response', required=True, help='inversion response (predicted data) file')
    parser.add_argument('-l', '--list', help='site list; adds raw data to the plots and names sites '
                                             'for formats that do not store names (e.g. WSINV3DMT)')
    parser.add_argument('-m', '--model', help='model file (listed in the summary)')
    parser.add_argument('-o', '--out-dir', default='Report', help='output folder (default: Report)')
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
    if os.path.exists(summary) and not args.force:
        print('autogen_report: {} already exists (use --force to overwrite)'.format(summary))
        return 1
    dataset = WSDS.Dataset(datafile=args.data, responsefile=args.response,
                           listfile=args.list or '', modelfile=args.model or '')
    if not dataset.rms:
        print('autogen_report: could not compute an RMS -- do the data and response files match?')
        return 1
    os.makedirs(args.out_dir, exist_ok=True)
    write_summary(dataset, summary, {'data': args.data, 'response': args.response,
                                     'list': args.list, 'model': args.model})
    pdfs = generate_misfit_plots(dataset, args.out_dir)
    print('Wrote {} and {} misfit PDF{} to {}'.format(summary, len(pdfs), '' if len(pdfs) == 1 else 's',
                                                      os.path.join(args.out_dir, 'misfit_curves')))
    return 0


if __name__ == '__main__':
    sys.exit(main())
