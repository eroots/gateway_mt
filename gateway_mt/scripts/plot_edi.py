"""Plot apparent resistivity / phase / tipper curves from EDI files.

Usage::

    plot_edi SITE.edi [SITE2.edi ...]   # show each EDI on screen
    plot_edi all                        # save a PNG for every .edi in the folder
    plot_edi all --format pdf -o figs   # ...as PDFs, into figs/
    plot_edi SITE.edi --save            # save instead of showing

Installed as the ``plot_edi`` console script.  Self-contained -- does not
use :mod:`gateway_mt.gplot`.
"""
import gateway_mt.utils as utils
import gateway_mt.data_structures as DS
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import argparse
import sys
import os


plot_options = {'marker_dict': {'Rxy':  'o',
								'Rxx':  'o',
								'Ryy':  'o',
								'Ryx':  'o',
								'Pxy':  'o',
								'Pxx':  'o',
								'Pyy':  'o',
								'Pyx':  'o',
								'TZXR': 'o',
								'TZXI': 'o',
								'TZYR': 'o',
								'TZYI': 'o'},
				'linestyle': '',
				'colour_dict': {'Rxy':  'b',
							    'Rxx':  'b',
							    'Ryy':  'r',
							    'Ryx':  'r',
							    'Pxy':  'b',
							    'Pxx':  'b',
							    'Pyy':  'r',
							    'Pyx':  'r',
							    'TZXR': 'b',
							    'TZXI': 'b',
							    'TZYR': 'r',
							    'TZYI': 'r'},
				'alpha': 	   {'Rxy':  1,
							    'Rxx':  0.5,
							    'Ryy':  0.5,
							    'Ryx':  1,
							    'Pxy':  1,
							    'Pxx':  0.5,
							    'Pyy':  0.5,
							    'Pyx':  1,
							    'TZXR': 1,
							    'TZXI': 1,
							    'TZYR': 1,
							    'TZYI': 1},
				'edgecolor': 'k',
				'markersize': 6,
				'edgewidth': 1,
				'x_limits': [1e-3, 1e4],
				'y_limits':     {'R': [1e0, 1e5],
								 'P': [0, 120],
								 'T': [-0.75, 0.75]}
								 }

def plot_site(site, save_fig=False, out_dir='.', fmt='png'):
	to_plot = []
	if 'ZXYR' in site.components:
		to_plot.append('Rxy')	
		to_plot.append('Ryx')
		to_plot.append('Rxx')
		to_plot.append('Ryy')
		# to_plot.append('PXX')
		to_plot.append('Pxy')
		to_plot.append('Pyx')
		# to_plot.append('PYY')
	if 'TZXR' in site.components:
		to_plot.append('TZXR')
		to_plot.append('TZYR')
		to_plot.append('TZXI')
		to_plot.append('TZYI')

	axes = {'rho': None,
			'phase': None,
			'tzr': None,
			'tzi': None}
	spec = {'rho': None,
			'phase': None,
			'tzr': None,
			'tzi': None}
	fig = plt.figure(figsize=(10, 12))
	bottom_plot = 'tzi'
	top_plot = 'rho'
	if len(to_plot) == 6:
		spec['rho'] = GridSpec(3, 3, figure=fig).new_subplotspec((0, 0), colspan=3, rowspan=2)
		spec['phase'] = GridSpec(3, 3, figure=fig).new_subplotspec((2, 0), colspan=3, rowspan=1)
		axes['rho'] = fig.add_subplot(spec['rho'])
		axes['phase'] = fig.add_subplot(spec['phase'])
		bottom_plot = 'phase'
	elif len(to_plot) == 4:
		top_plot = 'tzr'
		spec['tzr'] = GridSpec(2, 2, figure=fig).new_subplotspec((0, 0), colspan=2, rowspan=1)
		spec['tzi'] = GridSpec(2, 2, figure=fig).new_subplotspec((1, 0), colspan=2, rowspan=1)
		axes['tzr'] = fig.add_subplot(spec['tzr'])
		axes['tzi'] = fig.add_subplot(spec['tzi'])
	elif len(to_plot) > 6:
		spec['rho'] = GridSpec(6, 6, figure=fig).new_subplotspec((0, 0), colspan=6, rowspan=3)
		spec['phase'] = GridSpec(6, 6, figure=fig).new_subplotspec((3, 0), colspan=6, rowspan=1)
		spec['tzr'] = GridSpec(6, 6, figure=fig).new_subplotspec((4, 0), colspan=6, rowspan=1)
		spec['tzi'] = GridSpec(6, 6, figure=fig).new_subplotspec((5, 0), colspan=6, rowspan=1)
		axes['rho'] = fig.add_subplot(spec['rho'])
		axes['phase'] = fig.add_subplot(spec['phase'])
		axes['tzr'] = fig.add_subplot(spec['tzr'])
		axes['tzi'] = fig.add_subplot(spec['tzi'])

	for tp in to_plot:
		if tp.startswith('R'):
			rho, e, log_e = utils.compute_rho(site, calc_comp=tp[1:], errtype='errors')
			axes['rho'].errorbar((site.periods), (rho), xerr=None, yerr=e,
						         marker=plot_options['marker_dict'][tp],
						         linestyle=plot_options['linestyle'],
						         color=plot_options['colour_dict'][tp],
						         mec=plot_options['edgecolor'],
						         markersize=plot_options['markersize'],
						         mew=plot_options['edgewidth'],
						         alpha=plot_options['alpha'][tp],
						         label=tp)
			axes['rho'].set_ylabel('App. Rho (ohm-m)')
			axes['rho'].set_ylim(plot_options['y_limits']['R'])
			axes['rho'].set_xlim(plot_options['x_limits'])
			axes['rho'].set_xscale('log')
			axes['rho'].set_yscale('log')
		elif tp.startswith('P'):
			pha, e = utils.compute_phase(site, calc_comp=tp[1:], errtype='errors', wrap=True)
			axes['phase'].errorbar((site.periods), pha, xerr=None, yerr=e,
						           marker=plot_options['marker_dict'][tp],
						           linestyle=plot_options['linestyle'],
						           color=plot_options['colour_dict'][tp],
						           mec=plot_options['edgecolor'],
						           markersize=plot_options['markersize'],
						           mew=plot_options['edgewidth'],
						           alpha=plot_options['alpha'][tp],
						           label=tp)
			axes['phase'].set_ylabel('Phase (Degrees)')
			axes['phase'].set_ylim(plot_options['y_limits']['P'])
			axes['phase'].set_xlim(plot_options['x_limits'])
			axes['phase'].set_xscale('log')
		elif tp.startswith('T'):
			comp = 'tz{}'.format(tp[-1].lower())
			axes[comp].errorbar((site.periods), 
								 site.data[tp],
								 xerr=None, yerr=site.errors[tp],
						         marker=plot_options['marker_dict'][tp],
						         linestyle=plot_options['linestyle'],
						         color=plot_options['colour_dict'][tp],
						         mec=plot_options['edgecolor'],
						         markersize=plot_options['markersize'],
						         mew=plot_options['edgewidth'],
						         alpha=plot_options['alpha'][tp],
						         label=tp)
			axes[comp].set_xscale('log')
			axes[comp].grid(True)
			if comp.endswith('r'):
				axes[comp].set_ylabel(r'Re {Tz}')
			else:
				axes[comp].set_ylabel(r'Im {Tz}')
			axes[comp].set_ylim(plot_options['y_limits']['T'])
			axes[comp].set_xlim(plot_options['x_limits'])
	for ax_name, ax in axes.items():
		if ax:
			ax.legend(loc='upper left')
			if ax_name != bottom_plot:
				ax.set_xticklabels([])
	axes[bottom_plot].set_xlabel('Period (s)')
	axes[top_plot].set_title(site.name)
	if save_fig:
		plt.savefig(os.path.join(out_dir, '{}.{}'.format(site.name, fmt)))
		plt.close()
	else:
		plt.show()

def parse_args(argv=None):
	parser = argparse.ArgumentParser(
		prog='plot_edi',
		description='Plot apparent resistivity, phase and tipper from EDI files.')
	parser.add_argument('files', nargs='+', metavar='FILE',
						help="EDI file(s) to plot, or 'all' for every .edi in the current folder "
							 "('all' saves the figures instead of showing them)")
	parser.add_argument('--save', action='store_true',
						help='save the figures (SITE.png) instead of showing them')
	parser.add_argument('-o', '--out-dir', default='.',
						help='folder for saved figures (default: current folder)')
	parser.add_argument('--format', default='png', choices=('png', 'pdf'),
						help='file format for saved figures (default: png)')
	return parser.parse_args(argv)


def main(argv=None):
	args = parse_args(argv)
	save_fig = args.save
	files = args.files
	if len(files) == 1 and files[0].lower() == 'all':
		files = sorted(f for f in os.listdir() if f.endswith('.edi'))
		save_fig = True
		if not files:
			print('plot_edi: no .edi files in the current folder')
			return 1
	if save_fig:
		os.makedirs(args.out_dir, exist_ok=True)
	failed = 0
	for file in files:
		if not os.path.isfile(file):
			print('plot_edi: {} not found'.format(file))
			failed += 1
			continue
		try:
			data = DS.RawData(file)
			plot_site(data.sites[data.site_names[0]], save_fig=save_fig,
					  out_dir=args.out_dir, fmt=args.format)
		except Exception as err:  # e.g. WSFileError for an unreadable EDI
			print('plot_edi: skipped {} ({})'.format(file, getattr(err, 'message', err)))
			plt.close('all')
			failed += 1
	return 1 if failed else 0


if __name__ == '__main__':
	sys.exit(main())
