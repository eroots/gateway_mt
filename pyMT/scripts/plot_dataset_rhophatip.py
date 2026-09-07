import pyMT.utils as utils
from matplotlib.backends.backend_pdf import PdfPages
import pyMT.data_structures as DS
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import sys
import os
from copy import deepcopy
# from pyMT import gplot


plot_options = {'marker_dict': {'data': {'Rxy':  'o',
                                         'Rxx':  'v',
                                         'Ryy':  'v',
                                         'Ryx':  'o',
                                         'Pxy':  'o',
                                         'Pxx':  'v',
                                         'Pyy':  'v',
                                         'Pyx':  'o',
                                         'TZXR': 'o',
                                         'TZXI': 'o',
                                         'TZYR': 'o',
                                         'TZYI': 'o'},
                                 'resp': {'Rxy':  '',
                                          'Rxx':  '',
                                          'Ryy':  '',
                                          'Ryx':  '',
                                          'Pxy':  '',
                                          'Pxx':  '',
                                          'Pyy':  '',
                                          'Pyx':  '',
                                          'TZXR': '',
                                          'TZXI': '',
                                          'TZYR': '',
                                          'TZYI': ''}},

                'linestyle': {'data': {'Rxy':  '',
                                       'Rxx':  '',
                                       'Ryy':  '',
                                       'Ryx':  '',
                                       'Pxy':  '',
                                       'Pxx':  '',
                                       'Pyy':  '',
                                       'Pyx':  '',
                                       'TZXR': '',
                                       'TZXI': '',
                                       'TZYR': '',
                                       'TZYI': ''},
                              'resp': {'Rxy':  '-',
                                       'Rxx':  '--',
                                       'Ryy':  '-',
                                       'Ryx':  '--',
                                       'Pxy':  '-',
                                       'Pxx':  '--',
                                       'Pyy':  '--',
                                       'Pyx':  '-',
                                       'TZXR': '-',
                                       'TZXI': '--',
                                       'TZYR': '-',
                                       'TZYI': '--'}},
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
                'alpha':       {'Rxy':  0.8,
                                'Rxx':  0.5,
                                'Ryy':  0.5,
                                'Ryx':  0.8,
                                'Pxy':  0.8,
                                'Pxx':  0.5,
                                'Pyy':  0.5,
                                'Pyx':  0.8,
                                'TZXR': 0.8,
                                'TZXI': 0.8,
                                'TZYR': 0.8,
                                'TZYI': 0.8},
                'edgecolor': 'k',
                'markersize': 4,
                'edgewidth': 1,
                'x_limits': [1e-3, 1e4],
                'y_limits':     {'R': [1e-0, 10**4.5],
                                 'P': [-180, 180],
                                 'T': [-0.7, 0.7]}
                                 }

def plot_sites(sites, fig, resp_sites=None, data_type=['data'], save_fig=False, plots_per_row=2, use_legend=True):
    site_names = [site.name for site in sites['data']] # Assume at minimum the data has been passed in
    for ss, site_name in enumerate(site_names):
        to_plot = []
        # Assume 6 sites per page
        # Always plot tipper, just make it empty if its not there
        # if 'ZXYR' in site.components:
        #     to_plot.append('Rxy')   
        #     to_plot.append('Ryx')
        #     to_plot.append('Rxx')
        #     to_plot.append('Ryy')
        #     # to_plot.append('PXX')
        #     to_plot.append('Pxy')
        #     to_plot.append('Pyx')
        #     # to_plot.append('PYY')
        # if 'TZXR' in site.components:
        #     to_plot.append('TZXR')
        #     to_plot.append('TZYR')
        #     to_plot.append('TZXI')
        #     to_plot.append('TZYI')
        # if 'ZXYR' in site.components:
        to_plot.append('Rxy')   
        to_plot.append('Ryx')
        to_plot.append('Rxx')
        to_plot.append('Ryy')
        to_plot.append('Pxx')
        to_plot.append('Pxy')
        to_plot.append('Pyx')
        to_plot.append('Pyy')
        if 'TZXR' in sites['data'][ss].components:
            to_plot.append('TZXR')
            to_plot.append('TZYR')
            to_plot.append('TZXI')
            to_plot.append('TZYI')

        axes = {'rho': None,
                'phase': None,
                'tzr': None,
                'tzy': None}
        spec = {'rho': None,
                'phase': None,
                'tzr': None,
                'tzi': None}
        # fig = plt.figure(figsize=(10, 12))
        if ss < 2:
            row_number = 0
            col_number = ss*8 + ss
        else:
            row_number = 9
            col_number = (ss % plots_per_row) * 8 + (ss % plots_per_row)
        if not ('T' in [x[0] for x in to_plot]):
            bottom_plot = 'phase'
            top_plot = 'rho'
            spec['rho'] = GridSpec(17, 17, figure=fig).new_subplotspec((row_number, col_number), colspan=8, rowspan=5)
            spec['phase'] = GridSpec(17, 17, figure=fig).new_subplotspec((row_number+5, col_number), colspan=8, rowspan=3)
            # spec['tz'] = GridSpec(17, 17, figure=fig).new_subplotspec((row_number+6, col_number), colspan=8, rowspan=2)
            # spec['tzi'] = GridSpec(6, 6, figure=fig).new_subplotspec((5, 0), colspan=2, rowspan=1)
            axes['rho'] = fig.add_subplot(spec['rho'])
            axes['phase'] = fig.add_subplot(spec['phase'])
            # axes['tzr'] = fig.add_subplot(spec['tzr'])
            # axes['tzi'] = fig.add_subplot(spec['tzi'])
        elif len(to_plot) == 4:
            top_plot = 'tzr'
            spec['tzr'] = GridSpec(2, 2, figure=fig).new_subplotspec((0, 0), colspan=2, rowspan=1)
            spec['tzi'] = GridSpec(2, 2, figure=fig).new_subplotspec((1, 0), colspan=2, rowspan=1)
            spec['tzr'] = fig.add_subplot(spec['tzr'])
            spec['tzi'] = fig.add_subplot(spec['tzi'])
        else:
            bottom_plot = 'tz'
            top_plot = 'rho'
            spec['rho'] = GridSpec(17, 17, figure=fig).new_subplotspec((row_number, col_number), colspan=8, rowspan=4)
            spec['phase'] = GridSpec(17, 17, figure=fig).new_subplotspec((row_number+4, col_number), colspan=8, rowspan=2)
            spec['tz'] = GridSpec(17, 17, figure=fig).new_subplotspec((row_number+6, col_number), colspan=8, rowspan=2)
            # spec['tzi'] = GridSpec(6, 6, figure=fig).new_subplotspec((5, 0), colspan=2, rowspan=1)
            axes['rho'] = fig.add_subplot(spec['rho'])
            axes['phase'] = fig.add_subplot(spec['phase'])
            axes['tzr'] = fig.add_subplot(spec['tzr'])
            axes['tzi'] = fig.add_subplot(spec['tzi'])



        for dtype in data_type:
            site = sites[dtype][ss]
            for tp in to_plot:
                if dtype == 'data' and ss == 0:
                    label = tp
                else:
                    label = ''
                if tp.startswith('R'):
                    rho, e, log_e = utils.compute_rho(site, calc_comp=tp[1:], errtype='errors')
                    use_periods = deepcopy(site.periods)
                    idx = site.used_error['ZXYR'] == site.REMOVE_FLAG
                    use_errors = e
                    if np.any(idx):
                        for kk, flag in enumerate(idx):
                            if flag:
                                rho[kk] = np.nan
                                use_periods[kk] = np.nan 
                                use_errors[kk] = np.nan
                    if dtype == 'resp':
                        use_errors = None
                    axes['rho'].errorbar((use_periods), (rho), xerr=None, yerr=use_errors,
                                         marker=plot_options['marker_dict'][dtype][tp],
                                         linestyle=plot_options['linestyle'][dtype][tp],
                                         color=plot_options['colour_dict'][tp],
                                         mec=plot_options['edgecolor'],
                                         markersize=plot_options['markersize'],
                                         mew=plot_options['edgewidth'],
                                         alpha=plot_options['alpha'][tp],
                                         label=label)
                    axes['rho'].set_ylabel('App. Rho (ohm-m)')
                    axes['rho'].set_ylim(plot_options['y_limits']['R'])
                    axes['rho'].set_xlim(plot_options['x_limits'])
                    axes['rho'].set_xscale('log')
                    axes['rho'].set_yscale('log')
                elif tp.startswith('P'):
                    pha, e = utils.compute_phase(site, calc_comp=tp[1:], errtype='errors', wrap=False)
                    use_periods = deepcopy(site.periods)
                    idx = site.used_error['ZXYR'] == site.REMOVE_FLAG
                    use_errors = e
                    if np.any(idx):
                        for kk, flag in enumerate(idx):
                            if flag:
                                pha[kk] = np.nan
                                use_periods[kk] = np.nan
                                use_errors[kk] = np.nan
                    if dtype == 'resp':
                        use_errors = None
                    axes['phase'].errorbar((site.periods), pha, xerr=None, yerr=use_errors,
                                           marker=plot_options['marker_dict'][dtype][tp],
                                           linestyle=plot_options['linestyle'][dtype][tp],
                                           color=plot_options['colour_dict'][tp],
                                           mec=plot_options['edgecolor'],
                                           markersize=plot_options['markersize'],
                                           mew=plot_options['edgewidth'],
                                           alpha=plot_options['alpha'][tp],
                                           label=label)
                    axes['phase'].set_ylabel('Phase (Degrees)')
                    axes['phase'].set_ylim(plot_options['y_limits']['P'])
                    axes['phase'].set_xlim(plot_options['x_limits'])
                    axes['phase'].set_xscale('log')
                elif tp.startswith('T'):
                    # comp = 'tz{}'.format(tp[-1].lower())
                    comp = 'tz'
                    if tp in site.components:
                        use_periods = deepcopy(site.periods)
                        K = deepcopy(site.data[tp])
                        use_errors = deepcopy(site.used_error[tp])
                        idx = site.used_error['TZXR'] == site.REMOVE_FLAG
                        if np.any(idx):
                            for kk, flag in enumerate(idx):
                                if flag:
                                    # print(kk)
                                    K[kk] = np.nan
                                    use_periods[kk] = np.nan 
                                    use_errors[kk] = np.nan
                        if dtype == 'resp':
                            use_errors = None
                        axes[comp].errorbar((site.periods), 
                                             K,
                                             xerr=None, yerr=use_errors,
                                             marker=plot_options['marker_dict'][dtype][tp],
                                             linestyle=plot_options['linestyle'][dtype][tp],
                                             color=plot_options['colour_dict'][tp],
                                             mec=plot_options['edgecolor'],
                                             markersize=plot_options['markersize'],
                                             mew=plot_options['edgewidth'],
                                             alpha=plot_options['alpha'][tp],
                                             label=label)
                    axes[comp].set_xscale('log')
                    axes[comp].grid(True)
                    axes[comp].set_ylabel('Tz')
                    # if comp.endswith('r'):
                    #     axes[comp].set_ylabel(r'Re {Tz}')
                    # else:
                    #     axes[comp].set_ylabel(r'Im {Tz}')
                    axes[comp].set_ylim(plot_options['y_limits']['T'])
                    axes[comp].set_xlim(plot_options['x_limits'])
            for ax_name, ax in axes.items():
                if ax:
                    if ss == 0 and use_legend:
                        ax.legend(loc='upper left')
                    if ax_name != bottom_plot:
                        ax.set_xticklabels([])
                    if col_number > 0:
                        ax.set_yticklabels([])
                        ax.set_ylabel('')
        axes[bottom_plot].set_xlabel('Period (s)')
        axes[top_plot].set_title(site.name)
    if save_fig:
        plt.savefig(site.name)
        plt.close()
    else:
        return fig

# def main():
#     if sys.argv[1].lower() == 'all':
#         for site in os.listdir():
#             if site.endswith('.edi'):
#                 data = DS.RawData(site)
#                 plot_site(data.sites[data.site_names[0]], save_fig=True)
#     else:
#         data = DS.RawData(sys.argv[1])
#         plot_site(data.sites[data.site_names[0]], save_fig=False)


# if __name__ == '__main__':
#     main()
    # file = 'E:/Work/Regions/ATHA/Atha21_emtf/zxx/10124_2021-09-15-222830/13_flipHy_SL1.edi'
    
    # data = DS.RawData(file)
    




# listfiles = [r'G:/Other computers/My Computer/sync/Regions/ATHA/j2/all_wTHOT.lst']
# datafiles = [r'G:/Other computers/My Computer/sync/Regions/ATHA/gofem/model6/atha5_gofem-ZK_rm-sites.gdat']
# respfiles = [r'G:/Other computers/My Computer/sync/Regions/ATHA/gofem/model6/ZK2/inv_data_47.gdat']

# listfiles = ['G:/Other computers/My Computer/sync/Regions/ATHA/gofem/edi_interpolated/atha_mcarthur/mcarthur.lst']
# datafiles = ['G:/Other computers/My Computer/sync/Regions/ATHA/gofem/model6/wMcArthur/atha6-wMcArthur_ZK-MA-only.gdat']
# respfiles = ['G:/Other computers/My Computer/sync/Regions/ATHA/gofem/model6/wMcArthur/MA-only/ZK/inv2_data_13.gdat']
# listfiles = ['G:/Other computers/My Computer/sync/Regions/ATHA/gofem/edi_interpolated/atha_mcarthur/allall.lst']
# datafiles = ['G:/Other computers/My Computer/sync/Regions/ATHA/gofem/model6/wMcArthur/combined/atha-mcarthur_ZK-spliced.gdat']
# respfiles = ['G:/Other computers/My Computer/sync/Regions/ATHA/gofem/model6/wMcArthur/combined/atha_combined_ub2_7.gdat']
# listfiles = ['G:/Other computers/My Computer/sync/Regions/Yukon/edis/fix/edi_interpolated/region_A/region_A.lst',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/edis/fix/edi_interpolated/region_B/region_B_cull.lst',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/edis/fix/edi_interpolated/region_C/region_C_cull2.lst',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/edis/fix/edi_interpolated/region_D-expanded/region_D-expanded.lst']

# datafiles = ['G:/Other computers/My Computer/sync/Regions/Yukon/a_model1/region_A_ZK-fix.gdat',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/final_datasets/region_B/region_B_ZK.gdat',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/final_datasets/region_C/region-C_ZK.gdat',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/d_model1/region-D_ZK.gdat']

# respfiles = ['G:/Other computers/My Computer/sync/Regions/Yukon/a_model1/ZK/yuk_A-ZK_31.gdat',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/final_datasets/region_B/yuk_B-ZK2_47.gdat',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/final_datasets/region_C/yuk_C-ZK_final.gdat',
#              'G:/Other computers/My Computer/sync/Regions/Yukon/d_model1/ZK-cont/yuk-D_ZK_38.gdat']

listfiles = ['G:/Other computers/My Computer/sync/Regions/undercover/j2/und_cull4.lst']
datafiles = ['G:/Other computers/My Computer/sync/Regions/undercover/gofem/model3/und_cull4-Z.gdat']
respfiles = ['G:/Other computers/My Computer/sync/Regions/undercover/gofem/model3/Z/und_Z2_29.gdat']

# datasets = ['Region_A', 'Region_B', 'Region_C', 'Region_D']
datasets = ['und2025-Z']

out_path = 'G:/Other computers/My Computer/sync/Documents/undercover/'

for dd, name in enumerate(datasets):
    # if dd < 2:
    #     continue

    dataset = DS.Dataset(datafile=datafiles[dd], responsefile=respfiles[dd])
    out_name = '{}_misfit-plots.pdf'.format(name)

# listfile = r'E:/Work/sync/Regions/ATHA/j2/tie_line.lst'
# dataset = DS.Dataset(listfile=listfile)
# dataset = DS.Dataset(datafile=datafile, responsefile=respfile)



    sites_per_page = 4
    with PdfPages(''.join([out_path, out_name])) as pdf:
        for ii in range(0, dataset.data.NS, sites_per_page):
        # for ii in range(0, 24, sites_per_page):
            fig = plt.figure(figsize=(12, 16))
            site_names = dataset.data.site_names[ii:ii + sites_per_page]
            # data_sites = [dataset.raw_data.sites[site] for site in site_names]
            data_sites = [dataset.data.sites[site] for site in site_names]
            resp_sites = [dataset.response.sites[site] for site in site_names]
            sites = {'data': data_sites, 'resp': resp_sites}
            plot_sites(sites, fig=fig, data_type=['data', 'resp'], save_fig=False, plots_per_row=2, use_legend=True)
            pdf.savefig()
    plt.close()
