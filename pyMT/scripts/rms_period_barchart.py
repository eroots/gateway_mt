import numpy as np
import pyMT.data_structures as DS
import matplotlib.pyplot as plt
import pyMT.utils as utils


# listfile = 'E:/phd/NextCloud/data/Regions/MetalEarth/wst/j2/mantle/fullrun/wst_cullmantle.lst'
# data_file = 'E:/phd/NextCloud/data/Regions/MetalEarth/wst/fullmantle/cull/Z/ZK/wst_cullmantle3_LAMBERT_ZK_removed.dat'
# response_file = 'E:/phd/NextCloud/data/Regions/MetalEarth/wst/fullmantle/cull/Z/ZK/wstZK_lastIter.dat'
# save_path = 'E:/phd/NextCloud/Documents/ME_Transects/wst/rms_plots'
local_path = 'G:/Other Computers/My Computer/sync/Regions/'
save_path = 'G:/Other Computers/My Computer/sync/Documents/MSH/'
listfile = local_path + 'MSH/edis-all/empower/edi_interpolated/all_T1-rmdupes.lst'
data_file = local_path + 'MSH/model5/t1/mshT1_wTransect-ZK_rmdupes.gdat'
# data_file = local_path + 'MSH/edi_ffmt/mshT1-wTransectZK_rmdupes_unifiedError.dat'
response_file = local_path + 'MSH/model5/t1/ZK/finalize/inv_data_21.gdat'
response_file2 = local_path + 'MSH/model5/t1/ftest/inv_data_0-rmdupes.gdat'
response_file3 = local_path + 'MSH/model5/swap_data/t1/t1_swap_0.gdat'
# listfile = local_path + 'MSH/edis-all/empower/edi_interpolated/all_T2-wTransect-rmdupes.lst'
# # data_file = local_path + 'MSH/model5/t2/MSHT2_wTransect-ZK_rmdupes.gdat'
# data_file = local_path + 'MSH/edi_ffmt/mshT2-wTransectZK_rmdupes_unifiedError.dat'
# response_file = local_path + 'MSH/model5/t2/ZK/finalize/inv_data_21.gdat'
# response_file2 = local_path + 'MSH/model5/t2/ftest/inv_data_0.gdat'

# data_files = [data_file]
resp_files = [response_file, response_file2, response_file3]
NP = 22
NR = 12
all_prms = np.zeros((NP, len(resp_files)))
all_pcrms = np.zeros((NP, NR, len(resp_files)))
for ir, resp in enumerate(resp_files):

    dataset = DS.Dataset(listfile=listfile, datafile=data_file, responsefile=resp)

    # for ii, p in enumerate(dataset.data.periods):
    # p_rms = {site: {x for x in range(dataset.data.NP)} for site in dataset.data.site_names}
    all_rms = np.zeros((dataset.data.NP, dataset.data.NS, dataset.data.NR))
    for ii, site in enumerate(dataset.data.site_names):
        for jj, comp in enumerate(dataset.data.components):
            for kk, p in enumerate(dataset.data.periods):
                if dataset.data.sites[site].used_errors[comp][kk] != dataset.data.REMOVE_FLAG:
                    all_rms[kk][ii][jj] = (np.abs(dataset.response.sites[site].data[comp][kk] - dataset.data.sites[site].data[comp][kk]) / 
                                           dataset.data.sites[site].used_error[comp][kk])
                else:
                    all_rms[kk][ii][jj] = np.nan

    # Collapse down to the period RMS only
    pc_rms = np.sqrt(np.nanmean(all_rms**2, axis=1))
    p_rms = np.sqrt(np.nanmean(all_rms**2, axis=2))
    p_rms = np.sqrt(np.nanmean(p_rms**2, axis=1))
    all_prms[:, ir] = p_rms
    all_pcrms[:, :, ir] = pc_rms


# fig, axes = plt.subplots(nrows=4, ncols=2)

# idx1, idx2 = 0, 0
# for ii in range(int(dataset.data.NR)):
# plt.rcParams['axes.titlepad'] = -14
# for ii in range(8):
#     comp = dataset.data.components[ii]

#     # Plot both RMS values on  top of each for each component
#     # rms = np.sqrt(np.mean(all_pcrms[:, ii*2:(ii*2)+2, 0], axis=1))
#     # axes.plot(np.log10(dataset.data.periods), rms, 'k-', label='T1 Final Model')
#     # rms = np.sqrt(np.mean(all_pcrms[:, ii*2:2*ii+2, 1], axis=1))
#     # axes.plot(np.log10(dataset.data.periods), rms, 'b--', label='T1 Feature Test')
#     # axes[idx1, idx2].plot(np.log10(dataset.data.periods), all_pcrms[:, ii, 0], 'k-', label='T1 Final Model')
#     # axes[idx1, idx2].plot(np.log10(dataset.data.periods), all_pcrms[:, ii, 1], 'b--', label='T1 Final Model')
    
    
#     # Plot the difference between the RMS of each component by period
#     axes[idx1, idx2].plot(np.log10(dataset.data.periods), all_pcrms[:, ii, 1] - all_pcrms[:, ii, 0], 'k-')
#     axes[idx1, idx2].set_title(comp, y=1, pad=-14)
#     axes[idx1, idx2].set_ylim([-0.1, 0.75])
#     if idx2 == 0:
#         axes[idx1, idx2].set_ylabel('RMS', fontsize=14)
#     if idx1 == 3:
#         axes[idx1, idx2].set_xlabel('Log10 Period (s)', fontsize=14)
#     if idx2 > 0:
#         idx2 = 0
#         idx1 +=1
#     else:
#         idx2 += 1

# Plot as a single figure with the RMS for both on top of each other
fig, axes = plt.subplots(nrows=1, ncols=1)
plt.plot(np.log10(dataset.data.periods), all_prms[:,0], 'k-', label='T1 Final Model')
plt.plot(np.log10(dataset.data.periods), all_prms[:,1], 'b--', label='T1 Data + T2 Embedded C1')
plt.plot(np.log10(dataset.data.periods), all_prms[:,2], color='orange', linestyle='--', label='T1 Data + T2 Model')
twin = axes.twinx()
twin.plot(np.log10(dataset.data.periods), all_prms[:,2], color='orange', linestyle='--', label='T1 + Full T2')
axes.set_xlabel('Log10 Period (s)', fontsize=14)
axes.set_ylabel('nRMS Misfit', fontsize=14)
twin.tick_params(axis='y', labelcolor='orange')
axes.set_ylim([0.5, 1.8])
twin.set_ylim([2, 7])
axes.legend()
plt.show()


