"""Add raw-data tipper to a GoFEM data file, taking tipper errors from two other data files.

Personal script -- edit the paths at the top and run.
"""
import gateway_mt.data_structures as DS
import numpy as np

local_path = 'G:/Other computers/My Computer/sync/Regions/'

listfile = local_path + 'ATHA/gofem/edi_interpolated/atha_mcarthur/allall.lst'
data_1 = local_path + 'ATHA/gofem/model6/wMcArthur/atha6-wMcArthur-allfreq_Z.gdat'
data_2 = local_path + 'ATHA/gofem/model6/wMcArthur/atha6-wMcArthur_ZK-MA-only.gdat'
data_3 = local_path + 'ATHA/gofem/model6/atha5_gofem-ZK_rm-sites.gdat'
raw = DS.RawData(listfile=listfile)
d1 = DS.Data(datafile=data_1)
d2 = DS.Data(datafile=data_2)
d3 = DS.Data(datafile=data_3)

dK = raw.get_data(sites=d1.site_names,
				  periods=d1.periods,
				  components=d1.TIPPER_COMPONENTS)


for site in d1.site_names:
	# Copy the data from raw into the data structure
	d1.sites[site].add_component(data=dK[site].data,
								 errors=dK[site].errors,
								 errmap=dK[site].errmap)
	for comp in d1.TIPPER_COMPONENTS:
		idx = dK[site].errmap[comp] == d1.sites[site].NO_PERIOD_FLAG
		if np.any(idx):
			d1.sites[site].errors[comp][idx] = d1.REMOVE_FLAG
			d1.sites[site].used_errors[comp][idx] = d1.REMOVE_FLAG
	# Copy the errors from d2 into d1
	if site in d2.site_names:
		for ii, period in enumerate(d1.periods):
			for comp in d1.TIPPER_COMPONENTS:
				if period in d2.periods:
					idx = np.argmin(np.abs(period - np.array(d2.periods)))
					d1.sites[site].used_error[comp][ii] = d2.sites[site].used_error[comp][idx]
	# Copy the errors from d3 into d1					
	if site in d3.site_names:
		for ii, period in enumerate(d1.periods):
			for comp in d1.TIPPER_COMPONENTS:
				if period in d3.periods:
					idx = np.argmin(np.abs(period - np.array(d3.periods)))
					d1.sites[site].used_error[comp][ii] = d3.sites[site].used_error[comp][idx]

d1.inv_type = 5
d1.write(local_path + '/ATHA/gofem/model6/joint-datasets/atha-mcarthur_ZK-spliced')