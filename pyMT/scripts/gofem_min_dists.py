import numpy as np
import pyMT.data_structures as DS


csv_file = 'G:/Other computers/My Computer/sync/Regions/nevada/merlin/model6/receivers_projected.csv'
list_file = []
# list_file = 'G:/Other Computers/My Computer/sync/Regions/Yukon/edis/fix/edi_interpolated/region_B/region_B_cull.lst'
# list_file = r'G:/Other computers/My Computer/sync/Regions/Norway/GoFEM/data/ediNorway-12-2025/edi_flipped/allall.lst'
threshold = 300

if list_file:
	d = DS.RawData(list_file)
	locs = d.locations
	site_names = d.site_names
else:
	site_names = []
	x, y = [], []
	with open(csv_file, 'r') as f:
		lines = f.readlines()
	for line in lines:
		site_names.append(line.split()[1])
		x.append(float(line.split()[3]))
		y.append(float(line.split()[4]))
	locs = np.array((x, y)).T

# Change stuff if you're using the gofem location file instead

min_dist = []
sites = []
for ii, site in enumerate(site_names):
	dists = np.sqrt((locs[ii, 0] - locs[:, 0])**2 + (locs[ii, 1] - locs[:, 1])**2)
	idx = np.argpartition(dists, 1)[1]
	if dists[idx] < threshold:
		min_dist.append(dists[idx])
		sites.append([site, site_names[idx]])

for (d, s) in zip(min_dist, sites):
	print('{}, {}, {}'.format(round(d, 2), s[0], s[1]))

