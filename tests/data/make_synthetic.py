"""Generate the synthetic MT test dataset in ``tests/data/synthetic/``.

The dataset stands in for real survey data in the tests and docs. It is
built from analytic 1-D layered-earth responses, made more realistic by:

* a different conductor resistivity for each polarisation (so the xy and
  yx curves split, as over a 2-D structure),
* a geoelectric strike that rotates across the grid (non-zero diagonals),
* galvanic distortion (static shift) at a few sites,
* a tipper anomaly centred on the conductor,
* Gaussian noise at the size of the quoted errors.

Run from the repository root::

    python tests/data/make_synthetic.py

Output is deterministic (fixed random seed, no dates), so re-running it
should reproduce the committed files exactly.
"""

import os
import sys
from copy import deepcopy

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, os.pardir, os.pardir)))

import gateway_mt.data_structures as DS          # noqa: E402
from gateway_mt.modelling_1D import wait_1D       # noqa: E402

OUT = os.path.join(HERE, 'synthetic')
EDI_DIR = os.path.join(OUT, 'edis')
SEED = 20261029

MU0 = 4 * np.pi * 1e-7
FIELD_UNITS = 1 / (MU0 * 1e3)       # ohm -> mV/km/nT  (~795.77)

# Survey layout: a 6 x 5 grid of sites, ~5 km apart, centred here
CENTRE_LAT, CENTRE_LON = 49.0, -84.5
NX_SITES, NY_SITES, SPACING_KM = 6, 5, 5.0
PERIODS = np.logspace(-3, np.log10(5000), 40)
DATA_PERIOD_IDX = list(range(2, 38, 2))        # 18 periods for the inversion files
SHORT_SITES = ('SYN029', 'SYN030')             # AMT-style: no long periods
NO_TIPPER_SITES = ('SYN030',)
STATIC_SHIFT = {'SYN005': (1.6, 1.0), 'SYN012': (1.0, 0.55), 'SYN023': (0.7, 1.3)}

# Conductor anomaly (local km, x = north, y = east)
ANOMALY_X, ANOMALY_Y, ANOMALY_WIDTH = 2.0, -3.0, 8.0


# --------------------------------------------------------------------------- #
# Layout and 1-D structure
# --------------------------------------------------------------------------- #
def site_layout(rng):
    """Site names, local x/y (km), lat/long and elevation."""
    sites = []
    km_per_deg_lat = 111.2
    km_per_deg_lon = 111.2 * np.cos(np.deg2rad(CENTRE_LAT))
    n = 0
    for j in range(NY_SITES):
        for i in range(NX_SITES):
            n += 1
            x = (j - (NY_SITES - 1) / 2) * SPACING_KM + rng.uniform(-0.8, 0.8)
            y = (i - (NX_SITES - 1) / 2) * SPACING_KM + rng.uniform(-0.8, 0.8)
            sites.append({'name': 'SYN{:03d}'.format(n),
                          'x': x, 'y': y,
                          'lat': CENTRE_LAT + x / km_per_deg_lat,
                          'lon': CENTRE_LON + y / km_per_deg_lon,
                          'elev': 300 + 40 * np.sin(x / 6) + 25 * np.cos(y / 4)})
    return sites


def anomaly_weight(x, y):
    """0..1 weight of the conductor anomaly at local (x, y) km."""
    r2 = (x - ANOMALY_X) ** 2 + (y - ANOMALY_Y) ** 2
    return np.exp(-r2 / (2 * ANOMALY_WIDTH ** 2))


def layered_model(x, y, conductor_rho):
    """Layer resistivities and interface depths (m) under (x, y)."""
    w = anomaly_weight(x, y)
    surface_rho = 400 * 2 ** (np.sin(x / 7) + np.cos(y / 9))
    conductor_top = (18 - 9 * w) * 1e3
    rho = np.array([surface_rho, 3000, conductor_rho, 2000, 100])
    depths = np.array([0, 1e3, conductor_top, conductor_top + 5e3, 120e3])
    return rho, depths


def clean_response(site, periods):
    """Noise-free Z (2x2, ohm) and tipper (2,) at each period."""
    x, y = site['x'], site['y']
    w = anomaly_weight(x, y)
    rho_te = 50 * (1 - w) + 3 * w          # conductor seen strongly in xy...
    rho_tm = 8 * rho_te                    # ...and weakly in yx
    zs = []
    for rho_c in (rho_te, rho_tm):
        rho, depths = layered_model(x, y, rho_c)
        zs.append(wait_1D(periods, rho, np.diff(depths)))
    z_xy, z_yx = zs[0], -zs[1]

    strike = np.deg2rad(25 + 3 * x)        # rotates across the grid
    c, s = np.cos(strike), np.sin(strike)
    rot = np.array([[c, -s], [s, c]])
    gx, gy = STATIC_SHIFT.get(site['name'], (1.0, 1.0))
    distortion = np.diag([gx, gy])

    Z = np.zeros((len(periods), 2, 2), dtype=complex)
    for k in range(len(periods)):
        z_principal = np.array([[0, z_xy[k]], [z_yx[k], 0]])
        Z[k] = distortion @ rot @ z_principal @ rot.T

    # Tipper: points away from the anomaly, peaks near 10 s
    dx, dy = x - ANOMALY_X, y - ANOMALY_Y
    r = np.hypot(dx, dy) + 1e-6
    amp = 0.6 * (r / ANOMALY_WIDTH) * np.exp(-(r / ANOMALY_WIDTH) ** 2 / 2)
    logp = np.log10(periods)
    shape = np.exp(-(logp - 1) ** 2 / 1.5) + 0.4j * (logp - 1) * np.exp(-(logp - 1) ** 2 / 1.2)
    T = np.stack([amp * dx / r * shape, amp * dy / r * shape], axis=1)
    return Z, T


def errors_for(Z, T):
    """Standard deviations: 3% of sqrt|Zxy Zyx| for every Z element, 0.02 for tipper."""
    z_err = 0.03 * np.sqrt(np.abs(Z[:, 0, 1] * Z[:, 1, 0]))
    t_err = np.full(len(T), 0.02)
    return z_err, t_err


# --------------------------------------------------------------------------- #
# EDI writer (standard SEG layout, variances in the .VAR blocks)
# --------------------------------------------------------------------------- #
def dms(value):
    sign = '-' if value < 0 else ''
    value = abs(value)
    d = int(value)
    m = int((value - d) * 60)
    s = (value - d - m / 60) * 3600
    return '{}{:d}:{:02d}:{:05.2f}'.format(sign, d, m, s)


def write_block(f, name, values, fmt='{:>16.7E}'):
    f.write('>{} //{}\n'.format(name, len(values)))
    for i, v in enumerate(values):
        f.write(fmt.format(v))
        if (i + 1) % 6 == 0:
            f.write('\n')
    f.write('\n\n')


def write_edi(path, site, periods, Z, T, z_err, t_err, with_tipper=True):
    lat, lon = dms(site['lat']), dms(site['lon'])
    zf, ef = Z * FIELD_UNITS, z_err * FIELD_UNITS
    with open(path, 'w', newline='\n') as f:
        f.write('>HEAD\nDATAID="{}"\nACQBY="synthetic"\nFILEBY="gateway_mt make_synthetic.py"\n'
                'FILEDATE=01/01/26\nLAT={}\nLONG={}\nELEV={:.1f}\nSTDVERS="SEG 1.0"\n'
                'EMPTY=1.0E+32\n\n'.format(site['name'], lat, lon, site['elev']))
        f.write('>INFO\nMAXINFO=999\nSynthetic 1-D-based responses for gateway_mt tests.\n\n')
        f.write('>=DEFINEMEAS\nMAXCHAN=7\nMAXRUN=999\nMAXMEAS=9999\nUNITS=M\nREFTYPE=CART\n'
                'REFLAT={}\nREFLONG={}\nREFELEV={:.1f}\n\n'.format(lat, lon, site['elev']))
        f.write('>=MTSECT\nSECTID="{}"\nNFREQ={}\n\n'.format(site['name'], len(periods)))
        write_block(f, 'FREQ', 1 / periods)
        write_block(f, 'ZROT', np.zeros(len(periods)), fmt='{:>10.2f}')
        for (i, j), comp in zip(((0, 0), (0, 1), (1, 0), (1, 1)), ('ZXX', 'ZXY', 'ZYX', 'ZYY')):
            write_block(f, comp + 'R ROT=ZROT', zf[:, i, j].real)
            write_block(f, comp + 'I ROT=ZROT', zf[:, i, j].imag)
            write_block(f, comp + '.VAR ROT=ZROT', ef ** 2)
        if with_tipper:
            write_block(f, 'TROT.EXP', np.zeros(len(periods)), fmt='{:>10.2f}')
            for k, comp in enumerate(('TX', 'TY')):
                write_block(f, comp + 'R.EXP', T[:, k].real)
                write_block(f, comp + 'I.EXP', T[:, k].imag)
                write_block(f, comp + 'VAR.EXP', t_err ** 2)
        f.write('>END\n')


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #
def build_model(sites):
    """ModEM model: layered background plus a 3-D conductor under the anomaly."""
    core = [2000.0] * 16
    pads = [2000.0 * 1.5 ** k for k in range(1, 9)]
    cs = pads[::-1] + core + pads
    zcs = [50.0 * 1.2 ** k for k in range(35)]

    model = DS.Model()
    model.generate_dummy_model()
    model.xCS, model.yCS, model.zCS = list(cs), list(cs), zcs
    xc = (np.array(model.dx[:-1]) + np.array(model.dx[1:])) / 2 / 1e3
    yc = (np.array(model.dy[:-1]) + np.array(model.dy[1:])) / 2 / 1e3
    zc = (np.array(model.dz[:-1]) + np.array(model.dz[1:])) / 2
    vals = np.zeros((model.nx, model.ny, model.nz))
    for i, x in enumerate(xc):
        for j, y in enumerate(yc):
            rho, depths = layered_model(x, y, conductor_rho=50)
            layer = np.clip(np.searchsorted(depths, zc, side='right') - 1, 0, len(rho) - 1)
            vals[i, j, :] = rho[layer]
            w = anomaly_weight(x, y)
            top = depths[2]
            in_cond = (zc >= top) & (zc <= top + 5e3)
            vals[i, j, in_cond] = 50 * (1 - w) + 3 * w
    model.vals = vals
    model.rho_x = vals
    return model


# --------------------------------------------------------------------------- #
def main():
    rng = np.random.default_rng(SEED)
    os.makedirs(EDI_DIR, exist_ok=True)
    sites = site_layout(rng)

    clean = {}
    for site in sites:
        periods = PERIODS[:-6] if site['name'] in SHORT_SITES else PERIODS
        Z, T = clean_response(site, periods)
        z_err, t_err = errors_for(Z, T)
        noise_z = (rng.normal(size=Z.shape) + 1j * rng.normal(size=Z.shape)) * z_err[:, None, None]
        noise_t = (rng.normal(size=T.shape) + 1j * rng.normal(size=T.shape)) * t_err[:, None]
        write_edi(os.path.join(EDI_DIR, site['name'] + '.edi'), site, periods,
                  Z + noise_z, T + noise_t, z_err, t_err,
                  with_tipper=site['name'] not in NO_TIPPER_SITES)
        clean[site['name']] = site

    list_file = os.path.join(EDI_DIR, 'synthetic.lst')
    with open(list_file, 'w', newline='\n') as f:
        f.write('{}\n'.format(len(sites)))
        f.write('\n'.join(s['name'] for s in sites) + '\n')

    # Inversion-format data files, built with gateway_mt itself from the EDIs
    ds = DS.Dataset(listfile=list_file, datpath=EDI_DIR)
    ds.get_data_from_raw(periods=list(PERIODS[DATA_PERIOD_IDX]))
    ds.reset_dummy_periods()      # flag data periods outside a site's recorded band
    data = ds.data
    data.UTM_zone = '16U'
    data.origin = ds.raw_data.origin   # the MARE2DEM writer needs it set
    data.write(outfile=os.path.join(OUT, 'synthetic_modem.dat'), file_format='modem')
    data.write(outfile=os.path.join(OUT, 'synthetic_ws.data'), file_format='wsinv3dmt')
    data.write(outfile=os.path.join(OUT, 'synthetic_gofem.gdat'), file_format='gofem')
    data.write(outfile=os.path.join(OUT, 'synthetic_mare2dem.emdata'), file_format='mare2dem')

    # "Response": the noise-free values at the data periods, in the data's frame
    response = deepcopy(data)
    for name in response.site_names:
        Z, T = clean_response(clean[name], response.periods)
        site = response.sites[name]
        for (i, j), comp in zip(((0, 0), (0, 1), (1, 0), (1, 1)), ('ZXX', 'ZXY', 'ZYX', 'ZYY')):
            site.data[comp + 'R'] = Z[:, i, j].real
            site.data[comp + 'I'] = -Z[:, i, j].imag
        for k, comp in enumerate(('TZX', 'TZY')):
            if comp + 'R' in site.data:
                site.data[comp + 'R'] = T[:, k].real
                site.data[comp + 'I'] = -T[:, k].imag
    response.write(outfile=os.path.join(OUT, 'synthetic_modem_resp.dat'), file_format='modem')

    with open(os.path.join(OUT, 'receivers_projected.csv'), 'w', newline='\n') as f:
        for name, (x, y) in zip(data.site_names, data.locations):
            f.write('Dipole {} 1 {:.6e} {:.6e} {:.6e}\n'.format(
                name, x, y, -data.sites[name].locations['elev']))

    build_model(sites).write(os.path.join(OUT, 'synthetic.model'))

    with open(os.path.join(OUT, 'synthetic.pymt'), 'w', newline='\n') as f:
        f.write('% synthetic\n'
                'list edis/synthetic.lst\n'
                'data synthetic_modem.dat\n'
                'resp synthetic_modem_resp.dat\n'
                'model synthetic.model\n')
    print('Wrote synthetic dataset to {}'.format(OUT))


if __name__ == '__main__':
    main()
