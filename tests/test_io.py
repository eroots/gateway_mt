"""File input / output: EDIs, site lists, inversion data formats and models.

Round trips write a file with gateway_mt and read it back; cross-format tests
check that the same data read from different formats agrees; and the
synthetic response is checked against the analytic values it was built from.
"""
import importlib.util
import os
import shutil

import numpy as np
import pytest

import gateway_mt.data_structures as DS
import gateway_mt.IO as IO

from conftest import quiet, synthetic

LIST = synthetic('edis', 'synthetic.lst')


def _load_generator():
    """Import tests/data/make_synthetic.py as a module."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'make_synthetic.py')
    spec = importlib.util.spec_from_file_location('make_synthetic', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read(path, listfile='', file_format=None):
    with quiet():
        return DS.Data(datafile=path, listfile=listfile, file_format=file_format)


def assert_same_data(a, b, rtol=1e-6, atol=1e-10, b_errors='used_error', components=None):
    """Two Data objects hold the same sites, periods, values and (unflagged) errors.

    ``b_errors`` names the error attribute of ``b`` to compare with ``a.used_error``;
    ``components`` limits the comparison (default: all of ``a``'s).
    """
    assert list(a.site_names) == list(b.site_names)
    np.testing.assert_allclose(a.periods, b.periods, rtol=1e-4)
    for name in a.site_names:
        sa, sb = a.sites[name], b.sites[name]
        for comp in components or a.components:
            used = sa.used_error[comp] != sa.REMOVE_FLAG
            np.testing.assert_allclose(sb.data[comp][used], sa.data[comp][used], rtol=rtol, atol=atol,
                                       err_msg='{} {}'.format(name, comp))
            np.testing.assert_allclose(getattr(sb, b_errors)[comp][used], sa.used_error[comp][used],
                                       rtol=rtol, atol=atol, err_msg='{} {} errors'.format(name, comp))


# --------------------------------------------------------------------------- #
# Site lists and EDI files
# --------------------------------------------------------------------------- #
def test_read_site_list():
    names = IO.read_sites(LIST)
    assert len(names) == 30
    assert names[0] == 'SYN001' and names[-1] == 'SYN030'


def test_edis_read_as_raw_data():
    with quiet():
        raw = DS.RawData(listfile=LIST)
    assert raw.site_names == IO.read_sites(LIST)
    assert raw.sites['SYN001'].NP == 40
    assert raw.sites['SYN029'].NP == 34                         # AMT-style station
    assert 'TZXR' in raw.sites['SYN001'].components
    assert 'TZXR' not in raw.sites['SYN030'].components         # no tipper recorded


def test_edi_station_locations_match_the_survey_layout():
    """Lat / long parsed from the EDI DEFINEMEAS blocks equal the generator's station positions."""
    gen = _load_generator()
    layout = {s['name']: s for s in gen.site_layout(np.random.default_rng(gen.SEED))}
    with quiet():
        raw = DS.RawData(listfile=LIST)
    for name, site in raw.sites.items():
        assert site.locations['Lat'] == pytest.approx(layout[name]['lat'], abs=1e-5)
        assert site.locations['Long'] == pytest.approx(layout[name]['lon'], abs=1e-5)
        assert site.locations['elev'] == pytest.approx(layout[name]['elev'], abs=0.06)


def _edi_round_trip(tmp_path):
    with quiet():
        site = DS.RawData(listfile=LIST).sites['SYN005']
        IO.write_edi(site, str(tmp_path / 'SYN005.edi'))
        back = DS.RawData(listfile=str(tmp_path / 'SYN005.edi')).sites['SYN005']
    return site, back


def test_edi_write_read_round_trip_values(tmp_path):
    site, back = _edi_round_trip(tmp_path)
    np.testing.assert_allclose(back.periods, site.periods, rtol=1e-4)
    for comp in site.components:
        np.testing.assert_allclose(back.data[comp], site.data[comp], rtol=1e-6, atol=1e-12, err_msg=comp)


def test_edi_write_read_round_trip_errors(tmp_path):
    """The .VAR blocks hold variances, so errors survive a write / read."""
    site, back = _edi_round_trip(tmp_path)
    for comp in site.components:
        np.testing.assert_allclose(back.errors[comp], site.errors[comp], rtol=1e-6, err_msg=comp)


# --------------------------------------------------------------------------- #
# Inversion data formats
# --------------------------------------------------------------------------- #
@pytest.fixture(scope='module')
def modem_data():
    return _read(synthetic('synthetic_modem.dat'))


def test_modem_data_file(modem_data):
    assert modem_data.NS == 30 and modem_data.NP == 18
    assert modem_data.inv_type == 5                 # full impedance + tipper
    assert len(modem_data.components) == 12


@pytest.mark.parametrize('fmt, ext, needs_list', [('modem', '.dat', False),
                                                  ('wsinv3dmt', '.data', True),
                                                  ('gofem', '.gdat', False)])
def test_data_write_read_round_trip(tmp_path, modem_data, fmt, ext, needs_list):
    out = str(tmp_path / ('round_trip' + ext))
    if fmt == 'gofem':                               # GoFEM station positions live in a receiver file
        shutil.copy(synthetic('receivers_projected.csv'), tmp_path)
    with quiet():
        modem_data.write(outfile=out, file_format=fmt)
    back = _read(out, listfile=LIST if needs_list else '', file_format=fmt)
    assert_same_data(modem_data, back, rtol=1e-5)
    np.testing.assert_allclose(back.locations, modem_data.locations, atol=0.01)


@pytest.fixture
def mare2dem_round_trip(tmp_path, modem_data):
    out = str(tmp_path / 'round_trip.emdata')
    with quiet():
        modem_data.write(outfile=out, file_format='mare2dem')
        return _read(out, file_format='mare2dem')


def test_mare2dem_header_round_trip(mare2dem_round_trip, modem_data):
    """The UTM zone and origin survive a write / read (7-digit northings included)."""
    assert mare2dem_round_trip.UTM_zone == modem_data.UTM_zone
    np.testing.assert_allclose(mare2dem_round_trip.origin, modem_data.origin)


def test_mare2dem_write_read_round_trip(mare2dem_round_trip, modem_data):
    """MARE2DEM is 2-D, so only TZY is kept; everything else reads back (to its 6 significant figures)."""
    back = mare2dem_round_trip
    assert set(modem_data.components) - set(back.components) == {'TZXR', 'TZXI'}
    # 6 significant figures, and the error floors are recomputed from the rounded values
    assert_same_data(modem_data, back, rtol=5e-5, components=back.components)


def test_mare2dem_written_from_edis_alone(tmp_path):
    """Data made from EDIs has no UTM zone / origin; the writer uses 0 instead of failing."""
    out = str(tmp_path / 'from_edis.emdata')
    with quiet():
        data = DS.Dataset(listfile=synthetic('edis', 'synthetic.lst')).data
        data.write(outfile=out, file_format='mare2dem')
        back = _read(out, file_format='mare2dem')
    assert (back.NS, back.NP) == (data.NS, data.NP)


def test_gofem_file_agrees_with_modem(modem_data):
    """The same dataset stored as GoFEM reads back identical to the ModEM file."""
    gofem = _read(synthetic('synthetic_gofem.gdat'), file_format='gofem')
    assert_same_data(modem_data, gofem, rtol=1e-5)


def test_wsinv3dmt_file_agrees_with_modem(modem_data):
    """The same dataset stored as WSINV3DMT agrees with the ModEM file.

    WSINV3DMT keeps an error map (multipliers such as the x10 applied to
    outliers) next to the errors, and the reader applies it. ModEM files
    store only final errors, and the synthetic data was written without
    applying its map, so errors are compared where the map is 1.
    """
    ws = _read(synthetic('synthetic_ws.data'), listfile=LIST, file_format='wsinv3dmt')
    assert list(ws.site_names) == list(modem_data.site_names)
    n_mapped = 0
    for name in ws.site_names:
        a, b = modem_data.sites[name], ws.sites[name]
        for comp in ws.components:
            used = a.used_error[comp] != a.REMOVE_FLAG
            np.testing.assert_allclose(b.data[comp][used], a.data[comp][used], rtol=1e-5, atol=1e-10)
            unmapped = used & (b.errmap[comp] == 1)
            np.testing.assert_allclose(b.used_error[comp][unmapped], a.used_error[comp][unmapped], rtol=1e-5)
            n_mapped += np.sum(used & (b.errmap[comp] != 1))
    assert n_mapped > 0          # the file does carry an error map


def test_response_file_holds_the_analytic_responses():
    """synthetic_modem_resp.dat contains the noise-free responses make_synthetic.py computed."""
    gen = _load_generator()
    layout = {s['name']: s for s in gen.site_layout(np.random.default_rng(gen.SEED))}
    resp = _read(synthetic('synthetic_modem_resp.dat'))
    for name in ('SYN001', 'SYN015', 'SYN030'):
        site = resp.sites[name]
        Z, T = gen.clean_response(layout[name], site.periods)
        for comp, (i, j) in (('ZXY', (0, 1)), ('ZYX', (1, 0)), ('ZXX', (0, 0))):
            np.testing.assert_allclose(site.data[comp + 'R'], Z[:, i, j].real, rtol=1e-5, atol=1e-9)
            np.testing.assert_allclose(site.data[comp + 'I'], -Z[:, i, j].imag, rtol=1e-5, atol=1e-9)


# --------------------------------------------------------------------------- #
# FEMTIC
# --------------------------------------------------------------------------- #
def read_observe(path):
    """{block: {id: (x_km, y_km, rows)}} from a FEMTIC observe.dat; rows are
    (frequency, values..., errors...) per frequency."""
    n_values = {'MT': 16, 'PT': 8, 'VTF': 8}
    with open(path) as f:
        tok = f.read().split()
    out, i = {}, 0
    while tok[i] != 'END':
        block, n = tok[i], int(tok[i + 1])
        i += 2
        out[block] = {}
        for _ in range(n):
            sid, x, y, nf = int(tok[i]), float(tok[i + 2]), float(tok[i + 3]), int(tok[i + 4])
            i += 5
            width = 1 + n_values[block]
            out[block][sid] = (x, y, np.array(tok[i:i + nf * width], dtype=float).reshape(nf, width))
            i += nf * width
    return out


def _femtic_row(rows, period):
    hit = np.flatnonzero(np.isclose(rows[:, 0], 1 / period, rtol=1e-6))
    return rows[hit[0]] if hit.size else None


def test_femtic_file_holds_the_data_as_stored(tmp_path, modem_data):
    """Impedances and tippers go out as stored (ohm, exp(-i omega t): FEMTIC's convention),
    flagged points with a negative error (FEMTIC skips them), stations at the receiver-file
    positions in km, and the tipper stations with their own ids from 101."""
    out = str(tmp_path / 'observe.dat')
    receivers = IO.read_receivers(synthetic('receivers_projected.csv'))
    counts = IO.write_femtic(modem_data, out, locations=receivers)
    assert counts == {'MT': 30, 'VTF': 29}                       # SYN030 has no tipper
    obs = read_observe(out)
    assert sorted(obs['MT']) == list(range(1, 31)) and sorted(obs['VTF']) == list(range(101, 130))
    flagged = 0
    comps = [c + p for c in ('ZXX', 'ZXY', 'ZYX', 'ZYY') for p in ('R', 'I')]
    for ii, name in enumerate(modem_data.site_names):
        site = modem_data.sites[name]
        x, y, rows = obs['MT'][ii + 1]
        np.testing.assert_allclose((x, y), np.array(receivers[name][:2]) / 1000, atol=1e-6)
        for jj, period in enumerate(modem_data.periods):
            row = _femtic_row(rows, period)
            errors = np.array([site.used_error[c][jj] for c in comps])
            if row is None:                                    # dropped: nothing usable at this period
                assert (errors == modem_data.REMOVE_FLAG).all()
                flagged += 1
                continue
            for k, comp in enumerate(comps):
                if errors[k] == modem_data.REMOVE_FLAG:
                    assert row[9 + k] == -1
                else:
                    np.testing.assert_allclose(row[1 + k], site.data[comp][jj], rtol=1e-6, atol=1e-12)
                    np.testing.assert_allclose(row[9 + k], abs(errors[k]), rtol=1e-6)
    assert flagged                                             # SYN029 / SYN030 stop at ~470 s
    _x, _y, rows = obs['VTF'][101]
    row = _femtic_row(rows, modem_data.periods[0])
    site = modem_data.sites['SYN001']
    np.testing.assert_allclose(row[1:5], [site.data[c][0] for c in ('TZXR', 'TZXI', 'TZYR', 'TZYI')], rtol=1e-6)
    with open(tmp_path / 'observe_sites.csv') as f:
        lines = f.read().splitlines()
    assert lines[1].startswith('1,MT,SYN001,') and lines[31].startswith('101,VTF,SYN001,')


def test_femtic_phase_tensor_is_femtics_own(tmp_path, modem_data):
    """FEMTIC's phase tensor is X^-1 Y of its own impedance (ObservedDataStationPT.cpp), i.e. of
    the stored one; the stored phase tensor is its negative."""
    import copy
    data = copy.deepcopy(modem_data)
    with quiet():
        data.set_inversion_type(6)
    out = str(tmp_path / 'observe.dat')
    assert IO.write_femtic(data, out) == {'PT': 30}
    _x, _y, rows = read_observe(out)['PT'][1]
    site = data.sites['SYN001']
    for jj, period in enumerate(data.periods):
        z = np.array([[site.data['ZXXR'][jj] + 1j * site.data['ZXXI'][jj],
                       site.data['ZXYR'][jj] + 1j * site.data['ZXYI'][jj]],
                      [site.data['ZYXR'][jj] + 1j * site.data['ZYXI'][jj],
                       site.data['ZYYR'][jj] + 1j * site.data['ZYYI'][jj]]])
        row = _femtic_row(rows, period)
        np.testing.assert_allclose(row[1:5], np.linalg.solve(z.real, z.imag).ravel(), rtol=1e-5, atol=1e-6)
        assert (row[5:] > 0).all()


def test_femtic_through_write_data(tmp_path, modem_data):
    """file_format='femtic' uses the data's own locations."""
    out = str(tmp_path / 'observe.dat')
    with quiet():
        modem_data.write(outfile=out, file_format='femtic')
    x, y, _rows = read_observe(out)['MT'][1]
    np.testing.assert_allclose((x, y), modem_data.locations[0] / 1000, atol=1e-6)


def test_femtic_needs_every_site_located(tmp_path, modem_data):
    with pytest.raises(ValueError, match='SYN002'):
        IO.write_femtic(modem_data, str(tmp_path / 'observe.dat'), locations={'SYN001': (0, 0, 0)})


def test_write_list(tmp_path, modem_data):
    out = str(tmp_path / 'sites.lst')
    IO.write_list(modem_data, out)
    assert IO.read_sites(out) == list(modem_data.site_names)


# --------------------------------------------------------------------------- #
# Models
# --------------------------------------------------------------------------- #
@pytest.fixture(scope='module')
def model():
    with quiet():
        return DS.Model(synthetic('synthetic.model'))


def test_model_file(model):
    assert (model.nx, model.ny, model.nz) == (32, 32, 35)
    assert model.vals.shape == (32, 32, 35)
    assert len(model.dx) == model.nx + 1 and np.all(np.diff(model.dx) > 0)
    assert model.dz[0] == 0 and np.all(np.diff(model.dz) > 0)
    assert 1 < model.vals.min() and model.vals.max() < 1e4


def test_model_write_read_round_trip(tmp_path, model):
    out = str(tmp_path / 'round_trip.model')
    with quiet():
        model.write(out)
        back = DS.Model(out)
    for edges in ('dx', 'dy', 'dz'):
        np.testing.assert_allclose(getattr(back, edges), getattr(model, edges), rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(back.vals, model.vals, rtol=1e-5)
