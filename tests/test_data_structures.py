"""The core data structures: Dataset, Data, RawData, Model."""
from copy import deepcopy

import numpy as np
import pytest

import gateway_mt.data_structures as DS

from conftest import quiet, synthetic


def test_dataset_loads_all_parts(dataset):
    assert dataset.data.NS == dataset.response.NS == dataset.raw_data.NS == 30
    assert dataset.data.NP == 18
    assert dataset.site_names == dataset.data.site_names
    assert (dataset.model.nx, dataset.model.ny, dataset.model.nz) == (32, 32, 35)


def test_rms_is_about_one(dataset):
    """The synthetic data has noise at the size of its errors, so data vs. response RMS is ~1."""
    assert dataset.rms['Total'] == pytest.approx(1.0, abs=0.1)
    for name, station in dataset.rms['Station'].items():
        assert 0.6 < station['Total'] < 1.4, name
    for comp in dataset.data.components:
        assert 0.7 < dataset.rms['Component'][comp] < 1.3, comp


def test_tipper_only_data_loads_with_a_response(tmp_path):
    """SYN030 has no tipper, so in a tipper-only dataset all its points are flagged; its RMS is masked."""
    with quiet():
        for name, source in (('tip', 'synthetic_modem.dat'), ('tip_resp', 'synthetic_modem_resp.dat')):
            data = DS.Data(datafile=synthetic(source))
            data.inv_type = 3                               # tipper only
            data.write(str(tmp_path / name), file_format='modem')
        tipper = DS.Dataset(datafile=str(tmp_path / 'tip.dat'), responsefile=str(tmp_path / 'tip_resp.dat'))
    assert tipper.rms['Total'] == pytest.approx(1.0, abs=0.15)
    assert tipper.rms['Station']['SYN030']['Total'] is np.ma.masked
    assert 0.5 < tipper.rms['Station']['SYN001']['Total'] < 1.5


def test_periods_outside_a_stations_band_are_flagged(dataset):
    """SYN029 / SYN030 stop at ~470 s, so the two longest data periods are flagged; SYN030 has no tipper."""
    flag = dataset.data.REMOVE_FLAG
    long_periods = dataset.data.periods > 500
    for name in ('SYN029', 'SYN030'):
        errors = dataset.data.sites[name].used_error['ZXYR']
        assert np.all(errors[long_periods] == flag)
        assert not np.any(errors[~long_periods] == flag)
    assert np.all(dataset.data.sites['SYN030'].used_error['TZXR'] == flag)
    assert not np.any(dataset.data.sites['SYN001'].used_error['ZXYR'] == flag)


def test_raw_data_resampled_to_the_data_periods_reproduces_the_data_file(dataset):
    """get_data_from_raw on the EDIs gives back the values stored in the inversion data file."""
    with quiet():
        ds = DS.Dataset(listfile=synthetic('edis', 'synthetic.lst'))
        ds.get_data_from_raw(periods=list(dataset.data.periods))
    for name in dataset.data.site_names:
        stored, resampled = dataset.data.sites[name], ds.data.sites[name]
        for comp in dataset.data.components:
            used = stored.used_error[comp] != stored.REMOVE_FLAG
            np.testing.assert_allclose(resampled.data[comp][used], stored.data[comp][used],
                                       rtol=1e-6, atol=1e-12, err_msg='{} {}'.format(name, comp))


def test_data_locations_are_centred_on_the_survey(dataset):
    locs = dataset.data.locations
    assert locs.shape == (30, 2)
    centre = (locs.max(axis=0) + locs.min(axis=0)) / 2
    np.testing.assert_allclose(centre, 0, atol=1.0)
    assert np.ptp(locs[:, 0]) == pytest.approx(21e3, rel=0.1)   # 5 rows of stations, 5 km apart


def test_rotating_a_dataset_is_reversible(fresh_dataset):
    original = deepcopy({n: s.data for n, s in fresh_dataset.data.sites.items()})
    fresh_dataset.rotate_sites(azi=30)
    assert fresh_dataset.azimuth == 30
    assert not np.allclose(fresh_dataset.data.sites['SYN001'].data['ZXXR'], original['SYN001']['ZXXR'])
    fresh_dataset.rotate_sites(azi=0)
    for name, data in original.items():
        for comp in data:
            np.testing.assert_allclose(fresh_dataset.data.sites[name].data[comp], data[comp], atol=1e-12)


def test_remove_and_add_a_site(fresh_dataset, tmp_path):
    data = fresh_dataset.data
    removed = deepcopy(data.sites['SYN010'])
    data.remove_sites(sites=['SYN010'])
    assert data.NS == 29 and 'SYN010' not in data.site_names
    assert data.locations.shape[0] == 29
    data.add_site(removed)
    assert data.NS == 30 and 'SYN010' in data.site_names
    assert data.locations.shape[0] == 30
    with quiet():
        data.write(str(tmp_path / 'readded'), file_format='modem')
        back = DS.Data(datafile=str(tmp_path / 'readded.dat'))
    assert back.site_names == data.site_names
    np.testing.assert_allclose(back.locations, data.locations, atol=1e-6)


def test_phase_tensors_exist_for_every_period(dataset):
    site = dataset.raw_data.sites['SYN001']
    assert len(site.phase_tensors) == site.NP
    for pt in site.phase_tensors:
        assert pt.phi_max >= pt.phi_min


# --------------------------------------------------------------------------- #
# Model
# --------------------------------------------------------------------------- #
@pytest.fixture
def model():
    with quiet():
        return DS.Model(synthetic('synthetic.model'))


def test_model_mesh_edges_and_cell_sizes_agree(model):
    np.testing.assert_allclose(np.diff(model.dx), model.xCS)
    np.testing.assert_allclose(np.diff(model.dy), model.yCS)
    np.testing.assert_allclose(np.diff(model.dz), model.zCS)
    assert model.dx[0] == pytest.approx(-model.dx[-1])           # mesh centred on the survey


def test_deleting_and_inserting_mesh_lines(model):
    nx, vals = model.nx, model.vals.copy()
    model.dx_delete(0)
    assert model.nx == nx - 1 and model.vals.shape[0] == nx - 1
    np.testing.assert_allclose(model.vals, vals[1:])
    new_line = (model.dx[10] + model.dx[11]) / 2
    model.dx_insert(11, new_line)
    assert model.nx == nx and new_line in model.dx


def test_model_conductor_sits_under_the_anomaly(model):
    """The synthetic model's most conductive cells lie close to the tipper anomaly (2 km N, 3 km W)."""
    centres_x = (np.array(model.dx[:-1]) + np.array(model.dx[1:])) / 2
    centres_y = (np.array(model.dy[:-1]) + np.array(model.dy[1:])) / 2
    i, j, k = np.unravel_index(np.argmin(model.vals), model.vals.shape)
    assert abs(centres_x[i] - 2000) < 3000 and abs(centres_y[j] + 3000) < 3000
    assert 5e3 < (model.dz[k] + model.dz[k + 1]) / 2 < 25e3


# --------------------------------------------------------------------------- #
# Period tolerance
# --------------------------------------------------------------------------- #
def test_borrowed_periods_use_one_tolerance(dataset):
    """A period 5% from the nearest real one is borrowed at any period; 30% away is flagged."""
    site = dataset.raw_data.sites['SYN001']
    short, long = site.periods[2], site.periods[-3]               # one below 1 s, one above
    assert short < 1 < long
    resampled = site.get_data(periods=[short * 1.05, long * 1.05, short * 1.3, long * 1.3],
                              components=['ZXYR', 'ZXYI'])
    flagged = resampled.errmap['ZXYR'] == DS.Site.NO_PERIOD_FLAG
    assert list(flagged) == [False, True, False, True]            # sorted: short x1.05, x1.3, long x1.05, x1.3


def test_old_tolerance_names_still_work(fresh_dataset):
    raw = fresh_dataset.raw_data
    raw.low_tol = 0.05
    assert raw.tol == raw.high_tol == 0.05
    raw.remove_high_tol = 0.3
    assert raw.remove_tol == raw.remove_low_tol == 0.3
    with quiet():
        raw.get_data(sites=['SYN001'], periods=list(raw.sites['SYN001'].periods[:3]), lTol=0.02, hTol=0.15)
    assert raw.tol == 0.15                                         # the larger of the old pair
