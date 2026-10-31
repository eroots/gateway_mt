"""Model tools: grid conversion, transect sections and the ParaView export (no Qt needed)."""
import os

import numpy as np
import pytest

import gateway_mt.data_structures as DS
from gateway_mt import model_tools

from conftest import quiet, synthetic


@pytest.fixture
def model():
    with quiet():
        return DS.Model(synthetic('synthetic.model'))


def test_converts_a_model_to_a_vtk_grid(model):
    pytest.importorskip('pyvista')
    grid = model_tools.model_to_rectgrid(model)
    assert grid.n_cells == model.nx * model.ny * model.nz
    resistivity = grid.cell_data['Resistivity']               # log10(ohm-m)
    assert resistivity.min() == pytest.approx(np.log10(model.vals.min()))
    assert resistivity.max() == pytest.approx(np.log10(model.vals.max()))


def test_transect_path_follows_every_segment():
    """North-south segments used to collapse to a point; every segment is now sampled along its length."""
    qx, qy = model_tools.sample_path([0, 0, 10], [0, 10, 10], 20)      # north, then east
    assert qy.max() == 10 and len(np.unique(qy[qx == 0])) > 5          # the N-S leg has many points
    assert (qx[0], qy[0], qx[-1], qy[-1]) == (0, 0, 10, 10)
    steps = np.hypot(np.diff(qx), np.diff(qy))
    assert np.all(steps > 0)                                            # vertices not repeated


def test_stations_project_onto_the_transect():
    along, across = model_tools.project_onto_path([2, 11, 5], [1, 4, -3], [0, 10, 10], [0, 0, 10])
    np.testing.assert_allclose(along, [2, 14, 5])
    np.testing.assert_allclose(across, [1, 1, 3])
    assert model_tools.nearest_point(9, 9, [0, 10, 3], [0, 10, 3]) == (1, pytest.approx(np.sqrt(2)))


def test_section_samples_the_model_along_a_path(model):
    qx, qy, depths, values = model_tools.sample_section(model, [-10000, 10000], [0, 0], n_points=50)
    assert values.shape == (len(qx), model.nz + 1)
    np.testing.assert_allclose(depths, model.dz)
    assert np.log10(model.vals.min()) - 1e-9 <= values.min() and values.max() <= np.log10(model.vals.max()) + 1e-9
    # A ratio of two axes works too (Model Viewer's anisotropy views)
    values = model_tools.sample_section(model, [-10000, 10000], [0, 0], n_points=50, rho_axis='rho_x / vals')[3]
    np.testing.assert_allclose(values, 0, atol=1e-12)


def test_exports_to_paraview(tmp_path, model):
    """The model (with resolution) and stations are written in metres, as by the to_vtk command."""
    with quiet():
        dataset = DS.Dataset(listfile=synthetic('edis', 'synthetic.lst'), datafile=synthetic('synthetic_modem.dat'))
        west_edge, centre_east = model.dy[0], dataset.raw_data.origin[0]   # in m
        model.spatial_units = 'km'                                       # as shown in Model Viewer
        dataset.spatial_units = 'km'
        written = model_tools.export_to_paraview(str(tmp_path / 'out'), model, data=dataset.data,
                                                 raw_data=dataset.raw_data,
                                                 resolution=np.linspace(0.1, 1, model.vals.size).reshape(model.vals.shape))
    assert [os.path.basename(f) for f in written] == ['out_model.vtk', 'out_sites.vtk']
    text = (tmp_path / 'out_model.vtk').read_text()
    assert 'SCALARS Raw_Resolution' in text
    east = float(text.split('X_COORDINATES')[1].splitlines()[1].split()[0])
    assert east == pytest.approx(centre_east + west_edge, rel=1e-6)         # UTM easting, in m


def test_writes_a_section_to_csv(tmp_path):
    model_tools.write_section_csv(str(tmp_path / 'section.csv'), np.array([0., 3.]), np.array([0., 4.]),
                                  np.array([0., 1.]), np.array([[1., 2.], [3., 4.]]), 'km')
    rows = (tmp_path / 'section.csv').read_text().splitlines()
    assert rows[0].startswith('#') and rows[-1] == '3.0000,4.0000,5.0000,1.0000,4.00000'
