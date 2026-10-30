"""Plotting layer (gplot): sounding-curve plots and map views draw what the data says."""
import numpy as np
import pytest
from matplotlib.colors import Colormap
from matplotlib.figure import Figure

import gateway_mt.gplot as gplot
import gateway_mt.utils as utils
from gateway_mt.e_colours import colourmaps

PERIOD_IDX = 8          # ~6 s


GUI_COLOUR_MENU = ('bgy', 'viridis', 'jet', 'jet_plus', 'bwr', 'greys', 'turbo', 'turbo_capped',
                   'turbo_mod', 'twilight', 'twilight_shifted', 'colorwheel', 'hot')   # GUI_common.classes.ColourMenu


@pytest.mark.parametrize('name', sorted(set(colourmaps.COLOUR_MAPS + GUI_COLOUR_MENU + ('fire',))))
def test_colour_maps_resolve(name):
    """Every named colour map, the GUIs' colour menu, and an unprefixed colorcet name can be built."""
    cmap = colourmaps.get_cmap(name, N=32)
    assert isinstance(cmap, Colormap)
    assert len(cmap(0.5)) == 4             # an RGBA colour


@pytest.fixture
def six_sites(dataset):
    return dataset.get_sites(site_names=dataset.data.site_names[:6], dTypes='all')


def _curves(ax):
    """Plotted (x, y) arrays on an axis, ignoring empty legend-only lines."""
    return [(np.asarray(l.get_xdata(), float), np.asarray(l.get_ydata(), float))
            for l in ax.get_lines() if len(l.get_xdata())]


def test_sounding_curves_plot_the_apparent_resistivity(dataset, six_sites):
    """Data Plot's rho_xy curve for a site is compute_rho of that site, against log10(period)."""
    dpm = gplot.DataPlotManager(fig=Figure())
    dpm.components = ['RhoXY']
    dpm.sites = six_sites
    dpm.plot_data()
    assert len(dpm.axes) == 6
    assert dpm.axes[0].get_yscale() == 'log'
    curves = _curves(dpm.axes[0])           # raw data, inversion data, response
    assert [len(x) for x, _ in curves] == [40, 18, 18]
    for site, (x, y) in zip((dataset.data.sites['SYN001'], dataset.response.sites['SYN001']), curves[1:]):
        rho, _, _ = utils.compute_rho(site, calc_comp='xy')
        np.testing.assert_allclose(x, np.log10(site.periods))
        np.testing.assert_allclose(y, np.real(rho))


@pytest.mark.parametrize('components', [['ZXYR', 'ZXYI', 'ZYXR', 'ZYXI'], ['ZXXR', 'ZYYR'],
                                        ['PhaXY', 'PhaYX'], ['TZXR', 'TZYR'], ['RhoDet']])
def test_component_groups_plot(six_sites, components):
    dpm = gplot.DataPlotManager(fig=Figure())
    dpm.components = components
    dpm.sites = six_sites
    dpm.plot_data()
    assert len(dpm.axes) == 6
    assert all(_curves(ax) for ax in dpm.axes)


# --------------------------------------------------------------------------- #
# Map view
# --------------------------------------------------------------------------- #
def _map_view(dataset):
    fig = Figure()
    ax = fig.add_subplot(111)
    mv = gplot.MapView(fig=fig)
    mv.window['figure'], mv.window['axes'] = fig, [ax]
    for kind in ('data', 'raw_data', 'response'):
        mv.site_data[kind] = getattr(dataset, kind)
    mv.site_names = list(dataset.data.site_names)
    mv.dataset = None
    mv.interpolant = 'linear'
    mv.use_colourbar = False
    mv.site_locations['all'] = dataset.data.locations
    return mv, ax


def test_map_view_draws_an_induction_arrow_per_tipper_station(dataset):
    mv, ax = _map_view(dataset)
    mv.plot_induction_arrows(data_type=['data'], normalize=False, period_idx=PERIOD_IDX,
                             arrow_type=['R'], arrow_convention='parkinson')
    arrows = [c for c in ax.collections if type(c).__name__ == 'Quiver']
    assert len(arrows) == 1
    assert len(arrows[0].get_offsets()) == 29            # SYN030 has no tipper


def test_map_view_draws_a_phase_tensor_ellipse_per_station(dataset):
    mv, ax = _map_view(dataset)
    mv.plot_phase_tensor(data_type=['data'], normalize=True, fill_param='beta',
                         period_idx=PERIOD_IDX, cb_slot=0)
    assert len(ax.patches) == 30


@pytest.mark.parametrize('fill', ['rhoxy', 'phaxy', 'beta', 'tip'])
def test_plan_pseudosection(dataset, fill):
    mv, ax = _map_view(dataset)
    mv.plan_pseudosection(data_type=['data'], fill_param=fill, n_interp=30,
                          period_idx=PERIOD_IDX, cb_slot=0)
    # (matplotlib < 3.9 stores only the filled cells of a pcolor grid, so don't rely on its size)
    grid = np.ma.masked_invalid(np.asarray(ax.collections[0].get_array(), dtype=float))
    assert 0 < grid.count() <= 30 * 30
