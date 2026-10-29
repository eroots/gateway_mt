"""Matplotlib rendering layer for gateway_mt.

``gplot`` owns no data of its own.  It takes the in-memory objects built
by :mod:`gateway_mt.data_structures` (``Data``, ``RawData``, ``Response``,
``Model``, ``Site``, ``PhaseTensor``, ``CART``) plus a bare
:class:`matplotlib.figure.Figure` and draws onto it, calling *down* into
:mod:`gateway_mt.utils` (``compute_rho`` / ``compute_phase`` / ``compute_bost1D``,
``normalize*``, ``remove_outliers``, ``rotate_locs`` …) and
:mod:`gateway_mt.e_colours.colourmaps` for its colour maps.  It is driven *up*
by the GUIs (:mod:`gateway_mt.DataGUI.data_plot`, :mod:`gateway_mt.ModelGUI.model_viewer`)
and by several plotting scripts.

Two rendering classes, plus two small coordinate-formatter helpers:

``format_data_coords`` / ``format_model_coords``
    Callables assigned to ``ax.format_coord`` so the interactive toolbar
    shows period / frequency (sounding axes) or the value under the
    cursor (image / pcolor axes).

:class:`DataPlotManager`
    The multi-panel *sounding-curve* figure used by ``data_plot`` -- one
    subplot per site, plotting apparent resistivity / phase / impedance /
    tipper / phase-tensor components against period.  Behaviour is
    controlled by setting plain attributes (``components``, ``scale``,
    ``errors``, ``toggles`` …) then calling :meth:`~DataPlotManager.plot_data`;
    see the class docstring for the full list.

:class:`MapView`
    The plan-view *map* figure used by ``data_plot`` and ``model_viewer``
    -- site locations, induction arrows, phase-tensor ellipses / bars,
    interpolated pseudosections, RMS maps, and model plan / cross
    sections.  Also attribute-driven; ``self.dataset`` is expected to be
    assigned by the caller (it is not created in ``__init__``).

Notes
-----
* Methods decorated with ``@utils.enforce_input(...)``
  (``plot_phase_tensor``, ``plot_phase_bar``, ``plan_pseudosection``,
  ``plot_induction_arrows``, ``link_axes``, ``tensor_*_pseudosection``)
  coerce the listed arguments to the declared types but must be called
  with **keyword** arguments -- positional args past ``self`` break the
  wrapper, and any keyword not named in the decorator is dropped.
* Redrawing is coarse: most update paths clear whole axes and re-plot
  rather than mutating existing artists.
* Rotation of raw vs. inverted data through the plotting path has not
  been rigorously verified.
"""

import numpy as np
from matplotlib.figure import Figure
from scipy.interpolate import griddata
import matplotlib.patches as patches
from matplotlib.path import Path
from scipy.interpolate import RBFInterpolator
import gateway_mt.utils as utils
from gateway_mt.e_colours import colourmaps as cm
from copy import deepcopy
try:
    import naturalneighbor as nn
    has_nn = True
except ModuleNotFoundError:
    has_nn = False
    print('Natural Neighbor not found and will not be available')


def format_data_coords(x, y):
    """Cursor readout for a sounding-curve axis.

    Assigned to ``ax.format_coord``.  The x-axis carries
    :math:`\\log_{10}` period, so both the log value and the linear
    period / frequency are shown.

    Parameters
    ----------
    x, y : float
        Cursor position in axis (data) coordinates; ``x`` is
        :math:`\\log_{10}` period in seconds.

    Returns
    -------
    str
        Two-line string with the log-scale values and the linear
        period / frequency.
    """
    freq = np.log10(1 / 10**x)
    return '\n'.join(['Log10: period={}, frequency={}, y={}',
                      'period={}, frequency={}']).format(utils.truncate(x),
                                                         utils.truncate(freq),
                                                         utils.truncate(y),
                                                         utils.truncate(10**x),
                                                         utils.truncate(10**freq))


class format_model_coords(object):
    """Cursor readout for an image / ``pcolor`` axis.

    An instance is assigned to ``ax.format_coord``; calling it with the
    cursor position returns the label string, having looked up the cell
    the cursor falls in and the value stored there.

    Parameters
    ----------
    im : matplotlib.cm.ScalarMappable
        The image / mesh whose array holds the plotted values.
    X, Y : array_like
        1-D cell-edge (or cell-centre) coordinates along the x and y
        axes.  A model passes edges, so ``im`` has one fewer value per
        axis than ``X`` / ``Y``; a pseudosection passes centres and the
        sizes match.  Both cases are handled.
    x_label, y_label : str, optional
        Axis names used in the readout (default ``'y'``).
    use_log : bool, optional
        If ``True`` (default) the stored value is treated as
        :math:`\\log_{10}` and raised to the power of 10 before display.
    data_label : str, optional
        Name of the plotted quantity (``'Resistivity'`` by default);
        also selects the unit string shown.
    """

    def __init__(self, im, X, Y, x_label='y', y_label='y', use_log=True, data_label='Resistivity'):
        self.im = im
        self.x_label = x_label
        self.y_label = y_label
        self.X = X
        self.Y = Y
        self.use_log = use_log
        self.data_label = data_label

    def __call__(self, x, y):
        for ix, xx in enumerate(self.X):
            if xx > x:
                x_idx = min(ix, len(self.X) - 1) - 1
                break
            x_idx = len(self.X) - 2
        for iy, yy in enumerate(self.Y):
            if yy > y:
                y_idx = min(iy, len(self.Y) - 1) - 1
                break
            y_idx = len(self.Y) - 2
        vals = np.array(self.im.get_array())
        # Pseudosection images pass cell centres (sizes match X / Y); model
        # images pass cell edges, so the array is one smaller on each axis.
        if vals.size == self.X.size * self.Y.size:
            nx = self.X.size
            ny = self.Y.size
        else:
            nx = self.X.size - 1
            ny = self.Y.size - 1

        vals = np.reshape(vals, (ny, nx))[y_idx, x_idx]

        if self.use_log:
            vals = 10 ** vals
        if self.data_label.lower() == 'resistivity':
            self.data_units = 'Ωm'
        elif self.data_label.lower() in ('phase', 'skew', 'azimuth'):
            self.data_units = 'Degrees'
        elif self.data_label.lower() == 'depth':
            self.data_units = 'km'
        else:
            self.data_units = ''

        return '\t'.join(['{}: {:>4.4g} {}',
                          '{}: {:>4.4g} {}\n',
                          '{}: {} {}']).format(self.x_label, utils.truncate(self.X[x_idx]), 'km',
                                                        self.y_label, utils.truncate(self.Y[y_idx]), 'km',
                                                        self.data_label,
                                                        utils.truncate(vals), self.data_units)


class DataPlotManager(object):
    """Multi-panel sounding-curve figure: one subplot per site.

    Holds a :class:`matplotlib.figure.Figure` and renders MT response
    curves onto it -- apparent resistivity, phase, impedance, tipper,
    phase-tensor components or 1-D Bostick depth, plotted against
    :math:`\\log_{10}` period, with one subplot per site.  Used by
    :mod:`gateway_mt.DataGUI.data_plot` and by the ``plot_dataset`` /
    ``plot_misfit_curves`` / ``autogen_report`` scripts.

    Typical use: create with a figure, set the configuration attributes
    below, assign ``sites``, then call :meth:`plot_data`.

    Parameters
    ----------
    fig : matplotlib.figure.Figure, optional
        Figure to draw on.  If omitted a new one is created
        (:meth:`new_figure`).

    Attributes
    ----------
    sites : dict
        ``{data_type: [Site, ...]}`` for ``data_type`` in
        ``'raw_data'``, ``'data'``, ``'response'``, ``'smoothed_data'``
        -- typically the output of ``Dataset.get_sites(..., dTypes='all')``.
        The lists are parallel: index ``i`` is the same station in each.
    site1D : list
        Site objects for the optional 1-D overlay (``'1d'`` toggle).
    components : list of str
        Components to plot, e.g. ``['ZXYR', 'ZXYI']``, ``['RHOXY']``,
        ``['PHAXY']``, ``['TZXR']``, ``['PTXY']``, ``['BOSTXY']``,
        ``['BETA']``.  The first entry drives the axis labels / units.
    scale : {'sqrt(periods)', 'periods', 'log10', 'none'}
        Period weighting applied to impedance components before plotting.
    errors : {'mapped', 'raw', 'none'}
        Which error to draw as error bars (``used_error`` vs ``errors``).
    which_errors : list of str
        Data types for which error bars are drawn at all.
    toggles : dict
        ``{data_type: bool}`` -- whether each data type is drawn
        (includes the ``'1d'`` key).
    marker : dict
        ``{data_type: spec}`` where ``spec`` is ``'o'`` (points),
        ``'oo'`` (points with an edge), ``'-'`` / ``'--'`` (lines).
    show_outliers : bool
        If ``False``, y-limits are computed from an outlier-filtered copy
        of the data (:func:`gateway_mt.utils.remove_outliers`).
    plot_flagged_data : bool
        If ``False``, points whose error equals ``Site.REMOVE_FLAG`` are
        dropped before plotting.
    wrap_phase : int
        Passed through to :func:`gateway_mt.utils.compute_phase`.
    outlier_thresh : float
        Threshold for the outlier filter used when ``show_outliers`` is
        ``False``.
    min_ylim, max_ylim : float or None
        Hard y-limit overrides (``None`` = auto).
    pt_units : {'degrees', 'radians'}
        Units for phase-tensor components.
    tiling : list of int
        ``[nrows, ncols]`` subplot grid, recomputed by :meth:`gen_tiling`.
    use_designated_colours : bool
        If ``True``, colour each component by :attr:`designated_colours`
        rather than cycling :attr:`colour`.
    designated_colours, colour, y_labels, ax_lim_dict : dict / list
        Colour and label lookup tables.
    mec, markersize, edgewidth, edge_scale, linestyle, label_fontsize
        Matplotlib styling knobs.
    link_axes_bounds : bool
        If ``True``, :meth:`link_axes` shares one set of bounds across
        every subplot.
    artist_ref : dict
        ``{data_type: [artist, ...]}`` -- the artist returned for each
        axis by the most recent :meth:`plot_site`, kept for incremental
        redraws.
    fig : matplotlib.figure.Figure
    axes : list of matplotlib.axes.Axes
        The current subplots (populated by :meth:`new_figure` /
        :meth:`redraw_axes`).
    """

    def __init__(self, fig=None):
        self.link_axes_bounds = False
        self.axis_padding = 0.1
        self.use_designated_colours = True
        self.colour = ['darkgray', 'r', 'royalblue', 'g', 'm', 'y', 'lime', 'peru']
        self.designated_colours = {'xy': 'royalblue',
                                   'yx': 'red',
                                   'xx': 'lime',
                                   'yy': 'peru',
                                   'tzx': 'royalblue',
                                   'tzy': 'red',
                                   'aav': 'green',
                                   'gav': 'magenta',
                                   'det': 'darkgray',
                                   'ssq': 'orange'}
        self.sites = {'raw_data': [], 'data': [], 'response': [], 'smoothed_data': []}
        self.toggles = {'raw_data': True, 'data': True, 'response': True, '1d': False, 'smoothed_data': False}
        self.site1D = []
        self.marker = {'raw_data': 'o', 'data': 'oo', 'response': '-', '1d': '--', 'smoothed_data': '--'}
        self.errors = 'mapped'
        self.which_errors = ['data', 'raw_data']
        self.components = ['ZXYR']
        self.linestyle = '-'
        self.scale = 'sqrt(periods)'
        self.mec = 'k'
        self.markersize = 5
        self.edgewidth = 2
        self.edge_scale = 0.25
        self.label_fontsize = 14
        self.tiling = [0, 0]
        self.show_outliers = True
        self.plot_flagged_data = True
        self.wrap_phase = 1
        self.outlier_thresh = 2
        self.min_ylim = None
        self.max_ylim = None
        self.ax_lim_dict = {'rho': [0, 5], 'phase': [0, 120], 'impedance': [-1, 1],
                            'tipper': [-1, 1], 'beta': [-10, 10], 'azimuth': [-90, 90], 'skew': [-10, 10]}
        self.artist_ref = {'raw_data': [], 'data': [], 'response': [], '1d': [], 'smoothed_data': []}
        self.y_labels = {'rho': r'$Log_{10} Rho_{a}$',
                         'zxx': 'Impedance', 'zxy': 'Impedance', 'zyx': 'Impedance', 'zyy': 'Impedance',
                         'tzx': 'Magnitude', 'tzy': 'Magnitude',
                         'bos': r'$Log_{10} Rho_{a}$',
                         'ptx': 'Phase', 'pty': 'Phase', 'phi': 'Phase', 'pha': 'Phase',
                         'bet': 'Skew',
                         'azi': 'Azimuth'}
        self.pt_units = 'degrees'
        if fig is None:
            self.new_figure()
        else:
            self.fig = fig

    @property
    def site_names(self):
        """list of str: Station names, taken from the first non-empty
        ``sites`` list."""
        site_containers = self.sites.values()
        sites = next(site for site in site_containers if site != [])
        return [site.name for site in sites]

    @property
    def num_sites(self):
        """int: Largest ``sites`` list length -- i.e. the subplot count."""
        return max([len(data) for data in self.sites.values()])

    @property
    def units(self):
        """str: Unit string for the current component / scale, for axis
        labels (impedance in :math:`\\Omega`, optionally period-weighted;
        tipper unitless; resistivity in :math:`\\Omega`-m; phase / skew /
        azimuth in degrees)."""
        if self.components[0][0].upper() == 'Z':
            units = 'Ω'
            if self.scale.lower() == 'periods':
                units = ''.join(['s*', units])
            elif self.scale.lower() == 'sqrt(periods)':
                units = ''.join([r'$\sqrt{s}$*', units])
            elif self.scale.lower() == 'log10':
                units = ''.join([r'$Log_{10}$', units])
        elif self.components[0][0].upper() == 'T':
            units = 'Unitless'
        elif self.components[0][0].upper() == 'R' or self.components[0][0:4].upper() == 'BOST':
            units = r'${\Omega}$-m'
        elif self.components[0][0].upper() == 'P' or self.components[0][0:4].upper() in ('BETA', 'AZIM'):
            units = 'Degrees'
        
        return units

    def get_designated_colour(self, component, ii):
        """Colour for a component.

        Parameters
        ----------
        component : str
            Component name; matched case-insensitively against the keys
            of :attr:`designated_colours` (``'xy'``, ``'yx'``, ``'det'``,
            ``'ssq'`` …).
        ii : int
            Index into :attr:`colour`, used as the fallback.

        Returns
        -------
        str
            The designated colour if :attr:`use_designated_colours` is
            set and a key matches, otherwise ``colour[ii]``.
        """
        component = component.lower()
        if self.use_designated_colours:
            for key in self.designated_colours.keys():
                if key in component:
                    return self.designated_colours[key]
        return self.colour[ii]

    def new_figure(self):
        """Replace :attr:`fig` with a fresh empty ``Figure`` and clear
        :attr:`axes`."""
        self.fig = Figure()
        self.axes = []

    def redraw_axes(self):
        """Clear the figure and rebuild the subplot grid for the current
        :attr:`num_sites` (via :meth:`gen_tiling`)."""
        self.fig.clear()
        self.axes = []
        self.tiling = self.gen_tiling()
        for ii in range(self.num_sites):
            self.axes.append(self.fig.add_subplot(self.tiling[0], self.tiling[1], ii + 1))

    def gen_tiling(self):
        """Return ``[nrows, ncols]`` for :attr:`num_sites` subplots -- a
        single column for one or two sites, otherwise a near-square
        grid."""
        if self.num_sites < 3:
            tiling = [self.num_sites, 1]
        else:
            s1 = np.floor(np.sqrt(self.num_sites))
            s2 = np.ceil(self.num_sites / s1)
            tiling = [int(s1), int(s2)]
        return tiling

    def refresh(self):
        """Placeholder -- currently a no-op."""
        pass

    def replace_sites(self, sites_in, sites_out):
        """Swap plotted sites in place and redraw only the affected axes.

        Parameters
        ----------
        sites_in : dict
            ``{data_type: [Site, ...]}`` of replacement sites.
        sites_out : list of str
            Names of the currently-plotted sites to replace.

        Notes
        -----
        If the two collections differ in length the longer is truncated.
        Each changed axis is redrawn via :meth:`redraw_single_axis`
        rather than replotting the whole figure.
        """
        snames = self.site_names
        idx = []
        ns_in = max([len(sites) for sites in sites_in.values()])
        if ns_in > len(sites_out):
            sites_in = sites_in[:len(sites_out)]
        elif len(sites_out) > ns_in:
            sites_out = sites_out[:len(sites_in)]
        for dType, site_list in sites_in.items():
            jj = 0
            if sites_in[dType]:
                for ii, site in enumerate(snames):
                    if site in sites_out:
                        try:
                            self.sites[dType][ii] = sites_in[dType][jj]
                        except IndexError:
                            self.sites[dType].append(sites_in[dType][jj])
                        idx.append(ii)
                        jj += 1
            else:
                self.sites[dType] = []

        idx = list(set(idx))
        for ii in idx:
            self.redraw_single_axis(axnum=ii)

    def draw_all(self, sites=None):
        """Placeholder for a future incremental full redraw -- currently
        only stores ``sites`` if none are set."""
        if self.sites is None:
            self.sites = sites

    def redraw_single_axis(self, site_name='', axnum=None):
        """Clear and replot one subplot.

        Parameters
        ----------
        site_name : str, optional
            Station to redraw; its axis index is looked up in
            :attr:`site_names`.
        axnum : int, optional
            Axis index to redraw.  One of ``site_name`` / ``axnum`` must
            be given; ``axnum`` wins if both are.

        Notes
        -----
        Replots every enabled data type for that station, refreshes the
        legend (axis 0 only), bounds, labels and scales.
        """
        if (site_name == '' and axnum is None):
            print('You have to specify either the axis number or site name...')
            return
        elif axnum is None:
            axnum = self.site_names.index(site_name)
        self.axes[axnum].clear()
        site_name = self.site_names[axnum]
        Max = -99999
        Min = 99999
        for ii, dType in enumerate(['raw_data', 'data', 'response', '1d', 'smoothed_data']):
            if self.sites.get(dType):
                try:
                    site = self.sites[dType][axnum]
                except IndexError:
                    site = None
                if site is not None and self.toggles[dType]:
                    if dType != '1d':
                        self.axes[axnum], ma, mi, artist = self.plot_site(site, dType=dType,
                                                                          ax=self.axes[axnum])
                        Max = max(Max, ma)
                        Min = min(Min, mi)
                        if axnum >= len(self.artist_ref[dType]):
                            self.artist_ref[dType].append(artist)
                        else:
                            self.artist_ref[dType][axnum] = artist
                    else:
                        self.axes[axnum], ma, mi, artist = self.plot_site(self.site1D, dType='1d',
                                                                          ax=self.axes[axnum])
        if axnum == 0:
            self.set_legend()
        self.set_bounds(Max=Max, Min=Min, axnum=axnum)
        self.set_labels(axnum, site_name)
        self.set_axis_scales()

    def set_legend(self):
        """Draw the component legend on axis 0 (imaginary components get a
        triangle marker, everything else a circle)."""
        for ii, comp in enumerate(self.components):
            if 'i' in comp.lower():
                marker = 'v'
            else:
                marker = 'o'
            self.axes[0].plot([], [], color=self.get_designated_colour(comp, ii), label=comp, marker=marker)
        leg = self.axes[0].legend()
        leg.get_frame().set_alpha(0.4)

    def set_labels(self, axnum, site_name):
        """Title one subplot with the station name, and add a y-label on
        the leftmost column and an x-label on the bottom row (period, or
        depth for Bostick components)."""
        self.axes[axnum].set_title(site_name)
        rows = self.tiling[0]
        cols = self.tiling[1]
        if axnum % cols == 0:
            self.axes[axnum].set_ylabel('{} ({})'.format(
                self.y_labels[self.components[0][0:3].lower()], self.units), fontsize=self.label_fontsize)
        if axnum + 1 > cols * (rows - 1):
            if 'bost' in self.components[0].lower():
                self.axes[axnum].set_xlabel(r'$Log_{10}$ Depth (km)', fontsize=self.label_fontsize)
            else:
                self.axes[axnum].set_xlabel(r'$Log_{10}$ Period (s)', fontsize=self.label_fontsize)

    def plot_data(self, sites=None):
        """Render the whole figure -- the main entry point.

        Rebuilds the subplot grid if the site count changed, then for
        each data type (``raw_data`` first, so it plots underneath) calls
        :meth:`plot_site` for every station, and finally sets the bounds,
        legend and axis scales.

        Parameters
        ----------
        sites : dict, optional
            ``{data_type: [Site, ...]}``; only used if :attr:`sites` has
            not already been assigned.
        """
        if self.sites is None:
            self.sites = sites
        tiling = self.gen_tiling()
        self.tiling = tiling
        if self.fig is None:
            print('Opening new figure')
            self.new_figure()
        else:
            self.fig.clear()
        if len(self.fig.get_axes()) != self.num_sites:
            self.redraw_axes()
        Max = np.zeros([tiling[1] * tiling[0]])
        Min = np.zeros([tiling[1] * tiling[0]])
        if not self.components:
            self.components = ['ZXYR']
        last_ii = -1
        for jj, dType in enumerate(['raw_data', 'data', 'response', '1d', 'smoothed_data']):  # plot raw first, if avail
            if dType == '1d':
                pass
            else:
                for ii, site in enumerate(self.sites[dType]):
                    last_ii = max(last_ii, ii)
                    if site is not None and self.toggles[dType]:
                        self.axes[ii], ma, mi, artist = self.plot_site(site,
                                                                       dType=dType, ax=self.axes[ii])
                        Max[ii] = max(Max[ii], ma)
                        Min[ii] = min(Min[ii], mi)
                        if ii >= len(self.artist_ref[dType]):
                            self.artist_ref[dType].append(artist)
                        else:
                            self.artist_ref[dType][ii] = artist
                    if jj == 0:
                        self.set_labels(axnum=ii, site_name=site.name)
        self.set_bounds(Max=Max, Min=Min, axnum=list(range(last_ii + 1)))
        self.set_legend()
        self.set_axis_scales()

    @property
    def log10_yscale(self):
        """bool: Whether the y-axis should be log-scaled -- true for
        log10-scaled impedance and for any resistivity / Bostick
        component."""
        if (self.scale.lower() == 'log10' and self.components[0].lower().startswith('z')) or \
           (any(map(self.components[0].lower().__contains__, ('rho', 'bost')))):
           return True
        else:
           return False

    def set_bounds(self, Max, Min, axnum=None):
        """Set y-limits from the data range (padded by a quarter range,
        or a decade each side when log-scaled), honouring
        :attr:`min_ylim` / :attr:`max_ylim`.

        Parameters
        ----------
        Max, Min : float or sequence of float
            Data extremes -- per-axis sequences when ``axnum`` is a list,
            otherwise scalars.
        axnum : int or list of int, optional
            Axis / axes to bound.
        """
        try:
            for ax in axnum:
                if not self.min_ylim:
                    min_y = Min[ax]
                else:
                    min_y = self.min_ylim
                if not self.max_ylim:
                    max_y = Max[ax]
                else:
                    max_y = self.max_ylim
                axrange = max_y - min_y
                self.axes[ax].set_ylim([min_y - axrange / 4, max_y + axrange / 4])
        except TypeError:
            if not self.min_ylim:
                min_y = Min
            else:
                min_y = self.min_ylim
            if not self.max_ylim:
                max_y = Max
            else:
                max_y = self.max_ylim
            axrange = max_y - min_y
            if self.log10_yscale:
                y1 = 10**(np.log10(min_y) - 1)
                y2 = 10**(np.log10(max_y) + 1)
            else:
                y1 = min_y - axrange / 4
                y2 = max_y + axrange / 4
            try:
                self.axes[axnum].set_ylim([y1, y2])
            except ValueError:
                self.axes[axnum].set_ylim([-1, 1])
            except TypeError:
                self.axes[0].set_ylim([y1, y2])

    @utils.enforce_input(axnum=list, x_bounds=list, y_bounds=list)
    def link_axes(self, axnum, x_bounds=None, y_bounds=None):
        """Give every listed axis the same x / y limits.

        Call with keyword arguments (see the module note on
        ``enforce_input``).

        Parameters
        ----------
        axnum : list of int
            Axes to link; ``None`` means all.
        x_bounds, y_bounds : list of float, optional
            ``[low, high]`` limits; each defaults to the union of the
            current bounds across all axes.  ``y_bounds`` is raised to
            powers of 10 when the y-axis is log-scaled.
        """
        if axnum is None:
            axnum = range(0, len(self.axes))
        if not x_bounds:
            x_bounds = (min([ax.get_xbound()[0] for ax in self.axes]),
                        max([ax.get_xbound()[1] for ax in self.axes]))
        if not y_bounds:
            y_bounds = (min([ax.get_ybound()[0] for ax in self.axes]),
                        max([ax.get_ybound()[1] for ax in self.axes]))
        if self.log10_yscale:
                y_bounds = [10**x for x in y_bounds]
        for ax in axnum:
            self.axes[ax].set_xmargin(self.axis_padding)
            self.axes[ax].set_ymargin(self.axis_padding)
            self.axes[ax].set_ylim([y_bounds[0], y_bounds[1]])
            self.axes[ax].set_xlim([x_bounds[0], x_bounds[1]])

    def plot_site(self, site, dType='Data', ax=None, ):
        """Plot every configured component of one station onto one axis.

        For each entry in :attr:`components` the y-values are derived from
        ``site`` according to the component family -- apparent resistivity
        (:func:`gateway_mt.utils.compute_rho`), phase
        (:func:`gateway_mt.utils.compute_phase`), phase-tensor / invariant
        quantities (from ``site.phase_tensors`` / ``site.CART``), Bostick
        depth (:func:`gateway_mt.utils.compute_bost1D`) or a raw data component
        scaled by :attr:`scale` -- and drawn with ``ax.errorbar``.
        Imaginary components use a triangle marker.  A component with no
        usable data is annotated ``'No Data'``.

        Parameters
        ----------
        site : gateway_mt.data_structures.Site
            Station to plot.
        dType : str, optional
            Data type key (``'data'``, ``'raw_data'``, ``'response'``,
            ``'1d'``, ``'smoothed_data'``) -- selects the marker / line
            style from :attr:`marker` and whether error bars are drawn.
        ax : matplotlib.axes.Axes, optional
            Target axis.

        Returns
        -------
        ax : matplotlib.axes.Axes
            The same axis.
        max, min : float
            Data extremes across the plotted components (outlier-filtered
            when :attr:`show_outliers` is ``False``), used for y-limits.
        artist : matplotlib artist
            The artist from the last component drawn.
        """
        def pop_flagged_data(x, y, toplotErr, site, component):
            if component in site.components:
                e = site.used_error[component]
            else:
                e = site.used_error[site.components[0]]
            idx = []
            for ii, ie in enumerate(e):
                if ie == site.REMOVE_FLAG:
                    idx.append(ii)
            if idx:
                x = np.delete(x, idx)
                y = np.delete(y, idx)
                if toplotErr is not None:
                    toplotErr = np.delete(toplotErr, idx)
            return x, y, toplotErr
        response_types = ('response', '1d')
        ma = []
        mi = []
        linestyle = ''
        marker = self.marker[dType.lower()]
        edgewidth = 0
        if marker == 'oo':
            marker = 'o'
            edgewidth = self.edge_scale * self.markersize
        if marker == '-':
            marker = ''
            linestyle = self.linestyle
        if marker == '--':
            marker = ''
            linestyle = '--'
        if not self.components:
            self.components = [site.components[0]]
        if self.errors.lower() == 'raw':
            Err = site.errors
            errtype = 'errors'
        elif self.errors.lower() == 'mapped':
            Err = site.used_error
            errtype = 'used_error'
        elif self.errors.lower() == 'none':
            Err = None
            toplotErr = None
            errtype = 'none'
        if dType.lower() in response_types:
            Err = None
            toplotErr = None
            errtype = 'none'
        if dType.lower() not in self.which_errors:
            Err = None
            toplotErr = None
            errtype = 'none'
        for ii, comp in enumerate(self.components):
            
            try:
                if 'rho' in comp.lower():
                    toplot, e, log10_e = utils.compute_rho(site, calc_comp=comp, errtype=errtype)
                    # If there is dummy Z data then det, ssq, and gav calculations may give INFs.
                    if np.any(np.isinf(toplot)) or np.any(np.isnan(toplot)):
                        raise KeyError
                    ind = np.where(toplot == 0)[0]  # np.where returns a tuple, take first element
                    if ind.size != 0:
                        print('{} has bad points at periods {}'.format(site.name, site.periods[ind]))
                        print('Adjusting values, ignore them.')
                        toplot[ind] = np.max(toplot)
                    # toplot = np.log10(toplot)
                    if dType.lower() not in response_types and self.errors.lower() != 'none':
                        # toplotErr = log10_e
                        toplotErr = e
                elif 'pha' in comp.lower():
                    toplot, e = utils.compute_phase(site,
                                                    calc_comp=comp,
                                                    errtype=errtype,
                                                    wrap=self.wrap_phase)
                    if dType.lower() not in response_types and self.errors.lower() != 'none':
                        toplotErr = e
                elif 'pt' in comp.lower() or 'phi' in comp.lower() or 'beta' in comp.lower() or 'azimuth' in comp.lower():
                    # If PTs are actually the inverted data, take them directly from the site
                    if comp in site.components:
                        toplot = site.data[comp]
                        e = site.used_error[comp]
                    # Otherwise use the associated PT object
                    # If the site doesn't have all the impedance components, plot no data
                    
                    elif 'pt' in comp.lower():
                        toplot = np.array([getattr(site.phase_tensors[ii],
                                                   comp.upper())
                                           for ii in range(site.NP)])
                        e = np.array([getattr(site.phase_tensors[ii],
                                              comp.upper() + '_error')
                                      for ii in range(site.NP)])
                    else:
                        try:
                            toplot = np.array([getattr(site.phase_tensors[ii],
                                                       comp.lower())
                                               for ii in range(site.NP)])
                            e = np.array([getattr(site.phase_tensors[ii],
                                                  comp.lower() + '_error', 0) # Default error of 0 until I implement errors for invariants
                                          for ii in range(site.NP)])
                        except AttributeError:
                            # If the PT object doesn't exist, raise a KeyError (missing impedance components, so no data)
                            raise KeyError

                    # Convert to degrees
                    if dType.lower() not in response_types and self.errors.lower() != 'none':
                        toplotErr = e
                    else:
                        toplotErr = e * 0
                    if self.pt_units.lower() == 'degrees':
                        toplot = np.rad2deg(np.arctan(toplot))
                        try:
                            toplotErr = np.rad2deg(np.arctan(toplotErr))
                        except AttributeError:  # If plotting PTs, toplotErr will be None
                            toplotErr = None
                elif 'bost' in comp.lower():
                    # If the site doesn't have Z data, det calculations will fail
                    try:
                        toplot, depth = utils.compute_bost1D(site, comp=comp)[:2]
                        # toplot = np.log10(toplot)
                        toplotErr = None
                    except ValueError:
                        raise KeyError
                else:
                    toplot = site.data[comp]
                    if Err is not None:
                        toplotErr = Err[comp]
                    if not comp.lower().startswith('t'):
                        if self.scale == 'sqrt(periods)':
                            toplot = toplot * np.sqrt(site.periods)
                            if Err:
                                toplotErr = Err[comp] * np.sqrt(site.periods)
                        elif self.scale == 'periods':
                            toplot = toplot * site.periods
                            if Err:
                                toplotErr = Err[comp] * site.periods
                        elif self.scale.lower() == 'log10':
                            toplot = np.abs(toplot)
                            if Err:
                                toplotErr = np.abs(Err[comp])
                periods = site.periods
                if not self.plot_flagged_data:
                    if 'bost' in comp.lower():
                        depth, toplot, toplotErr = pop_flagged_data(depth, toplot, toplotErr, site, comp)
                    else:
                        periods, toplot, toplotErr = pop_flagged_data(periods, toplot, toplotErr, site, comp)
                    # depth = periods
                if 'i' in comp.lower():
                    if marker == 'o':
                        use_marker = 'v'
                        markersize = self.markersize * 1.25
                    else:
                        markersize = self.markersize
                        use_marker = marker
                else:
                    use_marker = marker
                    markersize = self.markersize
                    # marker = 'o'
                if 'bost' in comp.lower():
                    artist = ax.errorbar(np.log10(depth), toplot, xerr=None,
                                         yerr=None, marker=use_marker,
                                         linestyle=linestyle, color=self.get_designated_colour(comp, ii),
                                         mec=self.mec, markersize=markersize,
                                         mew=edgewidth, picker=3)
                                         # capsize=5)
                else:
                    artist = ax.errorbar(np.log10(periods), toplot, xerr=None,
                                         yerr=toplotErr, marker=use_marker,
                                         linestyle=linestyle, color=self.get_designated_colour(comp, ii),
                                         mec=self.mec, markersize=markersize,
                                         mew=edgewidth, picker=3)
                                         # capsize=5)
                if self.show_outliers and toplot.size != 0:
                    ma.append(max(toplot))
                    mi.append(min(toplot))
                else:
                    if ((self.toggles['raw_data'] and dType.lower() == 'raw_data') or \
                       (not self.toggles['raw_data'])) and (toplot.size != 0):
                        showdata = utils.remove_outliers(toplot, threshold=self.outlier_thresh)
                        ma.append(max(showdata))
                        mi.append(min(showdata))
                    else:
                        ma.append(0)
                        mi.append(0)
            except KeyError:
                artist = ax.text(0, 0, 'No Data')
                ma.append(0)
                mi.append(0)
        self.set_axis_scales()
        ax.format_coord = format_data_coords
        return ax, max(ma), min(mi), artist

    def set_axis_scales(self):
        """Switch every subplot's y-axis to log scale when
        :attr:`log10_yscale` is true (no-op otherwise)."""
        if self.log10_yscale:
            for axnum in range(len(self.site_names)):
                self.axes[axnum].set_yscale('log')


class MapView(object):
    """Plan-view map / pseudosection / model-slice figure.

    Renders station-located quantities onto a single map axis: site
    markers, induction arrows, phase-tensor ellipses and bars,
    interpolated plan pseudosections, RMS maps, and model plan / cross
    sections.  Used by :mod:`gateway_mt.DataGUI.data_plot`,
    :mod:`gateway_mt.ModelGUI.model_viewer` and the ``phase_tensor_map`` /
    ``plot_induction_map`` scripts.

    Like :class:`DataPlotManager` it is configured by setting attributes
    and then calling a ``plot_*`` method.  The data objects are supplied
    through :attr:`site_data`; several methods also expect the caller to
    have assigned :attr:`dataset` (a
    :class:`gateway_mt.data_structures.Dataset`), which ``__init__`` does *not*
    create.

    Parameters
    ----------
    figure : matplotlib.figure.Figure, optional
        Figure to draw on.  If omitted a new one is created.
    **kwargs
        Accepted for signature compatibility but ignored; use
        :meth:`set_kwargs` to bulk-set attributes.

    Attributes
    ----------
    COORD_SYSTEMS : tuple of str
        The recognised coordinate systems: ``'local'``, ``'utm'``,
        ``'latlong'``, ``'lambert'``.
    window : dict
        ``{'figure', 'axes', 'canvas', 'colorbar', 'colorbar_twin'}`` --
        the drawing surface.  ``window['axes']`` is a *list*; plotting
        methods use ``window['axes'][0]``.
    site_data : dict
        ``{'raw_data', 'data', 'response'}`` -> ``Data`` / ``RawData`` /
        ``Response`` objects.  Also reachable through the :attr:`data`,
        :attr:`raw_data` and :attr:`response` properties.
    dataset : gateway_mt.data_structures.Dataset
        Assigned by the caller; required by :meth:`set_spatial_units` and
        the RMS branch of :meth:`plot_locations`.
    site_names : list of str
        Stations to consider, in plot order.
    active_sites : list of str
        Highlighted subset; the rest are :attr:`generic_sites`.
    site_locations : dict
        ``{'generic', 'active', 'all'}`` -> ``(N, 2)`` arrays of
        ``[northing, easting]`` in the current coordinate system, filled
        by :meth:`set_locations` / :meth:`get_locations`.
    coordinate_system : str
        One of :attr:`COORD_SYSTEMS`; the setter reprojects the site
        locations (and :attr:`model`, if set).
    model : gateway_mt.data_structures.Model
        Model for :meth:`plot_plan_view` / :meth:`plot_x_slice` /
        :meth:`plot_y_slice`.
    interpolant : str
        Grid interpolant for :meth:`interpolate` -- ``'nearest'``,
        ``'linear'``, ``'cubic'``, ``'natural'`` (needs the optional
        ``naturalneighbor`` package; falls back to ``'linear'``), or
        ``'rbf-<kernel>'``.
    colourmap, invert_cmap, lut : str, bool, int
        Colour-map name, inversion flag and number of quantisation
        levels; combined by the :attr:`cmap` property.
    rho_cax, phase_cax, phase_split_cax, delta_cax, diff_cax, tipper_cax,
    skew_cax, model_cax, aniso_cax, depth_cax, rms_cax : list of float
        ``[min, max]`` colour-axis ranges per quantity.
    plot_rms, rms_plot_style : bool, list of str
        Whether :meth:`plot_locations` encodes station RMS, and how
        (``'size'``, ``'colour'``, ``'labels'``).
    min_pt_ratio, pt_ratio_cutoff, pt_scale, arc_scale, ellipse_VE,
    ellipse_linewidth, pt_delta_scale, pt_difference_definition,
    pt_rotation_axis
        Phase-tensor ellipse shape / scaling / differencing controls.
    induction_scale, induction_cutoff, induction_error_tol
        Induction-arrow scaling and rejection thresholds.
    padding_scale : float
        Fraction of the data extent added as a margin by
        :meth:`set_axis_limits` / :meth:`interpolate`.
    spatial_units : {'m', 'km'}
        Units for axis labels / readouts.
    site_marker, site_fill, site_interior, site_exterior, markersize,
    edgewidth, edge_scale, image_opacity, mesh, linewidth : misc
        Marker / mesh styling.
    cb_params : dict
        Colour-bar label / tick / rotation parameters used by
        :meth:`set_colorbar`.
    use_colourbar : bool
        Master switch for drawing colour bars.
    actors : dict
        Lists of artists added by each plot method, for later removal.
    """

    COORD_SYSTEMS = ('local', 'utm', 'latlong', 'lambert')

    def __init__(self, figure=None, **kwargs):
        self.window = {'figure': None, 'axes': None, 'canvas': None, 'colorbar': None, 'colorbar_twin': None}
        if figure:
            self.window = {'figure': figure, 'axes': figure.axes, 'canvas': figure.canvas, 'colorbar': None, 'colorbar_twin': None}
        self.axis_padding = 0.1
        self.fake_im = []
        self.colour = 'k'
        self.site_names = []
        self.site_data = {'raw_data': [], 'data': [], 'response': []}
        self._active_sites = []
        self.site_locations = {'generic': [], 'active': [], 'all': []}
        self.toggles = {'raw_data': False, 'data': False, 'response': False}
        self.actors = {'annotation': [], 'locations': [], 'pseudosection': [],
                       'model': [], 'arrows': [], 'ellipses': []}
        self.site_marker = 'o'
        self.site_fill = True
        self.site_interior = 'k'
        self.arrow_colours = {'raw_data': {'R': 'b', 'I': 'r'},
                              'data': {'R': 'b', 'I': 'r'},
                              'response': {'R': 'g', 'I': 'c'}}
        self.linestyle = '-'
        self.site_exterior = {'generic': 'k', 'active': 'k'}
        self.markersize = 5
        self.edgewidth = 2
        self.edge_scale = 0.255
        self.image_opacity = 0.33
        self._coordinate_system = 'local'
        self.artist_ref = {'raw_data': [], 'data': [], 'response': []}
        self.annotate_sites = 'active'
        self.colourmap = 'turbo'
        self.cb_params = {'rotation': 270,
                          'label_1': None,
                          'label_2': None,
                          'mappable_1': None,
                          'mappable_2': None,
                          'label_fontsize': 16,
                          'tick_fontsize':  14,
                          'labelpad': 20,
                          'ticks': None}
        self.invert_cmap = True
        self.has_nn = has_nn
        if self.has_nn:
            self.interpolant = 'natural'
        else:
            self.interpolant = 'nearest'
        self.rho_cax = [1, 5]
        self.phase_cax = [0, 90]
        self.phase_split_cax = [-40, 40]
        self.delta_cax = [0, 100]
        self.diff_cax = [-15, 15]
        self.tipper_cax = [0, 0.5]
        self.skew_cax = [-10, 10]
        self.model_cax = [1, 5]
        self.aniso_cax = [-2, 2]
        self.depth_cax = [0, 3.5]
        self.rms_cax = [0, 4]
        self.model = []
        self.padding_scale = 5
        self.plot_rms = False
        self.rms_plot_style = ['size']
        self.use_colourbar = True
        self.min_pt_ratio = 1 / 3
        self.pt_ratio_cutoff = 0
        self.pt_scale = 1
        self.arc_scale = 0.8
        self.pt_difference_definition = 'Booker'
        self.ellipse_VE = 1
        self.ellipse_linewidth = 1
        self.pt_delta_scale = 2
        self.pt_rotation_axis = 'x'
        self.induction_scale = 5
        self.induction_cutoff = 2
        self.induction_error_tol = 0.5
        self.rho_error_tol = 1
        self.phase_error_tol = 30
        self.include_outliers = True
        self.allowed_std = 3
        self.spatial_units = 'm'
        self.cb_fontsize = 14
        self.label_fontsize = 16
        self.tick_fontsize = 14
        self.mesh = False
        self.lut = 16
        self.linewidth = 0.005
        if figure is None:
            self.new_figure()
        else:
            self.window['figure'] = figure
        self.create_axes()

    @property
    def facecolour(self):
        """str: Marker fill -- :attr:`site_interior` when
        :attr:`site_fill`, otherwise ``'none'``."""
        if self.site_fill:
            return self.site_interior
        else:
            return 'none'

    @property
    def interpolant(self):
        """str: Grid interpolant name.  Setting ``'natural'`` when the
        optional ``naturalneighbor`` package is missing silently
        downgrades to ``'linear'``."""
        return self._interpolant

    @interpolant.setter
    def interpolant(self, interpolant):
        if interpolant == 'natural':
            if not has_nn:
                interpolant = 'linear'
        self._interpolant = interpolant

    @property
    def cmap(self):
        """matplotlib.colors.Colormap: :attr:`colourmap` resolved through
        :func:`gateway_mt.e_colours.colourmaps.get_cmap` with :attr:`lut` levels
        and :attr:`invert_cmap`."""
        return cm.get_cmap(self.colourmap, N=self.lut, invert=self.invert_cmap)

    @property
    def generic_sites(self):
        """list of str: Sites in :attr:`site_names` but not
        :attr:`active_sites` (unordered -- it is a set difference)."""
        return list(set(self.site_names) - set(self.active_sites))

    @property
    def data(self):
        """gateway_mt.data_structures.Data: shortcut for
        ``site_data['data']``."""
        return self.site_data['data']

    @data.setter
    def data(self, data):
        self.site_data['data'] = data

    @property
    def raw_data(self):
        """gateway_mt.data_structures.RawData: shortcut for
        ``site_data['raw_data']``."""
        return self.site_data['raw_data']

    @raw_data.setter
    def raw_data(self, raw_data):
        self.site_data['raw_data'] = raw_data

    @property
    def response(self):
        """gateway_mt.data_structures.Response: shortcut for
        ``site_data['response']``."""
        return self.site_data['response']

    @response.setter
    def response(self, response):
        self.site_data['response'] = response

    @property
    def active_sites(self):
        """list of str: Highlighted stations.  Assigning also refreshes
        :attr:`site_locations` via :meth:`set_locations`; a bare string
        is wrapped in a list."""
        return self._active_sites

    @active_sites.setter
    def active_sites(self, sites=None):
        if not sites:
            sites = []
        elif isinstance(sites, str):
            sites = [sites]
        self._active_sites = sites
        self.set_locations()

    @property
    def coordinate_system(self):
        """str: Active coordinate system (one of :attr:`COORD_SYSTEMS`).

        Assigning a new value re-derives every :attr:`site_locations`
        entry in that system and, if :attr:`model` is set, reprojects it
        too (``'latlong'`` model reprojection is not implemented and is
        skipped with a message).  A non-local system requires an
        initialised ``raw_data`` (see :meth:`verify_coordinate_system`).
        """
        return self._coordinate_system

    @coordinate_system.setter
    def coordinate_system(self, coordinate_system):
        coordinate_system = coordinate_system.lower()
        if coordinate_system in MapView.COORD_SYSTEMS and coordinate_system != self._coordinate_system:
            if self.verify_coordinate_system(coordinate_system):
                if coordinate_system in ('utm', 'lambert') and self._coordinate_system == 'local':
                    correction = utils.center_locs(np.fliplr(self.site_locations['all']))[1]
                else:
                    correction = np.array([0, 0])
                self._coordinate_system = coordinate_system
                generic_sites = list(set(self.site_names) - set(self.active_sites))
                self.site_locations['generic'] = self.get_locations(sites=generic_sites)
                self.site_locations['active'] = self.get_locations(sites=self.active_sites)
                self.site_locations['all'] = self.get_locations(sites=self.site_names)
                if self.model != []:
                    print('Correction factor is: {}'.format(correction))
                    station_center = np.array(utils.center_locs(np.fliplr(self.site_locations['all']))[1])
                    if coordinate_system == 'latlong':
                        print('Model to latlong not implemented')
                        return
                    self.model.project_model(system=coordinate_system, origin=station_center - correction)

    def verify_coordinate_system(self, coordinate_system):
        """Return ``True`` if ``coordinate_system`` can be used now.

        ``'local'`` is always fine; ``'utm'`` / ``'latlong'`` /
        ``'lambert'`` need an initialised ``raw_data`` (EDI locations).
        """
        if coordinate_system.lower() == 'local':
            return True
        elif coordinate_system.lower() in ('utm', 'latlong', 'lambert'):
            if self.site_data['raw_data'].initialized:
                return True
            else:
                return False

    def set_kwargs(self, kwargs):
        """Set any attribute already defined on the instance from
        ``kwargs``; unknown keys are ignored."""
        for key, value in kwargs.items():
            if key in self.__dict__.keys():
                setattr(self, key, value)

    def new_figure(self):
        """Create a fresh ``Figure`` and axis if ``window['figure']`` is
        empty."""
        if not self.window['figure']:
            self.window['figure'] = Figure()
            self.create_axes()

    def create_axes(self):
        """Populate ``window['axes']`` -- add a single ``111`` subplot,
        or adopt the figure's existing axes."""
        if not self.window['figure'].axes:
            self.window['axes'] = self.window['figure'].add_subplot(111)
        else:
            self.window['axes'] = self.window['figure'].axes

    def set_locations(self):
        """Recompute the ``'generic'`` / ``'active'`` / ``'all'`` entries
        of :attr:`site_locations` in the current coordinate system."""
        self.site_locations['generic'] = self.get_locations(sites=self.generic_sites)
        self.site_locations['active'] = self.get_locations(sites=self.active_sites)
        self.site_locations['all'] = self.get_locations(sites=self.site_names)

    def set_spatial_units(self, units):
        """Propagate ``'m'`` / ``'km'`` to :attr:`spatial_units`,
        :attr:`dataset` and the loaded ``raw_data`` / ``data`` objects
        (other values are ignored)."""
        if units.lower() in ('m', 'km'):
            self.spatial_units = units
            self.dataset.spatial_units = units
            if self.site_data['raw_data']:
                self.site_data['raw_data'].spatial_units = units
            if self.site_data['data']:
                self.site_data['data'].spatial_units = units

    def get_locations(self, sites=None, coordinate_system=None, inquire_only=False):
        """Return station locations as an ``(N, 2)`` ``[northing, easting]``
        array.

        Parameters
        ----------
        sites : list of str, optional
            Stations to return (default: all).
        coordinate_system : str, optional
            System to return them in; defaults to
            :attr:`coordinate_system`.  If given and ``inquire_only`` is
            ``False`` this also *sets* :attr:`coordinate_system` (which
            triggers reprojection).
        inquire_only : bool, optional
            If ``True``, use ``coordinate_system`` for this call only
            without changing the instance state.

        Returns
        -------
        numpy.ndarray
            Locations, rotated by the data azimuth if it differs from the
            stored one.
        """
        if not coordinate_system:
            coordinate_system = self.coordinate_system
        else:
            if not inquire_only:
                self.coordinate_system = coordinate_system
        if coordinate_system == 'local':
            azi = self.site_data['data'].azimuth
            check_azi = self.site_data['data'].check_azi()
            locs = self.site_data['data'].get_locs(site_list=sites, azi=azi)
        elif coordinate_system.lower() == 'utm':
            azi = self.site_data['raw_data'].azimuth
            check_azi = self.site_data['raw_data'].check_azi()
            locs = self.site_data['raw_data'].get_locs(sites=sites, mode='utm')
        elif coordinate_system.lower() == 'lambert':
            azi = self.site_data['raw_data'].azimuth
            check_azi = self.site_data['raw_data'].check_azi()
            locs = self.site_data['raw_data'].get_locs(sites=sites, mode='lambert')
        elif coordinate_system.lower() == 'latlong':
            azi = self.site_data['raw_data'].azimuth
            check_azi = self.site_data['raw_data'].check_azi()
            locs = self.site_data['raw_data'].get_locs(sites=sites, mode='latlong')
        if (azi % 360 != check_azi % 360):
            print('Rotating')
            locs = utils.rotate_locs(locs, azi)
        return locs

    def set_axis_labels(self, ax=None, xlabel=None, ylabel=None):
        """Label the map axes.  Explicit ``xlabel`` / ``ylabel`` win;
        otherwise Easting / Northing (projected systems) or Longitude /
        Latitude (``latlong``) are used, with the current spatial
        units."""
        if ax is None:
            ax = self.window['axes'][0]
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=self.label_fontsize)
        else:
            if self.coordinate_system.lower() in ('utm', 'local', 'lambert'):
                ax.set_xlabel('Easting ({})'.format(self.data.spatial_units), fontsize=self.label_fontsize)
            else:
                ax.set_xlabel(r'Longitude ($^{\circ}$)', fontsize=self.label_fontsize)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=self.label_fontsize)
        else:
            if self.coordinate_system.lower() in ('utm', 'local', 'lambert'):
                ax.set_ylabel('Northing ({})'.format(self.data.spatial_units), fontsize=self.label_fontsize)
            else:
                ax.set_ylabel(r'Latitude ($^{\circ}$)', fontsize=self.label_fontsize)
        ax.tick_params(axis='both', which='major', labelsize=self.tick_fontsize)
        
    def set_axis_limits(self, bounds=None, ax=None):
        """Set the map extent.

        Parameters
        ----------
        bounds : sequence of float, optional
            ``[xmin, xmax, ymin, ymax]``.  If omitted, the extent of
            ``site_locations['all']`` plus a ``1 / padding_scale``
            margin is used.
        ax : matplotlib.axes.Axes, optional
        """
        if ax is None:
            ax = self.window['axes'][0]
        if bounds is None:
            min_x, max_x = (min(self.site_locations['all'][:, 1]),
                            max(self.site_locations['all'][:, 1]))
            min_y, max_y = (min(self.site_locations['all'][:, 0]),
                            max(self.site_locations['all'][:, 0]))
            x_pad = (max_x - min_x) / self.padding_scale
            y_pad = (max_y - min_y) / self.padding_scale
            ax.set_xlim([min_x - x_pad, max_x + x_pad])
            ax.set_ylim([min_y - y_pad, max_y + y_pad])
        else:
            ax.set_xlim(bounds[:2])
            ax.set_ylim(bounds[2:])

    def set_colorbar(self, mappable, label='None', cb_slot=1):
        """Attach (or extend) the figure colour bar.

        Parameters
        ----------
        mappable : matplotlib.cm.ScalarMappable
            Mappable to draw the colour bar from.
        label : str, optional
            Colour-bar label.
        cb_slot : {1, 2}, optional
            ``1`` creates the primary colour bar; ``2`` adds a second
            scale -- a shared label if the two mappables span the same
            range, otherwise a twinned secondary axis.  ``> 2`` is not
            implemented.

        Notes
        -----
        The ``cb_slot == 2`` twin-axis branch has a known off-by-one in
        its tick labels.
        """
        if cb_slot == 1:
            self.cb_params['label_1'] = label
            self.cb_params['mappable_1'] = mappable
            self.window['colorbar'] = self.window['figure'].colorbar(mappable)
            self.window['colorbar'].set_label(label,
                                              rotation=self.cb_params['rotation'],
                                              labelpad=self.cb_params['labelpad'],
                                              fontsize=self.cb_params['label_fontsize'])
            self.window['colorbar'].ax.tick_params(labelsize=self.cb_params['tick_fontsize'])
        elif cb_slot == 2 and label != self.cb_params['label_1']:
            self.cb_params['label_2'] = label
            self.cb_params['mappable_2'] = mappable
            m1_min = np.min(self.cb_params['mappable_1'].get_array())
            m1_max = np.max(self.cb_params['mappable_1'].get_array())
            m2_min = np.min(self.cb_params['mappable_2'].get_array())
            m2_max = np.max(self.cb_params['mappable_2'].get_array())
            if ((m1_min == m2_min)
            and (m1_max == m2_max)):
                    # If the mappables are the same, just set one label
                # print('Mappables are equal')
                self.window['colorbar'].set_label('{}, {}'.format(self.cb_params['label_1'], self.cb_params['label_2']),
                                          rotation=self.cb_params['rotation'],
                                          labelpad=self.cb_params['labelpad'],
                                          fontsize=self.cb_params['label_fontsize'])
            else:
                # Known bug here where for some reason the colorbar adds one extra label at the end that is out of bounds
                # print('Mappables are not equal')
                # Otherwise set up a double colorbar
                # self.window['colorbar'] = self.window['colorbar'].ax
                pos = self.window['colorbar'].ax.get_position()
                map_1 = np.linspace(m1_min, m1_max, self.lut)
                map_2 = np.linspace(m2_min, m2_max, self.lut)
                self.window['colorbar'].ax.set_aspect('auto')
                def fwd(x):
                    # ret = np.interp(x, self.cb_params['mappable_1'].get_array(), self.cb_params['mappable_2'].get_array())
                    ret = np.interp(x, map_1, map_2)
                    ret[ret < m2_min] = m2_min
                    ret[ret > m2_max] = m2_max
                    # print([np.min(ret), np.max(ret)])
                    return ret
                def inv(x):
                    # ret = np.interp(x, self.cb_params['mappable_2'].get_array(), self.cb_params['mappable_1'].get_array())
                    ret = np.interp(x, map_2, map_1)
                    ret[ret < m1_min] = m1_min
                    ret[ret > m1_max] = m1_max
                    # print([np.min(ret), np.max(ret)])
                    return ret

                self.window['colorbar_twin'] = self.window['colorbar'].ax.secondary_yaxis(location='left',
                                                                                          functions=(fwd, inv))
                pos.x0 += 0.05
                pos.x1 += 0.05
                self.window['colorbar'].ax.set_position(pos)
                self.window['colorbar_twin'].set_position(pos)
                self.window['colorbar'].ax.set_ylabel(self.cb_params['label_1'], fontsize=self.cb_params['label_fontsize'])
                self.window['colorbar'].ax.yaxis.set_ticks_position('right')
                self.window['colorbar'].ax.yaxis.set_label_position('right')
                self.window['colorbar_twin'].set_ylabel(self.cb_params['label_2'], fontsize=self.cb_params['label_fontsize'])
                # self.window['colorbar_twin'].yaxis.set_ticks_position('left')
                # self.window['colorbar_twin'].yaxis.set_label_position('left')
                # self.window['colorbar_twin'].set_ylim([m2_min, m2_max])
                self.window['colorbar_twin'].tick_params(labelsize=self.cb_params['tick_fontsize'])
                self.window['colorbar_twin'].set_ylim([m2_min, m2_max])
        elif cb_slot > 2:
            print('Three colorbar slots not implemented yet')
            pass


    def interpolate(self, points, vals, n_interp):
        """Interpolate scattered station values onto a regular grid.

        Parameters
        ----------
        points : ndarray, shape (N, >=2)
            Station coordinates; only the first two columns are used
            (a third is passed to the natural-neighbour backend).
        vals : ndarray, shape (N,)
            Value at each station.
        n_interp : int
            Grid size per axis.

        Returns
        -------
        grid_vals : ndarray, shape (n_interp, n_interp)
        grid_x, grid_y : ndarray
            Meshgrid coordinates.
        x_flat, y_flat : ndarray
            The 1-D grid axes.

        Notes
        -----
        The grid spans the data extent plus a ``1 / padding_scale``
        margin.  The method follows :attr:`interpolant`:
        ``linear`` / ``cubic`` / ``nearest`` via
        :func:`scipy.interpolate.griddata`, ``natural`` via
        ``naturalneighbor`` when available, or ``rbf-<kernel>`` via
        :class:`scipy.interpolate.RBFInterpolator` (duplicate points are
        nudged apart first).
        """
        min_x = np.min(points[:, 0])
        max_x = np.max(points[:, 0])
        min_y = np.min(points[:, 1])
        max_y = np.max(points[:, 1])
        x_pad = (max_x - min_x) / self.padding_scale
        min_x, max_x = (min_x - x_pad, max_x + x_pad)
        y_pad = (max_y - min_y) / self.padding_scale
        min_y, max_y = (min_y - y_pad, max_y + y_pad)
        x_flat = np.linspace(min_x, max_x, n_interp)
        y_flat = np.linspace(min_y, max_y, n_interp)
        grid_x, grid_y = np.meshgrid(x_flat, y_flat)
        if self.interpolant.lower() in ('linear', 'cubic', 'nearest'):
            gx, gy = np.mgrid[min_x:max_x:n_interp*1j, min_y:max_y:n_interp*1j]
            grid_vals = griddata(points[:,:2], vals, (gx, gy), method=self.interpolant)
        elif self.interpolant.lower() == 'natural' and has_nn:
            grid_ranges = [[min_x, max_x, n_interp * 1j],
                           [min_y, max_y, n_interp * 1j],
                           [0, 1, 1]]
            grid_vals = np.squeeze(nn.griddata(points, vals, grid_ranges))                           
        elif 'rbf' in self.interpolant.lower():
            # Fix any duplicate points first
            # Rough hack, just add one centimeter to every point
            if len(np.unique(points[:, :2])) != len(points[:, :2]):
                for ii in range(len(points)):
                    points[ii, :2] += ii / 100
            grid = np.array([grid_x, grid_y])
            query_points = np.rollaxis(grid, 0, 3).reshape((grid.shape[1] * grid.shape[2], 2))
            grid_vals = RBFInterpolator(points[:, :2],
                                        vals,
                                        kernel=self.interpolant.split('-')[1].lower(),
                                        epsilon=1,
                                        degree=0)(query_points).reshape((grid.shape[1], grid.shape[2])).T
        
        return grid_vals, grid_x, grid_y, x_flat, y_flat

    def plot_locations(self, cb_slot=1):
        """Scatter the station markers onto the map.

        Generic and active stations are drawn separately.  When
        :attr:`plot_rms` is set, per-station total RMS (from
        ``dataset.rms``) is encoded according to :attr:`rms_plot_style`:
        ``'size'`` scales the markers, ``'colour'`` colour-maps them (and
        adds a colour bar in ``cb_slot``), ``'labels'`` annotates the
        value.

        Parameters
        ----------
        cb_slot : int, optional
            Colour-bar slot for the RMS colour scale (see
            :meth:`set_colorbar`).
        """
        marker_size = {'generic': [self.markersize ** 2], 'active': [self.markersize ** 2]}
        facecolour = {'generic': 'None', 'active': 'None'}
        edgecolour = deepcopy(self.site_exterior)
        edgewidth = self.markersize * self.edge_scale
        cmin, cmax = 1, 1
        if self.plot_rms:
            marker = self.site_marker
            generic_rms = np.array([self.dataset.rms['Station'][site]['Total']
                                    for site in self.generic_sites])
            active_rms = np.array([self.dataset.rms['Station'][site]['Total']
                                   for site in self.active_sites])
            if 'size' in self.rms_plot_style:
                marker_size['generic'] = generic_rms
                marker_size['active'] = active_rms
                marker_size['generic'] = (utils.normalize(marker_size['generic'],
                                                          lower=1, upper=2, explicit_bounds=True) *
                                          self.markersize) ** 2
                marker_size['active'] = (utils.normalize(marker_size['active'],
                                                         lower=1, upper=2, explicit_bounds=True) *
                                         self.markersize) ** 2
            if 'colour' in self.rms_plot_style:
                # cmin = np.min(np.concatenate((generic_rms, active_rms), axis=0))
                # cmax = np.max(np.concatenate((generic_rms, active_rms), axis=0))
                # facecolour['generic'] = utils.normalize(generic_rms,
                #                                         lower=1, upper=2, explicit_bounds=True)
                # facecolour['active'] = utils.normalize(active_rms,
                #                                         lower=1, upper=2, explicit_bounds=True)
                facecolour['generic'] = generic_rms
                facecolour['active'] = active_rms
                if 'size' in self.rms_plot_style:
                    edgecolour['generic'] = None
                    edgecolour['active'] = None
                # else:

        else:
            marker = self.site_marker
            facecolour['active'] = self.facecolour
            facecolour['generic'] = self.facecolour
            # edgewidth = 1
        if len(self.site_locations['generic']) > 0:
            try:
                self.window['axes'][0].scatter(self.site_locations['generic'][:, 1],
                                               self.site_locations['generic'][:, 0],
                                               marker=marker,
                                               s=marker_size['generic'],
                                               edgecolors=edgecolour['generic'],
                                               linewidths=edgewidth,
                                               # facecolors=facecolour,
                                               c=facecolour['generic'],
                                               zorder=9,
                                               cmap=self.cmap,
                                               vmin=self.rms_cax[0], vmax=self.rms_cax[1])
            except IndexError:
                pass

        if self.active_sites:
            self.window['axes'][0].scatter(self.site_locations['active'][:, 1],
                                           self.site_locations['active'][:, 0],
                                           marker=marker,
                                           s=marker_size['active'],
                                           edgecolors=edgecolour['active'],
                                           linewidths=edgewidth,
                                           # facecolors=facecolour,
                                           c=facecolour['active'],
                                           zorder=9,
                                           cmap=self.cmap,
                                           vmin=self.rms_cax[0], vmax=self.rms_cax[1])
        if self.plot_rms:
            if 'colour' in self.rms_plot_style and self.use_colourbar:
                fake_vals = np.linspace(self.rms_cax[0], self.rms_cax[1], self.site_locations['all'].shape[0])
                self.fake_im = self.window['axes'][0].scatter(self.site_locations['all'][:, 1],
                                                              self.site_locations['all'][:, 0],
                                                              c=fake_vals, cmap=self.cmap)
                self.fake_im.set_visible(False)
                label = 'R.M.S. Misfit'
                self.set_colorbar(mappable=self.fake_im, label=label, cb_slot=cb_slot)
                
                # self.window['colorbar'] = self.window['figure'].colorbar(mappable=self.fake_im)
                
                # self.window['colorbar'].set_label(label,
                                                  # rotation=270,
                                                  # labelpad=20,
                                                  # fontsize=18)
                self.window['colorbar'].set_ticks(np.linspace(self.rms_cax[0], self.rms_cax[1], min(self.lut, 24) + 1))
            if 'labels' in self.rms_plot_style:
                for ii, (xx, yy) in enumerate(self.site_locations['active']):
                    self.actors['locations'].append(self.window['axes'][0].annotate('{:5.2f}'.format(active_rms[ii]),
                                                                                     xy=(yy, xx),
                                                                                     xytext=(10, -4), textcoords='offset points'))
                for ii, (xx, yy) in enumerate(self.site_locations['generic']):
                    self.actors['locations'].append(self.window['axes'][0].annotate('{:5.2f}'.format(generic_rms[ii]),
                                                                                     xy=(yy, xx),
                                                                                     xytext=(10, -4), textcoords='offset points'))
        self.set_axis_labels()

    def plot_annotate(self):
        """Label station markers with their names -- the active subset for
        ``annotate_sites == 'active'``, every station for ``'all'``."""
        if self.annotate_sites == 'active':
            for ii, (xx, yy) in enumerate(self.site_locations['active']):
                self.actors['annotation'].append(self.window['axes'][0].annotate(self.active_sites[ii],
                                                                                 xy=(yy, xx),
                                                                                 xytext=(10, 4), textcoords='offset points'))
        elif self.annotate_sites == 'all':
            for ii, (xx, yy) in enumerate(self.site_locations['active']):
                    self.actors['annotation'].append(self.window['axes'][0].annotate(self.active_sites[ii],
                                                                                     xy=(yy, xx),
                                                                                     xytext=(10, 4), textcoords='offset points'))
            for ii, (xx, yy) in enumerate(self.site_locations['generic']):
                self.actors['annotation'].append(self.window['axes'][0].annotate(self.generic_sites[ii],
                                                                                 xy=(yy, xx),
                                                                                 xytext=(10, 4), textcoords='offset points'))
        self.set_axis_limits()

    @utils.enforce_input(data_type=list, normalize=bool, period_idx=int, arrow_type=list, arrow_convention=str)
    def plot_induction_arrows(self, data_type='data', normalize=True, period_idx=1, arrow_type=['R'], arrow_convention='parkinson'):
        """Draw tipper (induction) arrows at one period.

        Call with keyword arguments (see the module note on
        ``enforce_input``).

        Parameters
        ----------
        data_type : list of str
            Any of ``'data'``, ``'raw_data'``, ``'response'``; each is
            drawn in its :attr:`arrow_colours` colour.
        normalize : bool
            If ``True``, all arrows are drawn the same length.
        period_idx : int
            Period index.
        arrow_type : list of {'R', 'I'}
            Real and/or imaginary parts.
        arrow_convention : {'parkinson', 'weise'}
            Real-arrow sign convention (``'weise'`` reverses it).

        Notes
        -----
        Stations whose tipper is flagged at this period are skipped;
        arrows longer than :attr:`induction_cutoff` are zeroed.  A
        quiver key is added per drawn set.
        """
        x_max, x_min = (np.max(self.site_locations['all'][:, 0]), np.min(self.site_locations['all'][:, 0]))
        y_max, y_min = (np.max(self.site_locations['all'][:, 1]), np.min(self.site_locations['all'][:, 1]))
        max_length = np.max((x_max - x_min, y_max - y_min))
        num_keys = 0
        for R_or_L in arrow_type:
            for dType in data_type:
                X, Y = [], []
                colour = self.arrow_colours[dType][R_or_L]
                # if dType.lower() == 'data':
                #     colour = 'k'
                # elif dType.lower() == 'raw_data':
                #     colour = 'k'
                # elif dType.lower() == 'response':
                #     colour = 'r'
                idx = []
                for ii, site in enumerate(self.site_names):
                    site = self.site_data[dType].sites[site]
                    if np.all(np.array([site.used_error[comp][period_idx] for comp in site.TIPPER_COMPONENTS]) == site.REMOVE_FLAG):
                        continue
                    idx.append(ii)
                    if (R_or_L == 'R' and arrow_convention.lower() == 'weise'):
                        reverse = -1
                    else:
                        reverse = 1
                    if 'TZX' + R_or_L in site.components:
                        X.append(-site.data['TZX' + R_or_L][period_idx] * reverse)
                    else:
                        X.append(0)
                    if 'TZY' + R_or_L in site.components:
                        Y.append(-site.data['TZY' + R_or_L][period_idx] * reverse)
                    else:
                        Y.append(0)
                if idx:
                    arrows = np.transpose(np.array((X, Y)))
                    # arrows = utils.normalize_arrows(arrows)
                    lengths = np.sqrt(arrows[:, 0] ** 2 + arrows[:, 1] ** 2)
                    preserved_lengths = np.sqrt(arrows[:, 0] ** 2 + arrows[:, 1] ** 2)
                    largest_arrow = np.max(lengths)
                    lengths[lengths == 0] = 1
                    arrows[lengths > 1, 0] = 2 * arrows[lengths > 1, 0] / lengths[lengths > 1]
                    arrows[lengths > 1, 1] = 2 * arrows[lengths > 1, 1] / lengths[lengths > 1]
                    lengths = np.sqrt(arrows[:, 0] ** 2 + arrows[:, 1] ** 2)
                    cutoff_idx = lengths > self.induction_cutoff

                    arrows[cutoff_idx, :] = 0
                    lengths = np.sqrt(arrows[:, 0] ** 2 + arrows[:, 1] ** 2)
                    largest_arrow = np.max(lengths)

                    if normalize:
                        arrows = 0.5 * arrows / np.transpose(np.tile(lengths, [2, 1]))
                    arrows = arrows * self.induction_scale / (50) # * np.sqrt(lengths))
                    headwidth = max(3, 1)
                    # scale = self.induction_scale * max_length / 100
                    lengths = np.sqrt(arrows[:, 0] ** 2 + arrows[:, 1] ** 2)
                    largest_arrow = np.max(lengths)
                    width = 0.0025
                    quiv_handle = self.window['axes'][0].quiver(self.site_locations['all'][idx, 1],
                                                                self.site_locations['all'][idx, 0],
                                                                arrows[:, 1],
                                                                arrows[:, 0],
                                                                color=colour,
                                                                headwidth=headwidth,
                                                                width=width,
                                                                zorder=10,
                                                                scale=1,
                                                                scale_units=None)
                    self.set_axis_limits()
                else:
                    preserved_lengths = 0
                # If no arrows are drawn
                if not idx:
                    return
                # key_length_idx = np.argmin(np.abs(preserved_lengths - 0.5))
                key_length = 0.5 * self.induction_scale / 50
                horz_scale = (y_max - y_min) / 10
                vert_scale = (x_max - x_min) / 10
                label_pos_vert = x_min + vert_scale
                if R_or_L == 'R':
                    key_label = 'Real {} ({})'.format(dType, arrow_convention.capitalize())
                else:
                    key_label = 'Imaginary {}'.format(dType)
                x_offset = np.floor(num_keys / 2) * horz_scale * 5
                y_offset = (num_keys % 2) * (vert_scale / 2)
                qk = self.window['axes'][0].quiverkey(quiv_handle,
                                                      self.window['axes'][0].get_xlim()[0] + horz_scale * 4 + x_offset,
                                                      self.window['axes'][0].get_ylim()[0] + vert_scale - y_offset,
                                                      key_length, key_label + ', Length = 0.5',
                                                      labelpos='W',
                                                      coordinates='data')
                num_keys += 1

    def pt_fill_limits(self, fill_param, pt_type):
        """Return the ``(lower, upper)`` colour-axis range for a
        phase-tensor fill parameter.

        Parameters
        ----------
        fill_param : str
            The quantity used to colour the ellipses (``'beta'``,
            ``'alpha'``, ``'azimuth'``, ``'delta'``, ``'phi_max'``,
            ``'phi_split_pt'``, a ``'diff_*'`` residual, ...).
        pt_type : {'phi', 'phi_a', 'ua', 'va'}
            Tensor family; changes the range for the ``phi_*`` params.

        Returns
        -------
        tuple of float
        """
        if fill_param.startswith('diff') and 'delta' not in fill_param.lower():
            lower, upper = self.diff_cax
        elif fill_param in ['Lambda']:
            lower, upper = (0, 0.5)
        elif fill_param == 'beta':
            lower, upper = self.skew_cax
        elif fill_param == 'absbeta':
            lower, upper = [0, self.skew_cax[1]]
        elif fill_param in ['alpha', 'azimuth']:
            lower, upper = (-90, 90)
            # lower, upper = (0, 180)
        elif fill_param in ('delta', 'diff_delta'):
            lower, upper = (self.delta_cax[0], self.delta_cax[1])
        elif fill_param in ('phi_split', 'phi_split_z', 'phi_split_pt'):
            if fill_param == 'phi_split_pt':
                lower, upper = [0, self.phase_split_cax[1]]
            else:
                lower, upper = self.phase_split_cax
        elif fill_param.lower() == 'dimensionality':
            lower, upper = (1, 3)
        else:
            if pt_type.lower() == 'phi':
                if fill_param in ['phi_max', 'phi_min', 'det_phi', 'phi_1', 'phi_2', 'phi_3', 'phi_split', 'phi_split_z', 'phi_split_pt']:
                    lower, upper = self.phase_cax
            elif pt_type.lower() == 'phi_a':
                if fill_param in ['phi_max', 'phi_min', 'det_phi', 'phi_1', 'phi_2', 'phi_3', 'phi_split', 'phi_split_z', 'phi_split_pt']:
                    lower, upper = [-1*self.phase_cax[1], self.phase_cax[1]]
            elif pt_type.lower() == 'ua':
                if fill_param in ['phi_max', 'phi_min', 'det_phi', 'phi_1', 'phi_2', 'phi_3', 'phi_split', 'phi_split_z', 'phi_split_pt']:
                    lower, upper = self.rho_cax
            elif pt_type.lower() == 'va':
                if fill_param in ['phi_max', 'phi_min', 'det_phi', 'phi_1', 'phi_2', 'phi_3', 'phi_split', 'phi_split_z', 'phi_split_pt']:
                    lower, upper = [-1 * self.rho_cax[1], self.rho_cax[1]]
        return lower, upper

    def generate_rectangle(self, width, height, angle):
        """Return ``(x, y)`` vertex arrays for a ``width`` x ``height``
        rectangle (built from a squared sine/cosine parametrisation)
        rotated by ``angle`` radians."""
        t = np.arange(0, 2 * np.pi + np.pi / 30, np.pi / 30)
        x = width * (abs(np.cos(t)) * np.cos(t) + abs(np.sin(t)) * np.sin(t))
        y = height * (abs(np.cos(t)) * np.cos(t) - abs(np.sin(t)) * np.sin(t))
        R = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        x, y = np.matmul(R, np.array((x, y)))
        return x, y

    def generate_ellipse(self, phi):
        """Return ``(phi_x, phi_y)`` vertex arrays for the ellipse traced
        by the 2x2 tensor ``phi`` acting on the unit circle."""
        jx = np.cos(np.arange(0, 2 * np.pi + np.pi / 30, np.pi / 30))
        jy = np.sin(np.arange(0, 2 * np.pi + np.pi / 30, np.pi / 30))
        phi_x = phi[0, 0] * jx + phi[0, 1] * jy
        phi_y = phi[1, 0] * jx + phi[1, 1] * jy
        return phi_x, phi_y

    def generate_arc(self, phi, azi_1, azi_2):
        """Return ``(phi_x, phi_y, ok)`` for the arc of the ``phi``
        ellipse swept between azimuths ``azi_1`` and ``azi_2`` (radians).

        ``ok`` is ``False`` (and the coordinates ``None``) when the two
        azimuths are within 3 degrees of each other.  Used to draw the
        rotation arc for Castro-style phase-tensor differences.
        """
        if np.abs(azi_1 - azi_2) < np.pi/60:
            return None, None, False
        if azi_2 > azi_1:
            flip = -1
            jx = np.cos(np.arange(azi_1, azi_2 + np.pi/180, np.pi/180))
            jy = np.sin(np.arange(azi_1, azi_2 + np.pi/180, np.pi/180))
        else:
            flip = 1
            jx = np.cos(np.arange(azi_2, azi_1 + np.pi/180, np.pi/180))
            jy = np.sin(np.arange(azi_2, azi_1 + np.pi/180, np.pi/180))
        phi_x = phi[0, 0] * jx + phi[0, 1] * jy
        phi_y = phi[1, 0] * jx + phi[1, 1] * jy
        return phi_x[::flip], phi_y[::flip], True

    def get_tensor_params(self, pt_type, fill_param, site=None, tensor=None, period_idx=None):
        """Pull the ellipse geometry and fill value for one tensor.

        Parameters
        ----------
        pt_type : {'phi', 'phi_a', ...}
            ``'phi'`` uses ``site.phase_tensors``; ``'phi_a'`` and the
            anisotropy families (``'ua'``, ``'va'``, ...) use
            ``site.CART``.  For the latter the attribute names are
            prefixed with ``pt_type``.
        fill_param : str
            Attribute to return as the fill value.
        site : gateway_mt.data_structures.Site, optional
            Source of the tensor when ``tensor`` is not given.
        tensor : object, optional
            An explicit tensor (e.g. a residual) to read instead.
        period_idx : int, optional
            Period index into ``site``.

        Returns
        -------
        tuple
            ``(tensor, phi, phi_max, phi_min, azimuth, fill_val)``.
        """
        if pt_type == 'phi':
            if not tensor:
                tensor = site.phase_tensors[period_idx]
            phi = tensor.phi
            phi_max = tensor.phi_max
            phi_min = tensor.phi_min
            azimuth = tensor.azimuth
            fill_val = getattr(tensor, fill_param)
        elif pt_type == 'phi_a':
            if not tensor:
                tensor = site.CART[period_idx]
            phi = tensor.phi
            phi_max = tensor.phi_max
            phi_min = tensor.phi_min
            azimuth = tensor.azimuth
            fill_val = getattr(tensor, fill_param)
        else:
            if not tensor:
                tensor = site.CART[period_idx]
            phi = getattr(tensor, pt_type)
            phi_max = getattr(tensor, '_'.join([pt_type, 'phi_max']))
            phi_min = getattr(tensor, '_'.join([pt_type, 'phi_min']))
            azimuth = getattr(tensor, '_'.join([pt_type, 'azimuth']))
            fill_val = getattr(tensor, '_'.join([pt_type, fill_param]))
            if fill_param in ('phi_1', 'phi_2', 'phi_3', 'phi_min', 'phi_max', 'det_phi'):
                if pt_type.lower() == 'ua':
                    fill_val = np.log10(abs(fill_val))
                elif pt_type.lower() == 'va':
                    fill_val = np.sign(fill_val)*np.log10(abs(fill_val))
        return tensor, phi, phi_max, phi_min, azimuth, fill_val

    def resize_ellipse_aspect(self, azimuth, phi_max):
        """Return ``(phi_x, phi_y)`` for an ellipse of semi-major axis
        ``|phi_max|`` and axis ratio :attr:`min_pt_ratio`, oriented at
        ``azimuth`` -- used to stop very flat ellipses from vanishing."""
        R = np.array(((np.cos(azimuth), -np.sin(azimuth)),
                       (np.sin(azimuth), np.cos(azimuth))))
        max_element = np.abs(phi_max)
        minimal_phi = np.array(((max_element, 0), (0, max_element * self.min_pt_ratio)))
        phi_x, phi_y = self.generate_ellipse(minimal_phi)
        for kk in range(len(phi_x)):
            phi_new = np.matmul(R, np.array((phi_x[kk], phi_y[kk])))
            phi_x[kk], phi_y[kk] = phi_new
        return phi_x, phi_y

    def rescale_ellipse(self, phi_x, phi_y, rescale_factor):
        """Return ``[phi_x, phi_y]`` scaled uniformly by
        ``rescale_factor``."""
        return [rescale_factor * x for x in (phi_x, phi_y)]

    @utils.enforce_input(data_type=list, normalize=bool, fill_param=str, period_idx=int, pt_type=str, bostick_depth=float, cb_slot=int)
    def plot_phase_tensor(self, data_type='data', normalize=True, fill_param='Beta', period_idx=1, pt_type=None, bostick_depth=0, cb_slot=1):
        """Draw a phase-tensor ellipse per station at one period.

        Call with keyword arguments (see the module note on
        ``enforce_input``).

        Parameters
        ----------
        data_type : list of str
            One entry plots that tensor; two entries plot the
            ``data - response`` residual (per
            :attr:`pt_difference_definition`), optionally with a
            rotation arc for the ``'castro'`` definition.
        normalize : bool
            Scale every ellipse to the same size.
        fill_param : str
            Quantity used to colour the ellipses (see
            :meth:`pt_fill_limits`).
        period_idx : int
            Period index (ignored if ``bostick_depth`` is given, which
            picks the nearest Bostick depth instead).
        pt_type : str, optional
            Tensor family (default ``'phi'``).
        bostick_depth : float, optional
            Bostick depth (km) to plot at instead of ``period_idx``.
        cb_slot : int, optional
            Colour-bar slot; ``0`` suppresses it.

        Notes
        -----
        Stations with flagged / missing impedance, or an axis ratio
        below :attr:`pt_ratio_cutoff`, are drawn as plain markers.
        """
        ellipses = []
        fill_vals = []
        arcs = []
        arc_fills = []
        rescale_factor = 0
        if not pt_type:
            pt_type = 'phi'
        if fill_param != 'Lambda':
            fill_param = fill_param.lower()
        if len(data_type) == 2:
            data_type = ['data', 'response']            
        X_all, Y_all = self.site_locations['all'][:, 0], self.site_locations['all'][:, 1]
        scale = np.sqrt((np.max(X_all) - np.min(X_all)) ** 2 +
                        (np.max(Y_all) - np.min(Y_all)) ** 2)
        good_idx = []
        bad_idx = []
        for ii, site_name in enumerate(self.site_names):
            site = self.site_data[data_type[0]].sites[site_name]
            if len(data_type) == 1:
                if bostick_depth:
                    depths = utils.compute_bost1D(site=site, method='phase', comp='aav', filter_width=1)[1]
                    period_idx = np.argmin(abs(bostick_depth - depths))
                # Make sure the site actually has phase tensor data
                if (np.all(np.array([site.used_error[comp][period_idx] for comp in site.components]) == site.REMOVE_FLAG) or
                    (not site.phase_tensors) or site.phase_tensors[period_idx].remove_flag):
                    bad_idx.append(ii)
                    continue
                else:
                    tensor, phi, phi_max, phi_min, azimuth, fill_val = self.get_tensor_params(pt_type,
                                                                                              fill_param,
                                                                                              site=site,
                                                                                              period_idx=period_idx)
                
                # if ((phase_tensor.rhoxy_error / phase_tensor.rhoxy < self.rho_error_tol) and
                #     (phase_tensor.rhoyx_error / phase_tensor.rhoyx < self.rho_error_tol) and
                #     (phase_tensor.phasexy_error < self.phase_error_tol) and
                #     (phase_tensor.phaseyx_error < self.phase_error_tol)):
                if True:
                    cont = 1
                    phi_x, phi_y = self.generate_ellipse(phi)
                    radii = np.sqrt(phi_x ** 2 + phi_y ** 2)
                    if np.min(radii) < np.max(radii) * self.pt_ratio_cutoff:
                        cont = 0

                    if np.all(np.array([site.used_error[comp][period_idx] for comp in site.components]) == site.REMOVE_FLAG):
                        cont = 0
                    if np.min(radii) < np.max(radii) * self.min_pt_ratio:
                        phi_x, phi_y = self.resize_ellipse_aspect(azimuth, phi_max)

                    norm_x, norm_y = (phi_x, phi_y)
                else:
                    cont = 0
                    phi_x, phi_y = self.generate_ellipse(phi)
                    norm_x, norm_y = (phi_x, phi_y)
            else:
                if True:  # Buffer for error tolerances at some point...
                    # If one or both of the tensors is not defined (or the period is flagged), the subtraction should be removed
                    cont = 1
                    tensor1 = self.get_tensor_params(pt_type,
                                                     fill_param,
                                                     site=self.site_data[data_type[0]].sites[site_name],
                                                     period_idx=period_idx)[0]
                    tensor2 = self.get_tensor_params(pt_type,
                                                     fill_param,
                                                     site=self.site_data[data_type[1]].sites[site_name],
                                                     period_idx=period_idx)[0]
                    tensor1.difference_definition = self.pt_difference_definition
                    tensor2.difference_definition = self.pt_difference_definition
                    residual_tensor = tensor1 - tensor2
                    azi_1, azi_2 = tensor1.azimuth, tensor2.azimuth
                    # This makes sure that the ellipses are still generated using the proper methods, but the fill can be the differences
                    if pt_type == 'phi':
                        tensor, phi, phi_max, phi_min, azimuth, fill_val  = self.get_tensor_params(pt_type,
                                                                                                   'diff_'+fill_param, 
                                                                                                   tensor=residual_tensor)
                    else:
                        # If its not a phase tensor, I haven't defined fully what a subtraction is yet, so just do the default
                        tensor, phi, phi_max, phi_min, azimuth, fill_val  = self.get_tensor_params(pt_type,
                                                                                                   fill_param, 
                                                                                                   tensor=residual_tensor)
                    phi_x, phi_y = self.generate_ellipse(phi)
                    
                    # Rescale the ellipses so that those with small differences will vanish and large differences will get a bit bigger
                    # Delta ranges from 0 to 100, but is typically quite small, so maybe multiply it by 5?
                    if self.pt_difference_definition.lower() in ('booker', 'heise', 'hill', 'roots_1', 'roots_2'):
                        # Threw a sqrt on it so that extremely large delta's wouldn't grow too much bigger than everything else
                        rescale_factor = np.sqrt(np.max((np.min((self.pt_delta_scale, residual_tensor.delta * self.pt_delta_scale / 100)), 0.01)))
                        # phi_x = np.array([x*residual_tensor.diff_delta/100*3 for x in phi_x])
                        # phi_y = np.array([y*residual_tensor.diff_delta/100*3 for y in phi_y])
                    # elif self.pt_difference_definition.lower() in ('hill', 'castro', 'roots_1', 'roots_2'):
                    else:
                        # For these ones we don't want collapsing tensors
                        rescale_factor = 0
                    radii = np.sqrt(phi_x ** 2 + phi_y ** 2)
                    # if (np.min(radii) < np.max(radii) * self.pt_ratio_cutoff) or (tensor1.remove_flag or tensor2.remove_flag):
                    #     cont = 0
                    #     bad_idx.append(ii)
                    if np.min(radii) < np.max(radii) * self.min_pt_ratio:
                        phi_x, phi_y = self.resize_ellipse_aspect(azimuth, phi_max)                    
                    # if tensor1.remove_flag or tensor2.remove_flag:
                    #     cont = 0
                    #     bad_idx.append(ii)
                                       
                    # norm_x, norm_y = generate_ellipse(self.site_data[data_type[0]].sites[site_name].phase_tensors[period_idx].phi)
                    
            if cont:
                good_idx.append(ii)
                X, Y = X_all[ii], Y_all[ii]
                phi_x, phi_y = (1000 * phi_x / np.abs(phi_max),
                                1000 * phi_y / np.abs(phi_max))
                radius = np.max(np.sqrt(phi_x ** 2 + phi_y ** 2))
                # if radius > 1000:
                # phi_x, phi_y = [(self.pt_scale * scale / (radius * 100)) * x for x in (phi_x, phi_y)]
                # Fit the ellipse size to the size of the survey
                phi_x, phi_y = self.rescale_ellipse(phi_x=phi_x, phi_y=phi_y, rescale_factor=(self.pt_scale * scale / (radius * 100)))
                if rescale_factor:
                    # Rescale the ellipse again if you want individual ones to shrink or grow    
                    phi_x, phi_y = self.rescale_ellipse(phi_x=phi_x, phi_y=phi_y, rescale_factor=rescale_factor)
                ellipses.append([Y - phi_x, X - phi_y])
                # fill_vals.append(getattr(phase_tensor, fill_param))
                fill_vals.append(fill_val)
                if self.pt_difference_definition.lower() == 'castro' and len(data_type) == 2:
                    arc_x, arc_y, ret = self.generate_arc(phi, azi_1, azi_2)

                    # print([ii, arc_x, arc_y])
                    if ret:
                        arc_x, arc_y = (1000 * arc_x / np.abs(phi_max),
                                        1000 * arc_y / np.abs(phi_max))
                        arc_x, arc_y = [(self.arc_scale * scale / (radius * 100)) * x for x in (arc_x, arc_y)]
                        arcs.append(Path(np.array([Y - arc_x, X - arc_y]).T))
                        arc_fills.append([np.concatenate(([Y], Y - arc_x.T, [Y])), np.concatenate(([X], X - arc_y.T, [X]))])
                    else:
                        arcs.append(None)
                        arc_fills.append(None)
        fill_vals = np.array(fill_vals)
        if len(data_type) > 1:
            lower, upper = self.pt_fill_limits('diff_' + fill_param, pt_type)
        else:
            lower, upper = self.pt_fill_limits(fill_param, pt_type)
        # Alpha, beta, and therefore azimuth are already arctan'ed in data_structures
        if fill_param not in ('delta', 'Lambda', 'alpha', 'azimuth', 'beta', 'phi_split', 'phi_split_z', 'phi_split_pt', 'dimensionality') and pt_type in ('phi', 'phi_a'):
            fill_vals = np.rad2deg(np.arctan(fill_vals))
        if fill_param in ['alpha', 'azimuth', 'beta']:
            fill_vals = np.rad2deg(fill_vals)
                # fill_vals[fill_vals < 0] = 180 + fill_vals[fill_vals < 0]
        fill_vals[fill_vals > upper] = upper
        fill_vals[fill_vals < lower] = lower
        norm_vals = utils.normalize_range(fill_vals,
                                          lower_range=lower,
                                          upper_range=upper,
                                          lower_norm=0,
                                          upper_norm=1)
        for ii, ellipse in enumerate(ellipses):
            self.window['axes'][0].fill(ellipse[0], ellipse[1],
                                        color=self.cmap(norm_vals[ii]),
                                        zorder=3,
                                        edgecolor='k',
                                        linewidth=0)
            if self.ellipse_linewidth:
                self.window['axes'][0].plot(ellipse[0], ellipse[1],
                                            'k-', linewidth=self.ellipse_linewidth, zorder=3)
            if arcs:
                # print([ii, arcs[ii]])
                if arcs[ii]:
                    arc = patches.FancyArrowPatch(path=arcs[ii],
                                                  connectionstyle='arc3,rad=0.5',
                                                  arrowstyle='<-',
                                                  color='k', linewidth=1, #hatch='/',
                                                  mutation_scale=5,
                                                  zorder=5)
                    # arc_fill = patches.Polygon(arc_fills[ii], closed=True, facecolor='grey', zorder=4)
                    self.window['axes'][0].add_patch(arc)
                    self.window['axes'][0].fill(arc_fills[ii][0], arc_fills[ii][1],
                                                color='grey', zorder=4, linewidth=0)
        fake_vals = np.linspace(lower, upper, len(fill_vals))
        self.fake_im = self.window['axes'][0].scatter(self.site_locations['all'][good_idx, 1],
                                                      self.site_locations['all'][good_idx, 0],
                                                      c=fake_vals, cmap=self.cmap)
        # fake_im.colorbar()
        self.fake_im.set_visible(False)
        if self.use_colourbar and cb_slot:
            # self.window['colorbar'] = self.window['figure'].colorbar(mappable=self.fake_im)
            label = self.get_label(fill_param)
            self.set_colorbar(mappable=self.fake_im, label=label, cb_slot=cb_slot)
            # self.window['colorbar'].set_label(label,
            #                                   rotation=270,
            #                                   labelpad=20,
            #                                   fontsize=18)

        if bad_idx:
            self.window['axes'][0].scatter(self.site_locations['all'][bad_idx, 1],
                                           self.site_locations['all'][bad_idx, 0],
                                           marker=self.site_marker,
                                           s=self.markersize ** 2,
                                           edgecolors=self.site_exterior['generic'],
                                           linewidths=self.markersize * self.edge_scale,
                                           # facecolors=facecolour,
                                           c=self.site_interior,
                                           zorder=9)
        self.set_axis_limits()

    @utils.enforce_input(data_type=list, normalize=bool, fill_param=str, period_idx=int, pt_type=str, cb_slot=int)
    def plot_phase_bar(self, data_type='data', normalize=True, fill_param='Beta', period_idx=1, pt_type=None, cb_slot=2):
        """Draw a phase-tensor *bar* (rectangle) per station at one
        period -- an alternative to :meth:`plot_phase_tensor` where the
        rectangle's aspect encodes :math:`\\phi_{min} / \\phi_{max}` and
        its orientation the tensor azimuth.  Same keyword-only calling
        convention and ``fill_param`` / ``cb_slot`` meaning as
        :meth:`plot_phase_tensor`.
        """
        rectangles = []
        fill_vals = []
        if not pt_type:
            pt_type = 'phi'
        if fill_param != 'Lambda':
            fill_param = fill_param.lower()
        if len(data_type) == 2:
            data_type = ['data', 'response']
        X_all, Y_all = self.site_locations['all'][:, 0], self.site_locations['all'][:, 1]
        scale = np.sqrt((np.max(X_all) - np.min(X_all)) ** 2 +
                        (np.max(Y_all) - np.min(Y_all)) ** 2)
        good_idx = []
        bad_idx = []
        for ii, site_name in enumerate(self.site_names):

            site = self.site_data[data_type[0]].sites[site_name]
            if (np.all(np.array([site.used_error[comp][period_idx] for comp in site.components]) == site.REMOVE_FLAG) or
                    (not site.phase_tensors) or site.phase_tensors[period_idx].remove_flag):
                bad_idx.append(ii)
                continue
            else:
                good_idx.append(ii)
            tensor, phi, phi_max, phi_min, azimuth, fill_val = self.get_tensor_params(pt_type,
                                                                                      fill_param,
                                                                                      site=site,
                                                                                      period_idx=period_idx)
            xy = [Y_all[ii], X_all[ii]]
            # width = np.abs(phi_max / 5)
            if np.abs(phi_min / phi_max) < self.pt_ratio_cutoff:
                continue
            if np.abs(phi_min / phi_max) < self.min_pt_ratio:
                height = 2 * self.min_pt_ratio
            else:
                height = np.abs(2 * phi_min / phi_max)
            width = 0.2
            # height = 1
            width, height = [(self.pt_scale * scale) * x / 100 for x in (width, height)]
            xy[0] = xy[0] + height * (np.sin(azimuth)) / 2 - width * (np.cos(azimuth)) / 2
            xy[1] = xy[1] - height * (np.cos(azimuth)) / 2 - width * (np.sin(azimuth)) / 2
            rectangles.append((xy, width, height, np.rad2deg(azimuth)))
            #     ellipses.append([Y - phi_x, X - phi_y])
            fill_vals.append(fill_val)
        fill_vals = np.array(fill_vals)
        lower, upper = self.pt_fill_limits(fill_param, pt_type)
        # Alpha, beta, and therefore azimuth are already arctan'ed in data_structures
        if fill_param not in ('delta', 'Lambda', 'alpha', 'azimuth', 'beta', 'phi_split', 'phi_split_z', 'phi_split_pt', 'dimensionality') and pt_type in ('phi', 'phi_a'):
            fill_vals = np.rad2deg(np.arctan(fill_vals))
        if fill_param in ['alpha', 'azimuth', 'beta']:
            fill_vals = np.rad2deg(fill_vals)
        fill_vals[fill_vals > upper] = upper
        fill_vals[fill_vals < lower] = lower
        norm_vals = utils.normalize_range(fill_vals,
                                          lower_range=lower,
                                          upper_range=upper,
                                          lower_norm=0,
                                          upper_norm=1)
        for ii, rectangle in enumerate(rectangles):
            rect = patches.Rectangle(xy=rectangle[0],
                                     width=rectangle[1],
                                     height=rectangle[2],
                                     angle=rectangle[3],
                                     color=self.cmap(norm_vals[ii]),
                                     zorder=4)
            self.window['axes'][0].add_patch(rect)

        if self.use_colourbar:
            fake_vals = np.linspace(lower, upper, len(fill_vals))
            fake_im = self.window['axes'][0].scatter(self.site_locations['all'][good_idx, 1],
                                                     self.site_locations['all'][good_idx, 0],
                                                     c=fake_vals, cmap=self.cmap)
            # fake_im.colorbar()
            fake_im.set_visible(False)
            label = self.get_label(fill_param)
            self.set_colorbar(mappable=fake_im, label=label, cb_slot=cb_slot)
        self.set_axis_limits()

    @utils.enforce_input(data_type=list, fill_param=str, n_interp=int, period_idx=int, cb_slot=int)
    def plan_pseudosection(self, data_type='data', fill_param='rhoxy', n_interp=200, period_idx=0, cb_slot=1):
        """Interpolate a per-station scalar onto the map as a filled
        ``pcolor`` image at one period.

        Call with keyword arguments (see the module note on
        ``enforce_input``).

        Parameters
        ----------
        data_type : list of str
            One entry plots the field; two entries plot the
            ``response - data`` difference.
        fill_param : str
            Quantity to grid -- ``'rhoxy'`` / ``'rhoyx'`` (or ``'rho..-..'``
            for a component difference), ``'phaxy'``, ``'tip'``,
            ``'beta'`` / ``'absbeta'``, ``'azimuth'``, ``'bostxy'``,
            ``'pt_split'``, ``'phi_min'`` / ``'phi_max'``.
        n_interp : int
            Interpolation grid size per axis.
        period_idx : int
            Period index.
        cb_slot : int, optional
            Colour-bar slot; ``0`` suppresses it.

        Notes
        -----
        Only stations carrying the components the parameter needs are
        used; outliers are trimmed unless :attr:`include_outliers`.
        Interpolation follows :attr:`interpolant` (see
        :meth:`interpolate`).  The axis ``format_coord`` is set to read
        back the gridded value.
        """
        vals = []
        if len(data_type) == 2:
            data_type = ['data', 'response']
        for ii, dt in enumerate(data_type):
            data = self.site_data[dt]
            # temp_vals = []
            # good_idx = []
            if 'rho' in fill_param.lower():
                data_label = 'Resistivity'
                use_log = True
                use_sites = [site for site in data.site_names if set(('ZXYR', 'ZYXR')).issubset(set(data.sites[site].real_components))]
                if '-' in fill_param:
                    val = []
                    for site in use_sites:
                        val.append(utils.compute_rho(data.sites[site], calc_comp='xy', errtype='none')[0][period_idx] -
                                   utils.compute_rho(data.sites[site], calc_comp='yx', errtype='none')[0][period_idx])
                        val[-1] = np.sign(val[-1])*np.log10(np.abs(val[-1]))
                    vals.append(val)
                    # vals.append([np.log10(utils.compute_rho(site=data.sites[site], calc_comp='xy', errtype='none')[0][period_idx] - 
                    #                       utils.compute_rho(site=data.sites[site], calc_comp='yx', errtype='none')[0][period_idx]) for site in data.site_names])
                else:
                    vals.append([np.log10(utils.compute_rho(site=data.sites[site],
                                                            calc_comp=fill_param,
                                                            errtype='none')[0][period_idx]) for site in use_sites])
                # vals.append(temp_vals)
            elif 'pha' in fill_param.lower():
                data_label = 'Phase'
                use_log = False
                use_sites = [site for site in data.site_names if set(('ZXYR', 'ZYXR')).issubset(set(data.sites[site].real_components))]
                if '-' in fill_param:
                    vals.append([utils.compute_phase(data.sites[site], calc_comp='xy', errtype='none', wrap=1)[0][period_idx] -
                                 utils.compute_phase(data.sites[site], calc_comp='yx', errtype='none', wrap=1)[0][period_idx] for site in use_sites])
                else:
                    vals.append([utils.compute_phase(site=data.sites[site],
                                                     calc_comp=fill_param,
                                                     errtype='none',
                                                     wrap=1)[0][period_idx] for site in use_sites])
            elif 'tip' in fill_param.lower():
                data_label = 'Tipper Amplitude'
                use_log = False
                use_sites = [site for site in data.site_names if set(data.TIPPER_COMPONENTS).issubset(set(data.sites[site].real_components))]
                vals.append([np.sqrt((data.sites[site].data['TZXR'][period_idx] ** 2 +
                                     data.sites[site].data['TZYR'][period_idx] ** 2)) for site in use_sites])
            elif 'absbeta' in fill_param.lower():
                data_label = '|Skew|'
                use_log = False
                use_sites = [site for site in data.site_names if set(data.IMPEDANCE_COMPONENTS).issubset(set(data.sites[site].real_components))]
                vals.append([abs(np.rad2deg(data.sites[site].phase_tensors[period_idx].beta)) for site in use_sites])
            elif 'beta' in fill_param.lower():
                data_label = 'Skew'
                use_log = False
                use_sites = [site for site in data.site_names if set(data.IMPEDANCE_COMPONENTS).issubset(set(data.sites[site].real_components))]
                vals.append([np.rad2deg(data.sites[site].phase_tensors[period_idx].beta) for site in use_sites])
            elif 'azimuth' in fill_param.lower():
                data_label = 'Azimuth'
                use_log = False
                use_sites = [site for site in data.site_names if set(data.IMPEDANCE_COMPONENTS).issubset(set(data.sites[site].real_components))]
                vals.append([np.rad2deg(data.sites[site].phase_tensors[period_idx].azimuth) for site in use_sites])
            elif 'bost' in fill_param.lower():
                data_label = 'Depth (km)'
                use_log = True
                use_sites = [site for site in data.site_names if set(data.IMPEDANCE_COMPONENTS).issubset(set(data.sites[site].real_components))]
                # vals.append([np.log10(utils.compute_bost1D(data.sites[site], method='phase', comp=fill_param, filter_width=1)[1][period_idx]) for site in use_sites])
                if '-' in fill_param:
                    val = []
                    for site in use_sites:
                        val.append(utils.compute_bost1D(data.sites[site], method='bostick', comp='xy', filter_width=1)[1][period_idx] -
                                   utils.compute_bost1D(data.sites[site], method='bostick', comp='yx', filter_width=1)[1][period_idx])
                        val[-1] = np.sign(val[-1])*np.log10(np.abs(val[-1]))
                    vals.append(val)

                    # vals.append([(utils.compute_bost1D(data.sites[site], method='bostick', comp='xy', filter_width=1)[1][period_idx] -
                                  # utils.compute_bost1D(data.sites[site], method='bostick', comp='yx', filter_width=1)[1][period_idx]) for site in data.site_names])
                else:
                    vals.append([np.log10(utils.compute_bost1D(data.sites[site], method='bostick', comp=fill_param, filter_width=1)[1][period_idx]) for site in use_sites])
            elif 'pt_split' in fill_param.lower():
                data_label = 'Phase Split'
                use_log = False
                use_sites = [site for site in data.site_names if set(data.IMPEDANCE_COMPONENTS).issubset(set(data.sites[site].real_components))]
                vals.append([data.sites[site].phase_tensors[period_idx].phi_split_pt for site in use_sites])
            elif 'phi_min' in fill_param.lower():
                data_label = 'Phase'
                use_log = False
                use_sites = [site for site in data.site_names if set(data.IMPEDANCE_COMPONENTS).issubset(set(data.sites[site].real_components))]
                vals.append([np.rad2deg(np.arctan(data.sites[site].phase_tensors[period_idx].phi_min)) for site in use_sites])
            elif 'phi_max' in fill_param.lower():
                data_label = 'Phase'
                use_log = False
                use_sites = [site for site in data.site_names if set(data.IMPEDANCE_COMPONENTS).issubset(set(data.sites[site].real_components))]
                vals.append([np.rad2deg(np.arctan(data.sites[site].phase_tensors[period_idx].phi_max)) for site in use_sites])
        if len(vals) == 1:
            vals = np.array(vals[0])
            diff = False
        else:
            diff = True
            #  Should be response - data so a larger response gives a +ve percent difference
            if 'rho' in fill_param.lower() or 'bost' in fill_param.lower():
                # vals = 100 * (np.array(vals[1]) - np.array(vals[0])) / np.array(vals[0])
                vals = (np.array(vals[1])) - (np.array(vals[0]))
                use_log = False
            # elif 'pha' in fill_param.lower() or 'tip' in fill_param.lower():
            else:
                vals = (np.array(vals[1]) - np.array(vals[0]))  # / np.array(vals[0])
        if not self.include_outliers:
            val_mean = np.mean(vals)
            val_std = np.std(vals)
            good_idx = ((vals < (val_mean + self.allowed_std * val_std)) &
                        (vals > (val_mean - self.allowed_std * val_std)))
        elif 'tip' in fill_param.lower():
            good_idx = vals < self.induction_cutoff
        else:
            good_idx = abs(vals) < 1e10
        loc_idx = [ii for ii, site in enumerate(data.site_names) if site in use_sites]
        # Replace this with data.locations in case the 'data', 'response', and 'raw' have different site orders
        # use_locations = self.site_locations['all'][loc_idx, :]
        use_locations = data.locations[loc_idx, :]
        loc_x, loc_y = use_locations[good_idx, 1], use_locations[good_idx, 0]
        loc_z = np.zeros(loc_x.shape)
        points = np.transpose(np.array((loc_x, loc_y, loc_z)))
        # min_x, max_x = (min(loc_x), max(loc_x))
        # x_pad = (max_x - min_x) / self.padding_scale
        # min_x, max_x = (min_x - x_pad, max_x + x_pad)
        # min_y, max_y = (min(loc_y), max(loc_y))
        # y_pad = (max_y - min_y) / self.padding_scale
        # min_y, max_y = (min_y - y_pad, max_y + y_pad)
        # X = np.linspace(min_x, max_x, n_interp)
        # Y = np.linspace(min_y, max_y, n_interp)
        grid_vals, grid_x, grid_y, x_flat, y_flat = self.interpolate(points, vals[good_idx], n_interp)
        if diff and 'rho' in fill_param.lower():
            cax = self.diff_cax
            fill_param = 'Log10 Difference'
        elif (diff and 'pha' in fill_param.lower() or fill_param.lower() == ('phaxy-yx')) or \
             (diff and 'pt_split' in fill_param.lower()) or (diff and 'phi' in fill_param.lower()):
            cax = self.diff_cax
            fill_param = r'Difference ($^{\circ}$)'
        # Any other difference just use a generic label
        elif diff: # and 'tip' in fill_param.lower(): 
            cax = self.diff_cax
            fill_param = r'Difference'
        elif 'rho' in fill_param.lower():
            cax = self.rho_cax
        elif 'pha' in fill_param.lower() or 'phi' in fill_param.lower():
            cax = self.phase_cax
        elif 'tip' in fill_param.lower():
            cax = self.tipper_cax
        elif 'absbeta' in fill_param.lower():
            cax = [0, self.skew_cax[1]]
        elif 'beta' in fill_param.lower():
            cax = self.skew_cax
        elif 'azimuth' in fill_param.lower():
            cax = [-90, 90]
        elif 'bost' in fill_param.lower():
            cax = self.depth_cax
        elif 'pt_split' in fill_param.lower():
            cax = [0, self.phase_split_cax[1]]
        
        im = self.window['axes'][0].pcolor(grid_x, grid_y, grid_vals.T,
                                           cmap=self.cmap,
                                           vmin=cax[0], vmax=cax[1],
                                           zorder=0)
                                           # shading='auto')
        # if self.window['colorbar'] and self.use_colourbar:
        if cb_slot and self.use_colourbar:
            # self.window['colorbar'].remove()
            fake_vals = np.linspace(cax[0], cax[1], self.lut)
            fake_im = self.window['axes'][0].scatter(self.site_locations['all'][0:self.lut, 1],
                                                     self.site_locations['all'][0:self.lut, 0],
                                                     c=fake_vals, cmap=self.cmap)
            # fake_im.colorbar()
            fake_im.set_visible(False)
            self.set_colorbar(mappable=fake_im, label=data_label, cb_slot=cb_slot)
            # self.window['colorbar'] = self.window['figure'].colorbar(mappable=im)
            # self.window['colorbar'].set_label(fill_param,
            #                                   rotation=270,
            #                                   labelpad=20,
            #                                   fontsize=18)
        self.set_axis_limits()
        self.window['axes'][0].format_coord = format_model_coords(im,
                                              X=x_flat, Y=y_flat,
                                              x_label='Easting ({})'.format(self.spatial_units),
                                              y_label='Northing ({})'.format(self.spatial_units),
                                              data_label=data_label,
                                              use_log=use_log)

    def get_label(self, param):
        """Return the LaTeX axis / colour-bar label for a phase-tensor
        parameter name, or ``''`` if it is not recognised."""
        if param.lower() == 'alpha':
            label = r'$\alpha$'
        elif param.lower() == 'beta':
            label = r'$\beta$'
        elif param.lower() == 'phi_1':
            label = r'$\phi_1$'
        elif param.lower() == 'phi_2':
            label = r'$\phi_2$'
        elif param.lower() == 'phi_3':
            label = r'$\phi_3$'
        elif param.lower() == 'det_phi':
            label = r'$det \phi$'
        elif param.lower() == 'azimuth':
            label = r'Azimuth'
        elif param.lower() == 'phi_max':
            label = r'$\phi_{max}$'
        elif param.lower() == 'phi_min':
            label = r'$\phi_{min}$'
        elif param.lower() == 'phi_split_pt':
            label = r'$\phi_{max} - \phi_{min}$'
        else:
            label = ''
        return label

    def plot_plan_view(self, ax=None, z_slice=0, rho_axis='rho_x', cb_slot=0):
        """Plot a depth slice of :attr:`model` as :math:`\\log_{10}`
        resistivity.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
        z_slice : int
            Depth-layer index.
        rho_axis : str
            Model attribute to plot (``'rho_x'`` …); ``'a / b'`` plots
            the ratio of two attributes on the anisotropy colour axis.
        cb_slot : int, optional
            Colour-bar slot; ``0`` suppresses it.
        """
        rho_axis = rho_axis.lower()
        if not ax:
            ax = self.window['axes'][0]
        if self.mesh:
            edgecolor = 'k'
        else:
            edgecolor = None
        if '/' in rho_axis:
            rho_axis = rho_axis.split('/')
            v1 = getattr(self.model, rho_axis[0].strip())
            v2 = getattr(self.model, rho_axis[1].strip())
            vals = v1 / v2
            cax = self.aniso_cax
        else:
            vals = getattr(self.model, rho_axis, 'vals')
            cax = self.model_cax
        # debug_print(vals, 'test.txt')
        X = np.array(self.model.dy)
        Y = np.array(self.model.dx)
        im = ax.pcolormesh(X,
                           Y,
                           np.log10(vals[:, :, z_slice]),
                           cmap=self.cmap,
                           vmin=cax[0], vmax=cax[1],
                           zorder=0,
                           edgecolor=edgecolor,
                           linewidth=self.linewidth)
        self.set_axis_labels(ax=ax)
        # ax.set_ylabel('Northing (km)')
        # ax.set_title('Test')
        # ax.set_xlabel('Easting (km)')
        ax.format_coord = format_model_coords(im,
                                              X=X, Y=Y,
                                              x_label='Easting ({})'.format(self.spatial_units),
                                              y_label='Northing ({})'.format(self.spatial_units),
                                              data_label='Resistivity',
                                              use_log=True)
        if self.use_colourbar and cb_slot:
            self.set_colorbar(mappable=im, label=rho_axis, cb_slot=cb_slot)

    def plot_x_slice(self, ax=None, x_slice=0, rho_axis='rho_x'):
        """Plot a north-south (``x`` index fixed) vertical section of
        :attr:`model` as :math:`\\log_{10}` resistivity.  ``rho_axis``
        has the same meaning as in :meth:`plot_plan_view`."""
        rho_axis = rho_axis.lower()
        if not ax:
            ax = self.window['axes'][0]
        if self.mesh:
            edgecolor = 'k'
        else:
            edgecolor = None
        if '/' in rho_axis:
            rho_axis = rho_axis.split('/')
            v1 = getattr(self.model, rho_axis[0].strip())
            v2 = getattr(self.model, rho_axis[1].strip())
            vals = v1 / v2
            cax = self.aniso_cax
        else:
            vals = getattr(self.model, rho_axis, 'vals')
            cax = self.model_cax
        X = np.array(self.model.dy)
        Y = np.array(self.model.dz)
        im = ax.pcolormesh(X,
                           Y,
                           np.squeeze(np.log10(vals[x_slice, :, :])).T,
                           cmap=self.cmap,
                           vmin=cax[0], vmax=cax[1],
                           zorder=0,
                           edgecolor=edgecolor,
                           linewidth=self.linewidth)
        self.set_axis_labels(ax=ax, xlabel='Easting ({})'.format(self.spatial_units),
                                    ylabel='Depth ({})'.format(self.spatial_units))
        ax.format_coord = format_model_coords(im,
                                              X=X, Y=Y,
                                              x_label='Easting ({})'.format(self.spatial_units),
                                              y_label='Depth ({})'.format(self.spatial_units))

    def plot_y_slice(self, ax=None, y_slice=0, orientation='xz', rho_axis='rho_x'):
        """Plot an east-west (``y`` index fixed) vertical section of
        :attr:`model`.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
        y_slice : int
            Column index.
        orientation : {'xz', 'zx'}
            ``'xz'`` puts northing on the x-axis and depth (downward) on
            the y-axis; ``'zx'`` swaps them.
        rho_axis : str
            As in :meth:`plot_plan_view`.
        """
        rho_axis = rho_axis.lower()
        if not ax:
            ax = self.window['axes'][0]
        if self.mesh:
            edgecolor = 'k'
        else:
            edgecolor = None
        if '/' in rho_axis:
            rho_axis = rho_axis.split('/')
            v1 = getattr(self.model, rho_axis[0].strip())
            v2 = getattr(self.model, rho_axis[1].strip())
            vals = v1 / v2
            cax = self.aniso_cax
        else:
            vals = getattr(self.model, rho_axis, 'vals')
            cax = self.model_cax
        if orientation.lower() == 'xz':
            X = np.array(self.model.dx)
            Y = np.array(self.model.dz)
            to_plot = np.squeeze(np.log10(vals[:, y_slice, :])).T
            x_label = 'Northing ({})'.format(self.model.spatial_units)
            y_label = 'Depth ({})'.format(self.model.spatial_units)
        else:
            Y = np.array(self.model.dx)
            X = np.array(self.model.dz)
            to_plot = np.squeeze(np.log10(vals[:, y_slice, :]))
            y_label = 'Northing ({})'.format(self.model.spatial_units)
            x_label = 'Depth ({})'.format(self.model.spatial_units)
        im = ax.pcolormesh(X, Y, to_plot,
                           cmap=self.cmap,
                           vmin=cax[0], vmax=cax[1],
                           zorder=0,
                           edgecolor=edgecolor,
                           linewidth=self.linewidth)
        self.set_axis_labels(ax=ax, xlabel=x_label, ylabel=y_label)
        ax.format_coord = format_model_coords(im, X=X, Y=Y,
                                              x_label=x_label.split()[0],
                                              y_label=y_label.split()[0])

    def plot_image(self, image, extents, ax=None):
        """Draw a background raster (e.g. a georeferenced basemap) at
        ``extents`` = ``[x0, x1, y0, y1]`` with opacity
        :attr:`image_opacity`."""
        if not ax:
            ax = self.window['axes'][0]
        ax.imshow(image, extent=extents, alpha=self.image_opacity, zorder=2)

    @utils.enforce_input(sites=list, data_type=str, fill_param=str, periods=tuple, pt_type=str, x_axis=str, annotate_sites=bool)
    def tensor_ellipse_pseudosection(self, sites=None, data_type='data', fill_param='beta', periods=(), pt_type=None, x_axis='linear', annotate_sites=False):
        """Phase-tensor ellipse pseudosection: station distance across,
        :math:`\\log_{10}` period down.

        Call with keyword arguments (see the module note on
        ``enforce_input``).

        Parameters
        ----------
        sites : list of str
            Stations, in profile order (used as-is for ``x_axis='linear'``,
            otherwise sorted by the chosen coordinate).
        data_type : str
            ``'data'`` / ``'raw_data'`` / ``'response'``.
        fill_param : str
            Ellipse fill quantity (see :meth:`pt_fill_limits`).
        periods : tuple
            ``(lower, upper, skip)`` -- period window in seconds plus a
            stride (every ``skip + 1``-th period).
        pt_type : str, optional
            Tensor family (default ``'phi'``).
        x_axis : {'linear', 'x', 'y', 'lat', 'lon'}
            Horizontal coordinate: cumulative distance, or a projected /
            geographic axis.
        annotate_sites : bool, optional
            Currently unused.
        """
        ellipses = []
        fill_vals = []
        if not pt_type:
            pt_type = 'phi'
        if fill_param != 'Lambda':
            fill_param = fill_param.lower()

        if x_axis.lower() == 'x':
            X = np.array([self.site_data[data_type].sites[site_name].locations['X'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = 'Northing ({})'.format(self.site_data[data_type].spatial_units)
        elif x_axis.lower() == 'y':
            X = np.array([self.site_data[data_type].sites[site_name].locations['Y'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = 'Easting ({})'.format(self.site_data[data_type].spatial_units)
        elif x_axis.lower() in ('lat', 'latitude'):
            X = np.array([self.site_data[data_type].sites[site_name].locations['Lat'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = r'Latitude ($^{\circ}$)'
        elif x_axis.lower() in ('lon', 'long', 'longitude'):
            X = np.array([self.site_data[data_type].sites[site_name].locations['Long'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = r'Longitude ($^{\circ}$)'
        elif x_axis.lower() == 'linear':
            # If linear distance is used, assume the user has the sites in the desired order
            x = np.array([self.site_data[data_type].sites[site_name].locations['X'] for site_name in sites])
            y = np.array([self.site_data[data_type].sites[site_name].locations['Y'] for site_name in sites])
            X = utils.linear_distance(x, y)
            label = 'Distance ({})'.format(self.site_data[data_type].spatial_units)
        # scale = np.sqrt((np.max(X) - np.min(X)) ** 2 +
                        # (np.max(np.log10(periods[0])) - np.log10(np.min(periods[1]))) ** 2)
        x_scale = (np.max(X) - np.min(X))
        y_scale = (np.log10(periods[1]) - np.log10(periods[0]))
        X_all, Y_all = [], []
        for ii, site_name in enumerate(sites):
            site = self.site_data[data_type].sites[site_name]
            for ip, period in enumerate(site.periods):
                if (period >= periods[0] and period <= periods[1]) and (ip % (periods[2] + 1)) == 0:
                    X_all.append(X[ii])
                    Y_all.append(np.log10(period))
                    tensor, phi, phi_max, phi_min, azimuth, fill_val = self.get_tensor_params(pt_type,
                                                                                              fill_param,
                                                                                              site=site,
                                                                                              period_idx=ip)
                    phi_x, phi_y = self.generate_ellipse(phi)
                    radii = np.sqrt(phi_x ** 2 + phi_y ** 2)
                    if np.min(radii) < np.max(radii) * self.pt_ratio_cutoff:
                        X_all.pop()
                        Y_all.pop()
                        continue
                    if np.min(radii) < np.max(radii) * self.min_pt_ratio:
                        phi_x, phi_y = self.resize_ellipse_aspect(azimuth, phi_max)
                    current_x, current_y = X_all[-1], Y_all[-1]
                    phi_x, phi_y = (1000 * phi_x / np.abs(phi_max),
                                    1000 * phi_y / np.abs(phi_max))
                    radius = np.max(np.sqrt(phi_x ** 2 + phi_y ** 2))
                    # if radius > 1000:
                    phi_x = [self.ellipse_VE * (self.pt_scale * x_scale / (radius * 100)) * x for x in phi_x]
                    phi_y = [(self.pt_scale * y_scale / (radius * 100)) * y for y in phi_y]
                    # Note the +phi_x is to mirror the ellipse so that it matches the plan-view orientation
                    ellipses.append([current_x + phi_x, current_y - phi_y]) 
                    fill_vals.append(fill_val)
        fill_vals = np.array(fill_vals)
        lower, upper = self.pt_fill_limits(fill_param, pt_type)
        # Alpha, beta, and therefore azimuth are already arctan'ed in data_structures
        if fill_param not in ('delta', 'Lambda', 'alpha', 'azimuth', 'absbeta', 'beta', 'phi_split', 'phi_split_z', 'phi_split_pt', 'dimensionality') and pt_type in ('phi', 'phi_a'):
            fill_vals = np.rad2deg(np.arctan(fill_vals))
        if fill_param in ['alpha', 'azimuth', 'absbeta', 'beta']:
            fill_vals = np.rad2deg(fill_vals)
        if 'abs' in fill_param:
            fill_vals = np.abs(fill_vals)
        fill_vals[fill_vals > upper] = upper
        fill_vals[fill_vals < lower] = lower
        norm_vals = utils.normalize_range(fill_vals,
                                          lower_range=lower,
                                          upper_range=upper,
                                          lower_norm=0,
                                          upper_norm=1)
        for ii, ellipse in enumerate(ellipses):
            self.window['axes'][0].fill(ellipse[0], ellipse[1],
                                        color=self.cmap(norm_vals[ii]),
                                        zorder=3,
                                        edgecolor='k',
                                        linewidth=0)
            self.window['axes'][0].plot(ellipse[0], ellipse[1],
                                        'k-', linewidth=self.ellipse_linewidth, zorder=3)
        fake_vals = np.linspace(lower, upper, len(fill_vals))
        self.fake_im = self.window['axes'][0].scatter(X_all,
                                                      Y_all,
                                                      c=fake_vals, cmap=self.cmap)
        if x_axis.lower() == 'linear':
            self.window['axes'][0].set_xlabel('Distance ({})'.format(self.spatial_units), fontsize=self.label_fontsize)
        elif x_axis.lower() == 'x':
            self.window['axes'][0].set_xlabel('Northing ({})'.format(self.spatial_units), fontsize=self.label_fontsize)
        elif x_axis.lower() == 'y':
            self.window['axes'][0].set_xlabel('Easting ({})'.format(self.spatial_units), fontsize=self.label_fontsize)
        self.window['axes'][0].set_ylabel('Period (s)', fontsize=self.label_fontsize)

        self.fake_im.set_visible(False)
        if not self.window['colorbar'] and self.use_colourbar:
            # self.window['colorbar'] = self.window['figure'].colorbar(mappable=self.fake_im)
            label = self.get_label(fill_param)
            self.set_colorbar(mappable=self.fake_im, label=label)
            # self.window['colorbar'].set_label(label,
                                              # rotation=270,
                                              # labelpad=20,
                                              # fontsize=18)
        x_diff = abs(max(X_all) - min(X_all))
        y_diff = abs(max(Y_all) - min(Y_all))
        x_min, x_max = min(X_all), max(X_all)
        y_min, y_max = min(Y_all), max(Y_all)
        xlim = [x_min - x_diff/10, x_max + x_diff/10]
        ylim = [y_min - y_diff/10, y_max + y_diff/10]      
        aspect = self.ellipse_VE * (xlim[1] - xlim[0]) / (ylim[1] - ylim[0])
        self.set_axis_limits(bounds=xlim + ylim)
        self.window['axes'][0].set_aspect(aspect)
        self.window['axes'][0].invert_yaxis()

    @utils.enforce_input(sites=list, data_type=str, fill_param=str, periods=tuple, pt_type=str, x_axis=str)
    def tensor_bar_pseudosection(self, sites=None, data_type='data', fill_param='beta', periods=(), pt_type=None, x_axis='linear'):
        """Phase-tensor *bar* pseudosection -- the rectangle
        (:meth:`plot_phase_bar`-style) equivalent of
        :meth:`tensor_ellipse_pseudosection`, with the same keyword-only
        calling convention and ``sites`` / ``periods`` / ``x_axis``
        meaning.
        """
        rectangles = []
        fill_vals = []
        if not pt_type:
            pt_type = 'phi'
        if fill_param != 'Lambda':
            fill_param = fill_param.lower()

        if x_axis.lower() == 'x':
            X = np.array([self.site_data[data_type].sites[site_name].locations['X'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = 'Northing ({})'.format(self.site_data[data_type].spatial_units)
        elif x_axis.lower() == 'y':
            X = np.array([self.site_data[data_type].sites[site_name].locations['Y'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = 'Easting ({})'.format(self.site_data[data_type].spatial_units)
        elif x_axis.lower() in ('lat', 'latitude'):
            X = np.array([self.site_data[data_type].sites[site_name].locations['Lat'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = r'Latitude ($^{\circ}$)'
        elif x_axis.lower() in ('lon', 'long', 'longitude'):
            X = np.array([self.site_data[data_type].sites[site_name].locations['Long'] for site_name in sites])
            idx = np.argsort(X)
            X = X[idx]
            sites = [sites[c] for c in idx]
            label = r'Longitude ($^{\circ}$)'
        elif x_axis.lower() == 'linear':
            # If linear distance is used, assume the user has the sites in the desired order
            x = np.array([self.site_data[data_type].sites[site_name].locations['X'] for site_name in sites])
            y = np.array([self.site_data[data_type].sites[site_name].locations['Y'] for site_name in sites])
            X = utils.linear_distance(x, y)
            label = 'Distance ({})'.format(self.site_data[data_type].spatial_units)
        scale = np.sqrt((np.max(X) - np.min(X)) ** 2 +
                        (np.max(np.log10(periods[0])) - np.log10(np.min(periods[1]))) ** 2)
        x_scale = np.sqrt((np.max(X) - np.min(X)) ** 2)
        y_scale = np.sqrt((np.log10(periods[1]) - np.log10(periods[0])) ** 2)
        X_all, Y_all = [], []
        for ii, site_name in enumerate(sites):
            site = self.site_data[data_type].sites[site_name]
            for ip, period in enumerate(site.periods):
                if (period >= periods[0] and period <= periods[1]) and (ip % (periods[2] + 1)) == 0:
                    X_all.append(X[ii])
                    Y_all.append(np.log10(period))
                    tensor, phi, phi_max, phi_min, azimuth, fill_val = self.get_tensor_params(pt_type,
                                                                                              fill_param,
                                                                                              site=site,
                                                                                              period_idx=ip)
                    xy = [X_all[-1], Y_all[-1]]
                    # width = np.abs(phi_max / 5)
                    if np.abs(phi_min / phi_max) < self.pt_ratio_cutoff:
                        X_all.pop()
                        Y_all.pop()
                        continue
                    if np.abs(phi_min / phi_max) < self.min_pt_ratio:
                        height = 2 * self.min_pt_ratio
                    else:
                        height = np.abs(2 * phi_min / phi_max)
                    width = 0.2
                    height = height / 2
                    x, y = self.generate_rectangle(width, height, azimuth)
                    # # x *= (self.pt_scale * x_scale / (100))
                    # # y *= (self.pt_scale * y_scale / (100))
                    # phi_x = [self.ellipse_VE * (self.pt_scale * x_scale / (radius * 100)) * x for x in phi_x]
                    x = self.ellipse_VE * (self.pt_scale * x_scale / (100)) * x
                    y = (self.pt_scale * y_scale / (100)) * y
                    # Note the +phi_x is to mirror the ellipse so that it matches the plan-view orientation
                    rectangles.append([xy[0] + x, xy[1] - y])
                    fill_vals.append(fill_val)
        fill_vals = np.array(fill_vals)
        lower, upper = self.pt_fill_limits(fill_param, pt_type)
        # Alpha, beta, and therefore azimuth are already arctan'ed in data_structures
        if fill_param not in ('delta', 'Lambda', 'alpha', 'azimuth', 'beta', 'phi_split', 'phi_split_z', 'phi_split_pt', 'dimensionality') and pt_type in ('phi', 'phi_a'):
            fill_vals = np.rad2deg(np.arctan(fill_vals))
        if fill_param in ['alpha', 'azimuth', 'beta']:
            fill_vals = np.rad2deg(fill_vals)
        fill_vals[fill_vals > upper] = upper
        fill_vals[fill_vals < lower] = lower
        norm_vals = utils.normalize_range(fill_vals,
                                          lower_range=lower,
                                          upper_range=upper,
                                          lower_norm=0,
                                          upper_norm=1)
        for ii, rectangle in enumerate(rectangles):
            self.window['axes'][0].fill(rectangle[0], rectangle[1],
                                        color=self.cmap(norm_vals[ii]),
                                        zorder=4,
                                        edgecolor=None)
        self.window['axes'][0].invert_yaxis()
        x_diff = abs(max(X_all) - min(X_all))
        y_diff = abs(max(Y_all) - min(Y_all))
        x_min, x_max = min(X_all), max(X_all)
        y_min, y_max = min(Y_all), max(Y_all)
        xlim = [x_min - x_diff/10, x_max + x_diff/10]
        ylim = [y_min - y_diff/10, y_max + y_diff/10]
        aspect = self.ellipse_VE * (xlim[1] - xlim[0]) / (ylim[1] - ylim[0])
        self.set_axis_limits(bounds=xlim + ylim)
        
        self.window['axes'][0].set_aspect(aspect)
        self.window['axes'][0].invert_yaxis()
        # self.set_axis_limits(bounds=xlim + ylim)