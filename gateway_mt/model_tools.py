"""Working with resistivity models outside the GUIs.

Converting a model to a PyVista grid, sampling a vertical section along a
path, and exporting a model and its stations for ParaView. Model Viewer
(:mod:`gateway_mt.ModelGUI.model_viewer`) is built on these, but they need
no Qt, so scripts can use them too. To trim padding from a model first,
see :meth:`gateway_mt.data_structures.Model.trim`.

Paths are given as the easting and northing of their vertices, in the
model's units.
"""
from copy import deepcopy

import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI


def rho_values(model, rho_axis='rho_x'):
    """Resistivities of ``model`` from ``rho_axis``.

    ``rho_axis`` is a model attribute (``'rho_x'``, ``'rho_y'``,
    ``'rho_z'`` or ``'vals'``) or the ratio of two, written ``'a / b'``.
    """
    if '/' in rho_axis:
        a, b = (getattr(model, axis.strip()) for axis in rho_axis.split('/'))
        return a / b
    return getattr(model, rho_axis)


def model_to_rectgrid(model, resolution=None, rho_axis='rho_x'):
    """Convert a :class:`gateway_mt.data_structures.Model` to a PyVista grid.

    Parameters
    ----------
    model : gateway_mt.data_structures.Model
        The model to convert.  Its cell-edge coordinates ``dx`` / ``dy``
        / ``dz`` become the grid axes (``dz`` is negated so depth points
        down).
    resolution : gateway_mt.data_structures.Model, optional
        A resolution / sensitivity volume on the same mesh.  When given,
        two extra cell-data arrays are added: ``'RawResolution'``
        (:math:`\\log_{10}` of the raw values) and ``'Resolution'`` (the
        values divided by cell volume, clipped at the median, then
        standardised into roughly ``[0, 1]`` for use as a PyVista opacity
        field).  When omitted, ``'Resolution'`` is all ones.
    rho_axis : str, optional
        Where to take resistivities from (see :func:`rho_values`).

    Returns
    -------
    pyvista.RectilinearGrid
        Cell data ``'Resistivity'`` holds :math:`\\log_{10}\\rho` with the
        model's ``(nx, ny, nz)`` array transposed / flipped into PyVista's
        ``(y, x, -z)`` Fortran order.
    """
    import pyvista as pv
    vals = rho_values(model, rho_axis)
    grid = pv.RectilinearGrid(np.array(model.dy),
                              np.array(model.dx),
                              -np.array(model.dz[::-1]))
    vals = np.log10(np.swapaxes(np.flip(vals, 2), 0, 1)).flatten(order='F')
    grid.cell_data['Resistivity'] = vals
    if resolution:
        X, Y, Z = np.meshgrid(resolution.yCS, resolution.xCS, resolution.zCS)
        volume = np.swapaxes(np.flip(X * Y * Z, 2), 0, 1).flatten(order='F')
        resvals = np.swapaxes(np.flip(resolution.vals, 2), 0, 1).flatten(order='F')
        grid.cell_data['RawResolution'] = np.log10(resvals)
        resvals = resvals / (volume.ravel())# ** (2/3))
        resvals[resvals > np.median(resvals.ravel())] = np.median(resvals.ravel())
        resvals = resvals - np.mean(resvals)
        resvals = resvals / np.std(resvals)
        resvals = 0.5 + resvals * np.sqrt(0.5)
        resvals[resvals > 1] = 1
        grid.cell_data['Resolution'] = resvals
    else:
        grid.cell_data['Resolution'] = np.ones(vals.shape)
        grid.cell_data['Resolution'][0] = .99
    return grid


def sample_path(E, N, n_points):
    """Points along the polyline through ``(E, N)``, about ``n_points`` in all.

    The points are shared between the segments by length, follow each
    segment whatever its direction, and include every vertex once.

    Returns
    -------
    qx, qy : numpy.ndarray
    """
    E, N = np.asarray(E, dtype=float), np.asarray(N, dtype=float)
    lengths = np.hypot(np.diff(E), np.diff(N))
    total = np.sum(lengths)
    if total == 0:
        return E[:1], N[:1]
    counts = np.maximum(np.ceil(n_points * lengths / total).astype(int), 2)
    qx, qy = [E[:1]], [N[:1]]
    for ii, count in enumerate(counts):
        t = np.linspace(0, 1, count)[1:]
        qx.append(E[ii] + t * (E[ii + 1] - E[ii]))
        qy.append(N[ii] + t * (N[ii + 1] - N[ii]))
    return np.concatenate(qx), np.concatenate(qy)


def sample_section(model, E, N, n_points=200, rho_axis='rho_x'):
    """:math:`\\log_{10}\\rho` on a vertical section along a path.

    The path through ``(E, N)`` is sampled with :func:`sample_path`, and
    the model is interpolated linearly between cell centres at each of its
    depth mesh lines.

    Parameters
    ----------
    model : gateway_mt.data_structures.Model
    E, N : array_like
        Easting and northing of the path's vertices.
    n_points : int, optional
        About how many points to sample along the path.
    rho_axis : str, optional
        Where to take resistivities from (see :func:`rho_values`).

    Returns
    -------
    qx, qy : numpy.ndarray
        Easting and northing of the points along the path.
    depths : numpy.ndarray
        The model's depth mesh (``model.dz``).
    values : numpy.ndarray
        ``(len(qx), len(depths))`` :math:`\\log_{10}\\rho`.
    """
    qx, qy = sample_path(E, N, n_points)
    depths = np.array(model.dz)
    cell_N, cell_E, cell_z = model.cell_centers()
    vals = np.transpose(np.log10(rho_values(model, rho_axis)), [1, 0, 2])
    interpolator = RGI((cell_E, cell_N, cell_z), vals, bounds_error=False, fill_value=None)
    points = np.column_stack([np.repeat(qx, len(depths)),
                              np.repeat(qy, len(depths)),
                              np.tile(depths, len(qx))])
    values = interpolator(points).reshape(len(qx), len(depths))
    return qx, qy, depths, values


def nearest_point(x, y, px, py):
    """Index of, and distance to, the point of ``(px, py)`` nearest ``(x, y)``."""
    distances = np.hypot(np.asarray(px, dtype=float) - x, np.asarray(py, dtype=float) - y)
    idx = int(np.argmin(distances))
    return idx, distances[idx]


def project_onto_path(px, py, E, N):
    """Where points sit relative to the polyline through ``(E, N)``.

    Returns
    -------
    along : numpy.ndarray
        Distance along the path to the nearest point on it.
    across : numpy.ndarray
        Distance from the path.
    """
    px, py = np.atleast_1d(np.asarray(px, dtype=float)), np.atleast_1d(np.asarray(py, dtype=float))
    E, N = np.asarray(E, dtype=float), np.asarray(N, dtype=float)
    dE, dN = np.diff(E), np.diff(N)
    lengths = np.hypot(dE, dN)
    start = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    along = np.zeros(len(px))
    across = np.full(len(px), np.inf)
    for ii in range(len(lengths)):
        if lengths[ii] == 0:
            continue
        t = np.clip(((px - E[ii]) * dE[ii] + (py - N[ii]) * dN[ii]) / lengths[ii] ** 2, 0, 1)
        dist = np.hypot(px - (E[ii] + t * dE[ii]), py - (N[ii] + t * dN[ii]))
        closer = dist < across
        across[closer] = dist[closer]
        along[closer] = start[ii] + t[closer] * lengths[ii]
    return along, across


def export_to_paraview(outfile, model, data=None, raw_data=None, resolution=None,
                       utm_zone=None, use_elevation=False, local=True):
    """Write a model (and stations) for ParaView, as the ``to_vtk`` command does.

    Writes ``<outfile>_model.vtk`` and, if ``data`` has stations,
    ``<outfile>_sites.vtk``. The model and data are written in metres
    whatever units they are shown in.

    Parameters
    ----------
    outfile : str
        Base name for the files.
    model : gateway_mt.data_structures.Model
    data : gateway_mt.data_structures.Data, optional
    raw_data : gateway_mt.data_structures.RawData, optional
        If given (and ``local``), its station centre is used as the origin,
        so the files are in UTM coordinates. Otherwise the origin is 0, 0.
    resolution : numpy.ndarray, optional
        Resolution values on the model mesh, written alongside it.
    utm_zone : str, optional
        Written into the file headers ('999' if not known).
    use_elevation : bool, optional
        Place the stations at their elevations.
    local : bool, optional
        The model and data are in local coordinates (False if they are
        already in UTM etc., when no origin is added).

    Returns
    -------
    list of str
        The files written.
    """
    model = deepcopy(model)
    model.spatial_units = 'm'
    origin = (0, 0)
    if local and raw_data is not None and raw_data.initialized:
        raw_data = deepcopy(raw_data)
        raw_data.spatial_units = 'm'
        origin = raw_data.origin
    zone = utm_zone if utm_zone and str(utm_zone).lower() != 'undefined' else '999'
    if resolution is not None:
        model.resolution = resolution
    model.to_vtk(outfile, origin=origin, UTM=zone)
    written = ['{}_model.vtk'.format(outfile)]
    if data is not None and len(data.site_names):
        data = deepcopy(data)
        data.spatial_units = 'm'
        data.to_vtk(outfile, UTM=zone, origin=origin, use_elevation=use_elevation)
        written.append('{}_sites.vtk'.format(outfile))
    return written


def write_section_csv(path, easting, northing, depths, values, units):
    """Write a transect section to CSV: easting, northing, distance along the
    transect, depth and log10 resistivity, one line per point.

    ``values`` is ``(len(easting), len(depths))``, as from :func:`sample_section`.
    """
    distance = np.concatenate([[0], np.cumsum(np.hypot(np.diff(easting), np.diff(northing)))])
    with open(path, 'w') as f:
        f.write('# Easting ({0}), Northing ({0}), Distance ({0}), Depth ({0}), log10 Resistivity (ohm-m)\n'.format(units))
        for ix in range(len(easting)):
            for iz, depth in enumerate(depths):
                f.write('{:.4f},{:.4f},{:.4f},{:.4f},{:.5f}\n'.format(easting[ix], northing[ix], distance[ix],
                                                                      depth, values[ix, iz]))
