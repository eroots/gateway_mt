"""Stateless helper functions shared across :mod:`pyMT`.

This module is the package's utility layer.  It has no classes and does not
import any other :mod:`pyMT` module, so everything else -- :mod:`pyMT.data_structures`,
:mod:`pyMT.IO`, :mod:`pyMT.gplot`, the GUIs and the scripts -- is free to import
it.  Roughly, the functions fall into these groups:

``Container / type coercion``
    :func:`to_list`, :func:`list2np`, :func:`np2list`, :func:`flatten_list`,
    :func:`list_or_numpy`, :func:`enforce_input`, :func:`validate_input`.

``Numeric helpers``
    :func:`edge2center`, :func:`truncate`, :func:`normalize`,
    :func:`normalize_range`, :func:`percdiff`, :func:`closest_periods`,
    :func:`skin_depth`, :func:`is_all_empties`.

``Filtering``
    :func:`remove_outliers`, :func:`geotools_filter`.

``Site-location geometry``
    :func:`center_locs`, :func:`rotate_locs`, :func:`project_to_line`,
    :func:`linear_distance`.

``Coordinate conversion / map projection``
    :func:`parse_dms`, :func:`dms2dd`, :func:`dd2dms`, :func:`compute_longorigin`,
    :func:`geo2utm`, :func:`zone_convert`, :func:`project`, :func:`unproject`,
    :func:`generic_projection`, :func:`to_lambert`.

``MT transfer-function computations`` (operate on a
:class:`pyMT.data_structures.Site`)
    :func:`compute_MT_determinant`, :func:`compute_gav`, :func:`compute_aav`,
    :func:`compute_ssq`, :func:`compute_rho`, :func:`compute_phase`,
    :func:`compute_bost1D`, :func:`convert2impedance`, :func:`calculate_misfit`.

``Mesh / model``
    :func:`generate_zmesh`, :func:`generate_lateral_mesh`, :func:`regrid_model`.

``Filesystem``
    :func:`check_file`, :func:`fileparts`, :func:`check_extension`.
"""
import os
from copy import deepcopy
from math import log10, floor

import numpy as np
import pyproj
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.ndimage import median_filter


#: Magnetic permeability of free space (H/m).
MU = 4 * np.pi * 1e-7


# --------------------------------------------------------------------------- #
# Container / type coercion
# --------------------------------------------------------------------------- #
def to_list(obj):
    """Return ``obj`` as a list.

    Parameters
    ----------
    obj : object
        Any value.  Strings are treated as scalars (wrapped, not split).

    Returns
    -------
    list
        ``[obj]`` if ``obj`` is a scalar or string, otherwise ``list(obj)``.
    """
    if not hasattr(obj, '__iter__') or isinstance(obj, str):
        obj = [obj]
    else:
        obj = list(obj)
    return obj


def list2np(data):
    """Convert ``data`` to a :class:`numpy.ndarray`.

    Parameters
    ----------
    data : array_like

    Returns
    -------
    numpy.ndarray

    Raises
    ------
    ValueError
        Re-raised (after printing ``data``) if ``data`` cannot be turned into
        an array, e.g. a ragged nested sequence.
    """
    try:
        return np.array(data)
    except ValueError as e:
        print('Error here')
        print(data)
        print('Error over')
        raise(e)


def np2list(data):
    """Return the elements of ``data`` as a plain list.

    Parameters
    ----------
    data : iterable

    Returns
    -------
    list
    """
    return [x for x in data]


def flatten_list(List):
    """Flatten one level of nesting.

    Parameters
    ----------
    List : sequence of sequence
        A list of lists (or other iterables).

    Returns
    -------
    list
        The concatenated elements.  Elements are cast to ``float`` when every
        element allows it, otherwise they are returned unchanged.
    """
    try:
        ret = [float(point) for sublist in List for point in sublist]
    except ValueError:
        ret = [point for sublist in List for point in sublist]
    return ret


def list_or_numpy(func):
    """Decorate a 1-argument function so it accepts either a list or an array.

    If the wrapped function is called with a :class:`numpy.ndarray`, the array
    is converted to a list before the call and the result is converted back to
    an array afterwards.  Any other input type is passed through untouched.

    Parameters
    ----------
    func : callable
        A function taking a single sequence argument.

    Returns
    -------
    callable
    """
    def inner(obj):
        made_change = False
        if isinstance(obj, np.ndarray):
            obj = [x for x in obj]
            made_change = True
        retval = func(obj)
        if made_change:
            retval = np.array(retval)
        return retval
    return inner


def enforce_input(**types):
    """Decorator factory that coerces named keyword arguments to given types.

    Each keyword passed to :func:`enforce_input` names an argument of the
    decorated function and the type it should have.  When the function is
    called, any of those keyword arguments that is not already an instance of
    the requested type (and is not ``None``) is coerced: ``list`` via
    :func:`to_list`, :class:`numpy.ndarray` via :func:`numpy.array`, anything
    else by calling the type.

    Parameters
    ----------
    **types : type
        ``argument_name=required_type`` pairs.

    Returns
    -------
    callable
        A decorator.

    Examples
    --------
    >>> @enforce_input(periods=np.ndarray)
    ... def f(periods=None):
    ...     return type(periods).__name__
    >>> f(periods=[1, 2, 3])
    'ndarray'
    """
    def decorator(func):
        def new_func(*args, **kwargs):
            newargs = {}
            for k, t in types.items():
                arg = kwargs.get(k, None)
                if isinstance(arg, t) or arg is None:
                    newargs.update({k: arg})
                else:
                    if t == list:
                        newargs.update({k: to_list(kwargs[k])})
                    elif t is np.ndarray:
                        newargs.update({k: np.array(kwargs[k])})
                    else:
                        newargs.update({k: t(kwargs[k])})
            return func(*args, **newargs)
        return new_func
    return decorator


def validate_input(inval, expected_type):
    """Try to cast ``inval`` to ``expected_type``.

    Parameters
    ----------
    inval : object
        Value to convert (typically a string from a GUI field).
    expected_type : type
        Target type.

    Returns
    -------
    object or bool
        The converted value on success; ``False`` (after printing a message)
        if the conversion raises :class:`ValueError`.
    """
    try:
        retval = expected_type(inval)
    except ValueError:
        print('Conversion not possible. Please enter a valid input (of type {})'.format(expected_type))
        return False
    else:
        return retval


# --------------------------------------------------------------------------- #
# Numeric helpers
# --------------------------------------------------------------------------- #
def edge2center(x):
    """Return the midpoints of consecutive elements of ``x``.

    Parameters
    ----------
    x : array_like
        1-D array of ``N`` cell edges (e.g. mesh node positions).

    Returns
    -------
    numpy.ndarray
        Array of ``N - 1`` cell centres.
    """
    y = np.zeros(len(x) - 1)
    for ii in range(len(x) - 1):
        y[ii] = (x[ii] + x[ii + 1]) / 2
    return y


@list_or_numpy
def truncate(nums):
    """Round to five significant figures, element-wise.

    Parameters
    ----------
    nums : float or sequence of float or numpy.ndarray or None
        Value(s) to round.  Zeros are left as zero.  ``None`` returns ``None``.

    Returns
    -------
    same type as ``nums``
        The rounded value(s).  A list in, a list out; an array in, an array out
        (handled by the :func:`list_or_numpy` decorator).

    Notes
    -----
    The number of significant figures is fixed at 5.
    """
    sig = 5
    if nums is not None:
        if hasattr(nums, '__iter__'):
            retval = []
            for x in nums:
                if x == 0:
                    retval.append(0)
                else:
                    retval.append(round(x, -int(floor(log10(abs(x)))) + (sig - 1)))
        else:
            if nums == 0:
                retval = nums
            else:
                retval = round(nums, -int(floor(log10(abs(nums)))) + (sig - 1))
        return retval


def normalize(vals, lower=0, upper=1, explicit_bounds=False):
    """Min-max scale ``vals`` into ``[lower, upper]``.

    Parameters
    ----------
    vals : numpy.ndarray
        Values to scale.
    lower, upper : float, optional
        Target range (``explicit_bounds=False``) or the assumed input range
        (``explicit_bounds=True``).  Defaults ``0`` and ``1``.
    explicit_bounds : bool, optional
        If ``True``, treat ``lower``/``upper`` as the data range and simply map
        ``(vals - lower) / (upper - lower)``.  If ``False`` (default), use
        ``vals.min()``/``vals.max()`` as the data range and map onto
        ``[lower, upper]``.

    Returns
    -------
    numpy.ndarray
    """
    if explicit_bounds:
        norm_vals = (vals - lower) / (upper - lower)
    else:
        norm_vals = (vals - np.min(vals)) / \
                    (np.max(vals) - np.min(vals))
        norm_vals = norm_vals * (upper - lower) + lower
    return norm_vals


def normalize_range(vals, lower_range=0, upper_range=1, lower_norm=0, upper_norm=1):
    """Normalise ``vals`` against an explicitly given input range.

    ``lower_range`` and ``upper_range`` are appended to ``vals``, the whole lot
    is passed through :func:`normalize` onto ``[lower_norm, upper_norm]``, and
    the two appended points are dropped again.  This pins the normalisation to
    a known range rather than to the extremes of ``vals``.

    Parameters
    ----------
    vals : array_like
    lower_range, upper_range : float, optional
        The input range to normalise against.
    lower_norm, upper_norm : float, optional
        The output range.

    Returns
    -------
    numpy.ndarray
        Normalised values, same length as ``vals``.
    """
    vals = list(vals)
    vals.append(lower_range)
    vals.append(upper_range)
    norm_vals = normalize(np.array(vals), lower=lower_norm, upper=upper_norm)
    norm_vals = norm_vals[:-2]
    return norm_vals


def percdiff(val1, val2):
    """Symmetric percentage difference between ``val1`` and ``val2``.

    Parameters
    ----------
    val1, val2 : float or numpy.ndarray

    Returns
    -------
    float or numpy.ndarray
        ``|val1 - val2| / (0.5 * (val1 + val2))``.
    """
    return abs(val1 - val2) / (0.5 * (val1 + val2))


def closest_periods(available_p, wanted_p):
    """Snap each wanted period to the nearest available period.

    Parameters
    ----------
    available_p : sequence of float or numpy.ndarray
        Periods that actually exist in the data.
    wanted_p : sequence of float
        Requested periods.

    Returns
    -------
    list or numpy.ndarray
        Sorted, de-duplicated nearest matches.  Returned as an array when
        ``available_p`` is an array, otherwise a list.
    """
    retval = []
    for p in wanted_p:
        ind = np.argmin([abs(ap - p) for ap in available_p])
        retval.append(available_p[ind])
    retval = sorted(set(retval))
    if isinstance(available_p, np.ndarray):
        retval = np.array(retval)
    return retval


def skin_depth(resistivity, period):
    """Electromagnetic skin depth (metres).

    Parameters
    ----------
    resistivity : float or numpy.ndarray
        Resistivity in ohm-metres.
    period : float or numpy.ndarray
        Period in seconds.

    Returns
    -------
    float or numpy.ndarray
        ``500 * sqrt(resistivity * period)``.
    """
    return 500 * np.sqrt(resistivity * period)


def is_all_empties(dictionary):
    """Return ``True`` if every value in ``dictionary`` is empty/falsey.

    Parameters
    ----------
    dictionary : dict
        Values are expected to be iterables.

    Returns
    -------
    bool
    """
    empties = []
    for key, val in dictionary.items():
        if not any(True for ii in val):
            empties.append(key)
    return set(empties) == set(dictionary.keys())


# --------------------------------------------------------------------------- #
# Filtering
# --------------------------------------------------------------------------- #
def remove_outliers(data, size=7, threshold=1):
    """Replace outliers with a local median.

    Points whose fractional difference (:func:`percdiff`) from a running
    median filter exceeds ``threshold`` are replaced by the median value.

    Parameters
    ----------
    data : numpy.ndarray
        1-D data series.
    size : int, optional
        Width of the median filter window.  Default 7.
    threshold : float, optional
        Fractional-difference threshold above which a point is an outlier.
        Default 1.

    Returns
    -------
    numpy.ndarray
        Copy of ``data`` with outliers replaced.
    """
    med = median_filter(data, size=size)
    diff = percdiff(data, med)
    idx = abs(diff) > threshold
    output = np.array(data)
    output[idx] = med[idx]
    return output


def geotools_filter(x, y, fwidth=1, use_log=True):
    """Adaptive smoothing filter used for Niblett-Bostick transforms.

    Port of the GeoTools smoothing routine.  Each point is scored by how much
    it disagrees with its neighbours (step size and slope change); the output
    is a score-weighted triangular convolution of ``y``.  This suppresses
    isolated spikes while leaving smooth trends essentially untouched.

    Parameters
    ----------
    x : numpy.ndarray
        Abscissa (usually period).
    y : numpy.ndarray
        Ordinate to smooth (usually apparent resistivity or phase).
    fwidth : float, optional
        Filter width in normalised ``x`` units.  Default 1.
    use_log : bool, optional
        Work in log-log space (default).  Automatically disabled if ``x`` or
        ``y`` contain non-positive values that would make ``log10`` invalid.

    Returns
    -------
    numpy.ndarray
        Smoothed ``y``, same shape as the input.
    """
    RTD = 180 / np.pi
    DIFLIMIT = 0.1
    DSLLIMIT = 90.0
    if any(x < 0):
        logx = x
    else:
        logx = np.log10(x)
    if use_log:
        logy = np.log10(y)
        if any(np.isnan(logx)) or any(np.isnan(logy)):
            logx = x
            logy = y
            use_log = False
    else:
        logy = y
    hfwidth = fwidth / 2
    min_x = min(logx)
    min_y = min(logy)
    max_x = max(logx)
    max_y = max(logy)
    X = np.maximum(np.minimum((logx - min_x) / (max_x - min_x), 1), 0)
    Y = np.maximum(np.minimum((logy - min_y) / (max_y - min_y), 1), 0)

    ldif = np.zeros(X.shape)
    lslope = np.zeros(X.shape)
    rslope = np.zeros(X.shape)
    rdif = np.zeros(X.shape)
    ret_y = np.zeros(X.shape)
    score = np.zeros(X.shape)
    # Left and right differences
    for ii in range(len(X)):
        if ii > 0:
            ldif[ii] = Y[ii] - Y[ii - 1]
            if ldif[ii] == 0 and X[ii] == X[ii - 1]:
                lslope[ii] = 0
            else:
                lslope[ii] = RTD * np.arctan2(ldif[ii], X[ii] - X[ii - 1])

        if ii < len(X) - 1:
            rdif[ii] = Y[ii + 1] - Y[ii]
            if rdif[ii] == 0 and X[ii + 1] == X[ii]:
                rslope[ii] = 0
            else:
                rslope[ii] = RTD * np.arctan2(rdif[ii], X[ii + 1] - X[ii])
    ldif[0] = ldif[1]
    rdif[-1] = rdif[-2]
    lslope[0] = lslope[1]
    rslope[-1] = rslope[-2]

    # Generate score for each point
    for ii in range(len(X)):
        pen1 = 0
        pen2 = 0
        totaldif = abs(ldif[ii]) + abs(rdif[ii])
        if (((abs(ldif[ii]) > DIFLIMIT) and (abs(rdif[ii]) > DIFLIMIT) and
                (ldif[ii] * rdif[ii] < 0)) or (totaldif > 4 * DIFLIMIT)):
            pen1 = abs(ldif[ii]) - DIFLIMIT + abs(rdif[ii]) - DIFLIMIT
        dslope = rslope[ii] - lslope[ii]
        if (abs(dslope > DSLLIMIT)):
            pen2 = 0.01 * (abs(dslope) - DSLLIMIT)
        score[ii] = np.maximum(np.minimum(1 - pen1 - pen2, 1), 0.1)

    # Convolve filter with weighted, normalized Y data
    for ii in range(len(X)):
        point = 0
        WTSUM = 0
        for jj in range(len(X)):
            if (abs(X[ii] - X[jj]) <= fwidth):
                FILTER = max(min((hfwidth - abs(X[ii] - X[jj]) / hfwidth), 1), 0)
                WTSUM += FILTER * score[jj]
                point += FILTER * score[jj] * Y[jj]
        ret_y[ii] = (point / (max(WTSUM, 0.0001))) * (max_y - min_y) + min_y
    if use_log:
        ret_y = 10 ** ret_y
    return ret_y


# --------------------------------------------------------------------------- #
# Site-location geometry
# --------------------------------------------------------------------------- #
def center_locs(locations):
    """Recentre an ``(N, 2)`` array of locations on its bounding-box centre.

    Parameters
    ----------
    locations : numpy.ndarray
        Shape ``(N, 2)`` array of ``(x, y)`` positions.

    Returns
    -------
    centered : numpy.ndarray
        ``locations`` shifted so the bounding-box centre is at the origin.
    center : tuple of float
        The ``(x, y)`` centre that was subtracted.
    """
    X, Y = locations[:, 0], locations[:, 1]
    center = ((np.max(X) + np.min(X)) / 2, (np.max(Y) + np.min(Y)) / 2)
    X, Y = X - center[0], Y - center[1]
    return np.array([[x, y] for x, y in zip(X, Y)]), center


def rotate_locs(locs, azi=0):
    """Rotate ``(N, 2)`` locations about their bounding-box centre.

    Parameters
    ----------
    locs : numpy.ndarray
        Shape ``(N, 2)`` array of ``(x, y)`` positions.
    azi : float, optional
        Rotation angle in degrees, positive clockwise from north.  ``0``
        (default) returns ``locs`` unchanged.

    Returns
    -------
    numpy.ndarray
        Rotated locations, re-referenced to the original centre.
    """
    if azi == 0:
        return locs
    azi = np.deg2rad(azi)
    locs, center = center_locs(locs)
    R = np.array([[np.cos(azi), -np.sin(azi)], [np.sin(azi), np.cos(azi)]])
    # +ve azimuth moves clockwise from North
    rot_locs = np.matmul(locs, R.T)
    rot_locs += center
    return rot_locs


def project_to_line(locs, azi):
    """Orthogonally project points onto a line through the origin.

    Parameters
    ----------
    locs : numpy.ndarray
        Shape ``(N, 2)`` array of ``(x, y)`` positions.
    azi : float
        Azimuth of the line in degrees.

    Returns
    -------
    numpy.ndarray
        Shape ``(N, 2)`` array of the projected points.
    """
    a = np.array([0, 0])
    b = np.array([np.tan(np.deg2rad(azi)), 1])
    projected_locs = np.zeros(locs.shape)
    for ii in range(len(locs)):
        ap = locs[ii, :] - a
        ab = b - a
        projected_locs[ii, :] = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
    return projected_locs


def linear_distance(x, y):
    """Cumulative along-profile distance, snapped to the nearest node.

    Used to turn a set of 2-D site positions into a 1-D distance coordinate for
    pseudo-sections.

    Parameters
    ----------
    x, y : numpy.ndarray
        1-D arrays of site coordinates, in profile order.

    Returns
    -------
    numpy.ndarray
        For each site, the cumulative distance along the polyline to its
        nearest node.
    """
    nodes = np.array([x, y]).T
    linear_x = np.zeros(x.shape)
    linear_x[1:] = np.sqrt((x[1:] - x[:-1]) ** 2 +
                           (y[1:] - y[:-1]) ** 2)
    linear_x = np.cumsum(linear_x)
    linear_site = np.ones((x.shape))
    for ii, (x1, y1) in enumerate(zip(x, y)):
        dist = np.sum((nodes - np.array([x1, y1]).T) ** 2, axis=1)
        idx = np.argmin(dist)
        linear_site[ii] = linear_x[idx]
    return linear_site


# --------------------------------------------------------------------------- #
# Coordinate conversion / map projection
# --------------------------------------------------------------------------- #
#: Cache of :class:`pyproj.Proj` objects keyed by UTM zone number.
_projections = {}


def parse_dms(dms):
    """Split a ``"d:m:s"`` string into three floats.

    Parameters
    ----------
    dms : str
        Degrees, minutes and seconds separated by ``':'``.

    Returns
    -------
    tuple of float
        ``(degrees, minutes, seconds)``.

    Raises
    ------
    ValueError
        If any field is not numeric.
    """
    d, m, s = [float(x) for x in dms.split(':')]
    return d, m, s


def dms2dd(dms):
    """Convert a degrees:minutes:seconds string to decimal degrees.

    Parameters
    ----------
    dms : str
        ``"d:m:s"`` string.  A plain decimal-degree string is also accepted and
        returned as a float.

    Returns
    -------
    float
        Decimal degrees, sign taken from the degrees field.
    """
    try:
        d, m, s = parse_dms(dms)
        dd = abs(d) + abs(m) / 60 + abs(s) / 3600
        val = dd * np.sign(d)
    except ValueError:  # already decimal degrees
        return float(dms)
    return val


def dd2dms(dd):
    """Convert decimal degrees to a ``(degrees, minutes, seconds)`` tuple.

    Parameters
    ----------
    dd : float
        Decimal degrees.

    Returns
    -------
    tuple of float
        ``(degrees, minutes, seconds)``.  The sign is carried on the degrees
        component.
    """
    mult = -1 if dd < 0 else 1
    mnt, sec = divmod(abs(dd) * 3600, 60)
    deg, mnt = divmod(mnt, 60)
    return mult * deg, mnt, sec


# Backwards-compatible aliases for the older, near-duplicate implementations.
# ``dd_to_dms`` was byte-for-byte identical to :func:`dd2dms`.  ``dms_to_dd``
# differed only in its (buggy) handling of negative and zero values and had no
# in-tree callers; :func:`dms2dd` is the maintained version.
dms_to_dd = dms2dd
dd_to_dms = dd2dms


def compute_longorigin(Long=0):
    """Central meridian (origin longitude) of the UTM zone containing ``Long``.

    Parameters
    ----------
    Long : float, optional
        Longitude in degrees.  Default 0.

    Returns
    -------
    float
        The zone's central meridian longitude.
    """
    if Long == 0:
        return -3.0
    sixths = int(Long / 6)
    if sixths == 0:  # |Long| < 6
        return 3.0 if Long > 0 else -3.0
    return float(sixths * 6) + 3.0 * (1 if sixths > 0 else -1)


def geo2utm(Lat, Long, long_origin=999):
    """Convert geographic coordinates to UTM (WGS84), long-hand.

    A self-contained UTM forward transform, kept because :mod:`pyMT.IO` relies
    on its exact numerical behaviour.  For new code prefer :func:`project`.

    Parameters
    ----------
    Lat, Long : float
        Latitude and longitude in degrees.
    long_origin : float, optional
        Central meridian to use.  The sentinel ``999`` (default) means "compute
        it from ``Long``" via :func:`compute_longorigin`.

    Returns
    -------
    Easting, Northing : float
        UTM coordinates in metres.
    origin : float
        The central meridian that was used.
    """
    AA = 6378137
    ECC2 = 0.00669438
    K0 = 0.9996
    if long_origin == 999:
        origin = compute_longorigin(Long)
    else:
        origin = long_origin
    LatR, LongR = np.deg2rad([Lat, Long])
    eccPrimeSquared = ECC2 / (1 - ECC2)
    N = AA / (np.sqrt(1 - ECC2 * np.sin(LatR) ** 2.0))
    T = np.tan(LatR) ** 2
    C = eccPrimeSquared * np.cos(LatR) ** 2
    A = np.cos(LatR) * np.deg2rad((Long - origin))
    M = AA * ((1 -
               ECC2 / 4 -
               3 * ECC2 * ECC2 / 64 -
               5 * ECC2 * ECC2 * ECC2 / 256) * LatR -
              (3 * ECC2 / 8 +
               3 * ECC2 * ECC2 / 32 +
               45 * ECC2 * ECC2 * ECC2 / 1024) * np.sin(2 * LatR) +
              (15 * ECC2 * ECC2 / 256 + 45 * ECC2 * ECC2 * ECC2 / 1024) * np.sin(4 * LatR) -
              (35 * ECC2 * ECC2 * ECC2 / 3072) * np.sin(6 * LatR))
    Easting = K0 * N * \
        (A + (1 - T + C) * ((A ** 3) / 6) +
         (5 - 18 * T ** 3 + 72 * C - 58 * eccPrimeSquared) *
            (A ** 5) / 120) + 500000

    Northing = K0 * \
        (M + N * np.tan(LatR) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * ((A ** 4) / 24)) +
         (61 - 58 * T + T * T + 600 * C -
          330 * eccPrimeSquared) * ((A ** 6) / 720))
    return Easting, Northing, origin


def zone_convert(coordinates):
    """UTM zone number for a ``(longitude, latitude)`` pair.

    Includes the south-Norway and Svalbard zone exceptions.

    Parameters
    ----------
    coordinates : sequence of float
        ``(longitude, latitude)`` in degrees.

    Returns
    -------
    int
        UTM zone number.
    """
    if 56 <= coordinates[1] < 64 and 3 <= coordinates[0] < 12:
        return 32
    if 72 <= coordinates[1] < 84 and 0 <= coordinates[0] < 42:
        if coordinates[0] < 9:
            return 31
        elif coordinates[0] < 21:
            return 33
        elif coordinates[0] < 33:
            return 35
        return 37
    return int((coordinates[0] + 180) / 6) + 1


def generic_projection(x, y, from_proj=None, to_proj=None, utm_zone=None):
    """Transform coordinate arrays between named CRS shortcuts.

    Parameters
    ----------
    x, y : sequence of float
        Input coordinates.
    from_proj, to_proj : {'lambert', 'latlong', 'utm'}
        Source and target CRS.  ``'lambert'`` is EPSG:3979, ``'latlong'`` is
        EPSG:4326, ``'utm'`` is EPSG:4269.  Both must be given.
    utm_zone : object, optional
        Currently unused; reserved for a future zone-aware UTM option.

    Returns
    -------
    tuple of numpy.ndarray
        ``(out_x, out_y)``.  If ``from_proj``/``to_proj`` are missing this
        returns ``False``; if a name is unrecognised the inputs are returned
        unchanged.
    """
    if not (from_proj and to_proj):
        print('Have to define your projections...')
        return False
    proj_dict = {'lambert': 'epsg:3979',
                 'latlong': 'epsg:4326',
                 'utm': 'epsg:4269'}
    try:
        transformer = pyproj.Transformer.from_crs(proj_dict[from_proj], proj_dict[to_proj])
        out_x, out_y = np.zeros(len(x)), np.zeros(len(y))
        for ii, (xx, yy) in enumerate(zip(x, y)):
            out_x[ii], out_y[ii] = transformer.transform(xx, yy)
    except KeyError:
        print('Requested projection not recognized')
        return x, y
    return out_x, out_y


def project(coordinates, zone=None, letter=None):
    """Project a ``(longitude, latitude)`` pair to UTM using :mod:`pyproj`.

    Parameters
    ----------
    coordinates : sequence of float
        ``(longitude, latitude)`` in degrees.
    zone : int, optional
        UTM zone.  Auto-detected via :func:`zone_convert` if not given.
    letter : str, optional
        MGRS latitude band letter.  Auto-detected if not given.

    Returns
    -------
    zone : int
    letter : str
    x, y : float
        UTM easting/northing in metres (northing shifted by 10,000,000 for
        southern-hemisphere results).
    """
    def letter_convert(coordinates):
        return 'CDEFGHJKLMNPQRSTUVWXX'[int((coordinates[1] + 80) / 8)]
    if not (zone and letter):
        z = zone_convert(coordinates)
        L = letter_convert(coordinates)
    else:
        z = zone
        L = letter
    if z not in _projections:
        _projections[z] = pyproj.Proj(proj='utm', zone=z, ellps='WGS84')
    x, y = _projections[z](coordinates[0], coordinates[1])
    if y < 0:
        y += 10000000
    return z, L, x, y


def unproject(z, l, x, y):
    """Inverse of :func:`project`: UTM back to ``(longitude, latitude)``.

    Parameters
    ----------
    z : int
        UTM zone.
    l : str
        MGRS latitude band letter.  Letters before ``'N'`` are treated as
        southern hemisphere (the 10,000,000 m northing offset is removed).
    x, y : float
        UTM easting/northing in metres.

    Returns
    -------
    tuple of float
        ``(longitude, latitude)`` in degrees.
    """
    if z not in _projections:
        _projections[z] = pyproj.Proj(proj='utm', zone=z, ellps='WGS84')
    if l < 'N':
        y -= 10000000
    lng, lat = _projections[z](x, y, inverse=True)
    return (lng, lat)


def to_lambert(x, y):
    """Transform WGS84 lat/long arrays to the Canada LCC projection (EPSG:3979).

    Parameters
    ----------
    x : numpy.ndarray
        Latitudes.
    y : numpy.ndarray
        Longitudes.

    Returns
    -------
    lam_x, lam_y : numpy.ndarray
        Projected coordinates.
    """
    transformer = pyproj.Transformer.from_crs('epsg:4326', 'epsg:3979')
    lam_x = np.zeros(x.shape)
    lam_y = np.zeros(y.shape)
    for ii, (lat, lon) in enumerate(zip(x, y)):
        lam_x[ii], lam_y[ii] = transformer.transform(lat, lon)
    return lam_x, lam_y


# --------------------------------------------------------------------------- #
# MT transfer-function computations
#
# Each of these takes a pyMT.data_structures.Site and reads its ``data`` and
# ``used_error`` dictionaries (keys like 'ZXYR', 'ZXYI', ...).
# --------------------------------------------------------------------------- #
def compute_gav(site):
    """Geometric average of the off-diagonal impedances.

    Parameters
    ----------
    site : pyMT.data_structures.Site
        Must provide ``ZXY`` and ``ZYX`` data and errors.

    Returns
    -------
    gav : numpy.ndarray
        ``sqrt(Zxy * Zyx)`` (complex).
    gav_err : numpy.ndarray
        Propagated error estimate.

    Raises
    ------
    KeyError
        If a required component is missing.
    """
    try:
        gav = np.sqrt((site.data['ZXYR'] - 1j * site.data['ZXYI']) *
                      (site.data['ZYXR'] - 1j * site.data['ZYXI']))
        gav_err = np.sqrt(abs((site.used_error['ZXYR'] + 1j * site.used_error['ZXYI']) *
                              (site.used_error['ZYXR'] + 1j * site.used_error['ZYXI'])))
    except KeyError as e:
        print('Missing component needed for computation')
        raise e
    else:
        return gav, gav_err


def compute_aav(site):
    """Arithmetic average of the off-diagonal impedances.

    Parameters
    ----------
    site : pyMT.data_structures.Site
        Must provide ``ZXY`` and ``ZYX`` data and errors.

    Returns
    -------
    aav : numpy.ndarray
        ``(|Zxy| + |Zyx|) / 2`` (built component-wise, complex).
    aav_err : numpy.ndarray
        Propagated error estimate.

    Raises
    ------
    KeyError
        If a required component is missing.
    """
    try:
        aav = ((abs(site.data['ZXYR']) + 1j * abs(site.data['ZXYI'])) +
               (abs(site.data['ZYXR']) + 1j * abs(site.data['ZYXI']))) / 2
        aav_err = ((abs(site.used_error['ZXYR']) + abs(site.used_error['ZXYI'])) +
                   (abs(site.used_error['ZYXR']) + abs(site.used_error['ZYXI']))) / 4
    except KeyError as e:
        print('Missing component needed for computation')
        raise e
    else:
        return aav, aav_err


def compute_MT_determinant(site):
    """Determinant (Berdichevsky) invariant of the impedance tensor.

    Parameters
    ----------
    site : pyMT.data_structures.Site
        Must provide all four impedance components and their errors.

    Returns
    -------
    det : numpy.ndarray
        ``sqrt(Zxy * Zyx - Zxx * Zyy)`` (complex).
    det_err : numpy.ndarray
        Propagated error estimate.

    Raises
    ------
    KeyError
        If any impedance component is missing.
    """
    try:
        det = np.sqrt((site.data['ZXYR'] - 1j * site.data['ZXYI']) *
                      (site.data['ZYXR'] - 1j * site.data['ZYXI']) -
                      (site.data['ZXXR'] - 1j * site.data['ZXXI']) *
                      (site.data['ZYYR'] - 1j * site.data['ZYYI']))
        det_err = (site.used_error['ZXXR'] * (abs(site.data['ZYYR'] + 1j*site.data['ZYYI'])) +
                   site.used_error['ZYYR'] * (abs(site.data['ZXXR'] + 1j*site.data['ZXXI'])) +
                   site.used_error['ZXYR'] * (abs(site.data['ZYXR'] + 1j*site.data['ZYXI'])) +
                   site.used_error['ZYXR'] * (abs(site.data['ZXYR'] + 1j*site.data['ZXYI']))) / (2 * abs(det))
    except KeyError as e:
        print('Determinant cannot be computed unless all impedance components are available')
        raise e
    else:
        return det, det_err


def compute_ssq(site):
    """Sum-of-squares invariant of the impedance tensor.

    Parameters
    ----------
    site : pyMT.data_structures.Site
        Must provide all four impedance components and their errors.

    Returns
    -------
    ssq : numpy.ndarray
        ``sqrt((Zxx^2 + Zxy^2 + Zyx^2 + Zyy^2) / 2)`` (complex).
    det_err : numpy.ndarray
        Propagated error estimate (same expression as the determinant error,
        scaled by ``ssq``).

    Raises
    ------
    KeyError
        If any impedance component is missing.
    """
    try:
        zxy = (site.data['ZXYR'] - 1j * site.data['ZXYI'])
        zyx = (site.data['ZYXR'] - 1j * site.data['ZYXI'])
        zxx = (site.data['ZXXR'] - 1j * site.data['ZXXI'])
        zyy = (site.data['ZYYR'] - 1j * site.data['ZYYI'])
        ssq = np.sqrt((zxx * zxx + zxy * zxy + zyx * zyx + zyy * zyy) / 2)
        det_err = (site.used_error['ZXXR'] * (abs(site.data['ZYYR'] + 1j*site.data['ZYYI'])) +
                   site.used_error['ZYYR'] * (abs(site.data['ZXXR'] + 1j*site.data['ZXXI'])) +
                   site.used_error['ZXYR'] * (abs(site.data['ZYXR'] + 1j*site.data['ZYXI'])) +
                   site.used_error['ZYXR'] * (abs(site.data['ZXYR'] + 1j*site.data['ZXYI']))) / (2 * abs(ssq))
    except KeyError as e:
        print('SSQ cannot be computed unless all impedance components are available')
        raise e
    else:
        return ssq, det_err


def compute_rho(site, calc_comp=None, errtype='none'):
    """Apparent resistivity from impedance data for a given site.

    Parameters
    ----------
    site : pyMT.data_structures.Site
        Site holding the impedance data.
    calc_comp : str, optional
        Which apparent resistivity to compute.  One of
        ``{'xx', 'xy', 'yy', 'yx', 'det', 'gav', 'aav', 'ssq'}`` (default
        ``'det'``).  A leading ``'rho'`` is stripped, so ``'rhoxy'`` also works.
    errtype : str, optional
        Name of the :class:`~pyMT.data_structures.Site` error attribute to
        propagate (e.g. ``'used_error'``), or ``'none'`` (default) for no error.

    Returns
    -------
    rho_data : numpy.ndarray
        Apparent resistivity.
    rho_error : numpy.ndarray
        Linear apparent-resistivity error (zeros if ``errtype == 'none'``).
    rho_log10Err : numpy.ndarray
        Error in ``log10`` apparent resistivity.

    Notes
    -----
    Remove-flag entries in the source errors are preserved through to the
    output error arrays.
    """
    def compute_rho_error(site, errtyped, compd):
        zeR = getattr(site, errtyped)[''.join([compd, 'R'])]
        zeI = getattr(site, errtyped)[''.join([compd, 'I'])]
        z = site.data[''.join([compd, 'R'])] - 1j * site.data[''.join([compd, 'I'])]
        if len(zeR) == 0:
            zeR = zeI = z * 0
        C = 2 * site.periods / (MU * 2 * np.pi)
        rho_error = np.sqrt((C * np.real(z) * zeR) ** 2 + (C * np.imag(z) * zeI) ** 2)
        rho_log10Err = (rho_error * np.sqrt(2) /
                        ((C * np.real(z) ** 2 + C * np.imag(z) ** 2) * np.log(10)))
        idx = zeR == site.REMOVE_FLAG
        # Preserve remove flags
        rho_error[idx] = site.REMOVE_FLAG
        rho_log10Err[idx] = site.REMOVE_FLAG
        return rho_error, rho_log10Err

    COMPS = ('XX', 'XY',
             'YY', 'YX',
             'DET', 'GAV', 'AAV', 'SSQ')
    if not calc_comp:
        calc_comp = 'det'
    calc_comp = calc_comp.lower()
    if 'rho' in calc_comp:
        calc_comp = calc_comp.replace('rho', '')
    if calc_comp.upper() not in COMPS:
        print('{} not a valid option'.format(calc_comp))
        return
    if calc_comp.upper() not in COMPS[:4]:
        if calc_comp.lower() == 'det':
            det, err = compute_MT_determinant(site)
        elif calc_comp == 'gav':
            det, err = compute_gav(site)
        elif calc_comp == 'aav':
            det, err = compute_aav(site)
        elif calc_comp == 'ssq':
            det, err = compute_ssq(site)
        rho_data = det * np.conj(det) * site.periods / (MU * 2 * np.pi)
        if errtype.lower() != 'none':
            eXY, eXY_log10 = compute_rho_error(deepcopy(site), errtype, 'ZXY')
            eYX, eYX_log10 = compute_rho_error(site, errtype, 'ZYX')
            rho_error = np.max((eXY, eYX), axis=0)
            rho_log10Err = np.max((eXY_log10, eYX_log10), axis=0)
        else:
            rho_error = rho_data * 0
            rho_log10Err = rho_data * 0
    else:
        comp = ''.join(['Z', calc_comp.upper()])
        z = site.data[''.join([comp, 'R'])] - 1j * site.data[''.join([comp, 'I'])]
        rho_data = z * np.conj(z) * site.periods / (MU * 2 * np.pi)
        if errtype.lower() != 'none':
            rho_error, rho_log10Err = compute_rho_error(site, errtype, comp)
        else:
            rho_error = rho_log10Err = rho_data * 0
    return np.real(rho_data), np.real(rho_error), np.real(rho_log10Err)


def compute_phase(site, calc_comp=None, errtype=None, wrap=0):
    """Impedance phase for a given site.

    Parameters
    ----------
    site : pyMT.data_structures.Site
        Site holding the impedance data.
    calc_comp : str, optional
        Which phase to compute.  One of
        ``{'xx', 'xy', 'yy', 'yx', 'det', 'gav', 'aav', 'ssq'}`` (default
        ``'det'``).  A leading ``'pha'`` is stripped.
    errtype : str, optional
        Name of the :class:`~pyMT.data_structures.Site` error attribute to
        propagate, or ``None``/``'none'`` (default) for no error.
    wrap : int, optional
        If truthy, add 180 degrees to a ``'yx'`` phase so it plots in the first
        quadrant.  Default 0.

    Returns
    -------
    pha_data : numpy.ndarray
        Phase in degrees.
    pha_error : numpy.ndarray
        Phase error in degrees (zeros if no ``errtype``).

    Notes
    -----
    The phase error is approximated by averaging the phase over the eight
    combinations of +/- real and imaginary impedance error; this over-estimates
    the error when impedance errors are large.  Remove-flag entries are
    preserved.
    """
    def compute_phase_error(phase_data, z_real, z_imag, r_error, im_err):
        zE = np.zeros((len(z_real), 8))
        zE[:, 0] = np.angle(z_real - r_error + (z_imag + im_err), deg=True)
        zE[:, 1] = np.angle(z_real + r_error + (z_imag + im_err), deg=True)
        zE[:, 2] = np.angle(z_real + r_error + (z_imag - im_err), deg=True)
        zE[:, 3] = np.angle(z_real - r_error + (z_imag - im_err), deg=True)
        zE[:, 4] = np.angle(z_real + (z_imag - im_err), deg=True)
        zE[:, 5] = np.angle(z_real + (z_imag + im_err), deg=True)
        zE[:, 6] = np.angle(z_real + r_error + z_imag, deg=True)
        zE[:, 7] = np.angle(z_real - r_error + z_imag, deg=True)
        pha_error = np.mean(abs(phase_data[:, np.newaxis] - zE), 1)
        # Preserve remove flags
        idx = r_error == site.REMOVE_FLAG
        pha_error[idx] = site.REMOVE_FLAG
        return pha_error

    COMPS = ('XX', 'XY',
             'YY', 'YX',
             'DET', 'GAV', 'AAV', 'SSQ')
    if not errtype:
        errtype = 'none'
    if not calc_comp:
        calc_comp = 'det'
    calc_comp = calc_comp.lower()
    if 'pha' in calc_comp:
        calc_comp = calc_comp.replace('pha', '')
    if calc_comp.upper() not in COMPS:
        print('{} not a valid option'.format(calc_comp))
        return
    if calc_comp.upper() not in COMPS[:4]:
        if calc_comp == 'det':
            det, err = compute_MT_determinant(site)
        elif calc_comp == 'gav':
            det, err = compute_gav(site)
        elif calc_comp == 'aav':
            det, err = compute_aav(site)
        elif calc_comp == 'ssq':
            det, err = compute_ssq(site)
        pha_data = np.rad2deg(np.arctan2(np.imag(det), np.real(det)))
        if calc_comp == 'aav' or calc_comp == 'ssq':
            pha_data = pha_data % 90
        else:
            pha_data += 90

        if errtype.lower() != 'none':
            zR = site.data['ZXYR']
            zI = -1j*site.data['ZXYI']
            z = zR + zI
            zeR = getattr(site, errtype)['ZXYR']
            zeI = -1j*getattr(site, errtype)['ZXYI']
            eXY = compute_phase_error(np.angle(z, deg=True), zR, zI, zeR, zeI)
            zR = site.data['ZYXR']
            zI = -1j*site.data['ZYXI']
            z = zR + zI
            zeR = getattr(site, errtype)['ZYXR']
            zeI = -1j*getattr(site, errtype)['ZYXI']
            eYX = compute_phase_error(np.angle(z, deg=True), zR, zI, zeR, zeI)
            pha_error = np.max((eXY, eYX), axis=0)
        else:
            pha_error = 0 * pha_data
    else:
        comp = ''.join(['Z', calc_comp.upper()])
        zR = site.data[''.join([comp, 'R'])]
        zI = -1j * site.data[''.join([comp, 'I'])]
        z = zR + zI
        pha_data = np.angle(z, deg=True)
        pha_error = pha_data * 0
        if errtype.lower() != 'none':
            zeR = getattr(site, errtype)[''.join([comp, 'R'])]
            zeI = -1j * getattr(site, errtype)[''.join([comp, 'I'])]
            pha_error = compute_phase_error(pha_data, zR, zI, zeR, zeI)
        if 'yx' in calc_comp and wrap:
            pha_data += 180
    return pha_data, pha_error


def compute_bost1D(site, method='phase', comp=None, filter_width=1):
    """1-D Niblett-Bostick depth transform.

    Parameters
    ----------
    site : pyMT.data_structures.Site
        Site holding the impedance data.
    method : {'phase', 'bostick', 'phase-bostick'}, optional
        Transform variant.  Default ``'phase'``.
    comp : str, optional
        Component/invariant passed to :func:`compute_rho` / :func:`compute_phase`
        (default ``'det'``).  A leading ``'bost'`` is stripped.
    filter_width : float, optional
        Width passed to :func:`geotools_filter`.  Default 1.

    Returns
    -------
    bostick : numpy.ndarray
        Bostick resistivity.
    depth : numpy.ndarray
        Depth in km.
    rhofit : numpy.ndarray
        Smoothed apparent resistivity used in the transform.
    phase : numpy.ndarray
        Phase used in the transform.

    Notes
    -----
    Returns ``([0], [0], [0], [0])`` if the site has no valid apparent
    resistivity points.
    """
    if not comp:
        comp = 'det'
    if 'bost' in comp.lower():
        comp = comp.lower().replace('bost', '')
    rho = compute_rho(site, calc_comp=comp)[0]
    periods = site.periods
    if len(rho) != 0:
        rhofit = geotools_filter(periods, rho, fwidth=filter_width)
    else:
        print('No valid data points for bost1D.')
        return [0], [0], [0], [0]
    log_rho = np.log10(rhofit)
    log_freq = np.log10(1 / periods)
    depth = np.sqrt(rhofit * periods / (2 * np.pi * MU)) / 1000
    if method.lower() == 'bostick':
        slope = np.diff(log_rho) / abs(np.diff(log_freq))
        slope = list(slope)
        slope.append((log_rho[-1] - log_rho[-2]) / abs(log_freq[-1] - log_freq[-2]))
        slope = np.array(slope)
        slope[abs(slope) > 0.95] = np.sign(slope[abs(slope) > 0.95]) * 0.95
        bostick = rhofit * ((1 + slope) / (1 - slope))
        phase = compute_phase(site, calc_comp=comp)[0]
    elif method.lower() == 'phase-bostick':
        phase = compute_phase(site, calc_comp=comp)[0]
        phase = geotools_filter(periods, phase, use_log=False, fwidth=filter_width)
        slope = 1 - phase / 45
        bostick = rho * ((1 + slope) / (1 - slope))
    elif method.lower() == 'phase':
        phase = compute_phase(site, calc_comp=comp)[0]
        phase = geotools_filter(periods, phase, use_log=False, fwidth=filter_width)
        bostick = rho * ((np.pi / (2 * np.deg2rad(phase % 90))) - 1)
    return bostick, depth, rhofit, phase


def convert2impedance(rho, phase, periods, component):
    """Inverse of :func:`compute_rho` / :func:`compute_phase`.

    Parameters
    ----------
    rho : numpy.ndarray
        Apparent resistivity.
    phase : numpy.ndarray
        Phase in degrees.
    periods : numpy.ndarray
        Periods in seconds.
    component : str
        Impedance component name; only whether it contains ``'xy'`` or ``'yx'``
        matters (it sets the sign convention).

    Returns
    -------
    re, im : numpy.ndarray
        Real and imaginary impedance.
    """
    phase = np.deg2rad(phase)
    tp = np.tan(phase)
    re = np.sign(tp) * np.sqrt(5 * rho * MU / (periods * (1 + tp**2)))
    im = np.sign(tp) * np.sqrt(5 * rho * MU / (periods * (1 + 1 / (tp**2))))
    if 'xy' in component.lower():
        im = -im
    elif 'yx' in component.lower():
        re = -re
    return re, im


def calculate_misfit(data_site, response_site):
    """Per-component and per-period misfit between data and a model response.

    Parameters
    ----------
    data_site : pyMT.data_structures.Site
        The observed data (provides ``used_error``).
    response_site : pyMT.data_structures.Site
        The model response.  Its ``components`` define which are compared.

    Returns
    -------
    misfit : dict
        ``{component: array}`` of squared, error-normalised residuals per period.
    comp_misfit : dict
        The same, plus a ``'Total'`` key holding the per-period mean over the
        ``NR`` response components.

    Notes
    -----
    Returns ``None`` if the data and response have different numbers of periods.
    """
    components = response_site.components
    NP = len(data_site.periods)
    NR = len(response_site.components)
    misfit = {comp: np.zeros((NP)) for comp in components}
    comp_misfit = {comp: np.zeros((NP)) for comp in components}
    period_misfit = np.zeros((NP))
    if NP != len(response_site.periods):
        print('Number of periods in data and response not equal!')
        return
    for comp in components:
        if comp in data_site.components:
            try:
                misfit[comp] = (np.abs(response_site.data[comp] - data_site.data[comp]) /
                                data_site.used_error[comp]) ** 2
                comp_misfit[comp] = ((misfit[comp]))
                period_misfit += misfit[comp]
            except ValueError:
                # data/response shapes disagree for this component; skip it.
                continue
    period_misfit = (period_misfit / NR)
    comp_misfit.update({'Total': period_misfit})
    return misfit, comp_misfit


# --------------------------------------------------------------------------- #
# Mesh / model
# --------------------------------------------------------------------------- #
def generate_zmesh(min_depth=1, max_depth=500000, NZ=None):
    """Generate a logarithmically spaced vertical mesh.

    .. deprecated::
        Superseded by :meth:`pyMT.data_structures.Model.generate_zmesh`.  Kept
        for :mod:`pyMT.scripts.j2ws3d`.

    Parameters
    ----------
    min_depth, max_depth : float, optional
        Shallowest and deepest depth (m).  Defaults 1 and 500000.
    NZ : int or sequence of int
        Either a single layer count, or one count per decade between
        ``min_depth`` and ``max_depth``.

    Returns
    -------
    depths : list of float
        Layer boundary depths, with 0 prepended.
    zCS : list of float
        Layer thicknesses.
    ddz : numpy.ndarray
        Second difference of ``depths`` (a smoothness diagnostic).
    """
    print('This function is deprecated.\n')
    print('It functions as it used to, however, the function has been updated and included as a method within the pyMT Model object.')
    num_decade = int(np.ceil(np.log10(max_depth)) - floor(np.log10(min_depth)))
    try:
        if num_decade != len(NZ):
            print('Number of decades is: {}'.format(num_decade))
            print('Number of parameters (NZ) given is: {}'.format(len(NZ)))
            print('Length of NZ must be equal to the number of decades.')
            return
        decade = np.log10(min_depth)
        depths = []
        for n in NZ:
            dDecade = np.logspace(decade, min(np.floor(decade + 1), np.log10(max_depth)), int(n + 1))
            decade = floor(decade + 1)
            depths.append(dDecade)
        depths = flatten_list(depths)
        depths = np.array(np.unique(depths))
    except TypeError:
        depths = np.logspace(np.log10(min_depth), np.log10(max_depth), int(NZ))
    ddz = np.diff(np.diff(depths))
    if any(ddz < 0):
        print('Warning! Second derivative of depths is not always positive.')
    zCS = np.diff(depths)
    zCS = list(zCS)
    depths = list(depths)
    depths.insert(0, 0)
    zCS.insert(0, depths[1])
    return depths, zCS, ddz


def generate_lateral_mesh(site_locs, regular=True, min_x=None, model=None, max_x=None,
                          num_pads=None, pad_mult=None, DEBUG=True):
    """Generate a 1-D lateral (x or y) mesh from a set of station locations.

    Parameters
    ----------
    site_locs : numpy.ndarray
        1-D array of station coordinates along the axis being meshed.
    regular : bool, optional
        If ``True`` (default), build an evenly spaced mesh of pitch ``min_x``
        spanning the stations plus a 500 m margin.  If ``False``, use the
        station-aware "J2" mesher.
    min_x : float, optional
        Minimum / target cell size.
    max_x : float, optional
        Maximum gap allowed without a mesh node (J2 mesher only).
    model, num_pads, pad_mult : optional
        Currently unused; padding is applied elsewhere.  Retained for API
        compatibility.
    DEBUG : bool, optional
        Print progress from the J2 mesher.  Default ``True``.

    Returns
    -------
    xmesh : list or numpy.ndarray
        Mesh node positions.
    nmesh : int
        Number of nodes.
    """
    def regular_mesh(x_size, bounds):
        print(bounds, x_size)
        xmesh = list(np.arange(bounds[0], bounds[1], x_size))
        nmesh = len(xmesh)
        return xmesh, nmesh

    def j2_mesh(site_locs, min_x, max_x):
        MAX_X = 150
        xlocs = site_locs  # South-North
        xloc_sort = np.sort(xlocs)
        xoff = xloc_sort[0] + (xloc_sort[-1] - xloc_sort[0]) / 2
        xloc_sort = (xloc_sort - xoff).astype(int)
        x_seperation = np.diff(xloc_sort)
        max_sep = np.max(x_seperation)
        xmesh = np.zeros(MAX_X)
        avg_sep = np.mean(x_seperation)
        max_xmin = 3 * avg_sep
        if not min_x:
            min_x = max_xmin / 2
        if not max_x:
            max_x = min_x * 2
        if min_x > max_xmin:
            print('Minimum cell size shouldn\'t be more than {}'.format(max_xmin))
        dist_without_mesh = 0
        is_right_ofmesh = 0
        imesh = 1
        xmesh[0] = xloc_sort[0] - min_x * 0.5
        for ii in range(len(xloc_sort) - 1):
            # Look to station location on right
            dist_right = xloc_sort[ii + 1] - xloc_sort[ii]
            ifact = 0
            dist = 0
            while dist < dist_right:
                ifact += 1
                dist += min_x + min_x * ifact
            ifact -= ifact
            if ifact >= 2:
                for jj in range(1, max(ifact - 2, 1) + 1):
                    imesh += 1
                    xmesh[imesh - 1] = xmesh[imesh - 2] + min_x * jj
                imesh += 1
                xmesh[imesh - 1] = xloc_sort[ii] + (xloc_sort[ii + 1] - xloc_sort[ii]) / 2
                imesh += ifact - 1
                imesh_save = imesh
                xmesh[imesh - 1] = xloc_sort[ii + 1] - min_x / 2
                for jj in range(1, ifact - 1):
                    imesh -= 1
                    xmesh[imesh - 1] = xmesh[imesh] - min_x * jj
                imesh = imesh_save
                dist_without_mesh = 0
                is_right_ofmesh = ii + 1
            elif (dist_right >= min_x) and (dist_right <= 5 * min_x):
                if DEBUG:
                    print('Gap is small, splitting sites')
                imesh += 1
                xmesh[imesh - 1] = xloc_sort[ii] + (xloc_sort[ii + 1] - xloc_sort[ii]) / 2
                dist_without_mesh = 0
                is_right_ofmesh = ii + 1
            else:
                # no room for mesh?
                if DEBUG:
                    print('Sites too close, splitting mesh')
                dist_without_mesh = xloc_sort[ii + 1] - xmesh[imesh - 1]
            if dist_without_mesh > max_x:
                if DEBUG:
                    print("Gone too far without adding mesh")
                max_sep = 0
                for is_check in range(is_right_ofmesh, ii + 1):
                    check_sep = xloc_sort[is_check + 1] - xloc_sort[is_check]
                    if check_sep > max_sep:
                        is_save = is_check
                        max_sep = check_sep
                if imesh + 1 > MAX_X:
                    resp = input('Number of cells exceeds {}. Continue? (y/n)'.format(MAX_X))
                    if resp == 'n':
                        return
                imesh += 1
                xmesh[imesh - 1] = xloc_sort[is_save] + (xloc_sort[is_save + 1] - xloc_sort[is_save]) / 2
                dist_without_mesh = 0
                is_right_ofmesh = is_save + 1
        imesh += 1
        if imesh > MAX_X:
            resp = input('Number of cells exceeds {}. Continue? (y/n)'.format(MAX_X))
            if resp == 'n':
                return
        xmesh[imesh - 1] = xloc_sort[-1] + xloc_sort[-1] - xmesh[imesh - 2]
        nmesh = imesh
        return xmesh[:nmesh], nmesh

    if regular:
        min_x = max(min_x, 10)
        bounds = [np.min(site_locs) - 500, np.max(site_locs) + 500]
        xmesh, nmesh = regular_mesh(min_x, bounds)
        if nmesh < 2:
            xmesh.append(xmesh[0] + 1)
        print(xmesh, nmesh)
    else:
        xmesh, nmesh = j2_mesh(site_locs=site_locs, min_x=min_x, max_x=max_x)
    return xmesh, nmesh


def regrid_model(mod, new_x, new_y, new_z, rho_axis=None):
    """Nearest-neighbour resample of a model onto a new mesh.

    Parameters
    ----------
    mod : pyMT.data_structures.Model
        Source model.  Provides ``dx``/``dy``/``dz`` edges, ``vals`` (or
        ``rho_x``/``rho_y``/``rho_z``) and ``background_resistivity``.
    new_x, new_y, new_z : array_like
        New mesh cell edges.
    rho_axis : {'rho_x', 'rho_y', 'rho_z', None}, optional
        Which value array to resample.  ``None`` (default) uses ``mod.vals``.

    Returns
    -------
    numpy.ndarray
        Resampled values with shape ``(len(new_x) - 1, len(new_y) - 1,
        len(new_z) - 1)``.  Query points outside the source model take
        ``mod.background_resistivity``.
    """
    if not rho_axis:
        vals = mod.vals
    elif rho_axis.lower() == 'rho_x':
        vals = mod.rho_x
    elif rho_axis.lower() == 'rho_y':
        vals = mod.rho_y
    elif rho_axis.lower() == 'rho_z':
        vals = mod.rho_z
    else:
        vals = mod.vals
    x, y, z = (edge2center(arr) for arr in (mod.dx, mod.dy, mod.dz))
    X, Y, Z = (edge2center(arr) for arr in (new_x, new_y, new_z))
    X_grid, Y_grid, Z_grid = np.meshgrid(X, Y, Z)
    X_grid, Y_grid, Z_grid = (np.ravel(arr) for arr in (X_grid, Y_grid, Z_grid))
    interp = RGI((y, x, z),
                 np.transpose(vals, [1, 0, 2]),
                 method='nearest', bounds_error=False, fill_value=mod.background_resistivity)
    query_points = np.array((Y_grid, X_grid, Z_grid)).T
    new_vals = interp(query_points)
    new_vals = np.reshape(new_vals, [len(Y), len(X), len(Z)])
    new_vals = np.transpose(new_vals, [1, 0, 2])
    return new_vals


# --------------------------------------------------------------------------- #
# Filesystem
# --------------------------------------------------------------------------- #
def check_file(outfile):
    """Return ``True`` if ``outfile`` exists and can be opened for reading.

    Parameters
    ----------
    outfile : str
        Path to test.

    Returns
    -------
    bool
    """
    try:
        open(outfile, 'r')
    except FileNotFoundError:
        return False
    else:
        return True


def fileparts(filename):
    """Split a path into directory, basename and extension (MATLAB-style).

    Parameters
    ----------
    filename : str

    Returns
    -------
    path : str
        Directory portion.
    filename : str
        Basename (with extension).
    ext : str
        Extension, including the leading dot.
    """
    path = os.path.dirname(filename)
    filename = os.path.basename(filename)
    ext = os.path.splitext(filename)[1]
    return path, filename, ext


def check_extension(outfile, expected=''):
    """Ensure ``outfile`` ends with the ``expected`` extension.

    Parameters
    ----------
    outfile : str
        Filename to check.
    expected : str, optional
        Extension (without a dot) that ``outfile`` should have.

    Returns
    -------
    str
        ``outfile`` unchanged if its extension already contains ``expected``,
        otherwise ``outfile`` with ``.expected`` appended.
    """
    path, file, ext = fileparts(outfile)
    if expected in ext:
        return outfile
    else:
        return '.'.join([outfile, expected])


#: Misspelled original name, kept as an alias for existing callers.
check_extention = check_extension


@enforce_input(files=list)
def sort_files(files):
    """Classify a list of files by pyMT role, keyed by extension/content.

    Parameters
    ----------
    files : list of str
        Paths to classify.

    Returns
    -------
    dict
        Maps role names (``'model'``, ``'data'``, ``'response'``, ``'list'``,
        ``'resolution'``) to the last matching path.  A ``.dat``/``.data`` file
        whose first line contains ``'response'`` is classed as ``'response'``.
    """
    ret_dict = {}
    for file in files:
        name, ext = os.path.splitext(file)
        if ext in ('.rho', '.model', '.zani'):
            ret_dict.update({'model': file})
        elif ext in ('.dat', '.data'):
            with open(file, 'r') as f:
                line = f.readline()
            if 'response' in line:
                ret_dict.update({'response': file})
            else:
                ret_dict.update({'data': file})
        elif ext in ('.lst', '.list'):
            ret_dict.update({'list': file})
        elif ext in ('.reso', '.resolution'):
            ret_dict.update({'resolution': file})
    return ret_dict
