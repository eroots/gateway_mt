"""Small helpers for building test inputs."""
import numpy as np

import gateway_mt.data_structures as DS

MU0 = 4e-7 * np.pi
COMPONENTS = (('ZXX', (0, 0)), ('ZXY', (0, 1)), ('ZYX', (1, 0)), ('ZYY', (1, 1)))


def make_site(periods, Z, T=None, name='TEST', error=0.05):
    """A :class:`~gateway_mt.data_structures.Site` from complex impedances (and tipper).

    Parameters
    ----------
    periods : array (NP,)
    Z : complex array (NP, 2, 2)
        Impedance tensor (E/H, ohms) in the e^{+i w t} convention.
    T : complex array (NP, 2), optional
        Tipper (Tzx, Tzy).
    error : float
        Relative error on every impedance element, as a fraction of sqrt|Zxy Zyx|.
    """
    periods = np.asarray(periods, dtype=float)
    data, errors = {}, {}
    z_err = error * np.sqrt(np.abs(Z[:, 0, 1] * Z[:, 1, 0]))
    for comp, (i, j) in COMPONENTS:
        # gateway_mt stores the imaginary part with the opposite sign
        data[comp + 'R'] = np.real(Z[:, i, j]).copy()
        data[comp + 'I'] = -np.imag(Z[:, i, j]).copy()
        errors[comp + 'R'] = z_err.copy()
        errors[comp + 'I'] = z_err.copy()
    if T is not None:
        for k, comp in enumerate(('TZX', 'TZY')):
            data[comp + 'R'] = np.real(T[:, k]).copy()
            data[comp + 'I'] = -np.imag(T[:, k]).copy()
            errors[comp + 'R'] = np.full(len(periods), 0.02)
            errors[comp + 'I'] = np.full(len(periods), 0.02)
    return DS.Site(data=data, errors=errors, name=name, periods=periods.copy(), azimuth=0,
                   locations={'X': 0.0, 'Y': 0.0, 'Lat': 49.0, 'Long': -84.5, 'elev': 0.0})


def one_d_impedance(periods, rho, thicknesses):
    """Textbook layered-earth impedance (Ward & Hohmann recursion), independent of gateway_mt.

    ``thicknesses`` has one entry per layer except the bottom half-space.
    """
    rho = np.asarray(rho, dtype=float)
    out = []
    for period in np.asarray(periods, dtype=float):
        w = 2 * np.pi / period
        k = np.sqrt(1j * w * MU0 / rho)          # wavenumber in each layer
        intrinsic = 1j * w * MU0 / k              # intrinsic impedance of each layer
        Z = intrinsic[-1]
        for j in range(len(rho) - 2, -1, -1):
            t = np.tanh(k[j] * thicknesses[j])
            Z = intrinsic[j] * (Z + intrinsic[j] * t) / (intrinsic[j] + Z * t)
        out.append(Z)
    return np.array(out)


def one_d_tensor(Z):
    """The (NP, 2, 2) impedance tensor of a 1-D earth: Zxy = Z, Zyx = -Z, no diagonals."""
    tensor = np.zeros((len(Z), 2, 2), dtype=complex)
    tensor[:, 0, 1] = Z
    tensor[:, 1, 0] = -Z
    return tensor


def rotate_axes(Z, T, theta_deg):
    """Express an impedance tensor / tipper in axes rotated by ``theta_deg`` (one consistent rotation)."""
    th = np.deg2rad(theta_deg)
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, s], [-s, c]])
    Zr = np.einsum('ij,kjl,ml->kim', R, Z, R)
    Tr = np.einsum('kj,ij->ki', T, R) if T is not None else None
    return Zr, Tr
