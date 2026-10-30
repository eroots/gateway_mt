"""1-D magnetotelluric forward modelling and a simple 1-D inversion.

:func:`wait_1D` computes the impedance of a layered earth with Wait's
recursion. :func:`invert_1DMT` fits a smooth layered model to the average
apparent resistivity of a :class:`~gateway_mt.data_structures.Data` set,
using CMA-ES (the ``cma`` package) to minimise :func:`objective_function`.
The 1-D modelling window in the GUIs (:mod:`gateway_mt.GUI_common.windows`)
is built on these.

Impedances here are in ohms (E/H), with a positive imaginary part for a
half-space.
"""
import numpy as np
import cma
from gateway_mt.utils import compute_ssq, compute_MT_determinant


mu = 4 * np.pi * 1e-7

def compute_rho(Z, periods):
    """Apparent resistivity (ohm-m) of impedance ``Z`` (ohms) at ``periods`` (s)."""
    return np.real(Z * np.conj(Z) * periods / (mu * 2 * np.pi))

def compute_phase(Z):
    """Phase of ``Z`` in degrees."""
    return np.rad2deg(np.arctan2(np.imag(Z), np.real(Z)))

def layer_tops(dz):
    """Depth of the top of each layer, from the layer thicknesses ``dz``."""
    return np.concatenate([[0], np.cumsum(dz[:-1])])


def compute_det1D(Z):
    """Determinant impedance of a 1-D response, i.e. ``abs(Z)``.

    The determinant and SSQ invariants are equivalent in 1-D.
    """
    return np.sqrt((Z * np.conj(Z)))

def wait_1D(periods, model, dz):
    """Takes a model and periods and returns the 1-D impedance response.

    Parameters
    ----------
    periods : numpy.ndarray
        ``(NP,)`` periods (s) to generate the response at.
    model : array_like
        ``(M,)`` layer resistivities (ohm-m), top to bottom. The last layer
        is the basement half-space.
    dz : array_like
        ``(M,)`` layer thicknesses (m), one per layer. The last entry (the
        half-space) is not used, but must be present: with fewer than two
        entries the layer recursion is skipped, which is only right for a
        uniform half-space.

    Returns
    -------
    z_obs : numpy.ndarray
        ``(NP,)`` complex impedance (ohms).

    Examples
    --------
    A 1 km thick, 100 ohm-m layer over a 10 ohm-m half-space:

    >>> Z = wait_1D(np.logspace(-3, 3, 30), [100, 10], [1000, 0])
    """
    cond = 1 / np.array(model)
    omega = 2 * np.pi / periods
    z_obs = np.zeros(len(periods), dtype=np.complex128)
    for nfreq, w in enumerate(omega):
        prop_const = np.sqrt(1j*mu*cond[-1] * w)
        C = np.zeros(len(cond), dtype=np.complex128)
        C[-1] = 1 / prop_const
        if len(dz) > 1:
            for k in reversed(range(len(cond) - 1)):
                prop_layer = np.sqrt(1j*w*mu*cond[k])
                k1 = (C[k+1] * prop_layer + np.tanh(prop_layer * dz[k]))
                k2 = ((C[k+1] * prop_layer * np.tanh(prop_layer * dz[k])) + 1)
                C[k] = np.squeeze((1 / prop_layer) * (k1 / k2))
        z_obs[nfreq] = 1j * w * mu * C[0]

    return z_obs


def objective_function(model, *args):
    """Calculates the objective function for a given model and set of arguments.

    ``phi = |misfit - target_misfit| / regpar + roughness``, where ``misfit``
    is the RMS misfit of log10 apparent resistivity and ``roughness`` is
    ``sum(|m[i+1] - m[i]| ** model_norm)``.

    Parameters
    ----------
    model : numpy.ndarray
        log10 layer resistivities.
    *args
        ``(data_obs, err_obs, periods, dz, data_norm, model_norm, regpar,
        target_misfit, mantle_transition)``, as passed by
        :func:`invert_1DMT`. ``data_obs`` and ``err_obs`` are in log10
        apparent resistivity; ``data_norm`` and ``mantle_transition`` are
        not used here.

    Returns
    -------
    float
    """
    (data_obs, err_obs, periods, dz, data_norm, model_norm, regpar, target_misfit, mantle_transition) = args
    # rhoa_obs = compute_rho(compute_det1D(data_obs), periods)
    # model = 10**model
    # if mantle_transition:
    #     idx = np.argmin(410000 - np.array(dz))
    #     model[idx:] = np.log10(20)
    rhoa_mod = compute_rho(compute_det1D(wait_1D(periods, 10**model, dz)), periods)
    rhoa_mod = np.log10(rhoa_mod)
    # data_obs = np.log10(data_obs)
    # rhoa_mod = compute_rho(compute_det1D(z_mod), periods)
    # So far only rhoa modelling is done
    # pha_mod = compute_phase(compute_det1D(z_mod))
    reg_term = 0.
    for i in range(len(model)-1):
        reg_term += abs(model[i+1]-model[i])**model_norm

    misfit = np.sqrt(np.mean(((rhoa_mod - data_obs)**2) / err_obs**2))
    phi = (1/regpar) * np.abs(misfit - target_misfit) + reg_term

    return phi

def _invert_1DMT(rhoa_obs, err_obs, periods, dz, rho_initial,
                regpar=20, data_norm=2, model_norm=2, rho_min=1e-2,
                rho_max=1e6, maxiter=100, target_misfit=1.0, mantle_transition=False):

    x0 = np.log10(rho_initial)
    sigma0 = 0.2
    args = (rhoa_obs, err_obs, periods, dz, data_norm, model_norm, regpar, target_misfit, mantle_transition)
    options = {'bounds':[np.log10(rho_min), np.log10(rho_max)], 'maxiter':maxiter,
               'tolfunhist': 1e-2}
    if mantle_transition:
        # Fix the layers at and below the 410 km transition (and always the half-space) to 20 ohm-m
        below = np.nonzero(layer_tops(np.array(dz)) >= 410000)[0]
        idx = below[0] if len(below) else len(x0) - 1
        x0[idx:] = np.log10(20)
        options.update({'fixed_variables': {ii: np.log10(20) for ii in range(idx, len(x0))}})
    final_model, es = cma.fmin2(objective_function, x0, sigma0, options, args=args)
    final_model = 10**final_model
    z_best = wait_1D(periods, final_model, dz)
    rhoa_best = np.log10(compute_rho(compute_det1D(z_best), periods))
    rms_best = np.sqrt((np.abs(rhoa_best - rhoa_obs)**2 / err_obs**2).mean())

    return final_model, rhoa_best, z_best, rms_best

def invert_1DMT(data, rho_initial, dz,
                mode='ssq', uncertainty=0.1, regpar=0.1, data_norm=2,
                model_norm=2, rho_min=1e-2, rho_max=1e6,
                maxiter=100, target_misfit=1.0, mantle_transition=False):
    """Invert the average apparent resistivity of a data set for a 1-D model.

    Wraps ``_invert_1DMT`` to take a gateway_mt Data structure as its input:
    the data are the geometric mean over sites of the ``mode`` apparent
    resistivity (:meth:`~gateway_mt.data_structures.Data.average_rho`), with a
    uniform uncertainty in log10 units.

    Parameters
    ----------
    data : gateway_mt.data_structures.Data
        Data to invert.
    rho_initial : array_like
        Starting resistivity (ohm-m) for each layer.
    dz : array_like
        Layer thicknesses, as for :func:`wait_1D`.
    mode : str, optional
        Apparent resistivity to fit, passed to
        :func:`gateway_mt.utils.compute_rho` (default ``'ssq'``).
    uncertainty : float, optional
        Data uncertainty in log10(rho) (default 0.1).
    regpar : float, optional
        Trade-off parameter; smaller values weight the data fit more.
    data_norm : int, optional
        Not currently used (the misfit is always an RMS).
    model_norm : float, optional
        Exponent of the roughness term (default 2).
    rho_min, rho_max : float, optional
        Resistivity bounds (ohm-m).
    maxiter : int, optional
        Maximum CMA-ES iterations.
    target_misfit : float, optional
        RMS misfit to aim for (default 1.0).
    mantle_transition : bool, optional
        If True, fix the resistivity at 20 ohm-m in the layers whose tops are
        410 km or more down (and in the half-space).

    Returns
    -------
    final_model : numpy.ndarray
        Layer resistivities (ohm-m).
    rhoa_best : numpy.ndarray
        Apparent resistivity (ohm-m) of ``final_model``.
    z_best : numpy.ndarray
        Complex impedance (ohms) of ``final_model``.
    rms_best : float
        RMS misfit of ``rhoa_best``.
    """
    _, rhoa_obs = data.average_rho(comp=mode)
    rhoa_obs = np.log10(rhoa_obs)
    # rho_initial = np.array(list(rho_initial) + [hs])
    err_obs = np.ones(shape=rhoa_obs.shape) * uncertainty

    final_model, rhoa_best, z_best, rms_best = _invert_1DMT(rhoa_obs, err_obs, data.periods, dz, rho_initial,
                                                            regpar, data_norm, model_norm, rho_min,
                                                            rho_max, maxiter, target_misfit, mantle_transition)
    rhoa_best = 10**rhoa_best
    return final_model, rhoa_best, z_best, rms_best


