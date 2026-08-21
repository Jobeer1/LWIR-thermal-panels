"""
wave_physics/analytic_benchmarks.py — Closed-form textbook limits.

These functions provide *exact* analytic results used to validate any
numerical wave module (FDTD, CMT) and the multi-layer Fresnel layer that feeds
the cached full-wave path.  Nothing here is approximate; each routine is a
direct evaluation of a textbook formula (Born & Wolf; Harrington).

Contents
--------
* Single-interface Fresnel reflection/transmission (TE & TM, arbitrary angle,
  complex media).                       -> fresnel_amplitudes()
* Multi-layer planar stack R/T/A via the transfer-matrix (ABCD) method at any
  angle, TE or TM.                      -> multilayer_stack_rt()
* Rectangular TE10 cutoff.               -> rectangular_cutoff_wavelength_um()
* Circular TE11 cutoff.                  -> circular_cutoff_wavelength_um()
* Evanescent axial decay constant.       -> evanescent_decay_constant_um()
* Evanescent *power* transmission.       -> evanescent_power_transmission()
* Geometric cavity enhancement.          -> geometric_cavity_enhancement()

All conventions (exp(-i*omega*t), n~ = n + i*k, time-averaged Poynting flux)
are imported from :mod:`wave_physics.conventions`.
"""

from __future__ import annotations

import cmath
import math
from typing import List, Sequence, Tuple

import numpy as np

from . import conventions

# Cutoff constants (textbook first-mode values).
RECT_CUTOFF_CONSTANT = 2.0          # TE10  : lambda_c = 2 * min(W, D)
# TE11 of a circular guide has k_c * a = 1.84118 (a = radius).  With
# a = D/2, lambda_c = 2*pi*(D/2)/1.84118 = (pi/1.84118)*D ~= 1.7062 * D.
CIRC_TE11_CONSTANT = math.pi / 1.84118  # ~= 1.7062


# ---------------------------------------------------------------------------
# Single-interface Fresnel coefficients (TE / TM)
# ---------------------------------------------------------------------------

def fresnel_amplitudes(n1: complex, n2: complex, theta1_rad: float, pol: str = 'TE'
                       ) -> Tuple[complex, complex]:
    """Complex Fresnel amplitude coefficients (r, t) at a single interface.

    Parameters
    ----------
    n1, n2       : complex refractive indices of incident / transmitted media.
    theta1_rad   : angle of incidence (radians) in medium 1.
    pol          : 'TE' (s) or 'TM' (p).

    Returns
    -------
    (r, t) complex amplitude coefficients.  Power follows as R = |r|^2 and
    T = (Re(eta2)/Re(eta1)) * |t|^2.
    """
    if pol.upper() not in ('TE', 'TM'):
        raise ValueError("pol must be 'TE' or 'TM'")

    n1c = complex(n1)
    n2c = complex(n2)

    sin2 = n1c * math.sin(theta1_rad) / n2c
    cos1 = cmath.sqrt(complex(1.0 - math.sin(theta1_rad) ** 2))
    cos2 = cmath.sqrt(complex(1.0) - sin2 * sin2)

    if pol.upper() == 'TE':
        num_r = n1c * cos1 - n2c * cos2
        den_r = n1c * cos1 + n2c * cos2
        r = num_r / den_r
        t = 2.0 * n1c * cos1 / den_r
    else:  # TM
        num_r = n2c * cos1 - n1c * cos2
        den_r = n2c * cos1 + n1c * cos2
        r = num_r / den_r
        t = 2.0 * n1c * cos1 / den_r
    return r, t


def single_interface_rt(n1: complex, n2: complex, theta1_rad: float, pol: str = 'TE'
                        ) -> Tuple[float, float, float]:
    """Power R, T, A for a single interface (semi-infinite media).

    ``A`` is defined via energy balance as ``1 - R - T``; for a lossless
    interface ``A`` is ~0 (all energy reflected/transmitted).
    """
    r, t = fresnel_amplitudes(n1, n2, theta1_rad, pol)
    R = abs(r) ** 2
    cos1r = math.cos(theta1_rad)
    cos2r = _cos2_real(n1, n2, theta1_rad)
    eta1 = n1.real * cos1r   if pol.upper() == 'TE' else n1.real / max(cos1r, 1e-12)
    eta2 = n2.real * cos2r   if pol.upper() == 'TE' else n2.real / max(cos2r, 1e-12)
    T = (eta2 / eta1) * abs(t) ** 2 if eta1 > 1e-15 else 0.0
    T = max(0.0, min(T, 1.0))
    A = max(0.0, 1.0 - R - T)
    return R, T, A


def _cos2_real(n1: complex, n2: complex, theta1_rad: float) -> float:
    """Real part of cos(theta2) in medium 2 via Snell from medium 1."""
    sin2 = n1 * math.sin(theta1_rad) / n2
    arg = 1.0 - (sin2.real ** 2 - sin2.imag ** 2)
    return max(0.0, math.sqrt(arg)) if arg >= 0 else 0.0


# ---------------------------------------------------------------------------
# Multi-layer transfer-matrix (ABCD) method, arbitrary angle
# ---------------------------------------------------------------------------

def _layer_admittance(n: complex, cos_theta: complex, pol: str) -> complex:
    """Normalized optical admittance eta of a layer for the given polarisation."""
    c = complex(cos_theta)
    if pol.upper() == 'TE':
        return complex(n) * c
    else:  # TM (p): eta = n / cos(theta); diverges at grazing -> clamp
        if abs(c) < 1e-12:
            return complex(n) / complex(1e-12)
        return complex(n) / c



def multilayer_stack_rt(
    indices: Sequence[complex],
    thicknesses_um: Sequence[float],
    wavelength_um: float,
    theta0_rad: float = 0.0,
    pol: str = 'TE',
) -> Tuple[float, float, float]:
    """Power (R, T, A) for a planar multilayer stack at incidence angle theta0_rad.

    Parameters
    ----------
    indices         : complex indices of every layer, *including* the incident
                      medium (first) and the semi-infinite substrate (last).
    thicknesses_um  : physical thickness of every layer *except* the incident
                      medium and substrate (len == len(indices) - 2), in µm.
    wavelength_um   : free-space wavelength (µm).
    theta0_rad      : incidence angle in the first medium (radians).
    pol             : 'TE' or 'TM'.

    Returns
    -------
    (R, T, A) power coefficients with A = 1 - R - T (energy balance).
    """
    if len(indices) < 2:
        raise ValueError("Need at least incident medium and substrate.")
    if len(thicknesses_um) != len(indices) - 2:
        raise ValueError(
            "thicknesses_um length must equal len(indices) - 2.")

    n0 = complex(indices[0])
    ns = complex(indices[-1])

    # Propagate the complex sine through the stack (Snell invariant).
    sin0 = n0 * math.sin(theta0_rad)
    cos_thetas: List[complex] = [cmath.sqrt(complex(1.0) - (sin0 / n0) ** 2)]
    for n in indices[1:]:
        sin_i = sin0 / complex(n)
        cos_thetas.append(cmath.sqrt(complex(1.0) - sin_i * sin_i))

    k0 = 2.0 * math.pi / wavelength_um

    # ABCD product over the finite layers.
    A, B, C, D = complex(1.0), complex(0.0), complex(0.0), complex(1.0)
    for i, (n_i, d_i) in enumerate(zip(indices[1:-1], thicknesses_um)):
        theta_i = cos_thetas[i + 1]  # cos_thetas[0] is the incident medium
        eta_i = _layer_admittance(n_i, theta_i, pol)
        delta = k0 * complex(n_i) * theta_i * d_i
        cos_d = cmath.cos(delta)
        sin_d = cmath.sin(delta)
        M = np.array([[A, B], [C, D]], dtype=complex)
        # Layer characteristic matrix.  With the exp(-i*w*t) convention the
        # off-diagonal elements carry a NEGATIVE imaginary unit.
        m = np.array([[cos_d, -1j * sin_d / eta_i],
                      [-1j * eta_i * sin_d, cos_d]], dtype=complex)
        M = M @ m
        A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]

    eta0 = _layer_admittance(n0, cos_thetas[0], pol)
    etas = _layer_admittance(ns, cos_thetas[-1], pol)

    denom = eta0 * A + eta0 * etas * B + C + etas * D
    r = (eta0 * A + eta0 * etas * B - C - etas * D) / denom
    t = 2.0 * eta0 / denom

    R = abs(r) ** 2
    re0 = max(eta0.real, 1e-15)
    res = max(etas.real, 0.0)
    T = (res / re0) * abs(t) ** 2 if re0 > 0 else 0.0
    T = float(np.clip(T, 0.0, 1.0))
    R = float(np.clip(R, 0.0, 1.0))
    A = float(np.clip(1.0 - R - T, 0.0, 1.0))
    return R, T, A


# ---------------------------------------------------------------------------
# Waveguide cutoff frequencies (textbook)
# ---------------------------------------------------------------------------

def rectangular_cutoff_wavelength_um(width_um: float, depth_um: float) -> float:
    """TE10 cutoff wavelength of a rectangular guide: 2 * min(W, D)."""
    return RECT_CUTOFF_CONSTANT * min(width_um, depth_um)


def circular_cutoff_wavelength_um(diameter_um: float) -> float:
    """TE11 cutoff wavelength of a circular guide: 1.7062 * D."""
    return CIRC_TE11_CONSTANT * diameter_um


# ---------------------------------------------------------------------------
# Evanescent (sub-cutoff) axial decay
# ---------------------------------------------------------------------------

def evanescent_decay_constant_um(lambda_c_um: float, lambda_um: float) -> float:
    """Axial *field* decay constant kappa (per µm) of a sub-cutoff mode.

        kappa = 2*pi * sqrt(1/lambda_c**2 - 1/lambda**2)

    For lambda > lambda_c the axial propagation constant beta = i*kappa is
    purely imaginary, so the field decays as exp(-kappa*z) and the *power*
    as exp(-2*kappa*z).
    """
    if lambda_um <= lambda_c_um:
        # Propagating (or exactly at cutoff) -> no exponential decay.
        return 0.0
    lc = float(lambda_c_um)
    lam = float(lambda_um)
    return 2.0 * math.pi * math.sqrt(max(0.0, 1.0 / lc ** 2 - 1.0 / lam ** 2))


def evanescent_power_transmission(lambda_c_um: float, lambda_um: float,
                                  h_um: float) -> float:
    """Power that tunnels a length h_um in a sub-cutoff channel: T(h) = exp(-2*kappa*h)."""
    kappa = evanescent_decay_constant_um(lambda_c_um, lambda_um)
    return math.exp(-2.0 * kappa * h_um)


# ---------------------------------------------------------------------------
# Geometric cavity enhancement (used for Monte Carlo diagnostics)
# ---------------------------------------------------------------------------

def geometric_cavity_enhancement(area_walls_m2: float, area_base_m2: float,
                                 area_aperture_m2: float) -> float:
    """C_e = (A_walls + A_base) / A_aperture — the ray-geometric enhancement."""
    if area_aperture_m2 <= 0:
        return 0.0
    return (area_walls_m2 + area_base_m2) / area_aperture_m2
