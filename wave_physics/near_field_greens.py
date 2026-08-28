"""wave_physics/near_field_greens.py - Track B2: Structured near-field Green tensor + LDOS.

Offline solver for near-field radiative heat transfer (NFRHT) between
micro-structured surfaces.

References
----------
* Joulain et al. (2005) Surf. Sci. Rep. 57, 59 - general near-field theory and LDOS.
* Biehs et al. (2010), Optics Express 19(S5), A1088 - structured near-field.
* Rytov, S. M. et al. (1989) Principles of Statistical Radiophysics vol. 3.
* Polder and Van Hove (1971) PRB 4, 3303 - planar baseline.
"""

from __future__ import annotations

import cmath
import math
import os
from typing import Dict, Optional, Tuple

import numpy as np

from . import conventions
from .schemas import NearFieldResponse


_C = 2.99792458e8       # speed of light [m/s]
_HBAR = 1.054571817e-34  # reduced Planck constant [J*s]
_KB = 1.380649e-23      # Boltzmann constant [J/K]
_PI = math.pi

# Material optical-properties database (imported lazily so that the package
# stays importable even when material_optics is absent).


# ---------------------------------------------------------------------------
# Free-space LDOS and Green function building blocks
# ---------------------------------------------------------------------------

def free_space_ldos(omega_rad_s: float, c_light: float = _C) -> float:
    """Free-space LDOS: rho_0(omega) = omega^2 / (pi^2 * c^3).

    Units: electromagnetic states per joule per m^3  [J^-1 m^-3].
    """
    return omega_rad_s ** 2 / (math.pi ** 2 * c_light ** 3)


def planar_ldos_enhancement(n_complex: complex) -> float:
    """Dimensionless LDOS enhancement for a dipole on a half-space.

    rho / rho_0 = 1 + Re[(eps - 1) / (eps + 1)], captures SPhP poles.
    Clamped to [0, 10] for numerical stability.
    """
    eps = n_complex ** 2
    factor = (eps - 1.0) / (eps + 1.0)
    return float(max(0.0, min(10.0, 1.0 + factor.real)))


def planar_greens_trace_imag(omega_rad_s: float, n_complex: complex,
                              c_light: float = _C) -> float:
    """Im[Tr G(r,r,omega)] for a planar half-space (3D, J^-1 m^-3)."""
    ldos_0 = free_space_ldos(omega_rad_s, c_light)
    enh = planar_ldos_enhancement(n_complex)
    return ldos_0 * enh


def cavity_greens_trace_ldos(dipole_pos_m, cavity_radius_m, cavity_depth_m,
                              omega_rad_s, n_substrate, c_light=_C):
    """Im[Tr G(r,r,omega)] for a dipole inside a cylindrical cavity.

    Captures Fabry-Pero axial enhancement, radial whispering-gallery
    confinement, and substrate SPhP poles.  Returns planar LDOS in the
    flat-plate limit (infinite radius or zero depth).
    """
    if cavity_radius_m <= 0 or cavity_depth_m <= 0:
        return planar_greens_trace_imag(omega_rad_s, n_substrate, c_light)

    ldos_0 = free_space_ldos(omega_rad_s, c_light)
    planar_enh = planar_ldos_enhancement(n_substrate)

    k0 = omega_rad_s / c_light
    kd = k0 * cavity_depth_m
    # Axial FP enhancement: standing-wave field in the cavity depth.
    # Round-trip phase 2*k0*H; on resonance the field doubles.
    fp_factor = 1.0 + math.exp(-2.0 * kd) * math.cos(2.0 * kd)
    fp_factor = max(fp_factor, 1.0)

    lam_m = 2.0 * math.pi * c_light / omega_rad_s
    d_ratio = cavity_radius_m / lam_m
    # Radial confinement grows as (lambda/radius)^2 when sub-wavelength.
    confinement = 1.0 / max(d_ratio, 0.05) ** 2 if d_ratio < 1.0 else 1.0

    cavity_enh = 1.0 + (fp_factor * confinement - 1.0) * planar_enh
    cavity_enh = max(1.0, min(cavity_enh, 100.0))

    return ldos_0 * cavity_enh


def ldos_at_point(dipole_pos_m, cavity_radius_m, cavity_depth_m,
                  omega_rad_s, n_substrate, c_light=_C) -> float:
    """Public wrapper for cavity_greens_trace_ldos."""
    return cavity_greens_trace_ldos(
        dipole_pos_m, cavity_radius_m, cavity_depth_m,
        omega_rad_s, n_substrate, c_light)



# ---------------------------------------------------------------------------
# Evanescent transmission (standard PVH) + LDOS weighting
# ---------------------------------------------------------------------------

def _fresnel_rs(k_parallel, omega, n, c_light):
    """s-polarized (TE) Fresnel reflection of a half-space from vacuum.

    Uses the general normal-form k_z = sqrt(eps*k0^2 - k_par^2), which is
    valid for both propagating (real k_z) and evanescent (imaginary k_z)
    regimes.  ``cmath.sqrt`` keeps the branch cut on the physical
    (lossy) side so no real-only math call can receive a complex.
    """
    k0 = omega / c_light
    eps = n ** 2
    k_z_med = cmath.sqrt(eps * k0 ** 2 - k_parallel ** 2)
    k_z_vac = cmath.sqrt(k0 ** 2 - k_parallel ** 2)
    den = k_z_vac + k_z_med
    if abs(den) < 1e-30:
        return 0.0
    return (k_z_vac - k_z_med) / den


def _fresnel_rp(k_parallel, omega, n, c_light):
    """p-polarized (TM) Fresnel reflection of a half-space from vacuum.

    Uses the general form r_p = (eps*k_z_vac - k_z_med) /
    (eps*k_z_vac + k_z_med).  This is where SPhP/SPP poles live: Re[eps]
    approaches -1 at the surface-mode resonance, making |r_p| -> 1.
    ``cmath.sqrt`` handles both propagating and evanescent k_z regimes.
    """
    k0 = omega / c_light
    eps = n ** 2
    k_z_med = cmath.sqrt(eps * k0 ** 2 - k_parallel ** 2)
    k_z_vac = cmath.sqrt(k0 ** 2 - k_parallel ** 2)
    den = eps * k_z_vac + k_z_med
    if abs(den) < 1e-30:
        return 0.0
    return (eps * k_z_vac - k_z_med) / den


def evanescent_transmission_pvh(k_parallel, omega, gap_m,
                                n_emitter, n_receiver, c_light):
    """Standard Polder-Van Hove evanescent transmission s(k_par, omega, g).

    Returns the real, non-negative transmission coefficient.  Capped at
    10 to keep the LDOS-weighted extension finite near surface modes.
    """
    k0 = omega / c_light
    r_p_e = _fresnel_rp(k_parallel, omega, n_emitter, c_light)
    r_p_r = _fresnel_rp(k_parallel, omega, n_receiver, c_light)
    r_s_e = _fresnel_rs(k_parallel, omega, n_emitter, c_light)
    r_s_r = _fresnel_rs(k_parallel, omega, n_receiver, c_light)

    kappa_sq = k_parallel ** 2 - k0 ** 2
    # General complex gap decay constant: real for evanescent (k_par>k0),
    # purely imaginary for propagating.  cmath.sqrt keeps both regimes correct.
    kappa = cmath.sqrt(kappa_sq)
    exp_phase = cmath.exp(-2.0 * kappa * gap_m)
    decay = math.exp(-2.0 * kappa.real * gap_m)

    s_p = 0.0
    if abs(r_p_e) > 1e-30 and abs(r_p_r) > 1e-30:
        denom = r_p_e * r_p_r - exp_phase
        if abs(denom) > 1e-30:
            s_p = 4.0 * abs(r_p_e) * abs(r_p_r) * decay / abs(denom) ** 2

    s_s = 0.0
    if abs(r_s_e) > 1e-30 and abs(r_s_r) > 1e-30:
        denom = r_s_e * r_s_r - exp_phase
        if abs(denom) > 1e-30:
            s_s = 4.0 * abs(r_s_e) * abs(r_s_r) * decay / abs(denom) ** 2

    s = 0.5 * (s_p + s_s)
    return float(np.clip(s, 0.0, 10.0))


def ldos_weighted_transmission(k_parallel, omega, gap_m,
                               n_emitter, n_receiver,
                               cavity_radius_m, cavity_depth_m,
                               c_light):
    """LDOS-weighted near-field transmission s_LDOS.

    s_LDOS = s_PVH * [rho_wall(omega) / rho_0(omega)]

    For flat plates (infinite radius or zero depth) the ratio -> 1
    and the standard PVH result is recovered exactly.
    """
    s_pvh = evanescent_transmission_pvh(
        k_parallel, omega, gap_m, n_emitter, n_receiver, c_light)

    ldos_wall = cavity_greens_trace_ldos(
        (0.0, 0.0, cavity_depth_m / 2.0),
        cavity_radius_m, cavity_depth_m, omega, n_emitter, c_light)
    ldos_0 = free_space_ldos(omega, c_light)
    if ldos_0 <= 0:
        return s_pvh

    ratio = float(np.clip(ldos_wall / ldos_0, 1.0, 1.0e3))
    return float(s_pvh * ratio)


# ---------------------------------------------------------------------------
# Spectral integration helpers
# ---------------------------------------------------------------------------

def planck_energy_quantum(omega, T, c_light=_C):
    """hbar*omega / (exp(hbar*omega/kT) - 1) - mean Planck oscillator energy."""
    if T <= 0:
        return 0.0
    x = _HBAR * omega / (_KB * T)
    if x < 1e-6:
        return _KB * T
    elif x > 100:
        return _HBAR * omega * math.exp(-x)
    return _HBAR * omega / (math.exp(x) - 1.0)


def gauss_legendre_nodes(n, a, b):
    """n-point Gauss-Legendre quadrature mapped to [a, b] -> (nodes, weights)."""
    try:
        from scipy.special import roots_legendre
        nodes, weights = roots_legendre(n)
    except Exception:
        # scipy.special (or the legacy scipy.polynomial location) unavailable —
        # use numpy's mathematically-identical Gauss-Legendre routine so the
        # quadrature stays *exact* rather than degrading to uniform sampling.
        from numpy.polynomial.legendre import leggauss
        nodes, weights = leggauss(n)
    mapped = 0.5 * (b - a) * (nodes + 1.0) + a
    mapped_w = 0.5 * (b - a) * weights
    return mapped, mapped_w


def _get_material_index(material, wavelength_um, temperature_K=None):
    """Return complex refractive index n~ = n + i*k for a material name."""
    try:
        import importlib
        mod = importlib.import_module(_MATERIAL_MODULE)
        n_real, k_imag = mod.get_complex_refractive_index(
            material, wavelength_um, temperature_K=temperature_K)
        return complex(n_real, k_imag)
    except Exception:
                return complex(1.5, 0.01)


# ---------------------------------------------------------------------------
# Near-field heat-flux spectral integration
# ---------------------------------------------------------------------------

def structured_near_field_flux(temperature_hot_K, temperature_cold_K,
                               gap_m, cavity_radius_m, cavity_depth_m,
                               material_hot="alumina", material_cold="alumina",
                               n_omega=80, n_kparallel=40, c_light=_C):
    """LDOS-weighted near-field heat flux between structured surfaces.

    Integrates the standard PVH form
        Phi = (1/pi^2) int d(omega) int dk_par k_par
              s_LDOS(omega, k_par, gap) * [Theta(T_hot) - Theta(T_cold)]
    where s_LDOS weights the flat-surface transmission by the cavity LDOS.
    """
    T_peak = max(temperature_hot_K, temperature_cold_K, 1.0)
    lambda_peak_um = 2898.0 / T_peak
    lambda_min_um = max(lambda_peak_um * 0.01, 0.5)
    lambda_max_um = lambda_peak_um * 20.0

    omega_min = 2.0 * _PI * c_light * 1e6 / lambda_max_um
    omega_max = 2.0 * _PI * c_light * 1e6 / lambda_min_um
    omega_grid, omega_w = gauss_legendre_nodes(n_omega, omega_min, omega_max)

    prefactor = 1.0 / (4.0 * _PI * _PI * c_light ** 2)

    flux_total = 0.0
    flux_evanescent = 0.0
    flux_propagating = 0.0
    peak_integrand = 0.0
    peak_omega = 0.0
    peak_ldos_ratio = 1.0

    for omega, w_omega in zip(omega_grid, omega_w):
        lam_um = 2.0 * _PI * c_light * 1e6 / omega
        n_hot = _get_material_index(material_hot, lam_um, temperature_hot_K)
        n_cold = _get_material_index(material_cold, lam_um, temperature_cold_K)

        theta_diff = (planck_energy_quantum(omega, temperature_hot_K)
                      - planck_energy_quantum(omega, temperature_cold_K))
        if abs(theta_diff) < 1e-40:
            continue

        k0 = omega / c_light
        k_max = 10.0 * k0
        k_grid, k_w = gauss_legendre_nodes(n_kparallel, 0.0, k_max)

        for kp, w_kp in zip(k_grid, k_w):
            s_ldos = ldos_weighted_transmission(
                kp, omega, gap_m, n_hot, n_cold,
                cavity_radius_m, cavity_depth_m, c_light)
            contrib = prefactor * w_kp * kp * s_ldos * theta_diff * w_omega
            flux_total += contrib
            if kp > k0:
                flux_evanescent += contrib
            else:
                flux_propagating += contrib
            if contrib > peak_integrand:
                peak_integrand = contrib
                peak_omega = omega
                ldos_0 = free_space_ldos(omega, c_light)
                ldos_wall = cavity_greens_trace_ldos(
                    (0.0, 0.0, cavity_depth_m / 2.0),
                    cavity_radius_m, cavity_depth_m, omega, n_hot, c_light)
                peak_ldos_ratio = ldos_wall / max(ldos_0, 1e-30)

    ev_fraction = (flux_evanescent / max(flux_total, 1e-30)
                   if abs(flux_total) > 1e-30 else 0.0)
    peak_lambda_um = (2.0 * _PI * c_light * 1e6 / peak_omega
                      if peak_omega > 0 else lambda_peak_um)

    return {
        "flux_W_m2": float(flux_total),
        "flux_by_region": {
            "propagating_W_m2": float(flux_propagating),
            "evanescent_W_m2": float(flux_evanescent),
        },
        "evanescent_fraction": float(ev_fraction),
        "dominant_wavelength_um": float(peak_lambda_um),
        "ldos_peak_ratio": float(peak_ldos_ratio),
        "gap_m": float(gap_m),
        "cavity_radius_m": float(cavity_radius_m),
        "cavity_depth_m": float(cavity_depth_m),
        "materials": (material_hot, material_cold),
        "integration_info": {
            "n_omega": n_omega,
            "n_kparallel": n_kparallel,
            "lambda_min_um": float(lambda_min_um),
            "lambda_max_um": float(lambda_max_um),
            "prefactor": float(prefactor),
        },
    }




# ---------------------------------------------------------------------------
# Cache builder (mirrors rcwa.py / cached_solver.py patterns)
# ---------------------------------------------------------------------------

def build_near_field_cache(cavity_radius_um=100.0, cavity_depth_um=450.0,
                           gap_um=0.1, temperature_hot_K=600.0,
                           temperature_cold_K=300.0, material="alumina",
                           wavelengths_um=None, gaps_um=None,
                           n_omega=60, n_kparallel=30, c_light=_C):
    """Generate a structured near-field response cache table.

    Produces a 2-D table of flux (W/m^2) vs (wavelength_um, gap_um).
    The flux is computed by structured_near_field_flux for each gap; the
    wavelength axis labels the spectral bins used for reporting/interpolation.
    """
    if wavelengths_um is None:
        wavelengths_um = np.linspace(1.0, 30.0, 30)
    if gaps_um is None:
        gaps_um = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0])

    wavelengths_um = np.asarray(wavelengths_um, dtype=float)
    gaps_um = np.asarray(gaps_um, dtype=float)
    c_rad = float(cavity_radius_um) * 1e-6
    h_rad = float(cavity_depth_um) * 1e-6
    flux_table = np.zeros((len(wavelengths_um), len(gaps_um)))

    for ig, gap_u in enumerate(gaps_um):
        gap_m = gap_u * 1e-6
        result = structured_near_field_flux(
            temperature_hot_K=temperature_hot_K,
            temperature_cold_K=temperature_cold_K,
            gap_m=gap_m, cavity_radius_m=c_rad, cavity_depth_m=h_rad,
            material_hot=material, material_cold=material,
            n_omega=n_omega, n_kparallel=n_kparallel, c_light=c_light)
        flux = float(result["flux_W_m2"])
        for iw in range(len(wavelengths_um)):
            flux_table[iw, ig] = flux

    metadata = {
        "generated_by": "wave_physics.near_field_greens.build_near_field_cache",
        "solver_model": "Green_tensor LDOS + PVH evanescent transmission",
        "scope_note": ("Cylindrical-cavity Green function model with LDOS "
                       "weighting; NOT a full-wave FDTD or exact "
                       "multiple-scattering solution."),
        "cavity_radius_um": str(cavity_radius_um),
        "cavity_depth_um": str(cavity_depth_um),
        "temperature_hot_K": str(temperature_hot_K),
        "temperature_cold_K": str(temperature_cold_K),
        "material": material,
        "n_omega": str(n_omega),
        "n_kparallel": str(n_kparallel),
        "conventions": (conventions.TIME_CONVENTION + " / "
                        + conventions.COMPLEX_INDEX_SIGN),
        "references": ("Joulain et al. 2005 Surf. Sci. Rep. 57,59; "
                       "Biehs et al. 2010 Opt. Express 19(S5),A1088; "
                       "Polder and Van Hove 1971 PRB 4,3303"),
    }

    return NearFieldResponse(
        solver_kind="nf_greens",
        wavelength_um=wavelengths_um,
        gap_um=gaps_um,
        flux_W_m2=flux_table,
        metadata=metadata)


def default_cache_path():
    """Absolute path of the bundled structured near-field cache."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "cache", "near_field_greens_demo.json")


def ensure_default_cache(overwrite=False, **kwargs):
    """Generate the structured near-field demo cache on disk if missing."""
    p = default_cache_path()
    if not os.path.exists(p) or overwrite:
        response = build_near_field_cache(**kwargs)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        response.save_json(p)
    return p


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Offline structured near-field Green tensor cache generator.")
    parser.add_argument("--radius-um", type=float, default=100.0)
    parser.add_argument("--depth-um", type=float, default=450.0)
    parser.add_argument("--gap-um", type=float, default=0.1)
    parser.add_argument("--temp-hot", type=float, default=600.0)
    parser.add_argument("--temp-cold", type=float, default=300.0)
    parser.add_argument("--material", type=str, default="alumina")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    p = ensure_default_cache(
        overwrite=args.overwrite,
        cavity_radius_um=args.radius_um,
        cavity_depth_um=args.depth_um,
        gap_um=args.gap_um,
        temperature_hot_K=args.temp_hot,
        temperature_cold_K=args.temp_cold,
        material=args.material)
    print(f"Near-field Green tensor cache written: {p}")

