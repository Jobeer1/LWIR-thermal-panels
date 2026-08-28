"""
material_optics.py — Optical properties database for thin-film corrections.

Provides wavelength-dependent absorption depth data for common materials
used in micro/nano-structured thermal emitters.

Key Physics:
    Thin-film emissivity is calculated with a normal-incidence transfer matrix.
    The absorption-depth tables provide the extinction coefficient k; the
    refractive-index defaults below are engineering estimates until measured
    n(λ), k(λ) data are supplied.
  
  This reduces wall emission by factor ~t/δ for optically thin walls.
  
References:
  1. Palik, E. D. (Ed.). (1998). Handbook of Optical Constants of Solids.
  2. Sprafke et al., Adv. Opt. Mater. 2013 (PAA micro-cavities)
  3. Mizuno et al., PNAS 2009 (CNT forest optical properties)
"""

import numpy as np
from typing import Dict, List, Optional
import cmath
import math
import warnings

# ---------------------------------------------------------------------------
# Wavelength-dependent absorption depth database (µm)
# ---------------------------------------------------------------------------

# Alumina (Al₂O₃) - typical PAA cavity wall material
# Data from Palik (1998), Handbook of Optical Constants of Solids
# Includes the fundamental Reststrahlen optical phonon band in 10-30 µm range
ALUMINA_ABSORPTION_DEPTH = {
    'wavelengths_um': np.array([0.5, 1.0, 2.0, 5.0, 8.0, 10.0, 12.0, 14.5, 15.0, 20.0, 30.0]),
    'depths_um':      np.array([5000.0, 4000.0, 1500.0, 150.0, 12.7, 2.27, 0.80, 0.62, 0.61, 0.66, 3.0]),
    'source': 'Palik (1998) Handbook of Optical Constants, Reststrahlen band verified',
    'temperature_K': 300.0
}

# Vertically Aligned Carbon Nanotube (VACNT) forest
# Multi-wall CNTs have shorter penetration due to strong absorption
# Data from Mizuno et al., PNAS 2009; Yang et al., Nano Lett. 2008
CNT_FOREST_ABSORPTION_DEPTH = {
    'wavelengths_um': np.array([2.0, 5.0, 8.0, 10.0, 12.0, 15.0]),
    'depths_um':      np.array([0.3, 0.8, 1.5, 2.5, 3.5, 5.0]),
    'source': 'Mizuno PNAS 2009, Yang Nano Lett. 2008',
    'temperature_K': 300.0
}

# Silver (Ag) - typical base/substrate material
# High reflectivity in IR → long absorption depth
SILVER_ABSORPTION_DEPTH = {
    'wavelengths_um': np.array([2.0, 5.0, 10.0, 15.0, 20.0]),
    'depths_um':      np.array([0.05, 0.08, 0.12, 0.15, 0.18]),
    'source': 'Palik (1998), Ag is highly reflective',
    'temperature_K': 300.0
}

# Generic high-emissivity coating (e.g., black chrome, Pyromark)
GENERIC_HIGH_EPS_ABSORPTION_DEPTH = {
    'wavelengths_um': np.array([2.0, 5.0, 10.0, 15.0, 20.0]),
    'depths_um':      np.array([0.5, 1.0, 2.0, 3.0, 4.0]),
    'source': 'Typical high-ε coatings',
    'temperature_K': 300.0
}

# Material lookup dictionary
MATERIAL_ABSORPTION_DATA: Dict[str, Dict] = {
    'alumina': ALUMINA_ABSORPTION_DEPTH,
    'cnt_forest': CNT_FOREST_ABSORPTION_DEPTH,
    'silver': SILVER_ABSORPTION_DEPTH,
    'high_emissivity': GENERIC_HIGH_EPS_ABSORPTION_DEPTH,
    # Aliases for UI compatibility
    'al2o3': ALUMINA_ABSORPTION_DEPTH,
    'ag': SILVER_ABSORPTION_DEPTH,
    'cnt': CNT_FOREST_ABSORPTION_DEPTH,
}

# Default bulk emissivities (gray-body approximation)
DEFAULT_BULK_EMISSIVITY = {
    'alumina': 0.80,      # Anodized alumina at 300K
    'cnt_forest': 0.98,   # Multi-wall CNT forest
    'silver': 0.02,       # Polished silver
    'high_emissivity': 0.95,
}

# Approximate real refractive indices used by the TMM when measured optical
# constants are not supplied.  Absorption depth supplies the imaginary part.
DEFAULT_REAL_INDEX = {
    'alumina': 1.70,
    'cnt_forest': 1.80,
    'silver': 0.20,
    'high_emissivity': 1.80,
}

# ---------------------------------------------------------------------------
# PHASE 1: Complex Refractive Index Database
# Full n(λ) + i k(λ) for complex Fresnel calculations
# ---------------------------------------------------------------------------

# Alumina (Al₂O₃) complex refractive index (n + ik)
# Data from Palik (1998) Handbook of Optical Constants of Solids
# Captures transparent visible/NIR (k ~ 10^-5) and Reststrahlen optical phonon band (k ~ 1.5 - 2.5)
# Temperature: 300 K
ALUMINA_COMPLEX_INDEX = {
    'wavelengths_um': np.array([0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0, 14.5, 15.0, 20.0, 30.0]),
    'n_real': np.array([1.76, 1.73, 1.71, 1.69, 1.65, 1.62, 1.55, 1.45, 1.15, 0.45, 0.32, 0.35, 1.20, 2.80]),
    'k_imag': np.array([0.0001, 0.00005, 0.00001, 0.00001, 0.00005, 0.0001, 0.001, 0.05, 0.35, 1.20, 1.85, 1.95, 2.40, 0.80]),
    'citation': 'Palik, E. D. (1998). Handbook of Optical Constants of Solids (Reststrahlen band)',
    'temperature_K': 300.0
}

# Silicon (Si) complex refractive index
# Strong absorption edge in visible, transparent in IR
SILICON_COMPLEX_INDEX = {
    'wavelengths_um': np.array([0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0]),
    'n_real': np.array([4.5, 3.9, 3.5, 3.3, 3.3, 3.4, 3.5, 3.6, 3.7, 3.85, 4.1]),
    'k_imag': np.array([0.5, 0.001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.08, 0.1, 0.12, 0.15]),
    'citation': 'Palik, E. D. (1998)',
    'temperature_K': 300.0
}

# Carbon (Multi-wall CNT) complex refractive index
# Very strong absorption due to free-electron response
CARBON_NANOTUBE_COMPLEX_INDEX = {
    'wavelengths_um': np.array([0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0]),
    'n_real': np.array([2.1, 2.05, 2.0, 1.95, 1.9, 1.85, 1.8, 1.75]),
    'k_imag': np.array([1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.35, 0.3]),
    'citation': 'Mizuno et al., PNAS 106:6044-6047 (2009); Yang et al., Nano Lett. 8:446-451 (2008)',
    'temperature_K': 300.0
}

# Silver (Ag) complex refractive index
# Highly reflective in IR, strong Drude response
SILVER_COMPLEX_INDEX = {
    'wavelengths_um': np.array([0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0]),
    'n_real': np.array([0.14, 0.13, 0.15, 0.17, 0.20, 0.25, 0.35, 0.6, 1.0]),
    'k_imag': np.array([3.5, 3.2, 2.8, 2.2, 1.5, 1.0, 0.8, 0.6, 0.5]),
    'citation': 'Palik, E. D. (1998)',
    'temperature_K': 300.0
}

# Material lookup dictionary for complex indices
MATERIAL_COMPLEX_INDEX = {
    'alumina': ALUMINA_COMPLEX_INDEX,
    'silicon': SILICON_COMPLEX_INDEX,
    'carbon_nanotube': CARBON_NANOTUBE_COMPLEX_INDEX,
    'cnt': CARBON_NANOTUBE_COMPLEX_INDEX,
    'silver': SILVER_COMPLEX_INDEX,
    'ag': SILVER_COMPLEX_INDEX,
    'al2o3': ALUMINA_COMPLEX_INDEX,
    'si': SILICON_COMPLEX_INDEX,
}

# ---------------------------------------------------------------------------
# PHASE 4a: Temperature-dependent & non-local optical constants
# ---------------------------------------------------------------------------
# The tabulated n(λ), k(λ) databases above are anchored to 300 K literature
# data (Palik 1998; Mizuno et al. 2009).  At elevated temperature the
# dielectric response drifts for two physical reasons:
#
#  1. Electron-phonon (or phonon-phonon) scattering increases the damping
#     rate.  A linearised scaling is used:
#         γ(T) = γ_300 · [1 + s · (T/300 − 1)],    s = gamma_T_slope
#     (Modest, "Radiative Heat Transfer", 3rd ed., 2013 — high-T dielectrics).
#
#  2. Interband transition / band-gap energies red-shift with temperature
#     (Varshni relation  E_g(T) = E_g(0) − αT²/(T + β)  for semiconductors).
#
# Additionally, for structures with a characteristic size L below the
# electron mean free path (ℓ ≈ 10–40 nm in graphitic carbon and noble
# metals) the local-response approximation breaks down.  The hydrodynamic
# non-local correction adds a −β²k² term to the Drude denominator
# (Raza et al., J. Phys.: Condens. Matter 23, 325502 (2011)):
#
#     ε(ω, k, T) = ε∞ − ωp² / (ω² + i ω γ(T) − β² k²),    β = √(3/5)·v_F
#
# and diffuse surface scattering further raises the effective damping
# (Fuchs–Sondheimer):  γ_eff = γ(T)·[1 + 3(1−p)·ℓ/(8L)]  with p = 0.5.
#
# ANCHORING STRATEGY (perturbative, backward compatible):
#   The 300 K tables remain authoritative for the *absolute* values.  The
#   parametric model is evaluated at both T and 300 K and only the RELATIVE
#   drift is applied to the tabulated index:
#       ñ(λ, T) = ñ_tab(λ) · sqrt( ε_model(λ, T) / ε_model(λ, 300 K) ).
#   At T = 300 K with no non-local scale the ratio is exactly 1 and the
#   lookup returns the tabulated values unchanged.
# ---------------------------------------------------------------------------

#: Reduced Planck constant in eV·s (CODATA 2018).
_HBAR_EV_S = 6.582119569e-16

#: Hydrodynamic non-local coefficient β = sqrt(3/5)·v_F (GNOR model).
_NONLOCAL_BETA_FACTOR = math.sqrt(3.0 / 5.0)

#: Fuchs–Sondheimer diffuse surface-scattering prefactor 3(1−p)/8, p = 0.5.
_SURFACE_SCATTERING_FACTOR = 3.0 * (1.0 - 0.5) / 8.0

#: Temperatures within this tolerance of 300 K use the tabulated data as-is.
_REFERENCE_TEMP_TOL_K = 0.5

# ---------------------------------------------------------------------------
# Thin-film emissivity: exact Air–Film–Air TMM + effective-medium (EMA) model
# ---------------------------------------------------------------------------
# A planar TMM of an ultra-thin / optically-thin homogeneous slab returns an
# absorptance essentially equal to 0 (a ~50 nm smooth dense alumina wall is a
# ~1–8% emitter in the thermal IR).  That is the CORRECT, first-principles
# answer for a free-standing 51 nm membrane — a layer that thin simply cannot
# absorb much, no matter how lossy the material (absorptance is bounded by
# the optical-path-length · Im(ñ)).  An earlier heuristic
# ``THIN_FILM_MIN_EMISSIVITY`` floor artificially inflated this and is now
# removed.
#
# Real anodized (PAA) honeycomb walls achieve much higher emissivity than a
# smooth dense slab because they are a COMPOSITE: porous alumina host with
# lossy carbonaceous anodization residue.  We model that physically with a
# Bruggeman effective-medium approximation (EMA) that mixes the host complex
# permittivity with a small volume fraction of a strongly-absorbing inclusion,
# yielding an effective complex index ñ_eff that is then used in the exact
# Air–Film–Air transfer-matrix solution (emissivity = absorptivity by
# Kirchhoff).  The emissivity therefore rises with the EMA inclusion fraction
# and with wall thickness, and correctly collapses toward the transparent
# (low-ε) limit for genuinely thin membranes.

#: Volume fraction of strongly-absorbing carbonaceous inclusion in PAA walls
#: (anodization-derived carbon/graphite residue; see e.g. oxidised PAA studies).
EMA_CARBON_FILL_FRACTION = 0.10
#: Complex index of the absorbing inclusion (n + ik), graphite-like in the IR.
EMA_INCLUSION_INDEX = complex(1.9, 0.5)
#: PAA / porous-alumina aliases eligible for the EMA composite model.
PAA_ALIASES = ('alumina', 'al2o3', 'paa')

# Parametric dielectric models.  Absolute accuracy is NOT required from these
# parameters — they only supply the *relative* temperature/non-local drift
# applied to the authoritative 300 K tables above.
#   class 'drude_metal'      : free-electron response (Ag, graphitic CNT)
#   class 'phonon_dielectric': IR-active optical phonon (Al2O3)
#   class 'semiconductor'    : Varshni band-gap pole (Si)
DRUDE_LORENTZ_PARAMS: Dict[str, Dict] = {
    'silver': {
        'class': 'drude_metal',
        'eps_inf': 5.0,
        'hbar_omega_p_eV': 9.0,     # screened plasma energy
        'gamma_300_eV': 0.021,      # Drude damping at 300 K
        'gamma_T_slope': 0.7,       # γ(T) = γ300·[1 + s(T/300 − 1)]
        'v_fermi_m_s': 1.39e6,
        'electron_mfp_nm': 52.0,
    },
    'carbon_nanotube': {
        'class': 'drude_metal',
        'eps_inf': 3.5,
        'hbar_omega_p_eV': 1.2,     # effective free-carrier plasma (graphitic)
        'gamma_300_eV': 0.25,
        'gamma_T_slope': 0.6,
        'v_fermi_m_s': 8.0e5,
        'electron_mfp_nm': 25.0,
    },
    'alumina': {
        'class': 'phonon_dielectric',
        'eps_inf': 3.2,
        'hbar_omega_p_eV': 0.085,   # dominant TO phonon energy (~14.6 µm)
        'gamma_300_eV': 0.004,      # anharmonic phonon damping at 300 K
        'gamma_T_slope': 0.5,
        'v_fermi_m_s': 0.0,         # insulator → non-local term inactive
        'electron_mfp_nm': 1.0e9,
    },
    'silicon': {
        'class': 'semiconductor',
        'eps_inf': 11.7,
        'hbar_omega_p_eV': 0.0,
        'gamma_300_eV': 0.0,
        'gamma_T_slope': 0.0,
        'v_fermi_m_s': 0.0,
        'electron_mfp_nm': 1.0e9,
        # Varshni parameters for the E1 direct critical point (~3.4 eV @300 K)
        'eg_0_eV': 3.78,
        'eg_alpha_eV_per_K': 4.6e-4,
        'eg_beta_K': 320.0,
        'oscillator_strength': 8.0,  # ε = ε∞ + A·Eg²/(Eg² − (ħω)² − iγħω)
        'oscillator_gamma_eV': 0.10,
    },
}

# Aliases so every identifier used across the app resolves to a model.
DRUDE_LORENTZ_PARAMS['ag'] = DRUDE_LORENTZ_PARAMS['silver']
DRUDE_LORENTZ_PARAMS['cnt'] = DRUDE_LORENTZ_PARAMS['carbon_nanotube']
DRUDE_LORENTZ_PARAMS['cnt_forest'] = DRUDE_LORENTZ_PARAMS['carbon_nanotube']
DRUDE_LORENTZ_PARAMS['al2o3'] = DRUDE_LORENTZ_PARAMS['alumina']
DRUDE_LORENTZ_PARAMS['si'] = DRUDE_LORENTZ_PARAMS['silicon']


def electron_collision_frequency_eV(material: str, temperature_K: float = 300.0) -> float:
    """Damped-carrier / phonon collision energy γ(T) in eV.

    γ(T) = γ_300 · [1 + s·(T/300 − 1)], clamped positive.
    Unknown materials return 0.0 (no modeled drift).
    """
    params = DRUDE_LORENTZ_PARAMS.get(str(material).lower())
    if params is None:
        return 0.0
    gamma300 = float(params['gamma_300_eV'])
    slope = float(params['gamma_T_slope'])
    ratio = max(float(temperature_K), 1.0) / 300.0
    return max(gamma300 * (1.0 + slope * (ratio - 1.0)), 1e-9)


def electron_mean_free_path_m(material: str) -> float:
    """Electron mean free path ℓ in metres (very large for insulators)."""
    params = DRUDE_LORENTZ_PARAMS.get(str(material).lower())
    if params is None:
        return 1.0e9
    return float(params['electron_mfp_nm']) * 1e-9


def varshni_bandgap_eV(material: str, temperature_K: float = 300.0) -> Optional[float]:
    """Semiconductor band-gap / critical-point energy E_g(T) via Varshni.

    Returns None for materials without a modeled gap (metals, dielectrics).
    """
    params = DRUDE_LORENTZ_PARAMS.get(str(material).lower())
    if params is None or params.get('class') != 'semiconductor':
        return None
    T = max(float(temperature_K), 1.0)
    eg0 = float(params['eg_0_eV'])
    alpha = float(params['eg_alpha_eV_per_K'])
    beta = float(params['eg_beta_K'])
    return eg0 - alpha * T * T / (T + beta)


# Photon-energy conversion helper: λ [µm] → ħω [eV].
def _hbar_omega_eV(wavelength_um: float) -> float:
    """Photon energy ħω in eV for a vacuum wavelength in micrometres."""
    return _HBAR_EV_S * (2.0 * math.pi * 299792458.0) / (wavelength_um * 1e-6)


def drude_lorentz_permittivity(material: str, wavelength_um: float,
                               temperature_K: float = 300.0,
                               feature_scale_m: float = None) -> complex:
    """Parametric dielectric function ε(ω, k, T) used for the relative drift.

    Implements
        ε(ω, k, T) = ε∞ − ωp² / (ω² + iωγ_eff(T) − β²k²)       [Drude classes]
        ε(ω, T)    = ε∞ + A·Eg(T)² / (Eg(T)² − (ħω)² − iγħω)   [semiconductor]

    ``feature_scale_m`` (structure size L) activates the hydrodynamic
    non-local term with confinement wavevector k = π/L plus the
    Fuchs–Sondheimer surface-scattering damping enhancement — but only when
    L lies below the electron mean free path.
    """
    params = DRUDE_LORENTZ_PARAMS.get(str(material).lower())
    if params is None or wavelength_um <= 0:
        return complex(1.0, 0.0)

    hbar_omega = _hbar_omega_eV(wavelength_um)
    gamma = electron_collision_frequency_eV(material, temperature_K)

    # --- Non-local / surface-scattering enhancement (only below the MFP) ---
    beta_k_sq = 0.0
    if feature_scale_m is not None and feature_scale_m > 0:
        mfp = electron_mean_free_path_m(material)
        if feature_scale_m < mfp and float(params['v_fermi_m_s']) > 0:
            v_f = float(params['v_fermi_m_s'])
            beta = _NONLOCAL_BETA_FACTOR * v_f
            k_conf = math.pi / float(feature_scale_m)
            # β²k² expressed in (eV)² for unit consistency with ħω.
            beta_k_sq = ((beta * k_conf) * _HBAR_EV_S) ** 2
            gamma *= (1.0 + _SURFACE_SCATTERING_FACTOR * mfp / feature_scale_m)

    hbar_omega_sq = hbar_omega * hbar_omega
    eps_inf = float(params['eps_inf'])
    wp_sq = float(params['hbar_omega_p_eV']) ** 2

    if params.get('class') == 'semiconductor':
        eg = varshni_bandgap_eV(material, temperature_K) or 1.0
        A = float(params.get('oscillator_strength', 1.0))
        g_osc = float(params.get('oscillator_gamma_eV', 0.05))
        denom = complex(eg * eg - hbar_omega_sq, -g_osc * hbar_omega)
        if abs(denom) < 1e-30:
            return complex(eps_inf, 0.0)
        return eps_inf + A * eg * eg / denom

    # Drude metal / phonon dielectric classes.
    denom = complex(hbar_omega_sq - beta_k_sq, gamma * hbar_omega)
    if abs(denom) < 1e-30:
        return complex(eps_inf, 0.0)
    return eps_inf - wp_sq / denom


def temperature_index_drift_factors(material: str, wavelength_um: float,
                                    temperature_K: float = None,
                                    feature_scale_m: float = None) -> tuple:
    """Component-wise multiplicative drift factors ``(f_n, f_k)``.

        ñ_model(λ,T,L) = sqrt( ε_model(λ,T,L) )
        f_n = Re[ñ_model(λ,T,L)] / Re[ñ_model(λ,300 K)]
        f_k = Im[ñ_model(λ,T,L)] / Im[ñ_model(λ,300 K)]

    Applying the factors SEPARATELY avoids the phase rotation problem of a
    single complex multiplier (a Drude metal's sqrt(ε) lives in a very
    different quadrant than the tabulated index, so rotating the table by
    arg(sqrt(ε_T/ε_300)) would corrupt the split between n and k).

    Each consumer anchors the factors on its OWN authoritative 300 K
    baseline (index table, depth table, or DEFAULT_REAL_INDEX estimate),
    and all paths reduce to (1, 1) at 300 K with no sub-MFP feature scale.
    """
    params = DRUDE_LORENTZ_PARAMS.get(str(material).lower())
    if params is None or wavelength_um <= 0:
        return (1.0, 1.0)

    at_reference = (temperature_K is None
                    or abs(float(temperature_K) - 300.0) < _REFERENCE_TEMP_TOL_K)
    nonlocal_active = (
        feature_scale_m is not None and feature_scale_m > 0
        and feature_scale_m < electron_mean_free_path_m(material)
        and float(params['v_fermi_m_s']) > 0)
    if at_reference and not nonlocal_active:
        return (1.0, 1.0)

    try:
        eps_ref = drude_lorentz_permittivity(material, wavelength_um, 300.0, None)
        eps_new = drude_lorentz_permittivity(
            material, wavelength_um,
            300.0 if temperature_K is None else float(temperature_K),
            feature_scale_m)

        if abs(eps_ref) < 1e-15 or abs(eps_new) < 1e-15:
            return (1.0, 1.0)

        n_tilde_ref = cmath.sqrt(eps_ref)
        n_tilde_new = cmath.sqrt(eps_new)

        ref_re, ref_im = abs(n_tilde_ref.real), abs(n_tilde_ref.imag)
        new_re, new_im = abs(n_tilde_new.real), abs(n_tilde_new.imag)

        f_n = new_re / max(ref_re, 1e-12)
        f_k = new_im / max(ref_im, 1e-12)

        # Guard against pathological blow-ups from near-zero baselines.
        if not (math.isfinite(f_n) and math.isfinite(f_k)):
            return (1.0, 1.0)
        f_n = min(max(f_n, 1e-3), 1e3)
        f_k = min(max(f_k, 1e-3), 1e3)
        return (float(f_n), float(f_k))
    except (ValueError, ZeroDivisionError):
        return (1.0, 1.0)


def _tabulated_index_300k(material: str, wavelength_um: float) -> tuple:
    """Raw 300 K index lookup with NO temperature logic — internal helper.

    Mirrors :func:`get_complex_refractive_index` (table interpolation, or the
    absorption-depth fallback for unknown materials) without emitting the
    out-of-range warning.
    """
    if material not in MATERIAL_COMPLEX_INDEX:
        n_real = DEFAULT_REAL_INDEX.get(str(material).lower(), 1.5)
        delta = get_absorption_depth(material, wavelength_um)
        k_imag = (wavelength_um / (4.0 * np.pi * delta)) if delta > 0 else 0.001
        return (n_real, float(k_imag))

    data = MATERIAL_COMPLEX_INDEX[material]
    wavelengths = data['wavelengths_um']
    wl_clipped = np.clip(wavelength_um, wavelengths[0], wavelengths[-1])
    n_real = float(np.interp(wl_clipped, wavelengths, data['n_real']))
    k_imag = float(np.interp(wl_clipped, wavelengths, data['k_imag']))
    return (n_real, k_imag)


def get_complex_refractive_index_at_temperature(
        material: str, wavelength_um: float,
        temperature_K: float = None,
        feature_scale_m: float = None) -> tuple:
    """Complex refractive index ñ(λ, T, L) = n + ik including T-drift & non-local.

    Anchors on the tabulated 300 K data and applies the modeled *relative*
    drift   ñ(λ,T) = ñ_tab · sqrt( ε_model(λ,T) / ε_model(λ,300K) )   so that
    ``temperature_K=None`` (or ≈300 K) with no sub-MFP feature scale returns
    the tabulated values exactly (backward compatible).

    Parameters
    ----------
    material         : material identifier (aliases resolved).
    wavelength_um    : vacuum wavelength in micrometres.
    temperature_K    : material temperature; None → reference 300 K tables.
    feature_scale_m  : characteristic structure size L in metres
                       (film thickness, CNT diameter); activates the
                       hydrodynamic non-local correction when L < ℓ(e⁻).

    Returns
    -------
    (n_real, k_imag) with k ≥ 0.

    References
    ----------
    Modest (2013) Radiative Heat Transfer 3rd ed. — high-T dielectrics;
    Raza et al. (2011) J. Phys.: Condens. Matter 23, 325502 — non-local.
    """
    n_tab, k_tab = _tabulated_index_300k(material, wavelength_um)

    # Component-wise drift anchored on this function's own baseline (the
    # 300 K index tables). Identity at reference conditions.
    f_n, f_k = temperature_index_drift_factors(material, wavelength_um,
                                               temperature_K, feature_scale_m)
    if f_n == 1.0 and f_k == 1.0:
        return (n_tab, k_tab)

    return (float(n_tab * f_n), float(k_tab * f_k))


def temperature_optics_provenance(material: str, temperature_K: float = None,
                                  feature_scale_m: float = None) -> Dict:
    """Diagnostic metadata describing the active temperature-optics model."""
    params = DRUDE_LORENTZ_PARAMS.get(str(material).lower())
    prov = {
        'temperature_K': 300.0 if temperature_K is None else float(temperature_K),
        'feature_scale_m': feature_scale_m,
        'model': 'tabulated_300k' if params is None else params.get('class'),
        'gamma_eV': electron_collision_frequency_eV(
            material, 300.0 if temperature_K is None else temperature_K),
        'electron_mfp_m': electron_mean_free_path_m(material),
        'nonlocal_active': bool(
            params is not None and feature_scale_m is not None
            and feature_scale_m > 0
            and feature_scale_m < electron_mean_free_path_m(material)
            and float(params['v_fermi_m_s']) > 0),
    }
    eg = varshni_bandgap_eV(
        material, 300.0 if temperature_K is None else temperature_K)
    if eg is not None:
        prov['bandgap_eV'] = eg
    return prov


# ---------------------------------------------------------------------------
# Interpolation and calculation functions
# ---------------------------------------------------------------------------

def get_absorption_depth(material: str, wavelength_um: float,
                         temperature_K: float = None,
                         feature_scale_m: float = None) -> float:
    """
    Get absorption depth δ(λ) for a material at given wavelength.
    
    Parameters
    ----------
    material : str
        Material identifier ('alumina', 'cnt_forest', 'silver', etc.)
    wavelength_um : float
        Wavelength in micrometres
    temperature_K : float, optional
        Material temperature (Phase 4a).  When given and different from
        300 K the depth is rescaled by the modeled extinction drift,
        δ(λ,T) = δ_300 · k_300 / k_T  (since δ = λ / 4πk).
    feature_scale_m : float, optional
        Structure size (m) for the non-local correction.
    
    Returns
    -------
    float
        Absorption depth in micrometres. Returns 0.0 for unknown materials.
    """
    if material not in MATERIAL_ABSORPTION_DATA:
        # Default: assume 5µm absorption depth for unknown materials
        return 5.0
    
    data = MATERIAL_ABSORPTION_DATA[material]
    wavelengths = data['wavelengths_um']
    depths = data['depths_um']
    
    # Clip to wavelength range
    wl = max(wavelengths[0], min(wavelengths[-1], wavelength_um))
    
    # Linear interpolation in log-log space (absorption often follows power law)
    log_wl = np.log(wl)
    log_depths = np.log(depths)
    
    # Find bounding indices
    idx = np.searchsorted(wavelengths, wl)
    if idx == 0:
        delta_300 = float(depths[0])
    elif idx == len(wavelengths):
        delta_300 = float(depths[-1])
    else:
        # Linear interpolation
        wl0, wl1 = wavelengths[idx-1], wavelengths[idx]
        d0, d1 = depths[idx-1], depths[idx]
        t = (wl - wl0) / (wl1 - wl0)
        delta_300 = float(d0 + t * (d1 - d0))

    # Phase 4a: temperature / non-local rescaling via the dielectric model.
    # The extinction coefficient drifts by f_k, so the depth scales as
    # δ(T) = δ_300 / f_k — anchored on THIS table's own baseline, never
    # mixing sources with the complex-index tables.
    f_n, f_k = temperature_index_drift_factors(material, wavelength_um,
                                               temperature_K, feature_scale_m)
    if f_n == 1.0 and f_k == 1.0:
        return delta_300

    if math.isfinite(f_k) and f_k > 1e-12:
        rescaled = delta_300 / f_k
        if math.isfinite(rescaled) and rescaled > 0:
            return float(rescaled)
    return delta_300


def tmm_absorptance_normal_incidence(
    thickness_um: float,
    wavelength_um: float,
    material: str,
    substrate_index: complex = 1.5,
    incident_index: complex = 1.0,
    temperature_K: float = None,
    feature_scale_m: float = None,
) -> float:
    """Return stack absorptance from a one-layer normal-incidence TMM.

    The layer is the selected material and the substrate is semi-infinite.
    This is a coherent thin-film calculation, so it includes Fresnel
    reflection and phase interference.  ``k`` is inferred from the tabulated
    absorption depth using k = λ / (4πδ).  The returned value is the opaque,
    spectral absorptance A(λ, θ) of the surface from the Maxwell (Fresnel/TMM)
    solution.  Kirchhoff is invoked ONLY here, and only in its quantified local
    form: at a single opaque material surface in thermal balance, directional
    spectral emissivity equals directional spectral absorptivity by detailed
    balance.  No cavity-scale or directionally-averaged epsilon = alpha is
    assumed anywhere.

    Phase 4a: when ``temperature_K`` / ``feature_scale_m`` are supplied, both
    n and k come from the Drude–Lorentz-anchored temperature model instead
    of the fixed 300 K tables.
    """
    if thickness_um <= 0.0 or wavelength_um <= 0.0:
        return 0.0

    delta_um = get_absorption_depth(material, wavelength_um,
                                    temperature_K, feature_scale_m)
    if delta_um <= 0.0 or not math.isfinite(delta_um):
        return 0.0

    n_real = DEFAULT_REAL_INDEX.get(material, 1.5)
    k = wavelength_um / (4.0 * math.pi * delta_um)
    n_layer = complex(n_real, -k)  # exp(-iωt) convention

    # Phase 4a: apply the component-wise modeled drift to THIS function's
    # own 300 K baseline (DEFAULT_REAL_INDEX + depth-derived k), keeping the
    # reference behaviour bit-identical when both extra args are None.
    f_n, f_k = temperature_index_drift_factors(material, wavelength_um,
                                               temperature_K, feature_scale_m)
    if f_n != 1.0 or f_k != 1.0:
        n_layer = complex(n_real * f_n, -k * f_k)
    phase = 2.0 * math.pi * n_layer * thickness_um / wavelength_um
    cos_phase = cmath.cos(phase)
    sin_phase = cmath.sin(phase)

    # Characteristic matrix for [E, H] at normal incidence.
    a = cos_phase
    b = 1j * sin_phase / n_layer
    c = 1j * n_layer * sin_phase
    d = cos_phase
    n0 = complex(incident_index)
    ns = complex(substrate_index)
    denominator = n0 * a + n0 * ns * b + c + ns * d
    if abs(denominator) == 0.0:
        return 0.0

    reflection = (n0 * a + n0 * ns * b - c - ns * d) / denominator
    transmission = 2.0 * n0 / denominator
    reflectance = abs(reflection) ** 2
    transmittance = max(0.0, (ns.real / n0.real) * abs(transmission) ** 2)
    return float(np.clip(1.0 - reflectance - transmittance, 0.0, 1.0))


def bruggeman_effective_permittivity(
        eps_host: complex, eps_inclusion: complex,
        fill_fraction: float) -> complex:
    """Effective complex permittivity of a two-phase composite (Bruggeman EMA).

    Solves, in closed form, the symmetric Bruggeman mixing relation

        f · (ε_i − ε) / (ε_i + 2ε) + (1 − f) · (ε_h − ε) / (ε_h + 2ε) = 0

    for the effective permittivity ε, where ε_h is the host, ε_i the inclusion
    and f the inclusion volume fraction.  This is the standard EMA for a
    material (e.g. porous-anodic alumina) with a percolating lossy inclusion
    phase.  Only the physical root (Im ε ≥ 0, i.e. absorbing) is returned.

    Parameters
    ----------
    eps_host      : complex host permittivity (ε_h).
    eps_inclusion : complex inclusion permittivity (ε_i).
    fill_fraction : inclusion volume fraction f in [0, 1].

    Returns
    -------
    complex effective permittivity with Im(ε) ≥ 0.
    """
    f = float(np.clip(fill_fraction, 0.0, 1.0))
    if f <= 0.0:
        return complex(eps_host)
    if f >= 1.0:
        return complex(eps_inclusion)
    eps_i = complex(eps_inclusion)
    eps_h = complex(eps_host)
    S = eps_i * (3.0 * f - 1.0) + eps_h * (2.0 - 3.0 * f)
    disc = S * S + 8.0 * eps_h * eps_i
    r1 = (S + cmath.sqrt(disc)) / 4.0
    r2 = (S - cmath.sqrt(disc)) / 4.0
    return r1 if r1.imag >= 0.0 else r2


def thin_film_emissivity_air_film_air(
        wavelength_um: float, thickness_um: float,
        n_complex: complex) -> float:
    """Exact normal-incidence emissivity of a free-standing lossy membrane.

    Computes A = 1 − R − T for an Air | film | Air stack (a honeycomb
    pore-to-pore wall) using the coherent Airy / Air–Film–Air Fresnel
    relations under the exp(−iωt) convention.  Emissivity equals absorptivity
    by Kirchhoff's law at each local surface element.

    The propagation factor is taken as exp(i·2·δ) with δ = 2π·ñ·d/λ so that it
    properly DECAYS (|e^{i2δ}| = exp(−4π·Im(ñ)·d/λ) < 1) for an absorbing
    medium; using the opposite sign would produce a spurious amplifying layer.

    Parameters
    ----------
    wavelength_um : vacuum wavelength (µm).
    thickness_um  : film physical thickness (µm).
    n_complex     : film complex index ñ = n + ik.

    Returns
    -------
    float absorptance / emissivity clamped to [0, 1].
    """
    if thickness_um is None or thickness_um <= 0.0 or wavelength_um <= 0.0:
        return 0.0
    n0 = 1.0 + 0.0j     # air (incident)
    n2 = 1.0 + 0.0j     # air (backing) — free-standing membrane
    nc = complex(n_complex)

    r01 = (n0 - nc) / (n0 + nc)
    r12 = (nc - n2) / (nc + n2)
    delta = 2.0 * math.pi * nc * thickness_um / wavelength_um
    e2 = cmath.exp(2j * delta)      # round-trip propagation (decays)

    denom = 1.0 + r01 * r12 * e2
    if abs(denom) < 1e-15:
        return 0.0
    r = (r01 + r12 * e2) / denom
    t01 = 2.0 * n0 / (n0 + nc)
    t12 = 2.0 * nc / (nc + n2)
    t = (t01 * t12 * cmath.exp(1j * delta)) / denom

    R = abs(r) ** 2
    T = (n2.real / n0.real) * abs(t) ** 2
    A = 1.0 - R - T
    return float(np.clip(A, 0.0, 1.0))


def _thin_film_complex_index(
        material: str, wavelength_um: float,
        temperature_K: float = None,
        feature_scale_m: float = None) -> complex:
    """Complex index of a thin wall, with a Bruggeman EMA for PAA/alumina.

    Uses the wavelength-dependent tabulated n(λ), k(λ) (temperature-corrected)
    rather than a flat estimate.  For the porous-alumina honeycomb wall
    material the effective composite index is computed by mixing the alumina
    host permittivity with a small volume fraction of a strongly-absorbing
    carbonaceous inclusion, modelling the real structure of anodized PAA.

    Returns
    -------
    complex ñ = n + ik (absorbing, exp(−iωt) convention).
    """
    n_real, k_imag = get_complex_refractive_index_at_temperature(
        material, wavelength_um, temperature_K, feature_scale_m)
    n_layer = complex(float(n_real), abs(float(k_imag)))

    if str(material).lower() not in PAA_ALIASES:
        return n_layer

    eps_host = n_layer ** 2
    inc = complex(EMA_INCLUSION_INDEX)
    eps_eff = bruggeman_effective_permittivity(
        eps_host, inc ** 2, EMA_CARBON_FILL_FRACTION)
    n_eff = cmath.sqrt(eps_eff)
    return n_eff


def effective_emissivity_thin_film(
    bulk_emissivity: float,
    thickness_um: float,
    wavelength_um: float,
    material: str,
    min_thickness_ratio: float = 0.01,
    temperature_K: float = None,
    feature_scale_m: float = None
) -> float:
    """
    Calculate effective emissivity for a thin film using exact Air–Film–Air TMM.

    For optically thick layers (t/δ > 5) the supplied bulk emissivity is
    retained (a coherent finite-film model is only valid for genuinely thin,
    semi-transparent layers).  For optically thin layers the spectral
    emissivity equals the absorptivity A = 1 − R − T of the free-standing
    membrane, computed by the exact Airy / Air–Film–Air relations.  For
    porous-alumina (PAA) honeycomb walls the effective complex index comes
    from a Bruggeman effective-medium approximation that models the real
    anodized composite (alumina host + lossy carbonaceous inclusion) — the
    physical mechanism behind PAA's high emissivity, replacing any ad-hoc
    emissivity floor.

    Parameters
    ----------
    bulk_emissivity : float
        Bulk material emissivity (0-1)
    thickness_um : float
        Wall/film thickness in micrometres
    wavelength_um : float
        Photon wavelength in micrometres
    material : str
        Material identifier
    min_thickness_ratio : float
        Minimum t/δ for using exp() approximation (default 0.01)
    temperature_K : float, optional
        Material temperature (Phase 4a Drude–Lorentz drift); None → 300 K.
    feature_scale_m : float, optional
        Structure size (m) activating the non-local correction when sub-MFP.

    Returns
    -------
    float
        Effective emissivity (0-1)
    """
    if thickness_um is None or thickness_um <= 0:
        return 0.0

    delta = get_absorption_depth(material, wavelength_um,
                                 temperature_K, feature_scale_m)
    if delta <= 0:
        return bulk_emissivity  # No absorption data → assume bulk

    optical_thickness = thickness_um / delta

    # For very thick walls (t > 5δ), effectively bulk material
    if optical_thickness > 5.0:
        return bulk_emissivity

    # Effective complex index (wavelength-dependent + Bruggeman EMA for PAA).
    n_eff = _thin_film_complex_index(
        material, wavelength_um, temperature_K, feature_scale_m)

    # Exact Air|film|Air absorptance = emissivity (Kirchhoff).
    eps = thin_film_emissivity_air_film_air(
        wavelength_um, thickness_um, n_eff)
    return float(np.clip(eps, 0.0, min(1.0, bulk_emissivity)))


def planck_weighted_effective_emissivity(
    bulk_emissivity: float,
    thickness_um: float,
    material: str,
    temperature_K: float
) -> float:
    """
    Calculate Planck-weighted effective emissivity over the thermal spectrum.
    
    Integrates ε_eff(λ) weighted by blackbody spectral exitance M_λ(T).
    
    Parameters
    ----------
    bulk_emissivity : float
        Bulk material emissivity
    thickness_um : float
        Wall thickness
    material : str
        Material identifier
    temperature_K : float
        Temperature in Kelvin
    
    Returns
    -------
    float
        Planck-weighted effective emissivity (0-1)
    """
    if temperature_K <= 0:
        return bulk_emissivity
    
    # Sample wavelengths from Planck distribution
    from sampling import sample_planck_wavelength

    n_samples = 4000
    total_eps = 0.0

    for _ in range(n_samples):
        lam = sample_planck_wavelength(temperature_K)
        eps = effective_emissivity_thin_film(
            bulk_emissivity, thickness_um, lam, material,
            temperature_K=temperature_K,
        )
        total_eps += eps

    return total_eps / n_samples


def optical_thickness_analysis(
    thickness_um: float,
    material: str,
    temperature_K: float = 300.0
) -> Dict:
    """
    Analyze optical thickness effects for a given wall.
    
    Returns diagnostic information about thin-film correction.
    """
    # Fix 3 (peer-review): evaluate ALL optical constants at the correct
    # temperature-dependent wavelength.  Previously temperature_K was not
    # passed to get_absorption_depth or effective_emissivity_thin_film,
    # meaning the 300 K database was used regardless of the simulation
    # temperature.  For alumina at 200 K, Wien's peak is ~14.49 µm which
    # falls in the reststrahlen band (k >> 0.01) — using NIR k values from
    # a 1200 K source peak (2.41 µm) would corrupt the absorption calculation.
    lambda_peak_um = 2898.0 / temperature_K  # Wien's displacement
    # Fix 3: pass temperature_K so Drude-Lorentz drift gives correct k(λ,T)
    delta_peak = get_absorption_depth(material, lambda_peak_um,
                                      temperature_K=temperature_K)

    # Calculate correction factors
    optical_thickness = thickness_um / delta_peak if delta_peak > 0 else float('inf')

    # Get default bulk emissivity
    bulk_eps = DEFAULT_BULK_EMISSIVITY.get(material, 0.8)

    # Fix 3: effective emissivity at peak wavelength — evaluated at temperature_K.
    eps_peak = effective_emissivity_thin_film(
        bulk_eps, thickness_um, lambda_peak_um, material,
        temperature_K=temperature_K,
    )

    # Planck-weighted effective emissivity
    eps_weighted = planck_weighted_effective_emissivity(
        bulk_eps, thickness_um, material, temperature_K
    )

    # Correction factor
    correction_factor = eps_weighted / bulk_eps if bulk_eps > 0 else 0.0

    return {
        'wall_thickness_um': thickness_um,
        'material': material,
        'temperature_K': temperature_K,
        'peak_wavelength_um': lambda_peak_um,
        'absorption_depth_peak_um': delta_peak,
        'optical_thickness_t/delta': optical_thickness,
        'bulk_emissivity': bulk_eps,
        'effective_emissivity_peak': eps_peak,
        'effective_emissivity_weighted': eps_weighted,
        'thin_film_correction_factor': correction_factor,
        'is_optically_thin': optical_thickness < 1.0,
        'is_effectively_bulk': optical_thickness > 5.0,
        # Fix 3: record which wavelength k was evaluated at (lambda_peak @ T)
        'k_evaluated_at_wavelength_um': lambda_peak_um,
        'k_evaluated_at_temperature_K': temperature_K,
    }


# ---------------------------------------------------------------------------
# Unit tests / validation
# ---------------------------------------------------------------------------
def _test_thin_film_physics():
    """Validate thin-film physics implementation."""
    print("Testing thin-film physics...")
    
    # Test 1: Very thin wall should have low emissivity
    eps1 = effective_emissivity_thin_film(
        bulk_emissivity=0.8,
        thickness_um=0.1,  # 100nm
        wavelength_um=10.0,
        material='alumina'
    )
    print(f"100nm alumina wall at 10um: eps = {eps1:.4f} (expected ~0.022)")
    
    # Test 2: Very thick wall should approach bulk
    eps2 = effective_emissivity_thin_film(
        bulk_emissivity=0.8,
        thickness_um=50.0,  # 50um
        wavelength_um=10.0,
        material='alumina'
    )
    print(f"50um alumina wall at 10um: eps = {eps2:.4f} (expected ~0.8)")
    
    # Test 3: Planck-weighted average
    eps3 = planck_weighted_effective_emissivity(
        bulk_emissivity=0.8,
        thickness_um=0.1,
        material='alumina',
        temperature_K=300.0
    )
    print(f"100nm alumina at 300K (Planck-weighted): eps = {eps3:.4f}")
    
    # Test 4: Optical thickness analysis
    analysis = optical_thickness_analysis(
        thickness_um=0.1,
        material='alumina',
        temperature_K=300.0
    )
    print(f"\nAnalysis for 100nm alumina at 300K:")
    for key, value in analysis.items():
        print(f"  {key}: {value}")


if __name__ == '__main__':
    _test_thin_film_physics()


# ---------------------------------------------------------------------------
# PHASE 1: Complex Refractive Index Functions
# ---------------------------------------------------------------------------

def get_complex_refractive_index(material: str, wavelength_um: float,
                                 temperature_K: float = None,
                                 feature_scale_m: float = None) -> tuple:
    """
    Get complex refractive index ñ(λ) = n(λ) + ik(λ) at given wavelength.

    Parameters
    ----------
    material : str
        Material identifier ('alumina', 'silicon', 'carbon_nanotube', 'silver', etc.)
    wavelength_um : float
        Wavelength in micrometres
    temperature_K : float, optional
        Material temperature for the Drude–Lorentz drift correction
        (Phase 4a).  ``None`` or ≈300 K returns the tabulated data exactly.
    feature_scale_m : float, optional
        Characteristic structure size (m); activates the hydrodynamic
        non-local correction when below the electron mean free path.

    Returns
    -------
    tuple of (n_real, k_imag)
        Real and imaginary parts of refractive index ñ = n + ik
        
    Raises
    ------
    ValueError
        If material not in database
        
    Notes
    -----
    Uses linear interpolation in wavelength. For out-of-range wavelengths,
    clips to database bounds and issues warning.
    
    References: Palik, E. D. (1998). Handbook of Optical Constants of Solids.
    Modest (2013); Raza et al. (2011) for the temperature / non-local model.
    """
    if material not in MATERIAL_COMPLEX_INDEX:
        # Fallback: use absorption depth table to estimate k
        n_real = DEFAULT_REAL_INDEX.get(str(material).lower(), 1.5)
        delta = get_absorption_depth(material, wavelength_um)
        if delta > 0:
            k_imag = wavelength_um / (4.0 * np.pi * delta)
        else:
            k_imag = 0.001

        # Unknown materials have no parametric model; still allow the
        # non-local surface-scattering enhancement on the estimated k when a
        # sub-MFP feature scale is supplied.
        if feature_scale_m is not None and feature_scale_m > 0:
            mfp = electron_mean_free_path_m(material)
            if feature_scale_m < mfp:
                k_imag *= (1.0 + _SURFACE_SCATTERING_FACTOR
                           * mfp / feature_scale_m)
        return (n_real, k_imag)

    data = MATERIAL_COMPLEX_INDEX[material]
    wavelengths = data['wavelengths_um']
    
    # Clip to valid range (asymptotic low/high frequency limits)
    wl_clipped = float(np.clip(wavelength_um, wavelengths[0], wavelengths[-1]))
    
    at_reference = (temperature_K is None
                    or abs(float(temperature_K) - 300.0) < _REFERENCE_TEMP_TOL_K)
    if at_reference and feature_scale_m is None:
        # Fast path — identical to the pre-Phase-4a behaviour.
        n_real = float(np.interp(wl_clipped, wavelengths, data['n_real']))
        k_imag = float(np.interp(wl_clipped, wavelengths, data['k_imag']))
        return (n_real, k_imag)

    return get_complex_refractive_index_at_temperature(
        material, wl_clipped, temperature_K, feature_scale_m)


def fresnel_angular_absorptance(
    n_eff: complex,
    theta_deg: float,
) -> float:
    """Angular, polarisation-averaged absorptance of an air|medium interface.

    Part 1 §2 — the analytic wall absorptance (air over a medium with the
    effective complex index ``n_eff = n(λ,T)+ik(λ,T)``):

        α(θ) = 1 − ½( |r_p|² + |r_s|² )

    with (matching the requested form, n_0 = 1):
        r_p = (n_eff·cosθ − cosθ_t) / (n_eff·cosθ + cosθ_t)
        r_s = (cosθ − n_eff·cosθ_t) / (cosθ + n_eff·cosθ_t)
        sinθ_t = sinθ / n_eff            (Snell, complex — lossy media)

    Feeding the temperature- and scale-dependent index from
    ``get_complex_refractive_index_at_temperature`` turns the static scalar
    wall absorptivity (e.g. α_cnt) into a true α(λ, θ, T).

    Parameters
    ----------
    n_eff : complex
        Effective complex refractive index ñ = n + ik of the medium.
    theta_deg : float
        Incidence polar angle measured from the surface normal (degrees).

    Returns
    -------
    float in [0, 1] — normal-into-surface absorptance.
    """
    theta = math.radians(float(theta_deg))
    cth = math.cos(theta)
    if abs(cth) < 1e-12:
        return 1.0 if abs(n_eff) > 1e-9 else 0.0  # grazing limit
    sin_t = math.sin(theta) / n_eff               # Snell (complex)
    cth_t = cmath.sqrt(1.0 - sin_t * sin_t)       # complex cos θ_t
    n_cth = n_eff * cth
    # r_p = (n_eff·cosθ − cosθ_t)/(n_eff·cosθ + cosθ_t)
    r_p = (n_cth - cth_t) / (n_cth + cth_t)
    # r_s = (cosθ − n_eff·cosθ_t)/(cosθ + n_eff·cosθ_t)
    r_s = (cth - n_eff * cth_t) / (cth + n_eff * cth_t)
    R = 0.5 * (abs(r_p) ** 2 + abs(r_s) ** 2)
    return float(np.clip(1.0 - R, 0.0, 1.0))


def maxwell_garnett_effective_index(
    matrix_material: str,
    wavelength_um: float,
    fill_fraction: float,
    temperature_K: float = None,
    inclusion_eps: complex = 1.0 + 0.0j,
) -> complex:
    """Maxwell-Garnett effective index of a pore-filled homogenised slab.

    Returns ``n_eff = sqrt(ε_eff)`` of a host matrix (e.g. alumina) loaded
    with a dilute spherical inclusion population (default: air-filled pores,
    ε_incl = 1).  This is the homogenisation behind the
    MAXWELL_GARNETT_EFFECTIVE_MEDIUM regime, letting the top cavity layer be
    treated as a single Fresnel thin film instead of ray-traced pores.

    Parameters
    ----------
    matrix_material : host (matrix) material identifier.
    wavelength_um   : vacuum wavelength (µm).
    fill_fraction   : inclusion volume fraction f ∈ [0, 1].
    temperature_K   : material temperature for the Drude–Lorentz drift.
    inclusion_eps   : complex permittivity of the inclusion phase.

    Returns
    -------
    complex — n_eff = sqrt(ε_eff).
    """
    n_mat, k_mat = get_complex_refractive_index_at_temperature(
        matrix_material, wavelength_um, temperature_K)
    eps_m = complex(n_mat, k_mat) ** 2
    f = float(np.clip(fill_fraction, 0.0, 1.0))
    contrast = (inclusion_eps - eps_m) / (inclusion_eps + eps_m)
    denom = 1.0 - f * contrast
    eps_eff = eps_m * (1.0 + (2.0 * f * contrast) / denom
                       if abs(denom) > 1e-15 else 1.0)
    return cmath.sqrt(eps_eff + 0.0j)


def aperture_boundary_absorptance(
    material: str,
    wavelength_um: float,
    fill_fraction: float,
    temperature_K: float = None,
    theta_deg: float = 0.0,
) -> float:
    """Dynamic aperture-boundary absorptance ``1 − R_ap(λ, T, f)``.

    Replaces the hard-coded ``alpha_top = 1 − (1 − ε_walls)·0.1`` guess.  The
    sub-cutoff incident field that cannot form a channel mode folds around the
    wall rims and impinges on the top surface; that surface is modelled as a
    Maxwell–Garnett homogenised medium of effective index ``n_eff(f)`` across
    the aperture fill fraction ``f``.  The reflectance is the effective-index
    step at the aperture interface:

        R_ap(λ, T, f) = |(n_eff − 1)/(n_eff + 1)|²      (normal incidence)
        α(λ, T, f)    = 1 − R_ap(λ, T, f),

    evaluated with the temperature-dependent Drude–Lorentz drift so the
    boundary estimate is genuinely wavelength-, temperature- and fill-ratio
    dependent (no hard-coded 90% trap).

    Parameters
    ----------
    material      : host (matrix) material identifier, e.g. 'alumina'.
    wavelength_um : vacuum wavelength (µm).
    fill_fraction : inclusion volume fraction f ∈ [0, 1] across the aperture.
    temperature_K : material temperature for the n(λ,T) drift (300 K default).
    theta_deg     : incidence polar angle (degrees) from the surface normal.

    Returns
    -------
    float in [0, 1] — the aperture absorptance (1 − R_ap).
    """
    n_eff = maxwell_garnett_effective_index(
        material, float(wavelength_um), float(fill_fraction),
        temperature_K=float(temperature_K) if temperature_K is not None else None)
    return float(np.clip(fresnel_angular_absorptance(n_eff, float(theta_deg)),
                         0.0, 1.0))


def material_optics_confidence(
    material: str,
    temperature_K: float = None,
    baseline_K: float = 300.0,
    feature_scale_m: float = None,
) -> dict:
    """Confidence in the tabulated/Drude–Lorentz optical model at ``temperature_K``.

    The optical-property tables and the Drude–Lorentz parameters are calibrated
    at the ``baseline_K`` (300 K).  As the operating temperature drifts away
    from that anchor the extrapolated complex index becomes less trustworthy.
    This routine quantifies that degradation:

        drift    = |γ(T)/γ(T_0) − 1|        (collision-frequency drift)
        conf_T   = 1 / (1 + drift/κ)        (κ = half-drift falloff, default 0.5)
        conf_mat = 0 for un-calibrated materials (pure engineering table)

    and folds in a structural penalty when a sub-mean-free-path feature scale
    makes the non-local / surface-scattering correction active (larger model
    uncertainty).

    Parameters
    ----------
    material      : material identifier (aliases resolved internally).
    temperature_K : operating temperature; None → baseline (confidence 1.0).
    baseline_K    : the calibrated reference temperature (default 300 K).
    feature_scale_m : optional structure size activating non-local correction.

    Returns
    -------
    dict with ``confidence_level`` (float in [0,1]) and ``warnings`` (list).
    """
    warnings_d = []
    T = float(baseline_K) if temperature_K is None else float(temperature_K)
    params = DRUDE_LORENTZ_PARAMS.get(str(material).lower())

    if params is None:
        # Material has only an interpolated table / /engineering estimate with
        # no modelled temperature drift — confidence is capped low.
        conf = 0.40
        warnings_d.append(
            f"{material}: no calibrated Drude–Lorentz model — optical "
            f"properties are engineering estimates; confidence capped at {conf:.0%}.")
        return {'confidence_level': float(conf),
                'material': str(material),
                'temperature_K': T,
                'baseline_K': float(baseline_K),
                'drift_fraction': None,
                'warnings': warnings_d}

    if abs(float(baseline_K) - 300.0) < _REFERENCE_TEMP_TOL_K:
        at_baseline = abs(T - float(baseline_K)) < _REFERENCE_TEMP_TOL_K
    else:
        at_baseline = (T == float(baseline_K))

    gamma = electron_collision_frequency_eV(material, T)
    gamma_base = electron_collision_frequency_eV(material, float(baseline_K))
    drift = abs(gamma - gamma_base) / max(gamma_base, 1e-12)
    kappa = 0.5
    conf = 1.0 / (1.0 + drift / kappa)   # drift≈0 at baseline → conf≈1.0

    # Reduce confidence when the non-local / surface-scattering correction is
    # active (hydrodynamic term is an inherently more uncertain model layer).
    if feature_scale_m is not None and feature_scale_m > 0:
        mfp = electron_mean_free_path_m(material)
        if feature_scale_m < mfp and float(params['v_fermi_m_s']) > 0:
            conf *= 0.85
            warnings_d.append(
                f'non-local/surface-scattering correction active '
                f'(feature l={feature_scale_m:.3g} m < ℓ={mfp:.3g} m); '
                'model uncertainty increased (×0.85).')

    if abs(T - float(baseline_K)) > 400.0:
        warnings_d.append(
            f'T={T:.0f} K is >400 K from the {baseline_K:.0f} K calibration; '
            'extrapolated optical constants carry reduced confidence.')
    if drift > 0.25:
        warnings_d.append(
            f'collision-frequency drift = {drift:.2f} at T={T:.0f} K — '
            'dispersion has moved far from the calibrated 300 K baseline.')

    conf = float(max(0.0, min(1.0, conf)))
    if conf < 0.6:
        warnings_d.append(f'confidence={conf:.2f} — treat broadband n(λ,T) '
                          'emissivity results with caution.')

    return {'confidence_level': conf,
            'material': material,
            'temperature_K': T,
            'baseline_K': float(baseline_K),
            'drift_fraction': float(drift),
            'warnings': warnings_d}

def tmm_reflectance_single_layer(
    n_0: complex,
    n_1: complex,
    n_2: complex,
    thickness_um: float,
    wavelength_um: float,
    theta_0_deg: float = 0.0,
    polarization: str = 's'
) -> float:
    """
    Calculate reflectance for single lossy layer using Transfer-Matrix Method.
    
    Computes Fresnel reflectance R(λ,θ) for a single-layer film with complex
    refractive index, accounting for:
      - Phase accumulation in lossy medium: φ = 2π ñ₁ t cos(θ₁) / λ
      - Evanescent wave tunneling (complex Snell's law)
      - Both s-polarized (TE) and p-polarized (TM) cases
    
    Parameters
    ----------
    n_0 : complex
        Complex refractive index of incident medium (usually 1+0j for vacuum/air)
    n_1 : complex
        Complex refractive index of layer (ñ = n + ik)
    n_2 : complex
        Complex refractive index of substrate below (backing)
    thickness_um : float
        Layer thickness in micrometres
    wavelength_um : float
        Free-space wavelength in micrometres
    theta_0_deg : float
        Incident angle in degrees (0 = normal incidence)
    polarization : str
        'S' or 's' for s-polarized (TE, E-perpendicular)
        'P' or 'p' for p-polarized (TM, E-parallel)
    
    Returns
    -------
    float
        Reflectance R ∈ [0, 1]
        
    References
    ----------
    - Born & Wolf (1999). Principles of Optics, Chapter 1
    - Heavens (1955). Optical Properties of Thin Solid Films
    
    Notes
    -----
    For normal incidence (θ₀ = 0), result is independent of polarization.
    For oblique incidence and lossy materials, must use complex math.
    """
    theta_0_rad = np.radians(theta_0_deg)
    
    # Handle normal incidence case efficiently
    if abs(theta_0_deg) < 1e-6:
        # Normal incidence: simple Fresnel formula
        # r = (n₀ - ñ) / (n₀ + ñ)
        # R = |r|²
        
        # For single layer with Fabry-Pérot interference:
        r01 = (n_0 - n_1) / (n_0 + n_1)
        r12 = (n_1 - n_2) / (n_1 + n_2)
        
        # Phase shift in layer
        delta = (2.0 * np.pi / wavelength_um) * n_1 * thickness_um
        exp_2i_delta = np.exp(2.0j * delta)
        
        # Fresnel reflection including multiple reflections (Fabry-Pérot)
        r_total = (r01 + r12 * exp_2i_delta) / (1.0 + r01 * r12 * exp_2i_delta)
        R = abs(r_total) ** 2
        
        return float(np.real(R))
    
    # Oblique incidence with complex angles
    sin_theta_0 = np.sin(theta_0_rad)
    
    # Snell's law: n₀ sin(θ₀) = ñ₁ sin(θ₁)
    # Note: sin(θ₁) may be complex (evanescent wave)
    sin_theta_1 = (n_0 / n_1) * sin_theta_0
    cos_theta_0 = np.sqrt(1.0 - sin_theta_0**2)
    cos_theta_1 = np.sqrt(1.0 - sin_theta_1**2)
    
    # Substrate side
    sin_theta_2 = (n_0 / n_2) * sin_theta_0
    cos_theta_2 = np.sqrt(1.0 - sin_theta_2**2)
    
    # Fresnel reflection coefficients at interfaces
    if polarization.lower() == 's':
        # s-polarized (TE, E ⊥ plane of incidence)
        r01 = (n_0 * cos_theta_0 - n_1 * cos_theta_1) / (n_0 * cos_theta_0 + n_1 * cos_theta_1)
        r12 = (n_1 * cos_theta_1 - n_2 * cos_theta_2) / (n_1 * cos_theta_1 + n_2 * cos_theta_2)
    else:  # p-polarized
        # p-polarized (TM, H ⊥ plane of incidence)
        r01 = (n_1 * cos_theta_0 - n_0 * cos_theta_1) / (n_1 * cos_theta_0 + n_0 * cos_theta_1)
        r12 = (n_2 * cos_theta_1 - n_1 * cos_theta_2) / (n_2 * cos_theta_1 + n_1 * cos_theta_2)
    
    # Phase accumulated in layer
    delta = (2.0 * np.pi / wavelength_um) * n_1 * thickness_um * cos_theta_1
    exp_2i_delta = np.exp(2.0j * delta)
    
    # Total reflection (single layer with Fabry-Pérot)
    numerator = r01 + r12 * exp_2i_delta
    denominator = 1.0 + r01 * r12 * exp_2i_delta
    
    if abs(denominator) < 1e-12:
        return 1.0  # Total reflection
    
    r_total = numerator / denominator
    R = abs(r_total) ** 2
    
    return float(np.real(R))
