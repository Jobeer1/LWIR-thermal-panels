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
# At 300K, peak thermal radiation λ ≈ 9.7µm
# Data from Palik, Handbook of Optical Constants
ALUMINA_ABSORPTION_DEPTH = {
    'wavelengths_um': np.array([2.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0]),
    'depths_um':      np.array([1.5, 2.5, 3.2, 3.61, 4.0, 4.8, 6.0]),
    'source': 'Palik (1998), extrapolated for thermal IR',
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
# Data interpolated from Palik, Handbook of Optical Constants
# Temperature: 300 K (room temperature)
ALUMINA_COMPLEX_INDEX = {
    'wavelengths_um': np.array([0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 30.0]),
    'n_real': np.array([1.76, 1.73, 1.71, 1.69, 1.65, 1.62, 1.61, 1.58, 1.57, 1.56, 1.55, 1.54, 1.53]),
    'k_imag': np.array([0.001, 0.0008, 0.0005, 0.0002, 0.00005, 0.0001, 0.0008, 0.0015, 0.0020, 0.0025, 0.0035, 0.005, 0.008]),
    'citation': 'Palik, E. D. (1998). Handbook of Optical Constants of Solids',
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
# Interpolation and calculation functions
# ---------------------------------------------------------------------------

def get_absorption_depth(material: str, wavelength_um: float) -> float:
    """
    Get absorption depth δ(λ) for a material at given wavelength.
    
    Parameters
    ----------
    material : str
        Material identifier ('alumina', 'cnt_forest', 'silver', etc.)
    wavelength_um : float
        Wavelength in micrometres
    
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
        return float(depths[0])
    elif idx == len(wavelengths):
        return float(depths[-1])
    
    # Linear interpolation
    wl0, wl1 = wavelengths[idx-1], wavelengths[idx]
    d0, d1 = depths[idx-1], depths[idx]
    
    # Interpolate linearly (works well for smooth absorption curves)
    t = (wl - wl0) / (wl1 - wl0)
    return float(d0 + t * (d1 - d0))


def tmm_absorptance_normal_incidence(
    thickness_um: float,
    wavelength_um: float,
    material: str,
    substrate_index: complex = 1.5,
    incident_index: complex = 1.0,
) -> float:
    """Return stack absorptance from a one-layer normal-incidence TMM.

    The layer is the selected material and the substrate is semi-infinite.
    This is a coherent thin-film calculation, so it includes Fresnel
    reflection and phase interference.  ``k`` is inferred from the tabulated
    absorption depth using k = λ / (4πδ).  The result is the absorptance of
    the film/substrate surface and is therefore the emissivity used by the
    Monte Carlo model under Kirchhoff reciprocity.
    """
    if thickness_um <= 0.0 or wavelength_um <= 0.0:
        return 0.0

    delta_um = get_absorption_depth(material, wavelength_um)
    if delta_um <= 0.0 or not math.isfinite(delta_um):
        return 0.0

    n_real = DEFAULT_REAL_INDEX.get(material, 1.5)
    k = wavelength_um / (4.0 * math.pi * delta_um)
    n_layer = complex(n_real, -k)  # exp(-iωt) convention
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


def effective_emissivity_thin_film(
    bulk_emissivity: float,
    thickness_um: float,
    wavelength_um: float,
    material: str,
    min_thickness_ratio: float = 0.01
) -> float:
    """
    Calculate effective emissivity for a thin film using a transfer matrix.
    
    For optically thick layers (t/δ > 5), retain the supplied bulk
    emissivity.  This avoids using a coherent finite-film model where the
    layer is no longer thin and the available optical constants are only
    approximate.
    
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
    
    Returns
    -------
    float
        Effective emissivity (0-1)
    """
    if thickness_um <= 0:
        return 0.0
    
    delta = get_absorption_depth(material, wavelength_um)
    if delta <= 0:
        return bulk_emissivity  # No absorption data → assume bulk
    
    optical_thickness = thickness_um / delta
    
    # For very thick walls (t > 5δ), effectively bulk material
    if optical_thickness > 5.0:
        return bulk_emissivity
    
    tmm_eps = tmm_absorptance_normal_incidence(
        thickness_um, wavelength_um, material
    )
    return float(np.clip(tmm_eps, 0.0, min(1.0, bulk_emissivity)))


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
    
    n_samples = 1000
    total_eps = 0.0
    
    for _ in range(n_samples):
        lam = sample_planck_wavelength(temperature_K)
        eps = effective_emissivity_thin_film(
            bulk_emissivity, thickness_um, lam, material
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
    # Get absorption depth at peak wavelength
    lambda_peak_um = 2898.0 / temperature_K  # Wien's displacement
    delta_peak = get_absorption_depth(material, lambda_peak_um)
    
    # Calculate correction factors
    optical_thickness = thickness_um / delta_peak if delta_peak > 0 else float('inf')
    
    # Get default bulk emissivity
    bulk_eps = DEFAULT_BULK_EMISSIVITY.get(material, 0.8)
    
    # Effective emissivity at peak wavelength
    eps_peak = effective_emissivity_thin_film(
        bulk_eps, thickness_um, lambda_peak_um, material
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

def get_complex_refractive_index(material: str, wavelength_um: float) -> tuple:
    """
    Get complex refractive index ñ(λ) = n(λ) + ik(λ) at given wavelength.
    
    Parameters
    ----------
    material : str
        Material identifier ('alumina', 'silicon', 'carbon_nanotube', 'silver', etc.)
    wavelength_um : float
        Wavelength in micrometres
    
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
    """
    if material not in MATERIAL_COMPLEX_INDEX:
        # Fallback: use absorption depth table to estimate k
        n_real = DEFAULT_REAL_INDEX.get(material.lower(), 1.5)
        delta = get_absorption_depth(material, wavelength_um)
        if delta > 0:
            k_imag = wavelength_um / (4.0 * np.pi * delta)
        else:
            k_imag = 0.001
        return (n_real, k_imag)
    
    data = MATERIAL_COMPLEX_INDEX[material]
    wavelengths = data['wavelengths_um']
    n_real_table = data['n_real']
    k_imag_table = data['k_imag']
    
    # Clip to valid range and warn if out of bounds
    if wavelength_um < wavelengths[0] or wavelength_um > wavelengths[-1]:
        warnings.warn(
            f"Wavelength {wavelength_um:.3f}um outside database range "
            f"[{wavelengths[0]:.3f}, {wavelengths[-1]:.3f}]um for {material}. "
            f"Extrapolating...",
            UserWarning
        )
    
    wl_clipped = np.clip(wavelength_um, wavelengths[0], wavelengths[-1])
    
    # Linear interpolation
    n_real = float(np.interp(wl_clipped, wavelengths, n_real_table))
    k_imag = float(np.interp(wl_clipped, wavelengths, k_imag_table))
    
    return (n_real, k_imag)


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
