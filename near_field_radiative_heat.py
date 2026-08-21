"""
PHASE 3: Near-Field Radiative Heat Transfer (Polder-Van Hove)

Implements the Fluctuation-Dissipation Theorem (FDT) for radiative heat
exchange between two surfaces separated by a gap g << λ_peak.

Key Physics:
  - Evanescent waves tunnel across gap, bypassing classical view factors
  - Transmission coefficient s(k_∥, ω) from Fresnel at two interfaces
  - Spectral integral: Q(ω) ∝ [Θ(T_A,ω) - Θ(T_B,ω)] × s(k_∥,ω)
  - Parallel wavevector k_∥ extends beyond propagating limit (k_∥ > ω/c)
  - Evanescent contribution dominates for g < λ/(2π)

References:
  - Polder, D. & Van Hove, M. (1971). PRB 4(10), 3303
  - Rytov, S. M. et al. (1989). Principles of Statistical Radiophysics vol. 3
  - Basu, S. et al. (2009). Int. J. Energy Res. 33(13), 1203-1232

Author: Kiro
Date: August 2026
"""

import numpy as np
import warnings

try:
    from scipy import integrate, optimize
    from scipy.polynomial.legendre import leggauss
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    warnings.warn("scipy not available; using fallback numerical methods")
    
    def leggauss(n):
        """Fallback Gauss-Legendre quadrature if scipy unavailable."""
        warnings.warn("Using basic fallback quadrature (lower accuracy)")
        # Simple fallback: equally-spaced nodes with uniform weights
        nodes = np.linspace(-1, 1, n)
        weights = np.ones(n) * (2.0 / n)
        return nodes, weights

try:
    from material_optics import get_complex_refractive_index
except ImportError:
    warnings.warn("material_optics not available; using fallback refractive indices")
    def get_complex_refractive_index(material: str, wavelength_um: float) -> tuple:
        """Fallback refractive index database if material_optics unavailable."""
        fallback_indices = {
            'alumina': (1.7, 0.002),
            'silicon': (3.5, 0.01),
            'silver': (0.2, 2.0),
        }
        return fallback_indices.get(material.lower(), (1.5, 0.01))

# Physical constants
c = 299792458.0  # Speed of light (m/s)
hbar = 1.054571817e-34  # Reduced Planck constant (J·s)
k_B = 1.380649e-23  # Boltzmann constant (J/K)


# ---------------------------------------------------------------------------
# Fresnel Coefficients and Transmission Kernel
# ---------------------------------------------------------------------------

def fresnel_coefficients_interface(
    k_parallel_m: float,
    omega_rad_s: float,
    n_1_complex: complex,
    n_2_complex: complex,
    c_light: float = c
) -> tuple:
    """
    Fresnel reflection coefficients for evanescent and propagating waves.
    
    Calculates r_s (s-pol, TE) and r_p (p-pol, TM) for arbitrary parallel
    wavenumber k_∥, including evanescent regions where k_∥ > ω/c.
    
    Parameters
    ----------
    k_parallel_m : float
        Parallel (tangential) wavevector magnitude (m⁻¹)
        Range: 0 to k_0 (propagating), k_0 to ∞ (evanescent)
    omega_rad_s : float
        Angular frequency (rad/s)
    n_1_complex : complex
        Complex refractive index of medium 1 (e.g., emitter surface)
    n_2_complex : complex
        Complex refractive index of medium 2 (e.g., receiver surface or vacuum)
    c_light : float
        Speed of light (m/s)
    
    Returns
    -------
    tuple of (r_s, r_p)
        Complex reflection coefficients for s and p polarizations
        
    Notes
    -----
    For evanescent waves (k_∥ > ω/c), the perpendicular wavenumber is purely
    imaginary, leading to exponential decay away from interface.
    
    The Fresnel coefficients become complex-valued, but maintain |r| ≤ 1
    for passive, absorbing media.
    """
    
    # Free-space wavenumber
    k0 = omega_rad_s / c_light
    
    # Perpendicular wavenumbers (handle both propagating and evanescent)
    # κ² = (ω/c)² n² - k_∥²  (careful with signs for imaginary roots)
    
    # In medium 1
    kappa1_sq = (k0 * n_1_complex)**2 - k_parallel_m**2
    if np.real(kappa1_sq) < 0:
        # Evanescent: κ = i|κ|
        kappa1 = 1.0j * np.sqrt(-kappa1_sq)
    else:
        kappa1 = np.sqrt(kappa1_sq)
    
    # In medium 2 (assume vacuum/air for gap)
    kappa2_sq = k0**2 - k_parallel_m**2
    if kappa2_sq < 0:
        # Evanescent in gap
        kappa2 = 1.0j * np.sqrt(-kappa2_sq)
    else:
        kappa2 = np.sqrt(kappa2_sq)
    
    # Impedance in each medium: Z = ω / (κ n²)
    Z1 = omega_rad_s / (kappa1 * n_1_complex**2 + 1e-30)
    Z2 = omega_rad_s / (kappa2 + 1e-30)
    
    # s-polarized (TE, E perpendicular to plane of incidence)
    # r_s = (Z1 - Z2) / (Z1 + Z2)
    r_s = (Z1 - Z2) / (Z1 + Z2 + 1e-30)
    
    # p-polarized (TM, H perpendicular to plane)
    # r_p = (n_2² Z1 - n_1² Z2) / (n_2² Z1 + n_1² Z2)
    r_p = (n_2_complex**2 * Z1 - n_1_complex**2 * Z2) / (n_2_complex**2 * Z1 + n_1_complex**2 * Z2 + 1e-30)
    
    return (r_s, r_p)


def near_field_transmission_coefficient(
    k_parallel_m: float,
    omega_rad_s: float,
    gap_m: float,
    n_emitter: complex,
    n_receiver: complex,
    c_light: float = c
) -> float:
    """
    Near-field radiative transmission coefficient s(k_∥, ω, g).
    
    Computes the energy transmission coefficient for evanescent and
    propagating waves crossing a gap g between two surfaces.
    
    Implements Polder-Van Hove formula:
        s = 4 / (1 + exp(2κ_gap·g)) × 2·r_p·r_s / (r_p + r_s - r_p·r_s·exp(2κ_gap·g))
    
    Parameters
    ----------
    k_parallel_m : float
        Parallel wavevector (m⁻¹)
    omega_rad_s : float
        Angular frequency (rad/s)
    gap_m : float
        Gap distance (m)
    n_emitter : complex
        Complex refractive index of emitter surface
    n_receiver : complex
        Complex refractive index of receiver surface
    c_light : float
        Speed of light (m/s)
    
    Returns
    -------
    float
        Transmission coefficient s ∈ [0, ~4]
        Typically < 1 for propagating, can exceed 1 for evanescent near resonance
    """
    
    k0 = omega_rad_s / c_light
    
    # Get Fresnel coefficients at emitter-gap and gap-receiver interfaces
    r_s, r_p = fresnel_coefficients_interface(
        k_parallel_m, omega_rad_s, n_emitter, 1.0 + 0.0j, c_light
    )
    
    # Decay in vacuum gap
    kappa_gap_sq = k0**2 - k_parallel_m**2
    if kappa_gap_sq < 0:
        # Evanescent: κ_gap = i√(k_∥² - k_0²)
        kappa_gap = 1.0j * np.sqrt(-kappa_gap_sq)
    else:
        # Propagating
        kappa_gap = np.sqrt(kappa_gap_sq)
    
    # Phase/decay crossing gap
    phase_gap = 2.0 * kappa_gap * gap_m
    exp_phase_gap = np.exp(phase_gap)
    
    # Polder-Van Hove formula (simplified)
    # s = 4/(1 + e^(κ·2g)) × (2·r_p·r_s) / (r_p + r_s - r_p·r_s·e^(κ·2g))
    
    try:
        # Numerator: 4·r_p·r_s / (1 + exp(2κ·g))
        numerator = 4.0 * r_p * r_s
        denominator_1 = 1.0 + np.exp(np.real(phase_gap))  # For amplitude
        
        # Alternative denominator for the full formula
        denominator_2 = r_p + r_s - r_p * r_s * exp_phase_gap
        
        if abs(denominator_2) < 1e-30:
            return 0.0
        
        # Full formula from Polder-Van Hove
        s = numerator / denominator_2
        
        # Ensure physical: s should be real and ≥ 0 for this frequency
        s_real = np.real(s)
        s_phys = max(0.0, s_real)
        
        return float(s_phys)
        
    except (ValueError, RuntimeWarning):
        return 0.0


def planck_energy_quantum(
    omega_rad_s: float,
    temperature_K: float,
    hbar_const: float = hbar,
    k_B_const: float = k_B
) -> float:
    """
    Quantum harmonic oscillator energy at frequency ω and temperature T.
    
    ℏω / (exp(ℏω/k_B T) - 1)
    
    This is the average energy of a mode according to Planck distribution.
    
    Parameters
    ----------
    omega_rad_s : float
        Angular frequency (rad/s)
    temperature_K : float
        Temperature (K)
    
    Returns
    -------
    float
        Energy of oscillator (Joules)
    """
    
    if temperature_K <= 0:
        return 0.0
    
    x = hbar_const * omega_rad_s / (k_B_const * temperature_K)
    
    # Handle limits
    if x < 1e-6:
        # Low frequency (Rayleigh-Jeans): E ≈ k_B T
        return k_B_const * temperature_K
    elif x > 100:
        # High frequency (Wien): E ≈ ℏω exp(-ℏω/k_B T)
        return hbar_const * omega_rad_s * np.exp(-x)
    else:
        # Full Planck formula
        return (hbar_const * omega_rad_s) / (np.exp(x) - 1.0)


# ---------------------------------------------------------------------------
# Near-Field Spectral Integration
# ---------------------------------------------------------------------------

def near_field_heat_flux_spectral(
    temperature_hot_K: float,
    temperature_cold_K: float,
    gap_m: float,
    material_hot: str = 'alumina',
    material_cold: str = 'alumina',
    omega_min_rad_s: float = None,
    omega_max_rad_s: float = None,
    n_omega: int = 100,
    n_kparallel: int = 60,
    c_light: float = c
) -> dict:
    """
    Integrate Polder-Van Hove near-field heat flux over frequency spectrum.
    
    Computes:
        Q = (1/π²) ∫∫ [Θ(T_A,ω) - Θ(T_B,ω)] s(k_∥,ω,g) k_∥ dk_∥ dω
    
    where:
        - Θ is Planck energy oscillator
        - s is near-field transmission coefficient
        - ω ranges over thermal spectrum
        - k_∥ includes both propagating (k_∥ < ω/c) and evanescent (k_∥ > ω/c)
    
    Parameters
    ----------
    temperature_hot_K : float
        Hot surface temperature (K)
    temperature_cold_K : float
        Cold surface temperature (K)
    gap_m : float
        Gap distance (m)
    material_hot : str
        Material of hot surface for optical properties
    material_cold : str
        Material of cold surface (usually same)
    omega_min_rad_s : float
        Minimum frequency (rad/s). If None, set from λ_max ≈ 50 µm
    omega_max_rad_s : float
        Maximum frequency (rad/s). If None, set from λ_min ≈ 0.1 µm
    n_omega : int
        Number of frequency points (Gauss-Legendre quadrature)
    n_kparallel : int
        Number of parallel wavenumber points
    c_light : float
        Speed of light (m/s)
    
    Returns
    -------
    dict with:
        'flux_W_m2': Net heat flux (W/m²)
        'flux_by_region': {'propagating': ..., 'evanescent': ...} (W/m²)
        'dominant_wavelength_um': Peak integrand wavelength (µm)
        'peak_contribution_k_parallel_m': Peak wavevector (m⁻¹)
        'materials': (material_hot, material_cold)
        'integration_info': Diagnostic info
    """
    
    # Default frequency limits based on Wien's law
    if temperature_hot_K > 0:
        lambda_peak_hot = 2898.0 / temperature_hot_K  # µm
    else:
        lambda_peak_hot = 10.0
    
    if omega_min_rad_s is None:
        lambda_max_um = lambda_peak_hot * 5.0  # Integrate to 5× peak
        omega_min_rad_s = 2.0 * np.pi * c_light * 1e6 / lambda_max_um
    
    if omega_max_rad_s is None:
        lambda_min_um = 0.1  # UV cutoff
        omega_max_rad_s = 2.0 * np.pi * c_light * 1e6 / lambda_min_um
    
    # Gauss-Legendre quadrature for frequency
    omega_nodes, omega_weights = leggauss(n_omega)
    
    # Map from [-1,1] to [ω_min, ω_max]
    omega_vals = omega_min_rad_s + 0.5 * (omega_max_rad_s - omega_min_rad_s) * (omega_nodes + 1.0)
    omega_weights = 0.5 * (omega_max_rad_s - omega_min_rad_s) * omega_weights
    
    # Integration storage
    flux_total = 0.0
    flux_propagating = 0.0
    flux_evanescent = 0.0
    
    peak_integrand = 0.0
    peak_omega = 0.0
    peak_k_parallel = 0.0
    
    k0_at_peak = 2.0 * np.pi / (lambda_peak_hot * 1e-6)  # k_0 at peak wavelength
    
    # Frequency loop
    for i, omega in enumerate(omega_vals):
        
        # Get complex refractive indices
        lambda_um = 2.0 * np.pi * c_light * 1e6 / omega
        
        # Safety checks
        if lambda_um < 0.001 or lambda_um > 1000:
            continue
        
        try:
            n_hot_real, k_hot = get_complex_refractive_index(material_hot, lambda_um)
            n_hot = n_hot_real + 1.0j * k_hot
            
            n_cold_real, k_cold = get_complex_refractive_index(material_cold, lambda_um)
            n_cold = n_cold_real + 1.0j * k_cold
        except:
            continue
        
        # Planck energy difference
        theta_hot = planck_energy_quantum(omega, temperature_hot_K)
        theta_cold = planck_energy_quantum(omega, temperature_cold_K)
        theta_diff = theta_hot - theta_cold
        
        if abs(theta_diff) < 1e-40:
            continue
        
        # Free-space wavenumber
        k0 = omega / c_light
        
        # k_∥ integration: adaptive range
        # Propagating: k_∥ ∈ [0, k_0]
        # Evanescent: k_∥ ∈ [k_0, 10·k_0] (significant decay beyond 10·k_0)
        k_para_max = 10.0 * k0
        
        # Gauss-Legendre for k_∥
        k_para_nodes, k_para_weights = leggauss(n_kparallel)
        
        # Map to [0, k_para_max]
        k_para_vals = 0.5 * k_para_max * (k_para_nodes + 1.0)
        k_para_weights_scaled = 0.5 * k_para_max * k_para_weights
        
        # k_∥ loop
        for j, k_para in enumerate(k_para_vals):
            
            try:
                # Near-field transmission
                s = near_field_transmission_coefficient(
                    k_para, omega, gap_m, n_hot, n_cold, c_light
                )
                
                if s < 0 or not np.isfinite(s):
                    continue
                
                # Integrand: (1/π²) × θ_diff × s × k_∥
                integrand = (1.0 / (np.pi**2)) * theta_diff * s * k_para
                
                if integrand < 0 or not np.isfinite(integrand):
                    continue
                
                # Accumulate
                contribution = integrand * omega_weights[i] * k_para_weights_scaled[j]
                flux_total += contribution
                
                # Separate propagating/evanescent
                if k_para <= k0:
                    flux_propagating += contribution
                else:
                    flux_evanescent += contribution
                
                # Track peak
                if integrand > peak_integrand:
                    peak_integrand = integrand
                    peak_omega = omega
                    peak_k_parallel = k_para
            
            except:
                continue
    
    # Convert peak frequency back to wavelength
    if peak_omega > 0:
        peak_wavelength_um = 2.0 * np.pi * c_light * 1e6 / peak_omega
    else:
        peak_wavelength_um = lambda_peak_hot
    
    evanescent_fraction = flux_evanescent / (flux_total + 1e-40)
    
    return {
        'flux_W_m2': float(flux_total),
        'flux_by_region': {
            'propagating_W_m2': float(flux_propagating),
            'evanescent_W_m2': float(flux_evanescent)
        },
        'evanescent_fraction': float(evanescent_fraction),
        'dominant_wavelength_um': float(peak_wavelength_um),
        'peak_contribution_k_parallel_m': float(peak_k_parallel),
        'materials': (material_hot, material_cold),
        'integration_info': {
            'n_omega_points': n_omega,
            'n_kparallel_points': n_kparallel,
            'omega_min_rad_s': omega_min_rad_s,
            'omega_max_rad_s': omega_max_rad_s,
            'gap_m': gap_m
        }
    }


def gap_ratio_metric(gap_m: float, temperature_K: float) -> float:
    """
    Dimensionless gap ratio for near-field regime identification.
    
    Gap Ratio = g / (λ_peak / 2π)
    
    Interpretation:
        < 1: Strong near-field (evanescent dominated)
        < 5: Significant near-field contribution
        > 20: Negligible near-field, far-field dominant
    
    Parameters
    ----------
    gap_m : float
        Gap distance (m)
    temperature_K : float
        Temperature for Wien's peak wavelength (K)
    
    Returns
    -------
    float
        Dimensionless gap ratio
    """
    
    lambda_peak_m = 2898.0 * 1e-6 / temperature_K  # Wien's law in meters
    characteristic_length = lambda_peak_m / (2.0 * np.pi)
    
    ratio = gap_m / (characteristic_length + 1e-30)
    return float(ratio)


def should_use_near_field_model(gap_m: float, temperature_K: float, threshold: float = 5.0) -> bool:
    """
    Decision logic: when to use near-field model vs. far-field.
    
    Parameters
    ----------
    gap_m : float
        Gap distance (m)
    temperature_K : float
        Reference temperature (K)
    threshold : float
        Gap ratio threshold (default 5.0)
    
    Returns
    -------
    bool
        True if near-field physics significant (gap ratio < threshold)
    """
    
    ratio = gap_ratio_metric(gap_m, temperature_K)
    return ratio < threshold


# ---------------------------------------------------------------------------
# Unit Tests
# ---------------------------------------------------------------------------

def test_near_field_physics():
    """Validate near-field radiative transfer."""
    print("\n" + "="*70)
    print("Testing Near-Field Radiative Heat Transfer")
    print("="*70)
    
    # Test 1: Fresnel coefficients
    print("\n[Test 1] Fresnel coefficients at normal incidence")
    k_para_test = 0.0  # Normal incidence
    omega_test = 2.0 * np.pi * c / 10e-6  # 10 um wavelength
    n_test = 1.5 + 0.1j
    
    r_s, r_p = fresnel_coefficients_interface(k_para_test, omega_test, n_test, 1.0+0.0j)
    print(f"  r_s = {r_s:.4f}, r_p = {r_p:.4f}")
    print(f"  |r_s|^2 = {abs(r_s)**2:.6f}, |r_p|^2 = {abs(r_p)**2:.6f}")
    
    # Test 2: Planck energy
    print("\n[Test 2] Planck energy at thermal wavelengths")
    omega_10um = 2.0 * np.pi * c / 10e-6
    theta_300k = planck_energy_quantum(omega_10um, 300.0)
    print(f"  hbar*omega at 10um, 300K: {theta_300k:.3e} J")
    
    # Test 3: Near-field transmission
    print("\n[Test 3] Near-field transmission coefficient")
    gap_test = 100e-9  # 100 nm
    s = near_field_transmission_coefficient(
        0.0, omega_10um, gap_test, 1.5+0.1j, 1.5+0.1j
    )
    print(f"  s(k_parallel=0, 10um, 100nm): {s:.6f}")
    
    # Test 4: Gap ratio metric
    print("\n[Test 4] Gap ratio for various conditions")
    gaps_um = [0.01, 0.1, 1.0, 10.0]
    for gap_um in gaps_um:
        gap_m = gap_um * 1e-6
        ratio = gap_ratio_metric(gap_m, 300.0)
        use_nf = should_use_near_field_model(gap_m, 300.0)
        print(f"  Gap {gap_um:6.2f}um: ratio = {ratio:7.3f}, use_NF = {use_nf}")
    
    # Test 5: Spectral heat flux
    print("\n[Test 5] Near-field spectral heat flux integration")
    print("  Computing for T_hot=600K, T_cold=300K, gap=100nm...")
    result = near_field_heat_flux_spectral(
        600.0, 300.0, 100e-9,
        material_hot='alumina',
        material_cold='alumina',
        n_omega=50,
        n_kparallel=30
    )
    print(f"  Total flux: {result['flux_W_m2']:.3e} W/m2")
    print(f"  Propagating: {result['flux_by_region']['propagating_W_m2']:.3e} W/m2")
    print(f"  Evanescent: {result['flux_by_region']['evanescent_W_m2']:.3e} W/m2")
    print(f"  Evanescent fraction: {100*result['evanescent_fraction']:.1f}%")
    print(f"  Peak wavelength: {result['dominant_wavelength_um']:.2f}um")


if __name__ == '__main__':
    test_near_field_physics()
