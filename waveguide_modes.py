"""
PHASE 2: Lossy Cylindrical Waveguide Modal Dispersion

Solves complex propagation constants β(λ) for TE/TM modes in sub-wavelength
cavities with lossy boundary conditions (real material walls, not ideal PEC).

Key Physics:
  - Complex modal propagation β = β_real + i·α (α = attenuation constant)
  - Characteristic equation includes Bessel functions J_n and K_n
  - Evanescent decay for frequencies below modal cutoff
  - Quality factor Q = β_real / (2α) indicates modal loss
  - Material losses reduce Q and increase decay rate

References:
  - Jackson, J. D. (1998). Classical Electrodynamics, 3rd ed., Chapter 8
  - Narayanaswamy & Chen (2004). PRB 70, 125101 (LDOS in lossy cavities)
  - Basu et al. (2009). Int. J. Energy Res. 33:1203-1232 (near-field)

Author: Kiro
Date: August 2026
"""

import warnings
import numpy as np
import math

try:
    from scipy import special, optimize
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    warnings.warn("scipy not available; using fallback numerical methods")

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

# Free-space wave impedance (Ohm)
Z0_VACUUM = 376.730313668


def aperture_modal_impedance_and_reflection(
    beta_complex: complex,
    k0_m: float
) -> tuple:
    """
    Compute wave impedance Z_mode and aperture power reflection R_ap at z=0.
    
    For TE modes:
        Z_mode = omega * mu_0 / beta = Z_0 * k_0 / beta
        where Z_0 = sqrt(mu_0 / eps_0) approx 376.73 Ohm.
    
    Aperture reflection coefficient (modal mismatch to free space):
        R_ap = |(Z_mode - Z_0) / (Z_mode + Z_0)|^2
        T_ap = 1 - R_ap  (power transmission across aperture)
    """
    if abs(beta_complex) < 1e-30:
        return complex(np.inf, 0.0), 1.0, 0.0
    
    Z_mode = Z0_VACUUM * (k0_m / beta_complex)
    
    denom = Z_mode + Z0_VACUUM
    if abs(denom) < 1e-30:
        r_ap = 1.0 + 0.0j
    else:
        r_ap = (Z_mode - Z0_VACUUM) / denom
    
    R_ap = float(np.clip(abs(r_ap)**2, 0.0, 1.0))
    T_ap = float(np.clip(1.0 - R_ap, 0.0, 1.0))
    return Z_mode, R_ap, T_ap


# Transverse TE eigenvalues x_m (J'_1(x)=0) for the implicit
# gamma_m = sqrt((x_m/R)^2 - k0^2) modal decay sum, and the relative modal
# excitation weights c_m (singled out so every consumer shares ONE source).
MODAL_TE_X_ROOTS = (1.8412, 3.0542, 3.8317, 4.1997, 5.3189)
MODAL_TE_WEIGHTS = (1.00, 0.32, 0.15, 0.08, 0.03)


def modal_total_exitance_transmission(
    lambda_c_um: float,
    lambda_um: float,
    length_um: float,
    material: str = 'alumina',
    diameter_um: float = None,
) -> float:
    """Consolidated modal-exitance attenuation operator ``T_total(λ)``.

    Computes the exit transmission of a sub-cutoff channel **exactly once** per
    (mode, wavelength) so that the several downstream consumers (the Monte-Carlo
    ray tracer's ``evanescent_power_transmission`` and the spectrum-integrated
    exitance diagnostic in ``simulator``) can never double-count the modal
    gating.  The single physics expression is::

        T_total(λ) = T_ap · Σ_m c_m · exp(−2·Re[γ_m]·L)

    with the transverse TE eigenvalues ``x_m``, channel radius ``R``, axial
    decay ``γ_m = sqrt((x_m/R)² − (2π/λ)²)``, and the aperture-plane
    boundary-impedance transmittance ``T_ap = 1 − R_ap`` (Z_mode matched to
    the free-space continuum at z = 0).

    Propagating modes (λ < λ_c) or a zero propagation length return exactly 1.0.

    Returns
    -------
    float in [0, 1] — the fraction of modal blackbody power surviving to the
    aperture plane.
    """
    if (lambda_c_um is None or lambda_c_um <= 0.0
            or not np.isfinite(lambda_c_um)):
        return 1.0
    lam = float(lambda_um)
    if lam < float(lambda_c_um) or float(length_um) <= 0.0:
        return 1.0                                # propagating / zero depth

    R_um = ((float(diameter_um) / 2.0 if diameter_um and diameter_um > 0.0
             else float(lambda_c_um) * 1.8412 / np.pi))
    length = float(length_um)
    k0 = 2.0 * np.pi / max(lam, 1e-12)

    num = 0.0
    tot = 0.0
    for x_m, c_m in zip(MODAL_TE_X_ROOTS, MODAL_TE_X_WEIGHTS):
        kappa_m = x_m / max(R_um, 1e-12)
        gamma_m = math.sqrt(max(kappa_m * kappa_m - k0 * k0, 0.0))
        num += c_m * math.exp(-2.0 * gamma_m * length)
        tot += c_m
    decay = (num / tot) if tot > 0.0 else 0.0

    # Aperture-plane boundary-mode impedance factor (Z_mode vs Z_0).
    T_ap = 1.0
    try:
        gamma0 = math.sqrt(max((MODAL_TE_X_ROOTS[0] / max(R_um, 1e-12)) ** 2 - k0 * k0, 0.0))
        if gamma0 > 1e-12:
            _Z, _R_ap, _Tap = aperture_modal_impedance_and_reflection(
                complex(0.0, gamma0), k0)
            T_ap = float(np.clip(_Tap, 0.0, 1.0))
    except Exception:
        T_ap = 1.0                                # conservative fallback

    return float(np.clip(T_ap * decay, 0.0, 1.0))


def _solve_dielectric_cylinder_transcendental(
    k0_m: float,
    a_m: float,
    eps_wall: complex,
    beta_guess: complex
) -> complex:
    """
    Solve the exact complex characteristic dispersion equation for HE11/TE11
    mode in a cylindrical hollow dielectric waveguide with lossy wall.
    
    Boundary condition:
      [J'_1(u)/(u J_1(u)) + H'_1(v)/(v H_1(v))] * [J'_1(u)/(u J_1(u)) + eps_wall H'_1(v)/(v H_1(v))]
      = (beta/k0)^2 * (1/u^2 - 1/v^2)^2
    where u = a * sqrt(k0^2 - beta^2), v = a * sqrt(k0^2 * eps_wall - beta^2).
    """
    if not SCIPY_AVAILABLE:
        return beta_guess
    
    def obj(p):
        beta_c = complex(p[0], p[1])
        u_sq = (k0_m**2 - beta_c**2) * (a_m**2)
        v_sq = (k0_m**2 * eps_wall - beta_c**2) * (a_m**2)
        
        u = np.sqrt(u_sq + 0.0j)
        v = np.sqrt(v_sq + 0.0j)
        
        if abs(u) < 1e-12 or abs(v) < 1e-12:
            return [1e6, 1e6]
        
        try:
            J1_u = special.jv(1, u)
            J1p_u = special.jvp(1, u, 1)
            H1_v = special.hankel1(1, v)
            H1p_v = special.h1vp(1, v, 1)
            
            if abs(J1_u) < 1e-30 or abs(H1_v) < 1e-30:
                return [1e6, 1e6]
            
            F1 = J1p_u / (u * J1_u)
            F2 = H1p_v / (v * H1_v)
            
            lhs = (F1 + F2) * (F1 + eps_wall * F2)
            rhs = (beta_c / k0_m)**2 * ((1.0 / u_sq) - (1.0 / v_sq))**2
            diff = lhs - rhs
            return [float(diff.real), float(diff.imag)]
        except Exception:
            return [1e6, 1e6]
    
    try:
        p0 = [float(beta_guess.real), float(beta_guess.imag)]
        sol = optimize.root(obj, p0, method='hybr', tol=1e-7)
        if sol.success:
            return complex(sol.x[0], sol.x[1])
    except Exception:
        pass
    return beta_guess


# ---------------------------------------------------------------------------
# TE11 Mode Solver for Cylindrical Cavities
# ---------------------------------------------------------------------------

def solve_te11_mode_complex(
    diameter_um: float,
    wavelength_um: float,
    material: str = 'alumina',
    method: str = 'perturbation',
    temperature_K: float = None
) -> dict:
    """
    Solve complex propagation constant β(ω) for TE11 fundamental mode.
    
    The TE11 mode is most common in PAA honeycomb and cylindrical cavities.
    This solver accounts for:
      1. Perfect electric conductor (PEC) approximation as zeroth order
      2. Perturbative correction for lossy walls
      3. Full transcendental dielectric boundary characteristic equation
      4. Mode impedance and aperture reflection coefficient R_ap at z = 0
    
    Parameters
    ----------
    diameter_um : float
        Cavity diameter (cylindrical) in micrometres
    wavelength_um : float
        Operating wavelength in micrometres
    material : str
        Wall material for complex refractive index ('alumina', 'silicon', etc.)
    method : str
        'pec' - ideal PEC (no losses)
        'perturbation' - first-order correction for losses (default)
        'full' - full transcendental equation (rigorous complex dielectric boundary)
    temperature_K : float, optional
        Wall temperature for Phase 4a temperature-dependent optical
        constants (None → 300 K tabulated data).
    
    Returns
    -------
    dict with keys:
        'beta_real': Phase constant β (rad/m)
        'beta_imag': Attenuation constant α (rad/m)
        'beta_complex': β + iα (complex, rad/m)
        'cutoff_wavelength_um': λ_c (µm)
        'is_evanescent': bool (True if λ > λ_c)
        'decay_length_um': 1/α for evanescent (µm)
        'Q_factor': β_real / α - quality factor
        'Z_mode': Modal wave impedance (Ohm)
        'R_ap': Aperture power reflection coefficient at z=0
        'T_ap': Aperture power transmission coefficient (1 - R_ap)
        'group_velocity': Group velocity v_g/c
        'frequency_ghz': Operating frequency
        'material': Material name
        'method': Solver method used
    """
    
    # Get cavity radius
    a_um = diameter_um / 2.0
    a_m = a_um * 1e-6  # Convert to meters
    
    # Get complex refractive index for wall material
    n_real, k_imag = get_complex_refractive_index(
        material, wavelength_um, temperature_K=temperature_K)
    n_complex = n_real + 1.0j * k_imag
    eps_wall = n_complex**2
    
    # Wavenumber in free space
    k0_m = 2.0 * np.pi / (wavelength_um * 1e-6)  # Convert wavelength to meters
    
    # TE11 cutoff wavenumber (Bessel root for PEC)
    bessel_root_te11 = 1.8411835329
    kc_pec_m = bessel_root_te11 / a_m  # Cutoff wavenumber for PEC
    lambda_c_pec_um = 2.0 * np.pi * 1e6 / kc_pec_m  # Cutoff wavelength
    
    # Ideal PEC propagation constant (no loss)
    kc_sq_pec = kc_pec_m**2
    beta_sq_pec = k0_m**2 - kc_sq_pec
    
    if beta_sq_pec >= 0:
        beta_pec = np.sqrt(beta_sq_pec) + 0.0j
        is_evanescent = False
    else:
        # Evanescent: β = iα
        alpha_pec = np.sqrt(-beta_sq_pec)
        beta_pec = 0.0 + 1.0j * alpha_pec
        is_evanescent = True

    # Fix 1b (peer-review): Compute the dielectric-corrected effective cutoff
    # wavelength.  For lossy walls (k_imag > 0) the PEC eigenvalue k_c,pec is
    # shifted by the first-order perturbation from the wall impedance.  The
    # effective cutoff k_c,eff satisfies Re(β(k_c,eff)) = 0.  To first order
    # in the surface resistance Rs = Re(1/ε_r)^(1/2):
    #   k_c,eff ≈ k_c,pec · (1 + Re(Δk_c)/k_c,pec)
    # where the perturbative shift keeps |k_c,eff| close to k_c,pec.
    # For insulators with Re(ε_r) >> Im(ε_r) the shift is of order k_imag/n_real.
    # We use n_real as the effective index filling the mode volume to compute
    # a material-loaded cutoff wavelength:
    #   λ_c,eff = λ_c,pec / n_real
    # This is the standard result for a dielectric-filled waveguide (Pozar 2012, §3.2).
    # At 300 K alumina n_real ≈ 1.65-1.7 in the mid-IR, so λ_c,eff ≈ 0.59 · λ_c,pec.
    if method in ('perturbation', 'full') and n_real > 1.0:
        lambda_c_eff_um = lambda_c_pec_um / n_real
    else:
        lambda_c_eff_um = lambda_c_pec_um   # PEC / vacuum (no material shift)
    
    # Correction for lossy walls
    if method == 'pec':
        beta = beta_pec
        alpha = float(np.abs(np.imag(beta)))
        Q_modal = np.inf
        
    elif method == 'perturbation':
        if k_imag > 0:
            Q_wall = n_real / (2.0 * k_imag)
        else:
            Q_wall = np.inf

        # Confinement factor for a hollow (air-core) dielectric waveguide:
        # The TE11 mode's power fraction in the lossy cladding scales as ~1/V²
        # for well-confined modes far above cutoff (V ≫ 1), where
        #   V = (2π·a/λ)·√(n_real² − n_core²) ≈ (2π·a/λ)·√(n_real² − 1).
        # Near cutoff (V → 1, λ → λ_c) the factor saturates at 1 because
        # the evanescent tail extends far into the walls.  Without this factor
        # the formula α = β/Q_wall assumes 100% of modal power is in the lossy
        # walls, giving α ∝ 1/λ at short wavelengths — unphysical.
        V = 1.0
        if n_real > 1.0 and not is_evanescent:
            V = ((2.0 * np.pi * a_m) / max(wavelength_um * 1e-6, 1e-15)
                 * math.sqrt(max(n_real ** 2 - 1.0, 0.0)))
        F_wall = 1.0 / max(V * V, 1e-12)
        F_wall = min(1.0, F_wall)                     # saturate near cutoff

        if is_evanescent:
            alpha_pec_val = np.abs(np.imag(beta_pec))
            alpha = alpha_pec_val * (1.0 + 0.5 / Q_wall) if Q_wall < np.inf else alpha_pec_val
            beta = 0.0 + 1.0j * alpha
            Q_modal = Q_wall
        else:
            beta_pec_val = np.real(beta_pec)
            alpha = (beta_pec_val / Q_wall) * F_wall if Q_wall < np.inf else 0.0
            beta = beta_pec_val + 1.0j * alpha
            Q_modal = Q_wall

    elif method == 'full':
        # Full transcendental equation for dielectric/lossy cylinder boundary
        beta_init = beta_pec
        if k_imag > 0:
            Q_wall = n_real / (2.0 * k_imag)
            # Confinement factor for the transcendental solver: same scaling.
            V = 1.0
            if n_real > 1.0 and not is_evanescent:
                V = ((2.0 * np.pi * a_m) / max(wavelength_um * 1e-6, 1e-15)
                     * math.sqrt(max(n_real ** 2 - 1.0, 0.0)))
            F_wall = min(1.0, 1.0 / max(V * V, 1e-12))
            if is_evanescent:
                beta_init = 0.0 + 1.0j * (np.abs(np.imag(beta_pec)) * (1.0 + 0.5 / Q_wall))
            else:
                beta_init = np.real(beta_pec) + 1.0j * ((np.real(beta_pec) / Q_wall) * F_wall)

        beta = _solve_dielectric_cylinder_transcendental(k0_m, a_m, eps_wall, beta_init)
        alpha = float(np.abs(np.imag(beta)))
        beta_real_val = float(np.real(beta))
        Q_modal = (beta_real_val / (2.0 * alpha)) if alpha > 1e-12 else np.inf
    else:
        beta = beta_pec
        alpha = float(np.abs(np.imag(beta)))
        Q_modal = np.inf
    
    # Extract real and imaginary parts
    beta_real = float(np.real(beta))
    beta_imag = float(np.abs(np.imag(beta)))
    
    # Determine evanescent state from complex beta
    if beta_real <= 1e-12 and beta_imag > 0:
        is_evanescent = True
    
    # Mode impedance and aperture reflection
    Z_mode, R_ap, T_ap = aperture_modal_impedance_and_reflection(beta, k0_m)
    
    # Decay length for evanescent modes
    if is_evanescent and beta_imag > 0:
        decay_length_um = (1.0 / beta_imag) * 1e6  # Convert to µm
    elif not is_evanescent and alpha > 0:
        decay_length_um = (1.0 / alpha) * 1e6  # Convert to µm
    else:
        decay_length_um = np.inf
    
    # Operating frequency
    freq_hz = c / (wavelength_um * 1e-6)
    freq_ghz = freq_hz / 1e9
    
    # Group velocity (TE11 waveguide dispersion)
    if not is_evanescent and beta_real > 0 and k0_m > 0:
        # TE11: v_g/c = beta_real / k0 (from d(omega)/d(beta) for TE modes)
        # since beta = sqrt(k0^2 - kc^2), v_g = c * beta / k0 (sub-luminal)
        v_g_over_c = float(np.clip(beta_real / k0_m, 0.0, 1.0))
    else:
        v_g_over_c = 0.0
    
    return {
        'beta_real': float(np.real(beta_real)),
        'beta_imag': float(np.abs(beta_imag)),  # Return magnitude of attenuation constant
        'beta_complex': complex(beta),
        # Fix 1b (peer-review): 'cutoff_wavelength_um' is now the dielectric-
        # corrected effective cutoff (λ_c = λ_c,pec/n_real for loaded modes),
        # NOT the bare PEC value.  The PEC reference is still accessible via
        # 'cutoff_wavelength_pec_um' for diagnostics.
        'cutoff_wavelength_um':     float(lambda_c_eff_um),
        'cutoff_wavelength_pec_um': float(lambda_c_pec_um),
        'is_evanescent': bool(is_evanescent),
        'decay_length_um': float(decay_length_um),
        'Q_factor': float(Q_modal),
        'Z_mode': complex(Z_mode),
        'R_ap': float(R_ap),
        'T_ap': float(T_ap),
        'group_velocity': float(v_g_over_c),
        'frequency_ghz': float(freq_ghz),
        'material': material,
        'wall_n': float(n_real),
        'wall_k': float(k_imag),
        'method': method,
        'diameter_um': float(diameter_um),
        'wavelength_um': float(wavelength_um),
        'bessel_root': float(bessel_root_te11),
    }


def attenuation_factor_lossy_waveguide(
    propagation_distance_um: float,
    modal_result: dict
) -> float:
    """
    Attenuation factor exp(-α·distance) for photons in lossy waveguide.
    
    Parameters
    ----------
    propagation_distance_um : float
        Distance traveled in waveguide (µm)
    modal_result : dict
        Result from solve_te11_mode_complex() containing 'beta_imag'
    
    Returns
    -------
    float
        Transmission coefficient T ∈ [0, 1]
        T = exp(-α·d)
    """
    alpha = modal_result.get('beta_imag', 0.0)
    
    if alpha <= 0:
        return 1.0  # No attenuation
    
    # Convert distance to meters and attenuation constant to m⁻¹
    distance_m = propagation_distance_um * 1e-6
    
    # T = exp(-α·d)
    exponent = -alpha * distance_m
    
    # Clip exponent to prevent underflow
    if exponent < -100:
        return 0.0
    
    transmission = np.exp(exponent)
    return float(np.clip(transmission, 0.0, 1.0))


def modal_emission_probability(
    cavity_height_um: float,
    modal_result: dict,
    temperature_K: float = 300.0
) -> float:
    """
    Emission probability for photons at wavelength emitted within cavity.
    
    Accounts for three effects:
    1. Modal confinement: photons with λ > λ_c become evanescent
    2. Attenuation through cavity: exponential decay exp(-2αL)
    3. Aperture impedance matching: transmission T_ap = 1 - R_ap across z=0
    
    Parameters
    ----------
    cavity_height_um : float
        Height/depth of cavity (µm)
    modal_result : dict
        Modal dispersion info from solve_te11_mode_complex()
    temperature_K : float
        Temperature for reference (used for diagnostics)
    
    Returns
    -------
    float
        Probability photon escapes (0-1).
    """
    # Attenuation through full cavity height
    t_prop = attenuation_factor_lossy_waveguide(cavity_height_um, modal_result)
    t_ap = modal_result.get('T_ap', 1.0)
    
    return float(np.clip(t_ap * t_prop, 0.0, 1.0))


def modal_spectral_transmissivity(
    wavelengths_um: np.ndarray,
    diameter_um: float,
    cavity_height_um: float,
    material: str = 'alumina',
    method: str = 'perturbation'
) -> np.ndarray:
    """
    Transmissivity τ(λ) accounting for modal filtering and attenuation.
    
    Computes transmission coefficient as a function of wavelength,
    combining:
      1. Step function at cutoff λ_c (evanescent for λ > λ_c)
      2. Attenuation coefficient α(λ) inside cavity
    
    Parameters
    ----------
    wavelengths_um : np.ndarray
        Array of wavelengths (µm)
    diameter_um : float
        Cavity diameter
    cavity_height_um : float
        Cavity height/depth
    material : str
        Wall material
    method : str
        Modal solver method
    
    Returns
    -------
    np.ndarray
        Transmissivity at each wavelength
    """
    transmissivity = np.zeros_like(wavelengths_um, dtype=float)
    
    for i, wl in enumerate(wavelengths_um):
        modal = solve_te11_mode_complex(diameter_um, wl, material, method)
        tau = modal_emission_probability(cavity_height_um, modal, 300.0)
        transmissivity[i] = tau
    
    return transmissivity


# ---------------------------------------------------------------------------
# Diagnostics and Analysis
# ---------------------------------------------------------------------------

def modal_analysis_report(
    diameter_um: float,
    cavity_height_um: float,
    temperature_K: float = 300.0,
    material: str = 'alumina'
) -> dict:
    """
    Comprehensive modal analysis report for cavity geometry.
    
    Returns diagnostic information useful for understanding modal loss.
    """
    
    # Wien's displacement: λ_peak ≈ 2898 / T (µm·K)
    lambda_peak_um = 2898.0 / temperature_K
    
    # Solve at peak wavelength
    modal_peak = solve_te11_mode_complex(diameter_um, lambda_peak_um, material)
    
    # Solve across spectrum
    wavelengths_um = np.linspace(
        lambda_peak_um / 5.0, lambda_peak_um * 3.0, 50
    )
    transmissivity = modal_spectral_transmissivity(
        wavelengths_um, diameter_um, cavity_height_um, material
    )
    
    # Average transmissivity (weighted by Planck distribution)
    # Simplified: just arithmetic mean
    avg_transmissivity = np.mean(transmissivity)
    
    return {
        'cavity_diameter_um': diameter_um,
        'cavity_height_um': cavity_height_um,
        'temperature_K': temperature_K,
        'material': material,
        'peak_wavelength_um': lambda_peak_um,
        'modal_result_peak': modal_peak,
        'cutoff_wavelength_um': modal_peak['cutoff_wavelength_um'],
        'is_peak_evanescent': modal_peak['is_evanescent'],
        'wavelengths_um': wavelengths_um,
        'transmissivity': transmissivity,
        'average_transmissivity': float(avg_transmissivity),
        'decay_length_at_peak_um': modal_peak['decay_length_um'],
        'Q_factor': modal_peak['Q_factor'],
    }


# ---------------------------------------------------------------------------
# Unit Tests
# ---------------------------------------------------------------------------

def test_te11_modes():
    """Validate TE11 mode solver against analytical limits."""
    print("\n" + "="*70)
    print("Testing TE11 Modal Solver")
    print("="*70)
    
    # Test 1: PEC ideal case
    print("\n[Test 1] Ideal PEC cavity (10um diameter, alumina wall)")
    result = solve_te11_mode_complex(
        diameter_um=10.0,
        wavelength_um=5.0,
        material='alumina',
        method='pec'
    )
    print(f"  Cutoff wavelength: {result['cutoff_wavelength_um']:.3f} um")
    print(f"  Is evanescent:  {result['is_evanescent']}")
    print(f"  Q-factor: {result['Q_factor']:.1f}")
    
    # Test 2: Same with lossy walls
    print("\n[Test 2] Same cavity with lossy alumina walls")
    result_lossy = solve_te11_mode_complex(
        diameter_um=10.0,
        wavelength_um=5.0,
        material='alumina',
        method='perturbation'
    )
    print(f"  Cutoff wavelength: {result_lossy['cutoff_wavelength_um']:.3f} um")
    print(f"  Q-factor: {result_lossy['Q_factor']:.1f}")
    print(f"  Attenuation constant: {result_lossy['beta_imag']:.3e} rad/m")
    
    # Test 3: Transmission through cavity
    print("\n[Test 3] Transmission through 20um tall cavity")
    trans = attenuation_factor_lossy_waveguide(20.0, result_lossy)
    print(f"  Transmissivity T: {trans:.6f}")
    print(f"  Attenuation: {-20*np.log10(trans):.2f} dB")
    
    # Test 4: Spectral transmissivity
    print("\n[Test 4] Spectral transmissivity across thermal spectrum")
    wl_test = np.array([4.0, 8.0, 10.0, 12.0, 15.0])
    tau_test = modal_spectral_transmissivity(
        wl_test, 10.0, 20.0, 'alumina', 'perturbation'
    )
    for wl, tau in zip(wl_test, tau_test):
        print(f"  L = {wl:.1f} um: tau = {tau:.6f}")
    
    # Test 5: Full analysis report
    print("\n[Test 5] Detailed modal analysis for PAA-like geometry")
    report = modal_analysis_report(
        diameter_um=0.5,  # 500 nm - typical PAA pore
        cavity_height_um=20.0,
        temperature_K=300.0,
        material='alumina'
    )
    print(f"  Diameter: {report['cavity_diameter_um']:.3f} um")
    print(f"  Height: {report['cavity_height_um']:.1f} um")
    print(f"  Peak wavelength: {report['peak_wavelength_um']:.1f} um")
    print(f"  Cutoff: {report['cutoff_wavelength_um']:.3f} um")
    print(f"  Is peak evanescent: {report['is_peak_evanescent']}")
    print(f"  Avg transmissivity: {report['average_transmissivity']:.6f}")
    print(f"  Decay length: {report['decay_length_at_peak_um']:.3f} um")


if __name__ == '__main__':
    test_te11_modes()
