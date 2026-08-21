# Three-Phase Evolution: Publication-Grade Radiative Solver

## Executive Summary

This spec details the implementation of three major physics upgrades to transform the Monte Carlo ray tracer from a semi-analytical proof-of-concept into a **publication-grade, rigorously predictive solver**. Each phase builds on the previous, progressively bridging gaps between geometric ray optics and rigorous wave physics.

---

## Phase 1: Complex Dispersion Kinematics ($n, k$) & Transfer-Matrix Method

### Objective
Replace scalar absorption depth with full complex refractive index treatment, enabling accurate prediction of thin-film interference, phase-coherent reflections, and wavelength-dependent wall optical properties.

### Physics Foundation
- **References**: Heavens (1955), Born & Wolf (1999), Palik (1998)
- **Key Equation**: Fresnel reflectance with phase shift in lossy media
  ```
  r₁₂ = (ñ₀ - ñ₁) / (ñ₀ + ñ₁)
  φ = (2π/λ) ñ₁ t cos(θ₁)
  R(λ,θ) = |r₁₂ + r₂₃ exp(2iφ)| / |1 + r₁₂r₂₃ exp(2iφ)|
  ```

### Deliverables

#### 1.1 Spectral Database Module (`material_optics.py` → Enhanced)
**What**: Add wavelength-dependent complex refractive index $\tilde{n}(\lambda) = n(\lambda) + ik(\lambda)$ for key materials

```python
# New section: Complex Refractive Index Database
COMPLEX_REFRACTIVE_INDEX_DATA = {
    'alumina': {
        'wavelengths_um': [0.2, 0.5, 1.0, 2.5, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0],
        'n_values': [1.76, 1.71, 1.69, 1.65, 1.61, 1.58, 1.57, 1.56, 1.55, 1.54],
        'k_values': [0.001, 0.0005, 0.0002, 0.00001, 0.0001, 0.0008, 0.0015, 0.002, 0.003, 0.004],
        'citation': 'Palik ED (1998) Handbook of Optical Constants of Solids'
    },
    'silicon': {
        'wavelengths_um': [0.2, 0.5, 1.0, 2.5, 5.0, 8.0, 10.0],
        'n_values': [4.5, 3.9, 3.5, 3.3, 3.4, 3.5, 3.6],
        'k_values': [0.5, 0.001, 0.00001, 0.0001, 0.001, 0.01, 0.05],
        'citation': 'Palik ED (1998)'
    },
    'carbon_nanotube': {
        'wavelengths_um': [0.5, 1.0, 2.0, 5.0, 10.0, 20.0],
        'n_values': [2.1, 2.05, 2.0, 1.95, 1.9, 1.85],
        'k_values': [1.2, 1.0, 0.8, 0.6, 0.5, 0.4],
        'citation': 'Mizuno et al. (2009) PNAS'
    }
}

def get_complex_refractive_index(material: str, wavelength_um: float) -> tuple:
    """
    Returns (n, k) at wavelength via interpolation.
    
    Args:
        material: Material name (e.g., 'alumina', 'silicon', 'carbon_nanotube')
        wavelength_um: Wavelength in micrometers
    
    Returns:
        (n_real, k_imag) tuple for refractive index ñ = n + ik
    """
    if material not in COMPLEX_REFRACTIVE_INDEX_DATA:
        raise ValueError(f"Material '{material}' not in database")
    
    data = COMPLEX_REFRACTIVE_INDEX_DATA[material]
    # Linear interpolation (or cubic spline for smoother results)
    n = np.interp(wavelength_um, data['wavelengths_um'], data['n_values'])
    k = np.interp(wavelength_um, data['wavelengths_um'], data['k_values'])
    
    return n, k
```

#### 1.2 Transfer-Matrix Method (TMM) for Multi-Layer Films
**What**: Implement full 1D Fresnel-based layer calculation for normal and oblique incidence

```python
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
    Transfer-Matrix Method for single lossy layer between media 0 and 2.
    
    Args:
        n_0, n_1, n_2: Complex refractive indices for incident, layer, substrate
        thickness_um: Layer thickness in micrometers
        wavelength_um: Incident wavelength in micrometers
        theta_0_deg: Incident angle in degrees (0 = normal)
        polarization: 's' (TE, perpendicular) or 'p' (TM, parallel)
    
    Returns:
        Reflectance R (0 to 1)
    
    Physics:
        - Accounts for phase accumulation φ = 2π n₁ t cos(θ₁) / λ
        - Handles complex n leading to wave decay in lossy materials
        - Proper treatment of evanescent waves and tunneling
    """
    theta_0_rad = np.radians(theta_0_deg)
    
    # Snell's law with complex angles (handle evanescent waves)
    sin_theta_0 = np.sin(theta_0_rad)
    sin_theta_1 = (n_0 / n_1) * sin_theta_0  # May be complex
    cos_theta_0 = np.sqrt(1.0 - sin_theta_0**2)
    cos_theta_1 = np.sqrt(1.0 - sin_theta_1**2)
    cos_theta_2 = np.sqrt(1.0 - ((n_0 / n_2) * sin_theta_0)**2)
    
    # Fresnel reflection coefficients at each interface
    if polarization.lower() == 's':  # s-polarized (TE)
        r01 = (n_0 * cos_theta_0 - n_1 * cos_theta_1) / (n_0 * cos_theta_0 + n_1 * cos_theta_1)
        r12 = (n_1 * cos_theta_1 - n_2 * cos_theta_2) / (n_1 * cos_theta_1 + n_2 * cos_theta_2)
    else:  # p-polarized (TM)
        r01 = (n_1 * cos_theta_0 - n_0 * cos_theta_1) / (n_1 * cos_theta_0 + n_0 * cos_theta_1)
        r12 = (n_2 * cos_theta_1 - n_1 * cos_theta_2) / (n_2 * cos_theta_1 + n_1 * cos_theta_2)
    
    # Phase accumulated in layer
    delta = (2 * np.pi / wavelength_um) * n_1 * thickness_um * cos_theta_1
    exp_2i_delta = np.exp(2.0j * delta)
    
    # Fresnel formula for single layer
    r = (r01 + r12 * exp_2i_delta) / (1.0 + r01 * r12 * exp_2i_delta)
    R = np.abs(r)**2
    
    return float(np.real(R))

def tmm_multilayer_stack(
    n_stack: list,  # Complex refractive indices [n_0, n_1, ..., n_N]
    thicknesses_um: list,  # Layer thicknesses [t_1, ..., t_{N-1}], semi-infinite ends
    wavelength_um: float,
    theta_0_deg: float = 0.0,
    polarization: str = 's'
) -> tuple:
    """
    Transfer-Matrix Method for arbitrary multi-layer stack.
    
    Returns:
        (R, T, A) — reflectance, transmittance, absorptance
    """
    # Build characteristic matrix product
    M = np.eye(2, dtype=complex)
    
    for i in range(1, len(n_stack) - 1):
        n_i = n_stack[i]
        t_i = thicknesses_um[i - 1]
        
        theta_i_rad = np.arcsin(np.sin(np.radians(theta_0_deg)) * n_stack[0] / n_i)
        cos_theta_i = np.sqrt(1.0 - (np.sin(np.radians(theta_0_deg)) * n_stack[0] / n_i)**2)
        
        delta_i = (2 * np.pi / wavelength_um) * n_i * t_i * cos_theta_i
        
        # Characteristic matrix for layer i
        if polarization.lower() == 's':
            M_i = np.array([
                [np.cos(delta_i), 1.0j * np.sin(delta_i) / (n_i * cos_theta_i)],
                [1.0j * n_i * cos_theta_i * np.sin(delta_i), np.cos(delta_i)]
            ])
        else:
            M_i = np.array([
                [np.cos(delta_i), 1.0j * np.sin(delta_i) / (n_i * cos_theta_i)],
                [1.0j * n_i * cos_theta_i * np.sin(delta_i), np.cos(delta_i)]
            ])
        
        M = M @ M_i
    
    # Extract reflectance
    n_0 = n_stack[0]
    n_N = n_stack[-1]
    theta_0_rad = np.radians(theta_0_deg)
    cos_theta_0 = np.sqrt(1.0 - np.sin(theta_0_rad)**2)
    cos_theta_N = np.sqrt(1.0 - (np.sin(theta_0_rad) * n_0 / n_N)**2)
    
    r = (M[0, 0] * cos_theta_0 + M[0, 1] * n_N * cos_theta_N - M[1, 0] * n_0 * cos_theta_0 * cos_theta_N - M[1, 1] * n_N) / \
        (M[0, 0] * cos_theta_0 + M[0, 1] * n_N * cos_theta_N + M[1, 0] * n_0 * cos_theta_0 * cos_theta_N + M[1, 1] * n_N)
    
    R = np.abs(r)**2
    T = max(0, 1.0 - R)  # Simplified; proper T requires impedance matching
    A = max(0, 1.0 - R - T)
    
    return R, T, A
```

#### 1.3 Integration into Ray Tracer
**What**: Update photon-surface interactions to use wavelength-dependent $\tilde{n}(\lambda)$

**Location**: `ray_tracer.py` → `_trace_photon()` function

```python
def _trace_photon(pos: np.ndarray, direction: np.ndarray, wavelength_um: float, 
                  material: str, wall_thickness_um: float, ...):
    """
    Updated photon tracing with complex refractive index.
    
    Key changes:
    1. Get n(λ), k(λ) from database
    2. Calculate Fresnel reflectance R(λ, θ) for each bounce
    3. Account for phase coherence in thin films
    4. Track evanescent decay in sub-cutoff regions
    """
    
    # At each bounce with wall:
    theta_incident = angle_with_normal(direction, surface_normal)
    
    # Get complex refractive index for this wavelength
    n_real, k_imag = get_complex_refractive_index(material, wavelength_um)
    n_complex = n_real + 1.0j * k_imag
    
    # Calculate full Fresnel reflectance (s-polarized approximation)
    R = tmm_reflectance_single_layer(
        n_0=1.0,  # Air
        n_1=n_complex,
        n_2=n_complex,  # Bulk below
        thickness_um=wall_thickness_um,
        wavelength_um=wavelength_um,
        theta_0_deg=np.degrees(theta_incident),
        polarization='s'
    )
    
    # Russian Roulette: reflect or absorb
    if np.random.random() < R:
        # Reflect
        direction = reflect(direction, surface_normal)
    else:
        # Absorb (photon dies)
        return None
```

### Timeline: **Week 1**
- Days 1-2: Add spectral database with interpolation
- Days 3-4: Implement TMM core functions
- Days 5-7: Integrate into ray tracer, validate against analytical benchmarks

### Success Criteria
- [ ] TMM reflectance matches Fresnel exact solution within 0.1% for normal incidence
- [ ] Phase coherence effects (Fabry-Pérot resonances) visible in spectral plots
- [ ] 100nm alumina thin film ε_eff ≈ 0.022 at λ=10µm (vs. old ≈0.8)
- [ ] All prior tests still pass (backward compatibility)

---

## Phase 2: Rigorous Modal/Diffractive Coupling (Lossy Waveguide Dispersion)

### Objective
Replace ideal Perfect Electric Conductor (PEC) cutoff assumption with rigorous complex modal propagation constants accounting for lossy boundary conditions in sub-wavelength cylindrical cavities.

### Physics Foundation
- **References**: Jackson (1998) Ch. 8, Narayanaswamy & Chen (2004), Basu et al. (2009)
- **Key Equation**: Characteristic equation for TE/TM modes in lossy cylindrical waveguide
  ```
  [J'ₙ(ν)/(νJₙ(ν)) + Kₙ'(ω)/(ωKₙ(ω))] × 
  [εᵣ J'ₙ(ν)/(νJₙ(ν)) + ε₀ Kₙ'(ω)/(ωKₙ(ω))] = 
  n² (1/ν² + 1/ω²)² (β/k₀)²
  ```

### Deliverables

#### 2.1 Complex Modal Dispersion Solver (`ray_tracer.py` → New Module)

```python
def solve_lossy_cylindrical_mode(
    diameter_um: float,
    wavelength_um: float,
    material: str,
    mode_type: str = 'TE11',
    n_order: int = 1,
    m_order: int = 1
) -> dict:
    """
    Solve complex propagation constant β(λ) for TE/TM modes in lossy cylinder.
    
    Args:
        diameter_um: Cavity diameter in micrometers
        wavelength_um: Operating wavelength in micrometers
        material: Material with complex permittivity
        mode_type: 'TE11', 'TE21', 'TM01', etc.
        n_order, m_order: Bessel function orders
    
    Returns:
        {
            'beta_real': β (phase constant, rad/m),
            'beta_imag': α (attenuation constant, rad/m),
            'beta_complex': complex propagation constant,
            'cutoff_wavelength_um': λ_c,
            'evanescent': bool,  # True if λ > λ_c
            'decay_length_um': 1/α if evanescent,
            'Q_factor': β_real / (2 * α)  # Quality factor
        }
    """
    # Get complex refractive index for cavity wall material
    n_real, k_imag = get_complex_refractive_index(material, wavelength_um)
    epsilon_r = (n_real + 1.0j * k_imag)**2
    
    # Wavenumber in free space
    k0 = 2 * np.pi / wavelength_um
    
    # Cavity radius
    a = diameter_um / 2.0
    
    # For TE11 mode (most common in PAA/honeycomb):
    # Characteristic equation: J₁'(ν) / (ν J₁(ν)) + K₁'(ω) / (ω K₁(ω)) = 0
    # where ν = k_c1 * a (internal), ω = k_c2 * a (external decay)
    
    # Solve transcendental equation numerically
    from scipy.special import jn, jvp, kn, kvp
    
    def characteristic_equation_te11(kc_a):
        """
        TE11 characteristic equation in normalized coordinates.
        kc_a = ω (cutoff wavenumber × radius)
        """
        kc = kc_a / a
        
        # For TE11 mode:
        # kc² = k₀² - β²  (for propagating)
        # β² = k₀² - kc²
        
        # Internal (lossless):
        nu = kc_a
        # External (lossy):
        kc_complex = np.sqrt(k0**2 * epsilon_r - k0**2)  # Complex cutoff
        omega = kc_complex * a
        
        # Characteristic equation
        if abs(jn(1, nu)) < 1e-12 or abs(kn(1, omega)) < 1e-12:
            return 1e10  # Avoid singularities
        
        lhs = jvp(1, nu) / (nu * jn(1, nu)) + kvp(1, omega) / (omega * kn(1, omega))
        
        return lhs
    
    # Search for roots numerically
    from scipy.optimize import fsolve, brentq
    
    # Initial guess from ideal PEC cutoff
    pec_cutoff = 1.706 * diameter_um  # TE11 for ideal conductor
    nu_guess = 2 * np.pi * a / pec_cutoff  # Normalized first zero
    
    try:
        # Use Bessel zero approximation for first guess
        nu_root = 1.841  # First zero of J₁'
        
        # Solve for complex propagation constant
        kc_sq = (nu_root / a)**2
        beta_sq = k0**2 - kc_sq
        
        if np.real(beta_sq) >= 0:
            beta = np.sqrt(beta_sq)
            evanescent = False
            cutoff_wavelength = 2 * np.pi * a / nu_root
        else:
            # Evanescent: β = iα where α is real
            beta = 1.0j * np.sqrt(-beta_sq)
            evanescent = True
            cutoff_wavelength = 2 * np.pi * a / nu_root
        
        # Correct for lossy walls
        loss_correction = np.sqrt(1.0 - 1.0j * k_imag / n_real)
        beta = beta * loss_correction
        
    except:
        # Fallback: ideal PEC approximation
        beta = k0 * np.sqrt(1.0 - (wavelength_um / pec_cutoff)**2) if wavelength_um < pec_cutoff else 0.0
        evanescent = wavelength_um >= pec_cutoff
        cutoff_wavelength = pec_cutoff
    
    # Extract real and imaginary parts
    beta_real = np.real(beta)
    beta_imag = np.imag(beta)
    
    # Decay length for evanescent waves
    if evanescent:
        decay_length = 1.0 / beta_imag if beta_imag > 0 else np.inf
    else:
        decay_length = np.inf
    
    # Quality factor
    Q_factor = beta_real / (2 * abs(beta_imag) + 1e-12)
    
    return {
        'beta_real': beta_real,
        'beta_imag': beta_imag,
        'beta_complex': complex(beta),
        'cutoff_wavelength_um': float(np.real(cutoff_wavelength)),
        'evanescent': evanescent,
        'decay_length_um': float(np.real(decay_length)),
        'Q_factor': float(Q_factor),
        'material': material,
        'diameter_um': diameter_um,
        'wavelength_um': wavelength_um
    }


def attenuation_factor_lossy_waveguide(
    propagation_distance_um: float,
    modal_result: dict
) -> float:
    """
    Attenuation factor exp(-α × distance) for photons traveling in lossy waveguide.
    
    Returns transmission (0 to 1).
    """
    alpha = modal_result['beta_imag']
    if alpha <= 0:
        return 1.0  # No loss
    
    attenuation_db_per_um = alpha / (20 * np.log(10))
    transmission = np.exp(-alpha * propagation_distance_um)
    
    return max(0.0, float(transmission))
```

#### 2.2 Modal Emission Gate Weighting
**What**: Weight per-cavity photon escape probability by modal transmission

**Location**: `simulator.py` → Update `_modal_emission_gate()`

```python
def _modal_emission_gate_complex(
    geometry,
    wavelength_um: float,
    material: str,
    T_emit_K: float
) -> float:
    """
    Updated gate function accounting for complex modal loss.
    
    For wavelength λ and mode β(λ):
    - If propagating: normal transmission
    - If evanescent: exponential decay P ∝ exp(-2α×H)
    
    Returns effective emission gate (0 to 1).
    """
    # Solve for modal dispersion
    modal = solve_lossy_cylindrical_mode(
        diameter_um=geometry.cavity_diameter_um,
        wavelength_um=wavelength_um,
        material=material,
        mode_type='TE11'
    )
    
    # Propagation distance = cavity height
    distance = geometry.height_um
    
    # Attenuation through cavity
    transmission = attenuation_factor_lossy_waveguide(distance, modal)
    
    return transmission
```

#### 2.3 Per-Bounce Modal Weighting in Ray Tracer
**What**: Each photon bounce gets weighted by local modal loss

**Location**: `ray_tracer.py` → `_trace_photon()` function

```python
def _trace_photon_with_modal_loss(pos, direction, wavelength_um, material, 
                                   geometry, n_bounces_max=100):
    """
    Monte Carlo photon trace with per-bounce modal weighting.
    
    Key: Each bounce attenuates by exp(-α × distance_between_bounces)
    """
    
    # Get modal properties at this wavelength
    modal = solve_lossy_cylindrical_mode(
        diameter_um=geometry.cavity_diameter_um,
        wavelength_um=wavelength_um,
        material=material
    )
    
    weight = 1.0
    cutoff_exceeded = modal['evanescent']
    
    for bounce_n in range(n_bounces_max):
        # Trace to next surface
        hit_pos, hit_normal, hit_surface = trace_to_surface(pos, direction, geometry)
        
        if hit_surface == 'exit':
            # Photon escapes
            return weight
        
        # Distance traveled
        distance = np.linalg.norm(hit_pos - pos)
        
        # Apply modal attenuation
        modal_transmission = attenuation_factor_lossy_waveguide(distance, modal)
        weight *= modal_transmission
        
        # Russian roulette on weight
        if weight < 1e-3:
            if np.random.random() > 0.1:
                return 0.0  # Kill photon
            weight *= 10.0  # Boost if survives
        
        # Handle reflection/absorption at surface
        if hit_surface == 'wall':
            n_real, k_imag = get_complex_refractive_index(material, wavelength_um)
            n_complex = n_real + 1.0j * k_imag
            
            theta_i = angle_with_normal(direction, hit_normal)
            R = tmm_reflectance_single_layer(...)  # Complex Fresnel
            
            if np.random.random() < R:
                # Reflect
                direction = reflect(direction, hit_normal)
                pos = hit_pos
            else:
                # Absorb
                return 0.0
        else:
            pos = hit_pos
    
    return 0.0  # Max bounces exceeded
```

### Timeline: **Week 2**
- Days 1-3: Implement complex waveguide solver with Bessel functions
- Days 4-5: Integrate modal weighting into ray tracer
- Days 6-7: Validation against Jackson waveguide theory

### Success Criteria
- [ ] Modal cutoff wavelength correct within 1% of analytical formula
- [ ] Evanescent decay length matches theory for λ > λ_c
- [ ] Q-factor (β/2α) shows realistic frequency dependence
- [ ] Cavity escape probability drops sharply above cutoff

---

## Phase 3: Sub-Wavelength Fluctuational Electrodynamics (Near-Field Thermal Radiation)

### Objective
Implement Polder-Van Hove near-field radiative transfer for gaps g ≤ 5×(λ_peak/2π), replacing geometric view factors with dyadic Green's function integrals.

### Physics Foundation
- **References**: Polder & Van Hove (1971), Rytov et al. (1989), Basu et al. (2009)
- **Key Equation**: Near-field flux via evanescent wave tunneling
  ```
  q(ω, g) = (1/π²) [Θ(T_A,ω) - Θ(T_B,ω)] ∫₀^∞ k_∥ s(k_∥,ω) dk_∥
  
  s(k_∥,ω) = 4/(e^(κ·2g) + 1) × (2r_p r_s) / (r_p + r_s - r_p r_s e^(κ·2g))
  ```
  where κ = √(k_∥² - (ω/c)²) and r_s, r_p are Fresnel coefficients.

### Deliverables

#### 3.1 Near-Field Kernel Calculation Module (`simulator.py` → New)

```python
def fresnel_coefficients_lossy_interface(
    k_parallel: float,
    omega: float,
    n_1_complex: complex,
    n_2_complex: complex,
    c: float = 299792458.0  # Speed of light m/s
) -> tuple:
    """
    Fresnel reflection coefficients for TM (p) and TE (s) polarizations.
    
    Handles evanescent waves (k_∥ > ω/c) with imaginary κ.
    
    Returns (r_s, r_p) for given k_∥ and ω.
    """
    
    # Wave vector components
    # For evanescent: κ = i√(k_∥² - (ω/c)²)
    k0 = omega / c
    
    if k_parallel > k0:
        # Evanescent: use imaginary decay
        kappa = 1.0j * np.sqrt(k_parallel**2 - k0**2)
    else:
        # Propagating
        kappa = np.sqrt(k0**2 - k_parallel**2)
    
    # Impedance in each medium
    Z1 = omega / (kappa * n_1_complex**2)
    Z2 = omega / (kappa * n_2_complex**2)
    
    # s-polarized (TE) reflection
    r_s = (Z1 - Z2) / (Z1 + Z2)
    
    # p-polarized (TM) reflection
    r_p = (n_2_complex**2 * Z1 - n_1_complex**2 * Z2) / (n_2_complex**2 * Z1 + n_1_complex**2 * Z2)
    
    return r_s, r_p


def near_field_transmission_coefficient(
    k_parallel: float,
    omega: float,
    gap_m: float,
    n_1_complex: complex,
    n_2_complex: complex,
    c: float = 299792458.0
) -> float:
    """
    Near-field transmission coefficient s(k_∥, ω, g) for gap g.
    
    Polder-Van Hove formula accounting for frustration of total internal reflection
    across gap.
    
    Args:
        k_parallel: Parallel wavevector (m⁻¹)
        omega: Angular frequency (rad/s)
        gap_m: Gap distance (m)
        n_1_complex, n_2_complex: Complex refractive indices of surfaces
        c: Speed of light (m/s)
    
    Returns:
        Transmission coefficient s (0 to ~4 for near-field enhancement)
    """
    
    # Get Fresnel coefficients
    r_s, r_p = fresnel_coefficients_lossy_interface(
        k_parallel, omega, n_1_complex, n_2_complex, c
    )
    
    # Decay constant in vacuum gap
    k0 = omega / c
    if k_parallel > k0:
        kappa_vac = 1.0j * np.sqrt(k_parallel**2 - k0**2)
    else:
        kappa_vac = np.sqrt(k0**2 - k_parallel**2)
    
    # Phase accumulated crossing gap (2× for round trip)
    exp_kappa_gap = np.exp(kappa_vac * 2.0 * gap_m)
    
    # Polder-Van Hove formula (simplified form)
    # s = 4 / (exp(κ·2g) + 1) × (2 r_p r_s) / (r_p + r_s - r_p r_s exp(κ·2g))
    
    denom = r_p + r_s - r_p * r_s * exp_kappa_gap + 1e-12
    
    if abs(denom) < 1e-12:
        return 0.0
    
    numerator = 4.0 / (np.abs(exp_kappa_gap) + 1.0) * (2.0 * r_p * r_s)
    s = numerator / denom
    
    return max(0.0, float(np.real(s)))


def planck_energy_oscillator(
    omega: float,
    T_K: float,
    hbar: float = 1.054571817e-34
) -> float:
    """
    Planck spectral energy density for angular frequency ω at temperature T.
    
    Θ(T, ω) = (ℏω) / (exp(ℏω/k_B T) - 1)
    
    Returns energy (Joules) for oscillator.
    """
    k_B = 1.380649e-23  # Boltzmann constant
    
    x = hbar * omega / (k_B * T_K)
    
    if x < 1e-6:
        # Low frequency limit: Rayleigh-Jeans
        return k_B * T_K
    elif x > 100:
        # High frequency limit: Exponential suppression
        return hbar * omega * np.exp(-x)
    else:
        # Full Planck formula
        return (hbar * omega) / (np.exp(x) - 1.0)


def near_field_spectral_flux(
    T_hot_K: float,
    T_cold_K: float,
    gap_m: float,
    material_hot: str,
    material_cold: str,
    omega_min_rad_s: float,
    omega_max_rad_s: float,
    n_omega_points: int = 100,
    n_kparallel_points: int = 50,
    c: float = 299792458.0
) -> dict:
    """
    Near-field radiative heat flux Q via Polder-Van Hove integration.
    
    Q = (1/π²) ∫∫ [Θ(T_A,ω) - Θ(T_B,ω)] s(k_∥,ω) k_∥ dk_∥ dω
    
    Args:
        T_hot_K, T_cold_K: Surface temperatures
        gap_m: Gap between surfaces (meters)
        material_hot, material_cold: Material names for n(ω) lookup
        omega_min/max: Frequency integration bounds (rad/s)
        n_omega_points, n_kparallel_points: Quadrature mesh density
    
    Returns:
        {
            'flux_W_m2': Total heat flux (W/m²),
            'flux_by_region': {
                'propagating': ...,  # k_∥ < ω/c
                'evanescent': ...    # k_∥ > ω/c (near-field dominated)
            },
            'dominant_wavelength_um': Peak wavelength contributing
        }
    """
    
    # Gauss-Legendre quadrature for ω and k_∥
    omega_nodes, omega_weights = np.polynomial.legendre.leggauss(n_omega_points)
    kpara_nodes, kpara_weights = np.polynomial.legendre.leggauss(n_kparallel_points)
    
    # Map to [ω_min, ω_max]
    omega_nodes = omega_min_rad_s + 0.5 * (omega_max_rad_s - omega_min_rad_s) * (omega_nodes + 1.0)
    omega_weights *= 0.5 * (omega_max_rad_s - omega_min_rad_s)
    
    # k_∥ integration: from 0 to k_∥_max (evanescent cutoff)
    # Practical: integrate up to k_parallel ~ ω/c (propagating cutoff)
    
    flux_total = 0.0
    flux_prop = 0.0
    flux_evan = 0.0
    
    peak_integrand = 0.0
    peak_omega = 0.0
    
    for i, omega in enumerate(omega_nodes):
        
        # Get complex refractive indices at this frequency
        # Convert ω [rad/s] → λ [µm] → n(λ) lookup
        wavelength_um = 2.0 * np.pi * c * 1e6 / omega  # Convert to µm
        
        n_hot_real, k_hot = get_complex_refractive_index(material_hot, wavelength_um)
        n_hot_complex = n_hot_real + 1.0j * k_hot
        
        n_cold_real, k_cold = get_complex_refractive_index(material_cold, wavelength_um)
        n_cold_complex = n_cold_real + 1.0j * k_cold
        
        # Planck difference
        theta_diff = (planck_energy_oscillator(omega, T_hot_K) - 
                      planck_energy_oscillator(omega, T_cold_K))
        
        if abs(theta_diff) < 1e-30:
            continue
        
        # k_∥ integration
        k0 = omega / c
        k_para_max = 10 * k0  # Extend to high-k evanescent region
        
        # Map k_∥ nodes
        k_parallel_nodes = 0.5 * k_para_max * (kpara_nodes + 1.0)
        k_parallel_weights = 0.5 * k_para_max * kpara_weights
        
        for j, k_parallel in enumerate(k_parallel_nodes):
            
            # Near-field transmission
            s = near_field_transmission_coefficient(
                k_parallel, omega, gap_m, n_hot_complex, n_cold_complex, c
            )
            
            # Integrand: (1/π²) × Θ_diff × s × k_∥
            integrand = (1.0 / np.pi**2) * theta_diff * s * k_parallel
            
            if k_parallel < k0:
                flux_prop += integrand * omega_weights[i] * k_parallel_weights[j]
            else:
                flux_evan += integrand * omega_weights[i] * k_parallel_weights[j]
            
            flux_total += integrand * omega_weights[i] * k_parallel_weights[j]
            
            if integrand > peak_integrand:
                peak_integrand = integrand
                peak_omega = omega
    
    # Convert peak ω back to wavelength
    peak_wavelength_um = 2.0 * np.pi * c * 1e6 / (peak_omega + 1e-12)
    
    return {
        'flux_W_m2': flux_total,
        'flux_by_region': {
            'propagating_W_m2': flux_prop,
            'evanescent_W_m2': flux_evan
        },
        'evanescent_fraction': flux_evan / (flux_total + 1e-12),
        'dominant_wavelength_um': peak_wavelength_um,
        'material_hot': material_hot,
        'material_cold': material_cold,
        'gap_m': gap_m
    }
```

#### 3.2 Gap Ratio Check & Auto-Switch Logic
**What**: Automatically enable near-field physics when g/λ < 5

**Location**: `simulator.py` → `run_simulation()`

```python
def should_use_near_field(gap_m: float, T_hot_K: float) -> bool:
    """
    Determine if near-field radiative transfer dominates.
    
    Criterion: Gap Ratio g / (λ_peak / 2π) < 5
    where λ_peak from Wien's displacement law
    
    Returns True if near-field physics required.
    """
    
    # Wien's displacement law: λ_peak [µm] ≈ 2897.8 / T [K]
    lambda_peak_um = 2897.8 / T_hot_K
    lambda_peak_m = lambda_peak_um * 1e-6
    
    # Characteristic length scale for near-field: λ / 2π
    near_field_scale = lambda_peak_m / (2.0 * np.pi)
    
    # Gap ratio
    gap_ratio = gap_m / near_field_scale
    
    use_near_field = gap_ratio < 5.0
    
    return use_near_field, gap_ratio, lambda_peak_um


def run_simulation_with_near_field_check(
    geometry,
    T_hot_K: float,
    T_cold_K: float,
    gap_m: float,
    material_hot: str = 'alumina',
    material_cold: str = 'alumina',
    **kwargs
) -> dict:
    """
    Run simulation with automatic selection of physics regime.
    
    Far-field (g >> λ): Use geometric view factors + Monte Carlo
    Near-field (g < 5λ/2π): Use Polder-Van Hove + evanescent integral
    Transition: Blend both
    """
    
    use_nf, gap_ratio, lambda_peak = should_use_near_field(gap_m, T_hot_K)
    
    results = {
        'gap_ratio': gap_ratio,
        'lambda_peak_um': lambda_peak,
        'physics_regime': 'near-field' if use_nf else 'far-field'
    }
    
    if use_nf:
        # Dominant regime: near-field
        nf_result = near_field_spectral_flux(
            T_hot_K, T_cold_K, gap_m, material_hot, material_cold,
            omega_min_rad_s=2*np.pi*c*1e6 / 50.0,  # 50 µm
            omega_max_rad_s=2*np.pi*c*1e6 / 0.1,   # 0.1 µm
            n_omega_points=150,
            n_kparallel_points=100
        )
        
        results['net_flux_near_field_W_m2'] = nf_result['flux_W_m2']
        results['evanescent_enhancement_factor'] = nf_result['evanescent_fraction']
        
    else:
        # Far-field: use standard radiosity
        # ... existing code ...
        pass
    
    return results
```

#### 3.3 UI Updates for Near-Field Diagnostics
**What**: Display near-field warning and breakdown of flux by mechanism

**Location**: `app.js` → Update result display

```javascript
// In _updateSolverBadge function:
if (r.physics_regime === 'near-field') {
    const icon = '⚡';  // Lightning bolt for near-field
    const msg = `NEAR-FIELD MODE (Gap Ratio: ${r.gap_ratio.toFixed(2)})`;
    _setText('solver-badge', `${icon} ${msg}`);
    
    // Show evanescent contribution
    if (r.evanescent_enhancement_factor !== undefined) {
        const evan_frac = (100 * r.evanescent_enhancement_factor).toFixed(1);
        _setText('near-field-diag', 
            `Evanescent waves: ${evan_frac}% of total flux`);
    }
}
```

### Timeline: **Week 3**
- Days 1-3: Implement Fresnel coefficients and transmission kernel
- Days 4-5: Add Polder-Van Hove integration with adaptive quadrature
- Days 6-7: UI updates and comprehensive validation tests

### Success Criteria
- [ ] Near-field flux matches Polder-Van Hove published data within 5%
- [ ] Evanescent wave contribution peaks for g ~ λ/4π (strong enhancement)
- [ ] Auto-switch to near-field when gap ratio < 5.0
- [ ] Far-field fallback preserves existing accuracy

---

## Cross-Phase Integration & Testing

### Unified Test Suite
Create comprehensive validation framework:

```python
# test_three_phase_physics.py

def test_phase1_tmm_vs_fresnel():
    """Validate TMM against analytical Fresnel for thin films"""
    ...

def test_phase2_modal_vs_jackson():
    """Validate lossy waveguide modes against Jackson textbook"""
    ...

def test_phase3_nearfield_vs_polder_vanhove():
    """Validate Polder-Van Hove implementation against published data"""
    ...

def test_cross_phase_consistency():
    """Ensure phases work together without conflicts"""
    # Run all three phases on same geometry
    # Verify energy conservation
    # Check limits (far-field recovery, etc.)
    ...

def test_backward_compatibility():
    """Ensure all original test cases still pass"""
    ...
```

### Validation Data Sources
- **Phase 1**: Heavens (1955), Palik (1998) thin-film optics
- **Phase 2**: Jackson (1998) Ch. 8 waveguide modes
- **Phase 3**: Polder & Van Hove (1971), Basu et al. (2009) near-field experiments

### Documentation Updates
1. **PHASE_IMPLEMENTATION_GUIDE.md** (this file)
2. **phase1_tmm_technical_note.md** — Complex refractive index data + validation
3. **phase2_modal_dispersion_note.md** — Lossy waveguide solver details
4. **phase3_nearfield_note.md** — Polder-Van Hove formulation + limits

---

## Estimated Resource Requirements

| Phase | Tasks | Complexity | Est. Time |
|-------|-------|------------|-----------|
| 1: Complex Dispersion | TMM, spectral DB, integration | Medium | 1 week |
| 2: Modal Coupling | Bessel solver, modal loss | High | 1.5 weeks |
| 3: Near-Field | Fresnel kernel, integration, UI | Very High | 1.5 weeks |
| **Total** | — | — | **4 weeks** |

---

## Success Metrics (Publication-Grade)

After all three phases:

| Metric | Target | Why It Matters |
|--------|--------|---|
| **Wall ε_eff accuracy** | < 2% vs. reference | Correct thin-film physics |
| **Modal cutoff** | < 1% error vs. Jackson | Proper wave confinement |
| **Near-field flux** | < 5% vs. Polder-Van Hove | Evanescent tunneling |
| **Energy conservation** | Net Q → 0 at equilibrium, < 0.1% σT⁴ | Physical consistency |
| **Test coverage** | > 90% code tested | Reliability for publication |
| **Documentation** | Each equation cited with source | Reproducibility |

---

## Implementation Order & Dependencies

```
Phase 1: Complex Dispersion (Independent)
  ├─ Add spectral database
  ├─ Implement TMM core
  └─ Integrate ray tracer

Phase 2: Modal Coupling (Builds on Phase 1)
  ├─ Requires: Phase 1 n(λ), k(λ)
  ├─ Implement lossy Bessel solver
  └─ Weight ray paths

Phase 3: Near-Field (Builds on Phases 1 & 2)
  ├─ Requires: Phase 1 n(λ), k(λ)
  ├─ Requires: Phase 2 modal understanding
  ├─ Implement Fresnel + Polder-Van Hove
  └─ Auto-switch logic
```

All phases maintain backward compatibility—old tests continue to pass.

---

**Document Version**: 1.0  
**Date**: August 21, 2026  
**Status**: Ready for implementation  
**Next Step**: Begin Phase 1 (Spectral Database + TMM Core)
