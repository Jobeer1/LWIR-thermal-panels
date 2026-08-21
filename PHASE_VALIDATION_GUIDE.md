# Phase Implementation Validation Guide

## Overview

This guide validates the three-phase physics upgrade to the Monte Carlo ray tracer. Each phase includes test cases, success criteria, and expected outputs.

---

## Phase 1: Complex Dispersion & Transfer-Matrix Method

### Files
- **`material_optics.py`** (enhanced with Phase 1 additions)
- New databases: `COMPLEX_REFRACTIVE_INDEX_DATA`, `MATERIAL_COMPLEX_INDEX`
- New functions: `get_complex_refractive_index()`, `tmm_reflectance_single_layer()`

### Test Cases

#### Test 1.1: Complex Refractive Index Database
**Objective**: Verify complex index interpolation matches literature values

```python
# Expected output at key wavelengths
from material_optics import get_complex_refractive_index

# Alumina at 10 µm (thermal peak at 300K)
n, k = get_complex_refractive_index('alumina', 10.0)
print(f"n = {n:.4f} (expect ~1.57)")
print(f"k = {k:.6f} (expect ~0.002)")

# Silicon in visible
n, k = get_complex_refractive_index('silicon', 0.5)
print(f"n = {n:.2f} (expect ~3.9)")
print(f"k = {k:.4f} (expect ~0.001)")

# Carbon nanotube in IR
n, k = get_complex_refractive_index('carbon_nanotube', 10.0)
print(f"n = {n:.2f} (expect ~1.9)")
print(f"k = {k:.2f} (expect ~0.5)")
```

**Success Criteria**:
- ✓ All interpolated values within ±5% of literature (Palik 1998)
- ✓ No warnings for wavelengths within database range
- ✓ Graceful fallback for unknown materials

---

#### Test 1.2: Fresnel Reflectance (Single Layer)
**Objective**: Verify TMM against analytical Fresnel formula for normal incidence

```python
from material_optics import tmm_reflectance_single_layer
import numpy as np

# Test case: air-alumina interface at normal incidence
n_air = 1.0 + 0.0j
n_alumina = 1.57 + 0.002j
thickness = 0.1  # 100 nm
wavelength = 10.0  # µm

R = tmm_reflectance_single_layer(n_air, n_alumina, n_alumina, thickness, wavelength)

# Analytical Fresnel: R = |(n_1 - n_2) / (n_1 + n_2)|²
r_analytical = (1.0 - n_alumina) / (1.0 + n_alumina)
R_analytical = abs(r_analytical)**2

print(f"TMM result: R = {R:.6f}")
print(f"Fresnel analytical: R = {R_analytical:.6f}")
print(f"Error: {abs(R - R_analytical):.2e}")
```

**Success Criteria**:
- ✓ TMM ≤ 0.1% error vs. Fresnel for normal incidence
- ✓ Reflectance bounded [0, 1]
- ✓ Complex indices produce non-zero k dependence

---

#### Test 1.3: Thin-Film Correction Factor
**Objective**: Verify 100nm alumina wall shows ~35× reduction in emissivity

```python
from material_optics import effective_emissivity_thin_film, optical_thickness_analysis

# 100 nm alumina wall at 300K
thickness = 0.1  # µm
temperature = 300.0

analysis = optical_thickness_analysis(thickness, 'alumina', temperature)

print(f"Wall thickness: {analysis['wall_thickness_um']}µm")
print(f"Material: {analysis['material']}")
print(f"Peak wavelength: {analysis['peak_wavelength_um']:.2f}µm")
print(f"Absorption depth: {analysis['absorption_depth_peak_um']:.2f}µm")
print(f"Optical thickness t/δ: {analysis['optical_thickness_t/delta']:.4f}")
print(f"Bulk emissivity: {analysis['bulk_emissivity']:.3f}")
print(f"Thin-film ε_eff: {analysis['effective_emissivity_peak']:.4f}")
print(f"Correction factor: {analysis['thin_film_correction_factor']:.2f}")
```

**Expected Output**:
```
Wall thickness: 0.1µm
Material: alumina
Peak wavelength: 9.66µm
Absorption depth: 3.61µm
Optical thickness t/δ: 0.0277
Bulk emissivity: 0.800
Thin-film ε_eff: 0.0220
Correction factor: 0.0275 (i.e., ~1/36.4 reduction)
```

**Success Criteria**:
- ✓ Correction factor (ε_eff / ε_bulk) ≈ 1/35 ± 5%
- ✓ Optically thin detection working (t/δ < 1)
- ✓ Planck-weighted emissivity reasonable

---

#### Test 1.4: Phase Coherence Effects
**Objective**: Verify Fabry-Pérot interference visible in spectral reflectance

```python
from material_optics import tmm_reflectance_single_layer
import matplotlib.pyplot as plt

# Scan wavelength for fixed thickness (shows Fabry-Pérot fringes)
wavelengths = np.linspace(5, 20, 200)  # 5-20 µm
thicknesses = [0.1, 1.0, 5.0]  # 0.1, 1, 5 µm films

for t in thicknesses:
    reflectances = []
    for wl in wavelengths:
        R = tmm_reflectance_single_layer(
            1.0+0.0j, 1.57+0.002j, 1.57+0.002j, t, wl
        )
        reflectances.append(R)
    
    plt.plot(wavelengths, reflectances, label=f't={t:.1f}µm')

plt.xlabel('Wavelength (µm)')
plt.ylabel('Reflectance')
plt.legend()
plt.title('Fabry-Pérot Interference in Thin Films')
plt.show()
```

**Success Criteria**:
- ✓ Thicker films show more oscillations (higher finesse)
- ✓ Oscillations damped by absorption (k term)
- ✓ Envelopes match Beer-Lambert prediction

---

## Phase 2: Lossy Waveguide Modal Dispersion

### Files
- **`waveguide_modes.py`** (new module)
- New functions: `solve_te11_mode_complex()`, `attenuation_factor_lossy_waveguide()`, `modal_analysis_report()`

### Test Cases

#### Test 2.1: Cutoff Wavelength (PEC Limit)
**Objective**: Verify λ_c matches Jackson textbook for ideal conductor

```python
from waveguide_modes import solve_te11_mode_complex

# Test at wavelengths above and below cutoff
diameter = 10.0  # µm
wavelengths = [5.0, 8.5, 10.0, 15.0, 20.0]

for wl in wavelengths:
    result = solve_te11_mode_complex(
        diameter, wl, material='alumina', method='pec'
    )
    
    print(f"λ = {wl:.1f}µm:")
    print(f"  λ_c = {result['cutoff_wavelength_um']:.2f}µm")
    print(f"  Evanescent: {result['is_evanescent']}")
    print(f"  Q-factor: {result['Q_factor']:.1f}")
```

**Expected Output** (for 10µm diameter):
```
λ = 5.0µm: λ_c = 17.06µm, Evanescent: False, Q = inf
λ = 8.5µm: λ_c = 17.06µm, Evanescent: False, Q = inf
λ = 10.0µm: λ_c = 17.06µm, Evanescent: False, Q = inf
λ = 15.0µm: λ_c = 17.06µm, Evanescent: False, Q = inf
λ = 20.0µm: λ_c = 17.06µm, Evanescent: True, Q = inf
```

**Analytical Formula**: λ_c = 1.706 × diameter = 1.706 × 10 = 17.06 µm

**Success Criteria**:
- ✓ λ_c accurate within 0.1%
- ✓ Evanescent flag correct (λ > λ_c)
- ✓ PEC method gives infinite Q

---

#### Test 2.2: Lossy Wall Effects on Q-Factor
**Objective**: Verify material losses reduce Q-factor

```python
from waveguide_modes import solve_te11_mode_complex

# Same cavity, different materials
diameter = 10.0  # µm
wavelength = 10.0  # µm (propagating)

materials = ['alumina', 'silicon', 'carbon_nanotube']

for mat in materials:
    result = solve_te11_mode_complex(
        diameter, wavelength, material=mat, method='perturbation'
    )
    
    print(f"{mat}:")
    print(f"  Wall n: {result['wall_n']:.3f}")
    print(f"  Wall k: {result['wall_k']:.4f}")
    print(f"  Q-factor: {result['Q_factor']:.1f}")
    print(f"  Attenuation: {result['beta_imag']:.3e} rad/m")
```

**Expected Trend**:
- Alumina (low k): Q ~ 100-1000
- Silicon (medium k): Q ~ 10-100
- CNT (high k): Q ~ 1-10

**Success Criteria**:
- ✓ Q-factor decreases with increasing k
- ✓ Q ≈ n / (2k) relationship holds
- ✓ Attenuation constant > 0 for lossy materials

---

#### Test 2.3: Evanescent Decay Length
**Objective**: Verify sub-cutoff modes decay correctly

```python
from waveguide_modes import solve_te11_mode_complex, modal_analysis_report

# Small cavity (below cutoff at 10µm)
diameter = 5.0  # µm (λ_c = 8.53µm)
height = 20.0  # µm
wavelength = 15.0  # µm (sub-cutoff)

result = solve_te11_mode_complex(diameter, wavelength, 'alumina', 'perturbation')

print(f"Cavity diameter: {diameter}µm")
print(f"Wavelength: {wavelength}µm (sub-cutoff)")
print(f"Cutoff: {result['cutoff_wavelength_um']:.2f}µm")
print(f"Evanescent: {result['is_evanescent']}")
print(f"Decay length: {result['decay_length_um']:.3f}µm")
print(f"Height: {height}µm")
print(f"Decay ratio H/δ_ev: {height / result['decay_length_um']:.2f}")

# Expected: for sub-cutoff, photons decay rapidly
expected_survival = np.exp(-height / result['decay_length_um'])
print(f"Expected photon survival: {expected_survival:.3e}")
```

**Success Criteria**:
- ✓ Decay length calculated for evanescent modes
- ✓ δ_ev ~ λ/(2π) × √[(λ/λ_c)² - 1]⁻¹
- ✓ Photon survival drops exponentially with depth

---

#### Test 2.4: Modal Analysis Report
**Objective**: Generate comprehensive diagnostic for PAA-like geometry

```python
from waveguide_modes import modal_analysis_report

# Typical PAA geometry
report = modal_analysis_report(
    diameter_um=0.5,  # 500 nm
    cavity_height_um=20.0,  # 20 µm
    temperature_K=300.0,
    material='alumina'
)

print("Modal Analysis Report")
print("=" * 50)
for key, value in report.items():
    if key not in ['wavelengths_um', 'transmissivity']:
        print(f"{key}: {value}")

print(f"\nSpectral Analysis:")
print(f"  n_wavelengths: {len(report['wavelengths_um'])}")
print(f"  Avg transmissivity: {report['average_transmissivity']:.6f}")
```

**Expected Output**:
```
Modal Analysis Report
==================================================
cavity_diameter_um: 0.5
cavity_height_um: 20.0
temperature_K: 300.0
material: alumina
peak_wavelength_um: 9.66
cutoff_wavelength_um: 0.853
is_peak_evanescent: False
average_transmissivity: 0.999835
decay_length_at_peak_um: inf
Q_factor: 76.923
```

**Success Criteria**:
- ✓ Transmission across spectrum computed correctly
- ✓ Average transmissivity reasonable (not 0 or 1)
- ✓ Q-factor consistent with prior tests

---

## Phase 3: Near-Field Radiative Heat Transfer

### Files
- **`near_field_radiative_heat.py`** (new module)
- New functions: `fresnel_coefficients_interface()`, `near_field_transmission_coefficient()`, `near_field_heat_flux_spectral()`

### Test Cases

#### Test 3.1: Gap Ratio Metric
**Objective**: Verify dimensionless gap metric for regime identification

```python
from near_field_radiative_heat import gap_ratio_metric, should_use_near_field_model

# Test at different gap scales
gaps_um = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
temp = 300.0

print("Gap Ratio Analysis (T=300K, λ_peak≈9.66µm)")
print("=" * 60)
print(f"{'Gap (µm)':>10} | {'Gap Ratio':>12} | {'Near-Field?':>15}")
print("-" * 60)

for gap_um in gaps_um:
    gap_m = gap_um * 1e-6
    ratio = gap_ratio_metric(gap_m, temp)
    use_nf = should_use_near_field_model(gap_m, temp, threshold=5.0)
    print(f"{gap_um:10.3f} | {ratio:12.3f} | {str(use_nf):>15}")
```

**Expected Output**:
```
Gap Ratio Analysis (T=300K, λ_peak≈9.66µm)
============================================================
    Gap (µm) |     Gap Ratio |     Near-Field?
------------------------------------------------------------
       0.001 |         0.041 |            True   ← Strong NF
       0.010 |         0.410 |            True
       0.100 |         4.100 |            True
       1.000 |        41.000 |           False   ← Far-field
      10.000 |       410.000 |           False
     100.000 |      4100.000 |           False
```

**Success Criteria**:
- ✓ Threshold at gap ratio ≈ 5.0 (tunable)
- ✓ Sub-micron gaps always near-field
- ✓ Hundred-micron gaps always far-field

---

#### Test 3.2: Fresnel Coefficients
**Objective**: Verify complex Fresnel formulation handles evanescent waves

```python
from near_field_radiative_heat import fresnel_coefficients_interface
import numpy as np

# Test at normal incidence (k_∥ = 0)
omega = 2*np.pi * c / 10e-6  # 10 µm wavelength
n_test = 1.5 + 0.1j

r_s, r_p = fresnel_coefficients_interface(0.0, omega, n_test, 1.0+0.0j)

print(f"Normal incidence (10µm):")
print(f"  r_s = {r_s:.4f}")
print(f"  r_p = {r_p:.4f}")
print(f"  |r_s|² = {abs(r_s)**2:.6f}")
print(f"  |r_p|² = {abs(r_p)**2:.6f}")

# Analytical Fresnel at normal incidence
r_analytical = (1.0 - n_test) / (1.0 + n_test)
print(f"  Analytical: |r|² = {abs(r_analytical)**2:.6f}")
```

**Success Criteria**:
- ✓ Normal incidence matches Fresnel formula
- ✓ Oblique incidence includes evanescent decay
- ✓ Reflectances bounded [0, 1] for passive media

---

#### Test 3.3: Near-Field Transmission Coefficient
**Objective**: Verify transmission shows resonance near cutoff

```python
from near_field_radiative_heat import near_field_transmission_coefficient
import numpy as np

# Scan k_∥ at fixed frequency
omega = 2*np.pi * c / 10e-6
gap = 100e-9  # 100 nm
n = 1.5 + 0.1j

k0 = omega / c
k_parallel_range = np.linspace(0, 5*k0, 100)

transmission = []
for k_par in k_parallel_range:
    s = near_field_transmission_coefficient(k_par, omega, gap, n, n)
    transmission.append(s)

# Plot or analyze
prop_cutoff_idx = np.argmin(np.abs(k_parallel_range - k0))
print(f"Transmission at k_∥ = 0: {transmission[0]:.6f}")
print(f"Transmission at k_∥ = k_0 (cutoff): {transmission[prop_cutoff_idx]:.6f}")
print(f"Transmission at k_∥ = 5×k_0 (strong evanescent): {transmission[-1]:.6f}")
```

**Expected Trend**:
- Propagating region (k_∥ < k_0): Moderate transmission
- Near cutoff (k_∥ ≈ k_0): Peak transmission
- Evanescent (k_∥ >> k_0): Exponential decay

**Success Criteria**:
- ✓ Transmission maximum near k_∥ = k_0
- ✓ Smooth decay for large k_∥
- ✓ Physical bounds maintained

---

#### Test 3.4: Spectral Heat Flux Integration
**Objective**: Compute near-field heat transfer for realistic geometry

```python
from near_field_radiative_heat import near_field_heat_flux_spectral

# High-temperature difference, small gap
result = near_field_heat_flux_spectral(
    temperature_hot_K=600.0,
    temperature_cold_K=300.0,
    gap_m=100e-9,  # 100 nm
    material_hot='alumina',
    material_cold='alumina',
    n_omega=80,  # Coarse for speed
    n_kparallel=50
)

print("Near-Field Heat Flux Integration")
print("=" * 50)
print(f"T_hot:  {result['integration_info']['gap_m']*1e9:.1f} nm gap")
print(f"Total flux: {result['flux_W_m2']:.3e} W/m²")
print(f"  Propagating: {result['flux_by_region']['propagating_W_m2']:.3e} W/m²")
print(f"  Evanescent: {result['flux_by_region']['evanescent_W_m2']:.3e} W/m²")
print(f"  Evanescent fraction: {100*result['evanescent_fraction']:.1f}%")
print(f"Peak wavelength: {result['dominant_wavelength_um']:.2f}µm")

# Compare to blackbody Stefan-Boltzmann estimate
sigma = 5.67e-8  # W/(m²·K⁴)
Stefan_Boltzmann = sigma * (600**4 - 300**4)
print(f"\nStefan-Boltzmann (far-field): {Stefan_Boltzmann:.3e} W/m²")
print(f"Near-field / far-field ratio: {result['flux_W_m2'] / Stefan_Boltzmann:.2f}×")
```

**Expected Output** (for 100nm gap, ΔT=300K):
```
Near-Field Heat Flux Integration
==================================================
T_hot:  100.0 nm gap
Total flux: ~1e5 W/m² (orders of magnitude above SB)
  Propagating: ~1e4 W/m²
  Evanescent: ~1e5 W/m²
  Evanescent fraction: 90-95%
Peak wavelength: 9.66µm

Stefan-Boltzmann (far-field): ~6000 W/m²
Near-field / far-field ratio: ~15-50×
```

**Success Criteria**:
- ✓ Evanescent fraction significant (>50% for g < λ/4π)
- ✓ Total flux orders of magnitude above blackbody (expected)
- ✓ Integration stability (no NaNs or infinities)

---

## Integration Tests

### Test I.1: Ray Tracer with Phase 1 TMM
**Objective**: Verify ray tracer accepts complex Fresnel reflectances

```python
# In ray_tracer.py, update _trace_photon() to use:
# R_tmm = tmm_reflectance_single_layer(
#     n_0=1.0+0.0j,
#     n_1=n_material + 1.0j*k_material,
#     n_2=n_material + 1.0j*k_material,
#     thickness_um=wall_thickness,
#     wavelength_um=photon_wavelength,
#     theta_0_deg=incident_angle
# )
```

**Success Criteria**:
- ✓ Photon absorption probability matches R_tmm
- ✓ All existing Monte Carlo tests still pass
- ✓ Runtime unchanged (within 10%)

---

### Test I.2: Simulator with Phase 2 Modal Weighting
**Objective**: Verify cavity escape probability accounts for modal loss

```python
# In simulator.py, integrate:
# modal = solve_te11_mode_complex(diameter, wavelength, material)
# P_escape *= attenuation_factor_lossy_waveguide(height, modal)
```

**Success Criteria**:
- ✓ Sub-cutoff photons heavily attenuated
- ✓ Above-cutoff photons mostly escape
- ✓ Average emissivity reduction visible (~10-50%)

---

### Test I.3: Simulator with Phase 3 Near-Field Logic
**Objective**: Verify auto-switch to near-field model for small gaps

```python
# In simulator.py, before radiosity calculation:
# gap_ratio = gap_ratio_metric(gap, T_hot)
# if gap_ratio < 5.0:
#     Q_net = near_field_heat_flux_spectral(...)
# else:
#     Q_net = radiosity_far_field(...)
```

**Success Criteria**:
- ✓ Gap ratio correctly computed
- ✓ Near-field model activated for g < 5×λ/(2π)
- ✓ Results physically plausible (no divergence)

---

## Backward Compatibility Verification

### Test B.1: Existing Simulations
**Objective**: Ensure all prior test cases still pass

```bash
python test.py
python test_thin_film.py
python test_wave_benchmarks.py
```

**Success Criteria**:
- ✓ All tests pass
- ✓ Results match within ±1% (statistical variation)
- ✓ No deprecation warnings

---

## Performance Benchmarks

| Operation | Time | Requirement |
|-----------|------|-------------|
| Complex index lookup (1000×) | <10ms | <1s |
| TMM calculation (1000×) | <100ms | <1s |
| Modal solver (100×) | <500ms | <1s |
| Near-field integral (1×) | <5s | <30s |
| Full simulation (20k photons) | <30s | <60s |

---

## Expected Accuracy After All Phases

| Metric | Before Phases | After Phases | Target |
|--------|---------------|--------------|--------|
| Wall ε_eff error | 35× overestimation | <5% | <2% |
| Modal cutoff error | Ideal PEC only | <1% vs. Jackson | <0.5% |
| Near-field flux error | N/A | <10% vs. Polder-VH | <5% |
| Net flux at equilibrium | ±0.5-2% σT⁴ | <0.1% σT⁴ | <0.1% |

---

## Validation Output Template

When all tests pass, generate summary:

```
=============================================================
THREE-PHASE PHYSICS UPGRADE VALIDATION REPORT
=============================================================

Date: [date]
Tester: [name]
Simulator Version: [version]

PHASE 1: Complex Dispersion & TMM
  [✓] Complex refractive index database
  [✓] Fresnel reflectance calculations
  [✓] Thin-film correction factor
  [✓] Fabry-Pérot interference
  Status: VALIDATED

PHASE 2: Lossy Modal Dispersion
  [✓] TE11 cutoff wavelength
  [✓] Lossy wall Q-factor
  [✓] Evanescent decay
  [✓] Modal analysis report
  Status: VALIDATED

PHASE 3: Near-Field Radiative Transfer
  [✓] Gap ratio metric
  [✓] Fresnel coefficients (evanescent)
  [✓] Transmission coefficient
  [✓] Spectral heat flux integral
  Status: VALIDATED

INTEGRATION TESTS
  [✓] Ray tracer + Phase 1
  [✓] Simulator + Phase 2
  [✓] Simulator + Phase 3
  Status: VALIDATED

BACKWARD COMPATIBILITY
  [✓] Existing tests pass
  [✓] Results within 1%
  Status: CONFIRMED

PERFORMANCE
  [✓] All operations within timing budgets
  Status: ACCEPTABLE

OVERALL STATUS: ✓ READY FOR PUBLICATION
=============================================================
```

---

**Document Version**: 1.0  
**Status**: Ready for implementation  
**Last Updated**: August 21, 2026
