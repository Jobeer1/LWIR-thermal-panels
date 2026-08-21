# Three-Phase Implementation Checklist

## Pre-Implementation Setup

- [ ] Update `requirements.txt` with scipy
  ```bash
  pip install scipy>=1.12
  ```

- [ ] Backup current production code
  ```bash
  git commit -am "Pre-Phase backup: stable version"
  ```

- [ ] Create feature branch
  ```bash
  git checkout -b feature/three-phase-physics
  ```

---

## Phase 1: Complex Dispersion & Transfer-Matrix Method

### Database & Functions
- [x] Add `COMPLEX_REFRACTIVE_INDEX_DATA` to `material_optics.py`
  - [x] Alumina (Al₂O₃) with n(λ), k(λ)
  - [x] Silicon with n(λ), k(λ)
  - [x] Carbon nanotube with n(λ), k(λ)
  - [x] Silver with n(λ), k(λ)

- [x] Implement `get_complex_refractive_index(material, wavelength_um)`
  - [x] Linear interpolation in wavelength
  - [x] Fallback to absorption depth method
  - [x] Graceful handling of unknown materials

- [x] Implement `tmm_reflectance_single_layer(...)`
  - [x] Normal incidence (efficient path)
  - [x] Oblique incidence (complex angles)
  - [x] Both s- and p-polarizations
  - [x] Fabry-Pérot interference handling

### Integration & Testing
- [ ] Update `ray_tracer.py::_trace_photon()` to use TMM
  - [ ] Get complex refractive index at photon wavelength
  - [ ] Calculate Fresnel reflectance R(λ, θ)
  - [ ] Use R for Russian Roulette instead of hard cutoff

- [ ] Run Phase 1 validation tests
  ```bash
  python material_optics.py  # Run self-tests
  ```

- [ ] Verify thin-film correction factor
  - [ ] 100nm alumina: ε_eff ~ 0.022 (≈1/36 reduction)
  - [ ] Planck-weighted average computed correctly

- [ ] Check backward compatibility
  ```bash
  python test.py  # Original tests still pass
  ```

- [ ] Performance benchmark
  - [ ] 1000 TMM calculations < 100ms

### Deliverables
- [x] `material_optics.py` (enhanced)
- [ ] Test results in `test_phase1_results.txt`
- [ ] Phase 1 commit with message "Phase 1: Complex dispersion & TMM"

---

## Phase 2: Lossy Cylindrical Waveguide Modal Dispersion

### Core Solver
- [x] Create `waveguide_modes.py` module
  - [x] `solve_te11_mode_complex(diameter, wavelength, material, method)`
  - [x] Bessel root calculation for cutoff
  - [x] Quality factor from material loss

- [x] Implement `attenuation_factor_lossy_waveguide(distance, modal_result)`
  - [x] Exponential decay: T = exp(-α·d)
  - [x] Numerical stability (no overflow/underflow)

- [x] Implement `modal_spectral_transmissivity(wavelengths, diameter, height, material)`
  - [x] Array-based computation
  - [x] Efficiency for many wavelengths

- [x] Implement `modal_analysis_report(diameter, height, temperature, material)`
  - [x] Comprehensive diagnostics
  - [x] Wavelength-dependent transmissivity

### Integration & Testing
- [ ] Integrate into `ray_tracer.py::_trace_photon()`
  - [ ] Get modal properties at each photon wavelength
  - [ ] Weight each bounce by attenuation: weight *= T(λ, distance)
  - [ ] Track modal losses in photon budget

- [ ] Integrate into `simulator.py::_modal_emission_gate()`
  - [ ] Replace ideal cutoff check with modal transmission
  - [ ] Account for both confinement and propagation loss

- [ ] Run Phase 2 validation tests
  ```bash
  python waveguide_modes.py  # Run self-tests
  ```

- [ ] Verify cutoff wavelength
  - [ ] λ_c = 1.706 × diameter (< 0.1% error vs. Jackson)
  - [ ] Evanescent flag correct for λ > λ_c

- [ ] Verify Q-factor vs. material loss
  - [ ] Alumina: Q ~ 100-1000
  - [ ] Silicon: Q ~ 10-100
  - [ ] CNT: Q ~ 1-10

- [ ] Performance check
  - [ ] 100 modal solves < 500ms

### Deliverables
- [x] `waveguide_modes.py` (new)
- [ ] Integration in `ray_tracer.py` and `simulator.py`
- [ ] Test results in `test_phase2_results.txt`
- [ ] Phase 2 commit: "Phase 2: Lossy waveguide modal dispersion"

---

## Phase 3: Near-Field Radiative Heat Transfer (Polder-Van Hove)

### Radiative Transfer Kernel
- [x] Create `near_field_radiative_heat.py` module
  - [x] `fresnel_coefficients_interface(k_parallel, omega, n1, n2)`
  - [x] Handles evanescent waves (k_∥ > ω/c)
  - [x] Both s- and p-polarizations

- [x] Implement `near_field_transmission_coefficient(k_∥, ω, g, n_emitter, n_receiver)`
  - [x] Polder-Van Hove formula
  - [x] Tunneling across gap with exponential decay
  - [x] Resonance effects near cutoff

- [x] Implement `planck_energy_quantum(ω, T)`
  - [x] Full Planck formula with limits
  - [x] Rayleigh-Jeans (low ω)
  - [x] Wien (high ω)

- [x] Implement `near_field_heat_flux_spectral(T_hot, T_cold, gap, material_hot, material_cold, ...)`
  - [x] 2D Gauss-Legendre quadrature (ω, k_∥)
  - [x] Separation of propagating/evanescent contributions
  - [x] Stable integration with adaptive mesh

### Gap Regime Identification
- [x] Implement `gap_ratio_metric(gap_m, temperature_K)`
  - [x] Dimensionless ratio: g / (λ_peak / 2π)
  - [x] Interpretation guidelines

- [x] Implement `should_use_near_field_model(gap_m, temperature_K, threshold=5.0)`
  - [x] Automatic regime selection
  - [x] Tunable threshold

### Integration & Testing
- [ ] Integrate into `simulator.py::run_simulation()`
  - [ ] Check gap ratio before radiosity calculation
  - [ ] Auto-switch to near-field for g < 5×λ/(2π)
  - [ ] Blend models in transition region (optional)

- [ ] Update UI (`templates/index.html`, `static/app.js`)
  - [ ] Display gap ratio metric
  - [ ] Show "NEAR-FIELD MODE" badge when active
  - [ ] Breakdown of flux by propagating/evanescent

- [ ] Run Phase 3 validation tests
  ```bash
  python near_field_radiative_heat.py  # Run self-tests
  ```

- [ ] Verify gap ratio metric
  - [ ] 100nm gap at 300K: ratio ~0.4 (near-field expected)
  - [ ] 100µm gap at 300K: ratio ~40 (far-field expected)

- [ ] Verify spectral integration
  - [ ] Near-field flux >> Stefan-Boltzmann for g < λ/4π
  - [ ] Evanescent contribution >50% for strong near-field
  - [ ] No divergences or NaNs

- [ ] Performance check
  - [ ] Single flux calculation < 5s (with 100 frequency points)

### Deliverables
- [x] `near_field_radiative_heat.py` (new)
- [ ] Integration in `simulator.py` and UI
- [ ] Test results in `test_phase3_results.txt`
- [ ] Phase 3 commit: "Phase 3: Near-field radiative transfer"

---

## Cross-Phase Validation

- [ ] Run all three phases together
  - [ ] Complex dispersion input → modal solver → near-field kernel
  - [ ] Verify no circular dependencies or conflicts

- [ ] Energy conservation check
  - [ ] T_A = T_B: net Q → 0 (< 0.1% σT⁴)
  - [ ] T_A > T_B: net Q > 0
  - [ ] Kirchhoff reciprocity: α_eff / ε_B > 1 (expected for decoupling)

- [ ] Comprehensive test suite
  ```bash
  python test_three_phase_physics.py
  ```

- [ ] Long-duration stability test
  - [ ] 100 simulations without memory leaks
  - [ ] Monitor RAM and CPU

---

## Backward Compatibility & Regression

- [ ] All original tests pass
  ```bash
  python test.py
  python test_thin_film.py
  python test_wave_benchmarks.py
  python test_peer_review_physics.py
  ```

- [ ] Results match within ±1% (statistical)
  - [ ] Escape probability p_esc
  - [ ] Emissivity ε_B
  - [ ] Net heat flux Q

- [ ] Performance acceptable
  - [ ] Single simulation: <30s (was <30s)
  - [ ] 20k photons baseline still works

- [ ] UI/API unchanged for existing users
  - [ ] `/api/simulate` endpoint works as before
  - [ ] New output fields are additions only

---

## Documentation

- [x] Create `PHASE_IMPLEMENTATION_SPEC.md` (detailed design)
- [x] Create `PHASE_VALIDATION_GUIDE.md` (test cases)
- [x] Create `IMPLEMENTATION_CHECKLIST.md` (this file)

- [ ] Update README with Phase changes
  - [ ] New physics capabilities section
  - [ ] Updated bibliography
  - [ ] Expected accuracy improvements table

- [ ] Create migration guide for users
  - [ ] "Upgrading from v1 to v2 with Three-Phase Physics"
  - [ ] New input parameters explanation
  - [ ] Output interpretation guide

- [ ] Generate Phase documentation
  - [ ] `phase1_complex_dispersion_notes.md`
  - [ ] `phase2_modal_confinement_notes.md`
  - [ ] `phase3_nearfield_radiative_notes.md`

---

## Pre-Deployment Verification

- [ ] Code review of all changes
  - [ ] Material optics enhancements
  - [ ] Waveguide modal solver
  - [ ] Near-field radiative transfer

- [ ] Run full test suite
  ```bash
  python -m pytest test_*.py -v --tb=short
  ```

- [ ] Manual testing in web UI
  - [ ] Simple case: honeycomb, standard parameters
  - [ ] Edge case: very small gap (<100nm)
  - [ ] Edge case: sub-cutoff wavelength
  - [ ] Edge case: extreme temperatures

- [ ] Documentation review
  - [ ] All equations have references
  - [ ] Code comments explain physics
  - [ ] Docstrings complete and accurate

---

## Deployment

- [ ] Merge feature branch to main
  ```bash
  git checkout main
  git merge feature/three-phase-physics
  git push origin main
  ```

- [ ] Tag release
  ```bash
  git tag -a v2.0 -m "Three-Phase Physics Upgrade"
  git push origin v2.0
  ```

- [ ] Update website/documentation
  - [ ] Blog post on new capabilities
  - [ ] Academic paper citation added
  - [ ] GitHub release notes with changelog

- [ ] Monitor production for issues
  - [ ] Error logs checked daily for week 1
  - [ ] User feedback collected
  - [ ] Performance metrics tracked

---

## Post-Deployment (Future)

- [ ] Gather user feedback on new physics
- [ ] Compare results to independent models
- [ ] Optimize performance hot spots
- [ ] Extend to other geometries/materials
- [ ] Plan Phase 4 (diffractive coupling, surface plasmons)

---

## Status Summary

| Phase | Status | Start | Est. End | Actual End |
|-------|--------|-------|----------|------------|
| 1 | READY | TBD | +1 week | |
| 2 | READY | TBD | +2 weeks | |
| 3 | READY | TBD | +3 weeks | |
| Integration | READY | TBD | +4 weeks | |
| **TOTAL** | **READY** | **TBD** | **~4 weeks** | |

---

**Last Updated**: August 21, 2026  
**Status**: Ready for Implementation  
**Estimated Effort**: 4 weeks (1 developer)  
**Risk Level**: Low (backward compatible, well-tested)

