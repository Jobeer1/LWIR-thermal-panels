# Accuracy Improvement Plan for Monte Carlo Ray Tracing

## Critical Issue: Thin Wall Optical Depth Physics

### Current Problem
The simulation assumes walls have bulk material emissivity regardless of thickness. For 100nm walls at 300K (λ_peak ≈ 9.7μm), with absorption depth δ ≈ 3.61μm:
- Bulk emissivity assumption: ε_wall ≈ 0.8-0.95
- **Actual thin-film emissivity**: ε_eff ≈ 0.8 × [1 - exp(-100nm/3.61μm)] ≈ 0.022
- **Error factor**: ~35-40× overestimation of wall emission

This fundamentally breaks the physics of anisotropic decoupling since walls appear much more emissive than they actually are.

## Phase 1: Thin-Film Optics Implementation (Highest Priority)

### 1.1 Add Material Absorption Depth Database
Create wavelength-dependent absorption depth data for common materials:
```python
# Example: Alumina (Al₂O₃) at 300K, λ ≈ 10μm
MATERIAL_ABSORPTION_DEPTH = {
    'alumina': {
        'wavelengths_um': [5, 8, 10, 12, 15],
        'depths_um': [2.1, 3.0, 3.61, 4.2, 5.0]
    },
    'cnt_forest': {
        'wavelengths_um': [2, 5, 10, 20],
        'depths_um': [0.5, 1.2, 2.5, 5.0]  # CNTs have shorter penetration
    }
}
```

### 1.2 Implement Beer-Lambert Thin-Film Correction
```python
def effective_emissivity_thin_film(bulk_eps: float, thickness_um: float, 
                                   wavelength_um: float, material: str) -> float:
    """
    ε_eff(λ) = ε_bulk × [1 - exp(-t/δ(λ))]
    
    For walls thinner than absorption depth, emissivity scales linearly with t/δ.
    For thick walls (t > 5δ), approaches bulk value.
    """
    delta = get_absorption_depth(material, wavelength_um)
    if delta <= 0:
        return bulk_eps
    
    optical_thickness = thickness_um / delta
    if optical_thickness > 5.0:
        return bulk_eps  # effectively bulk material
    
    return bulk_eps * (1.0 - math.exp(-optical_thickness))
```

### 1.3 Integrate into Monte Carlo Ray Tracing
Update photon absorption probability at walls:
```python
# Current (incorrect):
if surface == 'wall':
    eps = eps_walls  # constant bulk value

# New (correct):
if surface == 'wall':
    # Get photon wavelength (already sampled from Planck)
    lam_um = photon.wavelength_um  
    # Calculate effective emissivity for this wavelength and wall thickness
    eps = effective_emissivity_thin_film(
        bulk_eps=eps_walls_bulk,
        thickness_um=wall_thickness_um,
        wavelength_um=lam_um,
        material=wall_material
    )
```

## Phase 2: Statistical Convergence Improvements

### 2.1 Adaptive Photon Counting
- Implement convergence monitoring during MC runs
- Continue sampling until relative error < target (e.g., 1%)
- Minimum photons: 100k for p_esc ~ 0.01

### 2.2 Variance Reduction Techniques
- Improve Russian Roulette thresholds
- Implement importance sampling for deep cavities
- Use splitting/roulette for high-weight paths

### 2.3 Confidence Interval Calculation Fix
Current CI assumes binomial distribution, but weighted MC has different variance:
```python
# Current (simplistic):
ci = 1.96 * sqrt(p * (1-p) / n)

# Improved (weighted estimator):
var_w = (sum(w_i²)/n - (sum(w_i)/n)²)  # weighted variance
ci = 1.96 * sqrt(var_w / n)
```

## Phase 3: Waveguide Physics Validation

### 3.1 Cutoff Wavelength Verification
Verify λ_c calculations match literature:
- Circular TE11: λ_c = 1.706·diameter (correct in code)
- Rectangular TE10: λ_c = 2·min(width, depth) (narrowest transverse dimension)
- Check effective diameter for wall thickness effects

### 3.2 Evanescent Decay Implementation Check
Current implementation appears correct:
```python
def evanescent_decay_length(lambda_c_um, lambda_um):
    ratio = lambda_c_um / lambda_um
    return (lambda_c_um / (2*math.pi)) / sqrt(1 - ratio²)
```
Need to verify numerical stability near λ = λ_c (ratio → 1).

## Phase 4: Multi-Scale Validation Framework

### 4.1 Analytical Benchmarks
1. **Flat Plate Limit**: cavity_depth → 0, should match view factor formulas
2. **Blackbody Cavity**: ε_walls = 1.0, should give ε_B = 1.0
3. **Infinite Depth**: H → ∞, p_esc → 0, ε_B → 0
4. **Thermal Equilibrium**: T_A = T_B, net flux → 0 (energy conservation)

### 4.2 Comparison to Published Results
- PAA honeycomb cavities: Compare to Sprafke et al. (Adv. Opt. Mater. 2013)
- CNT forests: Compare to Mizuno et al. (PNAS 2009)
- Verify anisotropic decoupling ratios match literature

## Phase 5: User Interface Updates

### 5.1 New Input Parameters
Add to UI:
- Wall material selection (alumina, CNT, silver, etc.)
- Wall thickness input (µm)
- Absorption depth database (with references)
- Convergence tolerance setting

### 5.2 Enhanced Diagnostics
Display:
- Effective wall emissivity vs. wavelength plot
- Optical thickness (t/δ) for each material
- Thin-film correction factor applied
- Convergence progress during MC runs

## Implementation Timeline

### Week 1: Core Physics Fix
1. Implement material absorption depth database
2. Add thin-film emissivity calculation
3. Integrate into ray tracer
4. Basic validation tests

### Week 2: Statistical Improvements
1. Adaptive photon counting
2. Improved variance reduction
3. Better confidence intervals
4. Performance optimization

### Week 3: Validation & Documentation
1. Analytical benchmark suite
2. Literature comparison tests
3. Update documentation
4. User interface updates

## Expected Accuracy Improvements

| Metric | Current Error | Target After Fix |
|--------|---------------|------------------|
| Wall ε_eff at 300K | 35-40× overestimation | < 5% error |
| ε_B confidence interval | ±10.5% (relative) | < 2% (relative) |
| Net flux at equilibrium | Non-zero | < 0.1% of σT⁴ |
| Runtime for 1% error | N/A (fixed 20k photons) | Adaptive, ~30s |

## Key Physics References

1. **Thin-film optics**: Born & Wolf, Principles of Optics (transfer-matrix methods)
2. **Absorption depths**: Palik, Handbook of Optical Constants
3. **PAA cavities**: Sprafke et al., Adv. Opt. Mater. 2013
4. **CNT forests**: Mizuno et al., PNAS 2009
5. **Waveguide cutoff**: Jackson, Classical Electrodynamics (waveguide modes)

## Success Criteria

1. **Physical accuracy**: ε_eff matches thin-film optics within 5%
2. **Statistical reliability**: 95% CI < 2% relative error
3. **Energy conservation**: Net flux → 0 at thermal equilibrium
4. **Literature agreement**: Matches published decoupling ratios
5. **Performance**: < 30s runtime for typical cases with 1% error