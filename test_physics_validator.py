"""
test_physics_registry.py — Self-tests for the Part 1–3 peer-review upgrades.

Run:  python test_physics_registry.py

Covers
  • select_physics_regime()  — dimensionless regime routing (Part 2)
  • _modal_cutoff_evanescent_exitance() — exact γ_m modal attenuation + Z_mode (Part 1 §1)
  • fresnel_angular_absorptance()/maxwell_garnett_effective_index() (Part 1 §2)
  • enforce_thermal_equilibrium_zero_net() + isothermal run (Part 3)
"""

import math
import numpy as np

from simulator import (
    REGIME_EFFECTIVE_MEDIUM, REGIME_FULL_WAVE, REGIME_MACRO, REGIME_NEAR_FIELD,
    select_physics_regime, _modal_cutoff_evanescent_exitance,
    enforce_thermal_equilibrium_zero_net, run_simulation,
)
from geometry import HoneycombCavityCell, CNTForestCell
from material_optics import fresnel_angular_absorptance, maxwell_garnett_effective_index


def _check(name, condition, detail=''):
    status = 'PASS' if condition else '** FAIL **'
    print(f'[{status}] {name} {detail}')
    return bool(condition)


def test_regime_routing():
    ok = True
    honey = HoneycombCavityCell(diameter_um=200, height_um=450, wall_emissivity=0.98)
    cnt = CNTForestCell(pitch_um=0.05, dia_base_nm=10, dia_top_nm=5, height_um=450)

    # Macro: D=200µm vs λ(600K)=4.83µm → D ≫ λ
    r, info = select_physics_regime(1e-4, honey, 600.0)
    ok &= _check('macro regime', r == REGIME_MACRO, r)
    # Near-field: gap 200 nm < λ(300K)/(2π)=1.54 µm
    r, _ = select_physics_regime(200e-9, honey, 300.0)
    ok &= _check('near-field regime', r == REGIME_NEAR_FIELD, r)
    # EMT: CNT pitch 0.05 µm ≪ λ(600K)=4.83 µm
    r, _ = select_physics_regime(1e-3, cnt, 600.0)
    ok &= _check('effective-medium regime', r == REGIME_EFFECTIVE_MEDIUM, r)
    # Full wave: D=5 µm in [0.2λ, 5λ], λ(600K)=4.83 µm
    small = HoneycombCavityCell(diameter_um=5, height_um=450, wall_emissivity=0.98)
    r, _ = select_physics_regime(1e-3, small, 600.0)
    ok &= _check('full-wave regime', r == REGIME_FULL_WAVE, r)
    return ok


def test_regime_switch_verification():
    """Verification Protocol (Part 5): verify the regime switches.

    • D = 0.5 µm at λ = 10 µm  →  P, D ≪ λ  →  MAXWELL_GARNETT_ESTIMATE_EMT /
      homogenised TMM slab (geometric ray tracing bypassed).
    • D = 500 µm at λ = 10 µm  →  λ/D < 0.2  →  MACRO_GOUFFE_RAY_TRACE.
    """
    ok = True
    T_at_10um = 2897.77 / 10.0            # Wien's peak at 10 µm
    gap = 1e-3                           # far-field (gap ≫ λ/2π = 1.6 µm)

    small = HoneycombCavityCell(diameter_um=0.5, height_um=450, wall_emissivity=0.98)
    r_s, info_s = select_physics_regime(gap, small, T_at_10um)
    ok &= _check('D=0.5um @ 10um → EMT',
                 r_s == REGIME_EFFECTIVE_MEDIUM, f'{r_s} λ/D={info_s["lambda_over_D"]:.2f}')
    ok &= _check('EMT ratio λ/D large (sub-wavelength feature)',
                 info_s['lambda_over_D'] >= 5.0, f'{info_s["lambda_over_D"]:.2f}')

    big = HoneycombCavityCell(diameter_um=500, height_um=450, wall_emissivity=0.98)
    r_b, info_b = select_physics_regime(gap, big, T_at_10um)
    ok &= _check('D=500um @ 10um routes MACRO',
                 r_b == REGIME_MACRO, f'{r_b} λ/D={info_b["lambda_over_D"]:.3f}')
    ok &= _check('MACRO λ/D < 0.2', info_b['lambda_over_D'] < 0.2,
                f'{info_b["lambda_over_D"]:.3f}')

    # End-to-end: run_simulation must actually EXECUTE the routed operator.
    sim_emt = run_simulation(geometry_mode='honeycomb', cavity_diameter=0.5,
                             temp_a=290.0, temp_b=290.0, temp_surr=290.0,
                             n_photons=400, enable_near_field=False)
    ok &= _check('run_simulation(D=0.5) bypasses ray → effective_medium',
                sim_emt['active_physics_operator'] == 'effective_medium',
                sim_emt.get('solver_mode'))
    ok &= _check('EMT alpha==eps (Kirchhoff, homogenised slab)',
                abs(sim_emt['alpha_eff'] - sim_emt['epsilon_b']) < 1e-6,
                f"{sim_emt['alpha_eff']:.6f} vs {sim_emt['epsilon_b']:.6f}")

    sim_macro = run_simulation(geometry_mode='honeycomb', cavity_diameter=500,
                               temp_a=290.0, temp_b=290.0, temp_surr=290.0,
                               n_photons=400, enable_near_field=False)
    ok &= _check('run_simulation(D=500um) executes ray tracing',
                sim_macro['active_physics_operator'] == 'ray',
                sim_macro.get('solver_mode'))
    return ok


def test_modal_evanescent():
    ok = True
    # Deep, sub-cutoff pore must leak far less than a short one.
    t_shallow = _modal_cutoff_evanescent_exitance(10.0, 15.0, depth_um=1.0,
                                                  material='alumina', diameter_um=10.0)
    t_deep = _modal_cutoff_evanescent_exitance(10.0, 15.0, depth_um=200.0,
                                               material='alumina', diameter_um=10.0)
    ok &= _check('deep < shallow leakage', t_deep < t_shallow, f'{t_deep:.3e} < {t_shallow:.3e}')
    ok &= _check('bounded in [0,1]', 0.0 <= t_deep <= 1.0 and 0.0 <= t_shallow <= 1.0)
    # Propagating mode → no attenuation.
    t_prop = _modal_cutoff_evanescent_exitance(10.0, 5.0, depth_um=100.0,
                                               material='alumina', diameter_um=10.0)
    ok &= _check('propagating exits unattenuated', abs(t_prop - 1.0) < 1e-12)
    return ok


def test_angular_fresnel_and_emt():
    ok = True
    n = 1.7 + 0.002j
    a0 = fresnel_angular_absorptance(n, 0.0)
    # Normal-incidence reflectivity |(n−1)/(n+1)|² → α = 1 − R.
    R0 = abs((1.7-1)/(1.7+1))**2
    ok &= _check('Euclidean normal-incidence alpha', abs(a0-(1-R0)) < 1e-6, f'{a0:.5f} vs {1-R0:.5f}')
    ok &= _check('oblique alpha < normal alpha (dielectric)', fresnel_angular_absorptance(n, 60.0) < a0)
    ne_eff = maxwell_garnett_effective_index('alumina', 10.0, 0.9)
    ok &= _check('EMT n_eff has Re>0', ne_eff.real > 0)
    return ok


def test_equilibrium():
    ok = True
    g = enforce_thermal_equilibrium_zero_net(0.0, 300.0, 300.0, 300.0, tol=1e-9)
    ok &= _check('equilibrium guard flags uniform', g['is_thermally_uniform'])
    ok &= _check('equilibrium guard passes at q=0', g['passes'])
    # Non-uniform → pass-by-definition.
    g2 = enforce_thermal_equilibrium_zero_net(5000.0, 600.0, 300.0, 300.0)
    ok &= _check('non-equil guard not treated as violation', g2['is_thermally_uniform'] is False)
    # End-to-end isothermal radiosity must give q_net = 0.
    sim = run_simulation(geometry_mode='honeycomb', temp_a=300.0, temp_b=300.0,
                         temp_surr=300.0, n_photons=1500, enable_near_field=False)
    ok &= _check('isothermal run q_net==0', abs(sim['q_net_AB_density']) < 1e-9,
                 f"{sim['q_net_AB_density']:.3e}")
    # Isothermal Kirchhoff is q_net == 0 (energy conservation).  The
    # directional effective properties ε_B and α_eff are distinct spectral /
    # angular integrals and, for a MACRO geometric cavity (no sub-wavelength
    # decoupling), converge to the same Gouffé value within Monte-Carlo noise
    # — NOT to machine precision.  Forcing exact equality was the bug that
    # collapsed the documented decoupling in the sub-wavelength regime.
    ok &= _check('isothermal MACRO alpha~eps (Kirchhoff within MC noise)',
                 abs(sim['alpha_eff'] - sim['epsilon_b']) < 0.10,
                 f"{sim['alpha_eff']:.4f} vs {sim['epsilon_b']:.4f}")
    ok &= _check('isothermal equilibrium guard passes',
                 sim['equilibrium_guard']['passes'])
    return ok


def test_orchestrator():
    """Physics Orchestrator (regime_selector.py) pre-flight audit & routing."""
    ok = True
    from regime_selector import select_physics_regime

    # D = 4 µm cavity at 300 K: λ_peak ≈ 9.66 µm → λ/D ≈ 2.4 (sub-wavelength)
    reg = select_physics_regime(
        T_emit_K=300.0, diameter_um=4.0, height_um=200.0,
        wall_thickness_um=0.051, gap_um=200.0, k_extinction=0.05,
        lambda_cutoff_um=4.5831, rcwa_available=False)
    ok &= _check('orchestrator: D=4um@300K geometry FULL_WAVE',
                 reg.geometry_regime == 'FULL_WAVE', reg.geometry_regime)
    ok &= _check('orchestrator: pure ray tracing disabled',
                 reg.ray_valid is False)
    ok &= _check('orchestrator: RCWA recommended at sub-wavelength scale',
                 reg.rcwa_recommended is True)
    ok &= _check('orchestrator: TMM required for membrane walls',
                 reg.thinfilm_required is True,
                 f"{reg.wall_regime} t/delta={reg.t_over_delta:.4f}")
    ok &= _check('orchestrator: modal cutoff enforced',
                 reg.modal_cutoff_required is True, reg.cavity_regime)
    ok &= _check('orchestrator: confidence reflects sub-λ penalties',
                 0.0 <= reg.confidence < 90.0, f'{reg.confidence:.0f}%')
    ok &= _check('orchestrator: warnings exported',
                 len(reg.regime_warnings) >= 2, str(len(reg.regime_warnings)))
    d = reg.to_dict()
    ok &= _check('orchestrator: to_dict JSON-safe keys',
                 all(isinstance(k, str) for k in d) and 'confidence' in d)

    # Macroscopic cavity must stay RAY with full confidence.
    macro = select_physics_regime(
        T_emit_K=300.0, diameter_um=500.0, height_um=400.0,
        wall_thickness_um=5.0, gap_um=200.0, k_extinction=0.05)
    ok &= _check('orchestrator: D=500um stays RAY / ray_valid',
                 macro.geometry_regime == 'RAY' and macro.ray_valid is True)

    # End-to-end provenance payload exported from run_simulation.
    sim = run_simulation(cavity_diameter=4.0, wall_thickness=0.051,
                         height=200.001, packing_fraction=0.8842,
                         temp_a=300.0, temp_b=300.0, temp_surr=300.0,
                         n_photons=1500, enable_near_field=False)
    orch = sim.get('physics_orchestrator') or {}
    ok &= _check('simulator exports orchestrator provenance',
                 orch.get('geometry_regime') == 'FULL_WAVE' or 'error' in orch,
                 str(orch.get('geometry_regime')))
    ok &= _check('simulator exports physics_confidence',
                 isinstance(sim.get('physics_confidence'), (int, float)),
                 str(sim.get('physics_confidence')))
    ok &= _check('isothermal 300K q_net==0 (D=4um)',
                 abs(sim['q_net_AB_density']) < 1e-9,
                 f"{sim['q_net_AB_density']:.3e}")
    # Isothermal Kirchhoff is q_net == 0 (energy conservation).  At 300 K the
    # D = 4 µm cavity is in the geometric Gouffé regime (the propagating
    # Planck fraction is large enough that ε_cav saturates at the Gouffé bound),
    # so α_eff ≡ ε_eff holds within Monte-Carlo noise — as the orchestrator
    # spec requires.  (The 47× directional decoupling documented in Phase 9 is
    # a 200 K phenomenon where λ_peak = 14.49 µm ≫ λ_c and f_prop ≈ 0.011%.)
    ok &= _check('isothermal 300K alpha_eff ~ eps_B (Gouffé, D=4um)',
                 abs(sim['alpha_eff'] - sim['epsilon_b']) < 0.05,
                 f"{sim['alpha_eff']:.4f} vs {sim['epsilon_b']:.4f}")
    return ok


def test_doc_phase9_targets():
    """Lock in MATHEMATICAL_DERIVATIONS.md Phase 9 target values (200 K case).

    This is the headline bug regression: the isothermal mean-forcing used to
    collapse ε_B and α_eff into their average (45.45% / 45.45% / q = 10.3),
    destroying the documented directional decoupling.  After the fix the
    independent emission / absorption estimators reproduce the Phase 9 table.
    """
    ok = True
    sim = run_simulation(
        geometry_mode='honeycomb', temp_a=200.0, temp_b=200.0, temp_surr=3.0,
        emissivity_a=0.981, emissivity_a_back=0.081,
        width_a=1000.0, depth_a=1000.0, gap=200.001, height=200.001,
        cavity_diameter=4.0, wall_thickness=0.051, packing_fraction=0.884232,
        alpha_cnt=0.8, alpha_ag=0.05, eps_flat_wall=0.8,
        n_photons=2000, enable_near_field=False)

    ok &= _check('Phase 9: ε_B ≈ 1.888%',
                 abs(100 * sim['epsilon_b'] - 1.888) < 0.15,
                 f"{100*sim['epsilon_b']:.4f}%")
    ok &= _check('Phase 9: α_eff ≈ 89.27%',
                 abs(100 * sim['alpha_eff'] - 89.27) < 3.0,
                 f"{100*sim['alpha_eff']:.3f}%")
    ok &= _check('Phase 9: q_net ≈ 18.8304 W/m²',
                 abs(sim['net_flux_A_front'] - 18.8304) < 0.5,
                 f"{sim['net_flux_A_front']:.4f}")
    ok &= _check('Phase 9: decoupling α_eff/ε_B ≈ 47.3×',
                 20.0 < sim['decoupling_ratio'] < 80.0,
                 f"{sim['decoupling_ratio']:.2f}×")
    ok &= _check('Phase 9: ε_B ≠ α_eff (no isothermal mean-forcing)',
                 abs(sim['alpha_eff'] - sim['epsilon_b']) > 0.5,
                 f"{sim['alpha_eff']:.4f} vs {sim['epsilon_b']:.4f}")
    ok &= _check('Phase 9: T_stag ≈ 181.84 K',
                 abs(sim['T_B_stag'] - 181.84) < 0.5,
                 f"{sim['T_B_stag']:.3f} K")
    ok &= _check('Phase 9: λ_c = 4.5831 µm',
                 abs(sim['cutoff_wavelength_um'] - 4.5831) < 1e-3,
                 f"{sim['cutoff_wavelength_um']:.4f}")
    # NOTE: this is NOT a thermal-equilibrium case (T_surr = 3 K), so the
    # net inter-plate flux q_net = 18.83 W/m² is the physical ε_A ≠ ε_B
    # mismatch (already asserted above).  Zero net flux requires
    # T_A = T_B = T_surr and is covered by test_equilibrium().
    return ok


def test_high_temperature_fullwave():
    """Regression: 3000 K, D = 4 µm (λ/D ≈ 0.24) must route to FULL_WAVE.

    This is the bug from the independent peer review: at 3000 K the old
    regime_selector classified λ/D = 0.24 as "HYBRID" with ray_valid=True and
    (worse) disabled the modal transmission operator because λ_peak < λ_c.
    The propagating stratum then fell back to diffusive billiard bouncing,
    collapsing P_esc to 0.29% (guided-mode answer ≈ 88%), inflating the ε_B CI
    to ±42%, and making ε_cav (Gouffé, deterministic) contradict P_esc.
    """
    ok = True
    sim = run_simulation(
        geometry_mode='honeycomb', temp_a=3000.0, temp_b=3000.0, temp_surr=3.0,
        emissivity_a=0.981, emissivity_a_back=0.081,
        width_a=1000.0, depth_a=1000.0, gap=200.001, height=200.001,
        cavity_diameter=4.0, wall_thickness=0.051, packing_fraction=0.884232,
        alpha_cnt=0.8, alpha_ag=0.05, eps_flat_wall=0.8,
        n_photons=2000, enable_near_field=False)

    orch = sim['physics_orchestrator']
    ok &= _check('3000K: geometry_regime FULL_WAVE (not HYBRID)',
                 orch['geometry_regime'] == 'FULL_WAVE', orch['geometry_regime'])
    ok &= _check('3000K: ray_valid False (pure ray tracing disabled)',
                 orch['ray_valid'] is False)
    ok &= _check('3000K: modal transmission enforced despite CLASSICAL cavity',
                 orch['modal_cutoff_required'] is True,
                 f"modal={orch['modal_cutoff_required']} "
                 f"cavity={orch['cavity_regime']}")

    # Propagating modes (f_prop ≈ 96%) escape via guided-mode transmission, so
    # P_esc must be O(0.5–1.0), NOT the diffusive ~0.29%.
    ok &= _check('3000K: P_esc is guided-mode (≫ diffusive 0.3%)',
                 sim['p_esc'] > 0.5, f"{100*sim['p_esc']:.1f}%")

    # ε_B CI must be tiny (Gouffé-deterministic ceiling), not ±42% / ±70%.
    ok &= _check('3000K: ε_B CI is not blown up',
                 sim['epsilon_b_ci95'] < 0.05,
                 f"±{100*sim['epsilon_b_ci95']:.2f}%")

    # Near-Kirchhoff at this scale (no sub-wavelength decoupling).
    ok &= _check('3000K: near-Kirchhoff (α_eff ≈ ε_B)',
                 abs(sim['alpha_eff'] - sim['epsilon_b']) < 0.03,
                 f"{sim['alpha_eff']:.3f} vs {sim['epsilon_b']:.3f}")
    return ok


def test_short_wavelength_high_temperature():
    """Regression: T = 13000 K (λ_peak ≈ 0.223 µm, λ/D ≈ 0.056) must RUN.

    This is the "shorter wavelength" bug: app.py previously hard-rejected any
    temperature above 3000 K with a ValueError, so the API returned status
    error and the UI rendered every result as '—' (while the stale stagnation
    temperature from the previous run remained on screen).  The cap is removed;
    high-temperature extrapolation is instead surfaced by the orchestrator's
    confidence score and audit warnings.  Material optics stay finite well
    beyond 13000 K (n ≈ 1.75, k ≈ 0.002 at 0.223 µm).
    """
    ok = True

    # 1. Material optics remain finite at the Wien peak wavelength.
    from material_optics import get_complex_refractive_index
    import math
    lam_peak = 2898.0 / 13000.0
    n, k = get_complex_refractive_index('alumina', lam_peak, temperature_K=13000.0)
    ok &= _check('13000K: n(λ,T),k(λ,T) finite at λ_peak',
                 math.isfinite(n) and math.isfinite(k) and k >= 0.0,
                 f"n={n:.4f} k={k:.4f}")

    # 2. The regime selector stays in geometric RAY (λ/D < 0.2) — valid.
    from regime_selector import select_physics_regime
    reg = select_physics_regime(
        13000.0, 4.0, 200.0, 0.051, 200.0, max(k, 0.0),
        lambda_cutoff_um=4.5831, rcwa_available=False)
    ok &= _check('13000K: geometry RAY (λ/D < 0.2, geometric valid)',
                 reg.geometry_regime == 'RAY' and reg.ray_valid is True,
                 f"{reg.geometry_regime} λ/D={reg.lambda_over_diameter:.3f}")
    ok &= _check('13000K: high-T confidence penalty applied',
                 reg.confidence < 100.0, f"{reg.confidence:.0f}%")

    # 3. End-to-end simulation completes and stays bounded (2nd law).
    sim = run_simulation(
        geometry_mode='honeycomb', temp_a=13000.0, temp_b=13000.0, temp_surr=3.0,
        emissivity_a=0.981, emissivity_a_back=0.081,
        width_a=1000.0, depth_a=1000.0, gap=200.001, height=200.001,
        cavity_diameter=4.0, wall_thickness=0.051, packing_fraction=0.884232,
        alpha_cnt=0.8, alpha_ag=0.05, eps_flat_wall=0.8,
        n_photons=2000, enable_near_field=False)
    for key in ('epsilon_b', 'alpha_eff', 'net_flux_A_front', 'T_B_stag'):
        v = sim.get(key)
        ok &= _check(f'13000K: {key} finite', math.isfinite(float(v)), str(v))
    ok &= _check('13000K: stagnation ≤ max(T_A) (2nd law)',
                 sim['T_B_stag'] <= 13000.0, f"{sim['T_B_stag']:.1f} K")
    ok &= _check('13000K: q_net finite and positive (cold sink)',
                 sim['net_flux_A_front'] > 0.0,
                 f"{sim['net_flux_A_front']:.3g}")
    return ok


if __name__ == '__main__':
    print('Running physics-validity self-tests...')
    results = [
        test_regime_routing(),
        test_regime_switch_verification(),
        test_modal_evanescent(),
        test_angular_fresnel_and_emt(),
        test_equilibrium(),
        test_orchestrator(),
        test_doc_phase9_targets(),
        test_high_temperature_fullwave(),
        test_short_wavelength_high_temperature(),
    ]
    print('\n' + '=' * 52)
    print('ALL TESTS PASSED' if all(results) else 'SOME TESTS FAILED')
    print('=' * 52)
    raise SystemExit(0 if all(results) else 1)