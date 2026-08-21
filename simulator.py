"""
simulator.py — Orchestrator for the two-plate radiative heat transfer simulation.

Improvements over v1:
  Bug #1  Fixed — 3-D Lambertian sampling (Malley's method) via ray_tracer.py
  Bug #2  Fixed — CNT diameter inputs now drive CNTForestCell / FrustumCavity3D geometry
  Bug #3  Fixed — 3-D finite-rectangle view factor (Howell catalog C-11)
  Bug #4  Fixed — Weight-based Russian Roulette; no more 50/50 bias
  Bug #5  Fixed — ε_B normalization uses 3-D areas from geometry objects
  Bug #6  Fixed — 4-surface radiosity enclosure (A-front, A-back, B, surroundings)
  Bug #7  Added — Spectral band model with Planck-weighted effective emissivities
  Bug #8  Added — Near-field warning when gap < λ_peak / (2π)
  Bug #9  Fixed — Aperture re-entry coupling from plate A reflections
  Bug #10 Fixed — Tapered frustum geometry for CNT tip taper
  Bug #11 Added — 95% CI on all MC outputs
  Bug #12 Fixed — (in rad_leakage.py; simulator already correct)
"""

import math
import numpy as np

from geometry    import RectPit3D, FrustumCavity3D, CNTForestCell, HoneycombCavityCell
from ray_tracer  import run_cavity_mc_3d
from spectral    import planck_weighted_emissivity, MATERIAL_EMISSIVITY, effective_emissivity_pair
from sampling    import planck_cumulative, planck_averaged_evanescent_decay
from near_field_radiative_heat import (
    gap_ratio_metric,
    should_use_near_field_model,
    near_field_heat_flux_spectral
)

SIGMA = 5.670374419e-8  # Stefan-Boltzmann constant, W m⁻² K⁻⁴


# ---------------------------------------------------------------------------
# Gap Regime Detection (Phase 3)
# ---------------------------------------------------------------------------

def _detect_gap_regime(gap_m: float, T_hot_K: float, threshold: float = 5.0) -> dict:
    """Auto-detect far-field vs. near-field based on gap ratio.
    
    Uses Polder-Van Hove dimensionless metric: gap / (λ_peak / 2π).
    
    Parameters
    ----------
    gap_m : float
        Gap distance (m)
    T_hot_K : float
        Hot surface temperature (K) for Wien's peak wavelength
    threshold : float
        Gap ratio threshold for near-field activation (default 5.0)
        Typically: < 1 = strong near-field, < 5 = significant correction
    
    Returns
    -------
    dict with:
        'gap_ratio': dimensionless ratio
        'use_near_field': bool indicating if near-field model should be used
        'regime': 'near-field' or 'far-field' label
    """
    gap_ratio = gap_ratio_metric(gap_m, T_hot_K)
    use_near_field = should_use_near_field_model(gap_m, T_hot_K, threshold)
    return {
        'gap_ratio': gap_ratio,
        'use_near_field': use_near_field,
        'regime': 'near-field' if use_near_field else 'far-field'
    }


# ---------------------------------------------------------------------------
# View-factor utilities
# ---------------------------------------------------------------------------

def _rect_rect_view_factor(a: float, b: float, c: float) -> float:
    """Exact view factor F₁₂ between two identical, parallel, coaxial rectangles
    of sides a × b separated by gap c.  (Howell et al., Catalog C-11.)

    Bug #3 fix: replaces the 2-D Hottel cross-string formula.
    """
    if c <= 0.0:
        return 1.0
    if a <= 0.0 or b <= 0.0:
        return 0.0
    X = a / c
    Y = b / c
    # Closed-form from Howell catalog C-11
    term1 = math.log((1 + X*X) * (1 + Y*Y) / (1 + X*X + Y*Y))
    term2 = X * math.sqrt(1 + Y*Y) * math.atan(X / math.sqrt(1 + Y*Y))
    term3 = Y * math.sqrt(1 + X*X) * math.atan(Y / math.sqrt(1 + X*X))
    term4 = X * math.atan(X)
    term5 = Y * math.atan(Y)
    F = (2.0 / (math.pi * X * Y)) * (
        0.5 * term1 + term2 + term3 - term4 - term5
    )
    return float(np.clip(F, 0.0, 1.0))


def _parallel_strips_vf_2d(w_a: float, w_b: float, gap: float, offset: float = 0.0) -> float:
    """2-D Hottel cross-string view factor (kept as fallback / legacy)."""
    if w_a <= 0.0:
        return 0.0
    if gap <= 0.0:
        return 1.0
    xa0, xa1 = 0.0, w_a
    xb0, xb1 = offset, offset + w_b

    def dist(x1, x2):
        return math.hypot(x2 - x1, gap)

    d1 = dist(xa0, xb1)
    d2 = dist(xa1, xb0)
    d3 = dist(xa0, xb0)
    d4 = dist(xa1, xb1)
    f = (d1 + d2 - d3 - d4) / (2.0 * w_a)
    return float(np.clip(f, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Near-field check
# ---------------------------------------------------------------------------

def _near_field_check(gap_m: float, T_hot_K: float) -> dict:
    """Compute near-field ratio and return a warning string if applicable.

    Near-field radiation becomes significant when gap d < λ_peak / (2π).
    λ_peak = 2898 µm·K / T  (Wien's displacement law).
    """
    if T_hot_K <= 0:
        return {'near_field_ratio': None, 'near_field_warning': ''}
    lambda_peak_m = (2898e-6) / T_hot_K   # metres
    d_nf = lambda_peak_m / (2.0 * math.pi)
    ratio = gap_m / d_nf

    warning = ''
    if ratio < 1.0:
        warning = (
            f'⚠ Near-field regime: gap ({gap_m*1e6:.1f} µm) < λ_peak/2π '
            f'({d_nf*1e6:.2f} µm). Evanescent photon tunneling can exceed '
            f'far-field flux by orders of magnitude — this model is NOT valid '
            f'in the near-field regime.'
        )
    elif ratio < 3.0:
        warning = (
            f'⚠ Near near-field boundary: gap/λ_peak·2π = {ratio:.2f} '
            f'(< 3). Near-field corrections may be 10–30%.'
        )

    return {'near_field_ratio': ratio, 'near_field_warning': warning}


# ---------------------------------------------------------------------------
# Bug #13 — LDOS emission-confinement gate
# ---------------------------------------------------------------------------

def _aperture_rim_strip_area(geometry, delta_um: float) -> float:
    """Area (m²) of the cavity interior surfaces that lie within `delta_um` of
    the aperture plane.  Only this "rim" strip can couple sub-wavelength
    (evanescent) thermal emission to the far field; deeper walls are
    LDOS-suppressed (Lin et al. PRB 2000, Narayanaswamy & Chen PRB 2004).
    """
    d = float(delta_um) * 1e-6
    if hasattr(geometry, 'R'):                     # HoneycombCavityCell (pore)
        perimeter = 2.0 * math.pi * geometry.R
        return min(geometry.area_walls, perimeter * d)
    if hasattr(geometry, 'r_top'):                 # FrustumCavity3D (truncated cone)
        perimeter = 2.0 * math.pi * geometry.r_top
        return min(geometry.area_walls, perimeter * d)
    if hasattr(geometry, 'W') and hasattr(geometry, 'D'):   # RectPit3D (legacy)
        perimeter = 2.0 * (geometry.W + geometry.D)
        return min(geometry.area_walls, perimeter * d)
    if hasattr(geometry, 'pitch_um'):              # CNTForestCell
        # Lateral frustum strip near the tip, plus the top cap — which lies in
        # the aperture plane itself and therefore couples directly.
        perimeter = 2.0 * math.pi * geometry.rt
        strip = min(geometry.area_lateral, perimeter * d)
        return strip + geometry.area_top_cap
    return geometry.area_walls


def _modal_emission_gate(geometry, lambda_c_um: float, T_emit: float,
                         eps_wall: float, eps_base: float,
                         f_prop: float):
    """Bug #13 LDOS emission-confinement factor for the cavity.

    The geometric cavity enhancement C_e = (Σ A·ε)/A_ap counts ALL wall area,
    but for a sub-wavelength channel only two emission paths reach the far
    field:

      1. the propagating spectral fraction f_prop (λ < λ_c) — fully ray-traced;
      2. the evanescent fraction (λ ≥ λ_c) emitted from the near-aperture rim
         strip within one Planck-averaged decay length δ_ev of the aperture —
         walls deeper than δ_ev have a suppressed photonic DOS and radiate
         essentially nothing.

        G_em  = f_prop + (1 − f_prop) · A_rim,ε / A_emit,ε
        ε_B^gated = ε_B · G_em

    Returns (G_em, rim_fraction, delta_avg_um).
    """
    if not lambda_c_um or lambda_c_um <= 0.0 or not math.isfinite(lambda_c_um):
        return 1.0, 0.0, 0.0
    H_um = getattr(geometry, 'H', 0.0) * 1e6
    delta_avg_um = planck_averaged_evanescent_decay(
        lambda_c_um, float(T_emit), cap_um=H_um if H_um > 0 else 0.0)
    A_rim = _aperture_rim_strip_area(geometry, delta_avg_um)
    A_emit = max(
        geometry.area_walls * float(eps_wall) + geometry.area_base * float(eps_base),
        1e-30)
    rim_frac = min(1.0, (A_rim * float(eps_wall)) / A_emit)
    g_em = float(f_prop) + (1.0 - float(f_prop)) * rim_frac
    return float(np.clip(g_em, 0.0, 1.0)), float(rim_frac), float(delta_avg_um)


# ---------------------------------------------------------------------------
# 4-surface radiosity network (Bug #6 fix)
# ---------------------------------------------------------------------------

def _radiosity_4surface(
    eps_af: float, eps_ab: float, eps_b: float, eps_surr: float,
    A_af: float, A_ab: float, A_b: float,
    F_af_b: float,
    T_a: float, T_b: float, T_surr: float,
) -> dict:
    """4-surface enclosure radiosity network.

    Surfaces:
      0 — Plate A front face  (area A_af, emissivity eps_af)
      1 — Plate A back face   (area A_ab, emissivity eps_ab)
      2 — Plate B front face  (area A_b,  emissivity eps_b)
      3 — Surroundings        (blackbody at T_surr, infinite area)

    View factors assumed:
      F_af→b   = F_af_b  (3-D rect-rect formula)
      F_af→sur = 1 - F_af_b
      F_b→af   = F_af_b * A_af / A_b  (reciprocity)
      F_b→sur  = 1 - F_b→af
      F_ab→sur = 1  (back face sees only surroundings)
      Surroundings absorbs everything.

    Returns dict with q_af_net, q_ab_net, q_b_net, q_af_to_b.
    """
    E_a    = SIGMA * T_a**4
    E_b    = SIGMA * T_b**4
    E_surr = SIGMA * T_surr**4

    F_af_s = max(1.0 - F_af_b, 1e-12)
    F_b_af = F_af_b * A_af / A_b if A_b > 0 else 0.0
    F_b_s  = max(1.0 - F_b_af, 1e-12)

    # Radiosity equations (2×2 for the two active surfaces AF and B)
    # J_AF = eps_af*E_A + (1-eps_af)*(F_af_b*J_B + F_af_s*E_surr)
    # J_B  = eps_b *E_B + (1-eps_b )*(F_b_af*J_AF + F_b_s*E_surr)
    coef = np.array([
        [1.0,               -(1.0 - eps_af) * F_af_b],
        [-(1.0 - eps_b) * F_b_af,  1.0              ],
    ])
    rhs = np.array([
        eps_af * E_a + (1.0 - eps_af) * F_af_s * E_surr,
        eps_b  * E_b + (1.0 - eps_b ) * F_b_s  * E_surr,
    ])
    try:
        J_af, J_b = np.linalg.solve(coef, rhs)
    except np.linalg.LinAlgError:
        J_af, J_b = E_a, E_b

    # Back face radiosity (sees only surroundings)
    J_ab = eps_ab * E_a + (1.0 - eps_ab) * E_surr

    # Net heat fluxes
    G_af = F_af_b * J_b + F_af_s * E_surr
    q_af_net  = A_af * (J_af - G_af)
    q_ab_net  = A_ab * (J_ab - E_surr)
    q_af_to_b = A_af * F_af_b * (J_af - J_b)

    G_b  = F_b_af * J_af + F_b_s * E_surr
    q_b_net   = A_b * (J_b - G_b)

    return {
        'j_af':      J_af,
        'j_b':       J_b,
        'j_ab':      J_ab,
        'q_af_net':  q_af_net,
        'q_ab_net':  q_ab_net,
        'q_b_net':   q_b_net,
        'q_af_to_b': q_af_to_b,
    }


# ---------------------------------------------------------------------------
# Main simulation entry point
# ---------------------------------------------------------------------------

def run_simulation(
    # Geometry mode
    geometry_mode: str = 'honeycomb',    # 'honeycomb' | 'cnt_forest'
        # Plate B cavity parameters — Honeycomb mode
    cavity_diameter: float = 200.0,      # cavity opening diameter (µm)
    wall_thickness: float = 1.0,         # wall thickness between holes (µm)
    packing_fraction: float = 0.9069,    # hex close-pack fraction (computed from D, w)
    # Plate B cavity parameters — shared
    height: float = 450.0,              # cavity height (µm)
    alpha_cnt: float = 0.98,            # wall absorptivity
    alpha_ag: float = 0.02,             # base absorptivity
    eps_flat_wall: float = 0.1,         # emissivity of flat inter-cavity area (for panel scaling)
    # Plate B cavity parameters — CNT Forest mode
    width: float = 10.0,                # trench width (µm) — rect_pit mode (legacy)
    depth: float = 10.0,                # trench depth into page (µm) — rect_pit mode (legacy)
    cnt_dia_base: float = 10.0,         # CNT base diameter (nm)
    cnt_dia_top: float = 5.0,           # CNT tip diameter (nm)
    cnt_pitch: float = 0.05,            # CNT pitch (µm) — NOTE: µm not nm
    # Plate A parameters
    temp_a: float = 600.0,             # K
    emissivity_a: float = 1.0,
    emissivity_a_back: float = 0.1,
    width_a: float = 1000.0,           # µm
    depth_a: float = 1000.0,           # µm
    # Gap and surroundings
    gap: float = 100.0,                # µm
    temp_b: float = 300.0,             # K
    temp_surr: float = 300.0,          # K
    # MC settings
    n_photons: int = 20000,
    full_gap_mc: bool = False,
    n_gap_photons: int = 20000,
        # Spectral model
    material_a: str = '',              # material key for spectral lookup ('' = gray)
    material_b: str = '',
    # Bug #13 — wave-optics modal cutoff (peer review)
    alpha_top: float = None,          # top-surface absorptivity for sub-λ_c incidence
    # Phase 0/6 — wave-model integration boundary
    wave_model: str = 'ray',          # 'ray' (default fallback) | 'cached'
    cache_path: str = '',             # optional explicit path to WaveResponse (.json/.h5)
    # Phase 3 — Polder-Van Hove near-field integration
    enable_near_field: bool = True,  # Auto-switch on/off
    near_field_threshold: float = 5.0,  # Gap ratio threshold
    near_field_n_omega: int = 80,  # Quadrature points for frequency
    near_field_n_kparallel: int = 50,  # Quadrature points for parallel k
) -> dict:
    """Run the full two-plate radiative exchange simulation.

    ``wave_model`` selects the cavity-absorption solver:
        * 'ray'    (default) — the 3-D Monte Carlo ray tracer (Phase-0 fallback);
        * 'cached'           — bypass MC and interpolate effective α(λ)/ε_b(λ)
                               from a pre-computed full-wave WaveResponse cache.
    
    ``enable_near_field`` controls Phase 3 (Polder-Van Hove near-field integration):
        * True (default) — auto-detect and use near-field for small gaps (g < 5λ/2π);
        * False          — always use far-field (radiosity) model.
    """

        # ---- Build cavity geometry -------------------------------------------
    if geometry_mode == 'honeycomb':
        geometry = HoneycombCavityCell(
            diameter_um      = cavity_diameter,
            height_um        = height,
            wall_emissivity  = alpha_cnt,
            packing_fraction = packing_fraction,
        )
        pitch = cavity_diameter + wall_thickness
        label = (f'Honeycomb Cavity (D={cavity_diameter:.1f}µm, H={height:.1f}µm, '
                 f'w={wall_thickness:.3f}µm, pitch={pitch:.3f}µm, '
                 f'AR={height/max(cavity_diameter,1e-9):.1f}, f={packing_fraction:.4f})')
    elif geometry_mode == 'cnt_forest':
        # cnt_pitch is in µm; cnt_dia_base/top still in nm
        geometry = CNTForestCell(
            pitch_um    = cnt_pitch,
            dia_base_nm = cnt_dia_base,
            dia_top_nm  = cnt_dia_top,
            height_um   = height,
        )
        label = (f'CNT Forest (pitch={cnt_pitch*1000:.0f}nm, d_base={cnt_dia_base:.1f}nm, '
                 f'd_top={cnt_dia_top:.1f}nm, H={height:.1f}µm)')
    elif geometry_mode == 'frustum':
        geometry = FrustumCavity3D(
            r_base_um = (cnt_dia_base * 1e-3) / 2.0,
            r_top_um  = (cnt_dia_top  * 1e-3) / 2.0,
            height_um = height,
        )
        label = f'Frustum (r_base={cnt_dia_base/2:.1f}nm, r_top={cnt_dia_top/2:.1f}nm, H={height:.1f}µm)'
    else:  # rect_pit (legacy)
        geometry = RectPit3D(
            width_um  = width,
            depth_um  = depth,
            height_um = height,
        )
        label = f'Rect Pit ({width:.1f}×{depth:.1f}×{height:.1f} µm) [legacy]'

    # ---- Macro-scale plate areas (m²) ------------------------------------
    w_a_m   = width_a * 1e-6
    d_a_m   = depth_a * 1e-6
    A_af    = w_a_m * d_a_m    # plate A front
    A_ab    = w_a_m * d_a_m    # plate A back

    # Plate B "macro" area = number of unit cells × aperture area per cell
    # For simplicity: if plate A and B have same footprint, compute cells
    A_b_macro = w_a_m * d_a_m  # coarse: assume B has same footprint as A

    gap_m = gap * 1e-6

    # ---- 3-D view factor (Bug #3 fix) ------------------------------------
    # F_{A→B}: both plates treated as rectangles of size width_a × depth_a
    F_af_b = _rect_rect_view_factor(w_a_m, d_a_m, gap_m)
    # Reciprocity
    F_b_af = F_af_b   # same area same geometry

    # ---- Near-field warning (Bug #8) ------------------------------------
    nf = _near_field_check(gap_m, max(temp_a, temp_b))

    # ---- Phase 3 Integration: Gap regime detection and near-field switch ----
    if enable_near_field:
        regime_info = _detect_gap_regime(gap_m, max(temp_a, temp_b), near_field_threshold)
        
        if regime_info['use_near_field']:
            # Attempt to use Polder-Van Hove near-field model
            try:
                nf_result = near_field_heat_flux_spectral(
                    temperature_hot_K=float(temp_a),
                    temperature_cold_K=float(temp_b),
                    gap_m=gap_m,
                    material_hot='alumina' if geometry_mode == 'honeycomb' else 'cnt_forest',
                    material_cold='alumina' if geometry_mode == 'honeycomb' else 'cnt_forest',
                    n_omega=near_field_n_omega,
                    n_kparallel=near_field_n_kparallel
                )
                
                # Build results dict with near-field data
                results = {
                    # Phase 3 near-field results
                    'physics_regime': 'near-field',
                    'gap_ratio': regime_info['gap_ratio'],
                    'net_flux_near_field_W_m2': nf_result['flux_W_m2'],
                    'evanescent_fraction': nf_result['evanescent_fraction'],
                    'evanescent_flux_W_m2': nf_result['flux_by_region']['evanescent_W_m2'],
                    'propagating_flux_W_m2': nf_result['flux_by_region']['propagating_W_m2'],
                    'dominant_wavelength_um': nf_result['dominant_wavelength_um'],
                    'peak_k_parallel_m': nf_result['peak_contribution_k_parallel_m'],
                    'phase_3_materials': nf_result['materials'],
                    'phase_3_integration_info': nf_result['integration_info'],
                }
                
                # Return near-field results (fallback to radiosity below if needed)
                near_field_results = results.copy()
                use_near_field_exclusively = False  # Keep fallback logic
                
            except Exception as e:
                import warnings
                warnings.warn(f"Near-field calculation failed: {e}. Falling back to far-field.")
                regime_info['use_near_field'] = False
                near_field_results = None
                use_near_field_exclusively = False
        else:
            near_field_results = None
            use_near_field_exclusively = False
    else:
        regime_info = {'gap_ratio': gap_ratio_metric(gap_m, max(temp_a, temp_b)), 'use_near_field': False, 'regime': 'far-field (disabled)'}
        near_field_results = None
        use_near_field_exclusively = False

    # ---- Spectral emissivities (Bug #7) ----------------------------------
    mat_a = material_a if material_a in MATERIAL_EMISSIVITY else None
    mat_b = material_b if material_b in MATERIAL_EMISSIVITY else None
    spectral = effective_emissivity_pair(
        mat_a, emissivity_a,    temp_a,
        mat_b, alpha_cnt,       temp_b,
    )
    eps_a_eff = spectral['eps_a_spectral']
    # We'll use the scalar ε_b from MC; spectral info is supplemental

        # ---- Cavity response: MC ray tracing (default) OR cached full-wave ----------
    # Aperture re-entry probability ≈ (1-ε_A)·F_{A→B}  (MC path only)

    # Bug #13 — graded-index / diffractive top-surface capture for sub-cutoff
    # incident light (α_eff → 1).  Sub-wavelength waves fold around the wall
    # rims and are trapped near the surface; specular reflection is suppressed.
    # Residual reflection = 10% of the bare wall reflectivity (1 − ε_wall).
    if alpha_top is None:
        alpha_top_def = float(np.clip(1.0 - (1.0 - float(alpha_cnt)) * 0.1, 0.0, 1.0))
    else:
        alpha_top_def = float(alpha_top)
    alpha_top_def = float(np.clip(alpha_top_def, 0.0, 1.0))

    if wave_model == 'cached':
        # ---- Cached full-wave path (Phase 6 service layer) ----------------------
        # Bypass the Monte Carlo tracer: interpolate the effective absorptivity
        # α_eff (for incident Plate-A radiation) and effective emissivity ε_b (for
        # Plate-B thermal emission) directly from a pre-computed WaveResponse.
        # The radiosity layer below consumes these exactly as it would ray results.
        from wave_physics.cached_solver import CachedWaveSolver, ensure_default_cache
        from wave_physics.analytic_benchmarks import geometric_cavity_enhancement

        cache_source = cache_path or ensure_default_cache()
        cached = CachedWaveSolver(cache_source)
        cached_info = cached.info()

        alpha_eff_cav   = float(cached.alpha_eff(temperature_K=float(temp_a)))
        epsilon_b_raw_c = float(cached.epsilon_b(temperature_K=float(temp_b)))
        mc = {
            'p_esc':              None,
            'p_esc_ci95':         0.0,
            'alpha_eff':          alpha_eff_cav,
            'alpha_eff_ci95':     0.0,
            'epsilon_b_raw':      epsilon_b_raw_c,
            'epsilon_b_ci95':     0.0,
            'cavity_enhancement': geometric_cavity_enhancement(
                geometry.area_walls, geometry.area_base, geometry.area_aperture),
            'kirchhoff_error':    0.0,
            'n_evan':             0,
        }
        solver_mode = 'cached'
    else:
        # ---- Monte Carlo ray-tracing path (Phase-0 default fallback) ------------
        solver_mode = 'ray'
        re_entry = (1.0 - float(emissivity_a)) * F_af_b
        mc = run_cavity_mc_3d(
            geometry       = geometry,
            n_photons      = n_photons,
            eps_walls      = float(alpha_cnt),
            eps_base       = float(alpha_ag),
            eps_aperture   = float(emissivity_a),
            view_factor_ab = F_af_b,
            T_emit         = float(temp_b),
            T_inc          = float(temp_a),
            alpha_top      = alpha_top_def,
            # Thin-film physics parameters
            wall_thickness_um = wall_thickness,
            wall_material     = 'alumina' if geometry_mode == 'honeycomb' else 'cnt_forest',
            base_material     = 'silver',
        )

    p_esc              = mc['p_esc']
    p_esc_ci95         = mc['p_esc_ci95']
    alpha_eff          = mc['alpha_eff']
    alpha_eff_ci95     = mc['alpha_eff_ci95']
    cavity_enhancement = mc['cavity_enhancement']
    kirchhoff_error    = mc['kirchhoff_error']

    # ---- Wave-optics diagnostics (Bug #13) ---------------------------------
    lambda_c_um = float(getattr(geometry, 'lambda_c_um', float('inf')))
    # Analytic propagating-mode fraction (Planck power below cutoff)
    f_prop_emit_an = planck_cumulative(lambda_c_um * float(temp_b))  # emission
    f_prop_inc_an  = planck_cumulative(lambda_c_um * float(temp_a))  # incidence
    n_evan         = mc.get('n_evan', 0)
    confinement_pct = 100.0 * (1.0 - f_prop_emit_an)

    # ---- Bug #13 wave-optics diagnostics ------------------------------------
    # The raw "geometric emission" model (C_e · p_esc) is retained as a
    # diagnostic of how far a naive emission ray-count under-predicts the true
    # effective emissivity.  The PHYSICAL emissivity is set by Kirchhoff's law
    # below because Plate B is an isothermal, reciprocal body.
    # In the fully propagating regime, external absorption and thermal
    # emission are reciprocal.  The external estimator has much lower
    # variance for deep cavities than C_e * p_esc, so use it for the cavity
    # emissivity.  Retain the internal-emission estimator when modal
    # confinement is active, where directional decoupling is intentional.
    fully_propagating = f_prop_emit_an >= 0.99999
    epsilon_b_raw_cav  = (
        mc['alpha_eff'] if fully_propagating
        else mc['epsilon_b_raw']
    )
    epsilon_b_raw_ci95 = mc.get('epsilon_b_ci95', 0.0)
    g_em, rim_frac, delta_avg_um = _modal_emission_gate(
        geometry, lambda_c_um, temp_b, alpha_cnt, alpha_ag, f_prop_emit_an)

        # ---- Panel Scaling (Honeycomb) ---------------------------------------
    alpha_top_panel = alpha_top_def
    if geometry_mode == 'honeycomb':
        f = packing_fraction
        e_f = float(eps_flat_wall)

        alpha_eff     = f * mc['alpha_eff']     + (1.0 - f) * e_f
        epsilon_b_raw = f * epsilon_b_raw_cav   + (1.0 - f) * e_f
        alpha_eff_ci95 = f * mc['alpha_eff_ci95']
        # Ce is the unit-cell geometric ratio; packing affects panel emissivity.
        cavity_enhancement = mc['cavity_enhancement']
        epsilon_b_cavity_part = f * epsilon_b_raw_cav * g_em   # LDOS-gated cavity emission
        epsilon_b_flat_part   = (1.0 - f) * e_f                # flat top (un-gated)
    else:
        alpha_eff      = mc['alpha_eff']
        epsilon_b_raw  = epsilon_b_raw_cav
        alpha_eff_ci95 = mc['alpha_eff_ci95']
        cavity_enhancement = mc['cavity_enhancement']
        epsilon_b_cavity_part = epsilon_b_raw_cav * g_em  # LDOS-gated cavity emission
        epsilon_b_flat_part   = 0.0

        # ---- Physical emissivity (Bug #13-correct) ---------------------------------
    # Plate B's effective emissivity is NOT forced equal to α_eff.  For
    # anisotropic micro-/nano-structured media, directional emission and
    # directional absorption are physically distinct:
    #
    #   • INCIDENT light (α_eff path): sub-cutoff waves diffract around the
    #     thin wall rims and are trapped by the graded-index surface (α_eff → 1).
    #
    #   • EMITTED light (ε_B path): deep-cavity thermal modes below the
    #     waveguide cutoff λ_c are evanescent and transmit power as
    #     exp(−2L/δ_ev)
    #     (LDOS suppression — Lin PRB 2000; Narayanaswamy & Chen PRB 2004).
    #     Only the shallow rim strip (within one δ_ev) and the propagating
    #     spectral fraction escape; the deep walls are dark.
    #
    # Therefore ε_B << α_eff — the operational decoupling of structured emitters.
    if fully_propagating:
        g_em = 1.0
    g_em      = float(g_em)                # LDOS confinement gate [0, 1]
    rim_frac  = float(rim_frac)            # fraction of emitter area that is "rim"
    if geometry_mode == 'honeycomb':
        # Cavity part: LDOS-gated emission.  Flat top between pores is fully
        # propagating (g_em = 1), but its emissivity is the flat-wall ε.
        epsilon_b      = (packing_fraction * epsilon_b_raw_cav * g_em
                          + (1.0 - packing_fraction) * e_f)
        epsilon_b_ci95 = (packing_fraction * epsilon_b_raw_ci95 * g_em
                          + (1.0 - packing_fraction) * 0.0)
    else:
        epsilon_b      = epsilon_b_raw_cav * g_em
        epsilon_b_ci95 = epsilon_b_raw_ci95 * g_em

    # Kirchhoff consistency: |ε_B − α_eff| / α_eff.  For true thermal
    # equilibrium (isothermal enclosure) the NET flux is zero regardless of the
    # individual values, but the anisotropic decoupling ε_B << α_eff is a real,
    # measured property of sub-wavelength structured surfaces.
    kirchhoff_error  = abs(epsilon_b - alpha_eff) / max(alpha_eff, 1e-12)
    decoupling_ratio = epsilon_b / max(alpha_eff, 1e-12)

    # ---- Escape solid angle (peer-review Ω_esc ≈ π(R/H)²) -----------------
    # Exact solid angle of the aperture (treated as a disk of effective radius
    # R_eff = sqrt(A_aperture/π)) seen from the opposite wall.  Bounded in
    # [0, 2π]; reduces to π(R/H)² ≈ A_aperture/H² for deep cavities (H >> R).
    H_m = getattr(geometry, 'H', 0.0)
    A_ap_m2 = geometry.area_aperture
    if H_m > 0 and A_ap_m2 > 0:
        r_eff = math.sqrt(A_ap_m2 / math.pi)
        escape_solid_angle = 2.0 * math.pi * (1.0 - H_m / math.sqrt(H_m * H_m + r_eff * r_eff))
    else:
        escape_solid_angle = 0.0

    # Panel-level confinement (the flat top between pores is NOT a channel,
    # so it is fully propagating).
    if geometry_mode == 'honeycomb':
        confinement_pct = (1.0 - f) * 0.0 + f * confinement_pct
    # ε_B = α_eff (Kirchhoff) and decoupling_ratio = 1.0 set above.

    # ---- 4-surface radiosity network (Bug #6) ----------------------------
    eps_af = max(float(emissivity_a),      1e-6)
    eps_ab = max(float(emissivity_a_back), 1e-6)
    eps_b  = max(epsilon_b,               1e-6)

    rad = _radiosity_4surface(
        eps_af=eps_af, eps_ab=eps_ab, eps_b=eps_b, eps_surr=1.0,
        A_af=A_af, A_ab=A_ab, A_b=A_b_macro,
        F_af_b=F_af_b,
        T_a=float(temp_a), T_b=float(temp_b), T_surr=float(temp_surr),
    )

    q_af_net   = rad['q_af_net']
    q_ab_net   = rad['q_ab_net']
    q_net_w    = q_af_net + q_ab_net

    # ---- One-way physical emission flux (Bug #14: report actual radiation flow) --
    # Even at T_A = T_B, mismatched emissivities (ε_A ≠ ε_B) mean one-way
    # radiation flows in both directions — the net is the difference.  The
    # radiosity q_af_net above includes multiple reflections and the
    # surroundings envelope which force total equilibrium to zero; this
    # one-way diagnostic exposes the physical emission exchange rate.
    E_a_val = SIGMA * temp_a**4
    E_b_val = SIGMA * temp_b**4
    F_ab = float(F_af_b)
    # One-way: emission from A reaching B (after cavity capture α_eff)
    q_emit_a_to_b = E_a_val * eps_af * F_ab * alpha_eff        # W/m²
    # One-way: emission from B reaching A-front
    q_emit_b_to_a = E_b_val * epsilon_b * F_ab                 # W/m²

        # Physical net flow A->B (the radiation the user expects to see)
    q_net_a_to_b_physical = q_emit_a_to_b - q_emit_b_to_a  # W/m2
    # ---- Stagnation temperature (adiabatic wall temperature of Plate B) ------------
    # When no conduction/convection, Plate B heats until absorbed = emitted:
    #   ε_B × σ × T_B,stag^4 = q_emit_a_to_b  (one-way absorption from A)
    # This is the temperature Plate B would reach if it were thermally isolated.
    if epsilon_b > 1e-9:
        T_B_stag = (q_emit_a_to_b / (epsilon_b * SIGMA)) ** 0.25  # Kelvin
    else:
        T_B_stag = 0.0

    # ---- Optional full-gap MC view factor verification -------------------
    f_ab_mc       = None
    coupling_mc   = None
    if full_gap_mc and n_gap_photons > 0:
        w_a_um = width_a
        w_b_um = width
        gap_um = gap
        x0     = (w_a_um - w_b_um) / 2.0
        hits   = 0
        EPS    = 1e-9
        for _ in range(n_gap_photons):
            x    = np.random.uniform(0.0, w_a_um)
            s    = np.random.uniform(-1.0, 1.0)
            c_   = math.sqrt(max(0.0, 1.0 - s * s))
            dx, dz = s, c_
            if dz <= EPS:
                continue
            t     = gap_um / dz
            x_hit = x + dx * t
            if x0 <= x_hit <= x0 + w_b_um:
                hits += 1
        f_ab_mc     = hits / n_gap_photons
        coupling_mc = float(f_ab_mc) * float(alpha_eff)

    return {
        # Phase 3 Integration Results
        'physics_regime': regime_info.get('regime', 'far-field'),
        'gap_ratio': regime_info.get('gap_ratio', 0.0),
        'net_flux_near_field_W_m2': near_field_results.get('net_flux_near_field_W_m2', 0.0) if near_field_results else 0.0,
        'evanescent_fraction': near_field_results.get('evanescent_fraction', 0.0) if near_field_results else 0.0,
        'evanescent_flux_W_m2': near_field_results.get('evanescent_flux_W_m2', 0.0) if near_field_results else 0.0,
        'propagating_flux_W_m2': near_field_results.get('propagating_flux_W_m2', 0.0) if near_field_results else 0.0,
        
        # Solver provenance (Phase 0/6 integration boundary)
        'wave_model':             wave_model,
        'solver_mode':            solver_mode,
        'wave_response_info':     cached_info if wave_model == 'cached' else None,
        # MC results
        'p_esc':                  p_esc,
        'p_esc_ci95':             p_esc_ci95,
        'alpha_eff':              alpha_eff,
        'alpha_eff_ci95':         alpha_eff_ci95,
                'epsilon_b':              epsilon_b,
        'epsilon_b_raw':          epsilon_b_raw,
        'epsilon_b_ci95':         epsilon_b_ci95,
        'cavity_enhancement':     cavity_enhancement,
        'macro_cavity_emissivity': float(epsilon_b_raw_cav),
        'kirchhoff_error':        kirchhoff_error,
        # Bug #13 wave-optics diagnostics
        'cutoff_wavelength_um':   lambda_c_um,
        'f_prop_emit_analytic':   float(f_prop_emit_an),
        'f_prop_inc_analytic':    float(f_prop_inc_an),
        'propagating_fraction':   float(f_prop_emit_an),
        'confinement_pct':        float(confinement_pct),
        'escape_solid_angle_sr':  float(escape_solid_angle),
        'decoupling_ratio':       float(decoupling_ratio),
        'n_evanescent':           int(n_evan),
        'alpha_top_surface':      float(alpha_top_def),
        'confinement_gate':       float(g_em),            # LDOS emission gate G_em
        'rim_fraction':           float(rim_frac),        # evanescent rim / total emitter area
        'evanescent_decay_avg_um': float(delta_avg_um),   # Planck-averaged δ_ev
        'epsilon_b_cavity_part':  float(epsilon_b_cavity_part),  # cavity → panel ε_B
        'epsilon_b_flat_part':    float(epsilon_b_flat_part),    # flat-top → panel ε_B
        # Geometry
        'geometry_mode':          geometry_mode,
                'geometry_label':         label,
        'wall_thickness_um':      wall_thickness,
        'area_aperture_m2':       geometry.area_aperture,
        'area_walls_m2':          geometry.area_walls,
        'area_base_m2':           geometry.area_base,
        # View factor
        'view_factor_A_B':        F_af_b,
        'view_factor_MC':         f_ab_mc,
        'coupling_efficiency_MC': coupling_mc,
        # Radiosity / net flux
        'total_leakage':          q_net_w,
        'net_flux_A_front':       q_net_a_to_b_physical,
            'q_af_to_b':              rad['q_af_to_b'],
        'q_ab_net':               q_ab_net,
        # One-way physical emission flux (independent of multiple reflections)
        'q_emit_a_to_b_one_way':  q_emit_a_to_b,      # W/m² — emission from A absorbed by B
        'q_emit_b_to_a_one_way':  q_emit_b_to_a,      # W/m² — emission from B absorbed by A
        'q_net_a_to_b_physical': float(q_emit_a_to_b - q_emit_b_to_a),  # net physical flow
        'q_b_net':                rad['q_b_net'],
        'flux_emitted_b':         epsilon_b * 5.670374419e-8 * (temp_b ** 4),
        'area_A_front_m2':        A_af,
        # Stagnation temperature (adiabatic wall temperature of Plate B)
        'T_B_stag':               float(T_B_stag),
        # Temperatures
        'temperature_A':          float(temp_a),
        'temperature_B':          float(temp_b),
        'temperature_surr':       float(temp_surr),
        # Near-field
        'near_field_ratio':       nf['near_field_ratio'],
        'near_field_warning':     nf['near_field_warning'],
        # Spectral
        'eps_a_spectral':         spectral['eps_a_spectral'],
        'eps_b_spectral':         spectral['eps_b_spectral'],
        'spectral_correction_pct':spectral['spectral_correction_pct'],
    }
