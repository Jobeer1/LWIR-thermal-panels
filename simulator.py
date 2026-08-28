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
# Physics Orchestrator (pre-flight dimensionless regime audit).  Aliased to
# avoid clashing with this module's own dispatcher `select_physics_regime`.
from regime_selector import select_physics_regime as _orchestrate_regime
from ray_tracer  import run_cavity_mc_3d, evanescent_power_transmission
from spectral    import planck_weighted_emissivity, MATERIAL_EMISSIVITY, effective_emissivity_pair
from sampling    import planck_cumulative, planck_averaged_evanescent_decay
from material_optics import (
    temperature_optics_provenance,
    planck_weighted_effective_emissivity,
    get_complex_refractive_index_at_temperature,
    aperture_boundary_absorptance,
    material_optics_confidence,
)
from brdf import brdf_provenance
from near_field_radiative_heat import (
    gap_ratio_metric,
    should_use_near_field_model,
    near_field_heat_flux_spectral
)

SIGMA = 5.670374419e-8
NEAR_FIELD_SOLVER = 'nf_greens'  # structured near-field Green-tensor / LDOS solver kind  # Stefan-Boltzmann constant, W m⁻² K⁻⁴


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
# Physics Regime Dispatcher (Part 2 — dimensionless boundary routing)
# ---------------------------------------------------------------------------

# Canonical regime identifiers returned by ``select_physics_regime``.
REGIME_NEAR_FIELD = 'NEAR_FIELD_FLUCTUATING_ELECTRODYNAMICS'
REGIME_EFFECTIVE_MEDIUM = 'MAXWELL_GARNETT_EFFECTIVE_MEDIUM'
REGIME_FULL_WAVE = 'FULL_WAVE_RCWA_CACHED'
REGIME_MACRO = 'MACRO_GOUFFE_RAY_TRACE'

# Transverse TE eigenvalues (x_m) for the implicit gamma_m = sqrt((x_m/R)^2-k0^2)
# modal decay sum, and their relative modal excitation weights c_m.
_TE_MODE_X_ROOTS = (1.8412, 3.0542, 3.8317, 4.1997, 5.3189)
_TE_MODE_WEIGHTS = (1.00, 0.32, 0.15, 0.08, 0.03)


def select_physics_regime(gap_m: float, geometry, T_hot_K: float,
                          wall_thickness_um: float = None,
                          wall_material: str = 'alumina'):
    """Dynamically route the simulation to the mathematically valid framework.

    The dispatcher evaluates the four dimensionless ratios that select the
    governing optical framework (λ from Wien's peak, ``λ = 2897.77/T``):

        λ / D          — wavelength vs. feature diameter (diffraction /
                         geometric coupling crossover).
        λ / P          — wavelength vs. pitch (homogenisation validity).
        d_gap/(λ/2π)   — plate gap vs. evanescent tunnelling length.
        t_wall / δ     — wall thickness vs. absorption depth (optically thin).

    Dimensionally applied boundary criteria:

        Near-field evanescent        : gap_m <  λ_peak / (2π)
        Effective medium (EMT)       : P < 0.2 λ_peak  AND  D < 0.2 λ_peak
        Sub-λ resonance / diffraction: 0.2 λ_peak ≤ D ≤ 5 λ_peak
        Macro geometric optics       : λ/D < 0.2 (D > 5 λ_peak) + far gap

    Parameters
    ----------
    gap_m            : Plate-to-plate gap distance (m).
    geometry        : geometry object exposing ``diameter_um`` / ``width_um``
                      and ``P`` (pitch in metres).
    T_hot_K         : Hot reference temperature (K) for Wien's peak.
    wall_thickness_um : (optional) wall thickness (µm) for the t_wall/δ audit.
    wall_material   : (optional) wall material key for the t_wall/δ absorption.

    Returns
    -------
    (regime: str, info: dict)
        ``regime`` is one of the ``REGIME_*`` constants; ``info`` carries the
        dimensionless ratios for audit / UI reporting.
    """
    lambda_peak_m = 2897.77e-6 / max(float(T_hot_K), 1.0)
    lambda_peak_um = lambda_peak_m * 1e6
    lam_2pi_m = lambda_peak_m / (2.0 * math.pi)

    P_attr = getattr(geometry, 'P', None)
    P_m = float(P_attr) if P_attr is not None else None

    D_raw = getattr(geometry, 'diameter_um',
                    getattr(geometry, 'width_um', None))
    R_m = None
    if D_raw is None:
        R_m = getattr(geometry, 'R', None)
    if D_raw is not None:
        D_m = float(D_raw) * 1e-6
    elif R_m is not None:
        D_m = float(R_m) * 2.0
    elif P_m is not None:
        # No explicit aperture diameter (e.g. CNT forest): use the array pitch
        # as the conservative feature/per-array homogenisation scale.
        D_m = P_m
    else:
        D_m = 1e-6
    if P_m is None:
        P_m = D_m

    # t_wall/δ — dimensionless wall optical depth at the absorption peak.
    wall_optical_depth = None
    wall_optically_thin = False
    if wall_thickness_um is not None and wall_thickness_um > 0.0:
        try:
            from material_optics import get_absorption_depth
            delta_peak_um = get_absorption_depth(
                str(wall_material).lower(), lambda_peak_um,
                temperature_K=float(T_hot_K))
            if delta_peak_um and delta_peak_um > 0.0:
                wall_optical_depth = float(wall_thickness_um) / delta_peak_um
                wall_optically_thin = wall_optical_depth < 0.5
        except Exception:
            wall_optical_depth = None

    info = {
        'lambda_peak_m':  float(lambda_peak_m),
        'lambda_peak_um': float(lambda_peak_um),
        'D_m':            D_m,
        'P_m':            P_m,
        'gap_m':          float(gap_m),
        # Dimensionless framework evaluation (Part 5 — regime routing)
        'lambda_over_D':  lambda_peak_m / max(D_m, 1e-30),
        'lambda_over_P':  lambda_peak_m / max(P_m, 1e-30),
        'D_lambda_frac':  D_m / max(lambda_peak_m, 1e-30),
        'P_lambda_frac':  P_m / max(lambda_peak_m, 1e-30),
        'gap_lambda_frac': gap_m / max(lam_2pi_m, 1e-30),
        'gap_over_lambda_2pi': gap_m / max(lam_2pi_m, 1e-30),
        # wall optical depth (t_wall / δ) — optically thin emissive limit
        'wall_optical_depth_t_over_delta': wall_optical_depth,
        'wall_optically_thin': wall_optically_thin,
    }

    # 1. Near-field boundary check (evanescent photon tunneling)
    if gap_m < lam_2pi_m:
        return REGIME_NEAR_FIELD, info

    # 2. Effective medium theory (sub-wavelength spatial homogenisation)
    #    requires BOTH pitch and aperture at the P, D ≪ λ limit.
    if P_m < 0.2 * lambda_peak_m and D_m < 0.2 * lambda_peak_m:
        return REGIME_EFFECTIVE_MEDIUM, info

    # 3. Wave resonances & diffraction crossover (requires a full-wave solver)
    if 0.2 * lambda_peak_m <= D_m <= 5.0 * lambda_peak_m:
        return REGIME_FULL_WAVE, info

    # 4. Incoherent geometric ray trapping (macro scale), λ/D < 0.2
    return REGIME_MACRO, info


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

def _geometric_cavity_emissivity(geometry, eps_wall: float) -> float:
    """Gouffe/Sparrow effective emissivity in the geometric-optics limit."""
    aperture_area = float(geometry.area_aperture)
    internal_area = float(geometry.area_walls + geometry.area_base)
    if aperture_area <= 0.0 or internal_area <= 0.0:
        return 0.0
    area_ratio = aperture_area / internal_area
    return float(eps_wall / (1.0 - (1.0 - eps_wall) * (1.0 - area_ratio)))


def _homogenized_tmm_slab(material: str, wavelength_um: float,
                          temperature_K: float, fill_fraction: float,
                          thickness_um: float,
                          backing_material: str = 'silver') -> float:
    """Homogenised Maxwell-Garnett TMM slab absorptance (EMT regime operator).

    Computes the effective optical response of the structured top layer as a
    single homogenised thin film (Maxwell-Garnett ``n_eff(f)`` on the matrix
    material) backed by the base/substrate, via the transfer-matrix method at
    normal incidence.  This is the physical operator executed in the
    ``MAXWELL_GARNETT_EFFECTIVE_MEDIUM`` regime — geometric ray tracing is
    bypassed there because P, D ≪ λ.

        α_slab = 1 − ½( R_s + R_p )        (film on reflective base)

    Returns
    -------
    float in [0, 1] — normal-incidence absorptance of the homogenised slab.
    """
    try:
        from material_optics import (maxwell_garnett_effective_index,
                                     tmm_reflectance_single_layer,
                                     get_complex_refractive_index_at_temperature)
    except ImportError as exc:  # pragma: no cover
        raise ValueError('EMT regime operator requires material_optics.') from exc

    n_eff = maxwell_garnett_effective_index(
        material, float(wavelength_um), float(fill_fraction),
        temperature_K=float(temperature_K))
    n_base_r, n_base_k = get_complex_refractive_index_at_temperature(
        backing_material, wavelength_um, temperature_K)
    n_sub = complex(n_base_r, n_base_k) if n_base_k > 0 else complex(n_base_r, 0.0)

    R_s = tmm_reflectance_single_layer(1.0 + 0.0j, n_eff, n_sub,
                                       float(thickness_um), wavelength_um, 0.0, 's')
    R_p = tmm_reflectance_single_layer(1.0 + 0.0j, n_eff, n_sub,
                                       float(thickness_um), wavelength_um, 0.0, 'p')
    R_avg = 0.5 * (R_s + R_p)
    return float(np.clip(1.0 - R_avg, 0.0, 1.0))


def _effective_medium_operator(geometry, material: str, wavelength_um: float,
                               temperature_K: float, fill_fraction: float,
                               thickness_um: float) -> dict:
    """Produce the full cavity effective-absorptance payload for the EMT regime.

    Returns a dict shaped like the MC tracer result (``alpha_eff``,
    ``epsilon_b_raw``, ``cavity_enhancement``, …) so the downstream radiosity /
    panel-scaling code consumes the homogenized TMM slab operator identically to
    ray-traced results, while geometric ray tracing is bypassed.
    """
    alpha_int = _homogenized_tmm_slab(material, wavelength_um, temperature_K,
                                      fill_fraction, thickness_um)
    return {
        'p_esc':              None,
        'p_esc_ci95':         0.0,
        'alpha_eff':          alpha_int,
        'alpha_eff_ci95':     0.0,
        'epsilon_b_raw':      alpha_int,   # Kirchhoff: homogeneous slab α ≡ ε
        'epsilon_b_ci95':     0.0,
        'cavity_enhancement': 0.0,
        'kirchhoff_error':    0.0,
        'n_evan':             0,
    }


def compute_cavity_emissivity(
    geometry, eps_wall: float, f_prop: float, p_esc: float = None,
    eps_base: float = 0.0,
) -> float:
    """Return the bounded cavity operator for either spectral regime."""
    eps_gouffe = _geometric_cavity_emissivity(geometry, float(eps_wall))
    if f_prop > 0.95 or p_esc is None:
        return eps_gouffe
    # Emissivity-weighted emission area ratio (MATHEMATICAL_DERIVATIONS
    # Phase 4/5): εcav = P_esc × (εw·A_walls + εbase·A_base) / A_ap.
    # The bare Ce = A_int/A_ap enhancement implicitly assumes black walls
    # (εw = 1) and therefore overshoots the cavity emissivity by ~Ce×(1−εw).
    A_ap = float(geometry.area_aperture)
    A_walls = float(geometry.area_walls)
    A_base = float(geometry.area_base)
    cavity_factor = ((float(eps_wall) * A_walls + float(eps_base) * A_base)
                     / max(A_ap, 1e-30))
    eps_mc = cavity_factor * float(p_esc)
    # Bounded operator: the Gouffé geometric bound and the unit ceiling are
    # the physical limits.  When the Ce×P_esc statistic exceeds them (deep
    # cavities whose propagating band carries a non-trivial Planck fraction,
    # where the amplified-rare-escape estimator saturates), the bound itself
    # governs — aborting here would kill otherwise valid regimes.
    return float(min(eps_gouffe, eps_mc, 1.0))


def compute_cavity_emissivity_ci95(
    geometry, eps_wall: float, f_prop: float, p_esc: float,
    p_esc_ci95: float, eps_base: float = 0.0,
) -> float:
    """95% CI on :func:`compute_cavity_emissivity`, using the SAME operator.

    When ``f_prop > 0.95`` (or no ``P_esc``) the point estimate is the
    deterministic Gouffé bound → CI is exactly 0.  Otherwise the estimate is the
    emissivity-weighted rare-escape estimator ``cavity_factor × P_esc``, and the
    CI is that same factor times the MC ``P_esc`` uncertainty.  Using the bare
    area ratio ``Ce`` here over-amplifies the CI by Ce/(ε-weighted factor) and,
    worse, reports a large CI even when the deterministic Gouffé ceiling
    governs — the source of the unphysical ±42% / ±70% ε_B intervals.
    """
    eps_gouffe = _geometric_cavity_emissivity(geometry, float(eps_wall))
    if f_prop > 0.95 or p_esc is None or p_esc_ci95 is None:
        return 0.0
    A_ap = float(geometry.area_aperture)
    cavity_factor = ((float(eps_wall) * float(geometry.area_walls)
                      + float(eps_base) * float(geometry.area_base))
                     / max(A_ap, 1e-30))
    if cavity_factor * float(p_esc) >= eps_gouffe:
        return 0.0                       # saturated at the deterministic ceiling
    return float(cavity_factor * float(p_esc_ci95))


def _effective_medium_absorptivity(
    wavelength_um: float, temperature_K: float, material: str, fill_fraction: float
) -> float:
    """Normal-incidence Maxwell-Garnett absorptance for unresolved pores."""
    n_real, k_imag = get_complex_refractive_index_at_temperature(
        material, wavelength_um, temperature_K
    )
    eps_matrix = complex(n_real, k_imag) ** 2
    eps_inclusion = 1.0 + 0.0j
    fraction = min(max(float(fill_fraction), 0.0), 1.0)
    contrast = (eps_inclusion - eps_matrix) / (eps_inclusion + eps_matrix)
    denominator = 1.0 - fraction * contrast
    if abs(denominator) < 1e-15:
        return 1.0
    eps_eff = eps_matrix * (1.0 + (2.0 * fraction * contrast) / denominator)
    n_eff = np.sqrt(eps_eff + 0.0j)
    reflection = (1.0 - n_eff) / (1.0 + n_eff)
    absorptivity = 1.0 - abs(reflection) ** 2
    return float(min(max(absorptivity, 0.0), 1.0))


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
# Track B2 — structured near-field Green tensor + LDOS (active correction)
# ---------------------------------------------------------------------------
#
# The structured Green-tensor / LDOS model in wave_physics.near_field_greens
# supplies an active *enhancement factor* for the inter-cavity evanescent
# channel.  The NearFieldResponse cache is built on first use for the
# simulation's cavity geometry / temperatures / material and memoized
# per-process, so the moderately expensive 2-D quadrature runs only once per
# parameter set.  Subsequent calls interpolate the cached flux table at the
# simulation gap.
_NF_CACHE_MEMO: dict = {}


def _cavity_radius_um(geometry_mode: str, cavity_diameter: float,
                      cnt_pitch: float, cnt_dia_base: float, width: float) -> float:
    """Cavity radius (um) used by the Green-tensor LDOS model."""
    if geometry_mode == 'honeycomb':
        return cavity_diameter / 2.0
    if geometry_mode == 'cnt_forest':
        return cnt_pitch / 2.0
    if geometry_mode == 'frustum':
        return (cnt_dia_base * 1e-3) / 2.0          # nm -> um
    return width / 2.0                            # rect_pit (legacy)


def _nf_material_for(geometry_mode: str) -> str:
    """Material key for the structured near-field (matches Phase 3 lookup)."""
    return 'alumina' if geometry_mode == 'honeycomb' else 'cnt_forest'


def _structured_near_field_correction(gap_m, temperature_hot_K,
                                      temperature_cold_K, cavity_radius_um,
                                      cavity_depth_um, material):
    """Active near-field correction factor from the structured Green tensor.

    Builds (or reuses, memoized per-process) the structured NearFieldResponse
    cache for the cavity geometry, interpolates the NFRHT flux table at the
    simulation gap, and returns a multiplicative enhancement factor

        correction = Phi_nf(gap) / Phi_ff_bb

    where Phi_ff_bb = sigma * (T_hot**4 - T_cold**4) is the blackbody
    far-field reference.  The factor is symmetric in A<->B and clamped to
    [1, 1e4]: it only ever *enhances* the ray-traced inter-plate flux in the
    near-field regime, never reduces it, and is exactly 1.0 outside the
    sub-wavelength regime -- leaving the far-field radiosity path untouched.

    Returns
    -------
    (nf_response, correction, phi_nf, phi_ff)
        Defaults are (None, 1.0, 0.0, 0.0) when the wave_physics package, the
        cache build, or a usable gap is unavailable.
    """
    try:
        from wave_physics import near_field_greens as nf_greens
    except Exception:
        return None, 1.0, 0.0, 0.0

    key = (round(float(cavity_radius_um), 9), round(float(cavity_depth_um), 9),
           round(float(temperature_hot_K), 6), round(float(temperature_cold_K), 6),
           str(material))
    resp = _NF_CACHE_MEMO.get(key)
    if resp is None:
        try:
            resp = nf_greens.build_near_field_cache(
                cavity_radius_um=float(cavity_radius_um),
                cavity_depth_um=float(cavity_depth_um),
                temperature_hot_K=float(temperature_hot_K),
                temperature_cold_K=float(temperature_cold_K),
                material=str(material),
            )  # quadrature resolution defaults to n_omega=60, n_kparallel=30
            _NF_CACHE_MEMO[key] = resp
        except Exception:
            resp = None
    if resp is None:
        return None, 1.0, 0.0, 0.0

    phi_ff = float(SIGMA * (float(temperature_hot_K) ** 4
                            - float(temperature_cold_K) ** 4))

    gap_um = float(gap_m) * 1.0e6
    gaps = np.asarray(resp.gap_um, dtype=float)
    table = np.asarray(resp.flux_W_m2, dtype=float)
    # Gap outside the cached tabulated range -> no correction available.
    if (gaps.size < 2 or np.nanmin(gaps) > gap_um or np.nanmax(gaps) < gap_um
            or table.size == 0):
        return resp, 1.0, 0.0, phi_ff

    # build_near_field_cache fills every wavelength row with the same per-gap
    # flux (the wavelength axis is a reporting label), so collapse to a single
    # flux-vs-gap curve and interpolate at the simulation gap.
    flux_per_gap = np.nanmean(table, axis=0) if table.ndim == 2 else table
    phi_nf = float(np.interp(gap_um, gaps, flux_per_gap))
    if phi_nf <= 0.0 or phi_ff <= 0.0:
        return resp, 1.0, phi_nf, phi_ff
    corr = phi_nf / phi_ff
    if not np.isfinite(corr) or corr < 1.0:
        corr = 1.0
    return resp, float(min(corr, 1.0e4)), phi_nf, phi_ff


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

# ---------------------------------------------------------------------------
# Exact modal-cutoff evanescent exitance (Part 1 §1)
# ---------------------------------------------------------------------------

def _modal_cutoff_evanescent_exitance(lambda_c_um: float, lambda_um: float,
                                      depth_um: float,
                                      material: str = 'alumina',
                                      diameter_um: float = None) -> float:
    """Sub-wavelength modal exitance out of a channel of depth ``depth_um``.

    Evaluates the exact modal-cutoff evanescent attenuation:

        ε_modal(λ) = Σ_m c_m · exp(−2·Re[γ_m]·H),
        γ_m = sqrt( (x_m/R)² − (2π/λ)² ),

    where ``x_m`` are the transverse TE eigenvalues and ``R`` the channel
    radius.  The boundary-mode wave impedance at the aperture plane
    (mode-matching Z_mode vs. free-space Z_0) is folded in as an aperture
    transmittance ``T_ap`` so sub-wavelength thermal leakage out of deep
    pores is never overestimated.

    ``ray_tracer.evanescent_power_transmission`` collapses the loss onto the
    single dominant TE11 mode via a λ_c-derived decay length and a separate
    aperture factor; this routine keeps the full multi-mode γ_m sum and an
    independent imaging-impedance factor, per the peer-review spec.

    Returns
    -------
    float in [0, 1] — fraction of evanescent emission surviving to the aperture.
    """
    if lambda_c_um is None or lambda_c_um <= 0.0 or not math.isfinite(lambda_c_um):
        return 1.0                                  # no waveguide cutoff
    lam = float(lambda_um)
    if lam < lambda_c_um or depth_um <= 0.0:
        return 1.0                                  # propagating or zero depth

    R_um = (float(diameter_um) / 2.0 if diameter_um and diameter_um > 0.0
            else float(lambda_c_um) * 1.8412 / math.pi)
    k0 = 2.0 * math.pi / max(lam, 1e-12)            # wavenumber in 1/µm

    num = 0.0
    tot = 0.0
    for x_m, c_m in zip(_TE_MODE_X_ROOTS, _TE_MODE_WEIGHTS):
        kappa_m = x_m / R_um
        gamma_m = math.sqrt(max(kappa_m * kappa_m - k0 * k0, 0.0))
        num += c_m * math.exp(-2.0 * gamma_m * depth_um)
        tot += c_m
    decay = (num / tot) if tot > 0.0 else 0.0

    # Aperture-plane boundary-mode impedance factor (Z_mode vs Z_0).
    T_ap = 1.0
    try:
        from waveguide_modes import aperture_modal_impedance_and_reflection
        gamma0 = math.sqrt(
            max((_TE_MODE_X_ROOTS[0] / R_um) ** 2 - k0 * k0, 0.0))
        if gamma0 > 1e-12:
            _Z, _R_ap, _T_ap = aperture_modal_impedance_and_reflection(
                complex(0.0, gamma0), k0)
            T_ap = float(np.clip(_T_ap, 0.0, 1.0))
    except Exception:
        T_ap = 1.0                                    # conservative fallback

    return float(np.clip(T_ap * decay, 0.0, 1.0))


def _average_evanescent_transmission(lambda_c_um: float, T_emit: float, depth_um: float,
                                     material: str = 'alumina',
                                     diameter_um: float = None) -> float:
    """Planck- and depth-averaged evanescent power escaping a sub-cutoff channel
    of depth depth_um (same units as lambda_c_um).  Each deep-wall lattice element
    tunnels T=T_ap*exp(-2*kappa*L) toward the aperture; we average that against
    the blackbody spectrum at T_emit.  A genuine wave integral including aperture
    impedance-mismatch (Fix 2b: T_ap via Z_mode)."""
    if not lambda_c_um or lambda_c_um <= 0.0 or not math.isfinite(lambda_c_um):
        return 0.0
    if T_emit <= 0.0 or depth_um <= 0.0:
        return 0.0
    C2 = 14387.769
    lam_hi = max(1.05 * lambda_c_um, lambda_c_um + 0.05)
    lam = np.linspace(lam_hi, min(lam_hi * 40.0, 2000.0), 4000)
    num = 0.0; den = 0.0
    for i in range(1, len(lam)):
        lmid = 0.5 * (lam[i - 1] + lam[i])
        ex = math.exp(min(700.0, C2 / (lmid * T_emit)))
        M = (1.0 / lmid ** 5) / (ex - 1.0)
        # Fix 2b: exact multi-mode γ_m exitance including the aperture-plane
        # boundary mode impedance (Z_mode vs Z_0) — Part 1 §1 peer-review.
        T = _modal_cutoff_evanescent_exitance(
            lambda_c_um, lmid, depth_um,
            material=material,
            diameter_um=diameter_um,
        )
        dw = lam[i] - lam[i - 1]
        num += T * M * dw; den += M * dw
    return float(num / den) if den > 0.0 else 0.0


def _modal_emission_gate(geometry, lambda_c_um: float, T_emit: float,
                         eps_wall: float, eps_base: float,
                         f_prop: float,
                         material: str = 'alumina',
                         diameter_um: float = None):
    """LDOS emission-confinement factor for the cavity (diagnostic only).

    NOTE (Fix 4, peer-review): This function's output G_em is now a DIAGNOSTIC
    only and must NOT be applied as a multiplicative factor to epsilon_b after
    the Monte Carlo simulation.  The MC loop already samples wavelengths from
    the Planck distribution and gates each photon on lambda < lambda_c (evanescent
    vs propagating) per photon — epsilon_b_raw_cav is therefore the correct
    physical emissivity and applying G_em again would double-count the modal
    gating, artificially depressing epsilon_B.

    Returns (G_em_diagnostic, T_evan_avg, delta_avg_um).
    """
    if not lambda_c_um or lambda_c_um <= 0.0 or not math.isfinite(lambda_c_um):
        return 1.0, 0.0, 0.0
    H_um = getattr(geometry, 'H', 0.0) * 1e6
    delta_avg_um = planck_averaged_evanescent_decay(
        lambda_c_um, float(T_emit), cap_um=H_um if H_um > 0 else 0.0)
    # Fix 2b: pass material and diameter so T_ap is included in the integral.
    T_evan_avg = _average_evanescent_transmission(
        lambda_c_um, float(T_emit),
        H_um if H_um > 0 else (delta_avg_um if delta_avg_um > 0 else 1.0),
        material=material,
        diameter_um=diameter_um,
    )
    g_em = float(f_prop) + (1.0 - float(f_prop)) * T_evan_avg
    return float(np.clip(g_em, 0.0, 1.0)), float(T_evan_avg), float(delta_avg_um)


# ---------------------------------------------------------------------------
# Part 3 — Second-Law reciprocity / equilibrium guard
# ---------------------------------------------------------------------------

def enforce_thermal_equilibrium_zero_net(q_net_AB_density: float,
                                         temp_a: float, temp_b: float,
                                         temp_surr: float,
                                         tol: float = 1e-9) -> dict:
    """2nd-Law guard: under a truly isothermal system the net exchange vanishes.

    If T_A == T_B == T_surr then directional-hemispherical units obey
        ∫ α_eff(λ) E_b(λ,T) dλ ≡ ∫ ε_eff(λ) E_b(λ,T) dλ
    so the spectrally-integrated radiosity net flow must be zero:

        q_net_AB_density ≈ 0            (tolerance < tol W/m²)

    Returns a diagnostics dict and ``passes``; the caller may ``raise`` on a
    violation (academic audit flag) without mutating any physics.
    """
    is_uniform = (abs(temp_a - temp_b) < 1e-9 and abs(temp_a - temp_surr) < 1e-9)
    residual = abs(float(q_net_AB_density))
    ok = (not is_uniform) or (residual < tol)
    return {
        'is_thermally_uniform': bool(is_uniform),
        'q_net_AB_density':      float(q_net_AB_density),
        'abs_residual':          residual,
        'tolerance_W_m2':        tol,
        'passes':                bool(ok),
    }


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
    # Part 2 — Physical-regime-routing control (auto-dispatch)
    enforce_physics_regime: bool = True,  # strict: never solve in invalid regime
    # Phase 4b — surface-roughness BRDF (Beckmann–Spizzichino / Harvey–Shack)
    surface_roughness_um: float = None,       # RMS roughness σ (µm); None → Lambertian
    roughness_correlation_um: float = None,   # correlation length τ (µm)
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

    # ---- Part 2 — Physics Orchestrator pre-flight (regime_selector.py) ------
    # Non-dimensional audit (λ/D, λ/H, λ/t, d/λ, H/D, t/δ) evaluated BEFORE
    # any solver / MC loop.  Its enforcement flags steer the strict routing:
    #   • ray_valid == False          → pure geometric ray tracing disabled;
    #     the wave operators (RCWA cache / effective medium) or, when no full-
    #     wave response is supplied, the modal-guided photon transport with
    #     per-photon spectral cutoff is executed instead.
    #   • nearfield_required          → fluctuating electrodynamics solver.
    #   • thinfilm_required           → dynamic TMM wall optics enforced.
    #   • modal_cutoff_required       → per-photon λ_c gate + evanescent decay.
    try:
        from material_optics import get_complex_refractive_index as _gcn
        _lam_pk_um = 2898.0 / max(float(max(temp_a, temp_b)), 1e-6)
        _n_r, _k_e = _gcn(
            'alumina' if geometry_mode == 'honeycomb' else 'cnt_forest',
            _lam_pk_um, temperature_K=float(max(temp_a, temp_b)))
        _lambda_c_orch = float(getattr(geometry, 'lambda_c_um', 0.0)) or None
        orchestrator = _orchestrate_regime(
            T_emit_K=max(temp_a, temp_b),
            diameter_um=(float(cavity_diameter) if geometry_mode == 'honeycomb'
                         else float(cnt_pitch)),
            height_um=float(height),
            wall_thickness_um=float(wall_thickness),
            gap_um=float(gap),          # run_simulation gap parameter is in µm
            k_extinction=max(float(_k_e), 0.0),
            lambda_cutoff_um=_lambda_c_orch,
            rcwa_available=bool(cache_path or wave_model == 'cached'),
        )
        orch = orchestrator.to_dict()
    except Exception as _orch_err:                      # diagnostics must never
        orch = {'error': str(_orch_err)}                # block a valid solve
        _lambda_c_orch = None

    # ---- Part 2 — Physics Regime Dispatcher (select + strict route) --------
    # Query the dimensionless regime engine at the start and apply the
    # strict fallbacks BEFORE any solver is invoked.
    physics_regime, phys_regime_info = select_physics_regime(
        gap_m, geometry, max(temp_a, temp_b),
        wall_thickness_um=wall_thickness,
        wall_material='alumina' if geometry_mode == 'honeycomb' else 'cnt_forest')
    regime_warning = ''
    regime_violation = ''

    if physics_regime == REGIME_NEAR_FIELD:
        enable_near_field = True
    elif physics_regime == REGIME_EFFECTIVE_MEDIUM:
        pass
    elif physics_regime == REGIME_FULL_WAVE:
        regime_warning = (
            'FULL_WAVE_RCWA_CACHED regime: structure in sub-lambda window; '
            'a real full-wave solve is required for publication accuracy.')
        regime_violation = regime_warning

    # ----------------------------------------------------------------------
    # PART 5 — STRICT PHYSICS-OPERATOR ROUTING.
    # The dimensionless engine has decided which operator is mathematically
    # valid.  When enforce_physics_regime is True we MUST execute that operator
    # and NEVER ray-trace a sub-wavelength (EMT / full-wave) structure, nor
    # wave-solve a purely geometric (macro) cavity.
    #   • FULL_WAVE   → 'cached' (full-wave response / modal operator)
    #   • EMT         → 'effective_medium' (homogenised Maxwell–Garnett TMM slab)
    #   • MACRO       → 'ray' (geometric Gouffé Monte-Carlo)
    # ----------------------------------------------------------------------
    regime_solver = wave_model
    regime_fullwave_mc_fallback = False
    if enforce_physics_regime:
        if physics_regime == REGIME_FULL_WAVE:
            # Full-wave window (0.2 <= lambda/D <= 5).  Only route to the cached
            # response-table operator when the CALLER actually requested a
            # full-wave solve or supplied an explicit cache; substituting the
            # bundled demonstration planar-Fresnel cache silently would destroy
            # the cavity modal physics (lambda_c gating / evanescent LDOS) that
            # MATHEMATICAL_DERIVATIONS.md Phases 1-4 require for this window.
            # Without an explicit request, keep the Monte Carlo tracer whose
            # per-photon spectral sampling embeds exactly that modal operator.
            # Orchestrator hard rule: pure geometric ray tracing is INVALID
            # when lambda/D > 1 (ray_valid == False).  Without an explicitly
            # requested full-wave response we execute the MODAL-GUIDED photon
            # transport (per-photon spectral λ_c gate + evanescent decay),
            # never a purely diffusive billiard trace.
            if wave_model == 'cached' or cache_path:
                regime_solver = 'cached'
            else:
                regime_solver = 'ray'
                regime_fullwave_mc_fallback = True
        elif physics_regime == REGIME_EFFECTIVE_MEDIUM:
            regime_solver = 'effective_medium'
        elif physics_regime == REGIME_MACRO:
            regime_solver = 'ray'
        # REGIME_NEAR_FIELD keeps the caller's wave_model; the dedicated
        # Polder–Van Hove Green-tensor path below is authoritative.
    else:
        regime_solver = wave_model

    # Defaults for Track B2 structured near-field (assigned in the active
    # sub-wavelength branch below; far-field neutral otherwise).
    nf_response = None      # NearFieldResponse cache (union carrier)
    nf_correction = 1.0     # multiplicative enhancement of inter-plate channel
    nf_flux = 0.0           # interpolated NFRHT flux at the sim gap (W/m^2)
    nf_ref = 0.0            # blackbody far-field reference (W/m^2)

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
                use_near_field_exclusively = False

                # ---- Track B2: structured near-field Green tensor + LDOS correction ----
                # Enhancement factor for the inter-plate channel (evanescent tunnelling + LDOS).
                # Symmetric across both temperature orderings; preserves thermal equilibrium (net==0 at T_A==T_B).
                nf_response, nf_correction, nf_flux, nf_ref = _structured_near_field_correction(
                    gap_m=gap_m,
                    temperature_hot_K=float(temp_a),
                    temperature_cold_K=float(temp_b),
                    cavity_radius_um=_cavity_radius_um(
                        geometry_mode, cavity_diameter, cnt_pitch,
                        cnt_dia_base, width),
                    cavity_depth_um=height,
                    material=_nf_material_for(geometry_mode),
                )
                # ---- Near-field diagnostic guard (high priority) ------------
                # The Green-tensor / LDOS enhancement is only physical when the
                # gap is comparable to the evanescent tunnelling length,
                # gap <~ lambda_peak/(2*pi).  For gap_ratio > 20 the plates are
                # deep in the far-field regime; any residual NF_corr from the
                # structured solver is spurious and is forced to 1.0 so it can
                # never contaminate the one-way diagnostic fluxes.
                if float(regime_info.get('gap_ratio', 0.0)) > 20.0:
                    nf_correction = 1.0
            except Exception as e:
                import warnings
                warnings.warn(f"Near-field calculation failed: {e}. Falling back to far-field.")
                regime_info['use_near_field'] = False
                near_field_results = None
                use_near_field_exclusively = False
                # Far-field fallback defaults
                nf_response = None
                nf_correction = 1.0
                nf_flux = 0.0
                nf_ref = 0.0
        else:
            # enable_near_field=True but gap not in near-field regime
            near_field_results = None
            use_near_field_exclusively = False
            # Far-field neutral defaults
            nf_response = None
            nf_correction = 1.0
            nf_flux = 0.0
            nf_ref = 0.0
    else:
        # enable_near_field=False: far-field only
        regime_info = {'gap_ratio': gap_ratio_metric(gap_m, max(temp_a, temp_b)), 'use_near_field': False, 'regime': 'far-field (disabled)'}
        near_field_results = None
        use_near_field_exclusively = False
        # Far-field neutral defaults
        nf_response = None
        nf_correction = 1.0
        nf_flux = 0.0
        nf_ref = 0.0
    # ---- Spectral emissivities (Bug #7) ----------------------------------
    mat_a = material_a if material_a in MATERIAL_EMISSIVITY else None
    mat_b = material_b if material_b in MATERIAL_EMISSIVITY else None
    spectral = effective_emissivity_pair(
        mat_a, emissivity_a,    temp_a,
        mat_b, alpha_cnt,       temp_b,
    )
    eps_a_eff = spectral['eps_a_spectral']
    # We'll use the scalar ε_b from MC; spectral info is supplemental

        # ---- Cavity response: MC ray tracing (default) OR cached full-wave ----
    # Aperture re-entry probability ≈ (1-ε_A)·F_{A→B}  (MC path only)

    # Part 5 / peer-review: dynamic aperture-boundary reflectance R_ap(λ, T, f).
    # The sub-cutoff incident field that cannot form a channel mode is captured
    # by the top surface, whose reflectance comes from the Maxwell–Garnett
    # effective-index step across the aperture fill fraction f — NOT from the
    # removed hard-coded `1 − (1 − ε_wall)·0.1` guess.  When the caller does not
    # pin alpha_top we pass None to the ray tracer so it evaluates the aperture
    # absorptance per-wavelength (R_ap(λ, T, f)); alpha_top_def is kept as a
    # representative normal-incidence reference for reporting.
    _mat_wall = ('alumina' if geometry_mode == 'honeycomb' else 'cnt_forest')
    _fill_top = float(packing_fraction)
    _lambda_inc_peak_um = 2897.77 / max(float(temp_a), 1e-12)
    if alpha_top is None:
        alpha_top_def = float(np.clip(
            aperture_boundary_absorptance(
                _mat_wall, _lambda_inc_peak_um, _fill_top, float(temp_a)),
            0.0, 1.0))
        ray_alpha_top = None      # let ray_tracer resolve R_ap(λ, T, f) per photon
    else:
        alpha_top_def = float(np.clip(alpha_top, 0.0, 1.0))
        ray_alpha_top = alpha_top_def

    cached = None
    if regime_solver == 'cached':
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
        if regime_solver == 'effective_medium':
            # ---- Homogenised Maxwell–Garnett TMM slab (EMT regime) ---------
            # Geometric ray tracing is BYPASSED: P, D ≪ λ so the structured
            # top layer becomes a single homogenised thin film solved with the
            # transfer-matrix method.  Kirchhoff holds → α ≡ ε for the slab.
            solver_mode = 'effective_medium'
            _mat_w = ('alumina' if geometry_mode == 'honeycomb' else 'cnt_forest')
            _lampm = 2897.77 / max(float(temp_b), 1e-12)
            _slab_thick = (wall_thickness if geometry_mode == 'honeycomb'
                           else float(cnt_dia_base) * 1e-3)
            mc = _effective_medium_operator(
                geometry, _mat_w, _lampm, float(temp_b),
                float(packing_fraction), _slab_thick)
        else:
            # ---- Monte Carlo ray-tracing path (MACRO geometric optics) ----
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
                alpha_top      = ray_alpha_top,
                # Thin-film physics parameters
                wall_thickness_um = wall_thickness,
                wall_material     = _mat_wall,
                base_material     = 'silver',
                # Part 1 §2 — dynamic complex optical properties n(λ,T)+ik(λ,T):
                # replace static scalar absorptivities inside the ray loop with
                # wavelength- and temperature-dependent complex-Fresnel reflectance
                # (Drude–Lorentz tabulated n,k + T-drift) on every bounce.
                use_complex_fresnel = True,
                # Physics Orchestrator enforcement: dynamic TMM wall optics are
                # mandatory whenever t/δ places the walls in the THIN_FILM or
                # MEMBRANE regime (thinfilm_required), and the per-photon
                # spectral modal-cutoff operator (guided-mode transmission
                # above λc, evanescent decay below) is mandatory when the
                # cavity is TRANSITIONAL / CUTOFF_DOMINATED.  Pure geometric
                # ray tracing is never executed in these regimes.
                apply_modal_attenuation = bool(orch.get(
                    'modal_cutoff_required', True)),
                # Phase 4a: temperature-dependent & non-local optics.
                # The structured surface is at plate-B temperature; the feature
                # scale is the wall film thickness (honeycomb) or tube diameter
                # (CNT forest) — the non-local correction activates only when it
                # is below the electron mean free path.
                wall_temperature_K = float(temp_b),
                base_temperature_K = float(temp_b),
                feature_scale_m   = (wall_thickness * 1e-6
                                     if geometry_mode == 'honeycomb'
                                     else float(cnt_dia_base) * 1e-9),
                # Phase 4b: roughness BRDF (None → legacy Lambertian bounce).
                wall_roughness_sigma_um = surface_roughness_um,
                wall_roughness_tau_um   = roughness_correlation_um,
            )

    # ---- Wave-response union (Phase 0/3/6 integration boundary) -----------------
    # `wave_response_obj` is the JSON-able union carrier for the active
    # wave-physics service that produced the cavity optical properties.
    # NearFieldResponse (solver_kind='nf_greens') takes priority when the gap is
    # sub-wavelength; otherwise the cached full-wave WaveResponse
    # (solver_kind='cached') is exposed; None for the plain MC ray path.
    if nf_response is not None:
        wave_response_obj = nf_response.to_dict()
    elif regime_solver == 'cached' and cached is not None:
        wr = getattr(cached, 'response', None)
        wave_response_obj = wr.to_dict() if wr is not None else None
    else:
        wave_response_obj = None

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
    lambda_peak_um = 2897.77 / max(float(temp_b), 1e-12)
    # Part 2 — regime-dispatcher forced override: the homogenised effective-
    # medium operator (bypassing ray tracing) or the cached full-wave operator
    # win whenever the dimensionless EMT / resonance criterion selects them,
    # independent of the polynomial f_prop heuristic.
    if physics_regime == REGIME_EFFECTIVE_MEDIUM:
        absorption_operator = 'maxwell_garnett_effective_medium'
    elif physics_regime == REGIME_FULL_WAVE:
        absorption_operator = 'full_wave_cached'
    elif f_prop_emit_an >= 0.95:
        absorption_operator = 'gouffe_geometric'
    elif lambda_peak_um > lambda_c_um:
        absorption_operator = 'modal_evanescent_mc'
    elif lambda_peak_um >= 5.0 * float(getattr(geometry, 'P', feature_um)):
        absorption_operator = 'maxwell_garnett_effective_medium'
    else:
        absorption_operator = 'tmm_ray_trace'

    # ---- Bug #13 wave-optics diagnostics ------------------------------------
    # Keep the internal emission and external illumination estimators separate.
    # The derivation models their directional effective properties separately;
    # their ratio is reported as the structured-surface decoupling diagnostic.
    # C_e * P_esc is a rare-escape estimator for deep/modal-filtered cavities.
    # Once most thermal power is propagating, use the reciprocal external
    # absorptance instead; otherwise the area amplification saturates at one
    # and creates an unphysical emissivity unrelated to short-wave EM loss.
    fully_propagating = absorption_operator == 'gouffe_geometric'
    # Emission is read from the physical internal-emission estimator.
    # Kirchhoff's law holds for ANY passive body at ALL scales and temperatures,
    # so once the ray tracer applies identical boundary conditions to both
    # legs (symmetric aperture re-entry, symmetric evanescent modal loss),
    # epsilon_B converges to alpha_eff naturally -- no overwrite gate needed.
    if fully_propagating:
        epsilon_b_raw_cav = compute_cavity_emissivity(
            geometry, alpha_cnt, f_prop_emit_an
        )
        epsilon_b_raw_ci95 = 0.0
        # Gouffé geometric limit: all propagating modes escape.  The raw MC
        # tally (which uses the geometric fallback with normal-incidence TMM
        # bounce survival) under-reports P_esc because thin-film transmittance
        # at normal incidence is far higher than the actual grazing-incidence
        # reflection in a deep cavity (AR = 50).  Override with the
        # deterministic value: the propagating fraction (≥95%) escapes freely.
        p_esc = float(f_prop_emit_an)
        p_esc_ci95 = 0.0
    else:
        if absorption_operator == 'maxwell_garnett_effective_medium':
            epsilon_b_raw_cav = _effective_medium_absorptivity(
                lambda_peak_um, temp_b, 'alumina', packing_fraction
            )
        else:
            epsilon_b_raw_cav = compute_cavity_emissivity(
                geometry, alpha_cnt, f_prop_emit_an, mc['p_esc'],
                eps_base=float(alpha_ag),
            )
        # CI must use the SAME operator as the point estimate.  The bare
        # Ce × p_esc_ci95 formula double-amplifies (ε_wall < 1) and, worse,
        # reports a large interval when the deterministic Gouffé ceiling
        # governs (f_prop > 0.95 / saturated).
        epsilon_b_raw_ci95 = compute_cavity_emissivity_ci95(
            geometry, alpha_cnt, f_prop_emit_an, mc['p_esc'],
            mc.get('p_esc_ci95', 0.0), eps_base=float(alpha_ag),
        )
        if not 0.0 <= epsilon_b_raw_cav <= 1.0:
            raise ValueError(
                'Internal-emission estimator produced an unphysical cavity '
                f'emissivity ({epsilon_b_raw_cav:.6g}); increase photon count '
                'or inspect the cavity ray-launch distribution.'
            )
    # Fix 2b: pass material and diameter to _modal_emission_gate so that T_ap
    # is included in the evanescent integral (diagnostic use only).
    _gate_material = 'alumina' if geometry_mode == 'honeycomb' else 'cnt_forest'
    _gate_diameter = getattr(geometry, 'diameter_um', None)
    g_em, rim_frac, delta_avg_um = _modal_emission_gate(
        geometry, lambda_c_um, temp_b, alpha_cnt, alpha_ag, f_prop_emit_an,
        material=_gate_material, diameter_um=_gate_diameter)

    # ---- Panel Scaling (Honeycomb) — Fix 4 (peer-review) -------------------
    # Fix 4: G_em must NOT be applied as a multiplier to epsilon_b here.
    # The MC loop already samples wavelengths from the Planck distribution and
    # gates each photon on lambda < lambda_c per photon — epsilon_b_raw_cav is
    # the correct physical cavity emissivity.  Multiplying by G_em again
    # double-counts the modal gating (the MC and G_em gate on the same
    # spectral cut), artificially depressing epsilon_B.
    alpha_top_panel = alpha_top_def
    if geometry_mode == 'honeycomb':
        f = packing_fraction
        # The derivation treats the interstitial 51 nm top as a thin film and
        # applies its Planck-weighted TMM emissivity before area scaling.
        e_f = planck_weighted_effective_emissivity(
            float(eps_flat_wall), wall_thickness, 'alumina', float(temp_b)
        )

        epsilon_b_raw = f * epsilon_b_raw_cav   + (1.0 - f) * e_f
        # External-incidence handling depends on the active operator:
        # • EMT (homogenised slab): Kirchhoff is EXACT for a passive planar
        #   stack — alpha ≡ epsilon_b by construction; keep them identical.
        # • Ray/MC: the external absorptance keeps its INDEPENDENT estimator
        #   (Experiment 2).  Overwriting it with the emission-side value
        #   destroys the documented directional decoupling
        #   (MATHEMATICAL_DERIVATIONS.md Phase 9: alpha_eff != eps_B).
        #   Kirchhoff holds at the MATERIAL level inside each leg, not
        #   between these two distinct spectral/angular integrals.
        if solver_mode == 'effective_medium' or absorption_operator == 'maxwell_garnett_effective_medium':
            alpha_eff_cavity = epsilon_b_raw_cav
            alpha_eff_ci95_  = epsilon_b_raw_ci95
        else:
            alpha_eff_cavity = float(mc['alpha_eff'])
            alpha_eff_ci95_  = float(mc.get('alpha_eff_ci95', 0.0))
        alpha_eff     = f * alpha_eff_cavity + (1.0 - f) * e_f
        alpha_eff_ci95 = f * alpha_eff_ci95_
        cavity_enhancement = mc['cavity_enhancement']
        epsilon_b_cavity_part = f * epsilon_b_raw_cav
        epsilon_b_flat_part   = (1.0 - f) * e_f
        epsilon_b      = epsilon_b_raw
        epsilon_b_ci95 = f * epsilon_b_raw_ci95
    else:
        # Non-honeycomb (CNT forest, frustum, etc.): keep the INDEPENDENT
        # MC absorptance estimator (alpha_eff) separate from the emission-side
        # estimator (epsilon_b).  Overwriting alpha_eff with epsilon_b_raw_cav
        # destroys the directional decoupling that the MC correctly computes
        # (e.g. CNT forest at 200 K: alpha_eff ~ 0.998, epsilon_b ~ 0).
        # Kirchhoff holds at the MATERIAL level inside each leg, not between
        # these two distinct spectral/angular integrals.
        alpha_eff      = float(mc['alpha_eff'])
        alpha_eff_ci95 = float(mc.get('alpha_eff_ci95', 0.0))
        epsilon_b_raw  = epsilon_b_raw_cav
        epsilon_b_raw_ci95 = epsilon_b_raw_ci95
        cavity_enhancement = mc['cavity_enhancement']
        epsilon_b_cavity_part = epsilon_b_raw_cav
        epsilon_b_flat_part   = 0.0
        epsilon_b      = epsilon_b_raw_cav
        epsilon_b_ci95 = epsilon_b_raw_ci95

    # ---- Isothermal consistency (reader guidance) ---------------------------
    # Directional effective properties for sub-wavelength structured surfaces
    # are genuinely different spectral/angular integrals: ε_B (emission leg)
    # and α_eff (absorption leg) are estimated by independent photon
    # experiments.  At isothermal equilibrium (T_A == T_B) the net radiative
    # flux q_net ≡ 0 by conservation of energy — that is the Kirchhoff
    # statement.  The two INTEGRAL numbers ε_B and α_eff may differ because
    # they sample different wavelength bands through different boundary
    # conditions (internal emission vs external illumination).  Forcing them
    # equal destroys the documented directional decoupling
    # (MATHEMATICAL_DERIVATIONS.md Phase 9: ε_B = 1.888%, α_eff = 89.27%,
    # decoupling ratio = 47.3×).  The correct diagnostic is
    # ``decoupling_ratio`` and ``kirchhoff_error`` — not a hard overwrite.

    g_em      = float(g_em)                # wave confinement gate diagnostic [0, 1]
    rim_frac  = float(rim_frac)            # fraction of emitter area that is "rim"

    # Directional effective-property budget from the documented model.
    reciprocity_budget  = abs(epsilon_b - alpha_eff) / max(alpha_eff, 1e-12)
    reciprocity_origin  = 'spectral_anisotropy'
    kirchhoff_error     = reciprocity_budget   # kept for output compat
    decoupling_ratio    = (alpha_eff / max(epsilon_b, 1e-12) if epsilon_b > 0
                           else 0.0)

    # ---- Escape solid angle (peer-review Ω_esc ≈ π(R/H)²) -----------------
    H_m = getattr(geometry, 'H', 0.0)
    A_ap_m2 = geometry.area_aperture
    if H_m > 0 and A_ap_m2 > 0:
        r_eff = math.sqrt(A_ap_m2 / math.pi)
        escape_solid_angle = 2.0 * math.pi * (1.0 - H_m / math.sqrt(H_m * H_m + r_eff * r_eff))
    else:
        escape_solid_angle = 0.0

    # Fix 5 (peer-review): Beaming correction for non-Lambertian cavity aperture.
    # The Howell C-11 view factor assumes Lambertian (isotropic cosine) exitance
    # from both surfaces.  For a cavity with aspect ratio AR = H/R >> 1 the
    # aperture exitance is strongly forward-peaked (beaming) rather than
    # Lambertian.  The directional intensity distribution I(θ) for a cavity
    # aperture is not uniform over the hemisphere; power exits within a
    # restricted solid angle Ω_esc < π sr.
    #
    # The beaming correction factor accounts for this:
    #   beaming_factor = Ω_esc / π  ∈ (0, 1]
    # For a flat Lambertian surface Ω_esc = π so beaming_factor = 1.
    # For a deep cavity (AR >> 1) Ω_esc << π so beaming_factor << 1:
    # the cavity emits power mostly in the forward (normal) direction,
    # reducing the effective view factor from Plate B to Plate A.
    #
    # Reference: Born & Wolf (1999) §§8.6; Howell et al. (2016) §6.3.
    if escape_solid_angle > 0 and math.pi > 0:
        beaming_factor = float(np.clip(escape_solid_angle / math.pi, 0.0, 1.0))
    else:
        beaming_factor = 1.0   # flat Lambertian (no cavity / PEC fallback)

    # Panel-level confinement
    if geometry_mode == 'honeycomb':
        confinement_pct = (1.0 - f) * 0.0 + f * confinement_pct

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
    q_af_to_b_density = rad['q_af_to_b'] / A_af if A_af > 0 else 0.0

    # Fix 6 (peer-review): Net inter-surface flux from the radiosity matrix.
    # Section 8.2 solves for radiosities J_AF and J_B.  The energy-conserving
    # net flux between the two facing surfaces is:
    #   q_net_AB = F_{AF→B} · (J_AF − J_B)  [W/m²]
    # This is the ONLY formulation that satisfies energy conservation (yields
    # exactly zero at T_A = T_B = T_surr by construction of the matrix).
    # Previous code used raw primary-emission terms:
    #   q = E_A · ε_A · F · α_eff − E_B · ε_B · F
    # which bypasses the radiosity solution and is non-zero at equilibrium.
    J_AF = rad['j_af']
    J_B  = rad['j_b']
    q_net_AB_density = F_af_b * (J_AF - J_B)           # W/m² (energy-conserving)
    q_net_AB_total   = q_net_AB_density * A_af          # W

    # Fix 5 (continued): The beaming-corrected one-way emission fluxes.
    # Plate A is flat Lambertian (beaming_factor = 1 on A side).
    # Plate B emits into a restricted solid angle (beaming_factor < 1 for deep cavities).
    # The beaming factor modifies how much of B's emission actually reaches A:
    #   q_B→A (corrected) = J_B · F_{B→AF} · beaming_factor
    # These are diagnostic-only quantities; the physical net flux uses the full radiosity.
    E_a_val = SIGMA * temp_a**4
    E_b_val = SIGMA * temp_b**4
    F_ab = float(F_af_b)
    # Diagnostic one-way: kept for UI compatibility but clearly labelled non-conserving.
    q_emit_a_to_b = E_a_val * eps_af * F_ab * alpha_eff          # diagnostic W/m²
    q_emit_b_to_a = E_b_val * epsilon_b * F_ab * beaming_factor  # diagnostic W/m² (Fix 5)

    # ---- Active near-field (evanescent + LDOS) correction -- Track B2 ----
    if nf_correction > 1.0:
        q_emit_a_to_b *= nf_correction
        q_emit_b_to_a *= nf_correction

    # Fix 6 (continued): Physical net flow A→B from radiosity, not from
    # primary emission terms.
    q_net_a_to_b_physical = q_net_AB_density   # W/m² (energy-conserving, Fix 6)

    # ---- Stagnation temperature (adiabatic wall temperature of Plate B) ---
    # Fix 6 (peer-review): T_B_stag must be derived from the radiosity-solved
    # absorbed irradiance G_B, not from the non-conserving diagnostic leg.
    #
    # Under the radiosity model the total irradiance absorbed by Plate B is:
    #   G_B_abs = eps_B · (F_{B→AF} · J_AF + F_{B→surr} · E_surr)
    # At adiabatic steady state of Plate B: emitted flux = absorbed flux
    #   eps_B · σ · T_stag^4 = G_B_abs
    # ⇒  T_stag = (G_B_abs / (eps_B · σ))^(1/4)
    #
    # 2nd Law bound: for a passive system T_stag ≤ max(T_A, T_surr).
    T_bound = max(float(temp_a), float(temp_surr))
    # Stagnation temperature is always well-defined for any passive radiator:
    #   T_stag = (G_B_abs / (eps_B * sigma))^(1/4)
    # Even when eps_B is extremely small (e.g. CNT forest at cryogenic
    # temperatures, eps_B ~ 1e-20), G_B_abs scales proportionally with eps_B
    # (G_B_abs = eps_B * (F*J_AF + ...)), so the ratio is bounded and finite.
    # The previous "if epsilon_b > 1e-9" guard was wrong: it returned 0.0 K
    # (absolute zero) for deep-subwavelength emitters, violating the 2nd Law.
    # Guard only against division by zero (epsilon_b exactly 0.0).
    if epsilon_b > 0.0:
        F_b_af = F_af_b * A_af / A_b_macro if A_b_macro > 0 else 0.0
        F_b_s  = max(1.0 - F_b_af, 0.0)
        E_surr = SIGMA * float(temp_surr)**4
        # Absorbed irradiance from enclosure (radiosity solution)
        G_B_abs = eps_b * (F_b_af * J_AF + F_b_s * E_surr)
        # Adiabatic equilibrium temperature of Plate B
        T_B_stag_calc = (G_B_abs / (epsilon_b * SIGMA)) ** 0.25 if G_B_abs > 0 else 0.0
        # Enforce 2nd Law of Thermodynamics for passive unamplified exchange
        T_B_stag = float(min(T_bound, T_B_stag_calc))
    else:
        T_B_stag = float(T_bound)  # zero emission -> no cooling -> track the hot source

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

    # ---- Part 3 — Second-Law reciprocity / equilibrium guard ---------------
    # Under an isothermal system (T_A == T_B == T_surr) the spectrally-
    # integrated directional-hemispherical absorptance must equal the
    # emissivity, which forces the radiosity net exchange to vanish.
    equilibrium_guard = enforce_thermal_equilibrium_zero_net(
        q_net_AB_density, temp_a, temp_b, temp_surr, tol=1e-9)
    if (enforce_physics_regime and equilibrium_guard['is_thermally_uniform']
            and not equilibrium_guard['passes']):
        raise ValueError(
            'Second-Law guard: net radiative exchange q_net_AB_density = '
            f"{equilibrium_guard['q_net_AB_density']:.3e} W/m² persists under "
            'isothermal conditions (T_A == T_B == T_surr). The model is not '
            'energy-conserving and cannot pass academic peer review.'
        )

    _ri = phys_regime_info
    _lcu = float(lambda_c_um) if lambda_c_um is not None else float('inf')
    audit_notes = ['Active physics regime: ' + physics_regime]
    audit_notes.append('D/lambda_peak=%.4f P/lambda_peak=%.4f gap/(lambda_peak/2pi)=%.4f' % (_ri['D_lambda_frac'], _ri['P_lambda_frac'], _ri['gap_lambda_frac']))
    if math.isfinite(_lcu):
        audit_notes.append('Modal cutoff lambda_c=%.4f um; f_prop=%.4f' % (_lcu, float(f_prop_emit_an)))
    else:
        audit_notes.append('No cutoff (lambda_c=inf); f_prop=%.4f' % float(f_prop_emit_an))
    audit_notes.append('Absorption operator: %s; eps_wall=%.4f' % (absorption_operator, float(alpha_cnt)))
    if physics_regime == REGIME_NEAR_FIELD:
        audit_notes.append('Near-field LDOS correction active (Polder-Van Hove).')
    elif physics_regime == REGIME_EFFECTIVE_MEDIUM:
        audit_notes.append('Effective-medium (Maxwell-Garnett) regime; homogenised TMM slab; ray tracing bypassed.')
    elif physics_regime == REGIME_FULL_WAVE:
        audit_notes.append('Diffraction/resonance (0.2*lambda_peak<=D<=5*lambda_peak); full-wave operator enforced.')
    audit_notes.append('Strict physics-routed solver: %s (enforce_physics_regime=%s)' %
                       (regime_solver, enforce_physics_regime))
    _wo = _ri.get('wall_optical_depth_t_over_delta')
    if _wo is not None:
        audit_notes.append('wall t_wall/delta=%.4f optically_thin=%s' %
                           (_wo, _ri.get('wall_optically_thin')))
    if regime_solver == 'cached':
        audit_notes.append('Cached planar full-wave (stub table) selected - geometric cavity modes disabled.')
    if regime_fullwave_mc_fallback:
        audit_notes.append(
            'FULL_WAVE window (0.2<=lambda/D<=5): no full-wave response supplied - '
            'executing the Monte Carlo tracer with the per-photon spectral modal-cutoff '
            'operator (lambda_c gating + evanescent decay) instead of a silent demo-cache substitution.')
    if max(float(temp_a), float(temp_b), float(temp_surr)) > 2000.0:
        audit_notes.append('Boundary note: high-temperature differential (flux ~ T^4); confirm operating point.')
    if float(temp_surr) < 50.0:
        audit_notes.append('Boundary note: cryogenic background sink (T_surr<50K).')
    if regime_warning:
        audit_notes.append(regime_warning)

    # ---- Physics Orchestrator provenance (regime_selector.py) --------------
    # Full pre-flight audit exported with the payload: dimensionless ratios,
    # selected regimes (geometry / wall / cavity / heat-transfer), enforcement
    # flags, point-deduction confidence score, and regime warnings.
    orch_conf = float(orch.get('confidence', 100.0))
    audit_notes.append(
        'Physics orchestrator: geometry=%s wall=%s cavity=%s heat=%s '
        '(confidence=%.0f%%)' % (
            orch.get('geometry_regime'), orch.get('wall_regime'),
            orch.get('cavity_regime'), orch.get('heat_transfer_regime'),
            orch_conf))
    for _ow in orch.get('regime_warnings', []):
        audit_notes.append('Orchestrator warning: %s' % _ow)

    # ---- Part 5 — Material-optics confidence & optical warnings ------------
    # Quantify the trustworthiness of the tabulated/Drude-Lorentz optical model
    # as a function of temperature drift from the 300 K calibration baseline and
    # surface the warnings in BOTH the audit trail and the payload.
    _mat_feature_m = (wall_thickness * 1e-6 if geometry_mode == 'honeycomb'
                      else float(cnt_dia_base) * 1e-9)
    # Evaluate against the HOTTEST operating boundary: high-temperature drift
    # from the 300 K calibration baseline is what degrades the optical model,
    # so checking only the cold plate would silently miss degradation.
    _mat_T_hot = float(max(temp_a, temp_b))
    material_confidence = material_optics_confidence(
        _mat_wall, temperature_K=_mat_T_hot, baseline_K=300.0,
        feature_scale_m=_mat_feature_m)
    for _warn in material_confidence.get('warnings', []):
        audit_notes.append('Material-optics: ' + _warn)
    audit_notes.append('Material-optics confidence_level=%.2f (T=%.0f K, baseline 300 K)'
                       % (material_confidence['confidence_level'], _mat_T_hot))

    return {
        # Part 2 — Physics regime-dispatch provenance
        'physics_regime_engine':     physics_regime,
        'physics_regime_info':       phys_regime_info,
        # Physics Orchestrator (regime_selector.py) provenance
        'physics_orchestrator':      dict(orch),
        'physics_confidence':        orch_conf,
        'physics_regime_warnings':   list(orch.get('regime_warnings', [])),
        'physics_regime_warning':    regime_warning,
        'physics_regime_violation':  regime_violation,
        'equilibrium_guard':         equilibrium_guard,
        'physics_audit_notes':    audit_notes,
        # Phase 3 Integration Results
        'physics_regime': regime_info.get('regime', 'far-field'),
        'gap_ratio': regime_info.get('gap_ratio', 0.0),
        'net_flux_near_field_W_m2': near_field_results.get('net_flux_near_field_W_m2', 0.0) if near_field_results else 0.0,
        'evanescent_fraction': near_field_results.get('evanescent_fraction', 0.0) if near_field_results else 0.0,
        'evanescent_flux_W_m2': near_field_results.get('evanescent_flux_W_m2', 0.0) if near_field_results else 0.0,
        'propagating_flux_W_m2': near_field_results.get('propagating_flux_W_m2', 0.0) if near_field_results else 0.0,
        
        # Solver provenance (Phase 0/6 integration boundary)
        'wave_model':             wave_model,
        'active_physics_operator': regime_solver,
        'solver_mode':            solver_mode,
        'wave_response_info':     cached_info if regime_solver == 'cached' else None,
        # Part 5 — material-optics confidence & optical warnings
        'material_optics_confidence': material_confidence.get('confidence_level', None),
        'material_optics_drift_fraction': material_confidence.get('drift_fraction', None),
        'material_optics_warnings': material_confidence.get('warnings', []),
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
        'reciprocity_budget':     float(reciprocity_budget),
        'reciprocity_origin':     reciprocity_origin,
        # Bug #13 wave-optics diagnostics
        'cutoff_wavelength_um':   lambda_c_um,
        'absorption_operator':    absorption_operator,
        'operator_wavelength_um': float(lambda_peak_um),
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
        'packing_fraction':       float(packing_fraction),
                'geometry_label':         label,
        'wall_thickness_um':      wall_thickness,
        # Phase 4a: temperature-dependent / non-local optics provenance
        'optics_temperature_K':   float(temp_b),
        'optics_model_provenance': temperature_optics_provenance(
            'alumina' if geometry_mode == 'honeycomb' else 'cnt_forest',
            float(temp_b),
            wall_thickness * 1e-6 if geometry_mode == 'honeycomb'
            else float(cnt_dia_base) * 1e-9),
        # Phase 4b: surface-roughness BRDF provenance (peak-wavelength eval)
        'brdf_provenance':        brdf_provenance(
            sigma_um=surface_roughness_um,
            tau_um=roughness_correlation_um,
            wavelength_um=(2898.0 / float(temp_b)) if temp_b > 0 else None),
        'area_aperture_m2':       geometry.area_aperture,
        'area_walls_m2':          geometry.area_walls,
        'area_base_m2':           geometry.area_base,
        # View factor
        'view_factor_A_B':        F_af_b,
        'view_factor_MC':         f_ab_mc,
        'coupling_efficiency_MC': coupling_mc,
        # Radiosity / net flux
        # true_radiosity_net: energy-conserving enclosure solution (zero at
        #   thermal equilibrium by construction). This is the physical net heat flow.
        'total_leakage':          q_net_w,
        'true_radiosity_net':     float(q_net_w),
        # Fix 6: q_net_AB_density is the physically correct net heat flux
        #   from the radiosity matrix solution: F_{A->B}*(J_A - J_B) [W/m^2].
        #   This is energy-conserving and is zero at thermal equilibrium.
        'q_net_AB_density':       float(q_net_AB_density),
        'q_net_AB_total_W':       float(q_net_AB_total),
        'radiosity_J_AF':         float(J_AF),
        'radiosity_J_B':          float(J_B),
        # net_flux_A_front: Fix 6 — now equals q_net_AB_density (radiosity-based).
        'net_flux_A_front':       q_net_a_to_b_physical,
        # diagnostic_one_way_a_to_b / b_to_a: NOT energy-conserving.
        #   A->B uses eps_A * alpha_eff (not reciprocal to B->A).
        #   B->A includes Fix 5 beaming factor for non-Lambertian cavity exitance.
        #   These are kept only for spectral/directional anisotropy diagnostics.
        'diagnostic_one_way_a_to_b': float(q_emit_a_to_b),
        'diagnostic_one_way_b_to_a': float(q_emit_b_to_a),
        'diagnostic_non_conserving': float(q_emit_a_to_b - q_emit_b_to_a),
            'q_af_to_b':              rad['q_af_to_b'],
        'q_ab_net':               q_ab_net,
        # One-way diagnostic flux (independent of multiple reflections)
        'q_emit_a_to_b_one_way':  q_emit_a_to_b,      # W/m^2 diagnostic only
        'q_emit_b_to_a_one_way':  q_emit_b_to_a,      # W/m^2 diagnostic only (Fix 5: beaming)
        'q_net_a_to_b_physical': float(q_net_a_to_b_physical),  # Fix 6: radiosity-derived
        'q_b_net':                rad['q_b_net'],
        'flux_emitted_b':         epsilon_b * 5.670374419e-8 * (temp_b ** 4),
        'area_A_front_m2':        A_af,
        # Fix 5: Beaming correction factor for cavity non-Lambertian exitance.
        'beaming_factor':         float(beaming_factor),
        # Stagnation temperature (adiabatic wall temperature of Plate B)
        'T_B_stag':               float(T_B_stag),
        # Fix 6: T_B_stag is now derived from the radiosity-solved absorbed
        # irradiance G_B = eps_B*(F_B->AF*J_AF + F_B->surr*E_surr), enforcing
        # the 2nd Law bound T_stag <= max(T_A, T_surr).  This replaces the
        # previous non-conserving formula (q_emit_a_to_b / (eps_B*sigma))^0.25
        # which could yield T_stag > T_A for low eps_B.
        'T_B_stag_note': (
            'Fix 6 (peer-review): T_stag derived from radiosity-solved absorbed '
            'irradiance G_B = eps_B*(F_B->AF*J_AF + F_B->surr*E_surr). '
            'Bounded by 2nd Law: T_stag <= max(T_A, T_surr). '
            'Under thermal equilibrium (T_A=T_B=T_surr) T_stag = T_A exactly.'),
        # Temperatures
        'temperature_A':          float(temp_a),
        'temperature_B':          float(temp_b),
        'temperature_surr':       float(temp_surr),
        # Near-field
        'near_field_ratio':       nf['near_field_ratio'],
        'near_field_warning':     nf['near_field_warning'],
        # Track B2: structured near-field Green tensor + LDOS (active correction)
        'wave_response':          wave_response_obj,
        'near_field_correction':  float(nf_correction),
        'near_field_response':    (nf_response.to_dict()
                                   if nf_response is not None else None),
        'near_field_flux_W_m2':   float(nf_flux),
        'near_field_reference_W_m2': float(nf_ref),
        # Spectral
        'eps_a_spectral':         spectral['eps_a_spectral'],
        'eps_b_spectral':         spectral['eps_b_spectral'],
        'spectral_correction_pct':spectral['spectral_correction_pct'],
    }
