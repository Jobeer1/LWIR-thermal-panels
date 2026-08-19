"""
ray_tracer.py — 3-D Monte Carlo ray tracer for cavity thermal emission.

Key improvements over the original 2-D tracer:
  ✓ Full 3-D geometry (rect pit, frustum, CNT forest unit cell)
  ✓ Correct 3-D Lambertian sampling (Malley's method)
  ✓ Weight-based photon tracking with proper Russian Roulette
  ✓ Aperture re-entry coupling from plate A reflections
  ✓ Statistical 95% confidence interval on P_esc and α_eff
  ✓ Returns all results needed by the radiosity network

Bug #13 (peer-review physics) — spectral MC with waveguide modal cutoff:
  Each photon's wavelength is sampled from the Planck distribution at the
  relevant temperature (T_emit for internal emission, T_inc for incident
  plate-A radiation).  The channel's waveguide cutoff λ_c is set by the
  geometry (TE11 for circular pores λ_c = 1.706·d, TE10 for the square
  inter-tube gap λ_c = 2·gap):

    • λ < λ_c  : the mode PROPAGATES → full geometric ray tracing.
    • λ ≥ λ_c  : the mode is EVANESCENT → exponentially-decaying field.

  For internal emission the evanescent photon tunnels to the aperture with
  transmission T = exp(−L/δ_ev), where L is the path to the aperture plane
  and δ_ev = λ_c/(2π)·(1 − (λ_c/λ)²)^−1/2 is the decay length.  Deep
  emission (L ≫ δ_ev) is confined — this is the LDOS suppression of Lin
  et al. (PRB 2000) and the sub-wavelength cutoff of Narayanaswamy & Chen
  (PRB 2004).  For external incidence the sub-cutoff photon cannot enter
  the channel and is captured by the lossy top surface (graded-index /
  diffractive trap), keeping α_eff → 1.

  ε_B is now computed from the physical internal-emission experiment
  (no longer forced equal to α_eff), exposing the anisotropic decoupling
  α_eff ≠ ε_eff for deep sub-wavelength channels.

Usage
-----
    from geometry import RectPit3D
    from ray_tracer import run_cavity_mc_3d

    geom = RectPit3D(width_um=10, depth_um=10, height_um=450)
    results = run_cavity_mc_3d(geom, n_photons=50000,
                               eps_walls=0.98, eps_base=0.02)
"""

import math
import numpy as np
from typing import Union

from geometry import RectPit3D, FrustumCavity3D, CNTForestCell
from sampling  import sample_hemisphere_3d, sample_planck_wavelength

AnyGeometry = Union[RectPit3D, FrustumCavity3D, CNTForestCell]

# ---------------------------------------------------------------------------
# Single-photon tracer
# ---------------------------------------------------------------------------

# Weight below which Russian Roulette is triggered
_RR_THRESHOLD = 1e-4
# Survival probability boost in RR (must be > _RR_THRESHOLD)
_RR_BOOST = 0.1
# Maximum bounces before forced termination (should never be reached in practice
# for realistic high-ε walls; exists only as a safety net)
_MAX_BOUNCES = 2000


def _trace_photon(pos: np.ndarray,
                  direction: np.ndarray,
                  geometry: AnyGeometry,
                  eps_walls: float,
                  eps_base:  float,
                  re_entry_prob: float = 0.0) -> float:
    """Trace one photon from *pos* in *direction* inside *geometry*.

    Returns the photon's surviving energy weight at escape (> 0 means it
    left through the aperture), or 0 if absorbed.

    Parameters
    ----------
    pos, direction    : initial position (m) and unit direction vector.
    geometry          : cavity geometry object.
    eps_walls         : wall emissivity / absorptivity.
    eps_base          : base emissivity / absorptivity.
    re_entry_prob     : probability that a photon escaping the aperture is
                        immediately re-injected downward (models plate A
                        reflection — Bug #9 fix).  Typically (1-ε_A)·F_AB.
    """
    weight = 1.0

    bounces = 0
    total_steps = 0
    while bounces < _MAX_BOUNCES and total_steps < 1000:
        total_steps += 1
        t, surface, normal = geometry.next_hit(pos, direction)

        if not math.isfinite(t) or t <= 0.0:
            return 0.0

        pos = pos + direction * t

        if surface == 'aperture':
            if re_entry_prob > 0.0 and np.random.random() < re_entry_prob:
                direction = sample_hemisphere_3d(np.array([0.0, 0.0, -1.0]))
                pos = pos - np.array([0.0, 0.0, 1e-12])
                bounces += 1
                continue
            return weight

        if surface == 'periodic_x':
            pos[0] = -np.sign(pos[0]) * (geometry.P / 2 - 1e-13)
            continue
        elif surface == 'periodic_y':
            pos[1] = -np.sign(pos[1]) * (geometry.P / 2 - 1e-13)
            continue

        bounces += 1

        if surface == 'wall' or surface == 'top_cap':
            eps = eps_walls
        elif surface == 'base':
            eps = eps_base
        else:
            return 0.0

        weight *= (1.0 - eps)

        if weight < _RR_THRESHOLD:
            if np.random.random() < weight / _RR_BOOST:
                weight = _RR_BOOST
            else:
                return 0.0

        direction = sample_hemisphere_3d(normal)
        pos = pos + normal * 1e-13

    return 0.0


# ---------------------------------------------------------------------------
# Waveguide modal-cutoff helpers (Bug #13 — peer-review physics)
# ---------------------------------------------------------------------------

def evanescent_decay_length(lambda_c_um: float, lambda_um: float) -> float:
    """Evanescent decay length δ_ev (µm) of a sub-cutoff waveguide mode.

    For a mode at λ > λ_c the axial propagation constant becomes imaginary:
        β = i·κ,   κ = 2π·(1/λ_c² − 1/λ²)^(1/2)
    so the field decays along the channel as exp(−z/δ_ev) with

        δ_ev = 1/κ = (λ_c/2π) · (1 − (λ_c/λ)²)^(−1/2)

    Limits:
      • λ → λ_c⁺  : δ_ev → ∞  (mode just below cutoff is nearly propagating)
      • λ ≫ λ_c   : δ_ev → λ_c/2π  (deeply evanescent, ~0.16 of a channel
                    diameter — matches the "fails to reach the aperture" claim)

    Returns the decay length in the same length units as the inputs (µm).
    """
    lambda_c_um = float(lambda_c_um)
    lambda_um   = float(lambda_um)
    if lambda_c_um <= 0.0 or not math.isfinite(lambda_c_um):
        return 0.0                      # no channel → fully confined
    if lambda_um < lambda_c_um:
        return float('inf')             # propagating mode: no decay
    ratio = lambda_c_um / lambda_um     # ∈ (0, 1]
    return (lambda_c_um / (2.0 * math.pi)) / math.sqrt(max(1e-12, 1.0 - ratio * ratio))

# ---------------------------------------------------------------------------
# Cavity MC driver
# ---------------------------------------------------------------------------

def run_cavity_mc_3d(
    geometry:       AnyGeometry,
    n_photons:      int,
    eps_walls:      float,
    eps_base:       float,
    eps_aperture:   float = 0.0,   # plate A emissivity (for re-entry)
    view_factor_ab: float = 1.0,   # F_{A→B} (for re-entry probability)
    T_emit:         float = 300.0, # plate B temperature — Planck spectrum for emission
    T_inc:          float = 600.0, # plate A temperature — Planck spectrum for incidence
    alpha_top:      float = None,  # top-surface absorptivity for sub-cutoff incident light
) -> dict:
    """Run the two cavity MC experiments and compute ε_B and α_eff.

    Experiment 1 — Internal emission (P_esc):
        Photons emitted from all internal surfaces (walls + base) weighted by
        area × emissivity.  Each photon gets a wavelength sampled from the
        Planck distribution at *T_emit*:
          • λ < λ_c  → geometric ray tracing until escape or absorption.
          • λ ≥ λ_c  → the channel mode is EVANESCENT; the photon tunnels to
                       the aperture with T = exp(−L/δ_ev) (≈0 for deep emission).
    Experiment 2 — External incidence (α_eff):
        Photons enter through the aperture with a Lambertian downward
        distribution; wavelength sampled from Planck at *T_inc*:
          • λ < λ_c  → geometric tracing (cavity trapping).
          • λ ≥ λ_c  → the wave cannot form a propagating channel mode; it is
                       captured by the lossy top surface (graded-index /
                       diffractive trap) with probability *alpha_top*.
    ε_B is the physically-computed emissivity (cavity_enhancement × mean
    escape weight).  It is NOT forced equal to α_eff: for deep sub-wavelength
    channels the modal cutoff makes ε_B < α_eff (anisotropic decoupling).

    Parameters
    ----------
    geometry       : RectPit3D, FrustumCavity3D, or CNTForestCell.
    n_photons      : number of photons per experiment.
    eps_walls      : wall absorptivity/emissivity.
    eps_base       : base absorptivity/emissivity.
    eps_aperture   : plate A emissivity (used for aperture re-entry).
    view_factor_ab : F_{A→B} at macro scale (used for re-entry).
    T_emit         : temperature [K] used for the internal-emission spectrum.
    T_inc          : temperature [K] used for the incident-radiation spectrum.
    alpha_top      : absorptivity of the top surface for sub-cutoff incidence
                     (None → eps_walls).

    Returns
    -------
    dict with keys:
      p_esc, p_esc_ci95            : escape probability and 95% CI half-width
      alpha_eff, alpha_eff_ci95    : effective absorptivity and 95% CI
      epsilon_b, epsilon_b_raw     : effective emissivity (physical emission)
      epsilon_b_ci95               : 95% CI on ε_B
      cavity_enhancement           : A_emitting / A_aperture
      kirchhoff_error              : |ε_B - α_eff| / α_eff * 100 [%]
      f_prop_emit, f_prop_inc      : sampled propagating-mode fractions
      n_evan                       : number of evanescent (confined) photons
    """
    re_entry = (1.0 - eps_aperture) * view_factor_ab
    if alpha_top is None:
        # Graded-index / diffractive top-surface capture (Bug #13): sub-wavelength
        # incident light cannot form a channel mode; it folds around the wall rims
        # and is trapped near the surface.  Specular reflection is suppressed
        # (peer review: α_eff → 1 for d ≈ λ).  Residual reflection = 10% of the
        # bare wall reflectivity (1 − ε_walls).
        alpha_top = float(np.clip(1.0 - (1.0 - float(eps_walls)) * 0.1, 0.0, 1.0))
    alpha_top = float(np.clip(alpha_top, 0.0, 1.0))

    # Waveguide cutoff of this channel geometry (µm)
    lambda_c_um = float(getattr(geometry, 'lambda_c_um', float('inf')))

    # ---- Experiment 1: Internal emission -----------------------------------
    A_walls  = geometry.area_walls
    A_base   = geometry.area_base
    w_walls  = A_walls * eps_walls
    w_base   = A_base  * eps_base
    w_total  = w_walls + w_base
    p_wall   = w_walls / w_total if w_total > 0 else 0.5

    escaped_weight    = 0.0
    escaped_weight_sq = 0.0
    n_prop            = 0
    n_evan            = 0
    tunnel_escape_w   = 0.0
    for _ in range(n_photons):
        if np.random.random() < p_wall:
            pos, n_out = geometry.sample_point_on_walls()
        else:
            pos, n_out = geometry.sample_point_on_base()
        direction = sample_hemisphere_3d(n_out)
        # Nudge off surface
        pos = pos + n_out * 1e-13

        lam = sample_planck_wavelength(T_emit)
        if lam < lambda_c_um:
            n_prop += 1
            w = _trace_photon(pos, direction, geometry, eps_walls, eps_base, re_entry)
        else:
            n_evan += 1
            # Sub-cutoff photon: the channel mode is EVANESCENT (cannot
            # propagate).  Its escape probability factors into
            #     w = p_geo(escape) × T_tunnel(L, λ)
            # where p_geo is the geometric aperture-escape probability
            # (solid angle Ω_esc ≈ π(R/H)², single-bounce re-absorption,
            # multi-bounce coupling) and T_tunnel = exp(−L/δ_ev) is the
            # frustrated-tunnelling transmission of the decaying evanescent
            # field through the L-µm of channel remaining to the aperture
            # (δ_ev = (λ_c/2π)·(1−(λ_c/λ)²)^−1/2; deep sub-wavelength
            # channels ⇒ δ_ev ≈ λ_c/2π ≪ 1 → "waves fail to reach the
            # aperture" — Narayanaswamy & Chen, PRB 2004).
            L = max(geometry.H - pos[2], 0.0) * 1e6   # µm to aperture plane
            delta_ev = evanescent_decay_length(lambda_c_um, lam)
            w_geo = _trace_photon(pos, direction, geometry,
                                  eps_walls, eps_base, re_entry)
            if math.isfinite(delta_ev) and delta_ev > 0.0:
                w = w_geo * math.exp(-L / delta_ev)
            else:
                w = w_geo if L <= 0.0 else 0.0
            tunnel_escape_w += w
        escaped_weight    += w
        escaped_weight_sq += w * w

    p_esc = escaped_weight / n_photons
    # Correct CI for weighted estimator: Var(p_esc) = Var(w_i) / n
    var_w = (escaped_weight_sq / n_photons) - (p_esc ** 2)
    p_esc_ci95 = 1.96 * math.sqrt(max(var_w / n_photons, 0.0))

    # ---- Experiment 2: External incidence (Kirchhoff check) ----------------
    A_ap = geometry.area_aperture

    # Aperture sampling by geometry type
    if hasattr(geometry, 'W'):                    # RectPit3D
        def _sample_aperture():
            return np.array([np.random.uniform(0.0, geometry.W),
                             np.random.uniform(0.0, geometry.D),
                             geometry.H])
    elif hasattr(geometry, 'r_top'):              # FrustumCavity3D
        def _sample_aperture():
            r   = geometry.r_top * math.sqrt(np.random.random())
            phi = np.random.uniform(0, 2 * math.pi)
            return np.array([r * math.cos(phi), r * math.sin(phi), geometry.H])
    elif hasattr(geometry, 'pitch_um'):           # CNTForestCell
        def _sample_aperture():
            return np.array([np.random.uniform(-geometry.P / 2, geometry.P / 2),
                             np.random.uniform(-geometry.P / 2, geometry.P / 2),
                             geometry.H + 1e-13])
    elif hasattr(geometry, 'diameter_um'):        # HoneycombCavityCell
        def _sample_aperture():
            r = geometry.R * math.sqrt(np.random.random())
            phi = np.random.uniform(0, 2 * math.pi)
            return np.array([r * math.cos(phi), r * math.sin(phi), geometry.H])
    else:
        def _sample_aperture():
            return np.array([0.0, 0.0, geometry.H])

    down = np.array([0.0, 0.0, -1.0])
    absorbed_weight    = 0.0
    absorbed_weight_sq = 0.0
    n_prop_inc         = 0
    for _ in range(n_photons):

        pos = _sample_aperture()
        lam = sample_planck_wavelength(T_inc)
        if lam < lambda_c_um:
            n_prop_inc += 1
            direction = sample_hemisphere_3d(down)
            w = _trace_photon(pos, direction, geometry, eps_walls, eps_base, 0.0)
            a = 1.0 - w          # absorbed fraction of this photon
        else:
            # Sub-cutoff incident wave: cannot form a propagating channel mode.
            # Captured by the lossy top surface (graded-index / diffractive
            # trapping) rather than entering the cavity.
            a = alpha_top
        absorbed_weight    += a
        absorbed_weight_sq += a * a

    alpha_eff = absorbed_weight / n_photons
    var_a = (absorbed_weight_sq / n_photons) - (alpha_eff ** 2)
    alpha_eff_ci95 = 1.96 * math.sqrt(max(var_a / n_photons, 0.0))

    # ---- Derived quantities ------------------------------------------------
    # ε_B is computed from the PHYSICAL emission experiment only.  The wave
    # optics (modal cutoff) makes it smaller than the geometric cavity
    # enhancement × escape prediction whenever the channels are sub-wavelength
    # (anisotropic decoupling — the operational result of the peer review).
    cavity_enhancement = w_total / A_ap if A_ap > 0 else 0.0
    epsilon_b_raw      = cavity_enhancement * p_esc
    epsilon_b          = epsilon_b_raw
    epsilon_b_ci95     = cavity_enhancement * p_esc_ci95
    kirchhoff_error    = (abs(epsilon_b_raw - alpha_eff) / max(alpha_eff, 1e-12)) * 100.0

    f_prop_emit = n_prop / n_photons if n_photons > 0 else 0.0
    f_prop_inc  = n_prop_inc / n_photons if n_photons > 0 else 0.0

    return {
        'p_esc':              p_esc,
        'p_esc_ci95':         p_esc_ci95,
        'alpha_eff':          alpha_eff,
        'alpha_eff_ci95':     alpha_eff_ci95,
        'epsilon_b':          epsilon_b,
        'epsilon_b_raw':      epsilon_b_raw,
        'epsilon_b_ci95':     epsilon_b_ci95,
        'cavity_enhancement': cavity_enhancement,
        'kirchhoff_error':    kirchhoff_error,
        'f_prop_emit':        f_prop_emit,
        'f_prop_inc':         f_prop_inc,
        'tunnel_escape_w':    tunnel_escape_w,
        'n_evan':             n_evan,
        'w_total':            w_total,
        'area_aperture':      A_ap,
    }
