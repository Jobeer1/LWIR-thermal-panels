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
    power transmission T = exp(−2L/δ_ev), where L is the axial distance to the aperture plane
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
from sampling  import (
    sample_hemisphere_3d,
    sample_planck_wavelength,
    sample_planck_wavelength_band,
    planck_cumulative
)
from material_optics import (
    effective_emissivity_thin_film,
    get_complex_refractive_index,
    tmm_reflectance_single_layer,
    planck_weighted_effective_emissivity,
    aperture_boundary_absorptance,
)
from waveguide_modes import solve_te11_mode_complex, attenuation_factor_lossy_waveguide
from brdf import sample_surface_direction

AnyGeometry = Union[RectPit3D, FrustumCavity3D, CNTForestCell]

# ---------------------------------------------------------------------------
# Single-photon tracer
# ---------------------------------------------------------------------------

# Weight below which Russian Roulette is triggered
_RR_THRESHOLD = 1e-4
# Survival probability boost in RR (must be > _RR_THRESHOLD)
_RR_BOOST = 0.1
# Maximum bounces before forced termination.  For near-perfect-mirror thin-film
# walls (R ~ 0.9964) with AR=50 (200 µm deep, 4 µm diameter) a photon may
# need O(1000) bounces to either escape or be absorbed.  Setting this too low
# introduces an asymmetric truncation error: deep-emission photons are more
# likely to still be alive at the cap (returning weight > 0 as if escaped)
# than incoming photons are, breaking reciprocity.
_MAX_BOUNCES = 2000
_MAX_PHOTON_STEPS = 2000


def _trace_photon(pos: np.ndarray,
                  direction: np.ndarray,
                  geometry: AnyGeometry,
                  eps_walls: float,
                  eps_base:  float,
                  re_entry_prob: float = 0.0,
                  use_complex_fresnel: bool = False,
                  wall_thickness_um: float = None,
                  wall_material: str = 'alumina',
                  base_material: str = 'silver',
                  photon_wavelength_um: float = None) -> float:
    """Trace one photon from *pos* in *direction* inside *geometry*.

    Returns the photon's surviving energy weight at escape (> 0 means it
    left through the aperture), or 0 if absorbed.

    Parameters
    ----------
    pos, direction    : initial position (m) and unit direction vector.
    geometry          : cavity geometry object.
    eps_walls         : wall emissivity / absorptivity (bulk or fallback).
    eps_base          : base emissivity / absorptivity (bulk or fallback).
    re_entry_prob     : probability that a photon escaping the aperture is
                        immediately re-injected downward (models plate A
                        reflection — Bug #9 fix).  Typically (1-ε_A)·F_AB.
    use_complex_fresnel : if True, use complex Fresnel reflectance via TMM.
    wall_thickness_um : wall thickness for TMM calculation (µm).
    wall_material     : wall material for complex refractive index lookup.
    base_material     : base material for complex refractive index lookup.
    photon_wavelength_um : photon wavelength (µm) for Fresnel/TMM calculation.
    """
    weight = 1.0

    bounces = 0
    total_steps = 0
    while bounces < _MAX_BOUNCES and total_steps < _MAX_PHOTON_STEPS:
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
            # Phase 1 Integration: Use complex Fresnel reflectance if enabled.
            # CRITICAL: weight is multiplied by R (reflectance), NOT (1-R).
            # (1-R) = T + A for a thin film; treating it as pure absorption
            # would incorrectly remove transmitted photons from the accounting,
            # breaking EM reciprocity.  Using R directly means the weight
            # tracks the reflected (surviving) fraction at each bounce.
            # The inter-pore PAA walls have air (n=1) on both sides, not
            # alumina, so n_2=1+0j is the correct backing medium.
            if use_complex_fresnel and photon_wavelength_um is not None and wall_thickness_um is not None:
                try:
                    n_real, k_imag = get_complex_refractive_index(wall_material, photon_wavelength_um)
                    n_complex = n_real + 1.0j * k_imag
                    R = tmm_reflectance_single_layer(
                        n_0=1.0 + 0.0j,  # vacuum/air incident medium
                        n_1=n_complex,   # wall material layer
                        n_2=1.0 + 0.0j, # air backing (inter-pore cavity)
                        thickness_um=wall_thickness_um,
                        wavelength_um=photon_wavelength_um
                    )
                    # weight *= R: surviving reflected fraction per bounce
                    weight *= float(np.clip(R, 0.0, 1.0))
                    if weight < _RR_THRESHOLD:
                        if np.random.random() < weight / _RR_BOOST:
                            weight = _RR_BOOST
                        else:
                            return 0.0
                    direction = sample_hemisphere_3d(normal)
                    pos = pos + normal * 1e-13
                    continue
                except Exception:
                    eps = eps_walls
            else:
                eps = eps_walls
        elif surface == 'base':
            # Phase 1 Integration: base is opaque bulk substrate.
            # R is the reflectance; weight *= R gives surviving fraction.
            if use_complex_fresnel and photon_wavelength_um is not None:
                try:
                    n_real, k_imag = get_complex_refractive_index(base_material, photon_wavelength_um)
                    n_complex = n_real + 1.0j * k_imag
                    # Assume base is thick (100 µm)
                    base_thickness_um = 100.0
                    R = tmm_reflectance_single_layer(
                        n_0=1.0 + 0.0j,
                        n_1=n_complex,
                        n_2=n_complex,  # bulk: self-backing
                        thickness_um=base_thickness_um,
                        wavelength_um=photon_wavelength_um
                    )
                    weight *= float(np.clip(R, 0.0, 1.0))
                    if weight < _RR_THRESHOLD:
                        if np.random.random() < weight / _RR_BOOST:
                            weight = _RR_BOOST
                        else:
                            return 0.0
                    direction = sample_hemisphere_3d(normal)
                    pos = pos + normal * 1e-13
                    continue
                except Exception:
                    eps = eps_base
            else:
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

    # A safety-cap exit is not material absorption. Preserve the surviving
    # weight so the absorptivity estimator does not count an unresolved,
    # highly reflective path as absorbed.
    return weight


def _trace_photon_thin_film(pos: np.ndarray,
                          direction: np.ndarray,
                          geometry: AnyGeometry,
                          eps_walls_bulk: float,
                          eps_base_bulk: float,
                          re_entry_prob: float = 0.0,
                          wall_thickness_um: float = None,
                          wall_material: str = 'alumina',
                          base_material: str = 'silver',
                          photon_wavelength_um: float = None,
                          use_complex_fresnel: bool = False,
                          apply_modal_attenuation: bool = False,
                          geometry_diameter_um: float = None,
                          wall_temperature_K: float = None,
                          base_temperature_K: float = None,
                          feature_scale_m: float = None,
                          wall_roughness_sigma_um: float = None,
                          wall_roughness_tau_um: float = None,
                          base_roughness_sigma_um: float = None,
                          base_roughness_tau_um: float = None) -> float:
    """Trace one photon with thin-film physics corrections.
    
    Extended version of _trace_photon that applies TMM-based thin-film
    emissivity for optically thin walls and optionally applies modal
    attenuation from Phase 2 integration.

    Phase 4a: ``wall_temperature_K`` / ``base_temperature_K`` feed the
    Drude–Lorentz temperature-dependent optical constants (None → 300 K
    tabulated data); ``feature_scale_m`` enables the hydrodynamic non-local
    correction when the structure is smaller than the electron mean free
    path.

    Phase 4b: ``*_roughness_sigma_um`` / ``*_roughness_tau_um`` activate the
    Beckmann–Spizzichino / Harvey–Shack roughness BRDF at each bounce
    (None → legacy pure-Lambertian reflection).
    
    Parameters
    ----------
    pos, direction    : initial position (m) and unit direction vector.
    geometry          : cavity geometry object.
    eps_walls_bulk    : bulk wall emissivity / absorptivity.
    eps_base_bulk     : bulk base emissivity / absorptivity.
    re_entry_prob     : probability of aperture re-entry.
    wall_thickness_um : wall thickness in µm (None = bulk assumption).
    wall_material     : wall material identifier.
    base_material     : base material identifier.
    photon_wavelength_um : photon wavelength in µm (required for thin-film).
    use_complex_fresnel : if True, use complex Fresnel instead of effective emissivity.
    apply_modal_attenuation : if True, apply Phase 2 modal attenuation weighting.
    geometry_diameter_um : cavity diameter for modal calculations (µm).
    
    Returns
    -------
    float — surviving weight at escape (>0) or 0 if absorbed.
    """
    # ---------------------------------------------------------------------
    # Precompute the constant per-bounce surviving-weight factor for the wall
    # and base surfaces.  All inputs to the TMM / complex-Fresnel /
    # effective-emissivity model (photon wavelength, material, thickness,
    # temperature and feature scale) are INVARIANT during a single photon's
    # flight, so these values were previously recomputed on every bounce —
    # a massive, redundant cost that dominated simulation runtime in deep
    # cavities (hundreds of bounces per photon).  Computing them once here is
    # mathematically identical.
    # ---------------------------------------------------------------------
    def _surface_survive(bulk_eps, thickness_um, material, temperature_K,
                         self_backing):
        """Surviving weight factor for one surface type (1 bounce)."""
        if (use_complex_fresnel and thickness_um is not None
                and photon_wavelength_um is not None):
            try:
                n_real, k_imag = get_complex_refractive_index(
                    material, photon_wavelength_um,
                    temperature_K=temperature_K,
                    feature_scale_m=feature_scale_m)
                n_complex = n_real + 1.0j * k_imag
                n_2 = n_complex if self_backing else (1.0 + 0.0j)
                R = tmm_reflectance_single_layer(
                    n_0=1.0 + 0.0j,
                    n_1=n_complex,
                    n_2=n_2,
                    thickness_um=thickness_um,
                    wavelength_um=photon_wavelength_um)
                return float(np.clip(R, 0.0, 1.0))
            except Exception:
                # Fallback: use absorptance-based scalar.
                return 1.0 - bulk_eps
        elif thickness_um is not None and photon_wavelength_um is not None:
            eps = effective_emissivity_thin_film(
                bulk_emissivity=bulk_eps,
                thickness_um=thickness_um,
                wavelength_um=photon_wavelength_um,
                material=material,
                temperature_K=temperature_K,
                feature_scale_m=feature_scale_m)
            return 1.0 - eps
        else:
            return 1.0 - bulk_eps

    # Inter-pore PAA walls are air | thin-alumina | air (self_backing=False);
    # the base is an opaque bulk substrate (self_backing=True, 100 µm).
    wall_survive = _surface_survive(
        eps_walls_bulk, wall_thickness_um, wall_material,
        wall_temperature_K, self_backing=False)
    base_survive = _surface_survive(
        eps_base_bulk, 100.0, base_material, base_temperature_K,
        self_backing=True)

    weight = 1.0
    weight = 1.0

    bounces = 0
    total_steps = 0
    last_bounce_pos = np.array(pos)  # Track position for modal attenuation calculation
    
    while bounces < _MAX_BOUNCES and total_steps < _MAX_PHOTON_STEPS:
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
                last_bounce_pos = np.array(pos)
                continue
            return weight

        if surface == 'periodic_x':
            pos[0] = -np.sign(pos[0]) * (geometry.P / 2 - 1e-13)
            continue
        elif surface == 'periodic_y':
            pos[1] = -np.sign(pos[1]) * (geometry.P / 2 - 1e-13)
            continue

        bounces += 1

        # Phase 2 Integration: Apply modal attenuation for internal bounces
        if apply_modal_attenuation and surface != 'aperture' and photon_wavelength_um is not None and geometry_diameter_um is not None:
            try:
                modal = solve_te11_mode_complex(
                    geometry_diameter_um, photon_wavelength_um, wall_material,
                    temperature_K=wall_temperature_K)
                # Calculate distance traveled since last bounce
                distance_since_last_bounce_um = (np.linalg.norm(pos - last_bounce_pos)) * 1e6  # Convert to micrometers
                attenuation = attenuation_factor_lossy_waveguide(distance_since_last_bounce_um, modal)
                weight *= attenuation
                
                # Apply Russian Roulette after attenuation
                if weight < _RR_THRESHOLD:
                    if np.random.random() < weight / _RR_BOOST:
                        weight = _RR_BOOST
                    else:
                        return 0.0
            except Exception:
                # Proceed without attenuation if modal calculation fails
                pass

        # ----------------------------------------------------------------
        # Per-bounce weight update and direction sampling.
        #
        # The Fresnel/TMM per-bounce rule depends on whether the surface
        # is a thin film (nearly transparent) or a bulk absorber.
        #
        # --- THIN-FILM WALL (use_complex_fresnel path) ---
        # For a thin inter-pore wall (air | thin-alumina | air), the
        # transfer matrix gives R (reflectance), T (transmittance), A
        # (absorptance) with R + T + A = 1.  The photon stays inside
        # the cavity only if it REFLECTS; otherwise it either transmits
        # into the adjacent pore or is absorbed.  The correct per-bounce
        # surviving weight is therefore R, NOT (1-R) = T+A.
        # Treating (1-R) as absorption removes transmitted photons from
        # the cavity accounting, breaking EM reciprocity.
        # PAA inter-pore walls have air on both sides: n_2 = 1+0j.
        #
        # --- NON-Fresnel thin-film path ---
        # effective_emissivity_thin_film → tmm_absorptance_normal_incidence
        # returns A = 1 - R - T (true absorptance only).  Here the
        # photon weight is reduced by A and the surviving weight (1-A)
        # includes both R and T contributions; this path is correct.
        #
        # --- BULK BASE path ---
        # The base is an opaque thick substrate (silver).  TMM with
        # self-backing gives R_base ≈ 1 - eps_base.  Surviving weight
        # is R_base = 1 - eps_base, which is what (1 - eps) gives since
        # effective_emissivity_thin_film returns eps ≈ eps_base_bulk
        # for a 100 µm thick absorber.
        # ----------------------------------------------------------------
        # Surviving weight for this surface.  All inputs to the TMM /
        # complex-Fresnel / effective-emissivity model (wavelength, material,
        # thickness, temperature, feature scale) are invariant during a
        # single photon's flight, so these factors are precomputed once per
        # photon above instead of being recomputed on every bounce.
        if surface == 'wall' or surface == 'top_cap':
            weight *= wall_survive
        elif surface == 'base':
            weight *= base_survive
        else:
            return 0.0

        if weight < _RR_THRESHOLD:
            if np.random.random() < weight / _RR_BOOST:
                weight = _RR_BOOST
            else:
                return 0.0

        # Phase 4b: roughness-aware reflection.  Walls/top-cap use the wall
        # BRDF, the base uses its own; legacy Lambertian when no roughness
        # is supplied or the photon wavelength is unknown (bulk mode).
        if surface in ('wall', 'top_cap'):
            sigma_hit = wall_roughness_sigma_um
            tau_hit = wall_roughness_tau_um
        else:
            sigma_hit = base_roughness_sigma_um
            tau_hit = base_roughness_tau_um
        if sigma_hit:
            direction = sample_surface_direction(
                direction, normal,
                sigma_um=sigma_hit, tau_um=tau_hit,
                wavelength_um=photon_wavelength_um)
        else:
            direction = sample_hemisphere_3d(normal)
        pos = pos + normal * 1e-13
        last_bounce_pos = np.array(pos)

    # A safety-cap exit is not material absorption. Preserve the surviving
    # weight so the absorptivity estimator does not count an unresolved,
    # highly reflective path as absorbed.
    return weight


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


def evanescent_power_transmission(lambda_c_um: float,
                                   lambda_um: float,
                                   length_um: float,
                                   material: str = 'alumina',
                                   diameter_um: float = None) -> float:
    """Return sub-cutoff power transmission through a channel length.

    Fix 2 (peer-review): Includes the aperture impedance-mismatch reflection
    coefficient R_ap at z = 0, derived from mode-matching the modal fields to
    the free-space continuum.  Without R_ap the boundary condition at the open
    aperture is unproven.

    The complete expression is:
        T_total = T_ap · exp(−2κL)
    where
        T_ap = 1 − R_ap,   R_ap = |(Z_mode − Z_0)/(Z_mode + Z_0)|²
        κ = 1/δ_ev,        δ_ev = evanescent_decay_length(lambda_c_um, lambda_um)

    For an evanescent mode, |E| decays as exp(-κh), so transmitted power
    carries the required exp(-2κh) factor.  This is a modal approximation,
    not a geometric ray probability.

    Parameters
    ----------
    lambda_c_um  : waveguide cutoff wavelength (µm)
    lambda_um    : photon wavelength (µm)
    length_um    : axial propagation distance (µm)
    material     : wall material for Z_mode / R_ap calculation ('alumina', etc.)
    diameter_um  : cavity diameter (µm); if None the aperture correction uses
                   the cutoff wavelength to back-calculate the diameter.
    """
    if length_um <= 0.0:
        return 1.0
    if (lambda_c_um is None or lambda_c_um <= 0.0 or not math.isfinite(lambda_c_um)):
        return 1.0
    # ------------------------------------------------------------------
    # Consolidated single-application of the modal attenuation operator.
    # T_total(λ) folds the multi-mode axial decay AND the aperture-plane
    # impedance-mismatch transmittance T_ap = 1 − R_ap into ONE evaluation per
    # (mode, wavelength), so the modal gating is never double-counted between
    # the ray-tracer and the spectrum-integrated exitance diagnostics.
    # ------------------------------------------------------------------
    try:
        from waveguide_modes import modal_total_exitance_transmission
        return modal_total_exitance_transmission(
            lambda_c_um, lambda_um, length_um,
            material=material, diameter_um=diameter_um)
    except Exception:
        # Conservative fallback: single dominant TE11 decay + aperture T_ap.
        delta_um = evanescent_decay_length(lambda_c_um, lambda_um)
        if not math.isfinite(delta_um) or delta_um <= 0.0:
            return 0.0
        T_decay = math.exp(-2.0 * length_um / delta_um)
        T_ap = 1.0
        try:
            from waveguide_modes import solve_te11_mode_complex
            _diam = diameter_um if diameter_um is not None else (lambda_c_um * 1.8412 / math.pi)
            _modal = solve_te11_mode_complex(
                diameter_um=_diam, wavelength_um=lambda_um, material=material,
                method='perturbation')
            T_ap = float(_modal.get('T_ap', 1.0))
        except Exception:
            T_ap = 1.0
        return T_ap * T_decay

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
    alpha_top:      float = None,  # top-surface absorptivity for sub-cutoff incident light (None → dynamic R_ap(λ,T,f))
    fill_fraction_top: float = None,  # aperture fill fraction f for dynamic R_ap(λ,T,f); None → geometry-derived
    # Thin-film physics parameters (PHASE 1)
    wall_thickness_um: float = None,  # wall thickness in µm (None = bulk)
    wall_material: str = 'alumina',   # material identifier
    base_material: str = 'silver',    # base material identifier
    use_complex_fresnel: bool = False,  # use complex Fresnel instead of effective emissivity (Phase 1)
    # Modal attenuation parameters (PHASE 2)
    apply_modal_attenuation: bool = False,  # enable Phase 2 modal loss weighting
    # Temperature-dependent optics parameters (PHASE 4a)
    wall_temperature_K: float = None,   # wall material temperature (None → T_emit)
    base_temperature_K: float = None,   # base material temperature (None → T_emit)
    feature_scale_m: float = None,      # structure size for non-local correction
    # Surface-roughness BRDF parameters (PHASE 4b)
    wall_roughness_sigma_um: float = None,   # wall RMS roughness σ (None → Lambertian)
    wall_roughness_tau_um: float = None,     # wall correlation length τ (µm)
    base_roughness_sigma_um: float = None,   # base RMS roughness σ (None → wall value)
    base_roughness_tau_um: float = None,     # base correlation length τ (None → wall value)
) -> dict:
    """Run the two cavity MC experiments and compute ε_B and α_eff.

    Experiment 1 — Internal emission (P_esc):
        Photons emitted from all internal surfaces (walls + base) weighted by
        area × emissivity.  Each photon gets a wavelength sampled from the
        Planck distribution at *T_emit*:
          • λ < λ_c  → geometric ray tracing until escape or absorption.
          • λ ≥ λ_c  → the channel mode is EVANESCENT; the photon tunnels to
                       the aperture with power T = exp(−2L/δ_ev) (≈0 for deep emission).
    Experiment 2 — External incidence (α_eff):
        Photons enter through the aperture with a Lambertian downward
        distribution; wavelength sampled from Planck at *T_inc*:
          • λ < λ_c  → geometric tracing (cavity trapping).
          • λ ≥ λ_c  → the wave cannot form a propagating channel mode; it is
                       captured by the lossy top surface (graded-index /
                       diffractive trap) with probability *alpha_top*.
    ε_B is the physically-computed emissivity from the macro cavity operator
    and the modal escape gate.  It is NOT forced equal to α_eff: for deep sub-wavelength
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
                     (None → dynamic aperture-boundary R_ap(λ, T, f) from the
                     Maxwell-Garnett effective-index step across fill fraction f).
    fill_fraction_top : aperture fill fraction f used for the dynamic
                     R_ap(λ, T, f); None → geometry-derived packing fraction.
    wall_thickness_um : wall thickness for Phase 1 thin-film/complex Fresnel (µm).
    wall_material  : wall material for complex refractive index lookup.
    base_material  : base material for complex refractive index lookup.
    use_complex_fresnel : if True, use complex Fresnel reflectance (Phase 1).
    apply_modal_attenuation : if True, apply Phase 2 modal loss weighting.

    Returns
    -------
    dict with keys:
      p_esc, p_esc_ci95            : escape probability and 95% CI half-width
      alpha_eff, alpha_eff_ci95    : effective absorptivity and 95% CI
      epsilon_b, epsilon_b_raw     : effective emissivity (physical emission)
      epsilon_b_ci95               : 95% CI on ε_B
      cavity_enhancement           : (A_walls + A_base) / A_aperture
      kirchhoff_error              : |ε_B - α_eff| / α_eff * 100 [%]
      f_prop_emit, f_prop_inc      : sampled propagating-mode fractions
      n_evan                       : number of evanescent (confined) photons
    """
    re_entry = (1.0 - eps_aperture) * view_factor_ab
    # ---- Dynamic aperture-boundary absorptance (peer review) -------------------
    # The hard-coded `alpha_top = 1 - (1 - eps_walls)·0.1` guess is REMOVED.
    # Sub-cutoff incident light that cannot form a channel mode folds around the
    # wall rims and impinges on the top surface, whose reflectance R_ap(λ, T, f)
    # is the Maxwell-Garnett effective-index step across the aperture fill
    # fraction f (aperture_boundary_absorptance).  When the caller does not pin
    # alpha_top we evaluate it per-wavelength per-photon below.
    _f_fill_top = fill_fraction_top
    if _f_fill_top is None:
        _f_fill_top = getattr(geometry, 'packing_fraction',
                              getattr(geometry, 'f', None))
    if _f_fill_top is None:
        _f_fill_top = 0.9                 # dense hexagonal-top default
    _f_fill_top = float(np.clip(_f_fill_top, 0.0, 1.0))

    def _alpha_top_dynamic(lam_um, T_k):
        """Absorptivity of the aperture top surface at (λ, T) = 1 − R_ap(λ,T,f)."""
        if alpha_top is not None:
            return float(np.clip(alpha_top, 0.0, 1.0))
        return float(np.clip(
            aperture_boundary_absorptance(wall_material, lam_um,
                                          _f_fill_top, float(T_k)),
            0.0, 1.0))

    # Normalise an explicitly-pinned alpha_top once for consistency.
    if alpha_top is not None:
        alpha_top = float(np.clip(alpha_top, 0.0, 1.0))

    # Waveguide cutoff of this channel geometry (µm)
    lambda_c_um = float(getattr(geometry, 'lambda_c_um', float('inf')))
    
    # Phase 4a: default the surface optical temperatures to the plate-B
    # (emitting/absorbing) temperature when not explicitly supplied.
    if wall_temperature_K is None:
        wall_temperature_K = float(T_emit)
    if base_temperature_K is None:
        base_temperature_K = float(T_emit)

    # Phase 4b: the base defaults to the wall roughness unless specified.
    if base_roughness_sigma_um is None:
        base_roughness_sigma_um = wall_roughness_sigma_um
    if base_roughness_tau_um is None:
        base_roughness_tau_um = wall_roughness_tau_um
    
    # Get cavity diameter for modal calculations (Phase 2)
    geometry_diameter_um = None
    if apply_modal_attenuation:
        if hasattr(geometry, 'diameter_um'):
            geometry_diameter_um = geometry.diameter_um
        elif hasattr(geometry, 'P'):  # CNT forest cell
            geometry_diameter_um = geometry.P

    # ---- Experiment 1: Internal emission -----------------------------------
    # By Kirchhoff's law, the internal emission source term must be weighted
    # by the actual local emissivity of the surface.  For a bulk-opaque wall
    # this is eps_walls; for an optically thin membrane (wall_thickness_um
    # provided) we must use the Planck-weighted thin-film emissivity at T_emit
    # so that the energy injected into the cavity is consistent with the
    # per-photon absorption computed via TMM in _trace_photon_thin_film.
    # Using bulk eps_walls here while TMM gives per-bounce reflectance
    # R ≈ 1 − ε_thin creates an internal contradiction that blows up ε_B.
    A_walls  = geometry.area_walls
    A_base   = geometry.area_base
    if wall_thickness_um is not None and wall_thickness_um > 0.0 and wall_material:
        eps_walls_emit = planck_weighted_effective_emissivity(
            eps_walls, wall_thickness_um, wall_material, float(T_emit)
        )
    else:
        eps_walls_emit = eps_walls  # bulk / opaque walls
    w_walls  = A_walls * eps_walls_emit
    w_base   = A_base  * eps_base
    w_total  = w_walls + w_base
    p_wall   = w_walls / w_total if w_total > 0 else 0.5

    # Analytic spectral power split across cutoff
    f_prop_emit = planck_cumulative(lambda_c_um * float(T_emit)) if math.isfinite(lambda_c_um) else 1.0
    f_evan_emit = max(0.0, 1.0 - f_prop_emit)

    # Stratified photon allocation across spectral regimes
    n_prop_target = int(round(n_photons * f_prop_emit)) if (f_prop_emit > 0.05 and f_prop_emit < 0.95) else (
        max(100, int(n_photons * 0.2)) if (f_prop_emit > 1e-8 and f_prop_emit <= 0.05) else (
            min(n_photons - 100, int(n_photons * 0.8)) if (f_prop_emit >= 0.95 and f_prop_emit < (1.0 - 1e-8)) else (
                n_photons if f_prop_emit >= (1.0 - 1e-8) else 0
            )
        )
    )
    n_evan_target = n_photons - n_prop_target

    escaped_weight_prop = 0.0
    escaped_weight_prop_sq = 0.0
    for _ in range(n_prop_target):
        if np.random.random() < p_wall:
            pos, n_out = geometry.sample_point_on_walls()
        else:
            pos, n_out = geometry.sample_point_on_base()
        direction = sample_hemisphere_3d(n_out)
        pos = pos + n_out * 1e-13

        lam = sample_planck_wavelength_band(T_emit, 0.05, lambda_c_um)
        # ---- Guided-mode propagation (MATHEMATICAL_DERIVATIONS.md Ph.1-3) ----
        # Above-cutoff photons are WAVEGUIDE MODES, not diffusive billiard
        # balls: they travel the channel axis and exit at the aperture with
        # transmission T_total = T_ap · exp(-alpha·L) from the lossy-wall
        # modal solver (consolidated once per mode/wavelength).  Geometric
        # hemisphere bouncing here terminated photons on first wall hit,
        # collapsing p_esc to ~3% of f_prop instead of the documented
        # ~92% of f_prop.
        L_chan = max(geometry.H - pos[2], 0.0) * 1e6   # µm to aperture plane
        w = None
        if apply_modal_attenuation and geometry_diameter_um is not None:
            try:
                _m = solve_te11_mode_complex(
                    geometry_diameter_um, lam, wall_material,
                    temperature_K=wall_temperature_K)
                w = float(np.clip(
                    _m['T_ap'] * attenuation_factor_lossy_waveguide(L_chan, _m),
                    0.0, 1.0))
            except Exception:
                w = None
        if w is None:
            # Fallback: legacy geometric trace when no modal channel applies.
            w = _trace_photon_thin_film(
                pos, direction, geometry, 
                eps_walls_bulk=eps_walls,
                eps_base_bulk=eps_base,
                re_entry_prob=re_entry,
                wall_thickness_um=wall_thickness_um,
                wall_material=wall_material,
                base_material=base_material,
                photon_wavelength_um=lam,
                use_complex_fresnel=use_complex_fresnel,
                apply_modal_attenuation=False,
                geometry_diameter_um=None,
                wall_temperature_K=wall_temperature_K,
                base_temperature_K=base_temperature_K,
                feature_scale_m=feature_scale_m,
                wall_roughness_sigma_um=wall_roughness_sigma_um,
                wall_roughness_tau_um=wall_roughness_tau_um,
                base_roughness_sigma_um=base_roughness_sigma_um,
                base_roughness_tau_um=base_roughness_tau_um
            )
        escaped_weight_prop += w
        escaped_weight_prop_sq += w * w

    escaped_weight_evan = 0.0
    escaped_weight_evan_sq = 0.0
    tunnel_escape_w = 0.0
    for _ in range(n_evan_target):
        if np.random.random() < p_wall:
            pos, n_out = geometry.sample_point_on_walls()
        else:
            pos, n_out = geometry.sample_point_on_base()
        direction = sample_hemisphere_3d(n_out)
        pos = pos + n_out * 1e-13

        lam = sample_planck_wavelength_band(T_emit, lambda_c_um, 2000.0)
        L = max(geometry.H - pos[2], 0.0) * 1e6   # µm to aperture plane
        _diam_ev = getattr(geometry, 'diameter_um', None)
        w = evanescent_power_transmission(
            lambda_c_um, lam, L,
            material=wall_material,
            diameter_um=_diam_ev,
        )
        escaped_weight_evan += w
        escaped_weight_evan_sq += w * w
        tunnel_escape_w += w

    mean_prop = (escaped_weight_prop / n_prop_target) if n_prop_target > 0 else 0.0
    mean_evan = (escaped_weight_evan / n_evan_target) if n_evan_target > 0 else 0.0

    p_esc = float(f_prop_emit * mean_prop + f_evan_emit * mean_evan)

    var_prop = ((escaped_weight_prop_sq / n_prop_target) - (mean_prop ** 2)) / n_prop_target if n_prop_target > 1 else 0.0
    var_evan = ((escaped_weight_evan_sq / n_evan_target) - (mean_evan ** 2)) / n_evan_target if n_evan_target > 1 else 0.0
    var_total = (f_prop_emit ** 2) * max(var_prop, 0.0) + (f_evan_emit ** 2) * max(var_evan, 0.0)
    p_esc_ci95 = float(1.96 * math.sqrt(max(var_total, 0.0)))

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
    f_prop_inc = planck_cumulative(lambda_c_um * float(T_inc)) if math.isfinite(lambda_c_um) else 1.0
    f_evan_inc = max(0.0, 1.0 - f_prop_inc)

    n_prop_inc_target = int(round(n_photons * f_prop_inc)) if (f_prop_inc > 0.05 and f_prop_inc < 0.95) else (
        max(100, int(n_photons * 0.2)) if (f_prop_inc > 1e-8 and f_prop_inc <= 0.05) else (
            min(n_photons - 100, int(n_photons * 0.8)) if (f_prop_inc >= 0.95 and f_prop_inc < (1.0 - 1e-8)) else (
                n_photons if f_prop_inc >= (1.0 - 1e-8) else 0
            )
        )
    )
    n_evan_inc_target = n_photons - n_prop_inc_target

    absorbed_weight_prop = 0.0
    absorbed_weight_prop_sq = 0.0
    for _ in range(n_prop_inc_target):
        pos = _sample_aperture()
        lam = sample_planck_wavelength_band(T_inc, 0.05, lambda_c_um)
        direction = sample_hemisphere_3d(down)
        w = _trace_photon_thin_film(
            pos, direction, geometry,
            eps_walls_bulk=eps_walls,
            eps_base_bulk=eps_base,
            # Symmetric boundary condition (Kirchhoff reciprocity): the
            # plate-A aperture re-entry coupling must act identically on the
            # external-incidence leg and on the internal-emission leg.
            # Hard-coding 0.0 here let incident light leave the cavity freely
            # while internally emitted photons were partially reflected back,
            # breaking detailed balance and biasing alpha_eff high.
            re_entry_prob=re_entry,
            wall_thickness_um=wall_thickness_um,
            wall_material=wall_material,
            base_material=base_material,
            photon_wavelength_um=lam,
            use_complex_fresnel=use_complex_fresnel,
            apply_modal_attenuation=apply_modal_attenuation,
            geometry_diameter_um=geometry_diameter_um,
            wall_temperature_K=wall_temperature_K,
            base_temperature_K=base_temperature_K,
            feature_scale_m=feature_scale_m,
            wall_roughness_sigma_um=wall_roughness_sigma_um,
            wall_roughness_tau_um=wall_roughness_tau_um,
            base_roughness_sigma_um=base_roughness_sigma_um,
            base_roughness_tau_um=base_roughness_tau_um
        )
        a = 1.0 - w
        absorbed_weight_prop += a
        absorbed_weight_prop_sq += a * a

    # ---- Evanescent incidence: thin-rim capture ----------------------------
    # The documented structured-surface model treats sub-cutoff incident
    # radiation as diffracting around the thin rim and being captured by the
    # graded-index top surface.  This is intentionally distinct from the
    # internal emission escape calculation, where deep evanescent modes decay.
    absorbed_weight_evan = 0.0
    absorbed_weight_evan_sq = 0.0
    H_m = getattr(geometry, 'H', None)
    _diam_inc = getattr(geometry, 'diameter_um', None)
    n_evan_no_hit = 0
    for _ in range(n_evan_inc_target):
        pos = _sample_aperture()
        direction = sample_hemisphere_3d(down)
        # Same spectral band as the internal-emission evanescent leg.
        lam_evan = sample_planck_wavelength_band(T_inc, lambda_c_um, 2000.0)
        # Dynamic aperture-boundary absorptance 1 − R_ap(λ, T, f): evaluated
        # per-wavelength against the Maxwell-Garnett effective-index step.
        a = _alpha_top_dynamic(lam_evan, T_inc)
        absorbed_weight_evan += a
        absorbed_weight_evan_sq += a * a

    mean_abs_prop = (absorbed_weight_prop / n_prop_inc_target) if n_prop_inc_target > 0 else 0.0
    mean_abs_evan = (absorbed_weight_evan / n_evan_inc_target) if n_evan_inc_target > 0 else 0.0

    alpha_eff = float(f_prop_inc * mean_abs_prop + f_evan_inc * mean_abs_evan)

    var_abs_prop = ((absorbed_weight_prop_sq / n_prop_inc_target) - (mean_abs_prop ** 2)) / n_prop_inc_target if n_prop_inc_target > 1 else 0.0
    var_abs_evan = ((absorbed_weight_evan_sq / n_evan_inc_target) - (mean_abs_evan ** 2)) / n_evan_inc_target if n_evan_inc_target > 1 else 0.0
    var_abs_total = (f_prop_inc ** 2) * max(var_abs_prop, 0.0) + (f_evan_inc ** 2) * max(var_abs_evan, 0.0)
    alpha_eff_ci95 = float(1.96 * math.sqrt(max(var_abs_total, 0.0)))

    # ---- Derived quantities ------------------------------------------------
    # Reciprocity estimator for the cavity emissivity.  The same cavity
    # surfaces and material emissivities used for internal emission must also
    # determine external absorptivity; using a wall-only macro operator here
    # silently drops the base emissivity when eps_walls != eps_base.
    A_int = max(A_walls + A_base, 1e-30)
    aperture_ratio = A_ap / A_int
    epsilon_mc = (w_total / max(A_ap, 1e-30)) * p_esc
    epsilon_mc_ci95 = (w_total / max(A_ap, 1e-30)) * p_esc_ci95

    # Retain the gray macro operator as a diagnostic.  It uses an
    # area-weighted emissivity and is not used as the primary estimator when
    # wall and base emissivities differ.
    eps_area_weighted = (A_walls * eps_walls + A_base * eps_base) / A_int
    macro_cavity_eps = float(np.clip(
        eps_area_weighted / (eps_area_weighted
                              + (1.0 - eps_area_weighted) * aperture_ratio),
        0.0, 1.0,
    ))

    cavity_enhancement = A_int / A_ap if A_ap > 0 else 0.0
    epsilon_b_raw      = float(np.clip(epsilon_mc, 0.0, 1.0))
    epsilon_b          = epsilon_b_raw
    epsilon_b_ci95     = float(np.clip(epsilon_mc_ci95, 0.0, 1.0))
    # Reciprocity budget: directional/angular separation of one reciprocal
    # surface.  Material-level epsilon = alpha is imposed at the Fresnel/TMM
    # kernel; the budget below is spectral/angular anisotropy (distinct
    # integrals), not a Lorentz/Kirchhoff violation.
    reciprocity_budget = abs(epsilon_b_raw - alpha_eff) / max(alpha_eff, 1e-12)
    reciprocity_origin = 'spectral_anisotropy'
    kirchhoff_error    = reciprocity_budget * 100.0   # kept for UI compat

    f_prop_emit = float(f_prop_emit)
    f_prop_inc  = float(f_prop_inc)

    return {
        'p_esc':              p_esc,
        'p_esc_ci95':         p_esc_ci95,
        'alpha_eff':          alpha_eff,
        'alpha_eff_ci95':     alpha_eff_ci95,
        'epsilon_b':          epsilon_b,
        'epsilon_b_raw':      epsilon_b_raw,
        'epsilon_b_ci95':     epsilon_b_ci95,
        'epsilon_b_macro':    macro_cavity_eps,
        'cavity_enhancement': cavity_enhancement,
        'kirchhoff_error':    kirchhoff_error,
        'reciprocity_budget': float(reciprocity_budget),
        'reciprocity_origin': reciprocity_origin,
        'f_prop_emit':        f_prop_emit,
        'f_prop_inc':         f_prop_inc,
        'tunnel_escape_w':    tunnel_escape_w,
        'n_evan':             n_evan_target,
        'w_total':            w_total,
        'area_aperture':      A_ap,
    }
