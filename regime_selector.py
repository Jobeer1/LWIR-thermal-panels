"""
regime_selector.py — Physics Orchestrator
=========================================

Pre-flight dimensionless-regime selector evaluated BEFORE any optical solver,
ray tracer, or radiosity calculation runs.  Computes non-dimensional geometric /
optical ratios (lambda/D, lambda/H, lambda/t_wall, d_gap/lambda, H/D, t/delta),
dynamically selects the valid physics regimes, deducts an audit-transparent
confidence score, and exposes strict solver-enforcement flags consumed by
``simulator.run_simulation``.
"""

from dataclasses import dataclass, field
from typing import List
import numpy as np


@dataclass
class PhysicsRegime:
    """Selected physics regimes + enforcement flags for one configuration."""
    # Dimensionless Metrics
    lambda_over_diameter: float
    lambda_over_height: float
    lambda_over_thickness: float
    gap_over_lambda: float
    height_over_diameter: float
    t_over_delta: float  # Wall thickness over skin depth

    # Selected Regimes
    geometry_regime: str       # "RAY", "FULL_WAVE", "EMT"
    wall_regime: str           # "BULK", "THIN_FILM", "MEMBRANE"
    cavity_regime: str         # "CLASSICAL", "TRANSITIONAL", "CUTOFF_DOMINATED"
    heat_transfer_regime: str  # "FAR_FIELD", "NEAR_FIELD"

    # Solver Enforcement Flags
    ray_valid: bool
    rcwa_recommended: bool
    thinfilm_required: bool
    modal_cutoff_required: bool
    nearfield_required: bool

    # Confidence Score & Diagnostics
    confidence: float
    regime_warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """JSON-able provenance payload for results / UI reports."""
        return {
            'dimensionless': {
                'lambda_over_diameter':  self.lambda_over_diameter,
                'lambda_over_height':    self.lambda_over_height,
                'lambda_over_thickness': self.lambda_over_thickness,
                'gap_over_lambda':       self.gap_over_lambda,
                'height_over_diameter':  self.height_over_diameter,
                't_over_delta':          self.t_over_delta,
            },
            'geometry_regime':      self.geometry_regime,
            'wall_regime':          self.wall_regime,
            'cavity_regime':        self.cavity_regime,
            'heat_transfer_regime': self.heat_transfer_regime,
            'ray_valid':            self.ray_valid,
            'rcwa_recommended':     self.rcwa_recommended,
            'thinfilm_required':    self.thinfilm_required,
            'modal_cutoff_required': self.modal_cutoff_required,
            'nearfield_required':   self.nearfield_required,
            'confidence':           self.confidence,
            'regime_warnings':      list(self.regime_warnings),
        }


def select_physics_regime(
    T_emit_K: float,
    diameter_um: float,
    height_um: float,
    wall_thickness_um: float,
    gap_um: float,
    k_extinction: float,
    lambda_cutoff_um: float = None,
    rcwa_available: bool = False
) -> PhysicsRegime:

    warnings = []
    confidence = 100.0

    # Wien's Displacement Law for peak emission wavelength
    lambda_peak_um = 2898.0 / max(float(T_emit_K), 1e-6)

    # 1. Calculate Dimensionless Ratios
    l_over_d = lambda_peak_um / max(diameter_um, 1e-6)
    l_over_h = lambda_peak_um / max(height_um, 1e-6)
    l_over_t = lambda_peak_um / max(wall_thickness_um, 1e-6)
    gap_over_l = gap_um / max(lambda_peak_um, 1e-6)
    h_over_d = height_um / max(diameter_um, 1e-6)

    # Optical skin depth delta = lambda / (4 * pi * k)
    skin_depth_um = (lambda_peak_um / (4.0 * np.pi * k_extinction)) if k_extinction > 1e-6 else 1e6
    t_over_delta = wall_thickness_um / skin_depth_um

    # 2. Geometry Regime & Ray Validity (matches simulator.py dispatcher).
    #    The three windows are:
    #      λ/D < 0.2         → RAY        (D > 5λ, geometric optics valid)
    #      0.2 ≤ λ/D ≤ 5.0   → FULL_WAVE  (diffraction / resonance; ray invalid)
    #      λ/D > 5.0         → EMT        (D ≪ λ, homogenised Maxwell-Garnett)
    if l_over_d < 0.2:
        geometry_regime = "RAY"
        ray_valid = True
        rcwa_rec = False
    elif 0.2 <= l_over_d <= 5.0:
        geometry_regime = "FULL_WAVE"
        ray_valid = False   # HARD ENFORCEMENT: pure geometric ray tracing invalid
        rcwa_rec = True
        confidence -= 25.0
        warnings.append(
            "Diffractive / sub-wavelength cavity (0.2 <= λ/D <= 5): geometric "
            "ray tracing disabled; routing to the waveguide modal solver (or "
            "full-wave RCWA when available).")
    else:
        geometry_regime = "EMT"
        ray_valid = False
        rcwa_rec = False
        confidence -= 25.0
        warnings.append(
            "Deep sub-wavelength cavity (λ/D > 5): homogenised "
            "Maxwell-Garnett effective medium; ray tracing bypassed.")

    # 3. Wall Regime Selection
    if t_over_delta > 3.0:
        wall_regime = "BULK"
        thinfilm_req = False
    elif 0.1 < t_over_delta <= 3.0:
        wall_regime = "THIN_FILM"
        thinfilm_req = True
    else:
        wall_regime = "MEMBRANE"
        thinfilm_req = True
        confidence -= 10.0
        warnings.append("Optically thin membrane limit (t <= 0.1δ): Beer-Lambert & substrate penetration dominates.")

    # 4. Cavity Cutoff Regime (diagnostic: λ_peak vs λ_c)
    if lambda_cutoff_um is not None and lambda_cutoff_um > 0:
        if lambda_peak_um < lambda_cutoff_um:
            cavity_regime = "CLASSICAL"
        elif 0.5 * lambda_cutoff_um <= lambda_peak_um <= 2.0 * lambda_cutoff_um:
            cavity_regime = "TRANSITIONAL"
        else:
            cavity_regime = "CUTOFF_DOMINATED"
            confidence -= 15.0
            warnings.append(
                "Cutoff dominated regime (λ > 2λ_c): evanescent modal exitance "
                "controls emission.")
    else:
        cavity_regime = "CLASSICAL"

    # The waveguide modal-transmission operator is required whenever the cavity
    # is a sub-wavelength waveguide (FULL_WAVE / EMT geometry): BOTH propagating
    # modes (guided T_ap·exp(-αL)) and evanescent modes (exp(-2L/δ)) need it.
    # "CLASSICAL" only means the Planck peak sits below λ_c — it does NOT make
    # geometric ray tracing valid for D ~ λ, so modal_req must follow the
    # GEOMETRY regime, not the cavity regime.  (At 3000 K, λ/D ≈ 0.24 → FULL_WAVE
    # yet λ_peak < λ_c; disabling the modal solver there collapsed P_esc to 0.29%
    # instead of the guided-mode ~92% and inflated the CI to ±42%.)
    modal_req = (geometry_regime in ("FULL_WAVE", "EMT"))

    # 5. Heat Transfer Regime (Near Field)
    nearfield_threshold = lambda_peak_um / (2.0 * np.pi)
    if gap_um < nearfield_threshold:
        heat_transfer_regime = "NEAR_FIELD"
        nearfield_req = True
    else:
        heat_transfer_regime = "FAR_FIELD"
        nearfield_req = False

    # 6. Additional Confidence Penalty Audits
    if T_emit_K > 1500.0:
        confidence -= 15.0
        warnings.append(f"High temperature ({T_emit_K} K): Dielectric function n(λ,T)+ik(λ,T) relies on extrapolation.")

    if rcwa_rec and not rcwa_available:
        confidence -= 20.0
        warnings.append(
            "Full-wave RCWA recommended by scale, but solver unavailable; "
            "falling back to the waveguide modal solver.")

    return PhysicsRegime(
        lambda_over_diameter=l_over_d,
        lambda_over_height=l_over_h,
        lambda_over_thickness=l_over_t,
        gap_over_lambda=gap_over_l,
        height_over_diameter=h_over_d,
        t_over_delta=t_over_delta,
        geometry_regime=geometry_regime,
        wall_regime=wall_regime,
        cavity_regime=cavity_regime,
        heat_transfer_regime=heat_transfer_regime,
        ray_valid=ray_valid,
        rcwa_recommended=rcwa_rec,
        thinfilm_required=thinfilm_req,
        modal_cutoff_required=modal_req,
        nearfield_required=nearfield_req,
        confidence=max(0.0, confidence),
        regime_warnings=warnings
    )


def _physics_regime_to_dict(self) -> Dict[str, Any]:
    """JSON-serialisable provenance payload for UI / reports."""
    return {
        'lambda_over_diameter':   float(self.lambda_over_diameter),
        'lambda_over_height':     float(self.lambda_over_height),
        'lambda_over_thickness':  float(self.lambda_over_thickness),
        'gap_over_lambda':        float(self.gap_over_lambda),
        'height_over_diameter':   float(self.height_over_diameter),
        't_over_delta':           float(self.t_over_delta),
        'geometry_regime':        self.geometry_regime,
        'wall_regime':            self.wall_regime,
        'cavity_regime':          self.cavity_regime,
        'heat_transfer_regime':   self.heat_transfer_regime,
        'ray_valid':              bool(self.ray_valid),
        'rcwa_recommended':       bool(self.rcwa_recommended),
        'thinfilm_required':      bool(self.thinfilm_required),
        'modal_cutoff_required':  bool(self.modal_cutoff_required),
        'nearfield_required':     bool(self.nearfield_required),
        'confidence':             float(self.confidence),
        'regime_warnings':        list(self.regime_warnings),
    }


PhysicsRegime.to_dict = _physics_regime_to_dict
