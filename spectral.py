"""
spectral.py — Spectral band model for gray-body approximation upgrade.

Provides a 5-band spectral model with Planck-weighted effective emissivities
for CNT forests and metallic (Ag) bases.

Literature sources for band emissivities:
  CNT forest: Mizuno et al. PNAS 2009; Yang et al. NL 2008
  Silver:     Palik, Handbook of Optical Constants, 1985
"""

import numpy as np
from sampling import planck_band_fraction

# ---------------------------------------------------------------------------
# Spectral band definitions  [lam_lo_um, lam_hi_um]
# ---------------------------------------------------------------------------
#
#  Band 1:  0.3 –  2   µm  (near-IR / visible)
#  Band 2:  2   –  5   µm  (mid-IR, CNT forest peak efficiency)
#  Band 3:  5   – 10   µm  (thermal IR, peak for T~300-600 K)
#  Band 4: 10   – 25   µm  (far-IR)
#  Band 5: 25   – 100  µm  (very far-IR / sub-mm)

BANDS = [
    (0.3,   2.0),
    (2.0,   5.0),
    (5.0,  10.0),
    (10.0, 25.0),
    (25.0, 100.0),
]

# ---------------------------------------------------------------------------
# Spectral emissivity data for common materials
# ---------------------------------------------------------------------------
# Each entry: list of (ε_band1, ε_band2, ε_band3, ε_band4, ε_band5)

MATERIAL_EMISSIVITY = {
    # VACNT forest (vertically aligned multi-wall CNT, ~10 µm pitch array)
    # Refs: Mizuno 2009 (ε > 0.98 in mid-IR), Yang 2008 (near-IR ~0.95)
    'cnt_forest': [0.95, 0.99, 0.99, 0.98, 0.97],

    # Silver (bulk, polished)
    # High reflectivity in IR; lower in visible/near-IR
    # Refs: Palik 1985
    'silver': [0.02, 0.02, 0.03, 0.03, 0.04],

    # Gold (polished)
    'gold': [0.02, 0.02, 0.03, 0.03, 0.04],

    # Oxidised silver / contaminated (common in ambient CNT experiments)
    'silver_oxidised': [0.05, 0.05, 0.08, 0.10, 0.12],

    # Generic metal (steel, etc.)
    'metal_generic': [0.15, 0.12, 0.10, 0.10, 0.10],

    # Gray body (flat spectrum — reduces to gray-body limit)
    'gray_0.98': [0.98, 0.98, 0.98, 0.98, 0.98],
    'gray_0.02': [0.02, 0.02, 0.02, 0.02, 0.02],
}


def planck_weighted_emissivity(emissivity_bands: list, T_K: float) -> float:
    """Compute the Planck-weighted total emissivity from per-band values.

    ε_total(T) = Σ_i ε_i · F(λ_{i,lo}, λ_{i,hi}, T)

    where F is the fractional blackbody function.

    Parameters
    ----------
    emissivity_bands : list of 5 floats, one per spectral band.
    T_K              : surface temperature in Kelvin.

    Returns
    -------
    float — Planck-weighted total emissivity in [0, 1].
    """
    if T_K <= 0:
        return float(np.mean(emissivity_bands))
    total = 0.0
    for (lam_lo, lam_hi), eps in zip(BANDS, emissivity_bands):
        f = planck_band_fraction(lam_lo, lam_hi, T_K)
        total += eps * f
    # Remaining fraction beyond last band (very far-IR) gets last band emissivity
    f_tail = 1.0 - sum(planck_band_fraction(lo, hi, T_K) for lo, hi in BANDS)
    total += emissivity_bands[-1] * max(0.0, f_tail)
    return float(np.clip(total, 0.0, 1.0))


def effective_emissivity_pair(
    mat_a: str | None, eps_a_scalar: float, T_a: float,
    mat_b: str | None, eps_b_scalar: float, T_b: float,
) -> dict:
    """Return Planck-weighted emissivities for plate A and B.

    If mat_* is a known material key, use literature spectral data.
    Otherwise fall back to the provided scalar (gray-body) value.

    Returns
    -------
    dict with keys: eps_a_spectral, eps_b_spectral, spectral_correction_pct
    """
    def _get(mat, scalar, T):
        if mat and mat in MATERIAL_EMISSIVITY:
            return planck_weighted_emissivity(MATERIAL_EMISSIVITY[mat], T)
        # Build a flat spectrum from the scalar
        return planck_weighted_emissivity([scalar] * 5, T)

    eps_a = _get(mat_a, eps_a_scalar, T_a)
    eps_b = _get(mat_b, eps_b_scalar, T_b)
    correction = abs(eps_a - eps_a_scalar) / max(eps_a_scalar, 1e-9) * 100.0

    return {
        'eps_a_spectral': eps_a,
        'eps_b_spectral': eps_b,
        'spectral_correction_pct': correction,
    }
