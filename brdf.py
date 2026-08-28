"""
brdf.py — Surface-roughness BRDF for the Monte Carlo ray tracer (Phase 4b).

Replaces the purely Lambertian (cosine-law) bounce assumption with a hybrid
specular / Beckmann-lobe / diffuse reflection model driven by the surface RMS
roughness σ and spatial correlation length τ.

Key physics
-----------
1. Beckmann–Spizzichino coherent (specular) fraction.  A wave reflected from
   a randomly rough surface retains a coherent specular component carrying
   the power fraction

       p_spec(λ, θᵢ) = exp[ −(4π σ cosθᵢ / λ)² ]

   (Beckmann & Spizzichino, *The Scattering of Electromagnetic Waves from
   Rough Surfaces*, 1987; Ogilvy 1991).  σ→0 or λ→∞ gives a perfect mirror;
   σ ≫ λ gives a fully incoherent reflection.

2. Total integrated scatter (TIS).  The incoherent complement is
       TIS = 1 − exp[−(4πσ/λ)²]
   which is the standard smooth-surface hemispherical scatter result
   (Stover, *Optical Scattering*, 2012).

3. Beckmann Gaussian-facet near-specular lobe.  The incoherent power that
   scatters forward is distributed around the mirror direction with a slope
   pdf p(s) ∝ exp(−s²/2m²), s = tanθ, with rms slope

       m = √2 · σ/τ

   for Gaussian height statistics with correlation length τ.

4. Harvey–Shack aureole.  The remaining incoherent power goes into a broad
   diffuse background characterised by the Harvey ABg surface transfer
   function  BSDF(θ) = A / (B² + θ²)^(g/2)  (Harvey 1976).  Here the
   background is sampled as Lambertian, with its weight growing as the
   surface roughens — the aureole limit of the Harvey–Shack model.

5. Energy conservation.  The Fresnel/TMM reflected power R is split between
   the three lobes *in the directional pdf only*; the photon weight carries
   the full R.  The weight-based absorption bookkeeping of the ray tracer is
   therefore untouched.

Conventions
-----------
* Lengths in micrometres (σ, τ, λ), consistent with the rest of the app.
* ``incident_dir`` points INTO the surface (dot(incident, normal) < 0);
  returned directions point AWAY (dot(dir, normal) ≥ 0).
* ``sigma_um=None`` or ≤ 0 selects the legacy pure-Lambertian path exactly
  (backward compatible).  σ → 0⁺ with the model enabled approaches a mirror.

References
----------
* Beckmann & Spizzichino (1987), Ch. 5 — coherent fraction, facet slopes.
* Harvey (1976); Stover (2012) — ABg BSDF, total integrated scatter.
* Siegel & Howell, *Thermal Radiation Heat Transfer* — mixed
  specular-diffuse reflection models for enclosure Monte Carlo.
"""

import math

import numpy as np

try:
    from sampling import sample_hemisphere_3d
except ImportError:  # pragma: no cover - standalone use
    raise ImportError('brdf.py requires sampling.py (Malley hemisphere sampler)')

# Slope-sampling guards: tan(θ) beyond ~15 is ≈ 86°, effectively grazing.
_MAX_SLOPE = 15.0
_MIN_SLOPE_RMS = 1e-4
_MAX_SLOPE_RMS = 5.0


# ---------------------------------------------------------------------------
# Analytic BRDF quantities
# ---------------------------------------------------------------------------

def beckmann_specular_fraction(sigma_um: float, wavelength_um: float,
                               cos_theta_i: float) -> float:
    """Coherent specular power fraction p_spec = exp[−(4πσcosθᵢ/λ)²].

    Parameters
    ----------
    sigma_um       : RMS surface roughness σ (µm).
    wavelength_um  : vacuum wavelength λ (µm).
    cos_theta_i    : cosine of the incidence angle from the surface normal.

    Returns
    -------
    float in [0, 1].  1.0 = perfect mirror, 0.0 = fully incoherent.
    """
    if sigma_um is None or sigma_um <= 0 or wavelength_um is None or wavelength_um <= 0:
        return 1.0
    g = 4.0 * math.pi * float(sigma_um) * max(0.0, min(1.0, float(cos_theta_i))) / float(wavelength_um)
    return math.exp(-g * g)


def total_integrated_scatter(sigma_um: float, wavelength_um: float) -> float:
    """Hemispherical incoherent scatter fraction TIS = 1 − exp[−(4πσ/λ)²]."""
    if sigma_um is None or sigma_um <= 0 or wavelength_um is None or wavelength_um <= 0:
        return 0.0
    g = 4.0 * math.pi * float(sigma_um) / float(wavelength_um)
    return 1.0 - math.exp(-g * g)


def rms_slope(sigma_um: float, tau_um: float) -> float:
    """RMS surface slope m = √2·σ/τ for Gaussian roughness statistics.

    Clamped to [_MIN_SLOPE_RMS, _MAX_SLOPE_RMS] for numerical sampling.
    """
    if tau_um is None or tau_um <= 0 or sigma_um is None or sigma_um <= 0:
        return _MAX_SLOPE_RMS
    m = math.sqrt(2.0) * float(sigma_um) / float(tau_um)
    return min(max(m, _MIN_SLOPE_RMS), _MAX_SLOPE_RMS)


def harvey_shack_bsdf(theta_rad: float, A: float = 1.0,
                      B: float = 0.005, g: float = 2.0) -> float:
    """Harvey ABg surface transfer function BSDF(θ) = A/(B² + θ²)^(g/2).

    Peak at θ = 0 (specular direction), power-law aureole decay.  Provided
    for BRDF evaluation / validation; the Monte Carlo sampler uses the
    equivalent mixture described in the module docstring.
    """
    t = abs(float(theta_rad))
    return A / ((B * B + t * t) ** (g / 2.0))


# ---------------------------------------------------------------------------
# Monte Carlo direction sampling
# ---------------------------------------------------------------------------

def _unit(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v)
    if norm < 1e-15:
        return np.array([0.0, 0.0, 1.0])
    return v / norm


def _mirror_direction(incident_dir: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Specular reflection d − 2(d·n)n, normalised."""
    d = _unit(incident_dir)
    n = _unit(normal)
    return _unit(d - 2.0 * float(np.dot(d, n)) * n)


def sample_beckmann_lobe(mirror_dir: np.ndarray, normal: np.ndarray,
                         slope_rms: float) -> np.ndarray:
    """Sample the Gaussian-facet lobe around the mirror direction.

    Slope s = tanθ ~ N(0, slope_rms²), azimuth uniform.  Falls back to the
    mirror direction if the facet tilt would push the ray below the surface.
    """
    m = min(max(float(slope_rms), _MIN_SLOPE_RMS), _MAX_SLOPE_RMS)
    m_dir = _unit(mirror_dir)
    n = _unit(normal)

    for _ in range(4):  # few retries for grazing facet draws
        z = float(np.random.normal(0.0, m))
        z = min(max(z, -_MAX_SLOPE), _MAX_SLOPE)
        theta = math.atan(z)
        phi = float(np.random.uniform(0.0, 2.0 * math.pi))

        # Local frame: mirror direction as +Z (Gram–Schmidt tangents).
        up = np.array([0.0, 1.0, 0.0]) if abs(m_dir[0]) < 0.9 else np.array([1.0, 0.0, 0.0])
        t = np.cross(m_dir, up)
        t /= np.linalg.norm(t)
        b = np.cross(m_dir, t)

        out = (math.cos(theta) * m_dir
               + math.sin(theta) * (math.cos(phi) * t + math.sin(phi) * b))
        out = _unit(out)
        if float(np.dot(out, n)) > 0.0:
            return out
    return m_dir


def sample_surface_direction(incident_dir: np.ndarray,
                             normal: np.ndarray,
                             sigma_um: float = None,
                             tau_um: float = None,
                             wavelength_um: float = None) -> np.ndarray:
    """Sample a reflected direction from the roughness BRDF mixture.

    Mixture (reflected-power pdf):
      1. coherent specular   — probability p_spec = exp[−(4πσcosθᵢ/λ)²]
                               → exact mirror direction;
      2. Beckmann facet lobe — probability (1−p_spec)·exp[−(4πσ/λ)²]
                               → Gaussian slope spread around the mirror;
      3. diffuse background  — remaining probability → Lambertian
                               (Harvey–Shack aureole limit).

    Parameters
    ----------
    incident_dir  : incoming photon direction (points into the surface).
    normal        : surface normal (points away from the surface).
    sigma_um      : RMS roughness σ (µm).  ``None``/≤0 → legacy Lambertian.
    tau_um        : roughness correlation length τ (µm); None/≤0 → incoherent
                    power goes fully diffuse (lobe disabled).
    wavelength_um : photon wavelength λ (µm).  ``None`` → legacy Lambertian
                    (the roughness effect is wavelength-dependent).

    Returns
    -------
    np.ndarray shape (3,), unit direction with dot(dir, normal) ≥ 0.
    """
    n = _unit(normal)

    # Legacy path — identical to the pre-Phase-4b tracer behaviour.
    if (sigma_um is None or sigma_um <= 0
            or wavelength_um is None or wavelength_um <= 0):
        return sample_hemisphere_3d(n)

    d = _unit(incident_dir)
    cos_i = -float(np.dot(d, n))
    if cos_i <= 0.0:
        cos_i = min(-cos_i, 1.0) if cos_i < 0 else 1.0
    cos_i = min(max(cos_i, 0.0), 1.0)

    mirror = _mirror_direction(d, n)
    p_spec = beckmann_specular_fraction(sigma_um, wavelength_um, cos_i)

    u = float(np.random.random())
    if u < p_spec:
        return mirror

    # Incoherent branch: Beckmann lobe vs diffuse aureole.
    f_lobe = 1.0 - total_integrated_scatter(sigma_um, wavelength_um)
    if tau_um is not None and tau_um > 0 and u < p_spec + (1.0 - p_spec) * f_lobe:
        return sample_beckmann_lobe(mirror, n, rms_slope(sigma_um, tau_um))
    return sample_hemisphere_3d(n)


def brdf_provenance(sigma_um: float = None, tau_um: float = None,
                    wavelength_um: float = None,
                    cos_theta_i: float = 1.0) -> dict:
    """Diagnostic metadata describing the active BRDF model."""
    active = bool(sigma_um is not None and sigma_um > 0
                  and wavelength_um is not None and wavelength_um > 0)
    return {
        'brdf_model': ('beckmann_harvey_shack' if active else 'lambertian'),
        'sigma_um': sigma_um,
        'tau_um': tau_um,
        'wavelength_um': wavelength_um,
        'specular_fraction': (beckmann_specular_fraction(
            sigma_um, wavelength_um, cos_theta_i) if active else 0.0),
        'total_integrated_scatter': (total_integrated_scatter(
            sigma_um, wavelength_um) if active else 0.0),
        'rms_slope': (rms_slope(sigma_um, tau_um)
                      if (active and tau_um) else None),
    }


if __name__ == '__main__':
    # Quick sanity print: specularity vs roughness at 10 µm.
    print('Beckmann specular fraction at lambda=10 um, normal incidence:')
    for sig in (0.01, 0.05, 0.1, 0.5, 1.0, 5.0):
        p = beckmann_specular_fraction(sig, 10.0, 1.0)
        print(f'  sigma={sig:5.2f} um  p_spec={p:.6f}  TIS={total_integrated_scatter(sig, 10.0):.6f}')

