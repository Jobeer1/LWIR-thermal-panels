"""
sampling.py — Monte Carlo direction sampling utilities.

Provides:
  sample_hemisphere_3d(normal)  : cosine-weighted (Lambertian) 3-D hemisphere
  rotate_to_normal(v, n)        : rotate a +Z-hemisphere sample onto arbitrary normal n
  planck_band_fraction(l1, l2, T): fractional blackbody power in wavelength band [l1, l2]
"""

import numpy as np

# ---------------------------------------------------------------------------
# 3-D Lambertian direction sampling
# ---------------------------------------------------------------------------

def sample_hemisphere_3d(normal: np.ndarray) -> np.ndarray:
    """Return a cosine-weighted (Lambertian) random unit direction in the
    hemisphere pointing away from *normal*.

    Uses Malley's method (Shirley & Morley 2003):
      - Sample a point uniformly on a unit disk in (x, y).
      - Lift to hemisphere: z = sqrt(1 - x² - y²).
      - Rotate the resulting +Z hemisphere direction to align with *normal*.

    This gives the correct 3-D Lambertian pdf:  p(θ) = cos θ / π.

    Parameters
    ----------
    normal : array-like of shape (3,)
        Outward surface normal (need not be unit length; will be normalised).

    Returns
    -------
    np.ndarray of shape (3,), unit vector.
    """
    # Uniform disk sample via rejection (fast for dense fills)
    while True:
        u1, u2 = np.random.uniform(-1.0, 1.0), np.random.uniform(-1.0, 1.0)
        if u1 * u1 + u2 * u2 <= 1.0:
            break
    x, y = u1, u2
    z = np.sqrt(max(0.0, 1.0 - x * x - y * y))
    v_local = np.array([x, y, z], dtype=float)
    return rotate_to_normal(v_local, np.asarray(normal, dtype=float))


def rotate_to_normal(v: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Rotate vector *v* (defined in the +Z hemisphere frame) so that +Z maps
    to *n*.  Uses a numerically stable Gram-Schmidt orthonormal basis (ONB).

    Parameters
    ----------
    v : np.ndarray (3,)  — vector in +Z frame, e.g. from Malley sampling.
    n : np.ndarray (3,)  — target normal (need not be unit).

    Returns
    -------
    np.ndarray (3,) unit vector.
    """
    n = n / np.linalg.norm(n)
    # Build a tangent basis (t, b) perpendicular to n
    if abs(n[0]) > 0.9:
        up = np.array([0.0, 1.0, 0.0])
    else:
        up = np.array([1.0, 0.0, 0.0])
    t = np.cross(n, up)
    t /= np.linalg.norm(t)
    b = np.cross(n, t)
    # Express v in the (t, b, n) frame
    return v[0] * t + v[1] * b + v[2] * n


# ---------------------------------------------------------------------------
# Legacy 2-D Lambertian (kept for backwards compatibility with rad leakage.py)
# ---------------------------------------------------------------------------

def sample_hemisphere_2d(normal_axis: str):
    """Correct 2-D Lambertian direction.  CDF inversion of p(θ)=½cos θ.

    Returns (dx, dz) unit vector.
    """
    s = np.random.uniform(-1.0, 1.0)   # sin θ = 2u − 1
    c = np.sqrt(max(0.0, 1.0 - s * s))
    if normal_axis == '+x':
        return c, s
    if normal_axis == '-x':
        return -c, s
    if normal_axis == '+z':
        return s, c
    if normal_axis == '-z':
        return s, -c
    raise ValueError(f'Unknown normal_axis: {normal_axis!r}')


# ---------------------------------------------------------------------------
# Planck band fractions  (for spectral model)
# ---------------------------------------------------------------------------

def planck_cumulative(lambdaT: float) -> float:
    """Cumulative blackbody fraction F(0→λT) evaluated at λT (µm·K).

    Standard Widger-Woodall series for the fractional blackbody function
    (Siegel & Howell, Appendix B). Returns value in [0, 1].
    """
    if lambdaT <= 0.0:
        return 0.0
    C2 = 14387.769  # µm·K  (second radiation constant hc/k_B)
    x = C2 / lambdaT
    total = 0.0
    for n in range(1, 512):
        en = np.exp(-n * x)
        term = en * (x * x * x / n + 3 * x * x / n**2 + 6 * x / n**3 + 6 / n**4)
        total += term
        if en < 1e-10:
            break
    return float(np.clip((15.0 / np.pi**4) * total, 0.0, 1.0))


def planck_band_fraction(lam1_um: float, lam2_um: float, T_K: float) -> float:
    """Fraction of blackbody power emitted between wavelengths lam1_um and
    lam2_um (both in micrometres) at temperature T_K [K].

    F(λ1→λ2, T) = F(0→λ2·T) − F(0→λ1·T).

    Returns a value in [0, 1].
    """
    lT1 = lam1_um * T_K
    lT2 = lam2_um * T_K
    return max(0.0, min(1.0, planck_cumulative(lT2) - planck_cumulative(lT1)))


def planck_averaged_evanescent_decay(lambda_c_um: float, T_K: float,
                                     cap_um: float = 0.0) -> float:
    """Planck-emissive-power-averaged evanescent decay length δ_ev (µm) of the
    sub-wavelength channel modes at temperature T_K.

        δ_ev(λ) = (λ_c / 2π) · (1 − (λ_c/λ)²)^(−1/2),    λ ∈ (λ_c, ∞)

    Right at cutoff the decay length diverges (the mode is nearly propagating);
    the average is therefore truncated to `cap_um` (typically the channel depth)
    so it stays finite.  Used to size the effective "rim" of cavity surface that
    couples sub-wavelength thermal emission to the aperture (LDOS confinement,
    Lin et al. PRB 2000 / Narayanaswamy & Chen PRB 2004).

    Returns 0.0 for degenerate input (no cutoff / cold emitter).
    """
    if not lambda_c_um or lambda_c_um <= 0.0 or not np.isfinite(lambda_c_um):
        return 0.0
    if T_K <= 0.0:
        return 0.0
    C2 = 14387.769  # µm·K  (second radiation constant hc/k_B)

    # Quadrature over the evanescent tail (skip λ_c exactly — singular).
    lam_hi = max(1.05 * lambda_c_um, lambda_c_um + 0.05)
    lam = np.linspace(lam_hi, min(lam_hi * 40.0, 2000.0), 20001)
    ratio = lambda_c_um / lam                      # ∈ (0, 1)
    delta = (lambda_c_um / (2.0 * np.pi)) / np.sqrt(np.maximum(1e-9, 1.0 - ratio * ratio))
    if cap_um and cap_um > 0.0:
        delta = np.minimum(delta, float(cap_um))
    # Blackbody spectral exitance M_λ ∝ λ⁻⁵ / (exp(C2/(λT)) − 1)
    with np.errstate(over='ignore', under='ignore', divide='ignore'):
        ex = np.exp(C2 / (lam * T_K))
        M = (1.0 / lam**5) / (ex - 1.0)
    M = np.maximum(M, 0.0)
    sE = float(np.sum(0.5 * (M[:-1] + M[1:]) * np.diff(lam)))
    if not (sE > 0.0) or not np.isfinite(sE):
        return 0.0
    sEd = float(np.sum(0.5 * (M[:-1] * delta[:-1] + M[1:] * delta[1:]) * np.diff(lam)))
    return float(sEd / sE)


_LAM_MIN_UM = 0.05
_LAM_MAX_UM = 2000.0
_N_GRID = 6000

# ---------------------------------------------------------------------------
# Pre-built log-wavelength grid (geometric-mean bin centres — uniform in ln λ).
# The CDF for a given temperature is cached so the expensive exp evaluation is
# done once per temperature instead of once per photon.
# ---------------------------------------------------------------------------
_PLANCK_EDGES = np.geomspace(_LAM_MIN_UM, _LAM_MAX_UM, _N_GRID + 1)
_PLANCK_MID = np.sqrt(_PLANCK_EDGES[:-1] * _PLANCK_EDGES[1:])   # geometric mid
_PLANCK_CDF_CACHE: dict = {}


def sample_planck_wavelength(T_K: float) -> float:
    """Sample a wavelength λ (µm) from the blackbody spectral exitance
    M_λ(λ, T) ∝ λ⁻⁵·(exp(C2/(λT)) − 1)⁻¹ at temperature T_K.

    Implemented by numerical inverse-CDF on a fixed log-spaced wavelength grid
    (0.05–2000 µm).  The per-temperature CDF is cached in a module-level dict so
    the cost is paid once per temperature — this keeps the Monte Carlo fast
    (≈1 ns per subsequent draw).

    Returns
    -------
    float — sampled wavelength in micrometres.
    """
    T_K = float(T_K)
    if T_K <= 0:
        raise ValueError('sample_planck_wavelength requires T_K > 0')

    # Cache key by rounded temperature (avoids spurious cache misses from
    # floating-point dust while keeping the binning negligible).
    key = round(T_K, 6)
    cdf = _PLANCK_CDF_CACHE.get(key)
    if cdf is None:
        C2 = 14387.769  # µm·K  (second radiation constant hc/k_B)
        lT = _PLANCK_MID * T_K
        # pdf per unit lnλ ∝ λ · M_λ(λ,T) = λ⁻⁴ / (exp(C2/λT) − 1)
        with np.errstate(over='ignore', under='ignore', divide='ignore'):
            ex = np.exp(C2 / lT)
            pdf = (1.0 / _PLANCK_MID**4) / (ex - 1.0)
        pdf = np.maximum(pdf, 0.0)
        s = float(pdf.sum())
        if not (s > 0) or not np.isfinite(s):
            _PLANCK_CDF_CACHE[key] = None   # remember the failure
        else:
            cdf = np.cumsum(pdf)
            cdf = cdf / cdf[-1]
            _PLANCK_CDF_CACHE[key] = cdf
        cdf = _PLANCK_CDF_CACHE.get(key)

    if cdf is None:
        # Degenerate fallback — Wien peak (deterministic).
        return 2898e-3 / T_K

    u = np.random.random()
    idx = int(np.searchsorted(cdf, u))
    if idx >= len(_PLANCK_MID):
        idx = len(_PLANCK_MID) - 1
    return float(_PLANCK_MID[min(max(idx, 0), len(_PLANCK_MID) - 1)])
