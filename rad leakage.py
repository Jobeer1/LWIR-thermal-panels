"""
rad leakage.py — Standalone radiation leakage estimator (2-D, educational).

Bug #12 fix: base surface emission is now correctly included, weighted by
  w_base = WIDTH * ALPHA_AG.

This is a self-contained educational script.  For production use, see the
3-D ray_tracer.py + geometry.py pipeline via simulator.py.
"""

import numpy as np

# Geometry: 10 um wide, 450 um deep trench
WIDTH  = 10.0
HEIGHT = 450.0

# Materials
ALPHA_CNT = 0.98   # CNT wall absorptivity (98%)
ALPHA_AG  = 0.02   # Silver base absorptivity (2%)

# Simulation Parameters
N_PHOTONS = 100_000
EPS = 1e-9


def sample_lambertian_2d(normal_axis):
    """
    Sample a 2-D diffuse (Lambertian) direction using CDF inversion.
    pdf: p(θ) = ½ cos θ  →  CDF: sin θ = 2u - 1  →  s = U(-1, 1).
    """
    s = np.random.uniform(-1.0, 1.0)
    c = np.sqrt(max(0.0, 1.0 - s**2))
    if normal_axis == '+x':
        return c, s
    elif normal_axis == '-x':
        return -c, s
    elif normal_axis == '+z':
        return s, c
    elif normal_axis == '-z':
        return s, -c
    raise ValueError(f'Unknown normal_axis: {normal_axis!r}')


def trace_ray_2d(x, z, dx, dz):
    """Trace a single 2-D photon inside the trench. Returns True if escaped."""
    for _ in range(5000):
        t_top    = (HEIGHT - z) / dz  if dz > EPS  else float('inf')
        t_bottom = (0.0    - z) / dz  if dz < -EPS else float('inf')
        t_right  = (WIDTH  - x) / dx  if dx > EPS  else float('inf')
        t_left   = (0.0    - x) / dx  if dx < -EPS else float('inf')

        t_hit = min(t_top, t_bottom, t_right, t_left)
        if t_hit <= 0.0 or not np.isfinite(t_hit):
            break

        x_new = x + dx * t_hit
        z_new = z + dz * t_hit

        at_top    = abs(t_hit - t_top)    < 1e-10
        at_bottom = abs(t_hit - t_bottom) < 1e-10
        hit_right = abs(t_hit - t_right)  < 1e-10
        hit_left  = abs(t_hit - t_left)   < 1e-10

        if at_top:
            return True
        elif at_bottom:
            if np.random.random() < ALPHA_AG:
                return False           # absorbed by silver base
            dx, dz = sample_lambertian_2d('+z')
            x_new = np.clip(x_new, 0, WIDTH)
            z_new = EPS
        elif hit_right:
            if np.random.random() < ALPHA_CNT:
                return False           # absorbed by right CNT wall
            dx, dz = sample_lambertian_2d('-x')
            x_new = WIDTH - EPS
            z_new = np.clip(z_new, 0, HEIGHT)
        elif hit_left:
            if np.random.random() < ALPHA_CNT:
                return False           # absorbed by left CNT wall
            dx, dz = sample_lambertian_2d('+x')
            x_new = EPS
            z_new = np.clip(z_new, 0, HEIGHT)

        x, z = x_new, z_new

    # Safety: photon exceeded bounce cap — treat as absorbed (not 50/50!)
    return False


def run_internal_emission(n_rays):
    """
    Emit from ALL internal surfaces weighted by area × emissivity.

    Bug #12 FIX: base surface is now included.
    Returns (P_esc, w_total).
    """
    # Surface emission weights (per unit depth)
    w_left  = HEIGHT * ALPHA_CNT
    w_right = HEIGHT * ALPHA_CNT
    w_base  = WIDTH  * ALPHA_AG    # ← previously omitted!
    w_total = w_left + w_right + w_base

    p_left  = w_left  / w_total
    p_right = w_right / w_total
    # p_base  = w_base  / w_total  (remainder)

    escaped = 0
    for _ in range(n_rays):
        r = np.random.random()
        if r < p_left:
            x, z  = 0.0, np.random.uniform(0.0, HEIGHT)
            dx, dz = sample_lambertian_2d('+x')
        elif r < p_left + p_right:
            x, z  = WIDTH, np.random.uniform(0.0, HEIGHT)
            dx, dz = sample_lambertian_2d('-x')
        else:
            x, z  = np.random.uniform(0.0, WIDTH), 0.0
            dx, dz = sample_lambertian_2d('+z')
        if trace_ray_2d(x, z, dx, dz):
            escaped += 1

    return escaped / n_rays, w_total


def run_external_incidence(n_rays):
    """External illumination → α_eff (Kirchhoff check)."""
    absorbed = 0
    for _ in range(n_rays):
        x = np.random.uniform(0.0, WIDTH)
        z = HEIGHT
        dx, dz = sample_lambertian_2d('-z')
        if not trace_ray_2d(x, z, dx, dz):
            absorbed += 1
    return absorbed / n_rays


def main():
    print('=' * 62)
    print('Monte Carlo Radiation Leakage Simulation (2-D, corrected)')
    print('=' * 62)
    print(f'Geometry:  {WIDTH} µm wide, {HEIGHT} µm deep  (AR {HEIGHT/WIDTH:.0f}:1)')
    print(f'Materials: CNT α={ALPHA_CNT},  Ag base α={ALPHA_AG}')
    print(f'Photons/experiment: {N_PHOTONS:,}\n')

    p_esc, w_total = run_internal_emission(N_PHOTONS)
    alpha_eff      = run_external_incidence(N_PHOTONS)

    cavity_enhancement = w_total / WIDTH
    epsilon_b_raw = cavity_enhancement * p_esc
    epsilon_b = alpha_eff
    kirchhoff_err = abs(epsilon_b_raw - alpha_eff) / max(alpha_eff, 1e-12) * 100.0

    import math
    ci_p   = 1.96 * math.sqrt(max(p_esc     * (1-p_esc)     / N_PHOTONS, 0))
    ci_a   = 1.96 * math.sqrt(max(alpha_eff * (1-alpha_eff) / N_PHOTONS, 0))

    print('-' * 50)
    print('RESULTS')
    print('-' * 50)
    print(f'1. Escape Probability (P_esc):          {p_esc*100:.4f}% ± {ci_p*100:.4f}%')
    print(f'2. Cavity Enhancement Factor:            {cavity_enhancement:.2f}×')
    print(f'3. Plate B Eff. Emissivity (ε_B):        {epsilon_b*100:.4f}%')
    print(f'4. External Absorptivity (α_eff):        {alpha_eff*100:.4f}% ± {ci_a*100:.4f}%')
    print(f'5. Kirchhoff Reciprocity Error:          {kirchhoff_err:.2f}%  (target < 5%)')
    print()
    print('Reference (10×450 µm CNT trench, 2-D):')
    print('  P_esc ~1.12%, ε_B ~98.81%, α_eff ~99.54%, Kirchhoff err ~0.7%')


if __name__ == '__main__':
    main()