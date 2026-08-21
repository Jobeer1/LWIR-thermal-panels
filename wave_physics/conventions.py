"""
wave_physics/conventions.py — Centralized physical conventions.

Every module in the wave-physics layer (and every external solver that feeds
data into it) MUST agree on the conventions declared here.  A mismatch in any
one of these silently corrupts phase, reflectance, absorptance, and emissivity
data.  Keeping them in one place makes the integration boundary explicit.

Declared conventions
--------------------
1. Time convention        : fields propagate as `exp(-i*w*t)`.  A plane wave in a
                            medium of complex index n~ = n + i*k therefore has
                            exp(i*(n~)*k0*z - i*w*t) = exp(-k*k0*z) * exp(i*n*k0*z-w*t)
                            i.e. a wave travelling toward +z *decays* if k > 0.

2. Complex index sign     : n~ = n + i*k   (n real part >= 0, k extinction >= 0).
                            This follows the optics literature (Born & Wolf,
                            Palik).  With the exp(-i*w*t) convention, k > 0 gives
                            decay of an absorbed wave in the +z direction.

3. Poynting flux         : time-averaged (real) power flow is
                            <S> = 0.5 * Re(E x H*).  On a planar interface the
                            transmitted *real* power balance uses Snell's law and
                            the real part of the propagation vector; reflectance
                            and transmittance are defined relative to this
                            time-averaged energy flux.

4. Energy balance        : for a passive, non-amplifying stack illuminated from
                            one side,
                                R + T + A = 1  (by direction),
                            and the reported residual is
                                delta_E = |1 - R - T - A|.
                            A solver that reports a large delta_E has a defect
                            (mesh, PML, material, or post-processing).

5. Polarization names    : 'TE' (s)  -> E perpendicular to plane of incidence.
                           'TM' (p)  -> E parallel   to plane of incidence.

6. Units                : all wavelengths in micrometres (um), all angles in
                           radians unless a function explicitly takes degrees.
                           Linear geometry is in SI metres in the existing
                           geometry.py objects; wave tables use um for the
                           spectral axis to match the rest of the app.
"""

from __future__ import annotations

import math
from typing import Tuple

# ---------------------------------------------------------------------------
# Conventions (single source of truth)
# ---------------------------------------------------------------------------

#: Harmonic time dependence.  Allowed values: '-i_omega_t' or '+i_omega_t'.
TIME_CONVENTION = '-i_omega_t'

#: Complex index-of-refraction sign convention.  Allowed: 'n_plus_ik' or 'n_minus_ik'.
COMPLEX_INDEX_SIGN = 'n_plus_ik'

#: Display / descriptive labels for the declared conventions.
TIME_CONVENTION_LABEL = r"exp(-i\omega t)" if TIME_CONVENTION == '-i_omega_t' else r"exp(+i\omega t)"
COMPLEX_INDEX_LABEL = r"\tilde{n} = n + i k"

#: Real scalar n for vacuum / free space.
VACUUM_INDEX = 1.0


def complex_index(n: float, k: float) -> complex:
    """Build the complex refractive index from real `n` and extinction `k`.

    Sign follows :data:`COMPLEX_INDEX_SIGN`.  With the '-i_omega_t' time
    convention, ``n + i*k`` gives a decaying wave for propagation toward +z.
    """
    if COMPLEX_INDEX_SIGN == 'n_plus_ik':
        return complex(n, k)
    return complex(n, -k)  # pragma: no cover - alternative convention


def time_averaged_poynting_label() -> str:
    """Human-readable description of the time-averaged Poynting flux."""
    return r"<S> = 0.5 * Re(E x H*)"  # time-averaged real power flow


# ---------------------------------------------------------------------------
# Energy-balance helpers
# ---------------------------------------------------------------------------

def energy_balance_error(reflectance: float,
                         transmittance: float,
                         absorptance: float) -> float:
    """delta_E = |1 - R - T - A| for a single-direction illumination.

    For a validated, passive stack this residual should be ~ 0 (<= numerical
    tolerance).  A large value flags an inconsistency in the solver or its
    post-processing.
    """
    return abs(1.0 - reflectance - transmittance - absorptance)


def check_energy_balance(reflectance: float,
                         transmittance: float,
                         absorptance: float,
                         tol: float = 1e-6) -> bool:
    """Return True when the stack conserves energy within `tol`."""
    return energy_balance_error(reflectance, transmittance, absorptance) <= tol


def real_snell_angle(n_incident: complex, n_transmitted: complex,
                     theta_incident_rad: float) -> Tuple[complex, bool]:
    """Propagation angle in the transmitted medium via Snell's law.

    Returns ``(sin_theta_t, is_propagating)`` where ``sin_theta_t`` is the
    (possibly complex) sine of the transmitted angle and ``is_propagating`` is
    False for total/internal reflection (real ``sin_theta_t > 1``).

    Uses the complex-index conventions above; the complex sine gives proper
    evanescent behavior (complex ``cos_theta_t``) for lossy media.
    """
    n1 = complex(n_incident)
    n2 = complex(n_transmitted)
    sin_t = n1 * math.sin(theta_incident_rad) / n2
    # Evanescent when the *real* media are non-absorbing and sin_t > 1.
    real_sin = abs(sin_t.real) if abs(n2.imag) < 1e-12 else float('nan')
    is_propagating = not (abs(n2.imag) < 1e-12 and real_sin > 1.0 + 1e-12)
    return sin_t, is_propagating
