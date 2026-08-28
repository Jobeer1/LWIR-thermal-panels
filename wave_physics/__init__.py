"""
wave_physics — Integration boundary for full-wave thermal-emission solvers.

This package defines the contracts, physical conventions, analytic benchmarks,
and cached-response service that let the Monte Carlo radiative-exchange app
consume full-wave (FDTD/CMT) results without scattering physics assumptions
through the rest of the codebase.

Main modules
------------
conventions        : centralized EM sign / unit conventions.
schemas            : versioned WaveResponse / NearFieldResponse dataclasses
                     and serialization.
analytic_benchmarks: closed-form textbook limits used for validation.
cached_solver      : CachedWaveSolver adapter that loads .json/.h5 tables.
near_field_greens  : Track B2 — structured near-field Green tensor + LDOS
                     solver (offline, cache-driven).

The default solver remains the Monte Carlo ray tracer
(wave_model == 'ray'); the cached path (wave_model == 'cached') is an
opt-in service that consumes pre-calculated full-wave data through the
WaveResponse contract.  Structured near-field corrections (Track B2)
are loaded via the NearFieldResponse cache when the gap falls below
lambda_peak / (2*pi).
"""

from . import conventions
from . import schemas
from . import analytic_benchmarks
from . import cached_solver

try:
    from . import near_field_greens
except ImportError:  # pragma: no cover - optional
    near_field_greens = None  # type: ignore

__version__ = '0.1.0'

__all__ = [
    'conventions',
    'schemas',
    'analytic_benchmarks',
    'cached_solver',
    'near_field_greens',
    '__version__',
]
