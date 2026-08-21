"""
wave_physics — Integration boundary for full-wave thermal-emission solvers.

This package defines the contracts, physical conventions, analytic benchmarks,
and cached-response service that let the Monte Carlo radiative-exchange app
consume full-wave (FDTD/CMT) results without scattering physics assumptions
through the rest of the codebase.

Main modules
------------
conventions          : centralized electromagnetic sign / unit conventions.
schemas              : versioned WaveResponse dataclass and serialization.
analytic_benchmarks  : closed-form textbook limits used for validation.
cached_solver        : CachedWaveSolver adapter that loads .json / .h5 tables.
make_cache           : builds a default analytic-benchmark response cache.

The default solver remains the Monte Carlo ray tracer (wave_model == 'ray');
the cached path (wave_model == 'cached') is an opt-in service that consumes
pre-calculated full-wave data through the WaveResponse contract.
"""

from . import conventions
from . import schemas
from . import analytic_benchmarks
from . import cached_solver

__version__ = '0.1.0'

__all__ = [
    'conventions',
    'schemas',
    'analytic_benchmarks',
    'cached_solver',
    '__version__',
]
