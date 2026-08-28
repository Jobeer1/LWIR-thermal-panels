"""
wave_physics/schemas.py — Versioned WaveResponse data contract.

The WaveResponse dataclass is the *integration boundary* between the existing
Monte Carlo / radiosity app and any full-wave solver (FDTD, CMT) or local cache
of pre-computed results.  It carries the spectral/angular arrays and provenance
that let simulator.py interpolate effective absorptivity / emissivity without
knowing which solver produced them.

Schema rules
------------
* All spectral axes are wavelengths in micrometres.
* Angular axes are theta (polar, from surface normal) and phi (azimuth), both
  in radians.
* reflectance / absorptance / transmittance arrays have the same shape as the
  broadcast grid of (theta, phi) for each wavelength — a shape
  ``(n_wavelength, n_theta, n_phi)``.  Scalar axes are broadcast (size 1).
* energy_balance_error[i,j,k] = |1 - R - T - A| at that grid point.
* solver_kind is one of 'ray', 'cached', 'cmt', 'fdtd', 'rcwa'.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, Optional

import numpy as np

from . import conventions

#: Current schema version.  Bump when the contract meaning changes.
SCHEMA_VERSION = '1.0'

#: Accepted solver kinds (provenance labels).
VALID_SOLVER_KINDS = ('ray', 'cached', 'cmt', 'fdtd', 'rcwa', 'nf_greens')

#: Names of all array fields (for serialization / validation).
ARRAY_FIELDS = (
    'wavelength_um', 'theta_rad', 'phi_rad',
    'reflectance', 'absorptance', 'transmittance',
    'energy_balance_error',
)


@dataclass
class WaveResponse:
    """Versioned container for solver R/T/A response tables."""

    version: str = SCHEMA_VERSION
    solver_kind: str = 'cached'
    metadata: Dict[str, object] = field(default_factory=dict)


    # ------------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Coerce list inputs to arrays and fill the energy-balance residual
        when it was not supplied by the source solver."""
        for name in ARRAY_FIELDS:
            v = getattr(self, name)
            if not isinstance(v, np.ndarray):
                setattr(self, name, np.asarray(v, dtype=float))
        if self.solver_kind not in VALID_SOLVER_KINDS:
            raise ValueError(
                f"Invalid solver_kind {self.solver_kind!r}; "
                f"expected one of {VALID_SOLVER_KINDS}.")

        if self.energy_balance_error.size == 0:
            if self.reflectance.size and self.absorptance.size:
                self.energy_balance_error = np.abs(
                    1.0 - self.reflectance - self.absorptance - self.transmittance)

    # ------------------------------------------------------------------ #

    def shape(self) -> tuple:
        """Return the broadcast response-array shape."""
        return tuple(self.reflectance.shape)

    def max_energy_balance_error(self) -> float:
        """Worst-case energy-balance residual over the whole table."""
        if self.energy_balance_error.size == 0:
            return 0.0
        return float(np.nanmax(self.energy_balance_error))

    def energy_balance_ok(self, tol: float = 1e-6) -> bool:
        """True if the entire table conserves energy within `tol`."""
        return self.max_energy_balance_error() <= tol

    # ------------------------------------------------------------------ #
    # Serialization
    # ------------------------------------------------------------------ #

    def to_dict(self) -> Dict[str, object]:
        """Convert to a plain (JSON-friendly) dict with list arrays."""
        d: Dict[str, object] = {
            'version': self.version,
            'solver_kind': self.solver_kind,
            'metadata': self.metadata,
            'conventions': {
                'time_convention': conventions.TIME_CONVENTION,
                'complex_index_sign': conventions.COMPLEX_INDEX_SIGN,
                'poynting_flux': conventions.time_averaged_poynting_label(),
                'energy_balance': r"delta_E = |1 - R - T - A|",
            },
        }
        for name in ARRAY_FIELDS:
            d[name] = np.asarray(getattr(self, name)).tolist()
        d['near_field'] = (self.near_field.to_dict()
                           if self.near_field is not None else None)
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, object]) -> "WaveResponse":
        """Build a WaveResponse from a dict (as produced by :meth:`to_dict`)."""
        kwargs: Dict[str, object] = {
            'version': str(d.get('version', SCHEMA_VERSION)),
            'solver_kind': str(d.get('solver_kind', 'cached')),
            'metadata': dict(d.get('metadata', {})),
        }
        for name in ARRAY_FIELDS:
            vals = d.get(name, [])
            kwargs[name] = np.asarray(vals, dtype=float) if vals is not None else np.array([])
        nf = d.get('near_field')
        if nf:
            # Reconstruct the NearFieldResponse carried in the near_field slot.
            kwargs['near_field'] = NearFieldResponse.from_dict(nf)  # type: ignore[arg-type]
        return cls(**kwargs)

    # Spectral / angular axes (µm, radians).  May be 1-D grids.
    wavelength_um: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    theta_rad: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    phi_rad: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Response arrays, shape (n_w, n_theta, n_phi), range [0, 1].
    reflectance: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    absorptance: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))
    transmittance: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    # Residual |1 - R - T - A| per grid point.
    energy_balance_error: np.ndarray = field(default_factory=lambda: np.array([], dtype=float))

    #: Optional structured near-field Green-tensor + LDOS response carried
    #: alongside the hemispherical R/T/A table.  This makes WaveResponse a
    #: union container: ray/cached/fdtt/rcwa/cmt responses fill the R/T/A axes,
    #: while a NearFieldResponse (solver_kind='nf_greens') is attached when the
    #: cavity gap is sub-wavelength (g < lambda_peak / 2*pi).
    near_field: Optional["NearFieldResponse"] = None

    # ------------------------------------------------------------------ #

    def save_json(self, path: str) -> None:
        """Write the response as UTF-8 JSON (default format for this repo)."""
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(self.to_dict(), fh, indent=2)

    def save_h5(self, path: str) -> None:
        """Write the response as an HDF5 file.

        Requires ``h5py``; raises ImportError otherwise so the repo still runs
        with only the (Flask, numpy) baseline dependencies.
        """
        try:
            import h5py
        except ImportError as exc:  # pragma: no cover - env dependent
            raise ImportError(
                "HDF5 serialization requires the optional 'h5py' package. "
                "Install it (pip install h5py) or use JSON instead.") from exc

        with h5py.File(path, 'w') as f:
            f.attrs['version'] = self.version
            f.attrs['solver_kind'] = self.solver_kind
            f.attrs['time_convention'] = conventions.TIME_CONVENTION
            f.attrs['complex_index_sign'] = conventions.COMPLEX_INDEX_SIGN
            for key, value in self.metadata.items():
                f.attrs[f'meta:{key}'] = str(value)
            for name in ARRAY_FIELDS:
                f.create_dataset(name, data=np.asarray(getattr(self, name)).astype(float))

    @classmethod
    def load_h5(cls, path: str) -> "WaveResponse":
        """Load a response written by :meth:`save_h5`."""
        try:
            import h5py
        except ImportError as exc:  # pragma: no cover - env dependent
            raise ImportError(
                "HDF5 loading requires the optional 'h5py' package.") from exc

        with h5py.File(path, 'r') as f:
            version = str(f.attrs.get('version', SCHEMA_VERSION))
            solver_kind = str(f.attrs.get('solver_kind', 'cached'))
            metadata = {
                str(k)[len('meta:'):]: str(v)
                for k, v in f.attrs.items() if str(k).startswith('meta:')
            }
            arrays = {name: np.asarray(f[name][...], dtype=float) for name in ARRAY_FIELDS}
        return cls(version=version, solver_kind=solver_kind, metadata=metadata, **arrays)

    def brief(self) -> str:
        """Short human-readable provenance summary."""
        grid = self.shape()
        return (f"WaveResponse v{self.version} solver={self.solver_kind} "
                f"grid={grid} max_δE={self.max_energy_balance_error():.2e} "
                f"λ∈[{np.min(self.wavelength_um):g},{np.max(self.wavelength_um):g}]µm")


# ---------------------------------------------------------------------------
# Module-level conveniences
# ---------------------------------------------------------------------------

def save_response(response: WaveResponse, path: str, fmt: str = 'auto') -> str:
    """Persist a WaveResponse, choosing format from extension if fmt='auto'."""
    fmt = fmt or 'auto'
    if fmt == 'auto':
        ext = os.path.splitext(path)[1].lower()
        fmt = 'h5' if ext in ('.h5', '.hdf5') else 'json'
    if fmt == 'h5':
        response.save_h5(path)
    else:
        response.save_json(path)
    return fmt


def load_response(path: str, fmt: str = 'auto') -> WaveResponse:
    """Load a WaveResponse from ``path`` (JSON or HDF5 by extension)."""
    if fmt == 'auto':
        ext = os.path.splitext(path)[1].lower()
        fmt = 'h5' if ext in ('.h5', '.hdf5') else 'json'
    if fmt == 'h5':
        return WaveResponse.load_h5(path)
    with open(path, 'r', encoding='utf-8') as fh:
        return WaveResponse.from_dict(json.load(fh))


# ---------------------------------------------------------------------------
# Near-field response schema (Track B2)
# ---------------------------------------------------------------------------

#: Solver kind for cached structured near-field Green's tensor results.
NEAR_FIELD_SOLVER_KIND = 'nf_greens'


@dataclass
class NearFieldResponse:
    """Versioned container for structured near-field heat-flux tables.

    The flux table has shape ``(n_wavelength, n_gap)`` and carries
    provenance metadata so the API/UI can display the solver mode,
    validity, and cache source.  Spectral axis is wavelength in
    micrometres; gap axis is gap distance in micrometres.

    Schema rules (mirrors WaveResponse):
    * All axes are wavelengths in micrometres and gaps in micrometres.
    * solver_kind must be 'nf_greens'.
    * metadata must include 'generated_by' and 'scope_note'.
    """

    version: str = SCHEMA_VERSION
    solver_kind: str = NEAR_FIELD_SOLVER_KIND
    metadata: Dict[str, object] = field(default_factory=dict)

    # 2-D spectral / spatial grid
    wavelength_um: object = field(default_factory=lambda: np.array([]))
    gap_um: object = field(default_factory=lambda: np.array([]))
    # flux table: (n_wavelength, n_gap)
    flux_W_m2: object = field(default_factory=lambda: np.array([]))

    # ------------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Coerce list inputs to arrays and validate solver_kind."""
        for name in ('wavelength_um', 'gap_um', 'flux_W_m2'):
            v = getattr(self, name)
            if not isinstance(v, np.ndarray):
                setattr(self, name, np.asarray(v, dtype=float))

        if self.solver_kind not in ('nf_greens',) + VALID_SOLVER_KINDS:
            raise ValueError(
                f"Invalid solver_kind {self.solver_kind!r}; "
                f"expected 'nf_greens' or one of {VALID_SOLVER_KINDS}.")

    # ------------------------------------------------------------------ #

    def shape(self) -> tuple:
        """Return the (n_wavelength, n_gap) shape of the flux table."""
        return tuple(self.flux_W_m2.shape)

    def max_abs_flux(self) -> float:
        """Maximum absolute flux over the whole table (W/m^2)."""
        if self.flux_W_m2.size == 0:
            return 0.0
        return float(np.nanmax(np.abs(self.flux_W_m2)))

    # ------------------------------------------------------------------ #
    # Serialization
    # ------------------------------------------------------------------ #

    def to_dict(self) -> dict:
        """Serialize to a plain dict (JSON-compatible)."""
        return {
            'version': self.version,
            'solver_kind': self.solver_kind,
            'metadata': self.metadata,
            'wavelength_um': np.asarray(self.wavelength_um).tolist(),
            'gap_um': np.asarray(self.gap_um).tolist(),
            'flux_W_m2': np.asarray(self.flux_W_m2).tolist(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "NearFieldResponse":
        """Construct from a dict (JSON round-trip)."""
        return cls(
            version=data.get('version', SCHEMA_VERSION),
            solver_kind=data.get('solver_kind', NEAR_FIELD_SOLVER_KIND),
            metadata=data.get('metadata', {}),
            wavelength_um=np.asarray(data.get('wavelength_um', []), dtype=float),
            gap_um=np.asarray(data.get('gap_um', []), dtype=float),
            flux_W_m2=np.asarray(data.get('flux_W_m2', []), dtype=float),
        )

    def save_json(self, path: str) -> None:
        """Write the near-field response as UTF-8 JSON."""
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(self.to_dict(), fh, indent=2)

    def brief(self) -> str:
        """Short human-readable provenance summary."""
        grid = self.shape()
        return (f"NearFieldResponse v{self.version} solver={self.solver_kind} "
                f"grid={grid} max_flux={self.max_abs_flux():.2e} W/m^2 "
                f"lambda[{np.min(self.wavelength_um):g},{np.max(self.wavelength_um):g}]um "
                                f"gap[{np.min(self.gap_um):g},{np.max(self.gap_um):g}]um")


def load_near_field_response(path: str) -> NearFieldResponse:
    """Load a NearFieldResponse from a JSON file path."""
    with open(path, 'r', encoding='utf-8') as fh:
        return NearFieldResponse.from_dict(json.load(fh))

