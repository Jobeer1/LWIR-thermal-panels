"""
wave_physics/cached_solver.py — Service-layer adapter for pre-computed R/T/A tables.

CachedWaveSolver loads a versioned :class:`WaveResponse` (JSON or HDF5) produced
by a full-wave solver (FDTD / CMT) and exposes the *effective* absorptivity and
emissivity that simulator.py needs, without running heavy electromagnetics in
the Flask request thread.

Responsibilities
----------------
* Load a WaveResponse from a local .json or .h5 file (or an in-memory object).
* Interpolate the hemisphere-integrated absorptance alpha(λ) and emissivity
  ε_b(λ) off the (wavelength, theta, phi) response grid.
* Optionally Planck-weight the spectrum at a given temperature to return a
  single effective value.
* Report provenance and energy-balance residuals so the API/UI can display the
  solver mode, validity, and cache source.

The Monte Carlo ray tracer remains the default fallback; this adapter is used
only when ``wave_model == 'cached'``.
"""

from __future__ import annotations

import os
from typing import Dict, Optional

import numpy as np

from . import conventions, analytic_benchmarks
from .schemas import WaveResponse, load_response

# ---------------------------------------------------------------------------
# Blackbody helpers (Planck weighting over a wavelength grid)
# ---------------------------------------------------------------------------

#: First radiation constant C1 = 2*pi*h*c^2  [W·µm^4/m²]
_C1 = 3.741771852e8
#: Second radiation constant C2 = h*c/kB      [µm·K]
_C2 = 1.438776877e4


def _blackbody_exitance(wavelength_um: np.ndarray, temperature_K: float) -> np.ndarray:
    """Planck spectral exitance M_lambda (W/m²/µm) on a wavelength grid."""
    lam = np.asarray(wavelength_um, dtype=float)
    lam = np.where(lam > 0, lam, np.nan)
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        m = _C1 / (lam ** 5 * (np.exp(_C2 / (lam * temperature_K)) - 1.0))
    return np.nan_to_num(m, nan=0.0, posinf=0.0, neginf=0.0)


def _hemisphere_average(values3d: np.ndarray,
                        theta_rad: np.ndarray,
                        phi_rad: np.ndarray) -> np.ndarray:
    """Solid-angle (Lambertian) average of a (nw, nt, np) response grid.

    Result is a length-(nw) array:  <f>(λ) = [∫∫ f cosθ sinθ dθ dφ] / π.
    """
    nw = values3d.shape[0]
    th = np.asarray(theta_rad, dtype=float)
    ph = np.asarray(phi_rad, dtype=float)

    nt = th.size
    np_ = ph.size

    # Weights: cos(theta)*sin(theta) dtheta * dphi, normalised by pi.
    dtheta = np.gradient(th) if nt > 1 else np.array([1.0])
    dphi = np.gradient(ph) if np_ > 1 else np.array([1.0])
    w_theta = np.cos(th) * np.sin(th) * dtheta          # (nt,)
    w_phi = dphi                                        # (np_,)

    # values3d: (nw, nt, np_)
    wt = w_theta.reshape(1, nt, 1)
    wf = w_phi.reshape(1, 1, np_)
    num = np.nansum(values3d * wt * wf, axis=(1, 2))    # (nw,)
    denom = np.sum(wt * wf)                              # scalar (~pi)
    return num / max(denom, 1e-12)


class CachedWaveSolver:
    """Adapter that interpolates effective optical properties from a response table."""

    def __init__(self, source: str = '',
                 response: Optional[WaveResponse] = None,
                 fmt: str = 'auto'):
        """Construct from a file path (:attr:`source`) or an in-memory response.

        Parameters
        ----------
        source   : path to a .json or .h5 WaveResponse file.  If empty, a
                   :attr:`response` object must be supplied.
        response : an already-loaded :class:`WaveResponse`.
        fmt      : 'auto' | 'json' | 'h5' (ignored when ``response`` is given).
        """
        if response is not None:
            self.response = response
            self.source = '<in-memory>'
        else:
            if not source:
                raise ValueError(
                    "CachedWaveSolver requires either 'source' or 'response'.")
            self.response = load_response(source, fmt=fmt)
            self.source = source

        self.response_shape = self.response.shape()
        self._build_interpolators()

    # ------------------------------------------------------------------ #

    def _build_interpolators(self) -> None:
        """Precompute hemisphere-averaged absorptance vs wavelength."""
        res = self.response
        A = np.asarray(res.absorptance, dtype=float)
        if A.ndim != 3:
            nw = res.wavelength_um.size
            n_resp = A.size
            if n_resp == 0:
                raise ValueError("WaveResponse has empty absorptance array.")
            if n_resp == nw:
                # Per-wavelength scalar: broadcast to (nw, 1, 1).
                A = A.reshape(nw, 1, 1)
            else:
                # Flat response (all entries identical): treat as a single
                # angular point per wavelength of a size-nw axis.
                if n_resp == 1:
                    A = np.full((nw, 1, 1), A.flat[0])
                else:
                    A = A.reshape(nw, 1, 1)
        self._alpha_spectrum = _hemisphere_average(A, res.theta_rad, res.phi_rad)
        self._eps_spectrum = self._alpha_spectrum.copy()
        self._wavelengths = np.asarray(res.wavelength_um, dtype=float)


    # ------------------------------------------------------------------ #

    def _interp1(self, xq: np.ndarray) -> np.ndarray:
        """Linear-interpolate the precomputed spectrum at query wavelengths (µm)."""
        x = self._wavelengths
        y = self._alpha_spectrum
        if x.size == 0:
            raise ValueError("No wavelength axis in cached response.")
        return np.interp(xq, x, y, left=y[0], right=y[-1])

    # ------------------------------------------------------------------ #

    def absorptance_spectrum(self, wavelengths_um: np.ndarray) -> np.ndarray:
        """Hemisphere-averaged absorptance alpha(λ) at the given wavelengths."""
        return np.clip(self._interp1(np.asarray(wavelengths_um, dtype=float)), 0.0, 1.0)

    def emissivity_spectrum(self, wavelengths_um: np.ndarray) -> np.ndarray:
        """Hemisphere-averaged emissivity ε(λ) (Kirchhoff basis, see notes)."""
        return np.clip(self._interp1(np.asarray(wavelengths_um, dtype=float)), 0.0, 1.0)

    # ------------------------------------------------------------------ #

    def alpha_eff(self, wavelength_um: Optional[float] = None,
                  temperature_K: Optional[float] = None) -> float:
        """Effective absorptivity of the structure.

        * If ``wavelength_um`` given      -> hemisphere-averaged α at that λ.
        * If ``temperature_K`` given      -> Planck-weighted effective α(T).
        * Otherwise                       -> plain mean over the cached spectrum.
        """
        if wavelength_um is not None:
            return float(self.absorptance_spectrum(np.array([wavelength_um]))[0])
        if temperature_K is not None and temperature_K > 0:
            return float(self._planck_weighted(self._alpha_spectrum, temperature_K))
        return float(np.nanmean(self._alpha_spectrum))

    def epsilon_b(self, wavelength_um: Optional[float] = None,
                  temperature_K: Optional[float] = None) -> float:
        """Effective emissivity of the structure (see class notes)."""
        if wavelength_um is not None:
            return float(self.emissivity_spectrum(np.array([wavelength_um]))[0])
        if temperature_K is not None and temperature_K > 0:
            return float(self._planck_weighted(self._eps_spectrum, temperature_K))
        return float(np.nanmean(self._eps_spectrum))

    # ------------------------------------------------------------------ #

    def _planck_weighted(self, spectrum: np.ndarray, temperature_K: float) -> float:
        """<f>_T = ∫ f(λ) M_λ(T) dλ / ∫ M_λ(T) dλ over the cached wavelength axis."""
        lam = self._wavelengths
        m = _blackbody_exitance(lam, temperature_K)
        denom = np.sum(m)
        if denom <= 0:
            return float(np.nanmean(np.asarray(spectrum)))
        spec = np.asarray(spectrum, dtype=float)
        return float(np.sum(spec * m) / denom)

    # ------------------------------------------------------------------ #

    def info(self) -> Dict[str, object]:
        """Provenance / validity diagnostics for UI display."""
        res = self.response
        return {
            'wave_model': 'cached',
            'solver_status': 'cached',
            'solver_kind': res.solver_kind,
            'source': self.source,
            'schema_version': res.version,
            'response_shape': self.response_shape,
            'wavelength_range_um': [float(np.min(self._wavelengths)),
                                    float(np.max(self._wavelengths))],
            'theta_range_rad': ([float(np.min(res.theta_rad)),
                                 float(np.max(res.theta_rad))]
                                if res.theta_rad.size else None),
            'energy_balance_error': float(res.max_energy_balance_error()),
            'energy_conservation_ok': bool(res.energy_balance_ok(tol=1e-6)),
            'metadata': dict(res.metadata),
        }



# ---------------------------------------------------------------------------
# Default (analytic) cache builder
# ---------------------------------------------------------------------------

def build_analytic_cache(
    wall_thickness_um: float = 1.0,
    wavelengths_um: Optional[np.ndarray] = None,
    theta_max_rad: float = np.pi / 2.0,
    n_theta: int = 9,
    n_phi: int = 4,
    wall_index: complex = 1.7 + 0.5j,
    substrate_index: complex = 3.5 + 0.05j,
    incident_index: complex = 1.0 + 0.0j,
    temperature_K: float = 300.0,
    geometry_diameter_um: float = 20.0,
    geometry_height_um: float = 450.0,
) -> WaveResponse:
    """Build a default full-wave-style cache from the analytic multilayer model.

    The cache models a single planar wall film on a lossy substrate (the
    honeycomb wall).  The spectral/angular absorptance is computed with the
    exact analytic multilayer Fresnel benchmark, giving a physically consistent,
    energy-conserving table that the cached simulator path can consume.

    This is a *demonstration* dataset (planar wall emissivity), NOT a full
    cavity modal FDTD solution.  Its provenance is recorded in the metadata so
    downstream consumers never mistake it for a true cavity-mode result.
    """
    if wavelengths_um is None:
        wavelengths_um = np.linspace(0.5, 30.0, 60)
    wavelengths_um = np.asarray(wavelengths_um, dtype=float)

    theta = np.linspace(0.0, theta_max_rad, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi)

    # indices array: incident | film | substrate
    indices = [incident_index, wall_index, substrate_index]

    nw, nt, np_ = len(wavelengths_um), len(theta), len(phi)
    A = np.zeros((nw, nt, np_))
    T = np.zeros((nw, nt, np_))
    R = np.zeros((nw, nt, np_))
    for i, lam in enumerate(wavelengths_um):
        for j, th in enumerate(theta):
            r_te, t_te, a_te = analytic_benchmarks.multilayer_stack_rt(
                indices, [wall_thickness_um], lam, th, 'TE')
            r_tm, t_tm, a_tm = analytic_benchmarks.multilayer_stack_rt(
                indices, [wall_thickness_um], lam, th, 'TM')
            # Unpolarized average.
            rr = 0.5 * (r_te + r_tm)
            tt = 0.5 * (t_te + t_tm)
            aa = 0.5 * (a_te + a_tm)
            R[i, j, :] = rr
            T[i, j, :] = tt
            A[i, j, :] = aa

    res = WaveResponse(
        solver_kind='cached',
        wavelength_um=wavelengths_um,
        theta_rad=theta,
        phi_rad=phi,
        reflectance=R,
        absorptance=A,
        transmittance=T,
        metadata={
            'generated_by': 'wave_physics.cached_solver.build_analytic_cache',
            'model': 'planar multilayer Fresnel (wall film on lossy substrate)',
            'note': 'Demonstration dataset; NOT a cavity-modal FDTD solution.',
            'wall_thickness_um': str(wall_thickness_um),
            'wall_index': str(wall_index),
            'substrate_index': str(substrate_index),
            'temperature_K': str(temperature_K),
            'geometry_diameter_um': str(geometry_diameter_um),
            'geometry_height_um': str(geometry_height_um),
            'conventions': conventions.TIME_CONVENTION + ' / ' + conventions.COMPLEX_INDEX_SIGN,
        },
    )
    return res


def default_cache_path() -> str:
    """Absolute path of the bundled default JSON cache."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, 'cache', 'default_wave_response.json')


def ensure_default_cache(overwrite: bool = False) -> str:
    """Make sure the default JSON cache exists on disk; return its path."""
    path = default_cache_path()
    if not os.path.exists(path) or overwrite:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        res = build_analytic_cache()
        res.save_json(path)
    return path


if __name__ == '__main__':
    # Regenerate the bundled default cache, e.g. to tweak default parameters.
    import sys
    path = ensure_default_cache(overwrite='--overwrite' in sys.argv)
    print(f'Default wave-response cache: {path}')


