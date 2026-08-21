#!/usr/bin/env python3
"""
test_wave_benchmarks.py — Continuous validation of the wave-physics layer.

Verifies that any new wave module matches textbook analytic limits:

* Single-interface and multi-layer Fresnel (TE/TM, arbitrary angle)
  against closed-form reflectivity and energy balance.
* Rectangular (2*min(W,D)) and circular (1.706*D) waveguide cutoffs.
* Analytic evanescent power decay  T(h) = exp(-2*kappa*h).
* The versioned WaveResponse schema (JSON + optional HDF5 round-trip).
* The CachedWaveSolver adapter consumes a response table cleanly.

Run:
    python test_wave_benchmarks.py          (unittest)
    python -m pytest test_wave_benchmarks.py
"""

import math
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from wave_physics import conventions, analytic_benchmarks as ab
from wave_physics.schemas import WaveResponse, save_response, load_response
from wave_physics.cached_solver import CachedWaveSolver, build_analytic_cache

try:
    import h5py  # noqa: F401
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False


def _closed_form_single_film(n0, n1, ns, wavelength, thickness):
    """Independent closed-form reflectivity for one film on a substrate
    (normal incidence, TE).  r = (r01 + r1s*e^{2i*beta}) / (1 + r01*r1s*e^{2i*beta})."""
    import cmath
    k0 = 2.0 * math.pi / wavelength
    r01 = (n0 - n1) / (n0 + n1)
    r1s = (n1 - ns) / (n1 + ns)
    beta = k0 * n1 * thickness
    r = (r01 + r1s * cmath.exp(2j * beta)) / (1.0 + r01 * r1s * cmath.exp(2j * beta))
    return abs(r) ** 2


class TestConventions(unittest.TestCase):
    def test_time_and_index_conventions(self):
        self.assertEqual(conventions.TIME_CONVENTION, '-i_omega_t')
        self.assertEqual(conventions.COMPLEX_INDEX_SIGN, 'n_plus_ik')
        self.assertEqual(conventions.complex_index(1.7, 0.5), complex(1.7, 0.5))

    def test_energy_balance_definition(self):
        # delta_E = |1 - R - T - A|
        self.assertAlmostEqual(conventions.energy_balance_error(0.3, 0.2, 0.5), 0.0)
        self.assertAlmostEqual(conventions.energy_balance_error(0.3, 0.2, 0.4), 0.1)
        self.assertFalse(conventions.check_energy_balance(0.3, 0.2, 0.4))
        self.assertTrue(conventions.check_energy_balance(0.3, 0.2, 0.5, tol=1e-12))


class TestFresnel(unittest.TestCase):
    def test_normal_incidence_air_glass(self):
        # Textbook: R = ((n1-n2)/(n1+n2))^2 with n1=1, n2=1.5
        R, T, A = ab.single_interface_rt(1.0, 1.5, 0.0, 'TE')
        self.assertAlmostEqual(R, ((1.0 - 1.5) / (1.0 + 1.5)) ** 2, places=10)
        self.assertAlmostEqual(T, 1.0 - R, places=10)
        self.assertAlmostEqual(A, 0.0, places=10)

    def test_tm_brewster_zero(self):
        # Lossless single interface: TM reflectance -> 0 at Brewster angle.
        th_b = math.atan(1.5 / 1.0)
        R, _, _ = ab.single_interface_rt(1.0, 1.5, th_b, 'TM')
        self.assertLess(R, 1e-12)


class TestMultilayer(unittest.TestCase):
    def setUp(self):
        self.inds = [1.0 + 0j, 1.7 + 0.5j, 3.5 + 0.05j]

    def test_normal_incidence_matches_closed_form(self):
        R, T, A = ab.multilayer_stack_rt(self.inds, [1.0], 5.0, 0.0, 'TE')
        Rref = _closed_form_single_film(1.0, 1.7 + 0.5j, 3.5 + 0.05j, 5.0, 1.0)
        self.assertAlmostEqual(R, Rref, places=9)
        self.assertAlmostEqual(R + T + A, 1.0, places=9)

    def test_zero_thickness_reduces_to_single_interface(self):
        # A zero-thickness film must reproduce the bare substrate interface.
        R, T, A = ab.multilayer_stack_rt(self.inds, [0.0], 5.0, 0.0, 'TE')
        R_if, T_if, _ = ab.single_interface_rt(1.0, self.inds[-1], 0.0, 'TE')
        self.assertAlmostEqual(R, R_if, places=10)
        self.assertAlmostEqual(R + T + A, 1.0, places=10)


    def test_lossless_slab_rt_equals_one(self):
        # air | glass d | air-matching (lossless), normal incidence.
        R, T, A = ab.multilayer_stack_rt([1.0, 1.5, 1.0], [1.0], 5.0, 0.0, 'TE')
        self.assertAlmostEqual(R + T, 1.0, places=9)
        self.assertAlmostEqual(A, 0.0, places=9)

    def test_lossy_slab_rt_plus_a_equals_one(self):
        R, T, A = ab.multilayer_stack_rt(self.inds, [1.0], 5.0, 0.0, 'TM')
        self.assertAlmostEqual(R + T + A, 1.0, places=9)
        self.assertGreater(A, 0.0)  # the lossy film absorbs

    def test_tm_and_te_differ_at_angle(self):
        R_te, _, _ = ab.multilayer_stack_rt([1.0, 1.5, 1.0], [1.0], 5.0, 1.0, 'TE')
        R_tm, _, _ = ab.multilayer_stack_rt([1.0, 1.5, 1.0], [1.0], 5.0, 1.0, 'TM')
        self.assertNotAlmostEqual(R_te, R_tm, places=3)


class TestWaveguideCutoff(unittest.TestCase):
    def test_rectangular_narrowest_dimension(self):
        self.assertAlmostEqual(ab.rectangular_cutoff_wavelength_um(100.0, 10.0), 20.0)
        self.assertAlmostEqual(ab.rectangular_cutoff_wavelength_um(10.0, 100.0), 20.0)

    def test_circular_te11(self):
        D = 20.0
        lam_c = ab.circular_cutoff_wavelength_um(D)
        self.assertAlmostEqual(lam_c / D, 1.7062, places=3)
        self.assertAlmostEqual(lam_c / D, math.pi / 1.84118, places=8)


class TestEvanescent(unittest.TestCase):
    def test_power_decay_matches_exp(self):
        lc, lam, h = 1.0, 2.0, 5.0
        kappa = ab.evanescent_decay_constant_um(lc, lam)
        expected = math.exp(-2.0 * kappa * h)
        self.assertAlmostEqual(
            ab.evanescent_power_transmission(lc, lam, h), expected, places=12)

    def test_formula_and_limits(self):
        lc, lam = 1.0, 2.0
        kappa = ab.evanescent_decay_constant_um(lc, lam)
        self.assertAlmostEqual(kappa, 2.0 * math.pi * math.sqrt(1 - 0.25), places=9)
        self.assertEqual(ab.evanescent_decay_constant_um(lc, 1.0), 0.0)   # propagating
        self.assertEqual(ab.evanescent_power_transmission(lc, 2.0, 0.0), 1.0)
        self.assertLess(ab.evanescent_power_transmission(lc, lam, 100.0), 1e-6)


class TestSchema(unittest.TestCase):
    def make_response(self):
        wl = np.linspace(1.0, 20.0, 5)
        th = np.array([0.0, 0.5, 1.0])
        a = np.full((5, 3, 1), 0.5)
        t = np.full((5, 3, 1), 0.3)
        r = 1.0 - a - t
        return WaveResponse(
            solver_kind='fdtd', wavelength_um=wl, theta_rad=th, phi_rad=[0.0],
            reflectance=r, absorptance=a, transmittance=t)

    def test_energy_balance_auto_filled(self):
        resp = self.make_response()
        self.assertAlmostEqual(resp.energy_balance_error.max(), 0.0, places=12)
        self.assertTrue(resp.energy_balance_ok(tol=1e-6))

    def test_json_round_trip(self):
        import json
        resp = self.make_response()
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'resp.json')
            resp.save_json(p)
            with open(p, encoding='utf-8') as fh:
                loaded = WaveResponse.from_dict(json.load(fh))
        np.testing.assert_allclose(loaded.absorptance, resp.absorptance)
        np.testing.assert_allclose(loaded.wavelength_um, resp.wavelength_um)
        self.assertEqual(loaded.solver_kind, 'fdtd')

    @unittest.skipUnless(HAS_H5PY, "h5py not installed")
    def test_h5_round_trip(self):
        resp = self.make_response()
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'resp.h5')
            save_response(resp, p)   # fmt auto by extension
            loaded = load_response(p)
        np.testing.assert_allclose(loaded.reflectance, resp.reflectance)
        self.assertEqual(loaded.solver_kind, 'fdtd')


class TestCachedSolver(unittest.TestCase):
    def test_values_in_unit_range(self):
        solver = CachedWaveSolver(response=build_analytic_cache())
        a = solver.alpha_eff(wavelength_um=5.0)
        self.assertGreaterEqual(a, 0.0)
        self.assertLessEqual(a, 1.0)
        self.assertGreater(solver.epsilon_b(temperature_K=300.0), 0.0)

    def test_info_reports_conservation(self):
        solver = CachedWaveSolver(response=build_analytic_cache())
        info = solver.info()
        self.assertEqual(info['solver_status'], 'cached')
        self.assertTrue(info['energy_conservation_ok'])

    def test_alpha_spectrum_shape(self):
        solver = CachedWaveSolver(response=build_analytic_cache())
        spec = solver.absorptance_spectrum(np.linspace(2.0, 15.0, 20))
        self.assertEqual(spec.shape, (20,))
        self.assertTrue(np.all((spec >= 0.0) & (spec <= 1.0)))


class TestAnalyticDefaultCacheFile(unittest.TestCase):
    def test_default_json_cache_loads(self):
        from wave_physics.cached_solver import ensure_default_cache
        path = ensure_default_cache()
        self.assertTrue(os.path.exists(path))
        solver = CachedWaveSolver(path)
        self.assertTrue(solver.info()['energy_conservation_ok'])


if __name__ == '__main__':
    unittest.main(verbosity=2)

