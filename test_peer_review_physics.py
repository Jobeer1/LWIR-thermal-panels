import math
import unittest

from geometry import RectPit3D
from material_optics import effective_emissivity_thin_film
from ray_tracer import evanescent_decay_length, evanescent_power_transmission


class PeerReviewPhysicsTests(unittest.TestCase):
    def test_tmm_thin_film_is_bounded_and_thickness_sensitive(self):
        thin = effective_emissivity_thin_film(
            bulk_emissivity=0.8,
            thickness_um=0.1,
            wavelength_um=10.0,
            material='alumina',
        )
        thick = effective_emissivity_thin_film(
            bulk_emissivity=0.8,
            thickness_um=10.0,
            wavelength_um=10.0,
            material='alumina',
        )
        self.assertGreaterEqual(thin, 0.0)
        self.assertLessEqual(thin, 0.8)
        self.assertLessEqual(thick, 0.8)
        self.assertNotEqual(thin, 0.8)

    def test_rectangular_cutoff_uses_narrowest_dimension(self):
        geometry = RectPit3D(width_um=100.0, depth_um=10.0, height_um=450.0)
        self.assertEqual(geometry.channel_cutoff_wavelength_um(), 20.0)

    def test_evanescent_transmission_uses_power_exponent(self):
        cutoff = 1.0
        wavelength = 2.0
        decay_length = evanescent_decay_length(cutoff, wavelength)
        transmission = evanescent_power_transmission(
            cutoff, wavelength, decay_length
        )
        self.assertAlmostEqual(transmission, math.exp(-2.0), places=12)
        self.assertEqual(evanescent_power_transmission(cutoff, wavelength, 0.0), 1.0)


if __name__ == '__main__':
    unittest.main()
