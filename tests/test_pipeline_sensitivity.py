import unittest
import numpy as np
from scripts.pipeline_sensitivity import paired_reduction_difference, profiled_template_amplitude


class TestPipelineSensitivity(unittest.TestCase):
    def test_difference_uncorrelated(self):
        d, s = paired_reduction_difference([5, 6], [3, 4], [2, 2], [3, 3])
        np.testing.assert_allclose(d, [2, 2])
        np.testing.assert_allclose(s, np.sqrt(13))

    def test_shared_noise_reduces_difference_uncertainty(self):
        _, independent = paired_reduction_difference([1, 2], [2, 3], [1, 1], [1, 1])
        _, shared = paired_reduction_difference([1, 2], [2, 3], [1, 1], [1, 1], 0.75)
        self.assertTrue(np.all(shared < independent))

    def test_profile_exact(self):
        t = np.array([-2., -1., 1., 2.])
        a, sigma = profiled_template_amplitude(100 + 3 * t, t, np.eye(4))
        self.assertAlmostEqual(a, 3.)
        self.assertAlmostEqual(sigma, 1 / np.sqrt(10))

    def test_offset_invariance(self):
        t = np.array([0., 1., 2., 4.])
        c = np.diag([1., 4., 9., 16.])
        a, _ = profiled_template_amplitude(50 + 2 * t, t, c)
        b, _ = profiled_template_amplitude(500 + 2 * t, t, c)
        self.assertAlmostEqual(a, b)

    def test_constant_template_rejected(self):
        with self.assertRaises(ValueError):
            profiled_template_amplitude([1, 2, 3], [5, 5, 5], np.eye(3))


if __name__ == "__main__":
    unittest.main()
