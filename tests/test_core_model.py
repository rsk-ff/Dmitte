import unittest

import numpy as np

from dmitte.calc_para import _logarithmic_mean
from dmitte.guassian import gaussian_puff_model
from dmitte.transfer_equotions import solve_con


class CoreModelRegressionTests(unittest.TestCase):
    def test_gaussian_accepts_scalar_time(self):
        x = np.array([1000.0, 3000.0])
        z = np.array([20.0, 20.0])
        result = gaussian_puff_model(
            x=x,
            y=0.0,
            z=z,
            t=1800.0,
            Q_total=3.7e15,
            UU=2.0,
            HEG=20.0,
            stability=1,
            t_release=3600.0,
            puff_num=10,
        )
        self.assertEqual(result.shape, x.shape)
        self.assertTrue(np.all(np.isfinite(result)))
        self.assertTrue(np.all(result >= 0))

    def test_logarithmic_mean_handles_equal_values_elementwise(self):
        a = np.array([1.0, 2.0, 4.0])
        b = np.array([1.0, 1.0, 2.0])
        result = _logarithmic_mean(a, b)

        self.assertTrue(np.all(np.isfinite(result)))
        self.assertAlmostEqual(result[0], 1.0)
        self.assertAlmostEqual(result[1], (2.0 - 1.0) / (np.log(2.0) - np.log(1.0)))

    def test_solve_con_does_not_mutate_transfer_rates(self):
        rates = np.zeros((3, 29), dtype=np.float64)
        rates[:, 17] = 0.693
        original = rates.copy()

        result = solve_con("hto", 1.0, rates, 0.0)

        np.testing.assert_array_equal(rates, original)
        self.assertEqual(result.shape, (4, 10))
        self.assertTrue(np.all(np.isfinite(result)))


if __name__ == "__main__":
    unittest.main()
