import unittest
from vcf import Vcf
import numpy as np


class TestVcf(unittest.TestCase):
    def test_convert_pval_to_neg_log10(self):
        self.assertEqual(0, Vcf.convert_pval_to_neg_log10(1))
        self.assertEqual(999, Vcf.convert_pval_to_neg_log10(0))

    def test_is_valid_float32(self):
        # small
        self.assertFalse(Vcf.is_float32_lossy(1e-37))
        self.assertFalse(Vcf.is_float32_lossy(np.finfo(np.float32).tiny))
        self.assertTrue(Vcf.is_float32_lossy(1e-50))
        self.assertFalse(None)
        self.assertFalse(Vcf.is_float32_lossy(-1e-37))
        self.assertTrue(Vcf.is_float32_lossy(-1e-50))

        # large
        self.assertTrue(Vcf.is_float32_lossy(999999999999999999999999999999999999999))
        self.assertTrue(Vcf.is_float32_lossy(-999999999999999999999999999999999999999))
        self.assertFalse(Vcf.is_float32_lossy(-99999999999999999999999999999999999999))
        self.assertFalse(Vcf.is_float32_lossy(99999999999999999999999999999999999999))


if __name__ == '__main__':
    unittest.main()
