import unittest
from vcf import Vcf
import numpy as np


class TestVcf(unittest.TestCase):
    def test_convert_pval_to_neg_log10(self):
        assert 0 == Vcf.convert_pval_to_neg_log10(1)
        assert 999 == Vcf.convert_pval_to_neg_log10(0)

    def test_is_valid_float32(self):
        # small
        assert Vcf.is_float32_lossy(1e-37) == False
        assert Vcf.is_float32_lossy(np.finfo(np.float32).tiny) == False
        assert Vcf.is_float32_lossy(1e-50) == True
        assert Vcf.is_float32_lossy(-1e-37) == False
        assert Vcf.is_float32_lossy(-1e-50) == True

        # large
        assert Vcf.is_float32_lossy(999999999999999999999999999999999999999) == True
        assert Vcf.is_float32_lossy(-999999999999999999999999999999999999999) == True
        assert Vcf.is_float32_lossy(-99999999999999999999999999999999999999) == False
        assert Vcf.is_float32_lossy(99999999999999999999999999999999999999) == False


if __name__ == '__main__':
    unittest.main()
