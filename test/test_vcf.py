import unittest
from vcf import Vcf


class TestVcf(unittest.TestCase):
    def test_convert_pval_to_neg_log10(self):
        self.assertEqual(0, Vcf.convert_pval_to_neg_log10(1))
        self.assertEqual(999, Vcf.convert_pval_to_neg_log10(0))


if __name__ == '__main__':
    unittest.main()
