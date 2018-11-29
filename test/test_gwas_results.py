import unittest
from gwas_results import GwasResult


class TestGwasResults(unittest.TestCase):
    def test_are_alleles_iupac(self):
        g = GwasResult('1', 100, 'A', 'T', 1, 0.1, 0.1, 100, 0.05, 0,  "snp")
        self.assertTrue(g.are_alleles_iupac())
        g = GwasResult('1', 100, 'A', 'wdeT', 1, 0.1, 0.1, 100, 0.05, 0, "snp")
        self.assertFalse(g.are_alleles_iupac())
