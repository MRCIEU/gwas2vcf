from gwas import Gwas
import pytest


def test_are_alleles_iupac():
    # should not raise exec
    g = Gwas('1', 100, 'A', 'T', 1, 0.1, 0.1, 100, 0.05, 0, "snp", None, None)
    g.check_alleles_iupac()

    with pytest.raises(AssertionError):
        g = Gwas('1', 100, 'A', 'wdeT', 1, 0.1, 0.1, 100, 0.05, 0, "snp", None, None)
        g.check_alleles_iupac()
