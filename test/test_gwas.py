from gwas import Gwas
import pytest
import pysam


def test_are_alleles_iupac():
    # should not raise exec
    g = Gwas('1', 100, 'A', 'T', 1, 0.1, 0.1, 100, 0.05, 0, "snp", None, None)
    g.check_alleles_iupac()

    with pytest.raises(AssertionError):
        g = Gwas('1', 100, 'A', 'wdeT', 1, 0.1, 0.1, 100, 0.05, 0, "snp", None, None)
        g.check_alleles_iupac()


def test_normalise():
    with pysam.FastaFile("test.fasta") as fasta:
        # lhs trim
        g = Gwas("1", 10, "GT", "GA", 0, 0, 0, 0, 0, "snp1", 0, 0, 0)
        g.normalise(fasta)
        assert g.chrom == "1"
        assert g.pos == 11
        assert g.ref == "T"
        assert g.alt == "A"
        g = Gwas("1", 10, "GT", "GAA", 0, 0, 0, 0, 0, "snp1", 0, 0, 0)
        g.normalise(fasta)
        assert g.chrom == "1"
        assert g.pos == 11
        assert g.ref == "T"
        assert g.alt == "AA"
        # rhs trim
        g = Gwas("1", 10, "GT", "AT", 0, 0, 0, 0, 0, "snp1", 0, 0, 0)
        g.normalise(fasta)
        assert g.chrom == "1"
        assert g.pos == 10
        assert g.ref == "G"
        assert g.alt == "A"
        g = Gwas("1", 9, "GGT", "AT", 0, 0, 0, 0, 0, "snp1", 0, 0, 0)
        g.normalise(fasta)
        assert g.chrom == "1"
        assert g.pos == 9
        assert g.ref == "GG"
        assert g.alt == "A"
