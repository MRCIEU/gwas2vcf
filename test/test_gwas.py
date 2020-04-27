from gwas import Gwas
import pytest
from vgraph import norm
import pysam


def test_are_alleles_iupac():
    g = Gwas('1', 1, 'A', 'T', None, None, None, None, None, None, None, None, None)
    g.check_alleles_iupac()

    with pytest.raises(AssertionError):
        g = Gwas('1', 1, 'A', 'wdeT', None, None, None, None, None, None, None, None, None)
        g.check_alleles_iupac()


def test_normalise():
    with pysam.FastaFile("test.fasta") as fasta:
        # SNV
        g = Gwas('test', 1, 'A', 'T', None, None, None, None, None, None, None, None, None)
        g.check_reference_allele(fasta)
        g.normalise(fasta)
        assert g.chrom == "test"
        assert g.pos == 1
        assert g.ref == "A"
        assert g.alt == "T"

        # left pad SNV
        g = Gwas('test', 10, 'ACACA', 'ACACT', None, None, None, None, None, None, None, None, None)
        g.check_reference_allele(fasta)
        g.normalise(fasta)
        assert g.chrom == "test"
        assert g.pos == 14
        assert g.ref == "A"
        assert g.alt == "T"

        # right pad SNV
        g = Gwas('test', 10, 'ACACA', 'TCACA', None, None, None, None, None, None, None, None, None)
        g.check_reference_allele(fasta)
        g.normalise(fasta)
        assert g.chrom == "test"
        assert g.pos == 10
        assert g.ref == "A"
        assert g.alt == "T"

        # left pad ins
        g = Gwas('test', 10, 'ACA', 'ACAGT', None, None, None, None, None, None, None, None, None)
        g.check_reference_allele(fasta)
        g.normalise(fasta)
        assert g.chrom == "test"
        assert g.pos == 12
        assert g.ref == "A"
        assert g.alt == "AGT"

        # left-align pad del
        g = Gwas('test', 10, 'ACACA', 'ACA', None, None, None, None, None, None, None, None, None)
        g.check_reference_allele(fasta)
        g.normalise(fasta)
        assert g.chrom == "test"
        assert g.pos == 9
        assert g.ref == "TAC"
        assert g.alt == "T"


def test_bioinformed_vgraph_normalise():
    """ Cython logic from bioinformed/vgraph """

    # get a reference sequence
    with pysam.FastaFile("test.fasta") as fasta:
        seq = fasta.fetch("test").upper()

        # SNV

        # use 0-based
        start = 0
        # use 1-based
        stop = 1
        ref = "A"
        alt = "T"

        trim_start, trim_stop, trim_alleles = norm.normalize_alleles(
            seq,
            start,
            stop,
            (ref, alt)
        )
        assert trim_alleles == ("A", "T")
        assert trim_start == 0
        assert trim_stop == 1

        # left pad SNV
        start = 9
        stop = 14
        ref = "ACACA"
        alt = "ACACT"

        trim_start, trim_stop, trim_alleles = norm.normalize_alleles(
            seq,
            start,
            stop,
            (ref, alt)
        )
        assert trim_alleles == ("A", "T")
        assert trim_start == 13
        assert trim_stop == 14

        # right pad SNV
        start = 9
        stop = 14
        ref = "ACACA"
        alt = "TCACA"

        trim_start, trim_stop, trim_alleles = norm.normalize_alleles(
            seq,
            start,
            stop,
            (ref, alt)
        )
        assert trim_alleles == ("A", "T")
        assert trim_start == 9
        assert trim_stop == 10

        # left pad ins
        start = 9
        stop = 12
        ref = "ACA"
        alt = "ACAGT"

        trim_start, trim_stop, trim_alleles = norm.normalize_alleles(
            seq,
            start,
            stop,
            (ref, alt)
        )

        assert trim_alleles == ("", "GT")
        assert trim_start == 12
        assert trim_stop == 12

        # left-align pad del
        start = 9
        stop = 14
        ref = "ACACA"
        alt = "ACA"

        trim_start, trim_stop, trim_alleles = norm.normalize_alleles(
            seq,
            start,
            stop,
            (ref, alt)
        )
        assert trim_alleles == ("AC", "")
        assert trim_start == 9
        assert trim_stop == 11
