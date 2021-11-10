import decimal
import os

import pysam
import pytest
from vgraph import norm

from gwas import Gwas
from pvalue_handler import PvalueHandler


def test_are_alleles_iupac():
    g = Gwas("test", 1, "A", "T", None, None, None, None, None, None, None, None, None)
    g.check_alleles_are_valid()

    with pytest.raises(AssertionError):
        g = Gwas(
            "test", 1, "A", "wdeT", None, None, None, None, None, None, None, None, None
        )
        g.check_alleles_are_valid()


def test_reverse_sign():
    g = Gwas("test", 1, "A", "T", 1, None, None, None, None, None, None, None, None)
    g.reverse_sign()
    assert g.chrom == "test"
    assert g.pos == 1
    assert g.ref == "T"
    assert g.alt == "A"
    assert g.b == -1


def test_update_dbsnp():
    with pysam.VariantFile(
        os.path.join(os.path.dirname(__file__), "dbsnp.vcf.gz")
    ) as dbsnp:
        g = Gwas("test", 1, "A", "T", 1, None, None, None, None, None, None, None, None)
        assert g.dbsnpid is None
        g.update_dbsnp(dbsnp)
        assert g.dbsnpid == "rs1234"
        g = Gwas("test", 2, "A", "T", 1, None, None, None, None, None, None, None, None)
        assert g.dbsnpid is None
        g.update_dbsnp(dbsnp)
        assert g.dbsnpid is None


def test_check_reference_allele():
    with pysam.FastaFile(
        os.path.join(os.path.dirname(__file__), "test.fasta")
    ) as fasta:
        g = Gwas("test", 1, "A", "T", 1, None, None, None, None, None, None, None, None)
        g.check_reference_allele(fasta)

        with pytest.raises(AssertionError):
            g = Gwas(
                "test", 1, "T", "A", 1, None, None, None, None, None, None, None, None
            )
            g.check_reference_allele(fasta)


def test_pvalues_precision_with_PvalueHandler():
    p_value_handler = PvalueHandler()
    p_value_strings = ["0.1", "0.00000000005", "7.5E-8", "8e-1000", "0.1E-10000"]
    gwas_list = []
    for p_value_string in p_value_strings:
        gwas_list.append(
            Gwas(
                "test",
                1,
                "A",
                "T",
                None,
                None,
                p_value_handler.neg_log_of_decimal(
                    p_value_handler.parse_string(p_value_string)
                ),
                None,
                None,
                None,
                None,
                None,
                None,
            )
        )

    # show that we get the number we expect when converted to a float
    assert gwas_list[0].nlog_pval == float(-decimal.Decimal(p_value_strings[0]).log10())
    assert gwas_list[1].nlog_pval == float(-decimal.Decimal(p_value_strings[1]).log10())
    assert gwas_list[2].nlog_pval == float(-decimal.Decimal(p_value_strings[2]).log10())
    assert gwas_list[3].nlog_pval == float(-decimal.Decimal(p_value_strings[3]).log10())
    assert gwas_list[4].nlog_pval == float(-decimal.Decimal(p_value_strings[4]).log10())

    # show that the show that the resulting flow is at least close to the exact decimal value
    epsilon = gwas_list[0].nlog_pval / 1000000.0
    assert (
        gwas_list[0].nlog_pval + epsilon > -decimal.Decimal(p_value_strings[0]).log10()
        and gwas_list[0].nlog_pval - epsilon
        < -decimal.Decimal(p_value_strings[0]).log10()
    )
    epsilon = gwas_list[1].nlog_pval / 1000000.0
    assert (
        gwas_list[1].nlog_pval + epsilon > -decimal.Decimal(p_value_strings[1]).log10()
        and gwas_list[1].nlog_pval - epsilon
        < -decimal.Decimal(p_value_strings[1]).log10()
    )
    epsilon = gwas_list[2].nlog_pval / 1000000.0
    assert (
        gwas_list[2].nlog_pval + epsilon > -decimal.Decimal(p_value_strings[2]).log10()
        and gwas_list[2].nlog_pval - epsilon
        < -decimal.Decimal(p_value_strings[2]).log10()
    )
    epsilon = gwas_list[3].nlog_pval / 1000000.0
    assert (
        gwas_list[3].nlog_pval + epsilon > -decimal.Decimal(p_value_strings[3]).log10()
        and gwas_list[3].nlog_pval - epsilon
        < -decimal.Decimal(p_value_strings[3]).log10()
    )
    epsilon = gwas_list[4].nlog_pval / 1000000.0
    assert (
        gwas_list[4].nlog_pval + epsilon > -decimal.Decimal(p_value_strings[4]).log10()
        and gwas_list[4].nlog_pval - epsilon
        < -decimal.Decimal(p_value_strings[4]).log10()
    )

    # Note:this is the problem we're trying to correct
    p = "7e-1000"
    pv = float(p)
    assert pv == 0.0

    # so that show that this problem is alleviated
    with pytest.raises(AssertionError):
        assert gwas_list[0].nlog_pval == 0.0
        assert gwas_list[1].nlog_pval == 0.0
        assert gwas_list[2].nlog_pval == 0.0
        assert gwas_list[3].nlog_pval == 0.0
        assert gwas_list[4].nlog_pval == 0.0


def test_normalise():
    with pysam.FastaFile(
        os.path.join(os.path.dirname(__file__), "test.fasta")
    ) as fasta:
        # SNV
        g = Gwas(
            "test", 1, "A", "T", None, None, None, None, None, None, None, None, None
        )
        g.check_reference_allele(fasta)
        g.normalise(fasta, padding=5)
        assert g.chrom == "test"
        assert g.pos == 1
        assert g.ref == "A"
        assert g.alt == "T"

        # left pad SNV
        g = Gwas(
            "test",
            10,
            "ACACA",
            "ACACT",
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        g.check_reference_allele(fasta)
        g.normalise(fasta, padding=5)
        assert g.chrom == "test"
        assert g.pos == 14
        assert g.ref == "A"
        assert g.alt == "T"

        # right pad SNV
        g = Gwas(
            "test",
            10,
            "ACACA",
            "TCACA",
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        g.check_reference_allele(fasta)
        g.normalise(fasta, padding=5)
        assert g.chrom == "test"
        assert g.pos == 10
        assert g.ref == "A"
        assert g.alt == "T"

        # left pad ins
        g = Gwas(
            "test",
            10,
            "ACA",
            "ACAGT",
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        g.check_reference_allele(fasta)
        g.normalise(fasta, padding=5)
        assert g.chrom == "test"
        assert g.pos == 12
        assert g.ref == "A"
        assert g.alt == "AGT"

        # left-align pad del
        g = Gwas(
            "test",
            10,
            "ACACA",
            "ACA",
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        g.check_reference_allele(fasta)
        g.normalise(fasta, padding=5)
        assert g.chrom == "test"
        assert g.pos == 9
        assert g.ref == "TAC"
        assert g.alt == "T"


def test_bioinformed_vgraph_normalise():
    """Cython logic from bioinformed/vgraph"""

    # get a reference sequence
    with pysam.FastaFile(
        os.path.join(os.path.dirname(__file__), "test.fasta")
    ) as fasta:
        seq = fasta.fetch("test").upper()

        # SNV

        # use 0-based
        start = 0
        # use 1-based
        stop = 1
        ref = "A"
        alt = "T"

        trim_start, trim_stop, trim_alleles = norm.normalize_alleles(
            seq, start, stop, (ref, alt)
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
            seq, start, stop, (ref, alt)
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
            seq, start, stop, (ref, alt)
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
            seq, start, stop, (ref, alt)
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
            seq, start, stop, (ref, alt)
        )
        assert trim_alleles == ("AC", "")
        assert trim_start == 9
        assert trim_stop == 11
