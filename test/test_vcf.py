from vcf import Vcf
import numpy as np


def test_convert_pval_to_neg_log10():
    assert Vcf.convert_pval_to_neg_log10(1) == 0
    assert Vcf.convert_pval_to_neg_log10(0) == 999


def test_is_valid_float32():
    # small
    assert not Vcf.is_float32_lossy(1e-37)
    assert not Vcf.is_float32_lossy(np.finfo(np.float32).tiny)
    assert Vcf.is_float32_lossy(1e-50)
    assert not Vcf.is_float32_lossy(-1e-37)
    assert Vcf.is_float32_lossy(-1e-50)

    # large
    assert Vcf.is_float32_lossy(999999999999999999999999999999999999999)
    assert Vcf.is_float32_lossy(-999999999999999999999999999999999999999)
    assert not Vcf.is_float32_lossy(-99999999999999999999999999999999999999)
    assert not Vcf.is_float32_lossy(99999999999999999999999999999999999999)
