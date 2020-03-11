from vcf import Vcf
import numpy as np
import pickle
from gwas import Gwas
import tempfile
from heapq import heappush, heappop


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


def test_gwas_unpickle():
    g1 = [Gwas("1", 101, "A", "T", 1, 0, 5e-8, 1000, 0.4, "rs1234", None, None, None),
          Gwas("1", 105, "A", "T", 1, 0, 5e-8, 1000, 0.4, "rs1234", None, None, None),
          Gwas("1", 102, "A", "T", 1, 0, 5e-8, 1000, 0.4, "rs1234", None, None, None)]
    g2 = []
    idx = []
    results = tempfile.TemporaryFile()
    for result in g1:
        heappush(idx, (result.pos, results.tell()))
        pickle.dump(result, results)

    while idx:
        pos = heappop(idx)
        results.seek(pos[1])
        result = pickle.load(results)
        g2.append(result)

    assert g2[0].pos == 101
    assert g2[1].pos == 102
    assert g2[2].pos == 105

    results.close()
