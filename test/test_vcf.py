from vcf import Vcf
import numpy as np
import pickle
import decimal
from gwas import Gwas
import tempfile
from heapq import heappush, heappop
from pvalueHandler import PvalueHandler

def test_convert_pval_to_neg_log10():
    p_value_handler = PvalueHandler()
    assert p_value_handler.neg_log_of_decimal(decimal.Decimal(1)) == 0
    assert p_value_handler.neg_log_of_decimal(decimal.Decimal(0)) == 999


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
    p_value_handler = PvalueHandler()
    precision_check = {"string":"1E-10000","nlog_value":10000}
    g1 = [Gwas("1", 101, "A", "T", 1, 0, 5e-8, 1000, 0.4, "rs1234", None, None, None),
          Gwas("1", 105, "A", "T", 1, 0, 5e-8, 1000, 0.4, "rs1234", None, None, None),
          Gwas("1", 102, "A", "T", 1, 0, 5e-8, 1000, 0.4, "rs1234", None, None, None),
          Gwas("1", 103, "A", "T", 1, 0, p_value_handler.neg_log_of_decimal(p_value_handler.parse_string(precision_check ["string"])),
               1000, 0.4, "rs1234", None, None, None)
          ]
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

    epsilon = g2[2].nlog_pval/1000000.0
    assert g2[0].pos == 101
    assert g2[1].pos == 102
    assert g2[2].pos == 103
    assert g2[2].nlog_pval > precision_check ["nlog_value"]- epsilon
    assert g2[2].nlog_pval < precision_check ["nlog_value"]+epsilon
    assert g2[3].pos == 105

    results.close()
