import unittest
from gwas_results import GwasResult


class TestGwasResults(unittest.TestCase):

    def test_read_from_text_file(self):
        harmonised = []
        for r in GwasResult.read_from_text_file('data/jointGwasMc_LDL.n100.txt', skip_n_rows=1):
            if not r.is_ref_allele_match_fasta('/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta'):
                r.reverse_sign()
                if not r.is_ref_allele_match_fasta('/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta'):
                    raise ValueError("Could not map record to FASTA {}".format(r))
            harmonised.append(r)


if __name__ == '__main__':
    unittest.main()
