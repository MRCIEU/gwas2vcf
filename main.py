import argparse
import logging
from gwas_results import GwasResult
from vcf import Vcf

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


def main():
    # TODO make skip lines optional
    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF')
    parser.add_argument('-o', dest='out', required=True, help='Path to output VCF')
    parser.add_argument('-g', dest='gwas', required=True, help='Path to GWAS summary stats')
    parser.add_argument('-f', dest='fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('-s', dest='skip', required=True, help='Number of rows to skip')
    args = parser.parse_args()

    # read in GWAS and harmonise alleles to reference fasta
    harmonised = []
    for r in GwasResult.read_from_text_file(args.gwas, skip_n_rows=int(args.skip)):
        if not r.is_ref_allele_match_fasta_htslib(args.fasta):
            r.reverse_sign()
            if not r.is_ref_allele_match_fasta_htslib(args.fasta):
                raise ValueError("Could not map record to FASTA {}".format(r))
        harmonised.append(r)

    # write to vcf
    Vcf.write_to_file(harmonised, args.out)


if __name__ == "__main__":
    main()
