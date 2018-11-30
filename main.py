import argparse
import logging
from gwas_results import GwasResult
from vcf import Vcf
import pysam
from harmonise import Harmonise
from datetime import datetime
import git
import sys
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


def main():
    repo = git.Repo(os.path.dirname(os.path.realpath(__file__)))
    sha = repo.head.object.hexsha

    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(sha))
    parser.add_argument('-o', dest='out', required=True, help='Path to output VCF')
    parser.add_argument('-g', dest='gwas', required=True, help='Path to GWAS tab-sep summary stats')
    parser.add_argument('-f', dest='fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('-b', dest='build', required=True, help='FASTA build assembly')
    parser.add_argument('-s', dest='skip', default=0, type=int, required=False, help='Number of rows to skip')
    parser.add_argument('-m', dest='max_missing', default=0.05, type=float, required=False,
                        help='Maximum fraction of missing variants permitted')
    parser.add_argument('-chrom_field', dest='chrom_field', type=int, required=True,
                        help='Column number for chromosome')
    parser.add_argument('-pos_field', dest='pos_field', type=int, required=True, help='Column number for chromosome')
    parser.add_argument('-ea_field', dest='ea_field', type=int, required=True, help='Column number for effect allele')
    parser.add_argument('-nea_field', dest='nea_field', type=int, required=True,
                        help='Column number for non-effect allele')
    parser.add_argument('-effect_field', dest='effect_field', type=int, required=True, help='Effect size field')
    parser.add_argument('-se_field', dest='se_field', type=int, required=True, help='SE field')
    parser.add_argument('-pval_field', dest='pval_field', type=int, required=True, help='P-Value field')
    parser.add_argument('-n_field', dest='n_field', type=int, required=False, help='Number of samples field')
    parser.add_argument('-dbsnp_field', dest='dbsnp_field', type=int, required=False, help='dbSNP identifier field')
    parser.add_argument('-ea_af_field', dest='ea_af_field', type=int, required=False,
                        help='Effect allele frequency field')
    parser.add_argument('-nea_af_field', dest='nea_af_field', type=int, required=False,
                        help='None effect allele frequency field')
    args = parser.parse_args()

    # read in GWAS and harmonise alleles to reference fasta
    gwas, total_variants = GwasResult.read_from_text_file(
        args.gwas,
        args.chrom_field,
        args.pos_field,
        args.ea_field,
        args.nea_field,
        args.effect_field,
        args.se_field,
        args.pval_field,
        dbsnp_field=args.dbsnp_field,
        n_field=args.n_field,
        ea_af_field=args.ea_af_field,
        nea_af_field=args.nea_af_field,
        skip_n_rows=args.skip)

    logging.info("Total variants: {}".format(total_variants))
    logging.info("Variants could not be read: {}".format(total_variants - len(gwas)))

    # print first lines for debugging
    for i in range(10):
        logging.info("Mapped line {}: {}".format(i, gwas[i]))

    with pysam.FastaFile(args.fasta) as fasta:

        # harmonise to FASTA
        harmonised, flipped_variants = Harmonise.align_gwas_to_fasta(gwas, fasta)

        logging.info("Variants harmonised: {}".format(len(harmonised)))
        logging.info("Variants discarded during harmonisation: {}".format(len(gwas) - len(harmonised)))
        logging.info("Alleles switched: {}".format(flipped_variants - (len(gwas) - len(harmonised))))

        # check number of skipped is acceptable
        logging.info("Skipped {} of {}".format(total_variants - len(harmonised), total_variants))
        if (total_variants - len(harmonised)) / total_variants > args.max_missing:
            raise RuntimeError("Too many sites skipped. The alleles must be on the forward strand.")

        # write to vcf
        Vcf.write_to_file(harmonised, args.out, fasta, args.build, {
            'gwas_harmonisation_command': ' '.join(sys.argv[1:]),
            'file_date': datetime.now().isoformat(),
            'total_variants': total_variants,
            'variants_not_read': total_variants - len(gwas),
            'harmonised_variants': len(harmonised),
            'variants_not_harmonised': len(gwas) - len(harmonised),
            'switched_alleles': flipped_variants - (len(gwas) - len(harmonised))
        })


if __name__ == "__main__":
    main()
