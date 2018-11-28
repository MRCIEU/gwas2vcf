import argparse
import logging
from gwas_results import GwasResult
from vcf import Vcf
import pysam
from harmonise import Harmonise

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


def main():
    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF')
    parser.add_argument('-o', dest='out', required=True, help='Path to output VCF')
    parser.add_argument('-g', dest='gwas', required=True, help='Path to GWAS tab-sep summary stats')
    parser.add_argument('-f', dest='fasta', required=True, help='Path to reference FASTA')
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
    parser.add_argument('-n0_field', dest='n0_field', type=int, required=False, help='N0 field')
    parser.add_argument('-dbsnp_field', dest='dbsnp_field', type=int, required=False, help='dbSNP identifier field')
    parser.add_argument('-n1_field', dest='n1_field', type=int, required=False, help='N1 field')
    parser.add_argument('-ea_af_field', dest='ea_af_field', type=int, required=False,
                        help='Effect allele frequency field')
    parser.add_argument('-nea_af_field', dest='nea_af_field', type=int, required=False,
                        help='None effect allele frequency field')
    args = parser.parse_args()

    # load fasta fai
    fasta = pysam.FastaFile(args.fasta)

    # read in GWAS and harmonise alleles to reference fasta
    gwas = GwasResult.read_from_text_file(
        args.gwas,
        args.chrom_field,
        args.pos_field,
        args.ea_field,
        args.nea_field,
        args.effect_field,
        args.se_field,
        args.pval_field,
        dbsnp_field=args.dbsnp_field,
        n0_field=args.n0_field,
        n1_field=args.n1_field,
        ea_af_field=args.ea_af_field,
        nea_af_field=args.nea_af_field,
        skip_n_rows=args.skip)

    # harmonise to FASTA
    harmonised, excluded_variants = Harmonise.align_gwas_to_fasta(gwas, fasta)

    # check number of skipped is acceptable
    logging.info("Skipped {} of {}".format(excluded_variants, len(gwas)))
    if excluded_variants / len(gwas) > args.max_missing:
        raise RuntimeError("Too many sites skipped.")

    # write to vcf
    Vcf.write_to_file(harmonised, args.out, fasta)

    fasta.close()


if __name__ == "__main__":
    main()
