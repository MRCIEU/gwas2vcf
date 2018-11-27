import argparse
import logging
from gwas_results import GwasResult
from vcf import Vcf

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


def main():
    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF')
    parser.add_argument('-o', dest='out', required=True, help='Path to output VCF')
    parser.add_argument('-g', dest='gwas', required=True, help='Path to GWAS summary stats')
    parser.add_argument('-f', dest='fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('-s', dest='skip', required=False, help='Number of rows to skip')
    parser.add_argument('-chrom_field', dest='chrom_field', required=False, help='Column number for chromosome')
    parser.add_argument('-pos_field', dest='pos_field', required=False, help='Column number for chromosome')
    parser.add_argument('-a1_field', dest='a1_field', required=True, help='Column number for allele 1')
    parser.add_argument('-a2_field', dest='a2_field', required=True, help='Column number for allele 2')
    parser.add_argument('-effect_field', dest='effect_field', required=True, help='Effect size field')
    parser.add_argument('-dbsnp_field', dest='dbsnp_field', required=False, help='dbSNP identifier field')
    parser.add_argument('-se_field', dest='se_field', required=True, help='SE field')
    parser.add_argument('-n0_field', dest='n0_field', required=False, help='N0 field')
    parser.add_argument('-n1_field', dest='n1_field', required=False, help='N1 field')
    parser.add_argument('-pval_field', dest='pval_field', required=True, help='P-Value field')
    parser.add_argument('-a2_af_field', dest='a2_af_field', required=True, help='Allele 2 frequency field')
    args = parser.parse_args()

    # count skipped lines
    excluded_variants = 0

    try:
        skip_lines = int(args.skip)
    except Exception:
        skip_lines = 0

    # read in GWAS and harmonise alleles to reference fasta
    harmonised = []
    orig = GwasResult.read_from_text_file(
        args.gwas,
        args.effect_field,
        args.se_field,
        args.pval_field,
        dbsnp_field=args.dbsnp_field,
        chrom_field=args.chrom_field,
        pos_field=args.pos_field,
        a1_field=args.a1_field,
        a2_field=args.a2_field,
        n0_field=args.n0_field,
        n1_field=args.n1_field,
        a2_af_field=args.a2_af_field,
        skip_n_rows=skip_lines)
    for r in orig:
        if not r.is_ref_allele_match_fasta_htslib(args.fasta):
            r.reverse_sign()
            if not r.is_ref_allele_match_fasta_htslib(args.fasta):
                logging.warning("Skipping record {}: could not match to FASTA".format(r))
                excluded_variants += 1
                continue
        harmonised.append(r)

    # check number of skipped is acceptable
    logging.info("Skipped {} of {}".format(excluded_variants, len(orig)))
    if excluded_variants / len(orig) > 0.05:
        raise RuntimeError("Too many sites skipped.")

    # write to vcf
    Vcf.write_to_file(harmonised, args.out)


if __name__ == "__main__":
    main()
