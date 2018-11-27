import argparse
import logging
from gwas_results import GwasResult
from vcf import Vcf
import pysam

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
    parser.add_argument('-a1_field', dest='a1_field', type=int, required=True, help='Column number for allele 1')
    parser.add_argument('-a2_field', dest='a2_field', type=int, required=True, help='Column number for allele 2')
    parser.add_argument('-effect_field', dest='effect_field', type=int, required=True, help='Effect size field')
    parser.add_argument('-se_field', dest='se_field', type=int, required=True, help='SE field')
    parser.add_argument('-pval_field', dest='pval_field', type=int, required=True, help='P-Value field')
    parser.add_argument('-n0_field', dest='n0_field', type=int, required=True, help='N0 field')
    parser.add_argument('-dbsnp_field', dest='dbsnp_field', type=int, required=False, help='dbSNP identifier field')
    parser.add_argument('-n1_field', dest='n1_field', type=int, required=False, help='N1 field')
    parser.add_argument('-a1_af_field', dest='a1_af_field', type=int, required=False, help='Allele 1 frequency field')
    parser.add_argument('-a2_af_field', dest='a2_af_field', type=int, required=False, help='Allele 2 frequency field')
    args = parser.parse_args()

    # load fasta fai
    fasta = pysam.FastaFile(args.fasta)

    # count skipped lines
    excluded_variants = 0

    # read in GWAS and harmonise alleles to reference fasta
    harmonised = []
    orig = GwasResult.read_from_text_file(
        args.gwas,
        args.chrom_field,
        args.pos_field,
        args.a1_field,
        args.a2_field,
        args.effect_field,
        args.se_field,
        args.pval_field,
        args.n0_field,
        dbsnp_field=args.dbsnp_field,
        n1_field=args.n1_field,
        a1_af_field=args.a1_af_field,
        a2_af_field=args.a2_af_field,
        skip_n_rows=args.skip)

    for r in orig:

        if not r.are_alleles_iupac():
            logging.warning("Skipping record {}: allele(s) are not standard IUPAC".format(r))
            excluded_variants += 1
            continue

        # get expected FASTA REF
        # TODO add support for indels
        expected_ref_allele = str(fasta.fetch(region="{}:{}-{}".format(r.chrom, r.pos, r.pos))).upper()

        if r.ref != expected_ref_allele:
            r.reverse_sign()
            if r.ref != expected_ref_allele:
                logging.warning("Skipping record {}: could not match to FASTA {}".format(r, expected_ref_allele))
                excluded_variants += 1
                continue

        harmonised.append(r)

    # check number of skipped is acceptable
    logging.info("Skipped {} of {}".format(excluded_variants, len(orig)))
    if excluded_variants / len(orig) > args.max_missing:
        raise RuntimeError("Too many sites skipped.")

    # write to vcf
    Vcf.write_to_file(harmonised, args.out)

    fasta.close()


if __name__ == "__main__":
    main()
