import argparse
import logging
from gwas_results import GwasResult
from vcf import Vcf
import pysam
from harmonise import Harmonise
import json
import os

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
    parser.add_argument('-n_field', dest='n_field', type=int, required=False, help='Number of samples field')
    parser.add_argument('-dbsnp_field', dest='dbsnp_field', type=int, required=False, help='dbSNP identifier field')
    parser.add_argument('-ea_af_field', dest='ea_af_field', type=int, required=False,
                        help='Effect allele frequency field')
    parser.add_argument('-nea_af_field', dest='nea_af_field', type=int, required=False,
                        help='None effect allele frequency field')
    args = parser.parse_args()

    # load fasta fai
    fasta = pysam.FastaFile(args.fasta)

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

    # harmonise to FASTA
    harmonised, flipped_variants = Harmonise.align_gwas_to_fasta(gwas, fasta)

    logging.info("Variants harmonised: {}".format(len(harmonised)))
    logging.info("Variants discarded during harmonisation: {}".format(len(gwas) - len(harmonised)))
    logging.info("Alleles switched: {}".format(flipped_variants - (len(gwas) - len(harmonised))))

    # check number of skipped is acceptable
    logging.info("Skipped {} of {}".format(total_variants - len(harmonised), total_variants))
    if (total_variants - len(harmonised)) / total_variants > args.max_missing:
        raise RuntimeError("Too many sites skipped.")

    # write to vcf
    Vcf.write_to_file(harmonised, args.out, fasta)

    # write metrics to json
    metrics = {
        'args': dict(),
        'counts': dict()
    }
    name, ext = os.path.splitext(args.out)

    # add args
    for arg in vars(args):
        metrics['args'][arg] = getattr(args, arg)

    # add counts
    metrics['counts']['total_variants'] = total_variants
    metrics['counts']['variants_not_read'] = total_variants - len(gwas)
    metrics['counts']['harmonised_variants'] = len(harmonised)
    metrics['counts']['variants_not_harmonised'] = len(gwas) - len(harmonised)
    metrics['counts']['switched_alleles'] = flipped_variants - (len(gwas) - len(harmonised))

    # write to file
    with open(name + ".json", 'w') as f:
        json.dump(metrics, f, ensure_ascii=False)

    fasta.close()


if __name__ == "__main__":
    main()
