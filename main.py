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
import json

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


def main():
    repo = git.Repo(os.path.dirname(os.path.realpath(__file__)))
    sha = repo.head.object.hexsha

    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF/BCF')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(sha))
    parser.add_argument('--out', dest='out', required=True, help='Path to output VCF/BCF')
    parser.add_argument('--data', dest='gwas', required=True, help='Path to GWAS summary stats')
    parser.add_argument('--ref', dest='fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('--json', dest='json', required=True, help='Path to parameters JSON')
    args = parser.parse_args()

    # load parameters from json
    with open(args.json) as f:
        j = json.load(f)

    # read in GWAS and harmonise alleles to reference fasta
    gwas, total_variants = GwasResult.read_from_text_file(
        args.gwas,
        j['chr_col'],
        j['pos_col'],
        j['ea_col'],
        j['oa_col'],
        j['beta_col'],
        j['se_col'],
        j['pval_col'],
        j['delimiter'],
        j['header'],
        ncase_field=j.get('ncase_col'),
        dbsnp_field=j.get('snp_col'),
        ea_af_field=j.get('eaf_col'),
        nea_af_field=j.get('oaf_col'),
        imp_z_field=j.get('imp_z_col'),
        imp_info_field=j.get('imp_info_col'),
        ncontrol_field=j.get('ncontrol_col')
    )

    logging.info("Total variants: {}".format(total_variants))
    logging.info("Variants could not be read: {}".format(total_variants - len(gwas)))

    # print first lines for debugging
    for i in range(10):
        try:
            logging.info("Mapped line {}: {}".format(i, gwas[i]))
        except IndexError:
            continue

    with pysam.FastaFile(args.fasta) as fasta:

        # harmonise to FASTA
        harmonised, flipped_variants = Harmonise.align_gwas_to_fasta(gwas, fasta)

        logging.info("Variants harmonised: {}".format(len(harmonised)))
        logging.info("Variants discarded during harmonisation: {}".format(len(gwas) - len(harmonised)))
        logging.info("Alleles switched: {}".format(flipped_variants - (len(gwas) - len(harmonised))))

        # number of skipped records
        logging.info("Skipped {} of {}".format(total_variants - len(harmonised), total_variants))

        params = {
            'gwas_harmonisation_command': ' '.join(sys.argv[1:]) + "; " + sha,
            'file_date': datetime.now().isoformat(),
            'counts.total_variants': total_variants,
            'counts.variants_not_read': total_variants - len(gwas),
            'counts.harmonised_variants': len(harmonised),
            'counts.variants_not_harmonised': len(gwas) - len(harmonised),
            'counts.switched_alleles': flipped_variants - (len(gwas) - len(harmonised))
        }

        if 'id' in j:
            params['gwas.id'] = j['id']

        if 'ncase_col' in j:
            params['gwas.type'] = 'case/control'
        else:
            params['gwas.type'] = 'continuous'

        # write to vcf
        Vcf.write_to_file(harmonised, args.out, fasta, j['build'], params)


if __name__ == "__main__":
    main()
