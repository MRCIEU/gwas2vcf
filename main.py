import argparse
import logging
import marshmallow
from gwas_results import GwasResult
from vcf import Vcf
import pysam
from harmonise import Harmonise
from datetime import datetime
import git
import os
import json
from param import Param
import sys


def main():
    repo = git.Repo(os.path.dirname(os.path.realpath(__file__)))
    sha = repo.head.object.hexsha

    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF/BCF')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(sha))
    parser.add_argument('--out', dest='out', required=True, help='Path to output VCF/BCF')
    parser.add_argument('--data', dest='gwas', required=True, help='Path to GWAS summary stats')
    parser.add_argument('--ref', dest='fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('--json', dest='json', required=True, help='Path to parameters JSON')
    parser.add_argument('--id', dest='id', default=None, required=False, help='GWAS identifier')
    parser.add_argument('--cohort_sample_size', type=int, dest='cohort_sample_size', required=False, default=None,
                        help='Total sample size of cohort')
    parser.add_argument('--cohort_frac_cases', type=float, dest='cohort_frac_cases', required=False, default=None,
                        help='Total fraction of cases in cohort')
    parser.add_argument("--log", dest="log", required=False, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging level")
    args = parser.parse_args()

    # set logging level
    if args.log:
        logging.basicConfig(level=getattr(logging, args.log), format='%(asctime)s %(levelname)s %(message)s')

    logging.info("GWAS Harmonisation {}".format(sha))

    # load parameters from json
    logging.info("Reading JSON parameters")
    try:
        schema = Param(strict=True)
        with open(args.json) as f:
            j = schema.load(json.load(f)).data
            logging.info("Parameters: {}".format(j))
    except json.decoder.JSONDecodeError as e:
        logging.error("Could not read json parameter file: {}".format(e))
        sys.exit()
    except marshmallow.exceptions.ValidationError as e:
        logging.error("Could not validate json parameter file: {}".format(e))
        sys.exit()

    # read in GWAS and harmonise alleles to reference fasta
    gwas, total_variants = GwasResult.read_from_file(
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
        ncontrol_field=j.get('ncontrol_col'),
        cohort_sample_size=args.cohort_sample_size,
        cohort_frac_cases=args.cohort_frac_cases
    )

    logging.info("Total variants: {}".format(total_variants))
    logging.info("Variants could not be read: {}".format(total_variants - len(gwas)))

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
            'counts.switched_alleles': flipped_variants - (len(gwas) - len(harmonised)),
            'cohort_sample_size': args.cohort_sample_size,
            'cohort_frac_cases': args.cohort_frac_cases,
            'gwas.id': args.id
        }

        if 'ncase_col' in j:
            params['gwas.type'] = 'case/control'
        else:
            params['gwas.type'] = 'continuous'

        # write to vcf
        Vcf.write_to_file(harmonised, args.out, fasta, j['build'], params)


if __name__ == "__main__":
    main()
