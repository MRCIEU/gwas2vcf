import argparse
import logging
import marshmallow
from gwas_results import GwasResult
from vcf import Vcf
import pysam
from harmonise import Harmonise
from datetime import datetime
import json
from param import Param
import sys


def main():
    version = "1.1.1"

    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF/BCF')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(version))
    parser.add_argument('--out', dest='out', required=True, help='Path to output VCF/BCF')
    parser.add_argument('--data', dest='gwas', required=True, help='Path to GWAS summary stats')
    parser.add_argument('--ref', dest='fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('--json', dest='json', required=True, help='Path to parameters JSON')
    parser.add_argument('--id', dest='id', required=True, help='Study identifier')
    parser.add_argument('--cohort_controls', type=int, dest='cohort_controls', required=False, default=None,
                        help='Total study number of controls (if case/control) or total sample size if continuous')
    parser.add_argument('--cohort_cases', type=int, dest='cohort_cases', required=False, default=None,
                        help='Total study number of cases')
    parser.add_argument('--rm_chr_prefix', dest='rm_chr_prefix', action='store_true', default=False, required=False,
                        help='Remove chr prefix from GWAS chromosome')
    parser.add_argument('--csi', dest='csi', action='store_true', default=False, required=False,
                        help='Default is to index tbi but use this flag to index csi')
    parser.add_argument("--log", dest="log", required=False, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging level")
    args = parser.parse_args()

    # set logging level
    if args.log:
        logging.basicConfig(level=getattr(logging, args.log), format='%(asctime)s %(levelname)s %(message)s')

    logging.info("GWAS2VCF {}".format(version))

    # check values are valid
    if args.cohort_cases is not None:
        if args.cohort_cases < 1:
            logging.error("Total study number of cases must be a positive number")
            sys.exit()

    if args.cohort_controls is not None:
        if args.cohort_controls < 1:
            logging.error("Total study number of controls must be a positive number")
            sys.exit()

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
        rm_chr_prefix=args.rm_chr_prefix
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

        file_metadata = {
            'gwas_harmonisation_command': ' '.join(sys.argv[1:]) + "; " + version,
            'file_date': datetime.now().isoformat()
        }

        sample_metadata = {
            'TotalVariants': total_variants,
            'VariantsNotRead': total_variants - len(gwas),
            'HarmonisedVariants': len(harmonised),
            'VariantsNotHarmonised': len(gwas) - len(harmonised),
            'SwitchedAlleles': flipped_variants - (len(gwas) - len(harmonised))
        }

        if args.cohort_controls is not None:
            sample_metadata['TotalControls'] = args.cohort_controls

        if args.cohort_cases is not None:
            sample_metadata['TotalCases'] = args.cohort_cases

        if 'ncase_col' in j or args.cohort_cases is not None:
            sample_metadata['StudyType'] = 'CaseControl'
        else:
            sample_metadata['StudyType'] = 'Continuous'

        # write to vcf
        Vcf.write_to_file(harmonised, args.out, fasta, j['build'], args.id, sample_metadata, file_metadata, args.csi)


if __name__ == "__main__":
    main()
