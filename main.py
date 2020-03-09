import argparse
import logging
import marshmallow
from gwas import Gwas
from vcf import Vcf
import pysam
from datetime import datetime
import json
from param import Param
import sys
import os


def main():
    version = "1.2.0"

    parser = argparse.ArgumentParser(description='Map GWAS summary statistics to VCF/BCF')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(version))
    parser.add_argument('--out', dest='out', required=False,
                        help='Path to output VCF/BCF. If not present then must be specified as \'out\' in json file')
    parser.add_argument('--data', dest='data', required=False,
                        help='Path to GWAS summary stats. If not present then must be specified as \'data\' in json file')
    parser.add_argument('--ref', dest='ref', required=True, help='Path to reference FASTA')
    parser.add_argument('--json', dest='json', required=True, help='Path to parameters JSON')
    parser.add_argument('--id', dest='id', required=False,
                        help='Study identifier. If not present then must be specified as \'id\' in json file')
    parser.add_argument('--cohort_controls', type=int, dest='cohort_controls', required=False, default=None,
                        help='Total study number of controls (if case/control) or total sample size if continuous. Overwrites value if present in json file.')
    parser.add_argument('--cohort_cases', type=int, dest='cohort_cases', required=False, default=None,
                        help='Total study number of cases. Overwrites value if present in json file.')
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

    logging.info("Gwas2VCF {}".format(version))
    logging.info("Arguments: {}".format(vars(args)))
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

    logging.info("Checking input arguments")
    if args.data is None:
        if 'data' in j.keys():
            vars(args)['data'] = j['data']
        else:
            logging.error("'data' filename not provided in arguments or json file")
            sys.exit()

    if args.out is None:
        if 'out' in j.keys():
            vars(args)['out'] = j['out']
        else:
            logging.error("out filename not provided in arguments or json file")
            sys.exit()

    if args.id is None:
        if 'id' in j.keys():
            vars(args)['id'] = j['id']
        else:
            logging.error("id not provided in arguments or json file")
            sys.exit()

    if args.cohort_cases is None and 'cohort_cases' in j.keys():
        vars(args)['cohort_cases'] = j['cohort_cases']

    if args.cohort_controls is None and 'cohort_controls' in j.keys():
        vars(args)['cohort_controls'] = j['cohort_controls']

    # check values are valid
    if args.cohort_cases is not None:
        if args.cohort_cases < 1:
            logging.error("Total study number of cases must be a positive number")
            sys.exit()

    if args.cohort_controls is not None:
        if args.cohort_controls < 1:
            logging.error("Total study number of controls must be a positive number")
            sys.exit()

    if not os.path.isfile(args.data):
        logging.error("{} file does not exist".format(args.data))
        sys.exit()

    if not os.path.isfile(args.ref):
        logging.error("{} file does not exist".format(args.ref))
        sys.exit()

    if not os.path.exists(os.path.dirname(args.out)):
        logging.error("{} output directory does not exist".format(args.out))
        sys.exit()

    # read in data
    # harmonise on-the-fly and write to pickle format
    # keep file index for each record and chromosome position to write out karyotypically sorted records later
    with pysam.FastaFile(args.ref) as fasta:
        gwas, idx, sample_metadata = Gwas.read_from_file(
            args.data,
            fasta,
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
            rsid_field=j.get('snp_col'),
            ea_af_field=j.get('eaf_col'),
            nea_af_field=j.get('oaf_col'),
            imp_z_field=j.get('imp_z_col'),
            imp_info_field=j.get('imp_info_col'),
            ncontrol_field=j.get('ncontrol_col'),
            rm_chr_prefix=args.rm_chr_prefix
        )

    # metadata
    file_metadata = {
        'Gwas2VCF_command': ' '.join(sys.argv[1:]) + "; " + version,
        'file_date': datetime.now().isoformat()
    }

    if args.cohort_controls is not None:
        sample_metadata['TotalControls'] = args.cohort_controls

    if args.cohort_cases is not None:
        sample_metadata['TotalCases'] = args.cohort_cases

    if 'ncase_col' in j or args.cohort_cases is not None:
        sample_metadata['StudyType'] = 'CaseControl'
    else:
        sample_metadata['StudyType'] = 'Continuous'

    # write to VCF
    # loop over sorted chromosome position and get record using random access
    Vcf.write_to_file(gwas, idx, args.out, fasta, j['build'], args.id, sample_metadata, file_metadata, args.csi)

    # close temp file to release disk space
    gwas.close()


if __name__ == "__main__":
    main()
