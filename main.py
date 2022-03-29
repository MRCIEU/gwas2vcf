import argparse
import json
import logging
import os
import sys
from datetime import datetime

import marshmallow
import pysam

from gwas import Gwas
from param import Param
from vcf import Vcf


def main():
    version = "1.4.3"

    parser = argparse.ArgumentParser(
        description="Map GWAS summary statistics to VCF/BCF"
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {version}"
    )
    parser.add_argument(
        "--out",
        dest="out",
        required=False,
        help="Path to output VCF/BCF. If not present then must be specified as 'out' in json file",
    )
    parser.add_argument(
        "--data",
        dest="data",
        required=False,
        help="Path to GWAS summary stats. If not present then must be specified as 'data' in json file",
    )
    parser.add_argument(
        "--ref", dest="ref", required=True, help="Path to reference FASTA"
    )
    parser.add_argument(
        "--dbsnp", dest="dbsnp", required=False, help="Path to reference dbSNP VCF"
    )
    parser.add_argument(
        "--json", dest="json", required=True, help="Path to parameters JSON"
    )
    parser.add_argument(
        "--id",
        dest="id",
        required=False,
        help="Study identifier. If not present then must be specified as 'id' in json file",
    )
    parser.add_argument(
        "--cohort_controls",
        type=int,
        dest="cohort_controls",
        required=False,
        default=None,
        help="Total study number of controls (if case/control) or total sample size if continuous. Overwrites value if present in json file.",
    )
    parser.add_argument(
        "--cohort_cases",
        type=int,
        dest="cohort_cases",
        required=False,
        default=None,
        help="Total study number of cases. Overwrites value if present in json file.",
    )
    parser.add_argument(
        "--csi",
        dest="csi",
        action="store_true",
        default=False,
        required=False,
        help="Default is to index tbi but use this flag to index csi",
    )
    parser.add_argument(
        "--log",
        dest="log",
        required=False,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level",
    )
    parser.add_argument(
        "--alias", dest="alias", required=False, help="Optional chromosome alias file"
    )

    args = parser.parse_args()

    # set logging level
    if args.log:
        logging.basicConfig(
            level=getattr(logging, args.log),
            format="%(asctime)s %(levelname)s %(message)s",
        )
    logging.info(f"Gwas2VCF {version}")
    logging.info(f"Arguments: {vars(args)}")
    logging.info("Reading JSON parameters")
    try:
        schema = Param(strict=True)
        with open(args.json, encoding='utf8') as json_fh:
            json_data = schema.load(json.load(json_fh)).data
            logging.info(f"Parameters: {json_data}")
    except json.decoder.JSONDecodeError as exception_name:
        logging.error(f"Could not read json parameter file: {exception_name}")
        sys.exit()
    except marshmallow.exceptions.ValidationError as exception_name:
        logging.error(f"Could not validate json parameter file: {exception_name}")
        sys.exit()

    logging.info("Checking input arguments")
    if args.data is None:
        if "data" in json_data.keys():
            vars(args)["data"] = json_data["data"]
        else:
            logging.error("'data' filename not provided in arguments or json file")
            sys.exit()

    if args.out is None:
        if "out" in json_data.keys():
            vars(args)["out"] = json_data["out"]
        else:
            logging.error("out filename not provided in arguments or json file")
            sys.exit()

    if args.id is None:
        if "id" in json_data.keys():
            vars(args)["id"] = json_data["id"]
        else:
            logging.error("id not provided in arguments or json file")
            sys.exit()

    if args.cohort_cases is None and "cohort_cases" in json_data.keys():
        vars(args)["cohort_cases"] = json_data["cohort_cases"]

    if args.cohort_controls is None and "cohort_controls" in json_data.keys():
        vars(args)["cohort_controls"] = json_data["cohort_controls"]

    # check values are valid
    if (args.cohort_cases is not None) and (args.cohort_cases < 1):
        logging.error("Total study number of cases must be a positive number")
        sys.exit()

    if (args.cohort_controls is not None) and (args.cohort_controls < 1):
        logging.error("Total study number of controls must be a positive number")
        sys.exit()

    if not os.path.isfile(args.data):
        logging.error(f"{args.data} file does not exist")
        sys.exit()

    if not os.path.isfile(args.ref):
        logging.error(f"{args.ref} file does not exist")
        sys.exit()

    if not os.path.exists(os.path.dirname(args.out)):
        logging.error(f"{args.out} output directory does not exist")
        sys.exit()

    dbsnp = pysam.VariantFile(args.dbsnp) if args.dbsnp is not None else None

    if args.alias is not None:
        alias = {}
        with open(args.alias, encoding='utf8') as alias_fh:
            for line in alias_fh:
                (key, val) = line.strip().split("\t")
                alias[key] = val
    else:
        alias = None

    # read in data
    # harmonise, left align and trim on-the-fly and write to pickle format
    # keep file index for each record and chromosome position to write out karyotypically sorted records later
    with pysam.FastaFile(args.ref) as fasta:
        gwas, idx, sample_metadata = Gwas.read_from_file(
            args.data,
            fasta,
            json_data["chr_col"],
            json_data["pos_col"],
            json_data["ea_col"],
            json_data["oa_col"],
            json_data["beta_col"],
            json_data["se_col"],
            json_data["pval_col"],
            json_data["delimiter"],
            json_data["header"],
            ncase_col_num=json_data.get("ncase_col"),
            rsid_col_num=json_data.get("snp_col"),
            ea_af_col_num=json_data.get("eaf_col"),
            nea_af_col_num=json_data.get("oaf_col"),
            imp_z_col_num=json_data.get("imp_z_col"),
            imp_info_col_num=json_data.get("imp_info_col"),
            ncontrol_col_num=json_data.get("ncontrol_col"),
            alias=alias,
            dbsnp=dbsnp,
        )

    if dbsnp is not None:
        dbsnp.close()

    # metadata
    file_metadata = {
        "Gwas2VCF_command": " ".join(sys.argv[1:]) + "; " + version,
        "file_date": datetime.now().isoformat(),
    }

    if args.cohort_controls is not None:
        sample_metadata["TotalControls"] = args.cohort_controls

    if args.cohort_cases is not None:
        sample_metadata["TotalCases"] = args.cohort_cases

    if "ncase_col" in json_data or args.cohort_cases is not None:
        sample_metadata["StudyType"] = "CaseControl"
    else:
        sample_metadata["StudyType"] = "Continuous"

    # write to VCF
    # loop over sorted chromosome position and get record using random access
    Vcf.write_to_file(
        gwas,
        idx,
        args.out,
        fasta,
        json_data["build"],
        args.id,
        sample_metadata,
        file_metadata,
        args.csi,
    )

    # close temp file to release disk space
    gwas.close()


if __name__ == "__main__":
    main()
