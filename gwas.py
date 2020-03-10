import logging
import gzip
import os
import tempfile
import pickle
import re
from heapq import heappush


class Gwas:

    def __init__(self, chrom, pos, ref, alt, b, se, pval, n, alt_freq, dbsnpid, ncase, imp_info, imp_z,
                 vcf_filter="PASS"):

        self.chrom = chrom
        self.pos = pos
        self.ref = str(ref).strip().upper()
        self.alt = str(alt).strip().upper()
        self.b = b
        self.se = se
        self.pval = pval
        self.alt_freq = alt_freq
        self.n = n
        self.dbsnpid = dbsnpid
        self.ncase = ncase
        self.imp_info = imp_info
        self.imp_z = imp_z
        self.vcf_filter = vcf_filter

    def reverse_sign(self):
        r = self.ref
        a = self.alt
        self.ref = a
        self.alt = r
        self.b = (self.b * -1)

        if self.imp_z is not None:
            self.imp_z = (self.imp_z * -1)

        try:
            self.alt_freq = (1 - self.alt_freq)
        except TypeError:
            self.alt_freq = None

    def check_reference_allele(self, fasta):
        assert self.ref == str(
            fasta.fetch(region="{}:{}-{}".format(
                self.chrom,
                self.pos,
                self.pos + len(self.ref) - 1
            ))
        ).upper()

    def check_alleles_iupac(self):
        for bp in self.alt:
            assert bp in {"A", "T", "G", "C"}
        for bp in self.ref:
            assert bp in {"A", "T", "G", "C"}

    def __str__(self):
        return str(self.__dict__)

    """ Function to read in GWAS data from plain text file """

    @staticmethod
    def read_from_file(
            path,
            fasta,
            chrom_field,
            pos_field,
            ea_field,
            nea_field,
            effect_field,
            se_field,
            pval_field,
            delimiter,
            header,
            ncase_field=None,
            rsid_field=None,
            ea_af_field=None,
            nea_af_field=None,
            imp_z_field=None,
            imp_info_field=None,
            ncontrol_field=None,
            rm_chr_prefix=False
    ):

        rsid_pattern = re.compile("^rs[0-9]*$")

        logging.info("Reading summary stats and mapping to FASTA: {}".format(path))
        logging.debug("File path: {}".format(path))
        logging.debug("CHR field: {}".format(chrom_field))
        logging.debug("POS field: {}".format(pos_field))
        logging.debug("EA field: {}".format(ea_field))
        logging.debug("NEA field: {}".format(nea_field))
        logging.debug("Effect field: {}".format(effect_field))
        logging.debug("SE field: {}".format(se_field))
        logging.debug("P fields: {}".format(pval_field))
        logging.debug("Delimiter: {}".format(delimiter))
        logging.debug("Header: {}".format(header))
        logging.debug("ncase Field: {}".format(ncase_field))
        logging.debug("dbsnp Field: {}".format(rsid_field))
        logging.debug("EA AF Field: {}".format(ea_af_field))
        logging.debug("NEA AF Field: {}".format(nea_af_field))
        logging.debug("IMP Z Score Field: {}".format(imp_z_field))
        logging.debug("IMP INFO Field: {}".format(imp_info_field))
        logging.debug("N Control Field: {}".format(ncontrol_field))

        metadata = {
            'TotalVariants': 0,
            'VariantsNotRead': 0,
            'HarmonisedVariants': 0,
            'VariantsNotHarmonised': 0,
            'SwitchedAlleles': 0
        }
        file_idx = dict()
        filename, file_extension = os.path.splitext(path)

        if file_extension == '.gz':
            logging.info("Reading gzip file")
            f = gzip.open(path, 'rt')
        else:
            logging.info("Reading plain text file")
            f = open(path, 'r')

        # skip header line (if present)
        if header:
            logging.info("Skipping header: {}".format(f.readline().strip()))

        # store results in a serialised temp file to reduce memory usage
        results = tempfile.TemporaryFile()

        for l in f:
            metadata['TotalVariants'] += 1
            s = l.strip().split(delimiter)

            logging.debug("Input row: {}".format(s))

            try:
                if rm_chr_prefix:
                    chrom = s[chrom_field].replace("chr", "")
                else:
                    chrom = s[chrom_field]
            except Exception as e:
                logging.warning("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                pos = int(float(s[pos_field]))  # float is for scientific notation
                assert pos > 0
            except Exception as e:
                logging.warning("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            ref = s[nea_field].upper()
            alt = s[ea_field].upper()

            try:
                b = float(s[effect_field])
            except Exception as e:
                logging.warning("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                se = float(s[se_field])
            except Exception as e:
                logging.warning("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                pval = float(s[pval_field])
            except Exception as e:
                logging.warning("Skipping line {}, {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                if ea_af_field is not None:
                    alt_freq = float(s[ea_af_field])
                elif nea_af_field is not None:
                    alt_freq = 1 - float(s[nea_af_field])
                else:
                    alt_freq = None
            except (IndexError, TypeError, ValueError) as e:
                logging.debug("Could not parse allele frequency: {}".format(e))
                alt_freq = None

            try:
                rsid = s[rsid_field]
                assert rsid_pattern.match(rsid)
            except (IndexError, TypeError, ValueError, AssertionError) as e:
                logging.debug("Could not parse dbsnp identifier: {}".format(e))
                rsid = None

            try:
                ncase = float(s[ncase_field])
            except (IndexError, TypeError, ValueError) as e:
                logging.debug("Could not parse number of cases: {}".format(e))
                ncase = None

            try:
                ncontrol = float(s[ncontrol_field])
            except (IndexError, TypeError, ValueError) as e:
                logging.debug("Could not parse number of controls: {}".format(e))
                ncontrol = None

            try:
                n = ncase + ncontrol
            except (IndexError, TypeError, ValueError) as e:
                logging.debug("Could not sum cases and controls: {}".format(e))
                n = ncontrol

            try:
                imp_info = float(s[imp_info_field])
            except (IndexError, TypeError, ValueError) as e:
                logging.debug("Could not parse imputation INFO: {}".format(e))
                imp_info = None

            try:
                imp_z = float(s[imp_z_field])
            except (IndexError, TypeError, ValueError) as e:
                logging.debug("Could not parse imputation Z score: {}".format(e))
                imp_z = None

            result = Gwas(
                chrom,
                pos,
                ref,
                alt,
                b,
                se,
                pval,
                n,
                alt_freq,
                rsid,
                ncase,
                imp_info,
                imp_z
            )

            logging.debug("Extracted row: {}".format(result))

            # check alleles
            try:
                result.check_alleles_iupac()
            except AssertionError as e:
                logging.warning("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            # harmonise alleles
            try:
                result.check_reference_allele(fasta)
            except AssertionError:
                try:
                    result.reverse_sign()
                    result.check_reference_allele(fasta)
                    metadata['SwitchedAlleles'] += 1
                except AssertionError as e:
                    logging.warning("Could not harmonise {}: {}".format(s, e))
                    metadata['VariantsNotHarmonised'] += 1
                    continue
            metadata['HarmonisedVariants'] += 1

            # keep file position for sorted recall later
            if result.chrom not in file_idx:
                file_idx[result.chrom] = []
            heappush(file_idx[result.chrom], (result.pos, results.tell()))
            pickle.dump(result, results)

        f.close()

        logging.info("Total variants: {}".format(metadata['TotalVariants']))
        logging.info("Variants could not be read: {}".format(metadata['VariantsNotRead']))
        logging.info("Variants harmonised: {}".format(metadata['HarmonisedVariants']))
        logging.info("Variants discarded during harmonisation: {}".format(metadata['VariantsNotHarmonised']))
        logging.info("Alleles switched: {}".format(metadata['SwitchedAlleles']))
        logging.info("Skipped {} of {}".format(
            metadata['VariantsNotRead'] + metadata['VariantsNotHarmonised'], metadata['TotalVariants'])
        )

        return results, file_idx, metadata
