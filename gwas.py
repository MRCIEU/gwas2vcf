import logging
import gzip
import os
import tempfile
import pickle
import re
from heapq import heappush
from vgraph import norm
from pvalue_handler import PvalueHandler


class Gwas:

    def __init__(self, chrom, pos, ref, alt, b, se, nlog_pval, n, alt_freq, dbsnpid, ncase, imp_info, imp_z,
                 vcf_filter="PASS"):

        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.b = b
        self.se = se
        self.nlog_pval = nlog_pval
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
        try:
            x = fasta.fetch(
                reference=self.chrom,
                start=self.pos - 1,
                end=self.pos + len(self.ref) - 1
            ).upper()
        except:
            assert 1 == 2
        assert self.ref == x

    def normalise(self, fasta, padding=100):
        # TODO handle padding edge cases
        # skip SNVs which do not need trimming
        if len(self.ref) < 2 and len(self.alt) < 2:
            return

        # zero based indexing
        pos0 = self.pos - 1
        # get reference sequence
        seq = fasta.fetch(reference=self.chrom, start=pos0 - padding, end=pos0 + padding).upper()
        # left-align and trim alleles
        start, stop, alleles = norm.normalize_alleles(
            seq,
            padding,
            padding + len(self.ref),
            (self.ref, self.alt)
        )
        # set trimmed alleles and new position
        self.ref = alleles[0]
        self.alt = alleles[1]
        self.pos = (pos0 - padding) + start + 1

        # add start base if lost during trimming
        if len(self.ref) == 0 or len(self.alt) == 0:
            # get distance from old and new positions
            dist = (self.pos - 1) - pos0
            # extract base from seq
            left_nucleotide = seq[(padding + dist) - 1: (padding + dist)]
            # set alleles and pos
            self.ref = left_nucleotide + self.ref
            self.alt = left_nucleotide + self.alt
            self.pos = self.pos - 1

    def update_dbsnp(self, dbsnp):
        if dbsnp is None:
            raise IOError("Could not read dbsnp file")
        self.dbsnpid = None
        for rec in dbsnp.fetch(contig=self.chrom, start=self.pos - 1, stop=self.pos):
            self.dbsnpid = rec.id
            break

    def check_alleles_are_vaild(self):
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
            alias=None,
            dbsnp=None
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

        # TODO use namedtuple
        metadata = {
            'TotalVariants': 0,
            'VariantsNotRead': 0,
            'HarmonisedVariants': 0,
            'VariantsNotHarmonised': 0,
            'SwitchedAlleles': 0,
            'NormalisedVariants': 0
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

        p_value_handler = PvalueHandler ()

        for l in f:
            metadata['TotalVariants'] += 1
            s = l.strip().split(delimiter)

            logging.debug("Input row: {}".format(s))

            try:
                if alias is not None:
                    if s[chrom_field] in alias:
                        chrom = alias[s[chrom_field]]
                    else:
                        chrom = s[chrom_field]
                else:
                    chrom = s[chrom_field]
            except Exception as e:
                logging.debug("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                pos = int(float(s[pos_field]))  # float is for scientific notation
                assert pos > 0
            except Exception as e:
                logging.debug("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            ref = str(s[nea_field]).strip().upper()
            alt = str(s[ea_field]).strip().upper()

            if ref == alt:
                logging.debug("Skipping: ref={} is the same as alt={}".format(ref, alt))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                b = float(s[effect_field])
            except Exception as e:
                logging.debug("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                se = float(s[se_field])
            except Exception as e:
                logging.debug("Skipping {}: {}".format(s, e))
                metadata['VariantsNotRead'] += 1
                continue

            try:
                pval = p_value_handler.parse_string(s[pval_field])
                nlog_pval = p_value_handler.neg_log_of_decimal(pval)
            except Exception as e:
                logging.debug("Skipping line {}, {}".format(s, e))
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
                nlog_pval,
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
                result.check_alleles_are_vaild()
            except AssertionError as e:
                logging.debug("Skipping {}: {}".format(s, e))
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
                    logging.debug("Could not harmonise {}: {}".format(s, e))
                    metadata['VariantsNotHarmonised'] += 1
                    continue
            metadata['HarmonisedVariants'] += 1

            # left align and trim variants
            if len(ref) > 1 and len(alt) > 1:
                try:
                    result.normalise(fasta)
                except Exception as e:
                    logging.debug("Could not normalise {}: {}".format(s, e))
                    metadata['VariantsNotHarmonised'] += 1
                    continue
                metadata['NormalisedVariants'] += 1

            # add or update dbSNP identifier
            if dbsnp is not None:
                result.update_dbsnp(dbsnp)

            # keep file position sorted by chromosome position for recall later
            if result.chrom not in file_idx:
                file_idx[result.chrom] = []
            heappush(file_idx[result.chrom], (result.pos, results.tell()))

            try:
                pickle.dump(result, results)
            except Exception as e:
                logging.error("Could not write to {}:".format(tempfile.gettempdir()), e)
                raise e

        f.close()

        logging.info("Total variants: {}".format(metadata['TotalVariants']))
        logging.info("Variants could not be read: {}".format(metadata['VariantsNotRead']))
        logging.info("Variants harmonised: {}".format(metadata['HarmonisedVariants']))
        logging.info("Variants discarded during harmonisation: {}".format(metadata['VariantsNotHarmonised']))
        logging.info("Alleles switched: {}".format(metadata['SwitchedAlleles']))
        logging.info("Normalised variants: {}".format(metadata['NormalisedVariants']))
        logging.info("Skipped {} of {}".format(
            metadata['VariantsNotRead'] + metadata['VariantsNotHarmonised'], metadata['TotalVariants'])
        )
        if (metadata['VariantsNotRead'] + metadata['VariantsNotHarmonised']) / metadata['TotalVariants'] > .2:
            logging.warning("More than 20% of variants not read or harmonised. Check your input")

        return results, file_idx, metadata
