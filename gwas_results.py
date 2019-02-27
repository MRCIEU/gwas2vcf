import logging
import gzip
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


class GwasResult:

    def __init__(self, chrom, pos, ref, alt, b, se, pval, n, alt_freq, dbsnpid, prop_cases, imp_info, imp_z,
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
        self.prop_cases = prop_cases
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

    def are_alleles_iupac(self):
        for bp in self.alt:
            if bp != 'A' and bp != 'T' and bp != 'C' and bp != 'G':
                return False

        for bp in self.ref:
            if bp != 'A' and bp != 'T' and bp != 'C' and bp != 'G':
                return False

        return True

    def __str__(self):
        return str(self.__dict__)

    @staticmethod
    def read_from_text_file(
            path,
            chrom_field,
            pos_field,
            ea_field,
            nea_field,
            effect_field,
            se_field,
            pval_field,
            dbsnp_field=None,
            n_field=None,
            ea_af_field=None,
            nea_af_field=None,
            imp_z_field=None,
            imp_info_field=None,
            prop_cases_field=None,
            skip_n_rows=0):

        logging.info("Reading summary stats and mapping to FASTA: {}".format(path))
        total_variants = 0
        filename, file_extension = os.path.splitext(path)

        if file_extension == '.gz':
            logging.info("Reading gzip file")
            f = gzip.open(path, 'rt')
        else:
            logging.info("Reading plain text file")
            f = open(path, 'r')

        results = []
        for n, l in enumerate(f):

            if n < skip_n_rows:
                logging.info("Skipping header: {}".format(l.strip()))
                continue

            total_variants += 1

            s = l.strip().split("\t")

            chrom = s[chrom_field]

            try:
                pos = int(float(s[pos_field]))  # float is for scientific notation
            except ValueError as e:
                logging.warning("Skipping {}: {}".format(s, e))
                continue

            ref = s[nea_field]
            alt = s[ea_field]
            b = float(s[effect_field])
            se = float(s[se_field])
            pval = float(s[pval_field])

            try:
                if ea_af_field is not None:
                    alt_freq = float(s[ea_af_field])
                elif nea_af_field is not None:
                    alt_freq = 1 - float(s[nea_af_field])
                else:
                    alt_freq = None
            except (IndexError, TypeError, ValueError):
                alt_freq = None

            try:
                dbsnpid = s[dbsnp_field]
            except (IndexError, TypeError, ValueError):
                dbsnpid = None

            try:
                n = float(s[n_field])
            except (IndexError, TypeError, ValueError):
                n = None

            try:
                prop_cases = float(s[prop_cases_field])
            except (IndexError, TypeError, ValueError):
                prop_cases = None

            try:
                imp_info = float(s[imp_info_field])
            except (IndexError, TypeError, ValueError):
                imp_info = None

            try:
                imp_z = float(s[imp_z_field])
            except (IndexError, TypeError, ValueError):
                imp_z = None

            result = GwasResult(
                chrom,
                pos,
                ref,
                alt,
                b,
                se,
                pval,
                n,
                alt_freq,
                dbsnpid,
                prop_cases,
                imp_info,
                imp_z
            )

            results.append(result)

        f.close()

        return results, total_variants
