from Bio.Seq import Seq
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


class GwasResult:

    # TODO set build from FASTA
    def __init__(self, chrom, pos, build, ref, alt, b, se, pval, n0, alt_freq, n1=None, dbsnpid=None,
                 vcf_filter="PASS"):

        self.chrom = chrom
        self.pos = pos
        self.build = build
        self.ref = str(ref).strip().upper()
        self.alt = str(alt).strip().upper()
        self.b = b
        self.se = se
        self.pval = pval
        self.n0 = n0
        self.alt_freq = alt_freq
        self.n1 = n1
        self.dbsnpid = dbsnpid
        self.vcf_filter = vcf_filter

    def reverse_sign(self):
        r = self.ref
        a = self.alt
        self.ref = a
        self.alt = r
        self.b = self.b * -1
        try:
            self.alt_freq = (1 - self.alt_freq)
        except TypeError:
            self.alt_freq = None

    # TODO only allow once per input
    # TODO safe guard palindromics
    def flip_strand(self):
        self.ref = Seq(self.ref).reverse_complement()
        self.alt = Seq(self.alt).reverse_complement()

    def are_alleles_iupac(self):
        for bp in self.alt:
            if bp != 'A' and bp != 'T' and bp != 'C' and bp != 'G':
                return False

        for bp in self.ref:
            if bp != 'A' and bp != 'T' and bp != 'C' and bp != 'G':
                return False

        return True

    def __str__(self):
        return f'{self.chrom}\t{self.pos}\t{self.ref}\t{self.alt}'

    @staticmethod
    def read_from_text_file(
            path,
            chrom_field,
            pos_field,
            a1_field,
            a2_field,
            effect_field,
            se_field,
            pval_field,
            n0_field,
            dbsnp_field=None,
            n1_field=None,
            a1_af_field=None,
            a2_af_field=None,
            skip_n_rows=0):

        logging.info("Reading summary stats and mapping to FASTA: {}".format(path))

        results = []
        with open(path, "r") as f:
            for n, l in enumerate(f):

                if n < skip_n_rows:
                    continue

                s = l.strip().split("\t")

                chrom = s[chrom_field]
                pos = int(s[pos_field])
                ref = s[a1_field]
                alt = s[a2_field]
                b = float(s[effect_field])
                se = float(s[se_field])
                pval = float(s[pval_field])
                n0 = float(s[n0_field])

                try:
                    if a2_af_field is not None:
                        alt_freq = float(s[a2_af_field])
                    elif a1_af_field is not None:
                        alt_freq = 1 - float(s[a1_af_field])
                    else:
                        alt_freq = None
                except (IndexError, TypeError, ValueError):
                    alt_freq = None

                try:
                    dbsnpid = s[dbsnp_field]
                except (IndexError, TypeError, ValueError):
                    dbsnpid = None

                try:
                    n1 = float(s[n1_field])
                except (IndexError, TypeError, ValueError):
                    n1 = None

                # TODO add builds
                result = GwasResult(
                    chrom,
                    pos,
                    "b37",
                    ref,
                    alt,
                    b,
                    se,
                    pval,
                    n0,
                    alt_freq,
                    n1=n1,
                    dbsnpid=dbsnpid
                )

                results.append(result)

        return results
