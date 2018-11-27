from pyfaidx import Fasta
import pysam
from Bio.Seq import Seq
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


class GwasResult:

    # TODO check alleles are IUPAC
    # TODO set build from FASTA
    def __init__(self, chrom, pos, build, ref, alt, dbsnpid=None, b=None, se=None, pval=None,
                 n0=None, n1=None, alt_freq=None, vcf_filter="PASS"):
        self.chrom = str(chrom)
        self.pos = int(pos)
        self.build = str(build)
        self.ref = str(ref).upper()
        self.alt = str(alt).upper()

        try:
            self.dbsnpid = str(dbsnpid)
        except (TypeError, ValueError):
            self.dbsnpid = None

        try:
            self.b = float(b)
        except (TypeError, ValueError):
            self.b = None

        try:
            self.se = float(se)
        except (TypeError, ValueError):
            self.se = None

        try:
            self.pval = float(pval)
        except (TypeError, ValueError):
            self.pval = None

        try:
            self.n0 = float(n0)
        except (TypeError, ValueError):
            self.n0 = None

        try:
            self.n1 = float(n1)
        except (TypeError, ValueError):
            self.n1 = None

        try:
            self.alt_freq = float(alt_freq)
        except (TypeError, ValueError):
            self.alt_freq = None

        self.vcf_filter = str(vcf_filter)

    def reverse_sign(self):
        r = self.ref
        a = self.alt
        self.ref = a
        self.alt = r
        self.b = self.b * -1

    # TODO only allow once per input
    # TODO safe guard palindromics
    def flip_strand(self):
        self.ref = Seq(self.ref).reverse_complement()
        self.alt = Seq(self.alt).reverse_complement()

    # TODO add support for indels
    def is_ref_allele_match_fasta(self, path):
        with Fasta(path) as fasta:
            seq = fasta.get_seq(self.chrom, self.pos, self.pos)
        return str(seq) == self.ref

    def is_ref_allele_match_fasta_htslib(self, path):
        file = pysam.FastaFile(path)
        seq = file.fetch(region="{}:{}-{}".format(self.chrom, self.pos, self.pos))
        file.close()
        return str(seq) == self.ref

    def __str__(self):
        return f'{self.chrom}\t{self.pos}\t{self.ref}\t{self.alt}'

    @staticmethod
    def read_from_text_file(path, skip_n_rows=0):
        logging.info("Reading sumamry stats and mapping to FASTA: {}".format(path))

        results = []
        with open(path, "r") as f:
            for n, l in enumerate(f):

                if n < skip_n_rows:
                    continue

                s = l.strip().split("\t")
                result = GwasResult(
                    chrom=s[1].split(":")[0].replace("chr", ""),
                    pos=int(s[1].split(":")[1]),
                    build="b37",
                    ref=s[3],
                    alt=s[4],
                    dbsnpid=s[2],
                    b=s[5],
                    se=s[6],
                    pval=s[8],
                    n1=int(round(float(s[7]))),
                    alt_freq=s[9]
                )

                results.append(result)

        return results
