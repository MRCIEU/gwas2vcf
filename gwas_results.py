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
        self.alt_freq = (1 - self.alt_freq)

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
    def read_from_text_file(
            path,
            effect_field,
            se_field,
            pval_field,
            dbsnp_field=None,
            chrom_field=None,
            pos_field=None,
            a1_field=None,
            a2_field=None,
            n0_field=None,
            n1_field=None,
            a2_af_field=None,
            skip_n_rows=0):

        logging.info("Reading sumamry stats and mapping to FASTA: {}".format(path))

        results = []
        with open(path, "r") as f:
            for n, l in enumerate(f):

                if n < skip_n_rows:
                    continue

                s = l.strip().split("\t")
                print(s)

                try:
                    chrom = s[int(chrom_field)]
                except IndexError:
                    chrom = None

                try:
                    pos = int(s[int(pos_field)])
                except Exception:
                    pos = None

                try:
                    ref = s[int(a1_field)]
                except Exception:
                    ref = None

                try:
                    alt = s[int(a2_field)]
                except Exception:
                    alt = None

                try:
                    dbsnpid = s[int(dbsnp_field)]
                except Exception:
                    dbsnpid = None

                try:
                    b = s[int(effect_field)]
                except Exception:
                    b = None

                try:
                    se = s[int(se_field)]
                except Exception:
                    se = None

                try:
                    pval = s[int(pval_field)]
                except Exception:
                    pval = None

                try:
                    n1 = float(s[int(n1_field)])
                except Exception:
                    n1 = None

                try:
                    n0 = float(s[int(n0_field)])
                except Exception:
                    n0 = None

                try:
                    alt_freq = float(s[int(a2_af_field)])
                except Exception:
                    alt_freq = None

                # TODO add builds
                result = GwasResult(
                    chrom=chrom,
                    pos=pos,
                    build="b37",
                    ref=ref,
                    alt=alt,
                    dbsnpid=dbsnpid,
                    b=b,
                    se=se,
                    pval=pval,
                    n1=n1,
                    n0=n0,
                    alt_freq=alt_freq
                )

                results.append(result)

        return results
