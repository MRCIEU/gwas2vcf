import pysam
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


class Vcf:

    # TODO set from fasta
    @staticmethod
    def get_contig_headers(path):
        contigs = []
        with open('test/data/contigs.vcf') as f:
            for l in f:
                contigs.append(l.strip())
            return contigs

    @staticmethod
    def write_to_file(gwas_results, path):
        logging.info("Writing to VCF: {}".format(path))

        header = pysam.VariantHeader()
        header.add_line(
            '##INFO=<ID=B,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele(s)">')
        header.add_line('##INFO=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">')
        header.add_line('##INFO=<ID=P,Number=A,Type=Float,Description="P-value for effect estimate">')
        header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">')
        header.add_line('##INFO=<ID=N1,Number=A,Type=Float,Description="Number of cases. 0 if continuous trait">')
        header.add_line(
            '##INFO=<ID=N0,Number=A,Type=Float,Description="Number of controls. Total sample size if continuous trait">')
        contigs = Vcf.get_contig_headers('test/data/contigs.vcf')
        for contig in contigs:
            header.add_line(contig)

        vcf = pysam.VariantFile(path, "w", header=header)

        for result in gwas_results:
            record = vcf.new_record()
            record.chrom = result.chrom
            record.pos = result.pos
            record.id = result.dbsnpid
            record.alleles = (result.ref, result.alt)
            record.filter.add(result.vcf_filter)
            record.info['B'] = result.b
            record.info['SE'] = result.se
            record.info['P'] = result.pval
            record.info['AF'] = result.alt_freq
            record.info['N1'] = result.n1
            record.info['N0'] = result.n0

            vcf.write(record)

        vcf.close()
