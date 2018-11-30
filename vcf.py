import pysam
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


class Vcf:

    @staticmethod
    def write_to_file(gwas_results, path, fasta, build, params=None):
        logging.info("Writing to VCF: {}".format(path))

        header = pysam.VariantHeader()
        header.add_line(
            '##INFO=<ID=B,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele(s)">')
        header.add_line('##INFO=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">')
        header.add_line('##INFO=<ID=PVAL,Number=A,Type=Float,Description="P-value for effect estimate">')
        header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">')
        header.add_line('##INFO=<ID=N,Number=A,Type=Float,Description="Sample size used to estimate genetic effect">')

        # add contig lengths
        assert len(fasta.references) == len(fasta.lengths)
        for n, contig in enumerate(fasta.references):
            header.add_line("##contig=<ID={},length={}, assembly={}>".format(contig, fasta.lengths[n], build))

        # add metadata
        if params is not None:
            for k in params:
                header.add_line('##{}={}'.format(k, params[k]))

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
            record.info['PVAL'] = result.pval
            record.info['AF'] = result.alt_freq
            record.info['N'] = result.n
            vcf.write(record)

        vcf.close()
