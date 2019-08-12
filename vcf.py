import pysam
import logging
import numpy as np


class Vcf:

    @staticmethod
    def convert_pval_to_neg_log10(p):
        # prevent negative 0 output
        if p == 1:
            return 0
        # prevent Inf output
        if p == 0:
            return 999
        return -np.log10(p)

    @staticmethod
    def write_to_file(gwas_results, path, fasta, build, trait_id, params=None):
        logging.info("Writing headers to BCF/VCF: {}".format(path))

        header = pysam.VariantHeader()
        header.add_line(
            '##FORMAT=<ID=EFFECT,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">')
        header.add_line('##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">')
        header.add_line('##FORMAT=<ID=L10PVAL,Number=A,Type=Float,Description="-log10 p-value for effect estimate">')
        header.add_line(
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in the association study">')
        header.add_line('##FORMAT=<ID=N,Number=A,Type=Float,Description="Sample size used to estimate genetic effect">')
        header.add_line(
            '##FORMAT=<ID=ZSCORE,Number=A,Type=Float,Description="Z-score provided if it was used to derive the EFFECT and SE fields">')
        header.add_line(
            '##FORMAT=<ID=SIMPINFO,Number=A,Type=Float,Description="Accuracy score of summary data imputation">')
        header.add_line(
            '##FORMAT=<ID=PROPCASES,Number=A,Type=Float,Description="Proportion of sample size that were cases in the GWAS">')
        header.samples.add(trait_id)

        # add contig lengths
        assert len(fasta.references) == len(fasta.lengths)
        for n, contig in enumerate(fasta.references):
            header.add_line("##contig=<ID={},length={}, assembly={}>".format(contig, fasta.lengths[n], build))

        # add metadata
        if params is not None:
            for k in params:
                header.add_line('##{}={}'.format(k, params[k]))

        vcf = pysam.VariantFile(path, "w", header=header)

        records = dict()
        for result in gwas_results:
            lpval = Vcf.convert_pval_to_neg_log10(result.pval)

            # check floats
            if result.b is not None and result.b > 0 and result.b < 1e-06:
                logging.warning(
                    "Effect field smaller than VCF specification. Expect loss of precision for: {}".format(
                        result.b)
                )
            if result.se is not None and result.se > 0 and result.se < 1e-06:
                # set SE to lowest possible value
                result.se = 1e-06
                logging.warning(
                    "Standard error field smaller than VCF specification. Expect loss of precision for: {}".format(
                        result.se)
                )
            if lpval is not None and lpval > 0 and lpval < 1e-06:
                logging.warning(
                    "-log10(pval) field smaller than VCF specification. Expect loss of precision for: {}".format(
                        lpval)
                )
            if result.alt_freq is not None and result.alt_freq > 0 and result.alt_freq < 1e-06:
                logging.warning(
                    "Allele frequency field smaller than VCF specification. Expect loss of precision for: {}".format(
                        result.alt_freq)
                )
            if result.n is not None and result.n > 0 and result.n < 1e-06:
                logging.warning(
                    "Sample size field smaller than VCF specification. Expect loss of precision for: {}".format(
                        result.n)
                )
            if result.imp_z is not None and result.imp_z > 0 and result.imp_z < 1e-06:
                logging.warning(
                    "Imputation Z score field smaller than VCF specification. Expect loss of precision for: {}".format(
                        result.imp_z)
                )
            if result.imp_info is not None and result.imp_info > 0 and result.imp_info < 1e-06:
                logging.warning(
                    "Imputation INFO field smaller than VCF specification. Expect loss of precision for: {}".format(
                        result.imp_info)
                )
            if result.prop_cases is not None and result.prop_cases > 0 and result.prop_cases < 1e-06:
                logging.warning(
                    "Proportion of cases field smaller than VCF specification. Expect loss of precision for: {}".format(
                        result.prop_cases)
                )

            record = vcf.new_record()
            record.chrom = result.chrom
            record.pos = result.pos
            record.id = result.dbsnpid
            record.alleles = (result.ref, result.alt)
            record.filter.add(result.vcf_filter)

            record.samples[trait_id]['EFFECT'] = result.b
            record.samples[trait_id]['SE'] = result.se
            record.samples[trait_id]['L10PVAL'] = lpval
            record.samples[trait_id]['AF'] = result.alt_freq
            record.samples[trait_id]['N'] = result.n
            record.samples[trait_id]['ZSCORE'] = result.imp_z
            record.samples[trait_id]['SIMPINFO'] = result.imp_info
            record.samples[trait_id]['PROPCASES'] = result.prop_cases

            # bank variants by chromosome
            if result.chrom not in records:
                records[result.chrom] = []

            records[result.chrom].append(record)

        # sort records & write records to VCF/BCF
        logging.info("Sorting records by chromosome and position")
        for contig in fasta.references:
            if contig in records:
                records[contig].sort(key=lambda x: x.pos)
                for record in records[contig]:
                    vcf.write(record)

        vcf.close()

        # index output file
        logging.info("Indexing output file")
        pysam.tabix_index(path, preset="vcf", force=True, csi=True)
