import pysam
import logging
import numpy as np
import pickle
from heapq import heappop

class Vcf:

    # no longer necessary
    # @staticmethod
    # def convert_pval_to_neg_log10(p):
    #     # prevent negative 0 output
    #     if p == 1:
    #         return 0
    #     # prevent Inf output
    #     if p == 0:
    #         return 999
    #     return -np.log10(p)

    @staticmethod
    def is_float32_lossy(f):
        if f == 0 or f == -0 or f is None or f == np.inf or f == -np.inf:
            return False

        # convert val to float32
        v = np.float32(f)

        return v == 0 or v == np.inf or v == -0 or v == -np.inf

    @staticmethod
    def remove_illegal_chars(s):
        if s is None:
            return None
        r = s.strip()
        r = r.replace(" ", "_")
        r = r.replace(";", "_")
        r = r.replace(":", "_")
        r = r.replace("=", "_")
        r = r.replace(",", "_")
        return r

    """ 
    Write GWAS file to VCF
    Expects an open file handle to a Pickle file of GWAS results & file index dict(chromosome[(position, offset)]) 
    """

    @staticmethod
    def write_to_file(gwas_file, gwas_idx, path, fasta, build, trait_id, sample_metadata=None, file_metadata=None,
                      csi=False):
        logging.info("Writing headers to BCF/VCF: {}".format(path))

        header = pysam.VariantHeader()

        # INFO
        header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')

        # FORMAT
        header.add_line(
            '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">')
        header.add_line('##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">')
        header.add_line('##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">')
        header.add_line(
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in the association study">')
        header.add_line(
            '##FORMAT=<ID=SS,Number=A,Type=Integer,Description="Sample size used to estimate genetic effect">')
        header.add_line(
            '##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score provided if it was used to derive the EFFECT and SE fields">')
        header.add_line('##FORMAT=<ID=SI,Number=A,Type=Float,Description="Accuracy score of summary data imputation">')
        header.add_line(
            '##FORMAT=<ID=NC,Number=A,Type=Integer,Description="Number of cases used to estimate genetic effect">')
        header.add_line(
            '##FORMAT=<ID=ID,Number=1,Type=String,Description="Study variant identifier">')

        # META
        header.add_line(
            '##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">')
        header.add_line(
            '##META=<ID=VariantsNotRead,Number=1,Type=Integer,Description="Number of variants that could not be read">')
        header.add_line(
            '##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">')
        header.add_line(
            '##META=<ID=VariantsNotHarmonised,Number=1,Type=Integer,Description="Total number of variants that could not be harmonised">')
        header.add_line(
            '##META=<ID=SwitchedAlleles,Number=1,Type=Integer,Description="Total number of variants strand switched">')
        header.add_line(
            '##META=<ID=TotalControls,Number=1,Type=Integer,Description="Total number of controls in the association study">')
        header.add_line(
            '##META=<ID=TotalCases,Number=1,Type=Integer,Description="Total number of cases in the association study">')
        header.add_line(
            '##META=<ID=StudyType,Number=1,Type=String,Description="Type of GWAS study [Continuous or CaseControl]">')

        # SAMPLES
        header.samples.add(trait_id)
        if file_metadata is not None:
            s = ""
            for k in sample_metadata:
                s += ",{}={}".format(k, sample_metadata[k])
            header.add_line('##SAMPLE=<ID={}{}>'.format(trait_id, s))

        # CONTIG
        assert len(fasta.references) == len(fasta.lengths)
        for n, contig in enumerate(fasta.references):
            header.add_line("##contig=<ID={},length={}, assembly={}>".format(contig, fasta.lengths[n], build))

        # add metadata
        if file_metadata is not None:
            for k in file_metadata:
                header.add_line('##{}={}'.format(k, file_metadata[k]))

        vcf = pysam.VariantFile(path, "w", header=header)

        # recall variant objects in chromosome position order
        logging.info("Writing variants to BCF/VCF: {}".format(path))
        for contig in fasta.references:
            if contig not in gwas_idx:
                continue

            while gwas_idx[contig]:
                chr_pos = heappop(gwas_idx[contig])

                # load GWAS result
                gwas_file.seek(chr_pos[1])
                result = pickle.load(gwas_file)

                result.nlog_pval = result.nlog_pval

                # check floats
                if Vcf.is_float32_lossy(result.b):
                    logging.warning(
                        "Effect field cannot fit into float32. Expect loss of precision for: {}".format(
                            result.b)
                    )
                if Vcf.is_float32_lossy(result.se):
                    result.se = np.float64(np.finfo(np.float32).tiny).item()
                    logging.warning(
                        "Standard error field cannot fit into float32. Expect loss of precision for: {}".format(
                            result.se)
                    )
                if Vcf.is_float32_lossy(result.nlog_pval):
                    logging.warning(
                        "-log10(pval) field cannot fit into float32. Expect loss of precision for: {}".format(
                            result.nlog_pval)
                    )
                if Vcf.is_float32_lossy(result.alt_freq):
                    logging.warning(
                        "Allele frequency field cannot fit into float32. Expect loss of precision for: {}".format(
                            result.alt_freq)
                    )
                if Vcf.is_float32_lossy(result.imp_z):
                    logging.warning(
                        "Imputation Z score field cannot fit into float32. Expect loss of precision for: {}".format(
                            result.imp_z)
                    )
                if Vcf.is_float32_lossy(result.imp_info):
                    logging.warning(
                        "Imputation INFO field cannot fit into float32. Expect loss of precision for: {}".format(
                            result.imp_info)
                    )

                record = vcf.new_record()
                record.chrom = result.chrom
                assert " " not in record.chrom
                record.pos = result.pos
                assert record.pos > 0
                record.id = Vcf.remove_illegal_chars(result.dbsnpid)
                record.alleles = (result.ref, result.alt)
                record.filter.add(result.vcf_filter)

                if result.alt_freq is not None:
                    record.info['AF'] = result.alt_freq

                if result.b is not None:
                    record.samples[trait_id]['ES'] = result.b
                if result.se is not None:
                    record.samples[trait_id]['SE'] = result.se
                if result.nlog_pval is not None:
                    record.samples[trait_id]['LP'] = result.nlog_pval
                if result.alt_freq is not None:
                    record.samples[trait_id]['AF'] = result.alt_freq
                if result.n is not None:
                    record.samples[trait_id]['SS'] = round(result.n)
                if result.imp_z is not None:
                    record.samples[trait_id]['EZ'] = result.imp_z
                if result.imp_info is not None:
                    record.samples[trait_id]['SI'] = result.imp_info
                if result.ncase is not None:
                    record.samples[trait_id]['NC'] = round(result.ncase)
                if result.dbsnpid is not None:
                    record.samples[trait_id]['ID'] = record.id

                # write to file
                vcf.write(record)

        vcf.close()

        # index output file
        logging.info("Indexing output file")
        pysam.tabix_index(path, preset="vcf", force=True, csi=csi)
