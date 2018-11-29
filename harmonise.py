import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')


class Harmonise:

    @staticmethod
    def align_gwas_to_fasta(variants, fasta):
        excluded_variants = 0
        flipped_variants = 0
        harmonised = []

        for variant in variants:

            if not variant.are_alleles_iupac():
                logging.warning("Skipping record {}: allele(s) are not standard IUPAC".format(variant))
                excluded_variants += 1
                continue

            # get expected FASTA REF
            expected_ref_allele = str(
                fasta.fetch(region="{}:{}-{}".format(
                    variant.chrom,
                    variant.pos,
                    variant.pos + (len(variant.ref) - 1)
                ))
            ).upper()

            if variant.ref != expected_ref_allele:
                variant.reverse_sign()
                flipped_variants += 1

                if len(variant.ref) != len(variant.alt):
                    # get expected FASTA REF
                    expected_ref_allele = str(
                        fasta.fetch(region="{}:{}-{}".format(
                            variant.chrom,
                            variant.pos,
                            variant.pos + (len(variant.ref) - 1)
                        ))
                    ).upper()

                if variant.ref != expected_ref_allele:
                    logging.warning("Skipping record {}: could not match alleles to FASTA".format(variant))
                    excluded_variants += 1
                    continue

            harmonised.append(variant)

        return harmonised, excluded_variants, flipped_variants
