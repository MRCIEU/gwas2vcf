import logging


class Harmonise:

    @staticmethod
    def align_gwas_to_fasta(variants, fasta):
        flipped_variants = 0
        harmonised = []

        for variant in variants:

            if not variant.are_alleles_iupac():
                logging.warning("Skipping record {}: allele(s) are not standard IUPAC".format(variant))
                continue

            # get expected FASTA REF
            try:
                expected_ref_allele = str(
                    fasta.fetch(region="{}:{}-{}".format(
                        variant.chrom,
                        variant.pos,
                        variant.pos + (len(variant.ref) - 1)
                    ))
                ).upper()
            except (TypeError, ValueError):
                logging.warning("Skipping record {}: problem getting ref allele".format(variant))
                continue

            if variant.ref != expected_ref_allele:
                variant.reverse_sign()
                flipped_variants += 1

                if len(variant.ref) != len(variant.alt):

                    # get expected FASTA REF
                    try:
                        expected_ref_allele = str(
                            fasta.fetch(region="{}:{}-{}".format(
                                variant.chrom,
                                variant.pos,
                                variant.pos + (len(variant.ref) - 1)
                            ))
                        ).upper()
                    except  (TypeError, ValueError):
                        logging.warning("Skipping record {}: problem getting ref allele".format(variant))
                        continue

                if variant.ref != expected_ref_allele:
                    logging.warning("Skipping record {}: could not match alleles to FASTA".format(variant))
                    continue

            harmonised.append(variant)

        return harmonised, flipped_variants
