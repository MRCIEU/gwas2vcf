#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
set -euo pipefail
cd $PBS_O_WORKDIR

module load languages/python-anaconda-5.0.1-3.6.3-phys015202
module load apps/bedtools2
module load apps/gatk-4.0.8.1
module load apps/bcftools-1.8

# params
gwas_path=""
vcf_path=""
fasta_path="/panfs/panasas01/sscm/ml18692/db/gatk/2.8/b37/human_g1k_v37.fasta"

# make VCF
python ~/apps/gwas_harmonisation/main.py \
-g "$gwas_path" \
-o "$vcf_path" \
-f "$fasta_path" \
-s  \
-chrom_field  \
-pos_field  \
-dbsnp_field  \
-ea_field  \
-nea_field  \
-ea_af_field  \
-effect_field  \
-se_field  \
-pval_field  \
-n_field 

# sort vcf
bedtools sort \
-i "$vcf_path" \
-faidx "$fasta_path".fai \
-header > $(echo "$vcf_path" | sed 's/.vcf/.sorted.vcf/g')

# validate vcf
gatk \
-T ValidateVariants \
-R "$fasta_path" \
-V $(echo "$vcf_path" | sed 's/.vcf/.sorted.vcf/g')

# combine multi allelics & output bcf
bcftools norm \
--check-ref e \
-f "$fasta_path" \
-m +any \
-Ob \
-o $(echo "$vcf_path" | sed 's/.vcf/.bcf/g') \
$(echo "$vcf_path" | sed 's/.vcf/.sorted.vcf/g')

# index bcf
bcftools index \
$(echo "$vcf_path" | sed 's/.vcf/.bcf/g')

