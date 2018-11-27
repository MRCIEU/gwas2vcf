#!/usr/bin/env bash
set -euo pipefail

g="test/data/jointGwasMc_LDL.split.n50k.txt"
f="/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta"
v="test/data/jointGwasMc_LDL.vcf"

# make VCF
/Users/ml/GitLab/gwas_harmonisation/venv/bin/python /Users/ml/GitLab/gwas_harmonisation/main.py \
-o "$v" \
-g "$g" \
-f "$f" \
-s 1 \
-chrom_field 2 \
-pos_field 3 \
-a1_field 5 \
-a2_field 6 \
-effect_field 7 \
-se_field 8 \
-n1_field 9 \
-pval_field 10 \
-dbsnp_field 4 \
-a2_af_field 11

# sort vcf
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools sort \
-i "$v" \
-faidx /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai \
-header > $(echo "$v" | sed 's/.vcf/.sorted.vcf/g')

# validate vcf
java -Xmx2g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V $(echo "$v" | sed 's/.vcf/.sorted.vcf/g')

# combine multi allelics
# TODO

# convert to bcf
# TODO