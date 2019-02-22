#!/usr/bin/env bash
set -euxo pipefail

# get data
curl -L https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz | \
gzip -dc > data/Meta-analysis_Locke_et_al_UKBiobank_2018_UPDATED.txt

g="data/Meta-analysis_Locke_et_al_UKBiobank_2018_UPDATED.txt"
f="/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta"
v="data/Meta-analysis_Locke_et_al_UKBiobank_2018_UPDATED.vcf"

# make VCF
/Users/ml/GitLab/gwas_harmonisation/venv/bin/python /Users/ml/GitLab/gwas_harmonisation/main.py \
-o "$v" \
-g "$g" \
-f "$f" \
-b "b37" \
-s 1 \
-chrom_field 0 \
-pos_field 1 \
-dbsnp_field 2 \
-ea_field 3 \
-nea_field 4 \
-ea_af_field 5 \
-effect_field 6 \
-se_field 7 \
-pval_field 8 \
-n_field 9

# sort vcf
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools sort \
-i "$v" \
-faidx "$f".fai \
-header > $(echo "$v" | sed 's/.vcf/.sorted.vcf/g')

# validate vcf
java -Xmx2g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R "$f" \
-V $(echo "$v" | sed 's/.vcf/.sorted.vcf/g')

# convertto bcf
bcftools view \
-Ob \
-o $(echo "$vcf_path" | sed 's/.vcf/.bcf/g') \
$(echo "$vcf_path" | sed 's/.vcf/.1kg.vcf/g')

# index bcf
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools index $(echo "$v" | sed 's/.vcf/.bcf/g')