#!/usr/bin/env bash
set -euxo pipefail

# get data
#curl -L http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz | \
#gzip -dc > data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt

# fix NL
#dos2unix data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt

g="data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt"
f="/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta"
v="data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.vcf"

# make VCF
/Users/ml/GitLab/gwas_harmonisation/venv/bin/python /Users/ml/GitLab/gwas_harmonisation/main.py \
-o "$v" \
-g "$g" \
-f "$f" \
-b "b37" \
-s 1 \
-dbsnp_field 1 \
-chrom_field 2 \
-pos_field 3 \
-ea_field 4 \
-nea_field 5 \
-ea_af_field 6 \
-effect_field 7 \
-se_field 8 \
-pval_field 9

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

# combine multi allelics & output bcf
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools norm \
--check-ref e \
-f "$f" \
-m +any \
-Ob \
-o $(echo "$v" | sed 's/.vcf/.bcf/g') \
$(echo "$v" | sed 's/.vcf/.sorted.vcf/g')

# index bcf
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools index $(echo "$v" | sed 's/.vcf/.bcf/g')