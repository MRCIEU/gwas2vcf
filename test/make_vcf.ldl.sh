#!/usr/bin/env bash
set -euxo pipefail

# get data
#curl -L http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz | \
#gzip -dc | \
#cut -s -f2- | \
#sed 's/^chr//g' | \
#tr ':' '\t' > test/data/jointGwasMc_LDL.txt

g="data/jointGwasMc_LDL.txt"
f="data/db/human/gatk/2.8/b37/human_g1k_v37.fasta"
v="data/jointGwasMc_LDL.vcf"

# make VCF
/Users/ml/GitLab/gwas_harmonisation/venv/bin/python /Users/ml/GitLab/gwas_harmonisation/main.py \
-o "$v" \
-g "$g" \
-f "$f" \
-s 1 \
-chrom_field 0 \
-pos_field 1 \
-dbsnp_field 2 \
-ea_field 3 \
-nea_field 4 \
-effect_field 5 \
-se_field 6 \
-n0_field 7 \
-pval_field 8 \
-ea_af_field 9

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