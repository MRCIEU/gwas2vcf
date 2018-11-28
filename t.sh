#!/usr/bin/env bash
set -euxo pipefail

g="test/data/jointGwasMc_LDL.txt"
f="/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta"
v="test/data/jointGwasMc_LDL.vcf"

# combine multi allelics
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools norm \
--check-ref e \
-f "$f" \
-m +any \
-Ob \
-o $(echo "$v" | sed 's/.vcf/.bcf/g') \
$(echo "$v" | sed 's/.vcf/.sorted.vcf/g')
