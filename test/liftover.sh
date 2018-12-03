#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
set -euo pipefail
cd $PBS_O_WORKDIR

module load apps/bcftools-1.8
module load apps/gatk-4.0.8.1

# params
bcf_path=""
target_fasta_path=""
chain_path=""

# convert to vcf
bcftools view "$bcf_path" > temp.vcf

# lift VCF to new build
gatk LiftoverVcf \
-I temp.vcf \
-O $(echo "$bcf_path" | sed 's/\.bcf/\.lifted\.vcf/g') \
-R "$target_fasta_path" \
--TAGS_TO_DROP B \
--TAGS_TO_DROP AF \
-C "$chain_path" \
--REJECT $(echo "$bcf_path" | sed 's/\.vcf/\.rejected\.vcf/g')

# convet to bcf
bcftools view $(echo "$bcf_path" | sed 's/\.bcf/\.lifted\.vcf/g') -Ob -o $(echo "$bcf_path" | sed 's/\.bcf/\.lifted\.bcf/g')
bcftools index $(echo "$bcf_path" | sed 's/\.bcf/\.lifted\.bcf/g')

rm temp.vcf
