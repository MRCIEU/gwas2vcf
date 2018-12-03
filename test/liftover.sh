#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
set -euo pipefail
cd $PBS_O_WORKDIR

module load apps/bcftools-1.8
module load apps/gatk-4.0.8.1

# params
bcf_path="HbA1c_METAL_European.bcf"
fasta_path="/panfs/panasas01/sscm/ml18692/db/gatk/2.8/b37/human_g1k_v37.fasta"
chain_path="/panfs/panasas01/sscm/ml18692/db/gatk/chains/b36tob37.chain"

# convert to vcf
bcftools view "$bcf_path" > temp.vcf

# lift VCF to new build
gatk LiftoverVcf \
-I temp.vcf \
-O $(echo "$bcf_path" | sed 's/\.bcf/\.lifted\.bcf/g') \
-R "$fasta_path" \
--TAGS_TO_DROP B \
--TAGS_TO_DROP AF \
-C "$chain_path" \
--REJECT $(echo "$bcf_path" | sed 's/\.vcf/\.rejected\.vcf/g') \
--CREATE_INDEX
