#!/usr/bin/env bash
set -euo pipefail

# get data
curl -L http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz | \
gzip -dc | \
cut -s -f2- | \
sed 's/^chr//g' | \
tr ':' '\t' | \
sed $'s/^SNP_hg19/chr\tpos/g' > example.txt

# harmonise against reference FASTA
python main.py \
--out example.bcf \
--data example.txt \
--ref /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--json example.json
