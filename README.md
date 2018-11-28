# GWAS Harmonisation

Tool to map GWAS summary statistics from tab-delim to VCF with on the fly harmonisation to a supplied reference FASTA.

### Installing

```
git clone https://ieugit-scmv-d0.epi.bris.ac.uk/ml18692/gwas_harmonisation.git
cd gwas_harmonisation
pip install -r ./requirements.txt
```

## Running the tests

```
python -m unittest discover gwas_harmonisation/test
```

## Usage

```
python gwas_harmonisation/main.py \
-o <output_vcf> \
-g <input_gwas> \
-f <fasta_path> \
-s <rows_to_skip> \
-chrom_field <n> \
-pos_field <n> \
-dbsnp_field <n> \
-a1_field <n> \
-a2_field <n> \
-a1_af_field <n> \
-effect_field <n> \
-se_field <n> \
-pval_field <n> \
-n0_field <n>
```

## Full workflow

```bash make_vcf.ldl.sh```