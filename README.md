# GWAS Harmonisation

Tool to map GWAS summary statistics from tab-delim to VCF with on the fly harmonisation to a supplied reference FASTA.

### Bug fixes ##

- 8c9a75b1 fixed an error where reading in GWAS files EA and NEA were incorrect leading wrong betas

### Installing

```
git clone https://ieugit-scmv-d0.epi.bris.ac.uk/ml18692/gwas_harmonisation.git
cd gwas_harmonisation
pip install --user -r ./requirements.txt
```

### Reference FASTA

```
# GRCh37/hg19/b37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai

# GRCh38/hg38/b38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
```

## Running the tests

Unit tests:

```
cd gwas_harmonisation
python -m unittest discover test
```

End-to-end tests:

```
# convert GIANT to BCF
bash test/make_vcf.bmi.sh

# convert GLGC to BCF
bash test/make_vcf.ldl.sh

# convert CVD to BCF
bash test/make_vcf.cvd.sh

# convert T2DM to BCF
bash test/make_vcf.t2dm.sh
```

## Usage

Column field numbers are 0-based

```
python gwas_harmonisation/main.py -h
  -h, --help                    show this help message and exit
  -v, --version                 show program's version number and exit
  -o OUT                        Path to output VCF
  -g GWAS                       Path to GWAS tab-sep summary stats
  -f FASTA                      Path to reference FASTA
  -b BUILD                      FASTA build assembly
  -s SKIP                       Number of rows to skip
  -m MAX_MISSING                Maximum fraction of missing variants permitted
  -chrom_field CHROM_FIELD      Column number for chromosome
  -pos_field POS_FIELD          Column number for chromosome
  -ea_field EA_FIELD            Column number for effect allele
  -nea_field NEA_FIELD          Column number for non-effect allele
  -effect_field EFFECT_FIELD    Effect size field
  -se_field SE_FIELD            SE field
  -pval_field PVAL_FIELD        P-Value field
  -n_field N_FIELD              Number of samples field
  -dbsnp_field DBSNP_FIELD      dbSNP identifier field
  -ea_af_field EA_AF_FIELD      Effect allele frequency field
  -nea_af_field NEA_AF_FIELD    None effect allele frequency field
```

## Full workflow on BC3

```
# edit as needed #
bash make_vcf.bc3.sh```
