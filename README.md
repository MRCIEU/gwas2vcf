# GWAS Harmonisation

Tool to map GWAS summary statistics from tab-delim to VCF with on the fly harmonisation to a supplied reference FASTA.

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
python -m unittest discover gwas_harmonisation/test
```

End-to-end tests:

```
# convert GIANT to BCF
bash test/make_vcf.bmi.sh

# convert GLGC to BCF
bash test/make_vcf.bmi.sh
```

## Usage

Column field numbers are 0-based

```
python gwas_harmonisation/main.py -h
```

## Full workflow

```bash make_vcf.ldl.sh```
