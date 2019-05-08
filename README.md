# GWAS Harmonisation

Tool to map GWAS summary statistics to VCF/BCF with on the fly harmonisation to a supplied reference FASTA.

### Install

```
# Src
git clone https://github.com/MRCIEU/gwas_harmonisation
cd gwas_harmonisation
```

```
# Native
pip install --user -r ./requirements.txt

# Docker
docker build -t gwas_harmonisation_wdl .
docker create -v /data/bgc/ref:/data/ref -v /data/bgc:/data -name gwas_harmonisation_wdl gwas_harmonisation_wdl 
```

### Reference FASTA

```
# GRCh36/hg18/b36
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b36/human_b36_both.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b36/human_b36_both.fasta.fai.gz

# GRCh37/hg19/b37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz

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
bash test/example.sh
```

## Usage

Column field numbers are 0-based

```
Map GWAS summary statistics to VCF/BCF

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit
  --out OUT      Path to output VCF/BCF
  --data GWAS    Path to GWAS summary stats
  --ref FASTA    Path to reference FASTA
  --json JSON    Path to parameters JSON
```

See param.py for JSON specification

## Dealing with missing variant frequency

Add variant frequency from 1000 genomes (or similar)

```
gatk VariantAnnotator \
-R ref.fasta \
-V harmonised.bcf \
-O harmonised_af.bcf \
--resource 1kg:1kg.vcf.gz \
-E 1kg.EUR_AF \
--resource-allele-concordance
``` 
