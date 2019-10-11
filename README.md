# GWAS to VCF harmonisation tool

Tool to map GWAS summary statistics to VCF/BCF with on the fly harmonisation to a supplied reference FASTA.

### Install

```
# src
git clone git@github.com:MRCIEU/gwas2vcf.git
cd gwas2vcf

# VirtualEnv
virtualenv venv
source ./venv/bin/activate
./venv/bin/pip install -r requirements.txt
./venv/bin/python main.py -h

# Docker
docker build -t gwas2vcf .
docker create -v /data:/data -name gwas2vcf gwas2vcf
docker run -it gwas2vcf:latest python main.py -h
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
cd gwas2vcf
python -m unittest discover test
```

End-to-end tests:

```
bash test/example.sh
```

## Usage

```
Map GWAS summary statistics to VCF/BCF

-h, --help            show this help message and exit
-v, --version         show program's version number and exit
--out OUT             Path to output VCF/BCF
--data GWAS           Path to GWAS summary stats
--ref FASTA           Path to reference FASTA
--json JSON           Path to parameters JSON
--id ID               Study identifier
--cohort_controls COHORT_CONTROLS
                      Total study number of controls (if case/control) or
                      total sample size if continuous
--cohort_cases COHORT_CASES
                      Total study number of cases
--rm_chr_prefix       Remove chr prefix from GWAS chromosome
--log {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                      Set the logging level
```

See param.py for JSON specification

## Combine multiallelics

Merge variants at single genetic position on to a single row

```
bcftools norm \
-f ref.fasta \
-m +any \
-O z \
-o norm.vcf.gz
```

## Validate VCF file

Check the file format is valid but ignore genotypes since these are missing

```
gatk ValidateVariants \
-V harmonised.vcf \
-R ref.fasta \
--dbsnp dbsnp.vcf \
--validation-type-to-exclude ALLELES
```

## Add variant frequency and dbSNP identifiers

Add variant frequency from 1000 genomes (or similar)

```
bcftools annotate \
-a 1kg.vcf.gz \
-c ID,AF \
-O z \
-o annotated.vcf.gz \
harmonised.vcf.gz
``` 

## Extract genome-wide significant variants to text file

```
bcftools query \
-i 'FORMAT/LP > 7.3' \
-f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP]\n' \
-o data.txt \
file.vcf.gz
```

## Merge multiple GWAS summary stats into a single file

Note: Merged GWAS BCFs are significantly slower to query; for best performance do not do this.

```
bcftools merge \
-O z \
-o merged.vcf.gz \
*.vcf.gz
```

## Known issues

VCF v4.2 cannot accommodate double precision floats. Decimals smaller than 1.1754944e-38 are rounded to 0. P values are encoded as -log10.
