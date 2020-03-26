# GWAS to VCF harmonisation tool

Tool to map GWAS summary statistics to VCF/BCF with on-the-fly harmonisation to a supplied reference FASTA

## Quick start

Use web interface [gwas2vcfweb](https://github.com/mrcieu/gwas2vcfweb)

## Run locally

Either run directly on a UNIX host or using Docker containerisation (recommended)

### Download

```sh
git clone git@github.com:MRCIEU/gwas2vcf.git
cd gwas2vcf
```

#### Native

```sh
virtualenv venv
source ./venv/bin/activate
./venv/bin/pip install -r requirements.txt
./venv/bin/python main.py -h
```

#### Docker

Pull existing image from DockerHub or build

```sh
docker pull mcgml/gwas2vcf
# OR
docker build -t gwas2vcf .
```

Run

```sh
docker run \
-v /path/to/fasta:/path/to/fasta \
-name gwas2vcf \
-it gwas2vcf:latest \
python main.py -h
```

### Reference FASTA

```sh
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

### Running the tests

Unit tests:

```sh
cd gwas2vcf
python -m pytest -v test
```

### Usage

```sh
usage: main.py [-h] [-v] [--out OUT] [--data DATA] --ref REF --json JSON [--id ID] [--cohort_controls COHORT_CONTROLS]
               [--cohort_cases COHORT_CASES] [--rm_chr_prefix] [--csi] [--log {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

Map GWAS summary statistics to VCF/BCF

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --out OUT             Path to output VCF/BCF. If not present then must be specified as 'out' in json file
  --data DATA           Path to GWAS summary stats. If not present then must be specified as 'data' in json file
  --ref REF             Path to reference FASTA
  --json JSON           Path to parameters JSON
  --id ID               Study identifier. If not present then must be specified as 'id' in json file
  --cohort_controls COHORT_CONTROLS
                        Total study number of controls (if case/control) or total sample size if continuous. Overwrites value if
                        present in json file.
  --cohort_cases COHORT_CASES
                        Total study number of cases. Overwrites value if present in json file.
  --rm_chr_prefix       Remove chr prefix from GWAS chromosome
  --csi                 Default is to index tbi but use this flag to index csi
  --log {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging level
```

Additional parameters are passed through a [JSON](https://www.w3schools.com/js/js_json_objects.asp) parameters file, see `param.py` for full details or example below.

```json
{
  "chr_col": 0,
  "pos_col": 1,
  "ea_col": 2,
  "oa_col": 3,
  "beta_col": 4,
  "se_col": 5,
  "pval_col": 6,
  "snp_col": 7,
  "eaf_col": 8,
  "imp_info_col": 9,
  "ncontrol_col": 10,
  "delimiter": "\t",
  "header": true,
  "build": "GRCh37"
}
```

### Example

See [gwas-vcf-performance](https://github.com/MRCIEU/gwas-vcf-performance/blob/master/workflow.Rmd) for a full implementation 

### Combine multiallelics

Merge variants at single genetic position on to a single row. This step is **highly** recommended to avoid duplicate RSIDs which is not supported by VCF

```sh
bcftools norm \
-f ref.fasta \
-m +any \
-O z \
-o norm.vcf.gz
```

### Validate VCF file

Check the file format is valid but ignore genotypes since there are none

```sh
gatk ValidateVariants \
-V harmonised.vcf \
-R ref.fasta \
--dbsnp dbsnp.vcf \
--validation-type-to-exclude ALLELES
```

### Add variant frequency and dbSNP identifiers

Add variant frequency from 1000 genomes (or similar)

```sh
bcftools annotate \
-a 1kg.vcf.gz \
-c ID,AF \
-O z \
-o annotated.vcf.gz \
harmonised.vcf.gz
```

### Extract genome-wide significant variants to text file

```sh
bcftools query \
-i 'FORMAT/LP > 7.3' \
-f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%ES\t%SE\t%LP]\n' \
-o data.txt \
file.vcf.gz
```

### Merge multiple GWAS summary stats into a single file

```sh
bcftools merge \
-O z \
-o merged.vcf.gz \
*.vcf.gz
```
