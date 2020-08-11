# Convert GWAS summary statistics to VCF

<!-- badges: start -->
[![Build Status](https://travis-ci.org/MRCIEU/gwas2vcf.svg?branch=master)](https://travis-ci.org/MRCIEU/gwas2vcf)
<!-- badges: end -->

Tool to map GWAS summary statistics to VCF/BCF with on-the-fly harmonisation to a supplied reference FASTA

Produces GWAS-VCF with version 1.0 of the [specification](https://github.com/MRCIEU/gwas-vcf-specification/releases/tag/1.0.0)

## Citation

```
The variant call format provides efficient and robust storage of GWAS summary statistics
Matthew S Lyon, Shea J Andrews, Benjamin L Elsworth, Tom R Gaunt, Gibran Hemani, Edoardo Marcora
bioRxiv 2020.05.29.115824; doi: https://doi.org/10.1101/2020.05.29.115824
```

## Quick start

Use web interface <http://vcf.mrcieu.ac.uk>

## Run locally

Either run directly on a UNIX host or using Docker containerisation (recommended)

### Download

```sh
git clone git@github.com:MRCIEU/gwas2vcf.git
cd gwas2vcf
```

#### Native

```sh
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
pip install git+git://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph
python main.py -h
```

#### Docker

Build docker image

```sh
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
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.fasta
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.fasta.fai

# GRCh37/hg19/b37
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai

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

Additional parameters are passed through a [JSON](https://www.w3schools.com/js/js_json_objects.asp) parameters file using ```--json <param.json>```, see `param.py` for full details and below example. Note that field columns start at 0.

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

## Working with GWAS-VCF

See below examples of working with GWAS-VCF. Let us know if you have other use cases through the issues page!

### Parsing libraries

See [R](https://github.com/mrcieu/gwasvcf) and [Python](https://github.com/mrcieu/pygwasvcf) libraries for reading GWAS summary statistics in GWAS-VCF

### Command-line manipulation

The following examples require:
- [bcftools](http://samtools.github.io/bcftools/)
- [gatk](https://gatk.broadinstitute.org/)
- [htslib](http://www.htslib.org/download/)

Please cite the relevant tool(s) if you use these examples.

#### Filter

##### Extract genome-wide significant variants

The LP field is -log10(P), 7.3 is approx 5e-8

```sh
bcftools filter \
-i 'FORMAT/LP > 7.3' \
-o output.vcf \
file.vcf.gz
```

##### Extract variants by gene

Requires annotation by Ensembl (see below)

```sh
bcftools filter \
-i 'INFO/ENSG_ID == "ENSG00000198670"' \
file.vcf.gz
```

##### Extract variants by pathway

Requires annotation by Reactome (see below)

```sh
bcftools filter \
-i 'INFO/Reactome_ID == "R-HSA-3000171"' \
file.vcf.gz
```

##### Select genome region for further analysis

 ```sh
bcftools filter \
-r 1:1000000-2000000 \
-o output.vcf.gz \
input.vcf.gz
```

#### Annotate

##### Add variant frequency

```sh
# download 1000 genomes phase 3 (hg19/GRCh37) allele frequencies and index
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi 

# annotate GWAS-VCF with 1kg allele frequencies
bcftools annotate \
-a ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz \
-c AF \
-O z \
-o output.vcf.gz \
input.vcf.gz
```

##### Add gene/pathway annotations

Download and merge input files

```sh
# download ensembl-to-position mapping and sort by ensgId
curl ftp://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz | \
gzip -dc | \
awk -F"\t|:|;|=" '$3=="gene" && $1 >= 1 && $1 <= 22 {print $11"\t"$1"\t"$4"\t"$5}' | \
sort -k1,1 > Ensembl2Position.sorted.txt

# download ensembl-to-pathway mapping and sort by ensgId
curl https://reactome.org/download/current/Ensembl2Reactome.txt | \
grep "Homo sapiens" | \
cut -s -f1,2 | \
sort -k1,1 > Ensembl2Reactome.sorted.txt

# merge tables by ensgId, sort, compress and index
join \
-t $'\t' \
--check-order \
Ensembl2Position.sorted.txt \
Ensembl2Reactome.sorted.txt | \
awk -F"\t" '{print $2"\t"$3-1"\t"$4"\t"$5}' | \
sort -k1,1V -k2,2n -k3,3n | \
bgzip -c > reactome.bed.gz
tabix -p bed reactome.bed.gz

# sort, compress and index ensembl-to-position mapping
awk -F"\t" '{print $2"\t"$3-1"\t"$4"\t"$1}' Ensembl2Position.sorted.txt | \
sort -k1,1V -k2,2n -k3,3n  | \
bgzip -c > ensembl.bed.gz
tabix -p bed ensembl.bed.gz
```

Annotate GWAS-VCF with gene ID

```sh
# annotate GWAS-VCF
bcftools annotate \
-a ensembl.bed.gz \
-c CHROM,FROM,TO,ENSG_ID \
-h <(echo '##INFO=<ID=ENSG_ID,Number=.,Type=String,Description="Ensembl gene ID">') \
-o output.vcf.gz \
-O z \
-l ENSG_ID:unique \
input.vcf.gz
```

Annotate GWAS-VCF with Reactome ID

```sh
# annotate GWAS-VCF
bcftools annotate \
-a reactome.bed.gz \
-c CHROM,FROM,TO,Reactome_ID \
-h <(echo '##INFO=<ID=Reactome_ID,Number=.,Type=String,Description="Reactome ID">') \
-o output.vcf.gz \
-O z \
-l Reactome_ID:unique \
input.vcf.gz
```

#### Convert

##### Export to NHGRI-EBI GWAS catalog format

```sh
# map to GWAS catalog format
bcftools query \
-e 'ID == "."' \
-f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t%AF\t[%ES\t%SE]\n' \
gwas.vcf.gz | \
awk 'BEGIN {print "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' > gwas.tsv

# validate file using [ss-validate](https://pypi.org/project/ss-validate)
ss-validate -f gwas.tsv
```

#### Liftover

##### Map genomic coordinates to another genome build (liftover)

This procedure requires a chain file which contains the chromosome base-position mapping between two genome builds

```sh
# download chain file
wget http://fileserve.mrcieu.ac.uk/ref/chains/b36tob37.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/b37tob36.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/b37tohg18.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/b37tohg19.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/hg18tob37.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/hg19toHg18.chain
```

```sh
# perform liftover
 gatk LiftoverVcf \
--INPUT input.vcf.gz \
--OUTPUT output.vcf.gz \
--REJECT rejected.vcf.gz \
--CHAIN file.chain \
--REFERENCE_SEQUENCE target.fasta \
--RECOVER_SWAPPED_REF_ALT false 
```

#### Merge

##### Combine multiple GWAS summary stats into a single file

This is useful for distributing QTL/molecular phenotype GWAS

```sh
bcftools merge \
-O z \
-o merged.vcf.gz \
*.vcf.gz
```

#### Validate

##### Check the file format is valid

```sh
gatk ValidateVariants \
-V input.vcf.gz \
-R ref.fasta \
--validation-type-to-exclude ALLELES
```