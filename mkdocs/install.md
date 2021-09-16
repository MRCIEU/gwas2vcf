# Install

## Quick start

Use web interface <a href="http://vcf.mrcieu.ac.uk" target="_blank">http://vcf.mrcieu.ac.uk</a>

## Run locally

Either run directly on a UNIX host or using Docker containerisation (recommended)

### Download

```sh
git clone https://github.com/MRCIEU/gwas2vcf.git
cd gwas2vcf
```

### Native

Requires Python v3.8

```sh
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
pip install git+git://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph
python main.py -h
```

### Docker

Pull image from DockerHub OR build image from source

```sh
# pull image from DockerHub
docker pull mrcieu/gwas2vcf

### OR ###

# build docker image from source
docker build -t gwas2vcf .
```

Run

```sh
docker run \
-v /path/to/fasta:/path/to/fasta \
--name gwas2vcf \
-it mrcieu/gwas2vcf:latest \
python main.py -h
```

### Usage

```sh
usage: main.py [-h] [-v] [--out OUT] [--data DATA] --ref REF [--dbsnp DBSNP] --json JSON [--id ID] [--cohort_controls COHORT_CONTROLS]
               [--cohort_cases COHORT_CASES] [--csi] [--log {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--alias ALIAS]

Map GWAS summary statistics to VCF/BCF

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program version number and exit
  --out OUT             Path to output VCF/BCF. If not present then must be specified as 'out' in json file
  --data DATA           Path to GWAS summary stats. If not present then must be specified as 'data' in json file
  --ref REF             Path to reference FASTA
  --dbsnp DBSNP         Path to reference dbSNP VCF
  --json JSON           Path to parameters JSON
  --id ID               Study identifier. If not present then must be specified as 'id' in json file
  --cohort_controls COHORT_CONTROLS
                        Total study number of controls (if case/control) or total sample size if continuous. Overwrites value if present in json
                        file.
  --cohort_cases COHORT_CASES
                        Total study number of cases. Overwrites value if present in json file.
  --csi                 Default is to index tbi but use this flag to index csi
  --log {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the logging level
  --alias ALIAS         Optional chromosome alias file
```

Additional parameters are passed through a [JSON](https://www.w3schools.com/js/js_json_objects.asp) parameters file using ```--json <param.json>```, see `param.py` for full details and below example. Note that field columns start at 0.

### Running the tests

Unit tests:

```sh
cd gwas2vcf
python -m pytest -v test
```

## Reference files

# FASTA

```sh
# GRCh36/hg18/b36
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.fasta.gz; gzip -d human_b36_both.fasta.gz
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.fasta.fai
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.dict

# GRCh37/hg19/b37
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.gz; gzip -d human_g1k_v37.fasta.gz
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.dict

# GRCh38/hg38/b38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
```

# dbSNP

```sh
# GRCh37/hg19/b37
wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.gz .
wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.gz.tbi .

# GRCh38/hg38/b38
wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz .
wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz.tbi .
```

Newer dbSNP builds can be obtained from the NCBI FTP but the VCF files have non-standard chromosome names which can be updated accordingly (thanks @darked89)

```
# download latest dbSNP VCF for hg38
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz

# define chromosome name mapping
echo -n "NC_000001.11    1
NC_000002.12    2
NC_000003.12    3
NC_000004.12    4
NC_000005.10    5
NC_000006.12    6
NC_000007.14    7
NC_000008.11    8
NC_000009.12    9
NC_000010.11    10
NC_000011.10    11
NC_000012.12    12
NC_000013.11    13
NC_000014.9     14
NC_000015.10    15
NC_000016.10    16
NC_000017.11    17
NC_000018.10    18
NC_000019.10    19
NC_000020.11    20
NC_000021.9     21
NC_000022.11    22
NC_000023.11    X
NC_000024.10    Y
NC_012920.1     MT
" > hg38_rename_chrom_names.tsv

# update chromosome names
bcftools annotate \
--rename-chrs hg38_rename_chrom_names.tsv \
--output-type z \
--output dbSNP_clean.vcf.gz GCF_000001405.39.gz

# index modified file
bcftools index dbSNP_clean.vcf.gz
```
