# Tutorial

Assuming the GWAS summary stats have a hg19/b37 chromosome name & position you can use these files:

* [Reference FASTA](install.md#fasta)
* [dbSNP VCF](install.md#dbsnp)
* [Alias File](https://raw.githubusercontent.com/MRCIEU/gwas2vcf/master/alias-b37.txt)

## Download GWAS

```sh
# obtain test gwas summary stats
wget https://raw.githubusercontent.com/MRCIEU/gwas2vcfweb/master/app/tests/data/example.1k.txt
```

## Create parameters file

```json
{
  "chr_col": 0,
  "pos_col": 1,
  "snp_col": 2,
  "ea_col": 3,
  "oa_col": 4,
  "beta_col": 5,
  "se_col": 6,
  "ncontrol_col": 7,
  "pval_col": 8,
  "eaf_col": 9,
  "delimiter": "\t",
  "header": true,
  "build": "GRCh37"
}
```

## Map GWAS summary stats to GWAS-VCF

```sh
SumStatsFile=/data/example.1k.txt
RefGenomeFile=/data/human_g1k_v37.fasta
ParamFile=/data/params.json
DbSnpVcfFile=/data/dbsnp.v153.b37.vcf.gz
VcfFileOutPath=/data/out.vcf
ID="test"

python /app/main.py \
--data ${SumStatsFile} \
--json ${ParamFile} \
--id ${ID} \
--ref ${RefGenomeFile} \
--dbsnp ${DbSnpVcfFile} \
--out ${VcfFileOutPath} \
--alias /app/alias-b37.txt
```

### Alias file

Some genome builds use the chr prefix on chromosome names i.e. ```chr1``` while others just use ```1```. This will cause issues if your GWAS summary statistics and FASTA you wish to map to have different chromosome names (although they use the same genome build).

One solution is to provide an alias file to map your GWAS summary stats chromosome name to another string. An example alias is provided in the repo ```alias-b37.txt``` and ```alias-hg38.txt```. The format is ```source-chr\tdest-chr``` one row per contig.