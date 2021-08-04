# Working with GWAS-VCF

See below examples of working with GWAS-VCF. Let us know if you have other use cases through the [issues](https://github.com/MRCIEU/gwas2vcf/issues) page!

## Data

Full GWAS summary statistics harmonised to GWAS-VCF are available on >14,000 datasets available from <a href="https://gwas.mrcieu.ac.uk/" target="_blank">https://gwas.mrcieu.ac.uk/</a>

## Parsing libraries

See <a href="https://github.com/mrcieu/gwasvcf" target="_blank">R</a> and <a href="https://github.com/mrcieu/pygwasvcf" target="_blank">Python3</a> libraries for reading GWAS summary statistics in GWAS-VCF

## Command-line manipulation

The following examples require:

* [bcftools](http://samtools.github.io/bcftools/)
* [gatk](https://gatk.broadinstitute.org/)
* [htslib](http://www.htslib.org/download/)

Please cite the relevant tool(s) if you use these examples.

### Filter

#### Extract genome-wide significant variants

The LP field is -log10(P), 7.3 is approx 5e-8

```sh
bcftools filter \
-i 'FORMAT/LP > 7.3' \
-o output.vcf \
file.vcf.gz
```

#### Extract variants by gene

Requires annotation by Ensembl (see below)

```sh
bcftools filter \
-i 'INFO/ENSG_ID == "ENSG00000198670"' \
file.vcf.gz
```

#### Extract variants by pathway

Requires annotation by Reactome (see below)

```sh
bcftools filter \
-i 'INFO/Reactome_ID == "R-HSA-3000171"' \
file.vcf.gz
```

#### Select genome region for further analysis

```sh
bcftools filter \
-r 1:1000000-2000000 \
-o output.vcf.gz \
input.vcf.gz
```

### Annotate

#### Add variant frequency

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

#### Add gene/pathway annotations

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

### Convert

#### Export to NHGRI-EBI GWAS catalog format

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

### Liftover

#### Map genomic coordinates to another genome build (liftover)

This procedure requires a chain file which contains the chromosome base-position mapping between two genome builds

```sh
# download chain file
wget http://fileserve.mrcieu.ac.uk/ref/chains/b36tob37.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/b37tob36.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/b37tohg18.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/b37tohg19.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/hg18tob37.chain
wget http://fileserve.mrcieu.ac.uk/ref/chains/hg19toHg18.chain
curl https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz | zcat | sed 's/chr//' > b37ToHg38.chain
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

### Merge

#### Combine multiple GWAS summary stats into a single file

This is useful for distributing QTL/molecular phenotype GWAS

```sh
bcftools merge \
-O z \
-o merged.vcf.gz \
*.vcf.gz
```

### Validate

#### Check the file format is valid

```sh
gatk ValidateVariants \
-V input.vcf.gz \
-R ref.fasta \
--validation-type-to-exclude ALLELES
```
