# Convert GWAS summary statistics to VCF

<!-- badges: start -->
![Build Status](https://github.com/MRCIEU/gwas2vcf/actions/workflows/test.yml/badge.svg)
[![DOI](https://img.shields.io/badge/doi-10.1186%2Fs13059--020--02248--0-blue)](https://doi.org/10.1186/s13059-020-02248-0)
<!-- badges: end -->

Tool to map GWAS summary statistics to VCF/BCF with on-the-fly harmonisation to a supplied reference FASTA

Produces GWAS-VCF with version 1.0 of the [specification](https://github.com/MRCIEU/gwas-vcf-specification/releases/tag/1.0.0)

## Documentation

Full documentation available from <https://mrcieu.github.io/gwas2vcf>

## GWAS summary statistics

Complete GWAS summary statistics harmonised to GWAS-VCF are available on >14,000 datasets available from the [OpenGWAS project](https://gwas.mrcieu.ac.uk)

## Tutorials

What can I do with GWAS-VCF?

- Command-line operations ([filter](https://mrcieu.github.io/gwas2vcf/downstream/#filter), [annotate](https://mrcieu.github.io/gwas2vcf/downstream/#annotate), [convert](https://mrcieu.github.io/gwas2vcf/downstream/#convert), [liftover](https://mrcieu.github.io/gwas2vcf/downstream/#liftover), [merge](https://mrcieu.github.io/gwas2vcf/downstream/#merge), [validate](https://mrcieu.github.io/gwas2vcf/downstream/#validate))
- Quality control ([QC](https://github.com/MRCIEU/opengwas-reports))
- Analysis workflows ([pathway](https://mrcieu.github.io/gwas2vcf/downstream/#extract-variants-by-pathway), [gene](https://mrcieu.github.io/gwas2vcf/downstream/#extract-variants-by-gene), [clumping](https://mrcieu.github.io/gwasglue/articles/finemapping.html), [finemapping](https://mrcieu.github.io/gwasglue/articles/finemapping.html), [LDSC](https://github.com/MRCIEU/gwas_processing), [colocalization](https://mrcieu.github.io/gwasglue/articles/colocalisation.html), [Mendelian randomization](https://mrcieu.github.io/gwasglue/articles/mr.html#using-gwas-vcf-files), [COJO](https://mrcieu.github.io/gwasglue/articles/cojo.html))
- File parsing ([R](https://github.com/MRCIEU/gwasvcf), [Python3](https://github.com/MRCIEU/pygwasvcf))

Let us know if you have other use cases through the issues page!

## Citation

Lyon M, Andrews S, Elsworth B, Gaunt T, Hemani G, Marcora E. The variant call format provides efficient and robust storage of GWAS summary statistics. Genome Biol 22, 32 (2021). <https://doi.org/10.1186/s13059-020-02248-0>
