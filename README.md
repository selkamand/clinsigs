# clinsigs

> [!WARNING]
> This pipeline is in development and not yet ready for use

## Problem Statement

This nextflow pipeline is inspired by the research pipelines [clinsigs](https://github.com/selkamand/clinsigs) and [sigscreen](https://github.com/selkamand/sigscreen) but is both more stable and more minimal so that it can one day be integrated in clinical pipelines.



## Design and Scope 

clinsigs is a nextflow pipeline for detecting mutational signatures in cancer. It requires you have at least one tumour biopsy sample characterised with any of the following technologies:

1. Whole-genome sequencing (tumour-normal)
2. RNAseq (total or polyA)
3. methylation arrays

clinsigs should be run after fundamental cancer analyses (alignment, variant calling, etc) are performed with oncoanalyser or similar tools. clinsigs seeks to infer as many characteristics of the tumour as possible from an single patient multiome.


> [!WARNING]
> Do NOT expect API stability. This tool runs a frankly silly combination of analyses, from driver mutation analysis, mutational signatures, various gene-dysfunction classifiers, circos visualisations, etc. As we figure out whats useful and whats not, a cleaner pipeline will be established in a separate repo.

If you're after cleaner implementations of these cancer analyses: check out the [sigverse](https://github.com/CCICB/sigverse) and [scarscape](https://github.com/CCICB/scarscape) that power some clinsigs analyses.

## Running the Pipeline

Start with a samplesheet (csv) with the following columns (include a header row!)

1. **sample**: sample identifier
2. **snv**: path to a somatic DNA snv vcf (ideally an oncoanalyser purple-annotated vcf) 
3. **sv**: path to a somatic DNA sv vcf (ideally an oncoanalyser purple-annotated vcf)
4. **cnv**: path a segment file (tsv)
5. **rnamut**: path to an RNA mutation vcf
6. **meth**: path to a methylation probe beta file (csv) with 2 cols (probe,beta)
7. **teal**: somatic telomere length estimated by teal
8. **bam**: somatic tumour bam (should be coord sorted and indexed with samtools index)

If your missing a modality, thats ok, just leave the field empty. Pipeline should dynamically adapt.

Once you have the samplesheet, you can run clinsigs as follows: 

```
nextflow run selkamand/clinsigs \
  -profile docker \
  -params params.yaml \
  --input samplesheet.csv
```

Due to the number of parameters required by clinsigs we recommend configuring the run with a param.yaml file. An example is provided [here](params.yaml). Expect to spend some time configuring your reference datasets. You can absolutely run this pipeline on a single-sample, but you get a lot more out if you load up reference matrices.


## What does clinsigs do?


### VCF Normalisation, Filtering, and Splitting by Clonality

1. **FILTER**: Filter snv, sv, and rnamut vcfs for PASS variants 
2. **Normalise** snv and rnamut VCFs: 
  - Left-align INDELS
  - Check REF alleles match reference genome
  - Split multiallelic sites into separate rows. 
  - We do NOT decompose complex variants (e.g. MNVs do NOT become consecutive SNVs). This is very important for downstream variant classification (e.g. doublet & indel signature analysis). 
3. For somatic SNVs, spawn 3 VCFs per sample - one containing only subclonal variants, one containing clonal variants, and one containing all variants. This allows mutation-signature analysis and other analyses to compare clonal vs subclonal profiles later in the pipeline.
4. Convert VCFs into tabular 1-based single position file format.

We then convert them to tabular, 1-based, file formats that feed into downstream analyses

### Signature Analysis (Reference signature fitting)

There are two steps in mutsig analysis. 
1. Classify mutations based on *classification schemes* and **tally** different types of mutations
2. **Fit** a model to tally matrices to predict which processes spawned those mutations.


### Signature Analysis (De novo discovery)

Traditional de novo signature analysis requires a large cohort. 

This pipeline however does try to flag samples with 'novel' signatures present. 

**Search strategy:**
1. Attempt to fit signatures while relaxing constraints that prevent overfitting.
2. Quantify goodness of fit. 
3. If mutation profile could not be reconstructed from known signatures, despite the sample having a decent number of mutations, then its profile may have derived from a mutagenic process we don't have a signature for (at least in databases used)



## Building the dockerfiles

From inside this directory navigate to `dockerfiles/<toolname>`.

Build local version for OSX

```{bash}
docker buildx build --platform linux/arm64 --load --tag selkamandcci/<toolname>:0.0.1 .
```

Build final version to push to dockerhub

```{bash}
docker buildx build --push --platform linux/amd64,linux/arm64 --tag selkamandcci/<toolname>:0.0.1 .
```

Repeat for other docker images (e.g `dockers/quast`)

## Testing the pipeline

```
nextflow run selkamand/clinsigs -profile docker,test
```

Replace docker with `singularity` / `apptainer` / etc. depending on your platform

