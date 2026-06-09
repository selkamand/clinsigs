#!/usr/bin/env nextflow
nextflow.enable.types = true

process BCFTOOLS_PASS {

    tag "${sample}"

    input:
    record(
        sample: String,
        vcf: Path
    )

    output:
    record(
        sample: sample,
        vcf: file("${vcf.baseName}.pass.vcf.gz"),
    )

    script:
    """
    set -euo pipefail

    bcftools view -f PASS ${vcf} > -Oz -o ${vcf.baseName}.pass.vcf.gz
    """
}


// Left-Normalise INDELS, split multiallelics into multiple lines 
// Importantly, does NOT split complex variants (like MNVs) into individual lines. 
process BCFTOOLS_NORMALISE_INDELS_AND_SPLIT_MULTIALLELICS {

    tag "${sample}"

    input:
    record(
        sample: String,
        vcf: Path
    )
    record(
        ref_fa: Path,
        reference_genome_fai: Path
    )

    output:
    record(
        sample: sample,
        vcf: file("${vcf.baseName}.normalise.vcf.gz"),
    )

    script:
    """
    set -euo pipefail

    bcftools norm \
        --multiallelic - \
        --check-ref e \
        --verbose 1 \
        --output-type z \
        --output "${vcf.baseName}.normalised.vcf.gz" \
        -f ${ref_fa} \
        ${vcf} 
    -f PASS ${vcf} > -Oz -o ${vcf.baseName}.pass.vcf.gz
    """
}

// Split a VCF into clonal vs subclonal based on value of INFO field 'SUBCL' (float).
// This is automatically annotated in purple SNV VCFs
process SPLIT_PURPLE_SNVS_BY_CLONALITY {
    tag "${sample}"

    input:
    record(
        sample: String,
        vcf: Path
    )

    output:
    record(
        sample: sample,
        all: file("${vcf.baseName}.all.vcf.gz"),
        subclonal: file("${vcf.baseName}.subclonal.vcf.gz"),
        clonal: file("${vcf.baseName}.clonal.vcf.gz"),
    )

    script:
    """
    set -euo pipefail
    
    # Get Subclonal / Clonal
    bcftools view -i 'INFO/SUBCL < 0.3' ${vcf} -Oz -o ${vcf.baseName}.clonal.vcf.gz
    bcftools view -e 'INFO/SUBCL < 0.3' ${vcf} -Oz -o ${vcf.baseName}.subclonal.vcf.gz
    bcftools view ${vcf} -Oz -o ${vcf.baseName}.all.vcf.gz
    """
}

// Convert a PURPLE SNV VCF to a TSV file. 
process PURPLE_SNV_VCF_TO_TSV {
    tag "${sample}"

    input:
    record(
        sample: String,
        vcf: Path
    )

    output:
    record(
        sample: sample,
        tsv: file("${vcf.baseName}.tsv"),
    )

    script:
    """
    set -euo pipefail

    bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/SUBCL\t%INFO/DP\t%INFO/PURPLE_AF\t%INFO/AF\n' \
    ${vcf} > ${vcf.baseName}.tsv
    """
}
