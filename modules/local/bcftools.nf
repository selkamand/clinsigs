#!/usr/bin/env nextflow
nextflow.enable.types = true

include { Clonality } from '../../types.nf'



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

    bcftools view -f PASS ${vcf} -Oz -o ${vcf.baseName}.pass.vcf.gz
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
        reference_genome_fasta: Path,
        reference_genome_fai: Path
    )

    output:
    record(
        sample: sample,
        vcf: file("${vcf.baseName}.normalised.vcf.gz"),
    )

    script:
    """
    set -euo pipefail

    bcftools norm \
        --multiallelic - \
        --check-ref e \
        --output-type z \
        --output "${vcf.baseName}.normalised.vcf.gz" \
        -f ${reference_genome_fasta} \
        ${vcf} 
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
        vcf: Path,
        clonality: Clonality
    )

    output:
    record(
        sample: sample,
        clonality: clonality,
        tsv: file("${vcf.baseName}.tsv"),
    )

    script:
    """
    set -euo pipefail

    # Note bcftools view is to isolate the tumour sample 
    # (sometimes VCFs will include germline and tumour sample)
    # Tumor sample name in VCF expected to match sample
    bcftools view -s ${sample} ${vcf} | \
    bcftools query \
        -HH -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/SUBCL\t%INFO/PURPLE_AF[\t%DP\t%AF]\n' \
        > "${vcf.baseName}.tsv"
    """
}


process PURPLE_SV_VCF_TO_BEDPE {
    tag "${sample}"

    input:
    record(
        sample: String,
        vcf: Path
    )

    output:
    record(
        sample: sample,
        bedpe: file("${vcf.baseName}.bedpe"),
    )

    script:
    """
    set -euo pipefail

    svcf -i ${vcf} --from purple --to bedpe > "${vcf.baseName}.bedpe"
    """
}
