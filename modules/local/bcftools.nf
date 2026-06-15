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
        vcf: file("${vcf.simpleName}.pass.vcf.gz"),
    )

    script:
    """
    set -euo pipefail

    bcftools view -f PASS ${vcf} -Oz -o ${vcf.simpleName}.pass.vcf.gz
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
        vcf: file("${vcf.simpleName}.normalised.vcf.gz"),
    )

    script:
    """
    set -euo pipefail

    bcftools norm \
        --multiallelic - \
        --check-ref e \
        --output-type z \
        --output "${vcf.simpleName}.normalised.vcf.gz" \
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
        all: file("${vcf.simpleName}.all.vcf.gz"),
        subclonal: file("${vcf.simpleName}.subclonal.vcf.gz"),
        clonal: file("${vcf.simpleName}.clonal.vcf.gz"),
    )

    script:
    """
    set -euo pipefail
    
    # Get Subclonal / Clonal
    bcftools view -i 'INFO/SUBCL < 0.3' ${vcf} -Oz -o ${vcf.simpleName}.clonal.vcf.gz
    bcftools view -e 'INFO/SUBCL < 0.3' ${vcf} -Oz -o ${vcf.simpleName}.subclonal.vcf.gz
    bcftools view ${vcf} -Oz -o ${vcf.simpleName}.all.vcf.gz
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
        tsv: file("${vcf.simpleName}.snv.${clonality}.tsv"),
    )

    script:
    """
    set -euo pipefail

    bioprep vcf -i ${vcf} --from purple --to tsv > ${vcf.simpleName}.snv.${clonality}.tsv
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
        bedpe: file("${vcf.simpleName}.breakpoints.bedpelike.tsv"),
    )

    script:
    """
    set -euo pipefail

    bioprep svcf -i ${vcf} --from purple --to bedpe > "${vcf.simpleName}.breakpoints.bedpelike.tsv"
    """
}

process PURPLE_SV_VCF_TO_BREAKEND_TSV {
    tag "${sample}"

    input:
    record(
        sample: String,
        vcf: Path
    )

    output:
    record(
        sample: sample,
        breakends_tsv: file("${vcf.simpleName}.breakends.tsv"),
    )

    script:
    """
    set -euo pipefail

    bioprep svcf -i ${vcf} --from purple --to breakend-tsv > "${vcf.simpleName}.breakends.tsv"
    """
}
