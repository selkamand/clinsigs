nextflow.enable.types = true

// include { Clonality ; Vartype } from '../../types.nf'

// TSV is purple copynumber segment file
process PREPROCESS_CNV_SEGMENTS {
    input:
    record(
        sample: String,
        tsv: Path
    )

    output:
    record(
        sample: sample,
        tsv: file("${tsv.baseName}.cleaned.segment.tsv"),
    )

    script:
    """
    cat ${tsv} > ${tsv.baseName}.cleaned.segment.tsv
    """
}
