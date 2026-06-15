nextflow.enable.types = true

// include { Clonality ; Vartype } from '../../types.nf'

// TSV is beta value TSV file
process CLEAN_METHYLATION_TSV {
    input:
    record(
        sample: String,
        tsv: Path
    )

    output:
    record(
        sample: sample,
        tsv: file("${tsv.baseName}.cleaned.meth.tsv"),
    )

    script:
    """
    cat ${tsv} > ${tsv.baseName}.cleaned.meth.tsv
    """
}
