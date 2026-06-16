nextflow.enable.types = true
include { Clonality ; Vartype ; Scheme } from '../../types.nf'


// Every tally process must output a tsv with the following columns:
// 1. Feature
// 2. Count

process SBS96 {
    input:
    record(
        sample: String,
        clonality: Clonality,
        vartype: Vartype,
        file: Path,
        scheme: Scheme
    )
    record(
        reference_genome_fasta: Path,
        reference_genome_fai: Path
    )

    output:
    record(
        sample: sample,
        tally: file("${sample}.tally.${clonality}.${scheme}.tsv"),
    )

    script:
    """
    bioprep tally --class SBS96 ${file} > ${sample}.tally.${clonality}.${scheme}.tsv
    """
}


process SBS6 {
    input:
    record(
        sample: String,
        clonality: Clonality,
        vartype: Vartype,
        file: Path,
        scheme: Scheme
    )
    record(
        reference_genome_fasta: Path,
        reference_genome_fai: Path
    )

    output:
    record(
        sample: sample,
        tally: file("${sample}.tally.${clonality}.${scheme}.tsv"),
    )

    script:
    """
    bioprep tally --class SBS6 ${file} > ${sample}.tally.${clonality}.${scheme}.tsv
    """
}
