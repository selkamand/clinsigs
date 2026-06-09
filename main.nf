nextflow.enable.dsl = 2
nextflow.enable.types = true

params {

    // Samplesheet
    input: Path
    outdir: String = "results"
}

include { BCFTOOLS_PASS ; BCFTOOLS_NORMALISE_INDELS_AND_SPLIT_MULTIALLELICS ; SPLIT_PURPLE_SNVS_BY_CLONALITY ; PURPLE_SNV_VCF_TO_TSV } from './modules/local/bcftools.nf'

record InputRecord {
    sample: String
    snv: Path?
    sv: Path?
    cnv: Path?
    rnamut: Path?
    meth: Path?
    teal: Path?
}

enum Vartype {
    SNV,
    SV,
    CNV,
    RNAmut,
    METH,
}
enum Clonality {
    all,
    clonal,
    subclonal,
}

record LongifiedRecord {
    sample: String
    vartype: Vartype
    clonality: Clonality
    file: Path
}

workflow {

    main:
    ch_samples = channel.of(params.input)
        .flatMap { csv -> csv.splitCsv(header: true) }
        .map { row ->

            // Pull required columns out of the row
            def sample = row.sample as String
            def snv = row.snv != null ? file(row.snv) : null
            def sv = row.sv != null ? file(row.sv) : null
            //def cnv = row.cnv != null ? file(row.cnv) : null
            def rnamut = row.rnamut != null ? file(row.rnamut) : null
            def meth = row.meth != null ? file(row.meth) : null
            def teal = row.teal != null ? file(row.teal) : null
            def bam = row.bam != null ? file(row.bam) : null

            if (!sample) {
                error("Missing sample column in manifest row: ${row}")
            }

            record(sample: sample, snv: snv, sv: sv, cnv: cnv, bam: bam, rnamut: rnamut, meth: meth, teal: teal) as InputRecord
        }

    publish:
    output1 = ch_samples
}

output {
    output1 {
        path { r -> "${params.outdir}/${r.sample}/files/" }
        mode 'copy'
    }
}
