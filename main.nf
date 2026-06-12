nextflow.enable.dsl = 2
nextflow.enable.types = true

params {

    // Samplesheet
    input: Path
    reference_genome_fasta: Path
    outdir: String = "results"
}

include {
    BCFTOOLS_PASS ;
    BCFTOOLS_PASS as BCFTOOLS_PASS_SVS ;
    BCFTOOLS_NORMALISE_INDELS_AND_SPLIT_MULTIALLELICS ;
    SPLIT_PURPLE_SNVS_BY_CLONALITY ;
    PURPLE_SNV_VCF_TO_TSV ;
    PURPLE_SV_VCF_TO_BEDPE
} from './modules/local/bcftools.nf'

include {
    Vartype ;
    Clonality
} from './types.nf'

record InputRecord {
    sample: String
    snv: Path?
    sv: Path?
    cnv: Path?
    rnamut: Path?
    meth: Path?
    teal: Path?
}

record LongifiedRecord {
    sample: String
    vartype: Vartype
    clonality: Clonality
    file: Path
}

record Reference {
    reference_genome_fasta: Path
    reference_genome_fai: Path
}

record ClonalityRecordWide {
    sample: String
    all: Path
    subclonal: Path
    clonal: Path
}
record ClonalityRecordLong {
    sample: String
    clonality: Clonality
    vcf: Path
}


workflow {

    main:

    // Create reference record
    def reference_genome_fai = file("${params.reference_genome_fasta}.fai")
    if (!reference_genome_fai.exists()) {
        error("Failed to find fai index for reference genome at ${reference_genome_fai}")
    }
    def ref = record(
        reference_genome_fasta: params.reference_genome_fasta,
        reference_genome_fai: reference_genome_fai,
    )

    // Parse Samplesheet
    ch_samples = channel.of(params.input)
        .flatMap { csv -> csv.splitCsv(header: true) }
        .map { row ->

            // Pull required columns out of the row
            def sample = row.sample as String
            def snv = row.snv == null | row.snv.isBlank() ? null : file(row.snv)
            def sv = row.sv == null | row.sv.isBlank() ? null : file(row.sv)
            def cnv = row.cnv == null | row.cnv.isBlank() ? null : file(row.cnv)
            def rnamut = row.rnamut == null | row.rnamut.isBlank() ? null : file(row.rnamut)
            def meth = row.meth == null | row.meth.isBlank() ? null : file(row.meth)
            def teal = row.teal == null | row.teal.isBlank() ? null : file(row.teal)
            def bam = row.bam == null | row.bam.isBlank() ? null : file(row.bam)

            if (!sample) {
                error("Missing sample column in manifest row: ${row}")
            }

            record(sample: sample, snv: snv, sv: sv, cnv: cnv, bam: bam, rnamut: rnamut, meth: meth, teal: teal)
        }


    // Split 
    ch_split = FILTER_NORMALISE_SPLITCLONES(ch_samples, ref)

    // Create a debug folder for our TSV descriptions of channels
    def debugfolder = file("${workflow.launchDir}/debug/")
    if (!debugfolder.exists()) {
        debugfolder.mkdirs()
    }
    def longified_files = file("${debugfolder}/inputs_longified.tsv")

    longified_files.text = "sample\tvartype\tclonality\tfile"
    ch_split.subscribe { r -> longified_files.append("\n${r.sample}\t${r.vartype}\t${r.clonality}\t${r.file.name}") }

    publish:
    preprocessed = ch_split
}

workflow FILTER_NORMALISE_SPLITCLONES {
    take:
    ch_samples: Channel<InputRecord>
    ref: Reference

    main:
    // Create a long sample X vcf channel for SNVs SVs and RNAmut VCFs to support filtering/normalisation
    ch_snv_vcfs = ch_samples
        .filter { r -> r.snv != null }
        .map { r -> record(sample: r.sample, vcf: r.snv) }
    ch_sv_vcfs = ch_samples
        .filter { r -> r.sv != null }
        .map { r -> record(sample: r.sample, vcf: r.sv) }
    ch_rnamut_vcfs = ch_samples
        .filter { r -> r.rnamut != null }
        .map { r -> record(sample: r.sample, vcf: r.rnamut) }

    // ch_cbv = ch_samples.map { r -> record(sample: r.sample, vcf: r.rnamut, vartype: Vartype.RNAmut) }

    // Filter and clonality split SNVS
    ch_snv_pass = BCFTOOLS_PASS(ch_snv_vcfs)
    ch_snv_norm = BCFTOOLS_NORMALISE_INDELS_AND_SPLIT_MULTIALLELICS(ch_snv_pass, ref)
    ch_snv_split = SPLIT_PURPLE_SNVS_BY_CLONALITY(ch_snv_norm)
    ch_snv_split_long = LONGIFY_CLONALITY(ch_snv_split)
    ch_snv_tsvs = PURPLE_SNV_VCF_TO_TSV(ch_snv_split_long)
    ch_snv_tsvs_with_vartype = ch_snv_tsvs.map { r ->
        record(
            sample: r.sample,
            vartype: Vartype.SNV,
            clonality: r.clonality,
            file: r.tsv,
        )
    }

    // Filter for PASS SVs and create breakpoint bedpe files + Breakend level TSVS
    ch_sv_pass = BCFTOOLS_PASS_SVS(ch_sv_vcfs)
    ch_sv_bedpe = PURPLE_SV_VCF_TO_BEDPE(ch_sv_pass)
    ch_sv_bedpe_with_vartype = ch_sv_bedpe.map { r ->
        record(
            sample: r.sample,
            vartype: Vartype.SV_BEDPE,
            clonality: Clonality.all,
            file: r.bedpe,
        )
    }

    // Pass along copynumber segment files (since we can't split by clonality)
    ch_cnv_segments_with_vartype = ch_samples
        .filter { r -> r.cnv != null }
        .map { r ->
            record(sample: r.sample, vartype: Vartype.CNV, clonality: Clonality.all, file: r.cnv)
        }

    // Pass along methylation files (since we can't split by clonality)
    ch_meth_with_vartype = ch_samples
        .filter { r -> r.meth != null }
        .map { r ->
            record(sample: r.sample, vartype: Vartype.METH, clonality: Clonality.all, file: r.meth)
        }

    // Combine different variant types into one very long channel of longified records
    ch_long = ch_snv_tsvs_with_vartype
        .mix(ch_sv_bedpe_with_vartype)
        .mix(ch_cnv_segments_with_vartype)
        .mix(ch_meth_with_vartype)

    emit:
    ch_long: Channel<LongifiedRecord>
}


workflow LONGIFY_CLONALITY {
    take:
    records: Channel<ClonalityRecordWide>

    main:
    longified = records.flatMap { r ->
        [
            record(sample: r.sample, clonality: Clonality.all, vcf: r.all),
            record(sample: r.sample, clonality: Clonality.clonal, vcf: r.clonal),
            record(sample: r.sample, clonality: Clonality.subclonal, vcf: r.subclonal),
        ]
    }

    emit:
    longified: Channel<ClonalityRecordLong>
}


output {
    preprocessed {
        path { r -> "${params.outdir}/${r.sample}/${r.vartype}/${r.clonality}/" }
        mode 'copy'
    }
}
