nextflow.enable.dsl=2

process PREPARE_SEEDS {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path manifest
    path purple
    path mask
    output:
    path 'seed_audit', emit: audit
    path 'ready/*_AA_CNV_SEEDS.bed', optional: true, emit: ready
    script:
    """
    prepare_seeds.py '${manifest}' '${purple}' '${mask}'
    """
}

process AMPLICON_SUITE {
    tag "${sample}"
    publishDir { "${params.outdir}/samples/${sample}" }, mode: 'copy'
    input:
    tuple val(sample), path('tumor.bam'), path('tumor.bam.bai'), path(seed)
    path reference
    output:
    path 'results/*'
    script:
    """
    mkdir -p results
    tar -xzf '${reference}'
    export AA_DATA_REPO="\$PWD"
    check_bam.sh tumor.bam GRCh38/GRCh38_noAlt.fa.fai

    python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py -s '${sample}' -t ${task.cpus} \\
        --bam tumor.bam --ref GRCh38 --cnv_bed '${seed}' \\
        --run_AA --run_AC --AA_solver clarabel --no_cstats --no_repo_check -o results
    report_context.py '${sample}'
    """
}

// Default AmpliconSuite seeding: CNVkit with the matched normal, then the
// wrapper's own cnv_prefilter and amplified_intervals filters. No other
// AmpliconSuite or AA options are changed from their defaults.
process AMPLICON_SUITE_DEFAULT {
    tag "${sample}"
    publishDir { "${params.outdir}/samples/${sample}" }, mode: 'copy'
    input:
    tuple val(sample), path('tumor.bam'), path('tumor.bam.bai'), path('normal.bam'), path('normal.bam.bai')
    path reference
    output:
    path 'results/*'
    script:
    """
    mkdir -p results
    tar -xzf '${reference}'
    export AA_DATA_REPO="\$PWD"
    check_bam.sh tumor.bam GRCh38/GRCh38_noAlt.fa.fai
    check_bam.sh normal.bam GRCh38/GRCh38_noAlt.fa.fai

    python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py -s '${sample}' -t ${task.cpus} \\
        --bam tumor.bam --normal_bam normal.bam --ref GRCh38 --run_AA --run_AC -o results
    report_context.py '${sample}'
    """
}

def validSample(row) {
    if (!row.sample || !(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
        error "Invalid sample ID: ${row.sample}"
    if (!row.bam?.endsWith('.bam') || !row.bai?.endsWith('.bai'))
        error "BAM and BAI are required for ${row.sample}"
    return row
}

workflow {
    if (!(params.mode in ['purple', 'default']))
        error "Unknown --mode ${params.mode}; use 'purple' or 'default'"
    if (!params.input || !params.reference || !params.outdir)
        error 'Set --input, --reference and --outdir'

    manifest = file(params.input, checkIfExists: true)
    rows = channel.of(manifest).splitCsv(header: true).map { row -> validSample(row) }
    reference = file(params.reference, checkIfExists: true)

    if (params.mode == 'default') {
        samples = rows.map { row ->
            if (!row.normal_bam?.endsWith('.bam') || !row.normal_bai?.endsWith('.bai'))
                error "normal_bam and normal_bai are required for ${row.sample} in default mode"
            tuple(row.sample, file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true),
                  file(row.normal_bam, checkIfExists: true), file(row.normal_bai, checkIfExists: true))
        }
        AMPLICON_SUITE_DEFAULT(samples, reference)
    } else {
        if (!params.purple || !params.mask)
            error 'Set --purple and --mask in purple mode'
        alignments = rows.map { row -> tuple(row.sample, row.bam, row.bai) }
        PREPARE_SEEDS(manifest, file(params.purple, checkIfExists: true),
                      file(params.mask, checkIfExists: true))
        seedBeds = PREPARE_SEEDS.out.ready.flatten().map { seed ->
            tuple(seed.name - '_AA_CNV_SEEDS.bed', seed)
        }
        samples = alignments.join(seedBeds, failOnDuplicate: true).map { sample, bam, bai, seed ->
            tuple(sample, file(bam, checkIfExists: true), file(bai, checkIfExists: true), seed)
        }
        AMPLICON_SUITE(samples, reference)
    }
}
