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
    tag sample
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
    samtools quickcheck -v tumor.bam
    test -s tumor.bam.bai

    samtools idxstats tumor.bam > bam.index.tsv
    awk 'NR==FNR { if (\$1 ~ /^chr([0-9]+|X|Y)\$/) ref[\$1]=\$2; next }
         \$1 in ref { if (\$2 != ref[\$1]) exit 1; delete ref[\$1] }
         END { if (length(ref)) exit 1 }' GRCh38/GRCh38_noAlt.fa.fai bam.index.tsv

    python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py -s '${sample}' -t ${task.cpus} \\
        --bam tumor.bam --ref GRCh38 --cnv_bed '${seed}' \\
        --run_AA --run_AC --AA_solver clarabel --no_cstats --no_repo_check -o results
    report_context.py '${sample}'
    """
}

workflow {
    if (!params.input || !params.purple || !params.mask || !params.reference || !params.outdir)
        error 'Set --input, --purple, --mask, --reference and --outdir'

    manifest = file(params.input, checkIfExists: true)
    alignments = channel.of(manifest).splitCsv(header: true).map { row ->
        if (!row.sample || !(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/))
            error "Invalid sample ID: ${row.sample}"
        if (!row.bam?.endsWith('.bam') || !row.bai?.endsWith('.bai'))
            error "BAM and BAI are required for ${row.sample}"
        tuple(row.sample, row.bam, row.bai)
    }

    PREPARE_SEEDS(manifest, file(params.purple, checkIfExists: true),
                  file(params.mask, checkIfExists: true))
    seedBeds = PREPARE_SEEDS.out.ready.flatten().map { seed ->
        tuple(seed.name - '_AA_CNV_SEEDS.bed', seed)
    }
    samples = alignments.join(seedBeds, failOnDuplicate: true).map { sample, bam, bai, seed ->
        tuple(sample, file(bam, checkIfExists: true), file(bai, checkIfExists: true), seed)
    }
    AMPLICON_SUITE(samples, file(params.reference, checkIfExists: true))
}
