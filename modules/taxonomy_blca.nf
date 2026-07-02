// BLCA taxonomy against a GTDB BLCA-formatted database.
// The ASV rep-seqs FASTA is split into chunks upstream (splitFasta); each chunk
// is classified in parallel here, then merged + parsed into tidy tables.
//
// BLCA is always run with -p 1 (the conda clustalo is built without OpenMP, so
// --threads>1 aborts); parallelism comes from the per-chunk tasks instead.

process blca_classify {
    conda (params.enable_conda ? "${projectDir}/envs/blca.yml" : null)
    cpus 1   // 1 cpu/chunk -> Nextflow runs ~max_cpus chunks in parallel

    input:
    path asv_chunk

    output:
    path "${asv_chunk.baseName}.blca.out", emit: blca_out

    script:
    // blca_main.py and parse_blca.py live in bin/ -> Nextflow puts them on PATH
    // for every executor (local/conda/docker/singularity).
    """
    blca_main.py \\
        -i ${asv_chunk} \\
        -r ${params.blca_tax} \\
        -q ${params.blca_db} \\
        --iset ${params.blca_minid} \\
        -n ${params.blca_nper} \\
        --seed ${params.random_seed} \\
        -p 1 \\
        -o ${asv_chunk.baseName}.blca.out
    """
}

process merge_blca {
    conda (params.enable_conda ? "${projectDir}/envs/blca.yml" : null)
    publishDir "${params.outdir}/taxonomy_blca", mode: params.publish_dir_mode
    cpus 1

    input:
    path blca_outs

    output:
    path "ASV_blca.out", emit: blca_out
    path "blca_taxonomy_table.csv", emit: tax_table
    path "blca_taxonomy_confidence.csv", emit: tax_conf

    script:
    """
    cat ${blca_outs} > ASV_blca.out
    parse_blca.py ASV_blca.out
    """
}
