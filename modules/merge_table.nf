// Final per-run ASV table: merge the DADA2 ASV-by-sample frequency table with a
// classifier's taxonomy (and its bootstrap confidence) into one wide TSV:
//   Feature ID | <count per sample...> | Domain..Species | confidence
// Confidence is a single column for NB (its bootstrap) or one column per rank
// for BLCA (--per-rank-confidence, from blca_taxonomy_confidence.csv). One table
// per ASV-based classifier; Emu/EMITS are read-level and have no ASVs, so they
// are not merged here. merge_asv_freq_taxonomy.py lives in bin/, so Nextflow puts
// it on PATH for every executor (local/conda/docker/singularity).

process merge_asv_table {
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    cpus 1

    input:
    tuple path(asv_freq), path(taxonomy), val(classifier), val(conf_flag)

    output:
    path "${classifier}_asv_table.tsv", emit: table

    script:
    """
    merge_asv_freq_taxonomy.py \\
        --asv-freq ${asv_freq} \\
        --taxonomy ${taxonomy} \\
        ${conf_flag} \\
        --out ${classifier}_asv_table.tsv
    """
}
