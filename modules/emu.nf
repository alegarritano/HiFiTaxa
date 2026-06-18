// Emu: EM-based species-level taxonomic profiling for full-length 16S reads
// (Curry et al. 2022, Nat. Methods). Operates on (optionally cutadapt-trimmed)
// per-sample FASTQ files and hits a custom GTDB-formatted Emu database built
// by bin/build_gtdb_emu_db.sh.
//
// Two processes:
//   emu_classify  : per-sample `emu abundance` -> <sample>_rel-abundance.tsv
//   emu_collate   : merges per-sample TSVs into species- and genus-level tables
//
// Unlike BLCA, Emu does NOT need DADA2 ASVs as input; it maps raw reads against
// the reference with minimap2 and runs an EM loop. The `--type` flag forwards
// to minimap2's preset; default is map-hifi (Emu 3.6+ ships a HiFi preset for
// PacBio CCS). Other valid options: map-ont, map-pb, sr, lr:hq, splice:hq.

process emu_classify {
    conda (params.enable_conda ? "${projectDir}/envs/emu.yml" : null)
    container params.emu_container
    publishDir "${params.outdir}/taxonomy_emu/per_sample", mode: params.publish_dir_mode
    cpus { params.emu_threads as int }
    label 'emu'

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    path("${sampleID}_rel-abundance.tsv"),      emit: abundance
    path("${sampleID}_read-assignment-distributions.tsv"), optional: true
    path("${sampleID}_emu_alignments.sam"),     optional: true

    script:
    def min_abund_opt = params.emu_min_abundance != null ? "--min-abundance ${params.emu_min_abundance}" : ""
    """
    emu abundance \\
        --db ${params.emu_db_dir} \\
        --type ${params.emu_type} \\
        --threads ${task.cpus} \\
        --keep-counts \\
        ${min_abund_opt} \\
        --output-dir . \\
        --output-basename ${sampleID} \\
        ${sampleFASTQ}
    """
}

process emu_collate {
    conda (params.enable_conda ? "${projectDir}/envs/emu.yml" : null)
    container params.emu_container
    publishDir "${params.outdir}/taxonomy_emu", mode: params.publish_dir_mode
    cpus 1
    label 'emu'

    input:
    path "*"                                                // staged per-sample _rel-abundance.tsv

    output:
    path "emu_species_table.tsv", emit: species_table
    path "emu_genus_table.tsv",   emit: genus_table

    script:
    """
    # Emu's combine-outputs takes a directory of per-sample _rel-abundance.tsv.
    mkdir -p stage
    cp *_rel-abundance.tsv stage/

    emu combine-outputs stage species
    emu combine-outputs stage genus

    mv stage/emu-combined-species.tsv emu_species_table.tsv
    mv stage/emu-combined-genus.tsv   emu_genus_table.tsv
    """
}
