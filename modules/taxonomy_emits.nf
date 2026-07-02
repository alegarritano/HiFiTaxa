// EMITS: EM-based read-level taxonomic profiling for full-length fungal ITS
// reads. The fungal-ITS analogue of the Emu 16S branch — used when marker==ITS.
//
// Like Emu, EMITS does NOT need DADA2 ASVs: it maps the (cutadapt-trimmed)
// per-sample reads against the UNITE reference with minimap2 and resolves
// abundances with an EM loop. minimap2 runs with the map-hifi preset and keeps
// multiple secondary hits (--secondary=yes -N 10 -p 0.95) so the EM step has the
// full set of near-best references to apportion each read across; `emits run`
// then collapses the PAF alignments to per-taxon abundance at species and genus
// rank against params.emits_db (db_unite/unite.fasta).
//
// emits_classify : per-sample minimap2 -> emits run -> <sample>_emits_{species,genus}.tsv

process emits_classify {
    conda (params.enable_conda ? "${projectDir}/envs/emits.yml" : null)
    container params.emits_container
    publishDir "${params.outdir}/taxonomy_emits/per_sample", mode: params.publish_dir_mode
    cpus { params.emits_threads as int }
    label 'emits'

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    path("${sampleID}_emits_species.tsv"), emit: species
    path("${sampleID}_emits_genus.tsv"),   emit: genus

    script:
    """
    minimap2 -cx map-hifi --secondary=yes -N 10 -p 0.95 -t ${task.cpus} \\
        ${params.emits_db} ${sampleFASTQ} > ${sampleID}.aln.paf

    # emits writes the FULL abundance table to --output; stdout is only a summary
    # (capturing stdout with '>' truncated the table to the top handful of taxa).
    emits run --input ${sampleID}.aln.paf --preset pacbio-hifi --rank species \\
        --output ${sampleID}_emits_species.tsv
    emits run --input ${sampleID}.aln.paf --preset pacbio-hifi --rank genus \\
        --output ${sampleID}_emits_genus.tsv
    """
}
