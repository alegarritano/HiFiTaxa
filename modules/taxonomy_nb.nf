// DADA2 R-based Naive-Bayes taxonomy, GTDB-only, two-step (genus + species).
// HiFiTaxa anchors BLCA, Emu, and NB on the same single GTDB release for direct
// comparability. assignTaxonomy with species-in-headers does not scale to full
// GTDB SSU r232 (bootstrap dilutes to Kingdom), so NB uses a genus-level
// reference for the bootstrap step + a species reference for exact-match
// addSpecies (see scripts/dada2_assign_tax.R).
//
// Inputs:
//   asv_seq_fasta : dada2_ASV.fasta from filter_dada2
//   asv_seq       : dada2-ccs_rep_filtered.qza
//   asv_freq      : dada2-ccs_table_filtered.qza
//   genus_db      : 6-rank (Kingdom..Genus) DADA2 reference for assignTaxonomy
//   species_db    : ">acc Genus species" DADA2 reference for addSpecies
//                   (both built by bin/build_gtdb_dada2_db.sh from the same
//                   BLCA-parsed GTDB that BLCA + Emu consume)
//   assign_script : bundled scripts/dada2_assign_tax.R
//
// Outputs:
//   best_tax.qza, best_taxonomy.tsv, best_taxonomy_withDB.tsv,
//   best_tax_merged_freq_tax.tsv, gtdb_nb.tsv

process nb_classify {
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "${params.outdir}/taxonomy_nb", pattern: 'best_tax*', mode: params.publish_dir_mode
    publishDir "${params.outdir}/taxonomy_nb",                       mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_seq_fasta
    path asv_seq
    path asv_freq
    path genus_db
    path species_db
    path assign_script

    output:
    path "best_tax.qza",                 emit: best_nb_tax_qza
    path "best_taxonomy.tsv",            emit: best_nb_tax
    path "best_taxonomy_withDB.tsv"
    path "best_tax_merged_freq_tax.tsv", emit: best_nb_tax_tsv
    path "gtdb_nb.tsv",                  emit: gtdb_nb_tsv

    script:
    """
    # QIIME2 image \$HOME is read-only; redirect caches DADA2 / R / QIIME2 use.
    export HOME="\$PWD/.dada2_home"
    mkdir -p "\$HOME"

    Rscript --vanilla ${assign_script} ${asv_seq_fasta} ${task.cpus} \\
        ${genus_db} ${species_db} ${params.nb_min_bootstrap} ${params.random_seed}

    qiime feature-table transpose --i-table ${asv_freq} \\
        --o-transposed-feature-table transposed-asv.qza

    qiime tools import --type "FeatureData[Taxonomy]" \\
        --input-format "TSVTaxonomyFormat" \\
        --input-path best_taxonomy.tsv --output-path best_tax.qza

    qiime metadata tabulate --m-input-file ${asv_seq} \\
        --m-input-file best_tax.qza \\
        --m-input-file transposed-asv.qza \\
        --o-visualization merged_freq_tax.qzv

    qiime tools export --input-path merged_freq_tax.qzv \\
        --output-path merged_freq_tax_tsv

    mv merged_freq_tax_tsv/metadata.tsv best_tax_merged_freq_tax.tsv
    """
}

// NB on a standalone ASV/sequence FASTA (the -entry taxonomy_only path).
// Same two-step classifier as nb_classify (genus assignTaxonomy + exact-match
// addSpecies) but WITHOUT the QIIME2 frequency-table merge: a standalone FASTA
// carries no per-sample ASV frequency table, so it just writes the per-sequence
// taxonomy. Outputs best_taxonomy.tsv, best_taxonomy_withDB.tsv, gtdb_nb.tsv.
process nb_classify_fasta {
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "${params.outdir}/taxonomy_nb", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_seq_fasta
    path genus_db
    path species_db
    path assign_script

    output:
    path "best_taxonomy.tsv",        emit: best_nb_tax
    path "best_taxonomy_withDB.tsv"
    path "gtdb_nb.tsv",              emit: gtdb_nb_tsv

    script:
    """
    # QIIME2 image \$HOME is read-only; redirect caches DADA2 / R use.
    export HOME="\$PWD/.dada2_home"
    mkdir -p "\$HOME"

    Rscript --vanilla ${assign_script} ${asv_seq_fasta} ${task.cpus} \\
        ${genus_db} ${species_db} ${params.nb_min_bootstrap} ${params.random_seed}
    """
}

// Single-step DADA2 NB for ITS (marker==ITS). Unlike the two-step GTDB path above
// (genus assignTaxonomy + exact-match addSpecies), ITS runs ONE assignTaxonomy()
// straight to Species (7 ranks: Kingdom..Species) against a UNITE reference that
// already carries the full species lineage in its headers — NO addSpecies, which
// collapses to ~0 species on the extracted ITS region. The R companion
// (scripts/dada2_assign_tax_singlestep.R) writes both outputs to the work dir:
//   best_taxonomy.tsv : QIIME2 TSVTaxonomyFormat (Feature ID, Taxon, Confidence)
//   nb_singlestep.tsv : per-ASV 7-rank table (Feature ID + Kingdom..Species)
//
// Inputs:
//   asv_seq_fasta : dada2_ASV.fasta from filter_dada2 (or a standalone query FASTA)
//   singlestep_db : 7-rank UNITE DADA2 reference (params.unite_dada2_singlestep_db,
//                   built by bin/build_unite_dada2_db.sh)
//   assign_script : bundled scripts/dada2_assign_tax_singlestep.R
process nb_classify_singlestep {
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "${params.outdir}/taxonomy_nb", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_seq_fasta
    path singlestep_db
    path assign_script

    output:
    path "best_taxonomy.tsv", emit: best_nb_tax
    path "nb_singlestep.tsv", emit: nb_tsv

    script:
    """
    # QIIME2 image \$HOME is read-only; redirect caches DADA2 / R use.
    export HOME="\$PWD/.dada2_home"
    mkdir -p "\$HOME"

    Rscript --vanilla ${assign_script} ${asv_seq_fasta} ${task.cpus} \\
        ${singlestep_db} ${params.nb_min_bootstrap} ${params.random_seed}
    """
}
