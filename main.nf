/*
===============================================================================
 HiFiTaxa
 PacBio HiFi 16S amplicon pipeline:
   QC -> primer trim -> ( QIIME2/DADA2 denoise -> ASV filter -> BLCA(GTDB) )
                        and / or
                      ( Emu (GTDB) species-level profiling on raw reads )

 Classifier choice is controlled by --classifier {blca|emu|nb|all}. BLCA
 stays the default; Emu is the EM-based species profiler (Curry et al. 2022,
 Nat. Methods); nb is DADA2's Naive-Bayes assignTaxonomy at genus level (Wang
 et al. 2007 AEM) plus exact-match addSpecies, running on DADA2 ASVs.

 Denoising logic adapted from PacBio HiFi-16S-workflow (DADA2 via QIIME2);
 BLCA taxonomy runs against a GTDB BLCA-formatted database; Emu runs against a
 GTDB-formatted Emu database built by bin/build_gtdb_emu_db.sh.
===============================================================================
*/
nextflow.enable.dsl=2

include { write_log } from './modules/utils'
include { QC_fastq; QC_fastq_post_trim; cutadapt; cutadapt_stats; length_filter } from './modules/qc'
include { prepare_qiime2_manifest; prepare_qiime2_manifest_skip_cutadapt; merge_manifest } from './modules/qc'
include { import_qiime2; demux_summarize } from './modules/qc'
include { learn_error; dada2_denoise; dada2_denoise_with_error_model } from './modules/dada2'
include { mergeASV; filter_dada2; dada2_qc } from './modules/dada2'
include { blca_classify; merge_blca } from './modules/taxonomy_blca'
include { emu_classify; emu_collate } from './modules/emu'
include { nb_classify; nb_classify_fasta; nb_classify_singlestep } from './modules/taxonomy_nb'
// --- marker==ITS read-prep + taxonomy modules ---
include { itsx_extract } from './modules/itsx'
include { emits_classify } from './modules/taxonomy_emits'

// ---- classifier selection ------------------------------------------------
// --classifier accepts a single name, a comma-separated list, or the shorthands
// 'all' (== blca,emu,nb) and 'both' (== blca,emu, kept for back-compat).
// The QIIME2 Naive-Bayes branch is canonically 'qiime2_nb'; 'nb', 'dada2_nb',
// 'qiime2-nb', 'dada2-nb' are accepted aliases (internally normalized to 'nb').
def VALID_CLASSIFIERS = ['blca', 'emu', 'nb']
def NB_ALIASES        = ['nb', 'dada2_nb', 'dada2-nb', 'qiime2_nb', 'qiime2-nb']
def _raw = params.classifier.toString().toLowerCase().trim()
def selected
if (_raw == 'all')       selected = ['blca', 'emu', 'nb']
else if (_raw == 'both') selected = ['blca', 'emu']
else                     selected = _raw.split(',').collect { it.trim() }
selected = selected.collect { it in NB_ALIASES ? 'nb' : it }
selected.each {
    if (!(it in VALID_CLASSIFIERS)) {
        exit 1, "ERROR: unknown classifier '${it}'; must be one of ${VALID_CLASSIFIERS} " +
                "(or NB aliases ${NB_ALIASES}), or shorthand 'all' / 'both'."
    }
}
def run_blca = 'blca' in selected
def run_emu  = 'emu'  in selected
def run_nb   = 'nb'   in selected

// ---- dynamic parameters --------------------------------------------------
if (params.input) {
    n_sample = file(params.input).countLines() - 1
    if (n_sample == 1) {
        dynamic_min_asv_totalfreq = 0
        dynamic_min_asv_sample = 0
        println("Only 1 sample. min_asv_sample and min_asv_totalfreq set to 0.")
    } else {
        dynamic_min_asv_totalfreq = params.min_asv_totalfreq
        dynamic_min_asv_sample = params.min_asv_sample
    }
} else {
    n_sample = 0
    dynamic_min_asv_totalfreq = 0
    dynamic_min_asv_sample = 0
}

if (params.skip_primer_trim) {
    dynamic_forward_primer = 'none'; dynamic_reverse_primer = 'none'; trim_cutadapt = "No"
} else {
    dynamic_forward_primer = params.forward_primer; dynamic_reverse_primer = params.reverse_primer; trim_cutadapt = "Yes"
}

// ---- marker awareness ----------------------------------------------------
// marker drives primers, length window, reference DB, read-prep, NB design and
// the read-level EM classifier. Validate up front; default '16S' keeps the
// classic GTDB behaviour, 'ITS' switches the whole fungal path on.
def MARKER = params.marker.toString().toUpperCase().trim()
if (!(MARKER in ['16S', 'ITS'])) {
    exit 1, "ERROR: unknown --marker '${params.marker}'; must be '16S' or 'ITS'."
}
def is_its = (MARKER == 'ITS')
def reference_name = is_its ? 'UNITE (db_unite)' : 'GTDB (db)'
def read_em_name   = is_its ? 'EMITS' : 'Emu'
def nb_design_name = is_its ? 'single-step (7-rank assignTaxonomy, no addSpecies)'
                            : 'two-step (genus assignTaxonomy + exact-match addSpecies)'
def read_prep_name = is_its ? 'cutadapt -> itsxrust ITS extraction -> DADA2'
                            : 'cutadapt -> DADA2'

log_text = """
  HiFiTaxa pipeline
  =======================
  Samples in TSV:            $n_sample
  Marker:                    $MARKER
  Reference:                 $reference_name
  Read prep:                 $read_prep_name
  Filter reads above Q:      $params.filterQ
  Trim primers (cutadapt):   $trim_cutadapt
  Forward primer:            $params.forward_primer
  Reverse primer:            $params.reverse_primer
  Classifier(s):             ${selected.join(',')}
  NB design:                 $nb_design_name
  Read-level EM classifier:  $read_em_name
  --- DADA2 (needed for BLCA + NB) ---
  DADA2 min/max len:         $params.min_len / $params.max_len
  DADA2 maxEE / minQ:        $params.max_ee / $params.minQ
  DADA2 pooling:             $params.pooling_method
  Min ASV total freq:        $dynamic_min_asv_totalfreq
  Min ASV samples:           $dynamic_min_asv_sample
  --- BLCA branch ---
  GTDB BLCA db:              $params.blca_db
  GTDB BLCA taxonomy:        $params.blca_tax
  BLCA chunk size (ASVs):    $params.blca_chunk_size
  --- Emu branch ---
  Emu DB dir:                $params.emu_db_dir
  Emu minimap2 preset:       $params.emu_type
  Emu threads/sample:        $params.emu_threads
  Emu read length filter:    $params.min_len-$params.max_len bp (matches DADA2)
  Emu min-abundance:         ${params.emu_min_abundance != null ? params.emu_min_abundance : 'default'}
  --- NB branch ---
  NB DADA2 two-step (genus assignTaxonomy + exact-match addSpecies, bootstrap=$params.nb_min_bootstrap):
    GTDB genus ref   : $params.gtdb_dada2_genus_db
    GTDB species ref : $params.gtdb_dada2_species_db
"""

workflow {
    if (!params.input)    { exit 1, "ERROR: --input <samples.tsv> is required" }
    if (!params.metadata) { exit 1, "ERROR: --metadata <metadata.tsv> is required" }

    log.info(log_text)
    write_log(log_text)

    sample_file = channel.fromPath(params.input)
        .splitCsv(header: ['sample', 'fastq'], skip: 1, sep: "\t")
        // absolute paths used as-is; relative paths (e.g. the bundled example
        // manifest) resolve against the pipeline directory so they work on any clone.
        .map{ row -> tuple(row.sample, file(row.fastq.startsWith('/') ? row.fastq : "${projectDir}/${row.fastq}")) }
    metadata_file = channel.fromPath(params.metadata)

    // ---- QC + primer trim (always; both classifier branches consume the trimmed reads) ----
    QC_fastq(sample_file)

    // `reads_for_emu` holds per-sample (sampleID, FASTQ) tuples ready for Emu.
    // For BLCA, we additionally need a QIIME2 manifest + import; that work runs
    // only when the BLCA branch is selected.
    def reads_for_emu
    def qiime2_manifest
    def filtered_fastq_files

    if (params.skip_primer_trim) {
        reads_for_emu        = QC_fastq.out.filtered_fastq
        // marker==ITS still gets the itsxrust ITS extraction before import (only the
        // cutadapt primer trim is skipped); 16S is unchanged. EMITS also runs on the
        // itsxrust-extracted reads.
        def reads_for_import_skip
        if (is_its) {
            itsx_extract(QC_fastq.out.filtered_fastq)
            reads_for_import_skip = itsx_extract.out.fastq
            filtered_fastq_files  = itsx_extract.out.fastq.map { sid, fq -> fq }
            reads_for_emu         = itsx_extract.out.fastq   // EMITS on itsxrust-extracted reads (default)
        } else {
            reads_for_import_skip = QC_fastq.out.filtered_fastq
            filtered_fastq_files  = QC_fastq.out.filtered_fastq_files
        }
        if (run_blca) {
            prepare_qiime2_manifest_skip_cutadapt(reads_for_import_skip.collect(), metadata_file)
            qiime2_manifest = prepare_qiime2_manifest_skip_cutadapt.out.sample_trimmed_file.flatten()
            import_qiime2(qiime2_manifest, filtered_fastq_files.collect())
        }
    } else {
        cutadapt(QC_fastq.out.filtered_fastq, dynamic_forward_primer, dynamic_reverse_primer)
        QC_fastq_post_trim(cutadapt.out.cutadapt_fastq)
        // primer-removal stats, printed to the log right after cutadapt
        cutadapt_stats(cutadapt.out.summary_tocollect.collect())
        cutadapt_stats.out.stats.splitText().view { "[primer-removal] " + it.trim() }
        // Emu (16S) profiles the cutadapt-trimmed reads directly. EMITS (ITS) profiles
        // the itsxrust ITS-extracted reads (reads_for_emu is re-pointed to itsx_extract
        // inside the is_its block below), so the read-level fungal profiler classifies the
        // same ITS span as the ASV path — the configuration used for the reported results.
        reads_for_emu        = cutadapt.out.cutadapt_fastq

        // marker==ITS read prep: pull the full ITS span out of the trimmed reads with
        // itsxrust BEFORE the QIIME2 import, so DADA2 denoises ITS-only sequences, and
        // feed the same itsxrust reads to EMITS. 16S keeps the cutadapt -> import path.
        def reads_for_import
        if (is_its) {
            itsx_extract(cutadapt.out.cutadapt_fastq)
            reads_for_import     = itsx_extract.out.fastq
            filtered_fastq_files = itsx_extract.out.fastq.map { sid, fq -> fq }
            reads_for_emu        = itsx_extract.out.fastq   // EMITS classifies the itsxrust ITS-extracted reads (default)
        } else {
            reads_for_import     = cutadapt.out.cutadapt_fastq
            filtered_fastq_files = cutadapt.out.cutadapt_fastq_files
        }
        if (run_blca) {
            prepare_qiime2_manifest(reads_for_import.collect(), metadata_file)
            qiime2_manifest = prepare_qiime2_manifest.out.sample_trimmed_file.flatten()
            import_qiime2(qiime2_manifest, filtered_fastq_files.collect())
        }
    }

    // DADA2 + filter run when EITHER BLCA or NB is requested (both consume ASVs).
    if (run_blca || run_nb) {
        demux_summarize(import_qiime2.out)

        if (params.learn_error_sample) {
            learn_error(params.learn_error_sample, params.learnError_script)
            dada2_denoise_with_error_model(import_qiime2.out, params.dadaCCS_script, params.minQ, learn_error.out.dada2_error_model)
            mergeASV(dada2_denoise_with_error_model.out.asv_seq.collect(),
                     dada2_denoise_with_error_model.out.asv_freq.collect(),
                     dada2_denoise_with_error_model.out.asv_stats.collect())
        } else {
            dada2_denoise(import_qiime2.out, params.dadaCCS_script, params.minQ)
            mergeASV(dada2_denoise.out.asv_seq.collect(),
                     dada2_denoise.out.asv_freq.collect(),
                     dada2_denoise.out.asv_stats.collect())
        }

        filter_dada2(mergeASV.out.asv_freq, mergeASV.out.asv_seq, dynamic_min_asv_totalfreq, dynamic_min_asv_sample)
        dada2_qc(mergeASV.out.asv_stats, filter_dada2.out.asv_freq, metadata_file)
    }

    // ---------------- BLCA branch (ASV chunks -> BLCA/GTDB with bootstrap) ----------------
    if (run_blca) {
        // Chunk size: user value, else AUTO ~= one chunk per usable core.
        asv_chunks = filter_dada2.out.asv_seq_fasta.flatMap { f ->
            int total = f.countFasta()
            int cs = params.blca_chunk_size ? (params.blca_chunk_size as int)
                     : Math.max(1, (int) Math.ceil(total / (params.max_cpus as double)))
            f.splitFasta(by: cs, file: true)
        }
        blca_classify(asv_chunks)
        merge_blca(blca_classify.out.blca_out.collect())
    }

    // ---------------- NB branch ------------------------------------------------
    // 16S: DADA2 two-step (genus assignTaxonomy bootstrap + exact-match addSpecies)
    //      against the same full GTDB release as BLCA/Emu (scripts/dada2_assign_tax.R).
    // ITS: DADA2 single-step (one 7-rank assignTaxonomy straight to Species, NO
    //      addSpecies) against the UNITE reference (scripts/dada2_assign_tax_singlestep.R).
    if (run_nb) {
        if (is_its) {
            singlestep_ch = channel.fromPath(params.unite_dada2_singlestep_db,   checkIfExists: true)
            r_script_ss   = channel.fromPath(params.dadaAssignSinglestep_script, checkIfExists: true)
            nb_classify_singlestep(filter_dada2.out.asv_seq_fasta,
                                   singlestep_ch, r_script_ss)
        } else {
            genus_ch   = channel.fromPath(params.gtdb_dada2_genus_db,   checkIfExists: true)
            species_ch = channel.fromPath(params.gtdb_dada2_species_db, checkIfExists: true)
            r_script   = channel.fromPath(params.dadaAssign_script,     checkIfExists: true)
            nb_classify(filter_dada2.out.asv_seq_fasta,
                        filter_dada2.out.asv_seq,
                        filter_dada2.out.asv_freq,
                        genus_ch, species_ch, r_script)
        }
    }

    // ---------------- Read-level EM branch (length-filtered reads) ----------------
    // Length-filter to the same min_len/max_len window DADA2 enforces.
    //   16S -> Emu   (cutadapt-trimmed reads;        minimap2 + EM vs the GTDB-Emu DB)
    //   ITS -> EMITS (itsxrust ITS-extracted reads;  minimap2 map-hifi + emits run vs the UNITE FASTA)
    if (run_emu) {
        length_filter(reads_for_emu, params.min_len, params.max_len)
        if (is_its) {
            emits_classify(length_filter.out.filtered)
        } else {
            emu_classify(length_filter.out.filtered)
            emu_collate(emu_classify.out.abundance.collect())
        }
    }
}

// Taxonomy-only entry: skip QC/denoise, classify an existing ASV/sequence fasta.
// Honours --classifier for BLCA, Emu, and/or NB. NB here runs the two-step
// classifier (genus assignTaxonomy + exact-match addSpecies) on the FASTA and
// writes the per-sequence taxonomy; it skips the ASV-by-sample frequency merge
// (a standalone FASTA carries no frequency table).
//   nextflow run main.nf -entry taxonomy_only -profile standard \
//       --asv_fasta ASV.fasta --classifier all
workflow taxonomy_only {
    if (!params.asv_fasta) exit 1, "ERROR: -entry taxonomy_only requires --asv_fasta <ASV.fasta>"

    // parse --classifier; 'all' = blca+emu+nb, 'both' = blca+emu (legacy)
    def _craw = params.classifier.toString().toLowerCase().trim()
    def sel
    if (_craw == 'all')       sel = ['blca', 'emu', 'nb']
    else if (_craw == 'both') sel = ['blca', 'emu']
    else                      sel = _craw.split(',').collect { it.trim() }
    def nb_aliases = ['nb', 'dada2_nb', 'dada2-nb', 'qiime2_nb', 'qiime2-nb']
    sel = sel.collect { it in nb_aliases ? 'nb' : it }
    def to_blca = 'blca' in sel
    def to_emu  = 'emu'  in sel
    def to_nb   = 'nb'   in sel
    if (!to_blca && !to_emu && !to_nb) {
        to_blca = true
        log.info("No blca/emu/nb in --classifier; defaulting to BLCA.")
    }

    if (to_blca) {
        log.info("Taxonomy-only: BLCA on ${params.asv_fasta} vs ${params.blca_db}")
        asv_chunks = channel.fromPath(params.asv_fasta, checkIfExists: true).flatMap { f ->
            int total = f.countFasta()
            int cs = params.blca_chunk_size ? (params.blca_chunk_size as int)
                     : Math.max(1, (int) Math.ceil(total / (params.max_cpus as double)))
            f.splitFasta(by: cs, file: true)
        }
        blca_classify(asv_chunks)
        merge_blca(blca_classify.out.blca_out.collect())
    }

    if (to_emu) {
        log.info("Taxonomy-only: Emu on ${params.asv_fasta} vs ${params.emu_db_dir}")
        emu_in = channel.fromPath(params.asv_fasta, checkIfExists: true).map { f -> tuple(f.baseName, f) }
        emu_classify(emu_in)
        emu_collate(emu_classify.out.abundance.collect())
    }

    if (to_nb) {
        log.info("Taxonomy-only: NB on ${params.asv_fasta} (genus assignTaxonomy + addSpecies)")
        genus_ch   = channel.fromPath(params.gtdb_dada2_genus_db,   checkIfExists: true)
        species_ch = channel.fromPath(params.gtdb_dada2_species_db, checkIfExists: true)
        r_script   = channel.fromPath(params.dadaAssign_script,     checkIfExists: true)
        nb_classify_fasta(channel.fromPath(params.asv_fasta, checkIfExists: true),
                          genus_ch, species_ch, r_script)
    }
}
