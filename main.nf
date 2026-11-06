/*
===============================================================================
 HiFiTaxa
 PacBio HiFi amplicon pipeline for the 16S rRNA gene (default) and fungal ITS:
   16S: QC -> cutadapt primer trim -> ( QIIME2/DADA2 denoise -> ASV filter -> BLCA and/or
        DADA2 naive Bayes against GTDB ) and/or ( Emu species-level profiling of the QC reads )
   ITS: QC -> itsxrust ITS extraction -> ( DADA2 denoise -> BLCA and/or single-step naive
        Bayes against UNITE ) and/or ( EMITS profiling of the extracted reads )

 Classifier choice is controlled by --classifier {blca|emu|emits|nb|all}. BLCA is
 the 16S default and EMITS the ITS default; Emu is the EM-based species profiler
 (Curry et al. 2022, Nat. Methods); nb is DADA2's naive Bayes assignTaxonomy
 (Wang et al. 2007 AEM) at genus level plus exact-match addSpecies for 16S, and a
 single-step seven-rank assignTaxonomy for ITS, running on DADA2 ASVs.

 Denoising logic adapted from PacBio HiFi-16S-workflow (DADA2 via QIIME2);
 BLCA taxonomy runs against a GTDB or UNITE BLCA-formatted database; Emu runs
 against a GTDB-formatted Emu database built by bin/build_gtdb_emu_db.sh.
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
include { merge_asv_table } from './modules/merge_table'
// --- marker==ITS read-prep + taxonomy modules ---
include { itsx_extract } from './modules/itsx'
include { emits_classify } from './modules/taxonomy_emits'

// ---- marker awareness (resolved before classifier so defaults can depend on it) ----
// marker drives primers, length window, reference DB, read-prep, NB design and the
// read-level EM classifier. Default '16S' keeps the classic GTDB behaviour; 'ITS'
// switches the whole fungal path on.
def MARKER = params.marker.toString().toUpperCase().trim()
if (!(MARKER in ['16S', 'ITS'])) {
    exit 1, "ERROR: unknown --marker '${params.marker}'; must be '16S' or 'ITS'."
}
def is_its = (MARKER == 'ITS')

// ---- classifier selection ------------------------------------------------
// --classifier accepts a single name, a comma-separated list, or the shorthands
// 'all' (== blca,emu,nb) and 'both' (== blca,emu, kept for back-compat). The
// read-level EM classifier is 'emu' for 16S and 'emits' for ITS (same slot);
// 'emits' is an accepted alias. The QIIME2 Naive-Bayes branch is canonically
// 'qiime2_nb'; 'nb', 'dada2_nb', 'qiime2-nb', 'dada2-nb' are aliases (-> 'nb').
// Default is marker-aware: 16S -> blca (ASV path); ITS -> emits (read-level,
// paired with itsxrust). The DADA2/ASV steps run only when blca or nb is asked.
def VALID_CLASSIFIERS = ['blca', 'emu', 'nb']
def NB_ALIASES        = ['nb', 'dada2_nb', 'dada2-nb', 'qiime2_nb', 'qiime2-nb']
def EMITS_ALIASES     = ['emits']
def _raw = (params.classifier ?: (is_its ? 'emits' : 'blca')).toString().toLowerCase().trim()
def selected
if (_raw == 'all')       selected = ['blca', 'emu', 'nb']
else if (_raw == 'both') selected = ['blca', 'emu']
else                     selected = _raw.split(',').collect { it.trim() }
selected = selected.collect { it in NB_ALIASES ? 'nb' : (it in EMITS_ALIASES ? 'emu' : it) }
selected.each {
    if (!(it in VALID_CLASSIFIERS)) {
        exit 1, "ERROR: unknown classifier '${it}'; must be one of blca, ${is_its ? 'emits' : 'emu'}, nb " +
                "(or shorthand 'all' / 'both')."
    }
}
def run_blca = 'blca' in selected
def run_emu  = 'emu'  in selected
def run_nb   = 'nb'   in selected
def selected_display = selected.collect { (it == 'emu' && is_its) ? 'emits' : it }

// ---- reference database pre-flight -----------------------------------------
// `nextflow run` on its own builds no reference database (the launcher,
// bin/run_pipeline.py, does that on first use). So before any process is
// submitted, check that every reference file the selected classifiers need for
// this marker exists and is not empty, and stop with the build command for each
// missing one, rather than failing hours into the run. --check_db prints the
// same report and exits without running the pipeline. The tests mirror the
// launcher preflights and bin/doctor.sh: file present and non-empty, and for
// BLCA also the BLAST index that makeblastdb writes next to the FASTA. Under
// -preview the result is reported but not enforced (CI compiles the workflow
// without any database).
def _rel = { p ->
    def s = p.toString(); def root = "${projectDir}/"
    s.startsWith(root) ? s.substring(root.length()) : s
}
def _size = { long n ->
    n >= 1e9 ? String.format('%.1f GB', n / 1e9) : n >= 1e6 ? String.format('%.1f MB', n / 1e6)
             : n >= 1e3 ? String.format('%.0f kB', n / 1e3) : "${n} B"
}
// One reference file: OK when present and non-empty (a dangling symlink is MISSING).
def _file_item = { String cls, String what, String path, String fix ->
    def f = new File(path)
    def status = !f.isFile() ? 'MISSING' : f.length() == 0 ? 'EMPTY' : 'OK'
    [cls: cls, what: what, path: _rel(path), status: status, note: status == 'OK' ? _size(f.length()) : '', fix: fix]
}
// The BLAST index next to a BLCA FASTA: makeblastdb -out <fasta> writes <fasta>.n??
// (one volume) or <fasta>.NN.n?? volumes plus a <fasta>.nal alias for very large sets.
def _blast_item = { String cls, String fasta, String fix ->
    def f = new File(fasta).getAbsoluteFile()
    def dir = f.getParentFile()
    def idx = (dir?.isDirectory() ? dir.listFiles().toList() : []).findAll {
        it.name.startsWith(f.name + '.') && it.name ==~ /.*\.n[a-z]{2}$/ && it.isFile() && it.length() > 0
    }
    def ext = idx.collect { it.name.substring(f.name.length()) }.sort()
    def ok  = ('.nal' in ext) || ['.nhr', '.nin', '.nsq'].every { it in ext }
    [cls: cls, what: 'BLAST index', path: _rel(fasta) + '.n*',
     status: ok ? 'OK' : (idx ? 'INCOMPLETE' : 'MISSING'),
     note: ok ? "${idx.size()} files" : (idx ? "found only ${ext.join(' ')}" : ''), fix: fix]
}

def db_items = []
def db_howto
if (is_its) {
    def udir     = _rel(params.unite_db_dir)
    def fix_blca = "bash bin/build_unite_blca_db.sh <unite_release.fasta> ${udir} 0"
    def fix_nb   = "bash bin/build_unite_dada2_db.sh ${_rel(params.blca_db)} ${_rel(params.blca_tax)} ${_rel(params.unite_dada2_singlestep_db)}"
    def fix_emits = "bash bin/build_emits_db.sh <unite_release.fasta> ${udir}"
    if (run_blca) {
        db_items << _file_item('blca', 'UNITE BLCA FASTA',    params.blca_db.toString(),  fix_blca)
        db_items << _file_item('blca', 'UNITE BLCA taxonomy', params.blca_tax.toString(), fix_blca)
        db_items << _blast_item('blca', params.blca_db.toString(), fix_blca)
    }
    if (run_nb)  db_items << _file_item('nb',    'UNITE single-step NB reference', params.unite_dada2_singlestep_db.toString(), fix_nb)
    if (run_emu) db_items << _file_item('emits', 'EMITS minimap2 target, raw UNITE FASTA',  params.emits_db.toString(), fix_emits)
    db_howto = "Run the build commands from ${projectDir}. <unite_release.fasta> is the UNITE general FASTA release for Fungi\n" +
               "(sh_general_release_dynamic_*.fasta) downloaded from https://unite.ut.ee/repository.php. The launcher builds all three\n" +
               "references from that file in one go: python bin/run_pipeline.py --marker ITS --unite-fasta <unite_release.fasta> ...\n" +
               "See docs/databases.md."
} else {
    def rel  = (params.gtdb_release ?: '').toString() ==~ /\d+/ ? params.gtdb_release.toString() : '232'
    def dirs = [_rel(params.gtdb_db_dir), _rel(params.nb_db_dir), _rel(params.emu_db_dir)]
    def fix_gtdb = "bash bin/build_gtdb_db.sh ${rel} ${dirs[0]}" + ((dirs[1] != 'db_nb' || dirs[2] != 'db_emu') ? " ${dirs[1]} ${dirs[2]}" : '')
    if (run_blca) {
        db_items << _file_item('blca', 'GTDB BLCA FASTA',    params.blca_db.toString(),  fix_gtdb)
        db_items << _file_item('blca', 'GTDB BLCA taxonomy', params.blca_tax.toString(), fix_gtdb)
        db_items << _blast_item('blca', params.blca_db.toString(), fix_gtdb)
    }
    if (run_emu) {
        db_items << _file_item('emu', 'Emu database taxonomy',  "${params.emu_db_dir}/taxonomy.tsv",        fix_gtdb)
        db_items << _file_item('emu', 'Emu database sequences', "${params.emu_db_dir}/species_taxid.fasta", fix_gtdb)
    }
    if (run_nb) {
        db_items << _file_item('nb', 'DADA2 NB genus reference',   params.gtdb_dada2_genus_db.toString(),   fix_gtdb)
        db_items << _file_item('nb', 'DADA2 NB species reference', params.gtdb_dada2_species_db.toString(), fix_gtdb)
    }
    db_howto = "Run the build command from ${projectDir}; it downloads GTDB r${rel} once and asks which of BLCA, NB and Emu to format.\n" +
               "The launcher builds whatever is missing without those questions: python bin/run_pipeline.py --test --classifier ${selected_display.join(',')}\n" +
               "See docs/databases.md."
}
def db_missing = db_items.findAll { it.status != 'OK' }
def db_title   = "Reference database check: marker ${MARKER}, classifier(s) ${selected_display.join(',')}, paths relative to ${projectDir}"
def db_found   = db_items.collect { String.format('  %-6s %-10s %-44s %s%s', it.cls, it.status, it.path, it.what, it.note ? " (${it.note})" : '') }
def db_fixes   = db_missing.collect { String.format('  %-6s %-10s %s (%s)  ->  %s', it.cls, it.status, it.path, it.what, it.fix) }

if (params.check_db) {
    log.info(db_title + '\n' + db_found.join('\n'))
    if (db_missing) {
        log.info("Missing or incomplete, with the command that builds each:\n" + db_fixes.join('\n') + '\n' + db_howto)
        exit 1, "--check_db: ${db_missing.size()} of ${db_items.size()} reference file(s) missing or empty; nothing was run."
    }
    log.info("All ${db_items.size()} reference file(s) present. --check_db set, so exiting without running the pipeline.")
    exit 0
}
if (db_missing) {
    def db_msg = "ERROR: reference database check failed: ${db_missing.size()} of ${db_items.size()} reference file(s) missing or empty " +
                 "for marker ${MARKER}, classifier(s) ${selected_display.join(',')} (paths relative to ${projectDir}).\n" +
                 db_fixes.join('\n') + '\n' + db_howto + '\n' +
                 "Calling nextflow run directly does not build the databases; nothing has been run. " +
                 "Add --check_db to the same command to repeat this check on its own."
    if (workflow.preview) log.warn(db_msg + "\n(-preview: reported only, not enforced)")
    else                  error(db_msg)
}

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

// ---- marker-derived display names + which branches actually run ----------
def reference_name = is_its ? 'UNITE (db_unite)' : 'GTDB (db)'
def runs_dada2     = (run_blca || run_nb)   // the DADA2/ASV path runs ONLY for BLCA or NB
def trim_display   = is_its ? 'No (itsxrust trims the SSU/LSU flanks)' : trim_cutadapt
def read_prep_name = is_its
    ? ('itsxrust ITS extraction (no cutadapt)' + (runs_dada2 ? ' -> DADA2 for ASVs' : '; read-level EMITS only, no DADA2'))
    : ('cutadapt primer trim' + (runs_dada2 ? ' -> DADA2 for ASVs' : '') + (run_emu ? '; Emu on QC reads' : ''))

// Only print the sections for branches that will actually run, marker-aware.
def dada2_block = runs_dada2 ? """
  --- DADA2 (ASVs for BLCA / NB) ---
  DADA2 min/max len:         $params.min_len / $params.max_len
  DADA2 maxEE / minQ:        $params.max_ee / $params.minQ
  DADA2 pooling:             $params.pooling_method
  Min ASV total freq:        $dynamic_min_asv_totalfreq
  Min ASV samples:           $dynamic_min_asv_sample""" : ''

def blca_block = run_blca ? """
  --- BLCA branch ---
  BLCA db:                   $params.blca_db
  BLCA taxonomy:             $params.blca_tax
  BLCA chunk size (ASVs):    $params.blca_chunk_size""" : ''

def readem_block = run_emu ? (is_its ? """
  --- EMITS branch (read-level, on itsxrust reads) ---
  EMITS UNITE target:        $params.emits_db
  EMITS minimap2 preset:     map-hifi
  EMITS threads/sample:      $params.emits_threads
  EMITS read length filter:  $params.min_len-$params.max_len bp""" : """
  --- Emu branch (read-level, on QC reads) ---
  Emu DB dir:                $params.emu_db_dir
  Emu minimap2 preset:       $params.emu_type
  Emu threads/sample:        $params.emu_threads
  Emu read length filter:    $params.min_len-$params.max_len bp
  Emu min-abundance:         ${params.emu_min_abundance != null ? params.emu_min_abundance : 'default'}""") : ''

def nb_block = run_nb ? (is_its ? """
  --- NB branch (single-step 7-rank assignTaxonomy) ---
  UNITE 7-rank ref:          $params.unite_dada2_singlestep_db""" : """
  --- NB branch (two-step genus assignTaxonomy + addSpecies, bootstrap=$params.nb_min_bootstrap) ---
    GTDB genus ref   : $params.gtdb_dada2_genus_db
    GTDB species ref : $params.gtdb_dada2_species_db""") : ''

log_text = """
  HiFiTaxa pipeline
  =======================
  Samples in TSV:            $n_sample
  Marker:                    $MARKER
  Reference:                 $reference_name
  Read prep:                 $read_prep_name
  Filter reads above Q:      $params.filterQ
  Trim primers (cutadapt):   $trim_display
  Forward primer:            $params.forward_primer
  Reverse primer:            $params.reverse_primer
  Classifier(s):             ${selected_display.join(',')}${dada2_block}${blca_block}${readem_block}${nb_block}
"""

// ---- citations for the tools a run used -------------------------------------
// At the end of a successful run the pipeline prints which papers to cite for the
// tools and references that this marker and classifier set actually used, and
// writes the same list to <outdir>/CITATIONS.txt. Short forms only; the full
// author lists and DOIs are in CITATION.cff and README.md (References and what
// to cite). Keep the three files in step when a reference changes.
def REFS = [
    hifitaxa: "N. Garritano A, Modolon F, Warren CR (2026) HiFiTaxa: a reproducible analysis pipeline for full length 16S rRNA gene and ITS amplicon sequencing from PacBio reads. bioRxiv. doi:10.64898/2026.11.03.763332",
    nextflow: "Di Tommaso P et al. (2017) Nextflow enables reproducible computational workflows. Nat Biotechnol 35:316-319. doi:10.1038/nbt.3820",
    seqkit:   "Shen W et al. (2016) SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE 11:e0163962. doi:10.1371/journal.pone.0163962",
    cutadapt: "Martin M (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17:10-12. doi:10.14806/ej.17.1.200",
    qiime2:   "Bolyen E et al. (2019) Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37:852-857. doi:10.1038/s41587-019-0209-9",
    dada2:    "Callahan BJ et al. (2016) DADA2: high-resolution sample inference from Illumina amplicon data. Nat Methods 13:581-583. doi:10.1038/nmeth.3869",
    dada2_ccs: "Callahan BJ et al. (2019) High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide resolution. Nucleic Acids Res 47:e103. doi:10.1093/nar/gkz569",
    rdp_nb:   "Wang Q et al. (2007) Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol 73:5261-5267. doi:10.1128/AEM.00062-07",
    blca:     "Gao X et al. (2017) A Bayesian taxonomic classification method for 16S rRNA gene sequences with improved species-level accuracy. BMC Bioinformatics 18:247. doi:10.1186/s12859-017-1670-4",
    emu:      "Curry KD et al. (2022) Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data. Nat Methods 19:845-853. doi:10.1038/s41592-022-01520-4",
    itsxrust: "O'Brien A et al. (2026) ITSxRust: ITS region extraction with partial-chain recovery and structured diagnostics for long-read amplicon sequencing. Methods Ecol Evol 17:2843-2853. doi:10.1111/2041-210X.70393",
    emits:    "O'Brien A et al. (2026) EMITS: expectation-maximization abundance estimation for fungal ITS communities from long-read sequencing. Methods Ecol Evol 17:3233-3244. doi:10.1111/2041-210X.70418",
    gtdb:     "Parks DH et al. (2026) GTDB release 10: a complete and systematic taxonomy for 715 230 bacterial and 17 245 archaeal genomes. Nucleic Acids Res 54:D743-D754. doi:10.1093/nar/gkaf1040",
    unite:    ["Abarenkov K et al. (2024) The UNITE database for molecular identification and taxonomic communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered. Nucleic Acids Res 52:D791-D797. doi:10.1093/nar/gkad1039",
               "and the UNITE release you built from, e.g. Abarenkov K et al. (2025) UNITE general FASTA release for Fungi, version 10.0 (19 Feb 2025). UNITE Community. doi:10.15156/BIO/3301229"],
]
// its: marker ITS. reads: raw reads went through QC (false for -entry taxonomy_only).
// trim: cutadapt ran (16S without --skip_primer_trim). asvs: QIIME 2 + DADA2 denoising
// ran (BLCA or NB on ASVs). blca / nb / em: those classifiers ran; em is Emu for 16S
// and EMITS for ITS. NB is DADA2's assignTaxonomy, so it cites DADA2 and Wang et al.
// Denoising cites both DADA2 papers: the 2016 method and the 2019 PacBio CCS workflow
// that qiime dada2 denoise-ccs runs (for ITS as well as 16S).
def citation_text = { boolean its, boolean reads, boolean trim, boolean asvs, boolean blca, boolean nb, boolean em ->
    def rows = []
    // the pipeline itself, cited on its label line with no role text
    rows << ['HiFiTaxa', REFS.hifitaxa, []]
    rows << ['Nextflow', 'workflow engine', REFS.nextflow]
    if (reads)                 rows << ['seqkit',   'read QC, quality and length filtering', REFS.seqkit]
    if (reads && its)          rows << ['ITSxRust', 'ITS extraction from the QC reads', REFS.itsxrust]
    if (reads && !its && trim) rows << ['cutadapt', 'primer trimming', REFS.cutadapt]
    if (asvs)                  rows << ['QIIME 2',  'read import, DADA2 denoising wrapper and ASV tables', REFS.qiime2]
    if (asvs || nb)            rows << ['DADA2',    (asvs ? 'ASV denoising of the CCS reads (denoise-ccs)' : '') + (asvs && nb ? ' and ' : '') +
                                                    (nb ? (its ? 'the NB classifier (assignTaxonomy)' : 'the NB classifier (assignTaxonomy, addSpecies)') : ''),
                                                    asvs ? [REFS.dada2, REFS.dada2_ccs] : REFS.dada2]
    if (blca)                  rows << ['BLCA',     'ASV taxonomy with per-rank bootstrap confidence', REFS.blca]
    if (nb)                    rows << ['NB',       'the naive Bayes algorithm that DADA2 assignTaxonomy implements', REFS.rdp_nb]
    if (em && !its)            rows << ['Emu',      'read-level species profiling', REFS.emu]
    if (em && its)             rows << ['EMITS',    'read-level ITS profiling', REFS.emits]
    rows << [its ? 'UNITE' : 'GTDB', 'reference sequences and taxonomy', its ? REFS.unite : REFS.gtdb]
    // the third element is one reference or a list of them, each on its own indented line
    rows.collect { r -> String.format('  %-9s %s', r[0], r[1]) + [r[2]].flatten().collect { "\n            ${it}" }.join('') }.join('\n')
}
// Set by whichever workflow entry runs; read by workflow.onComplete at the end of the file.
def cite = [text: null]

workflow {
    if (!params.input)    { exit 1, "ERROR: --input <samples.tsv> is required" }
    if (!params.metadata) { exit 1, "ERROR: --metadata <metadata.tsv> is required" }

    log.info(log_text)
    write_log(log_text)
    cite.text = citation_text(is_its, true, !params.skip_primer_trim, run_blca || run_nb, run_blca, run_nb, run_emu)

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
            itsx_extract(QC_fastq.out.filtered_fastq, file(params.itsx_hmm))
            reads_for_import_skip = itsx_extract.out.fastq
            filtered_fastq_files  = itsx_extract.out.fastq.map { sid, fq -> fq }
            reads_for_emu         = itsx_extract.out.fastq   // EMITS on itsxrust-extracted reads (default)
        } else {
            reads_for_import_skip = QC_fastq.out.filtered_fastq
            filtered_fastq_files  = QC_fastq.out.filtered_fastq_files
        }
        if (run_blca || run_nb) {
            prepare_qiime2_manifest_skip_cutadapt(reads_for_import_skip.collect(), metadata_file)
            qiime2_manifest = prepare_qiime2_manifest_skip_cutadapt.out.sample_trimmed_file.flatten()
            import_qiime2(qiime2_manifest, filtered_fastq_files.collect())
        }
    } else {
        // Primer handling is marker-specific.
        //  - ITS runs NO cutadapt: itsxrust extracts the ITS1-5.8S-ITS2 span from the QC
        //    reads and trims the conserved SSU/LSU flanks that carry the primers (1391F in
        //    the SSU, ITS4ngsUni in the LSU), so both DADA2 (ASVs) and EMITS work on
        //    primer-free ITS reads without a cutadapt step.
        //  - 16S keeps cutadapt, but only on the DADA2/ASV path (there is no ITS
        //    extraction to remove primers before denoising); Emu runs on the QC reads
        //    directly, since minimap2 soft-clips the primer ends.
        def reads_for_import
        if (is_its) {
            itsx_extract(QC_fastq.out.filtered_fastq, file(params.itsx_hmm))
            reads_for_import     = itsx_extract.out.fastq
            filtered_fastq_files = itsx_extract.out.fastq.map { sid, fq -> fq }
            reads_for_emu        = itsx_extract.out.fastq            // EMITS on the itsxrust ITS-extracted reads
        } else {
            cutadapt(QC_fastq.out.filtered_fastq, dynamic_forward_primer, dynamic_reverse_primer)
            QC_fastq_post_trim(cutadapt.out.cutadapt_fastq)
            cutadapt_stats(cutadapt.out.summary_tocollect.collect())
            cutadapt_stats.out.stats.splitText().view { "[primer-removal] " + it.trim() }
            reads_for_import     = cutadapt.out.cutadapt_fastq
            filtered_fastq_files = cutadapt.out.cutadapt_fastq_files
            reads_for_emu        = QC_fastq.out.filtered_fastq       // Emu on the QC reads (no cutadapt)
        }
        if (run_blca || run_nb) {
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

    // ---------------- Final per-classifier ASV table (counts + taxonomy + confidence) ----------------
    // One wide TSV per ASV-based classifier: Feature ID, per-sample counts, the 7 ranks,
    // and confidence (per-rank for BLCA, single bootstrap for NB). Emu/EMITS are read-level
    // (no ASVs), so they are not merged. Published to $outdir/<classifier>_asv_table.tsv.
    merge_inputs = channel.empty()
    if (run_blca) {
        merge_inputs = merge_inputs.mix(
            filter_dada2.out.asv_freq_tsv
                .combine(merge_blca.out.tax_conf)
                .map { freq, tax -> tuple(freq, tax, 'blca', '--per-rank-confidence') })
    }
    if (run_nb) {
        nb_tax = is_its ? nb_classify_singlestep.out.best_nb_tax : nb_classify.out.best_nb_tax
        merge_inputs = merge_inputs.mix(
            filter_dada2.out.asv_freq_tsv
                .combine(nb_tax)
                .map { freq, tax -> tuple(freq, tax, 'nb', '--keep-confidence') })
    }
    merge_asv_table(merge_inputs)
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
    cite.text = citation_text(is_its, false, false, false, to_blca, to_nb, to_emu)

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

// ---- end of run: what to cite ------------------------------------------------
// Printed after a successful run (under -preview as the list the command would
// use) and written to <outdir>/CITATIONS.txt next to parameters.txt. Nothing is
// printed when the run failed, or when it stopped before the workflow started
// (--check_db, a failed reference database check).
workflow.onComplete {
    if (!workflow.success || !cite.text) return
    def header = workflow.preview
        ? "Tools this command would use; please cite them when you publish (-preview, nothing was run):"
        : "Please cite the tools this run used (full references in CITATION.cff and README.md):"
    def body = "\n  ${header}\n${cite.text}"
    if (workflow.preview) { log.info(body); return }
    try {
        def out = file(params.outdir).toAbsolutePath()
        out.mkdirs()
        out.resolve('CITATIONS.txt').text = "${header}\n${cite.text}\n"
        log.info(body + "\n  This list is also in ${out.resolve('CITATIONS.txt')}")
    } catch (Exception e) {
        log.info(body)
        log.warn("Could not write ${params.outdir}/CITATIONS.txt: ${e.message}")
    }
}
