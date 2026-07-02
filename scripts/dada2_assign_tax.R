# HiFiTaxa — DADA2 taxonomy, GTDB-only, two-step (genus NB + exact-match species).
#
# Background
# ----------
# PacBio's HiFi-16S-workflow runs DADA2 assignTaxonomy() with the SPECIES-level
# lineage in the reference headers. That works against their small (~62k-seq)
# pre-built GTDB file, but it does NOT scale to the full GTDB SSU r232 reference
# (957,312 seqs / ~81,774 unique species lineages): with that many near-identical
# terminal labels the Wang/RDP bootstrap cannot commit to a pick and collapses to
# Kingdom-only assignments (verified empirically; ~Desktop/nb_scale_test/RESULTS.md).
#
# This is the documented reason the standard DADA2 SILVA/RDP training sets are
# built at GENUS level for assignTaxonomy(), with species assigned separately by
# exact matching via addSpecies(). HiFiTaxa therefore uses that canonical two-step
# design so NB runs against the SAME full GTDB release as BLCA and Emu:
#
#   1. assignTaxonomy() against a GENUS-level (6-rank) reference  -> Kingdom..Genus
#      with per-rank bootstrap support (minBoot, outputBootstraps).
#   2. addSpecies() against a species-assignment reference        -> Species by
#      exact (100% identity, tryRC) matching, only where the genus is consistent.
#
# Note: species via exact matching is high-precision but low-recall on full-length
# reads (only ASVs that exactly match a GTDB genome's 16S get a species call). The
# reported Confidence is the NB bootstrap at the deepest assigned rank (genus
# level); the species overlay is exact-match (no bootstrap).
#
# Args: <asv_fasta> <threads> <genus_db> <species_db> [minBoot=80]

library(dada2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4 || length(args) > 6) {
  stop(
    "Usage: Rscript dada2_assign_tax.R <asv_fasta> <threads> <genus_db> <species_db> [minBoot=80] [seed=42]",
    call. = FALSE
  )
}

seqs_path   <- args[1]
threads     <- as.numeric(args[2])
genus_db    <- args[3]
species_db  <- args[4]
minBoot_num <- if (length(args) >= 5) as.numeric(args[5]) else 80
seed_num    <- if (length(args) >= 6) as.integer(args[6]) else 42L

# Seed the RNG so the assignTaxonomy Wang/RDP bootstrap is reproducible across
# runs. Best-effort under multithread (DADA2 forks workers); set.seed fixes the
# sampling stream assignTaxonomy draws from.
set.seed(seed_num)

RANKS    <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
NB_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")  # from assignTaxonomy

seqs   <- getSequences(seqs_path)
otu_id <- names(seqs)
if (is.null(otu_id)) otu_id <- paste0("seq", seq_along(seqs))

# ----- Step 1: genus-level Naive-Bayes (assignTaxonomy, genus reference) ------
nb <- assignTaxonomy(seqs,
  refFasta = genus_db, minBoot = minBoot_num,
  multithread = threads, outputBootstraps = TRUE
)
tax  <- nb$tax    # character matrix, columns Kingdom..Genus (NA below minBoot)
boot <- nb$boot   # integer matrix, bootstrap support 0-100 per NB rank

# ----- Step 2: species by exact match (REFERENCE-CHUNKED addSpecies) -----------
# addSpecies appends a "Species" column (epithet only) where an exact-match
# reference of the SAME genus as the assignTaxonomy call exists; tryRC handles
# orientation. Two scale guards are required at full-GTDB reference size:
#
#  (a) Non-ACGT queries: addSpecies aborts the ENTIRE batch if ANY single query
#      contains a non-ACGT base ("Non-ACGT characters present in the query
#      sequences"). One ambiguous ASV would silently drop the species call for
#      ALL ASVs, so we run addSpecies only on the ACGT-clean ASVs; any ASV with a
#      non-ACGT base gets Species = NA individually instead of poisoning the batch.
#
#  (b) Reference scale: a SINGLE addSpecies() call against a full GTDB-size
#      species reference (~10^6 seqs) SILENTLY under-detects exact matches -- a
#      match that an ~80k-seq reference chunk finds is missed when the same seqs
#      sit in the whole reference (verified on the r232 holdout: single-shot
#      committed 0/300 reads that a chunked reference committed 298/300), and it
#      can also hit R's vector limit. So we split the reference into chunks, run
#      addSpecies(allowMultiple=TRUE) against each, union the per-query epithet
#      sets, and keep the species only where exactly ONE survives across all
#      chunks. That reproduces whole-reference addSpecies(allowMultiple=FALSE)
#      without the scale failure.
SPECIES_REF_CHUNK_SIZE <- 80000L

split_reference_chunks <- function(ref_path, chunk_size, out_dir) {
  ref <- Biostrings::readDNAStringSet(ref_path)   # headers preserved verbatim
  n <- length(ref)
  if (n == 0L) return(character(0))
  grp   <- (seq_len(n) - 1L) %/% chunk_size
  files <- character(0)
  for (g in sort(unique(grp))) {
    fn <- file.path(out_dir, sprintf("refchunk_%04d.fa", g))
    Biostrings::writeXStringSet(ref[grp == g], filepath = fn)
    files <- c(files, fn)
  }
  files
}

species_epithet <- rep(NA_character_, length(seqs))
seq_chr  <- as.character(seqs)
is_clean <- !grepl("[^ACGTacgt]", seq_chr)
if (sum(!is_clean) > 0) {
  message(sprintf(
    "addSpecies: %d/%d ASVs contain non-ACGT bases; running species on the %d ACGT-clean ASVs only (Species=NA for the rest).",
    sum(!is_clean), length(seqs), sum(is_clean)))
}
if (any(is_clean)) {
  # addSpecies matches by sequence; dedup the clean ASVs and map results back.
  clean_seq <- seq_chr[is_clean]
  uniq      <- !duplicated(clean_seq)
  tax_in    <- tax[is_clean, , drop = FALSE][uniq, , drop = FALSE]
  rownames(tax_in) <- clean_seq[uniq]

  chunk_dir <- file.path(tempdir(), "hifitax_refchunks")
  unlink(chunk_dir, recursive = TRUE)
  dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
  chunk_files <- tryCatch(
    split_reference_chunks(species_db, SPECIES_REF_CHUNK_SIZE, chunk_dir),
    error = function(e) { message("reference split failed (", conditionMessage(e), ")"); character(0) })
  message(sprintf("addSpecies: species reference split into %d chunk(s) of <=%d seqs; %d unique clean ASVs.",
                  length(chunk_files), SPECIES_REF_CHUNK_SIZE, nrow(tax_in)))

  spec_sets <- vector("list", nrow(tax_in))
  for (cf in chunk_files) {
    sp <- tryCatch(
      addSpecies(tax_in, refFasta = cf, tryRC = TRUE,
                 allowMultiple = TRUE, verbose = FALSE),
      error = function(e) {
        message("addSpecies failed on ", basename(cf), " (", conditionMessage(e), "); skipping chunk.")
        NULL
      })
    if (is.null(sp) || !("Species" %in% colnames(sp))) next
    spv <- as.character(sp[, "Species"])
    for (i in which(!is.na(spv) & nzchar(spv)))
      spec_sets[[i]] <- union(spec_sets[[i]], strsplit(spv[i], "/", fixed = TRUE)[[1]])
  }
  unlink(chunk_dir, recursive = TRUE)

  # exactly one species across all chunks -> assign; >1 (ambiguous) or 0 -> NA.
  sp_uniq <- vapply(spec_sets,
    function(s) if (length(s) == 1L) s else NA_character_, character(1))
  names(sp_uniq) <- rownames(tax_in)
  species_epithet[is_clean] <- unname(sp_uniq[clean_seq])
}

# ----- Normalise into per-rank data frames ------------------------------------
tax_df  <- as.data.frame(tax,  stringsAsFactors = FALSE)
boot_df <- as.data.frame(boot, stringsAsFactors = FALSE)

# Ensure all NB rank columns exist (assignTaxonomy may return fewer if td<6)
for (r in NB_RANKS) if (!(r %in% colnames(tax_df)))  tax_df[[r]]  <- NA_character_
for (r in NB_RANKS) if (!(r %in% colnames(boot_df))) boot_df[[r]] <- NA_integer_
tax_df  <- tax_df[,  NB_RANKS, drop = FALSE]
boot_df <- boot_df[, NB_RANKS, drop = FALSE]

# Build the binomial species string: "Genus epithet" where both present.
genus_vec  <- tax_df[["Genus"]]
binomial   <- ifelse(!is.na(species_epithet) & !is.na(genus_vec),
                     paste(genus_vec, species_epithet), NA_character_)

# ----- gtdb_nb.tsv : per-rank detail (Kingdom..Genus NB + exact-match Species) -
# Columns: Feature.ID, Kingdom..Genus, Species, boot.Kingdom..boot.Genus, Assignment.
# Species has no bootstrap column (it is exact-match, not NB).
gtdb_out <- data.frame(
  "Feature.ID" = otu_id,
  tax_df,
  Species = binomial,
  setNames(boot_df, paste0("boot.", NB_RANKS)),
  Assignment = "GTDB r232 (genus NB + exact-match species)",
  check.names = FALSE, stringsAsFactors = FALSE
)
write.table(gtdb_out, "gtdb_nb.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE)

# ----- best_taxonomy.tsv : QIIME2 TSVTaxonomyFormat (Feature ID, Taxon, Confidence)
# Taxon = "d__..; p__..; ...; s__..", Unclassified where NA.
# Confidence = NB bootstrap at the FIRST unclassified NB rank, or the Genus
# bootstrap if classified all the way to genus (mirrors PacBio's confidence
# semantics, evaluated over the 6 NB ranks; species is an exact-match overlay).
prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")

to_save  <- data.frame()
to_save2 <- data.frame()
for (i in seq_along(seqs)) {
  nb_vals <- as.character(unlist(tax_df[i, NB_RANKS]))   # length 6
  # Confidence from NB bootstrap
  na_nb <- which(is.na(nb_vals))
  if (length(na_nb) > 0) {
    conf <- boot_df[i, min(na_nb)]
  } else {
    conf <- boot_df[i, length(NB_RANKS)]   # fully classified to Genus -> Genus boot
  }
  # Full 7-level lineage (NB ranks + exact-match species)
  lineage_vals <- c(nb_vals, binomial[i])
  lineage_vals[is.na(lineage_vals)] <- "Unclassified"
  taxon <- paste(paste0(prefixes, lineage_vals), collapse = "; ")

  to_save  <- rbind(to_save,  data.frame(
    "Feature ID" = otu_id[i], "Taxon" = taxon, "Confidence" = conf,
    check.names = FALSE, stringsAsFactors = FALSE))
  to_save2 <- rbind(to_save2, data.frame(
    "Feature ID" = otu_id[i], "Taxon" = taxon, "Confidence" = conf,
    "Assignment Database" = "GTDB r232",
    check.names = FALSE, stringsAsFactors = FALSE))
}

write.table(to_save, "best_taxonomy.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE,
  col.names = c("Feature ID", "Taxon", "Confidence"))

write.table(to_save2, "best_taxonomy_withDB.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE,
  col.names = c("Feature ID", "Taxon", "Confidence", "Assignment Database"))
