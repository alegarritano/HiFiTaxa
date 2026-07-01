# HiFiTaxa — DADA2 taxonomy, single-step (7-rank assignTaxonomy, NO addSpecies).
#
# Background
# ----------
# This is the ITS-appropriate Naive-Bayes design. The two-step GTDB classifier
# (scripts/dada2_assign_tax.R) does a genus-level assignTaxonomy() bootstrap and
# overlays species by exact-match addSpecies(). That works for 16S/GTDB but
# collapses to ~0 species on ITS/UNITE (the addSpecies exact-match overlay finds
# almost no full-length hits on the extracted ITS region), so for ITS we run a
# SINGLE assignTaxonomy() call straight to species against a 7-rank
# (Kingdom..Species) reference, with the species lineage already in the headers.
#
# Unlike full GTDB SSU r232 (957k seqs / ~82k species lineages) — where so many
# near-identical terminal labels dilute the Wang/RDP bootstrap to Kingdom-only —
# UNITE is small enough that a single-step assignTaxonomy() commits to
# species-level picks. minBoot defaults to 80, matching the two-step path.
#
# Mirrors the R logic in
#   HiFiTaxa_Fungi/bin/classify_nb_singlestep.sh
# (assignTaxonomy, outputBootstraps=TRUE, taxLevels=Kingdom..Species), wired to
# the HiFiTaxa 4-arg contract and fixed output filenames.
#
# Usage: Rscript dada2_assign_tax_singlestep.R <query.fasta> <cpus> <ref_db> <minBoot>
#   <query.fasta> : ASV / query sequences (FASTA)
#   <cpus>        : multithread count
#   <ref_db>      : 7-rank DADA2 reference (Kingdom..Species in headers), e.g.
#                   db_unite/unite_full_singlestep_ref.fa.gz
#   <minBoot>     : assignTaxonomy bootstrap cutoff (default 80 if omitted)
#
# Writes:
#   best_taxonomy.tsv : QIIME2 TSVTaxonomyFormat (Feature ID, Taxon, Confidence),
#                       7-rank lineage "k__..; p__..; ...; s__.."
#   nb_singlestep.tsv : per-ASV 7-rank table (Feature ID + Kingdom..Species)

library(dada2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 4) {
  stop(
    "Usage: Rscript dada2_assign_tax_singlestep.R <query.fasta> <cpus> <ref_db> <minBoot>",
    call. = FALSE
  )
}

seqs_path   <- args[1]
threads     <- as.numeric(args[2])
ref_db      <- args[3]
minBoot_num <- if (length(args) == 4) as.numeric(args[4]) else 80

RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

seqs   <- getSequences(seqs_path)
otu_id <- names(seqs)
if (is.null(otu_id)) otu_id <- paste0("seq", seq_along(seqs))

# ----- Single-step 7-rank Naive-Bayes (assignTaxonomy straight to species) -----
nb <- assignTaxonomy(seqs,
  refFasta = ref_db, minBoot = minBoot_num,
  multithread = threads, outputBootstraps = TRUE,
  taxLevels = RANKS
)
tax  <- nb$tax    # character matrix, columns Kingdom..Species (NA below minBoot)
boot <- nb$boot   # integer matrix, bootstrap support 0-100 per rank

# ----- Normalise into per-rank data frames ------------------------------------
tax_df  <- as.data.frame(tax,  stringsAsFactors = FALSE)
boot_df <- as.data.frame(boot, stringsAsFactors = FALSE)

# Ensure all rank columns exist (assignTaxonomy may return fewer if td<7).
for (r in RANKS) if (!(r %in% colnames(tax_df)))  tax_df[[r]]  <- NA_character_
for (r in RANKS) if (!(r %in% colnames(boot_df))) boot_df[[r]] <- NA_integer_
tax_df  <- tax_df[,  RANKS, drop = FALSE]
boot_df <- boot_df[, RANKS, drop = FALSE]

# ----- nb_singlestep.tsv : per-ASV 7-rank table -------------------------------
nb_out <- data.frame(
  "Feature ID" = otu_id,
  tax_df,
  check.names = FALSE, stringsAsFactors = FALSE
)
write.table(nb_out, "nb_singlestep.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE)

# ----- best_taxonomy.tsv : QIIME2 TSVTaxonomyFormat (Feature ID, Taxon, Confidence)
# Taxon = "k__..; p__..; ...; s__..", Unclassified where NA.
# Confidence = NB bootstrap at the FIRST unclassified rank, or the Species
# bootstrap if classified all the way to species.
prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")

to_save <- data.frame()
for (i in seq_along(seqs)) {
  rank_vals <- as.character(unlist(tax_df[i, RANKS]))   # length 7
  # Confidence from NB bootstrap
  na_rk <- which(is.na(rank_vals))
  if (length(na_rk) > 0) {
    conf <- boot_df[i, min(na_rk)]
  } else {
    conf <- boot_df[i, length(RANKS)]   # fully classified to Species -> Species boot
  }
  lineage_vals <- rank_vals
  lineage_vals[is.na(lineage_vals)] <- "Unclassified"
  taxon <- paste(paste0(prefixes, lineage_vals), collapse = "; ")

  to_save <- rbind(to_save, data.frame(
    "Feature ID" = otu_id[i], "Taxon" = taxon, "Confidence" = conf,
    check.names = FALSE, stringsAsFactors = FALSE))
}

write.table(to_save, "best_taxonomy.tsv",
  quote = FALSE, sep = "\t", row.names = FALSE,
  col.names = c("Feature ID", "Taxon", "Confidence"))
