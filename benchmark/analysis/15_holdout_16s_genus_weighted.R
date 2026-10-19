#!/usr/bin/env Rscript
# =====================================================================
# REVISED 2026-10: decision 5 (new script). 16S holdout, genus rank, two views.
#
#   Per sequence      every test sequence counts once (the value in holdout_accuracy.csv).
#   Per genus (equal) accuracy is computed inside each genus and the genera are then
#                     averaged with equal weight, so a genus with 20,000 test sequences
#                     counts as much as a genus with one.
#
# Scored subset (same as 01_holdout_accuracy.R at genus rank): test sequences whose
# true genus is present in the 90% training fraction. A genus enters the equal weight
# mean when it has at least one such test sequence.
#
# The script also measures how much of the per sequence difference between BLCA and NB
# comes from the Escherichia group: Escherichia and the two genera that GTDB r232 split
# from it (ECMA0423, G047199095).
#
# Outputs:
#   tables/holdout_16s_genus_weighted.csv          the summary table
#   tables/holdout_16s_genus_weighted_ecoli.csv    the three Escherichia group genera
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
ECOLI_GROUP <- c("Escherichia","ECMA0423","G047199095")

hd   <- file.path(BR, "benchmark_holdout_16s_seed42"); rdir <- file.path(hd, "results")
truth <- read_lineage_tsv(file.path(hd,"test_10.taxonomy"))
train <- read_lineage_tsv(file.path(hd,"reference_90.taxonomy"))
train_genus <- unique(na.omit(train$Genus))
bl <- read_blca(rdir); nb <- read_nb(rdir)

# same sequence set and order as 01_holdout_accuracy.R
ids <- intersect(truth$id, bl$tab$id)
T  <- truth   %>% filter(id %in% ids) %>% arrange(id)
B  <- bl$tab  %>% filter(id %in% ids) %>% arrange(id)
Bc <- bl$conf %>% filter(id %in% ids) %>% arrange(id)
N  <- nb      %>% filter(id %in% ids) %>% arrange(id)
stopifnot(all(T$id==B$id), all(T$id==Bc$id), all(T$id==N$id))
cc <- suppressWarnings(as.numeric(Bc$Genus_conf)); cc[is.na(cc)] <- 0
B80 <- B$Genus; B80[cc < 80] <- NA                      # BLCA@80 at genus rank

keep <- !is.na(T$Genus) & (T$Genus %in% train_genus)   # genus present in training
d <- tibble(genus = T$Genus[keep],
            BLCA      = !is.na(B$Genus[keep]) & B$Genus[keep] == T$Genus[keep],
            `BLCA@80` = !is.na(B80[keep])     & B80[keep]     == T$Genus[keep],
            NB        = !is.na(N$Genus[keep]) & N$Genus[keep] == T$Genus[keep])
CL <- c("BLCA","BLCA@80","NB")

per_genus <- d %>% group_by(genus) %>%
  summarise(n = n(), across(all_of(CL), sum, .names = "correct_{.col}"), .groups = "drop")
for (cl in CL) per_genus[[paste0("acc_", cl)]] <- 100 * per_genus[[paste0("correct_", cl)]] / per_genus$n

# one row of the summary table for a set of genera and a weighting
total_diff <- sum(per_genus$correct_BLCA) - sum(per_genus$correct_NB)   # per sequence BLCA minus NB
summarise_set <- function(pg, subset_label, weighting) {
  if (weighting == "Per sequence") {
    acc <- sapply(CL, function(cl) 100 * sum(pg[[paste0("correct_", cl)]]) / sum(pg$n))
  } else {
    acc <- sapply(CL, function(cl) mean(pg[[paste0("acc_", cl)]]))
  }
  diff_n <- sum(pg$correct_BLCA) - sum(pg$correct_NB)
  tibble(genera = subset_label, weighting = weighting,
         n_genera = nrow(pg), n_sequences = sum(pg$n),
         BLCA = acc[["BLCA"]], `BLCA@80` = acc[["BLCA@80"]], NB = acc[["NB"]],
         BLCA_minus_NB_points = acc[["BLCA"]] - acc[["NB"]],
         correct_BLCA = sum(pg$correct_BLCA), `correct_BLCA@80` = sum(pg$`correct_BLCA@80`),
         correct_NB = sum(pg$correct_NB),
         BLCA_minus_NB_sequences = diff_n,
         share_of_per_sequence_difference_pct = 100 * diff_n / total_diff)
}
eco   <- per_genus %>% filter(genus %in% ECOLI_GROUP)
other <- per_genus %>% filter(!genus %in% ECOLI_GROUP)
tab <- bind_rows(
  summarise_set(per_genus, "All genera", "Per sequence"),
  summarise_set(per_genus, "All genera", "Per genus (equal weight)"),
  summarise_set(eco,   "Escherichia group (Escherichia, ECMA0423, G047199095)", "Per sequence"),
  summarise_set(other, "All genera except the Escherichia group", "Per sequence"),
  summarise_set(other, "All genera except the Escherichia group", "Per genus (equal weight)"))

dir.create(file.path(OUT,"tables"), showWarnings = FALSE, recursive = TRUE)
write_csv(tab, file.path(OUT,"tables","holdout_16s_genus_weighted.csv"))
write_csv(eco %>% arrange(desc(n)) %>%
            mutate(BLCA_minus_NB_sequences = correct_BLCA - correct_NB,
                   share_of_per_sequence_difference_pct = 100 * BLCA_minus_NB_sequences / total_diff),
          file.path(OUT,"tables","holdout_16s_genus_weighted_ecoli.csv"))

# consistency check against the per sequence value written by 01_holdout_accuracy.R
f01 <- file.path(OUT,"tables","holdout_accuracy.csv")
if (file.exists(f01)) {
  a <- read_csv(f01, show_col_types = FALSE, progress = FALSE) %>% filter(marker == "16s", Rank == "Genus")
  for (cl in CL) stopifnot(a$n[a$classifier == cl] == sum(per_genus$n),
                           a$correct[a$classifier == cl] == sum(per_genus[[paste0("correct_", cl)]]))
  cat("Per sequence counts agree with holdout_accuracy.csv (16s, Genus).\n")
}

options(width = 220)
cat("\n== 16S holdout, genus rank: per sequence and per genus (equal weight) accuracy (%) ==\n")
tab %>% mutate(across(c(BLCA, `BLCA@80`, NB, BLCA_minus_NB_points, share_of_per_sequence_difference_pct), ~round(., 2))) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\n== Escherichia group genera ==\n")
eco %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nWrote tables/holdout_16s_genus_weighted.csv and tables/holdout_16s_genus_weighted_ecoli.csv\n")
