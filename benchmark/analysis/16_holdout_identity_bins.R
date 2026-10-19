#!/usr/bin/env Rscript
# =====================================================================
# Holdout accuracy stratified by the identity of each test sequence to its
# best hit in the 90% training fraction (identity_<marker>.tsv written by
# identity.pbs or identity_resume.pbs). Same scoring as 01_holdout_accuracy.R (rescuable
# subset, BLCA@80 mask, ITS placeholder exclusion), split into identity bins.
# Usage: HIFITAXA_DATA=<data> HIFITAXA_OUT=<out> IDENT_DIR=<dir with identity_*.tsv>
#        Rscript 16_holdout_identity_bins.R
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
IDENT_DIR <- Sys.getenv("IDENT_DIR", unset = .SCRIPT_DIR)
BINS <- c("100", "99 to <100", "97 to <99", "<97", "no hit")
bin_of <- function(p) ifelse(is.na(p), "no hit", ifelse(p >= 100, "100", ifelse(p >= 99, "99 to <100", ifelse(p >= 97, "97 to <99", "<97"))))

score_bins <- function(marker) {
  hd <- file.path(BR, sprintf("benchmark_holdout_%s_seed42", marker)); rdir <- file.path(hd, "results")
  truth <- read_lineage_tsv(file.path(hd, "test_10.taxonomy"))
  ref <- read_lineage_tsv(file.path(hd, "reference_90.taxonomy"))
  train <- setNames(lapply(RANKS, function(r) unique(na.omit(ref[[r]]))), RANKS)
  bl <- read_blca(rdir, path_conf = NULL); nb <- read_nb(rdir)
  idt <- read_tsv(file.path(IDENT_DIR, sprintf("identity_%s.tsv", marker)), show_col_types = FALSE, progress = FALSE,
                  na = c("NA", "")) %>% mutate(pident = suppressWarnings(as.numeric(pident)))
  ids <- intersect(truth$id, bl$tab$id)
  T <- truth %>% filter(id %in% ids) %>% arrange(id); B <- bl$tab %>% filter(id %in% ids) %>% arrange(id)
  Bc <- bl$conf %>% filter(id %in% ids) %>% arrange(id); N <- nb %>% filter(id %in% ids) %>% arrange(id)
  stopifnot(all(T$id == B$id), all(T$id == Bc$id), all(T$id == N$id))
  bin <- bin_of(idt$pident[match(T$id, idt$id)])
  excl <- setNames(lapply(RANKS, function(r) rep(FALSE, nrow(T))), RANKS)
  if (marker == "its") { excl$Species <- !is.na(T$Species) & is_placeholder_sp(T$Species)
                         excl$Genus <- !is.na(T$Genus) & is_placeholder_genus(T$Genus) }
  Bmask <- B
  for (r in RANKS) { cc <- suppressWarnings(as.numeric(Bc[[paste0(r, "_conf")]])); cc[is.na(cc)] <- 0; Bmask[[r]][cc < 80] <- NA }
  bind_rows(lapply(RANKS, function(r) {
    tv <- T[[r]]; keep <- !is.na(tv) & (tv %in% train[[r]]) & !excl[[r]]
    bind_rows(lapply(BINS, function(bn) {
      k <- keep & bin == bn; denom <- sum(k)
      bind_rows(lapply(c("BLCA", "BLCA@80", "NB"), function(cl) {
        pv <- switch(cl, "BLCA" = B[[r]], "BLCA@80" = Bmask[[r]], "NB" = N[[r]])
        tibble(marker = marker, Rank = r, identity_bin = bn, classifier = cl, n = denom,
               accuracy = ifelse(denom > 0, 100 * sum(k & !is.na(pv) & pv == tv) / denom, NA_real_),
               pct_incorrect = ifelse(denom > 0, 100 * sum(k & !is.na(pv) & pv != tv) / denom, NA_real_),
               pct_unassigned = ifelse(denom > 0, 100 * sum(k & is.na(pv)) / denom, NA_real_))
      }))
    }))
  })) %>% mutate(identity_bin = factor(identity_bin, levels = BINS))
}
res <- bind_rows(score_bins("16s"), score_bins("its"))
write_csv(res, file.path(OUT, "tables", "holdout_accuracy_by_identity.csv"))
cat("\n== Share of test sequences per identity bin (rescuable at species) ==\n")
res %>% filter(Rank == "Species", classifier == "BLCA") %>% group_by(marker) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>% select(marker, identity_bin, n, pct) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\n== Accuracy (%) by identity bin, genus and species ==\n")
res %>% filter(Rank %in% c("Genus", "Species")) %>% mutate(accuracy = round(accuracy, 1)) %>%
  select(marker, Rank, identity_bin, classifier, n, accuracy) %>%
  pivot_wider(names_from = classifier, values_from = accuracy) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nWrote tables/holdout_accuracy_by_identity.csv\n")
