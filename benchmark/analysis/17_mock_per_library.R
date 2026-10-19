#!/usr/bin/env Rscript
# =====================================================================
# ATCC MSA-1003 mock: genus level composition PER LIBRARY for every classifier
# (Supplementary Table S13). Same folding and the same taxon rows as Table S3
# (04_mock_composition.R): the 18 mock genera, Unassigned (ASV classifiers only),
# the foreign genera that any classifier reports above FP_THRESH of the pooled
# reads, and one row for all other foreign genera. For BLCA, BLCA@80 and NB the
# value is the share of the library's ASV reads; for Emu it is Emu's per library
# relative abundance (its tables are renormalised over assigned reads, so Emu
# has no Unassigned row). mean_of_libraries is the plain mean over the eight
# libraries; Table S3 pools the reads of the ASV classifiers instead, so its
# values differ slightly from that mean for BLCA, BLCA@80 and NB.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
dir <- file.path(BR, "results_atcc_16s")
# per library ASV counts
ft <- read_tsv(file.path(dir, "dada2", "feature_table.tsv"), comment = "# Constructed", show_col_types = FALSE, progress = FALSE)
names(ft)[1] <- "id"; libs <- sort(setdiff(names(ft), "id"))
long_counts <- ft %>% pivot_longer(all_of(libs), names_to = "library", values_to = "count")
asv_comp <- function(tab, name) {
  g <- tibble(id = tab$id, taxon = norm_genus(tab$Genus))
  long_counts %>% inner_join(g, by = "id") %>% mutate(taxon = ifelse(is.na(taxon), "Unassigned", taxon)) %>%
    group_by(library, taxon) %>% summarise(count = sum(count), .groups = "drop") %>%
    group_by(library) %>% mutate(relab = 100 * count / sum(count)) %>% ungroup() %>%
    transmute(classifier = name, library, taxon, relab)
}
bl <- read_blca(dir)$tab; b8 <- read_blca_at80(dir); nb <- read_nb(dir)
# Emu per library (columns of emu_genus_table.tsv; NA = not reported = 0)
e <- read_tsv(file.path(dir, "taxonomy_emu", "emu_genus_table.tsv"), show_col_types = FALSE, progress = FALSE)
esamp <- setdiff(names(e), c("species", "genus", "family", "order", "class", "phylum", "superkingdom"))
emu <- e %>% transmute(taxon = norm_genus(clean(str_replace(genus, "^genus:", ""))), across(all_of(esamp))) %>%
  filter(!is.na(taxon)) %>% pivot_longer(all_of(esamp), names_to = "library", values_to = "relab") %>%
  mutate(relab = ifelse(is.na(relab), 0, relab)) %>% group_by(library, taxon) %>%
  summarise(relab = sum(relab), .groups = "drop") %>% group_by(library) %>%
  mutate(relab = 100 * relab / sum(relab)) %>% ungroup() %>% transmute(classifier = "Emu", library, taxon, relab)
comp <- bind_rows(asv_comp(bl, "BLCA"), asv_comp(b8, "BLCA@80"), emu, asv_comp(nb, "NB"))
# row set as in Table S3: mock genera, Unassigned, foreign genera above FP_THRESH of the pooled composition
pooled <- read_csv(file.path(OUT, "tables", "atcc_composition_genus.csv"), show_col_types = FALSE, progress = FALSE)
listed_foreign <- setdiff(pooled$taxon, c(atcc_genus, "Unassigned"))
comp <- comp %>% mutate(row = case_when(taxon %in% atcc_genus ~ taxon, taxon == "Unassigned" ~ "Unassigned",
                                        taxon %in% listed_foreign ~ taxon, TRUE ~ "Other foreign genera")) %>%
  group_by(classifier, library, row) %>% summarise(relab = sum(relab), .groups = "drop")
expected <- tibble(row = names(ATCC_TRUTH), expected = 100 * as.numeric(ATCC_TRUTH)) %>%
  mutate(row = norm_genus(str_split_fixed(row, " ", 2)[, 1])) %>% group_by(row) %>% summarise(expected = sum(expected), .groups = "drop")
row_order <- c(expected %>% arrange(desc(expected), row) %>% pull(row), "Unassigned", sort(listed_foreign), "Other foreign genera")
wide <- comp %>% pivot_wider(names_from = library, values_from = relab, values_fill = 0) %>%
  right_join(tidyr::expand_grid(classifier = c("BLCA", "BLCA@80", "Emu", "NB"), row = row_order), by = c("classifier", "row")) %>%
  mutate(across(all_of(libs), ~ replace_na(.x, 0))) %>%
  left_join(expected, by = "row") %>% mutate(expected = replace_na(expected, 0)) %>%
  mutate(mean_of_libraries = rowMeans(across(all_of(libs)))) %>%
  mutate(classifier = factor(classifier, levels = c("BLCA", "BLCA@80", "Emu", "NB")), row = factor(row, levels = row_order)) %>%
  arrange(classifier, row) %>% select(classifier, taxon = row, expected, all_of(libs), mean_of_libraries) %>%
  filter(!(classifier == "Emu" & taxon == "Unassigned"))
write_csv(wide, file.path(OUT, "tables", "atcc_composition_genus_per_library.csv"))
cat("Wrote tables/atcc_composition_genus_per_library.csv:", nrow(wide), "rows\n")
print(as.data.frame(wide %>% filter(taxon %in% c("Escherichia", "Porphyromonas", "Mycobacterium", "Unassigned")) %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
