#!/usr/bin/env Rscript
# =====================================================================
# Read-weighted mock accuracy (supports the manuscript's fungal-mock prose:
# "% of reads in a correct mock genus", "% reads to a correct mock species",
# and the foreign-read fraction). Complements the taxon-count detection (02).
#
# For each classifier x rank:
#   reads_correct  = reads whose folded taxon is an EXPECTED mock taxon
#   reads_foreign  = reads whose folded taxon is outside every mock genus
#                    (excludes unclassified/placeholder = "unassigned")
#   accuracy       = reads_correct / total classified reads
#   foreign_frac   = reads_foreign / total classified reads
# BLCA/NB weight ASVs by feature_table counts; Emu/EMITS use their read counts.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
# reads per folded taxon for a classifier (returns tibble taxon, reads)
asv_reads <- function(ids, taxa, cnt, rank) {
  tb <- tibble(id=ids, taxon=if(rank=="Species") norm_species(taxa) else norm_genus(taxa)) %>%
        inner_join(cnt, by="id")
  tb %>% group_by(taxon) %>% summarise(reads=sum(count),.groups="drop")
}
emu_reads <- function(dir, rank, marker) {
  if (marker=="16s") {
    e <- read_emu(dir, tolower(rank))   # relab; scale to pseudo-reads by *1e6 (ratios only)
    e %>% transmute(taxon=if(rank=="Species") norm_species(taxon) else norm_genus(taxon), reads=relab*1e6) %>%
      group_by(taxon) %>% summarise(reads=sum(reads),.groups="drop")
  } else {
    e <- read_emits(dir, tolower(rank))
    e %>% transmute(taxon=if(rank=="Species") norm_species(taxon) else norm_genus(taxon), reads=relab*1e6) %>%
      group_by(taxon) %>% summarise(reads=sum(reads),.groups="drop")
  }
}

score_rw <- function(tb, exp_set, exp_genus, rank) {
  tb <- tb %>% filter(!is.na(taxon))
  total <- sum(tb$reads)
  correct <- sum(tb$reads[tb$taxon %in% exp_set])
  # foreign = folded GENUS not in the mock genus set (species rows -> take genus token)
  gtok <- if (rank=="Species") norm_genus(str_split_fixed(tb$taxon," ",2)[,1]) else tb$taxon
  foreign <- sum(tb$reads[!(gtok %in% exp_genus)])
  tibble(rank=rank, reads_total=round(total), accuracy=100*correct/total, foreign_frac=100*foreign/total)
}

run_rw <- function(marker) {
  if (marker=="16s") { dir<-file.path(BR,"results_atcc_16s"); es<-atcc_species; eg<-atcc_genus; rl<-"Emu" }
  else                { dir<-file.path(BR,"results_fungi_its"); es<-fung_species; eg<-fung_genus; rl<-"EMITS" }
  cnt<-read_feature_table(dir); bl<-read_blca(dir)$tab; nb<-read_nb(dir)
  bind_rows(
    score_rw(asv_reads(bl$id,bl$Genus,cnt,"Genus"),   eg,eg,"Genus")  %>% mutate(classifier="BLCA"),
    score_rw(asv_reads(bl$id,bl$Species,cnt,"Species"),es,eg,"Species")%>% mutate(classifier="BLCA"),
    score_rw(asv_reads(nb$id,nb$Genus,cnt,"Genus"),   eg,eg,"Genus")  %>% mutate(classifier="NB"),
    score_rw(asv_reads(nb$id,nb$Species,cnt,"Species"),es,eg,"Species")%>% mutate(classifier="NB"),
    score_rw(emu_reads(dir,"Genus",marker),   eg,eg,"Genus")  %>% mutate(classifier=rl),
    score_rw(emu_reads(dir,"Species",marker), es,eg,"Species")%>% mutate(classifier=rl)
  ) %>% mutate(marker=ifelse(marker=="16s","16S (ATCC MSA-1003)","Fungal ITS (103-sp mock)"))
}

rw <- bind_rows(run_rw("16s"), run_rw("its"))
dir.create(file.path(OUT,"tables"),showWarnings=FALSE,recursive=TRUE)
write_csv(rw, file.path(OUT,"tables","mock_readweighted_accuracy.csv"))

cat("\n== Read-weighted mock accuracy (% reads) ==\n")
rw %>% transmute(marker,classifier,rank,
     accuracy=round(accuracy,1), foreign_pct=round(foreign_frac,1)) %>%
  as.data.frame() %>% print(row.names=FALSE)
cat("\nWrote tables/mock_readweighted_accuracy.csv\n")
