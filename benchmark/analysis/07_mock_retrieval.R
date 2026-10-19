#!/usr/bin/env Rscript
# =====================================================================
# Mock RETRIEVAL barplots, 16S and ITS. Follows manuscript Fig. 2a/b (16S)
# and 2e/f (ITS):
#   (a/e) # non-redundant taxa reported per classifier vs Expected, genus+species
#   (b/f) detection split into Detected / Missed / False pos., dodged by
#         classifier (BLCA, BLCA@80, Emu|EMITS, NB), genus+species
# Read-weighted (feature table / Emu / EMITS). Manuscript classifier colours.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
detect <- function(relab, exp_set, rank) {
  d <- relab %>% mutate(taxon=if(rank=="Species") norm_species(taxon) else norm_genus(taxon)) %>%
       filter(!is.na(taxon)) %>% group_by(taxon) %>% summarise(relab=sum(relab),.groups="drop")
  det <- length(intersect(exp_set, d$taxon[d$relab>0]))
  fo <- d %>% filter(!(taxon %in% exp_set), relab>FP_THRESH)
  fo <- if (rank=="Species") fo %>% filter(!is_placeholder_sp(taxon)) else
        fo %>% filter(!str_detect(taxon,"(ceae$| sp$|^Fungi|Incertae|_fam_|_gen_)"))
  tibble(rank=rank, Detected=det, Missed=length(exp_set)-det, `False pos.`=nrow(fo),
         nonred=det+nrow(fo), n_expected=length(exp_set))
}

marker_detect <- function(marker) {
  if (marker=="16s") { dir<-file.path(BR,"results_atcc_16s"); eg<-atcc_genus; es<-atcc_species; rl<-"Emu" }
  else                { dir<-file.path(BR,"results_fungi_its"); eg<-fung_genus; es<-fung_species; rl<-"EMITS" }
  cnt<-read_feature_table(dir); bl<-read_blca(dir)$tab; b8<-read_blca_at80(dir); nb<-read_nb(dir)
  rd <- function(tab,col) relab_from_counts(tab$id, tab[[col]], cnt)
  em <- function(rk) if(marker=="16s") read_emu(dir,tolower(rk)) else read_emits(dir,tolower(rk))
  rows <- bind_rows(
    detect(rd(bl,"Genus"),    eg,"Genus")   %>% mutate(classifier="BLCA"),
    detect(rd(bl,"Species"),  es,"Species") %>% mutate(classifier="BLCA"),
    detect(rd(b8,"Genus"),    eg,"Genus")   %>% mutate(classifier="BLCA@80"),
    detect(rd(b8,"Species"),  es,"Species") %>% mutate(classifier="BLCA@80"),
    detect(em("Genus"),       eg,"Genus")   %>% mutate(classifier=rl),
    detect(em("Species"),     es,"Species") %>% mutate(classifier=rl),
    detect(rd(nb,"Genus"),    eg,"Genus")   %>% mutate(classifier="NB"),
    detect(rd(nb,"Species"),  es,"Species") %>% mutate(classifier="NB"))
  rows$marker <- marker; rows$readlevel <- rl; rows
}

d16 <- marker_detect("16s"); dIT <- marker_detect("its")
# Output names carry the threshold when FP_THRESH is not the manuscript default of 0.001, so that
# the Supplementary Table S11 counts (FP_THRESH=0 and FP_THRESH=0.0001) do not overwrite Table S2.
THR_TAG <- if (FP_THRESH == 0.001) "" else if (FP_THRESH == 0) "_threshold_any_read" else sprintf("_threshold_%gpct", 100 * FP_THRESH)
write_csv(bind_rows(d16,dIT) %>% select(marker,classifier,rank,Detected,Missed,`False pos.`,nonred,n_expected),
          file.path(OUT,"tables",paste0("mock_retrieval",THR_TAG,".csv")))

# --- panel builders (manuscript style) ---
mlevels <- function(rl) c("Expected","BLCA","BLCA@80",rl,"NB")
panel_nonred <- function(dd, rl, tag) {
  exp <- dd %>% distinct(rank,n_expected) %>% mutate(classifier="Expected",nonred=n_expected)
  d <- bind_rows(dd %>% select(rank,classifier,nonred), exp %>% select(rank,classifier,nonred))
  d$classifier <- factor(d$classifier, levels=mlevels(rl)); d$rank<-factor(d$rank,levels=c("Genus","Species"))
  ggplot(d, aes(classifier,nonred,fill=classifier)) +
    geom_col(width=0.8, colour="black", linewidth=0.2) +
    geom_text(aes(label=nonred), vjust=-0.3, size=3) +
    facet_wrap(~rank, scales="free_y") +
    scale_fill_manual(values=CLF_COL, guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
    labs(x=NULL, y="# non-redundant taxa", title=tag) +
    theme_ms(12) + theme(axis.text.x=element_text(angle=25,hjust=1))
}
panel_detect <- function(dd, rl, tag) {
  d <- dd %>% select(rank,classifier,Detected,Missed,`False pos.`) %>%
       pivot_longer(c(Detected,Missed,`False pos.`),names_to="cat",values_to="n")
  d$cat<-factor(d$cat,levels=c("Detected","Missed","False pos."))
  d$classifier<-factor(d$classifier,levels=c("BLCA","BLCA@80",rl,"NB")); d$rank<-factor(d$rank,levels=c("Genus","Species"))
  ggplot(d, aes(cat,n,fill=classifier)) +
    geom_col(position=position_dodge(0.85), width=0.8) +
    geom_text(aes(label=n), position=position_dodge(0.85), vjust=-0.3, size=2.6) +
    facet_wrap(~rank, scales="free_y") +
    scale_fill_manual(values=CLF_COL, name=NULL) +
    scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
    labs(x=NULL, y="# taxa", title=tag) +
    theme_ms(12) + theme(legend.position="top")
}

fig16 <- panel_nonred(d16,"Emu","a   16S mock: non-redundant taxa") /
         panel_detect(d16,"Emu","b   16S mock: detection vs expected")
figIT <- panel_nonred(dIT,"EMITS","a   ITS mock: non-redundant taxa") /
         panel_detect(dIT,"EMITS","b   ITS mock: detection vs expected")
ggsave(file.path(OUT,"figures",paste0("fig_mock_retrieval_16s",THR_TAG,".png")), fig16, width=9, height=8, dpi=200)
ggsave(file.path(OUT,"figures",paste0("fig_mock_retrieval_16s",THR_TAG,".pdf")), fig16, width=9, height=8)
ggsave(file.path(OUT,"figures",paste0("fig_mock_retrieval_its",THR_TAG,".png")), figIT, width=9, height=8, dpi=200)
ggsave(file.path(OUT,"figures",paste0("fig_mock_retrieval_its",THR_TAG,".pdf")), figIT, width=9, height=8)

cat("\n== Mock retrieval ==\n")
bind_rows(d16,dIT) %>% transmute(marker,classifier,rank,
   s=sprintf("det %d / miss %d / FP %d (nonred %d/%d)",Detected,Missed,`False pos.`,nonred,n_expected)) %>%
  as.data.frame() %>% print(row.names=FALSE)
cat(sprintf("\nWrote fig_mock_retrieval_{16s,its}%s.{png,pdf} + tables/mock_retrieval%s.csv (FP_THRESH = %g)\n", THR_TAG, THR_TAG, FP_THRESH))
