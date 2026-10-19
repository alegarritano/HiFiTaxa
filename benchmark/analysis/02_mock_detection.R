#!/usr/bin/env Rscript
# =====================================================================
# Mock community detection (corrected data) — Fig. 2b (16S) / 2f (ITS)
# Now with REAL read abundances for every route:
#   16S BLCA/NB : ASV read counts from dada2/feature_table.tsv
#   16S Emu     : emu_{genus,species}_table.tsv
#   ITS BLCA/NB : ASV read counts from dada2/feature_table.tsv (now shipped)
#   ITS EMITS   : full per-sample EMITS abundance table (now complete)
# recovered / missed vs the mock truth; false positive = foreign taxon
# carrying > 0.1% of reads (placeholders excluded). Synonym folding as manuscript.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
score_detect <- function(pred, exp_set, rank) {
  pred <- pred %>% mutate(taxon=if(rank=="Species") norm_species(taxon) else norm_genus(taxon)) %>%
          filter(!is.na(taxon)) %>% group_by(taxon) %>% summarise(relab=sum(relab),.groups="drop")
  detected  <- pred$taxon[pred$relab>0]
  recovered <- intersect(exp_set, detected); missed <- setdiff(exp_set, detected)
  foreign <- pred %>% filter(!(taxon %in% exp_set), relab>FP_THRESH)
  if (rank=="Species") foreign <- foreign %>% filter(!is_placeholder_sp(taxon)) else
                       foreign <- foreign %>% filter(!str_detect(taxon,"(ceae$| sp$|^Fungi|Incertae|_fam_|_gen_)"))
  # precompute the names BEFORE the tibble: inside tibble() the column `missed=length(missed)`
  # shadows the character vector, so `paste(sort(missed))` there would stringify the count.
  missed_str <- paste(sort(missed), collapse="; ")
  tibble(rank=rank, n_expected=length(exp_set), recovered=length(recovered), missed=length(missed),
         false_positive=nrow(foreign), missed_taxa=missed_str,
         fp_taxa=paste(sort(foreign$taxon),collapse="; "),
         fp_relab=paste(sprintf("%s(%.1f%%)",foreign$taxon,100*foreign$relab),collapse="; "))
}

run_marker <- function(marker) {
  if (marker=="16s") {
    dir<-file.path(BR,"results_atcc_16s"); cnt<-read_feature_table(dir)
    bl<-read_blca(dir)$tab; nb<-read_nb(dir)
    res<-bind_rows(
      score_detect(relab_from_counts(bl$id,norm_genus(bl$Genus),cnt),   atcc_genus,  "Genus")  %>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(bl$id,norm_species(bl$Species),cnt),atcc_species,"Species")%>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(nb$id,norm_genus(nb$Genus),cnt),   atcc_genus,  "Genus")  %>% mutate(classifier="NB"),
      score_detect(relab_from_counts(nb$id,norm_species(nb$Species),cnt),atcc_species,"Species")%>% mutate(classifier="NB"),
      score_detect(read_emu(dir,"genus"),                              atcc_genus,  "Genus")  %>% mutate(classifier="Emu"),
      score_detect(read_emu(dir,"species"),                            atcc_species,"Species") %>% mutate(classifier="Emu"))
    res$marker<-"16S (ATCC MSA-1003)"; res$readlevel<-res$classifier=="Emu"
  } else {
    dir<-file.path(BR,"results_fungi_its"); cnt<-read_feature_table(dir)
    bl<-read_blca(dir)$tab; nb<-read_nb(dir)
    res<-bind_rows(
      score_detect(relab_from_counts(bl$id,norm_genus(bl$Genus),cnt),   fung_genus,  "Genus")  %>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(bl$id,norm_species(bl$Species),cnt),fung_species,"Species")%>% mutate(classifier="BLCA"),
      score_detect(relab_from_counts(nb$id,norm_genus(nb$Genus),cnt),   fung_genus,  "Genus")  %>% mutate(classifier="NB"),
      score_detect(relab_from_counts(nb$id,norm_species(nb$Species),cnt),fung_species,"Species")%>% mutate(classifier="NB"),
      score_detect(read_emits(dir,"genus"),                            fung_genus,  "Genus")  %>% mutate(classifier="EMITS"),
      score_detect(read_emits(dir,"species"),                          fung_species,"Species") %>% mutate(classifier="EMITS"))
    res$marker<-"Fungal ITS (103-sp mock)"; res$readlevel<-res$classifier=="EMITS"
  }
  res
}

det <- bind_rows(run_marker("16s"), run_marker("its"))
dir.create(file.path(OUT,"tables"),showWarnings=FALSE,recursive=TRUE)
write_csv(det%>%select(marker,classifier,rank,n_expected,recovered,missed,false_positive,missed_taxa,fp_relab),
          file.path(OUT,"tables",paste0("mock_detection", if (FP_THRESH == 0.001) "" else sprintf("_threshold_%g", FP_THRESH), ".csv")))

cat("\n== Mock detection (recovered/expected, +FP) ==\n")
det%>%transmute(marker,classifier,rank,s=sprintf("%2d/%d rec, %d FP",recovered,n_expected,false_positive))%>%
  as.data.frame()%>%print(row.names=FALSE)

pl<-det%>%select(marker,classifier,rank,recovered,missed,false_positive)%>%
  pivot_longer(c(recovered,missed,false_positive),names_to="cat",values_to="n")
pl$cat<-factor(pl$cat,levels=c("recovered","missed","false_positive"),labels=c("Recovered","Missed","False positive"))
pl$classifier<-factor(pl$classifier,levels=c("BLCA","NB","Emu","EMITS")); pl$rank<-factor(pl$rank,levels=c("Genus","Species"))
expl<-det%>%distinct(marker,rank,n_expected)
p<-ggplot(pl,aes(classifier,n,fill=cat))+geom_col(width=0.72,position=position_stack(reverse=TRUE))+
  geom_hline(data=expl,aes(yintercept=n_expected),linetype="dashed",colour="grey30")+
  facet_grid(rank~marker,scales="free")+
  scale_fill_manual(values=c("Recovered"="#4daf4a","Missed"="#bdbdbd","False positive"="#e41a1c"),name=NULL)+
  labs(x=NULL,y="Number of taxa",title="Mock community detection by classifier",
       subtitle="Dashed = expected taxa; false positives (>0.1% reads) stack on top. Read counts from feature table / Emu / EMITS.")+
  theme_bw(base_size=12)+theme(legend.position="top",panel.grid.minor=element_blank(),strip.text=element_text(face="bold"))
ggsave(file.path(OUT,"figures","fig_mock_detection.png"),p,width=9,height=7.5,dpi=200)
ggsave(file.path(OUT,"figures","fig_mock_detection.pdf"),p,width=9,height=7.5)
cat("\nWrote figures/fig_mock_detection.{png,pdf} + tables/mock_detection.csv\n")
