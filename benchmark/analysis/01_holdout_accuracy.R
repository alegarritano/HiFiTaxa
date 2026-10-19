#!/usr/bin/env Rscript
# =====================================================================
# Holdout accuracy + orphan novelty backoff (corrected full holdouts)
# Reproduces manuscript Fig. 2c/d (16S GTDB) and 2g/h (ITS UNITE).
#
# Rescuable subset (manuscript definition): a test sequence counts at rank r
# only if its TRUE taxon at r is present in the 90% training fraction
# (reference_90.taxonomy). Accuracy at r = correct / rescuable-at-r.
#
# Orphan novelty backoff:
#   orphan species (species absent from training, genus present) -> score genus
#   orphan genera  (genus absent from training,  family present) -> score family
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
# unique taxa present in the 90% training fraction, per rank
train_rank_sets <- function(ref90_path) {
  ref <- read_lineage_tsv(ref90_path)
  setNames(lapply(RANKS, function(r) unique(na.omit(ref[[r]]))), RANKS)
}

score_marker <- function(marker) {
  hd   <- file.path(BR, sprintf("benchmark_holdout_%s_seed42", marker))
  rdir <- file.path(hd, "results")
  truth <- read_lineage_tsv(file.path(hd,"test_10.taxonomy"))
  train <- train_rank_sets(file.path(hd,"reference_90.taxonomy"))
  bl <- read_blca(rdir, path_conf=NULL); nb <- read_nb(rdir)

  ids <- intersect(truth$id, bl$tab$id)
  T  <- truth   %>% filter(id %in% ids) %>% arrange(id)
  B  <- bl$tab  %>% filter(id %in% ids) %>% arrange(id)
  Bc <- bl$conf %>% filter(id %in% ids) %>% arrange(id)
  N  <- nb      %>% filter(id %in% ids) %>% arrange(id)
  stopifnot(all(T$id==B$id), all(T$id==Bc$id), all(T$id==N$id))

  # REVISED 2026-10: decision 4 (fungal ITS holdout ONLY; the 16S holdout is unchanged because GTDB
  # placeholder names are genome defined species). Test sequences whose UNITE species label is a
  # placeholder ("Tomentella sp", "Fungi sp", "GS01 phy Incertae sedis sp") are excluded from every
  # species level calculation, and those whose genus label is of the "Incertae sedis" kind
  # ("Agaricales_gen_Incertae_sedis") from every genus level calculation, including the two
  # absent taxon (orphan) analyses. excl[[r]] is all FALSE for the other ranks and for 16S.
  excl <- setNames(lapply(RANKS, function(r) rep(FALSE, nrow(T))), RANKS)
  if (marker=="its") {
    excl$Species <- !is.na(T$Species) & is_placeholder_sp(T$Species)
    excl$Genus   <- !is.na(T$Genus)   & is_placeholder_genus(T$Genus)
  }

  # BLCA@80: keep a rank's call only if its bootstrap >= 80
  Bmask <- B
  for (r in RANKS) { cc<-suppressWarnings(as.numeric(Bc[[paste0(r,"_conf")]])); cc[is.na(cc)]<-0
                     Bmask[[r]][cc<80] <- NA }

  # per-rank accuracy + decomposition (correct/incorrect/unassigned) on the rescuable subset
  acc <- bind_rows(lapply(RANKS, function(r){
    tv <- T[[r]]
    keep <- !is.na(tv) & (tv %in% train[[r]]) & !excl[[r]]   # REVISED 2026-10: decision 4 (& !excl)
    denom <- sum(keep)
    bind_rows(lapply(c("BLCA","BLCA@80","NB"), function(cl){
      pv <- switch(cl,"BLCA"=B[[r]],"BLCA@80"=Bmask[[r]],"NB"=N[[r]])
      cc <- sum(keep & !is.na(pv) & pv==tv); ic <- sum(keep & !is.na(pv) & pv!=tv); un <- sum(keep & is.na(pv))
      tibble(marker=marker, Rank=r, classifier=cl, n=denom, correct=cc,
             accuracy=ifelse(denom>0,100*cc/denom,NA_real_),
             pct_incorrect=ifelse(denom>0,100*ic/denom,NA_real_),
             pct_unassigned=ifelse(denom>0,100*un/denom,NA_real_))
    }))
  }))

  # BLCA@80 coverage / precision on the rescuable subset
  covprec <- bind_rows(lapply(RANKS, function(r){
    tv<-T[[r]]; keep<-!is.na(tv)&(tv %in% train[[r]])&!excl[[r]]; pv<-Bmask[[r]]; kept<-keep&!is.na(pv)   # REVISED 2026-10: decision 4 (&!excl)
    tibble(marker=marker, Rank=r, coverage=100*sum(kept)/sum(keep),
           precision=ifelse(sum(kept)>0,100*sum(kept&pv==tv)/sum(kept),NA_real_))
  }))

  # novelty backoff. species-orphan: genus absent? no -> score genus. genus-orphan: family present -> score family.
  # REVISED 2026-10: decision 4. The three definitions below drop placeholder species labels and
  # "Incertae sedis" genus labels (ITS only; excl is all FALSE for 16S). The *_all versions keep the
  # old definitions so that the old and new denominators can be reported side by side.
  sp_orphan_all <- !is.na(T$Species) & !(T$Species %in% train$Species)
  g_present_all <- !is.na(T$Genus)   &  (T$Genus   %in% train$Genus)
  g_orphan_all  <- !is.na(T$Genus)   & !(T$Genus   %in% train$Genus)
  sp_orphan <- sp_orphan_all & !excl$Species
  g_present <- g_present_all & !excl$Genus
  g_orphan  <- g_orphan_all  & !excl$Genus
  f_present <- !is.na(T$Family)  &  (T$Family  %in% train$Family)
  backoff <- bind_rows(
    bind_rows(lapply(c("BLCA","BLCA@80","NB"), function(cl){
      idx <- sp_orphan & g_present; pv<-switch(cl,"BLCA"=B$Genus,"BLCA@80"=Bmask$Genus,"NB"=N$Genus)[idx]; tv<-T$Genus[idx]; n<-length(tv)
      tibble(marker=marker,classifier=cl,level="Orphan species -> genus",n=n,
             Correct=100*sum(!is.na(pv)&pv==tv)/n,Incorrect=100*sum(!is.na(pv)&pv!=tv)/n,Unassigned=100*sum(is.na(pv))/n)})),
    bind_rows(lapply(c("BLCA","BLCA@80","NB"), function(cl){
      idx <- g_orphan & f_present; pv<-switch(cl,"BLCA"=B$Family,"BLCA@80"=Bmask$Family,"NB"=N$Family)[idx]; tv<-T$Family[idx]; n<-length(tv)
      tibble(marker=marker,classifier=cl,level="Orphan genus -> family",n=n,
             Correct=100*sum(!is.na(pv)&pv==tv)/n,Incorrect=100*sum(!is.na(pv)&pv!=tv)/n,Unassigned=100*sum(is.na(pv))/n)})))
  # REVISED 2026-10: decision 4. Old and new denominators of every affected calculation.
  f_present <- !is.na(T$Family) & (T$Family %in% train$Family)
  sp_resc_all <- !is.na(T$Species) & (T$Species %in% train$Species)
  denom <- tibble(marker=marker,
    quantity=c("Test sequences scored",
               "Test sequences with a placeholder species label",
               "Test sequences with an Incertae sedis genus label",
               "Species accuracy denominator (species present in training)",
               "Genus accuracy denominator (genus present in training)",
               "Species absent, genus present (scored at genus)",
               "Genus absent, family present (scored at family)"),
    n_before=c(length(ids), NA, NA, sum(sp_resc_all), sum(g_present_all),
               sum(sp_orphan_all & g_present_all), sum(g_orphan_all & f_present)),
    n_after =c(length(ids), sum(excl$Species), sum(excl$Genus), sum(sp_resc_all & !excl$Species), sum(g_present),
               sum(sp_orphan & g_present), sum(g_orphan & f_present)))
  list(acc=acc, covprec=covprec, backoff=backoff, n_scored=length(ids),
       n_sp_orphan=sum(sp_orphan), n_g_orphan=sum(g_orphan), denom=denom)
}

s16 <- score_marker("16s"); sits <- score_marker("its")
message(sprintf("16S scored %d (sp-orphan %d, genus-orphan %d) | ITS scored %d (sp-orphan %d, genus-orphan %d)",
  s16$n_scored,s16$n_sp_orphan,s16$n_g_orphan, sits$n_scored,sits$n_sp_orphan,sits$n_g_orphan))

acc<-bind_rows(s16$acc,sits$acc); covprec<-bind_rows(s16$covprec,sits$covprec); backoff<-bind_rows(s16$backoff,sits$backoff)
mlab<-c("16s"="16S (GTDB r232 holdout)","its"="Fungal ITS (UNITE holdout)")
acc$marker_lab<-mlab[acc$marker]; backoff$marker_lab<-mlab[backoff$marker]
acc$Rank<-factor(acc$Rank,levels=RANKS); acc$classifier<-factor(acc$classifier,levels=c("BLCA","BLCA@80","NB"))

dir.create(file.path(OUT,"tables"),showWarnings=FALSE,recursive=TRUE)
write_csv(acc,     file.path(OUT,"tables","holdout_accuracy.csv"))
write_csv(covprec, file.path(OUT,"tables","holdout_blca80_coverage_precision.csv"))
write_csv(backoff, file.path(OUT,"tables","holdout_orphan_backoff.csv"))
# REVISED 2026-10: decision 4. denominators before and after the placeholder exclusion (ITS only)
write_csv(sits$denom, file.path(OUT,"tables","holdout_its_placeholder_denominators.csv"))
cat("\n== ITS holdout denominators before and after excluding placeholder labels ==\n")
print(as.data.frame(sits$denom), row.names=FALSE)

pal<-c("BLCA"="#1b9e77","BLCA@80"="#7570b3","NB"="#d95f02")
p_acc<-ggplot(acc,aes(Rank,accuracy,colour=classifier,group=classifier))+
  geom_line(linewidth=0.9)+geom_point(size=2.2)+facet_wrap(~marker_lab)+
  scale_colour_manual(values=pal,name="Classifier")+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))+
  labs(x=NULL,y="Overall accuracy (%)",title="Per-rank accuracy on the leave-10%-out reference holdout",
       subtitle="Rescuable subset (true taxon at that rank present in the 90% training fraction)")+
  theme_bw(base_size=12)+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="top",panel.grid.minor=element_blank())
bl<-backoff%>%pivot_longer(c(Correct,Incorrect,Unassigned),names_to="outcome",values_to="pct")
bl$outcome<-factor(bl$outcome,levels=c("Correct","Incorrect","Unassigned")); bl$classifier<-factor(bl$classifier,levels=c("BLCA","BLCA@80","NB"))
p_bo<-ggplot(bl,aes(classifier,pct,fill=outcome))+geom_col(width=0.7)+facet_grid(level~marker_lab)+
  scale_fill_manual(values=c("Correct"="#4daf4a","Incorrect"="#e41a1c","Unassigned"="#999999"),name=NULL)+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25))+
  labs(x=NULL,y="% of orphan queries",title="Novelty backoff: nearest-rank call for taxa absent from the reference")+
  theme_bw(base_size=12)+theme(legend.position="top",panel.grid.minor=element_blank(),strip.text.y=element_text(size=8))
fig<-p_acc/p_bo+patchwork::plot_annotation(tag_levels="a")+patchwork::plot_layout(heights=c(1,1.2))
ggsave(file.path(OUT,"figures","fig_holdout_accuracy.png"),fig,width=10,height=10,dpi=200)
ggsave(file.path(OUT,"figures","fig_holdout_accuracy.pdf"),fig,width=10,height=10)

cat("\n== Holdout accuracy (%) — rescuable subset ==\n")
acc%>%select(marker,Rank,classifier,accuracy)%>%pivot_wider(names_from=classifier,values_from=accuracy)%>%
  mutate(across(where(is.numeric),~round(.,1)))%>%as.data.frame()%>%print(row.names=FALSE)
cat("\n== Orphan backoff (% correct genus/family) ==\n")
backoff%>%select(marker,level,classifier,Correct)%>%mutate(Correct=round(Correct,1))%>%
  pivot_wider(names_from=classifier,values_from=Correct)%>%as.data.frame()%>%print(row.names=FALSE)
cat("\nWrote figures/fig_holdout_accuracy.{png,pdf} + tables/holdout_*.csv\n")
