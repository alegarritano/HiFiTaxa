#!/usr/bin/env Rscript
# =====================================================================
# Clade exclusion (corrected data): recovery on the FULL reference vs the
# reference with the mock species DEPLETED. Tests genus fallback. Real read
# abundances (feature_table / Emu / full EMITS). Fig-equivalent to the paper's
# clade-exclusion (species collapses; a good classifier holds genus).
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
# folded reported relative abundance per classifier x rank for one result dir
reported_relab <- function(dir, marker) {
  cnt <- read_feature_table(dir)
  bl <- read_blca(dir)$tab; nb <- read_nb(dir)
  out <- list(
    `BLCA|Genus`   = relab_from_counts(bl$id, norm_genus(bl$Genus),   cnt),
    `BLCA|Species` = relab_from_counts(bl$id, norm_species(bl$Species),cnt),
    `NB|Genus`     = relab_from_counts(nb$id, norm_genus(nb$Genus),   cnt),
    `NB|Species`   = relab_from_counts(nb$id, norm_species(nb$Species),cnt))
  if (marker=="16s") {
    out[["Emu|Genus"]]   <- read_emu(dir,"genus")   %>% mutate(taxon=norm_genus(taxon))
    out[["Emu|Species"]] <- read_emu(dir,"species") %>% mutate(taxon=norm_species(taxon))
  } else {
    out[["EMITS|Genus"]]   <- read_emits(dir,"genus")   %>% mutate(taxon=norm_genus(taxon))
    out[["EMITS|Species"]] <- read_emits(dir,"species") %>% mutate(taxon=norm_species(taxon))
  }
  out
}

# --- depletable-subset restriction (local fix; remove after a corrected re-run) -------
# deplete_reference.py removes references by EXACT name match (after GTDB _<UPPER>
# folding only) and does NOT apply SPECIES_ALIAS/GENUS_ALIAS. Four mock members whose
# GTDB/UNITE labels diverge from the mock's names were therefore never removed from the
# "depleted" reference, so their reference was still present and they cannot be validly
# clade-excluded:
#   16S: Escherichia coli (leaked via Escherichia boydii + ECMA0423 sp047199055) and
#        Schaalia odontolytica (leaked via Pauljensenia dentiphila).
#   ITS: Cortinarius ochraceobrunneus (leaked via Phlegmacium ochraceobrunneum).
# We evaluate clade exclusion only over the mock members that were actually depletable.
# After re-running the depletion with benchmark/mock_exclude/*_corrected.txt, delete this
# block (and the two setdiff() lines below) to score the full mock.
UNDEPLETABLE <- list(
  `16s`=list(Genus=c("Escherichia","Schaalia"),
             Species=c("Escherichia coli","Schaalia odontolytica")),
  `its`=list(Genus=c("Cortinarius"),
             Species=c("Cortinarius ochraceobrunneus")))

score_pair <- function(marker, full_dir, excl_dir) {
  eg <- if(marker=="16s") atcc_genus else fung_genus; es <- if(marker=="16s") atcc_species else fung_species
  eg <- setdiff(eg, UNDEPLETABLE[[marker]]$Genus)      # drop taxa the depletion step could not remove
  es <- setdiff(es, UNDEPLETABLE[[marker]]$Species)
  rf <- reported_relab(full_dir,marker); re <- reported_relab(excl_dir,marker)
  bind_rows(lapply(names(rf), function(k){
    p<-str_split_fixed(k,"\\|",2); cl<-p[1]; rk<-p[2]; exp<-if(rk=="Genus") eg else es
    det<-function(tb) length(intersect(exp, tb$taxon[tb$relab>0]))
    tibble(marker=marker, classifier=cl, rank=rk,
           full=det(rf[[k]]), excluded=det(re[[k]]), n_expected=length(exp))
  }))
}

res <- bind_rows(
  score_pair("16s", file.path(BR,"results_atcc_16s"),  file.path(BR,"results_cladeexcl_16s")),
  score_pair("its", file.path(BR,"results_fungi_its"), file.path(BR,"results_cladeexcl_its")))
res$marker <- c("16s"="16S (ATCC MSA-1003)","its"="Fungal ITS (103-sp mock)")[res$marker]
dir.create(file.path(OUT,"tables"),showWarnings=FALSE,recursive=TRUE)
write_csv(res, file.path(OUT,"tables","clade_exclusion.csv"))

cat("\n== Clade exclusion: recovered (full -> species-depleted / expected) ==\n")
res%>%transmute(marker,classifier,rank,ch=sprintf("%d -> %d / %d",full,excluded,n_expected))%>%
  as.data.frame()%>%print(row.names=FALSE)

pl<-res%>%pivot_longer(c(full,excluded),names_to="db",values_to="recovered")
pl$db<-factor(pl$db,levels=c("full","excluded"),labels=c("Full reference","Mock species removed"))
pl$rank<-factor(pl$rank,levels=c("Genus","Species")); pl$classifier<-factor(pl$classifier,levels=c("BLCA","NB","Emu","EMITS"))
expl<-res%>%distinct(marker,rank,n_expected)
p<-ggplot(pl,aes(classifier,recovered,fill=db))+geom_col(position=position_dodge(width=0.75),width=0.7)+
  geom_hline(data=expl,aes(yintercept=n_expected),linetype="dashed",colour="grey30")+
  facet_grid(rank~marker,scales="free")+
  scale_fill_manual(values=c("Full reference"="#377eb8","Mock species removed"="#ff7f00"),name=NULL)+
  labs(x=NULL,y="Mock taxa recovered",title="Clade exclusion: recovery when the mock species is removed from the reference",
       subtitle="Species recovery collapses (removed from DB); genus recovery holds if the classifier backs off gracefully. Dashed = expected.")+
  theme_bw(base_size=12)+theme(legend.position="top",panel.grid.minor=element_blank(),strip.text=element_text(face="bold"))
ggsave(file.path(OUT,"figures","fig_clade_exclusion.png"),p,width=9,height=7.5,dpi=200)
ggsave(file.path(OUT,"figures","fig_clade_exclusion.pdf"),p,width=9,height=7.5)
cat("\nWrote figures/fig_clade_exclusion.{png,pdf} + tables/clade_exclusion.csv\n")
