#!/usr/bin/env Rscript
# =====================================================================
# Mock community TAXAPLOTS (composition), 16S and fungal ITS.
# Follows the manuscript pattern (Suppl. S1 / composition panels): stacked
# read-weighted relative-abundance bars, x = Expected + classifiers, top-14
# taxa in the Okabe-Ito palette + Other (light grey) + Unassigned (dark grey),
# genus and species panels. 16S classifiers BLCA/BLCA@80/Emu/NB; ITS
# BLCA/BLCA@80/NB/EMITS. ITS "Expected" is the equal-volume pool reference
# (each of the 103 mock species weighted equally, collapsed to the rank).
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
build_comp <- function(methods, expected=NULL, ntop=14) {
  all <- bind_rows(lapply(names(methods), function(m) methods[[m]] %>% mutate(method=m)))
  if (!is.null(expected)) all <- bind_rows(expected %>% mutate(method="Expected"), all)
  rankbase <- if (!is.null(expected)) expected else
              all %>% group_by(taxon) %>% summarise(relab=max(relab),.groups="drop")
  rankbase <- rankbase %>% filter(!taxon %in% c("Unassigned","Other")) %>% arrange(desc(relab))
  top <- head(rankbase$taxon, ntop)
  all <- all %>% mutate(grp=ifelse(taxon %in% top, taxon,
                                   ifelse(taxon=="Unassigned","Unassigned","Other"))) %>%
    group_by(method,grp) %>% summarise(relab=sum(relab),.groups="drop")
  lev <- c(top,"Other","Unassigned")
  all$grp <- factor(all$grp, levels=lev)
  cols <- setNames(c(TAXA_PAL[seq_along(top)],OTHER_COL,UNASSIGNED_COL), lev)
  list(df=all, cols=cols)
}

taxaplot <- function(bc, method_levels, title, ylab) {
  d <- bc$df; d$method <- factor(d$method, levels=method_levels)
  ggplot(d, aes(method, relab, fill=grp)) +
    geom_col(width=0.72, colour="white", linewidth=0.2) +
    scale_fill_manual(values=bc$cols, name=NULL) +
    scale_y_continuous(labels=percent_format(accuracy=1), expand=expansion(mult=c(0,0.01))) +
    labs(x=NULL, y=ylab, title=title) +
    guides(fill=guide_legend(ncol=3, byrow=TRUE)) +
    theme_ms(12) + theme(legend.position="bottom", legend.text=element_text(face="italic",size=8),
                         legend.key.size=unit(10,"pt"), axis.text.x=element_text(size=10))
}

## ------------------- 16S mock -------------------
d16 <- file.path(BR,"results_atcc_16s"); c16 <- read_feature_table(d16)
bl<-read_blca(d16)$tab; b8<-read_blca_at80(d16); nb<-read_nb(d16)
exp_g <- tibble(species=names(ATCC_TRUTH),relab=as.numeric(ATCC_TRUTH)) %>%
  mutate(taxon=norm_genus(str_split_fixed(norm_species(species)," ",2)[,1])) %>%
  group_by(taxon) %>% summarise(relab=sum(relab),.groups="drop")
exp_s <- tibble(taxon=norm_species(names(ATCC_TRUTH)),relab=as.numeric(ATCC_TRUTH))
m16g <- list(BLCA=comp_reads(bl$id,norm_genus(bl$Genus),c16),
             `BLCA@80`=comp_reads(b8$id,norm_genus(b8$Genus),c16),
             Emu=read_emu(d16,"genus")%>%mutate(taxon=norm_genus(taxon))%>%group_by(taxon)%>%summarise(relab=sum(relab),.groups="drop"),
             NB=comp_reads(nb$id,norm_genus(nb$Genus),c16))
m16s <- list(BLCA=comp_reads(bl$id,norm_species(bl$Species),c16),
             `BLCA@80`=comp_reads(b8$id,norm_species(b8$Species),c16),
             Emu=read_emu(d16,"species")%>%mutate(taxon=norm_species(taxon))%>%group_by(taxon)%>%summarise(relab=sum(relab),.groups="drop"),
             NB=comp_reads(nb$id,norm_species(nb$Species),c16))
lev16 <- c("Expected","BLCA","BLCA@80","Emu","NB")
pA <- taxaplot(build_comp(m16g,exp_g), lev16, "a   16S mock: genus",  "Mean relative abundance")
pB <- taxaplot(build_comp(m16s,exp_s), lev16, "b   16S mock: species","Mean relative abundance")
f16 <- pA + pB
ggsave(file.path(OUT,"figures","fig_mock_taxaplot_16s.png"), f16, width=13, height=7, dpi=200)
ggsave(file.path(OUT,"figures","fig_mock_taxaplot_16s.pdf"), f16, width=13, height=7)

## ------------------- ITS mock -------------------
# The 103-species fungal mock is EQUAL-VOLUME pooled, so there is no abundance
# ground truth. We therefore drop the (meaningless) equal-'Expected' column and
# show each classifier's observed composition, ranked by observed abundance so
# the taxa the classifiers actually report get coloured (not an arbitrary set).
dI <- file.path(BR,"results_fungi_its"); cI <- read_feature_table(dI)
blI<-read_blca(dI)$tab; b8I<-read_blca_at80(dI); nbI<-read_nb(dI)
mItg <- list(BLCA=comp_reads(blI$id,norm_genus(blI$Genus),cI),
             `BLCA@80`=comp_reads(b8I$id,norm_genus(b8I$Genus),cI),
             NB=comp_reads(nbI$id,norm_genus(nbI$Genus),cI),
             EMITS=read_emits(dI,"genus")%>%mutate(taxon=norm_genus(taxon))%>%group_by(taxon)%>%summarise(relab=sum(relab),.groups="drop"))
mIts <- list(BLCA=comp_reads(blI$id,norm_species(blI$Species),cI),
             `BLCA@80`=comp_reads(b8I$id,norm_species(b8I$Species),cI),
             NB=comp_reads(nbI$id,norm_species(nbI$Species),cI),
             EMITS=read_emits(dI,"species")%>%mutate(taxon=norm_species(taxon))%>%group_by(taxon)%>%summarise(relab=sum(relab),.groups="drop"))
levI <- c("BLCA","BLCA@80","NB","EMITS")
pE <- taxaplot(build_comp(mItg), levI, "a   ITS mock: genus",  "Mean relative abundance")
pF <- taxaplot(build_comp(mIts), levI, "b   ITS mock: species","Mean relative abundance")
fIT <- pE + pF
ggsave(file.path(OUT,"figures","fig_mock_taxaplot_its.png"), fIT, width=13, height=7.5, dpi=200)
ggsave(file.path(OUT,"figures","fig_mock_taxaplot_its.pdf"), fIT, width=13, height=7.5)

cat("Wrote fig_mock_taxaplot_16s.{png,pdf} and fig_mock_taxaplot_its.{png,pdf}\n")
