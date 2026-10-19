#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Figure S4: holdout call decomposition, 16S and fungal ITS.
# Companion to main Fig. 2b/2e (which show accuracy alone). For every rank and
# classifier the calls on the rescuable subset are split into Correct (= the
# per rank accuracy), Incorrect (a wrong taxon) and Unassigned (no call), and
# stacked to 100 percent. This separates a classifier that ERRS from one that
# ABSTAINS: at species BLCA is mostly Correct with a large Incorrect fraction,
# whereas NB is rarely Incorrect but leaves nearly half of sequences Unassigned.
#
# Design: the three outcomes are drawn as SIDE BY SIDE (dodged) bars at each rank,
# so each outcome can be compared directly across ranks and classifiers. The three
# still sum to 100 percent at every rank. Fixed left to right order (Correct,
# Incorrect, Unassigned) gives a position encoding independent of colour. Labels
# are selective: values below 1 percent are suppressed rather than printing zeros.
# Palette Correct/Incorrect/Unassigned is CVD validated (see 00_common.R note).
# Numbers are in supplementary table holdout_accuracy.csv.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
acc <- read_csv(file.path(OUT,"tables","holdout_accuracy.csv"),show_col_types=FALSE,progress=FALSE)

panel <- function(mk, title, leg=TRUE) {
  d <- acc %>% filter(marker==mk) %>%
       transmute(Rank, classifier,
                 Correct=accuracy, Incorrect=pct_incorrect, Unassigned=pct_unassigned) %>%
       pivot_longer(c(Correct,Incorrect,Unassigned), names_to="cat", values_to="pct")
  d$Rank <- factor(d$Rank, levels=RANKS)
  if (mk=="its") levels(d$Rank)[levels(d$Rank)=="Domain"] <- "Kingdom"
  # left to right dodge order: Correct, Incorrect, Unassigned
  d$cat <- factor(d$cat, levels=c("Correct","Incorrect","Unassigned"))
  d$classifier <- factor(d$classifier, levels=c("BLCA","BLCA@80","NB"))

  ggplot(d, aes(Rank, pct, fill=cat)) +
    geom_col(position=position_dodge(0.8), width=0.75) +
    geom_text(aes(label=ifelse(pct >= 1, round(pct), "")),   # suppress near-zero labels
              position=position_dodge(0.8), angle=90, hjust=-0.2, vjust=0.5,
              size=7/.pt, colour="grey20") +
    facet_wrap(~classifier) +
    scale_fill_manual(values=CAT_COL, name=NULL) +
    scale_y_continuous(limits=c(0,116), breaks=seq(0,100,25),
                       expand=expansion(mult=c(0,0.02))) +
    labs(x=NULL, y="% of rescuable test sequences", title=title) +
    theme_ms(12) +
    theme(legend.position=if(leg) "top" else "none", legend.text=element_text(size=10),
          panel.grid=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
}

fig <- panel("16s","a   16S holdout (GTDB r232)", TRUE) /
       panel("its","b   Fungal ITS holdout (UNITE)", FALSE)
ggsave(file.path(OUT,"figures","fig_holdout_decomposition.png"), fig, width=10, height=9, dpi=200)
ggsave(file.path(OUT,"figures","fig_holdout_decomposition.pdf"), fig, width=10, height=9)

cat("\n== Holdout call decomposition (% of rescuable) ==\n")
acc %>% transmute(marker,Rank,classifier,
                  Correct=round(accuracy,1),Incorrect=round(pct_incorrect,1),
                  Unassigned=round(pct_unassigned,1)) %>%
  filter(Rank %in% c("Genus","Species")) %>% as.data.frame() %>% print(row.names=FALSE)
cat("\nWrote figures/fig_holdout_decomposition.{png,pdf}\n")
