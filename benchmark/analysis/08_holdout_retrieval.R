#!/usr/bin/env Rscript
# =====================================================================
# Holdout RETRIEVAL barplots, 16S and ITS. Follows manuscript Fig. 2c/d (16S)
# and 2g/h (ITS):
#   (c/g) per-rank overall accuracy, dodged bars (BLCA, BLCA@80, NB), value labels
#   (d/h) orphan novelty backoff, stacked Correct/Incorrect/Unassigned, per
#         classifier, faceted "Orphan species (genus call)" | "Orphan genera (family call)"
# Reads the verified tables written by 01_holdout_accuracy.R.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
acc <- read_csv(file.path(OUT,"tables","holdout_accuracy.csv"), show_col_types=FALSE, progress=FALSE)
bo  <- read_csv(file.path(OUT,"tables","holdout_orphan_backoff.csv"), show_col_types=FALSE, progress=FALSE)
acc$Rank <- factor(acc$Rank, levels=RANKS)
acc$classifier <- factor(acc$classifier, levels=c("BLCA","BLCA@80","NB"))

rank_lab <- function(marker) if(marker=="its") c("Domain"="Kingdom") else c()  # ITS uses Kingdom for Domain

panel_acc <- function(mk, tag) {
  d <- acc %>% filter(marker==mk)
  if (mk=="its") levels(d$Rank)[levels(d$Rank)=="Domain"] <- "Kingdom"
  ggplot(d, aes(Rank, accuracy, fill=classifier)) +
    geom_col(position=position_dodge(0.8), width=0.75, colour="black", linewidth=0.2) +
    geom_text(aes(label=round(accuracy)), position=position_dodge(0.8),
              angle=90, hjust=-0.12, size=3.4) +
    scale_fill_manual(values=CLF_COL, name=NULL) +
    scale_y_continuous(limits=c(0,113), breaks=seq(0,100,25), oob=scales::oob_keep) +
    labs(x=NULL, y="Accuracy (% correct)", title=tag) +
    theme_ms(12) + theme(legend.position="top", axis.text.x=element_text(angle=30,hjust=1))
}
panel_backoff <- function(mk, tag) {
  d <- bo %>% filter(marker==mk) %>%
       pivot_longer(c(Correct,Incorrect,Unassigned),names_to="cat",values_to="pct")
  d$cat <- factor(d$cat, levels=c("Unassigned","Incorrect","Correct"))
  d$classifier <- factor(d$classifier, levels=c("BLCA","BLCA@80","NB"))
  d$level <- factor(d$level, levels=c("Orphan species -> genus","Orphan genus -> family"),
                    labels=c("Orphan species (genus call)","Orphan genera (family call)"))
  ggplot(d, aes(classifier, pct, fill=cat)) +
    geom_col(width=0.7) + facet_wrap(~level) +
    scale_fill_manual(values=CAT_COL, name=NULL,
                      breaks=c("Correct","Incorrect","Unassigned")) +
    scale_y_continuous(expand=expansion(mult=c(0,0.02))) +
    labs(x=NULL, y="% of orphan sequences", title=tag) +
    theme_ms(12) + theme(legend.position="top")
}

fig16 <- panel_acc("16s","a   16S holdout: per-rank accuracy") /
         panel_backoff("16s","b   16S holdout: novelty backoff")
figIT <- panel_acc("its","a   ITS holdout: per-rank accuracy") /
         panel_backoff("its","b   ITS holdout: novelty backoff")
ggsave(file.path(OUT,"figures","fig_holdout_retrieval_16s.png"), fig16, width=9, height=8, dpi=200)
ggsave(file.path(OUT,"figures","fig_holdout_retrieval_16s.pdf"), fig16, width=9, height=8)
ggsave(file.path(OUT,"figures","fig_holdout_retrieval_its.png"), figIT, width=9, height=8, dpi=200)
ggsave(file.path(OUT,"figures","fig_holdout_retrieval_its.pdf"), figIT, width=9, height=8)
cat("Wrote fig_holdout_retrieval_{16s,its}.{png,pdf}\n")
