#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Figure S6: do BLCA and NB agree, and who is right?
# Alluvial (Sankey) diagram of the holdout, following each test sequence
# from Domain to Species. At every rank each sequence is placed in one of five
# outcomes from comparing the two classifiers to the GTDB/UNITE truth:
#   Both correct            BLCA = NB = truth
#   BLCA correct only       BLCA = truth, NB wrong or unassigned
#   NB correct only         NB = truth, BLCA wrong or unassigned
#   Agree, both wrong       BLCA = NB but neither matches truth (a shared error)
#   Both wrong / unassigned neither matches truth (and they do not agree wrongly)
# Ribbon width is the number of sequences taking that path, so the figure shows
# AGREEMENT (concordant vs split) and CORRECTNESS (who wins when they split) and
# how both degrade toward the species rank.
#
# "Both correct" (green) is anchored at the BASE of every stack. Because that band
# dominates, the y-axis is BROKEN low, through the always-both-correct region
# (16S is >=47.6% both-correct at every rank, ITS >=32.8%), so the green bulk is
# compressed and the smaller outcomes stacked above it are enlarged. The break is
# drawn with ggbreak; ggalluvial's ribbons survive it intact.
#
# Restricted to the species-rescuable subset (truth species present in the 90%
# training fraction), so a correct call is achievable at every rank shown.
# BLCA@80, Emu and EMITS are omitted: this compares the two distinct ASV
# classifiers, and a holdout has no reads for the read level profilers.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
suppressPackageStartupMessages({
  library(ggalluvial); library(ggbreak); library(ggplotify); library(cowplot)
})

CATS <- c("Both correct","BLCA correct only","NB correct only",
          "Agree, both wrong","Both wrong / unassigned")
# green = jointly right, blue = BLCA-only, pink = NB-only, yellow = shared error,
# navy = total failure / abstain
CAT_PAL <- c("Both correct"="#06d6a0","BLCA correct only"="#118ab2",
             "NB correct only"="#ef476f","Agree, both wrong"="#ffd166",
             "Both wrong / unassigned"="#073b4c")

## ---- per-rank BLCA-vs-NB outcome for each rescuable test sequence ----
agree_paths <- function(marker) {
  hd <- file.path(BR, sprintf("benchmark_holdout_%s_seed42", marker)); rd <- file.path(hd,"results")
  truth <- read_lineage_tsv(file.path(hd,"test_10.taxonomy"))
  train <- read_lineage_tsv(file.path(hd,"reference_90.taxonomy"))
  bl <- read_blca(rd)$tab; nb <- read_nb(rd)
  ids <- Reduce(intersect, list(truth$id, bl$id, nb$id))
  o <- function(x) x[match(ids,x$id),]
  T <- o(truth); B <- o(bl); N <- o(nb)
  train_sp <- unique(na.omit(train$Species))
  keep <- !is.na(T$Species) & T$Species %in% train_sp           # species-rescuable subset
  # REVISED 2026-10: decision 4. ITS only: placeholder species labels are not scored at species level.
  if (marker=="its") keep <- keep & !is_placeholder_sp(T$Species)
  T<-T[keep,]; B<-B[keep,]; N<-N[keep,]
  cat_rank <- function(r) {
    tv<-T[[r]]; bv<-B[[r]]; nv<-N[[r]]
    bc <- !is.na(bv) & bv==tv; nc <- !is.na(nv) & nv==tv
    agr <- !is.na(bv) & !is.na(nv) & bv==nv
    out <- rep("Both wrong / unassigned", length(tv))
    out[agr & !bc & !nc] <- "Agree, both wrong"
    out[bc & !nc] <- "BLCA correct only"
    out[nc & !bc] <- "NB correct only"
    out[bc & nc]  <- "Both correct"
    out
  }
  df <- as.data.frame(lapply(RANKS, cat_rank), stringsAsFactors=FALSE); names(df) <- RANKS
  if (marker=="its") names(df)[1] <- "Kingdom"
  df %>% count(across(everything()), name="n")   # aggregate identical paths
}

## ---- long (lodes) form + total count ----
to_lodes <- function(marker) {
  agg <- agree_paths(marker); axes <- setdiff(names(agg), "n")
  n_total <- sum(agg$n); agg$alluvium <- seq_len(nrow(agg))
  lodes <- agg %>% pivot_longer(all_of(axes), names_to="rank", values_to="outcome")
  lodes$rank <- factor(lodes$rank, levels=axes)
  # ggalluvial default reverse=TRUE places the FIRST level at the TOP, so rev()
  # the levels to anchor "Both correct" at the BOTTOM. Legend order is restored
  # by scale_fill_manual(breaks=CATS) below.
  lodes$outcome <- factor(lodes$outcome, levels=rev(CATS))
  list(lodes=lodes, n_total=n_total)
}

## ---- clean white theme (no grid, no grey) ----
theme_sankey <- theme_ms(12) +
  theme(legend.position="none", panel.grid=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x=element_text(face="bold"), plot.title=element_text(face="bold", size=12, hjust=0),
        panel.background=element_rect(fill="white", colour=NA),
        plot.background =element_rect(fill="white", colour=NA),
        axis.text.y.right=element_blank(), axis.ticks.y.right=element_blank(),
        axis.line.y.right=element_blank())

base_plot <- function(marker) {
  d <- to_lodes(marker)
  ggplot(d$lodes, aes(x=rank, stratum=outcome, alluvium=alluvium, y=n, fill=outcome)) +
    geom_flow(width=1/3, alpha=0.45) +                 # NB: colour=NA here suppresses the ribbons
    geom_stratum(width=1/3, colour="grey40", linewidth=0.2) +
    scale_fill_manual(values=CAT_PAL, breaks=CATS, name=NULL) +
    labs(x=NULL, y="% of test sequences") +
    guides(fill=guide_legend(nrow=1))
}

## ---- one panel with a low broken y-axis ----
# brk_pct: percent range removed (must stay below the both-correct floor); SCALES
# gives the expanded upper region this many times the height of the base strip.
SCALES <- 8
make_panel <- function(marker, brk_pct, title) {
  n <- to_lodes(marker)$n_total
  lo <- brk_pct[1]/100*n; hi <- brk_pct[2]/100*n; pos <- seq(0,100,25)/100*n
  base_plot(marker) + ggtitle(title) + theme_sankey +
    scale_y_continuous(name="% of test sequences", breaks=pos, labels=seq(0,100,25),
                       expand=expansion(mult=c(0,0.02))) +
    scale_y_break(c(lo,hi), scales=SCALES, space=0.12, expand=FALSE)
}

n16 <- to_lodes("16s")$n_total; nit <- to_lodes("its")$n_total
pa <- make_panel("16s", c(3,44), sprintf("a   16S holdout (GTDB r232), BLCA vs NB, n = %s species-rescuable sequences", format(n16, big.mark=",")))
pb <- make_panel("its", c(3,29), sprintf("b   Fungal ITS holdout (UNITE), BLCA vs NB, n = %s species-rescuable sequences", format(nit, big.mark=",")))

# one shared legend on top; ggbreak plots must be flattened to grobs (patchwork "/" misbehaves)
leg <- cowplot::get_legend(base_plot("16s") + theme_sankey +
         theme(legend.position="top", legend.text=element_text(size=11)) +
         guides(fill=guide_legend(nrow=1)))
ga <- ggplotify::as.ggplot(print(pa)); gb <- ggplotify::as.ggplot(print(pb))
fig <- cowplot::plot_grid(leg, ga, gb, ncol=1, rel_heights=c(0.07, 1, 1))
ggsave(file.path(OUT,"figures","fig_agreement_sankey.png"), fig, width=12, height=11, dpi=200, bg="white")
ggsave(file.path(OUT,"figures","fig_agreement_sankey.pdf"), fig, width=12, height=11, bg="white")

cat("\n== BLCA vs NB agreement by rank (16S, % of species-rescuable) ==\n")
d16 <- agree_paths("16s")
for (r in c("Family","Genus","Species")) {
  s <- d16 %>% group_by(.data[[r]]) %>% summarise(n=sum(n),.groups="drop") %>% mutate(pct=round(100*n/sum(n),1))
  cat("\n--",r,"--\n"); print(as.data.frame(s), row.names=FALSE)
}
cat("\nWrote figures/fig_agreement_sankey.{png,pdf}\n")
