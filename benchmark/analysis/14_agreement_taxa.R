#!/usr/bin/env Rscript
# =====================================================================
# Which 16S GENERA do BLCA and NB agree-and-get-right, and which do they contest?
# Companion to Fig S6 (which aggregates the outcomes per rank): here the same
# outcomes are resolved PER TAXON at the genus rank, so the figure names the
# genera that drive the BLCA-vs-NB divergence rather than just its size.
#
# The holdout is split into two populations, reported separately:
#   IN REFERENCE  the true species IS in the 90% training fraction (n = 90,688).
#                 The classifier has seen the species; the genus is trivially present.
#   ORPHAN        the true species is NOT in the reference (n = 5,043) -> a novel
#                 species. Of these, 4,130 still have their GENUS in the reference,
#                 so a correct genus call is achievable and the genus comparison is
#                 fair; the other 913 are novel at genus too (no correct genus call
#                 exists), so they are kept in the TABLE but excluded from the FIGURE.
#
# Per sequence the BLCA vs NB genus call is scored against the truth as one of:
#   Both correct / BLCA correct only / NB correct only /
#   Agree, both wrong (same wrong genus) / Disagree, both wrong (different wrong
#   genera) / Unassigned (>=1 classifier abstained, neither correct).
# AGREEMENT = both classifiers correct. DISAGREEMENT = 1 - agreement.
#
# Outputs, per population:
#   figures/fig_agreement_taxa_16s_<pop>.{png,pdf}  top TOPN genera by disagreement, n>=MINN
#   tables/agreement_genus_16s_<pop>.csv            EVERY genus (no floor, n>=1)
# A minimum-n floor is required for the figure, not cosmetic: ranking all genera by
# disagreement with no floor returns 100 singleton genera tied at 100% and drops
# Escherichia off the list entirely. The CSVs carry the unfiltered data.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
CATS <- c("Both correct","BLCA correct only","NB correct only",
          "Agree, both wrong","Disagree, both wrong","Unassigned")
CAT_PAL <- c("Both correct"="#06d6a0","BLCA correct only"="#118ab2",
             "NB correct only"="#ef476f","Agree, both wrong"="#ffd166",
             "Disagree, both wrong"="#073b4c","Unassigned"="#9aa0a6")
MINN <- 5       # figure only: a genus needs this many test sequences to be ranked
TOPN <- 100     # figure only: show the TOPN most-disagreeing genera
LBL  <- 8       # genus label size in pt (fixed; page height is derived from it)

## ---- per-sequence genus outcome across the WHOLE 16S holdout ----
hd <- file.path(BR,"benchmark_holdout_16s_seed42"); rd <- file.path(hd,"results")
truth <- read_lineage_tsv(file.path(hd,"test_10.taxonomy"))
train <- read_lineage_tsv(file.path(hd,"reference_90.taxonomy"))
bl <- read_blca(rd)$tab; nb <- read_nb(rd)
ids <- Reduce(intersect, list(truth$id, bl$id, nb$id))
o <- function(x) x[match(ids,x$id),]; T<-o(truth); B<-o(bl); N<-o(nb)

tr_sp <- unique(na.omit(train$Species)); tr_ge <- unique(na.omit(train$Genus))
sp_in_ref <- !is.na(T$Species) & T$Species %in% tr_sp    # species seen in training?
ge_in_ref <- !is.na(T$Genus)   & T$Genus   %in% tr_ge    # genus seen in training?

tv<-T$Genus; bv<-B$Genus; nv<-N$Genus
bc<-!is.na(bv)&bv==tv; nc<-!is.na(nv)&nv==tv             # each classifier correct?
same<-!is.na(bv)&!is.na(nv)&bv==nv                       # both assigned the SAME genus
out<-rep("Unassigned",length(tv))                        # residual default: a call was missing
out[!bc&!nc & !is.na(bv)&!is.na(nv)&bv!=nv]<-"Disagree, both wrong"
out[same & !bc & !nc]<-"Agree, both wrong"
out[bc&!nc]<-"BLCA correct only"; out[nc&!bc]<-"NB correct only"; out[bc&nc]<-"Both correct"

d <- tibble(taxon=tv, outcome=out, agree=bc&nc,
            population = ifelse(sp_in_ref, "in_reference", "orphan"),
            genus_in_ref = ge_in_ref) %>% filter(!is.na(taxon))

## ---- per-genus summary (no floor: this is what the CSV reports) ----
genus_stats <- function(dd) {
  dd %>% group_by(taxon) %>%
    summarise(n = n(),
              genus_in_reference = first(genus_in_ref),
              both_correct       = sum(outcome=="Both correct"),
              blca_correct_only  = sum(outcome=="BLCA correct only"),
              nb_correct_only    = sum(outcome=="NB correct only"),
              agree_both_wrong   = sum(outcome=="Agree, both wrong"),
              disagree_both_wrong= sum(outcome=="Disagree, both wrong"),
              unassigned         = sum(outcome=="Unassigned"),
              .groups="drop") %>%
    mutate(agreement_pct    = round(100*both_correct/n, 2),   # both classifiers correct
           disagreement_pct = round(100 - 100*both_correct/n, 2)) %>%
    arrange(desc(disagreement_pct), desc(n))
}

## ---- one figure for one population ----
build_fig <- function(dd, pop, title, note) {
  per <- genus_stats(dd) %>%
    mutate(agree_lbl = ifelse(both_correct==n, 100L,
                              pmin(99L, as.integer(round(agreement_pct)))))  # 100 only if perfect
  scored <- per %>% filter(n>=MINN)
  # only genera that actually show disagreement; an all-green bar says nothing in a
  # "most disagreeing" chart (it matters for orphans, where most genera are perfect)
  sel   <- scored %>% filter(disagreement_pct > 0) %>%
           arrange(desc(disagreement_pct), desc(n)) %>% slice_head(n=TOPN)
  n_zero <- sum(scored$disagreement_pct == 0)          # perfect genera, omitted
  n_dis  <- nrow(scored) - n_zero                      # genera with any disagreement
  ny  <- nrow(sel)                                     # of those, how many are shown

  plotdat <- dd %>% filter(taxon %in% sel$taxon) %>%
    count(taxon, outcome, name="k") %>%
    group_by(taxon) %>% mutate(pct=100*k/sum(k)) %>% ungroup()
  plotdat$outcome <- factor(plotdat$outcome, levels=rev(CATS))       # Both correct at left
  ord <- sel %>% arrange(disagreement_pct) %>% pull(taxon)           # highest disagreement at top
  plotdat$taxon <- factor(plotdat$taxon, levels=ord)
  nlab <- sel %>% mutate(taxon=factor(taxon, levels=ord))

  PAGE_W <- 8.27                       # A4 width
  CHROME <- 1.55                       # title + 2-line subtitle + legend + x-axis + margins
  nlbl   <- LBL * 0.82
  row_pt <- LBL / 0.95                 # pitch that keeps an LBL-pt label clear of its neighbours
  PAGE_H <- CHROME + ny * row_pt / 72
  ytop   <- ny + 0.5

  p <- ggplot(plotdat, aes(pct, taxon, fill=outcome)) +
    geom_col(width=0.9, colour="white", linewidth=0.1) +
    geom_text(data=nlab, aes(x=103, y=taxon, label=format(n, big.mark=",")),
              inherit.aes=FALSE, hjust=0, size=nlbl/.pt, colour="grey35") +
    geom_text(data=nlab, aes(x=118, y=taxon, label=paste0(agree_lbl, "%")),
              inherit.aes=FALSE, hjust=0, size=nlbl/.pt, colour="grey35") +
    annotate("text", x=103, y=ytop, label="n", hjust=0, vjust=0, size=nlbl/.pt,
             fontface="bold", colour="grey30") +
    annotate("text", x=118, y=ytop, label="agree", hjust=0, vjust=0, size=nlbl/.pt,
             fontface="bold", colour="grey30") +
    coord_cartesian(clip="off") +
    scale_fill_manual(values=CAT_PAL, breaks=CATS, name=NULL, drop=FALSE) +
    scale_x_continuous("% of the genus's test sequences", limits=c(0,122),
                       breaks=seq(0,100,25), expand=expansion(mult=c(0,0))) +
    scale_y_discrete(expand=expansion(add=0.6)) +
    labs(y=NULL, title=title,
         subtitle=sprintf("%s\nAgreement = both classifiers correct. Bars: n, %% agreement. Showing %d of %d disagreeing genera (%d with %d+ seqs; %d perfect).",
                          note, ny, n_dis, nrow(scored), MINN, n_zero)) +
    guides(fill=guide_legend(nrow=1)) +
    theme_ms(12) +
    theme(legend.position="top", legend.text=element_text(size=8),
          legend.margin=margin(1,1,1,1), legend.box.spacing=unit(2,"pt"),
          legend.key.size=unit(9,"pt"),
          panel.grid=element_blank(), axis.ticks.y=element_blank(),
          plot.title=element_text(size=12.5),
          plot.subtitle=element_text(size=8, colour="grey30"),
          plot.margin=margin(4, 8, 4, 4),
          axis.text.y=element_text(face="italic", size=LBL, margin=margin(r=1)))

  f <- file.path(OUT,"figures",sprintf("fig_agreement_taxa_16s_%s", pop))
  ggsave(paste0(f,".png"), p, width=PAGE_W, height=PAGE_H, dpi=300, bg="white")
  ggsave(paste0(f,".pdf"), p, width=PAGE_W, height=PAGE_H, bg="white")
  cat(sprintf("[%s] %d seqs | %d genera n>=%d (%d with zero disagreement) | showing %d (%.1f%%-%.1f%% disagreement) | page %.2f x %.2f in\n",
              pop, nrow(dd), nrow(scored), MINN, n_zero, ny,
              min(sel$disagreement_pct), max(sel$disagreement_pct), PAGE_W, PAGE_H))
  invisible(sel)
}

dir.create(file.path(OUT,"tables"), showWarnings=FALSE, recursive=TRUE)

## ---- population 1: species present in the reference ----
d_ref <- d %>% filter(population=="in_reference")
write_csv(genus_stats(d_ref), file.path(OUT,"tables","agreement_genus_16s_in_reference.csv"))
build_fig(d_ref, "in_reference",
          "16S, species in reference: the genera BLCA and NB disagree on most",
          sprintf("Test sequences whose true species IS in the 90%% reference (n = %s), highest disagreement at top.",
                  format(nrow(d_ref), big.mark=",")))

## ---- population 2: orphans (novel species) ----
d_orp     <- d %>% filter(population=="orphan")                    # all 5,043 -> table
d_orp_fig <- d_orp %>% filter(genus_in_ref)                        # 4,130 -> figure (genus scoreable)
write_csv(genus_stats(d_orp), file.path(OUT,"tables","agreement_genus_16s_orphan.csv"))
build_fig(d_orp_fig, "orphan",
          "16S orphans (novel species): the genera BLCA and NB disagree on most",
          sprintf("Test sequences whose true species is NOT in the reference but whose GENUS is (n = %s), highest disagreement at top.",
                  format(nrow(d_orp_fig), big.mark=",")))

cat(sprintf("\nOrphan table also carries the %d sequences whose genus is absent from the reference (genus_in_reference=FALSE).\n",
            sum(!d_orp$genus_in_ref)))
cat("Wrote figures/fig_agreement_taxa_16s_{in_reference,orphan}.{png,pdf}\n")
cat("Wrote tables/agreement_genus_16s_{in_reference,orphan}.csv (every genus, no n floor)\n")
