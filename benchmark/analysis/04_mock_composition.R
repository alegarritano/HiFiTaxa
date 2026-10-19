#!/usr/bin/env Rscript
# =====================================================================
# ATCC mock — per-genus relative-abundance heatmap (manuscript Fig. S1c style).
# Each cell shows the OBSERVED read-weighted relative abundance (%), coloured by
# its deviation from the expected staggered composition (blue = under, red =
# over, white = match). Columns: Expected + BLCA, BLCA@80, Emu, NB. Rows: the 18
# mock genera (by expected abundance), then Unassigned, then any foreign genus a
# method reports above 0.1%. The stacked composition taxaplots are in
# 06_mock_taxaplots.R.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
dir <- file.path(BR,"results_atcc_16s"); cnt <- read_feature_table(dir)
bl <- read_blca(dir)$tab; b8 <- read_blca_at80(dir); nb <- read_nb(dir)

exp_g <- tibble(species=names(ATCC_TRUTH), relab=as.numeric(ATCC_TRUTH)) %>%
  mutate(taxon=norm_genus(str_split_fixed(norm_species(species)," ",2)[,1])) %>%
  group_by(taxon) %>% summarise(Expected=sum(relab),.groups="drop")

# observed genus relative abundance per method (comp_reads folds NA -> Unassigned)
g_blca <- comp_reads(bl$id, norm_genus(bl$Genus), cnt) %>% rename(BLCA=relab)
g_b80  <- comp_reads(b8$id, norm_genus(b8$Genus), cnt) %>% rename(`BLCA@80`=relab)
g_nb   <- comp_reads(nb$id, norm_genus(nb$Genus), cnt) %>% rename(NB=relab)
g_emu  <- read_emu(dir,"genus") %>% mutate(taxon=norm_genus(taxon)) %>%
          group_by(taxon) %>% summarise(Emu=sum(relab),.groups="drop")

comp <- Reduce(function(a,b) full_join(a,b,by="taxon"),
               list(exp_g, g_blca, g_b80, g_emu, g_nb)) %>%
        mutate(across(-taxon, ~replace_na(.,0)))

# row set + order: mock genera (by expected) -> Unassigned -> foreign (>0.1% any method)
is_mock <- comp$taxon %in% atcc_genus
maxobs  <- pmax(comp$BLCA, comp$`BLCA@80`, comp$Emu, comp$NB)
comp <- comp %>% filter(is_mock | taxon=="Unassigned" | maxobs>0.001)
mock_ord <- comp %>% filter(taxon %in% atcc_genus) %>% arrange(desc(Expected)) %>% pull(taxon)
foreign_ord <- comp %>% filter(!taxon %in% atcc_genus, taxon!="Unassigned") %>%
               mutate(m=pmax(BLCA,`BLCA@80`,Emu,NB)) %>% arrange(desc(m)) %>% pull(taxon)
row_ord <- c(mock_ord, "Unassigned", foreign_ord)
row_ord <- row_ord[row_ord %in% comp$taxon]
write_csv(comp %>% mutate(across(c(Expected,BLCA,`BLCA@80`,Emu,NB),~round(100*.,3))),
          file.path(OUT,"tables","atcc_composition_genus.csv"))

hm <- comp %>% pivot_longer(c(Expected,BLCA,`BLCA@80`,Emu,NB),
                            names_to="method", values_to="obs") %>%
  left_join(comp %>% select(taxon, Exp=Expected), by="taxon") %>%
  mutate(dev = ifelse(method=="Expected", 0, obs - Exp))
hm$method <- factor(hm$method, levels=c("Expected","BLCA","BLCA@80","Emu","NB"))
hm$taxon  <- factor(hm$taxon,  levels=rev(row_ord))

p <- ggplot(hm, aes(method, taxon, fill=dev)) +
  geom_tile(colour="grey88") +
  geom_text(aes(label=sprintf("%.1f", 100*obs)), size=2.7) +
  scale_fill_gradient2(low="#2166ac", mid="white", high="#b2182b", midpoint=0,
                       limits=c(-0.2,0.2), oob=scales::squish,
                       labels=percent_format(accuracy=1),
                       name="observed\nminus expected\n(white = match)") +
  scale_x_discrete(position="top") +
  labs(x=NULL, y=NULL,
       title="ATCC MSA-1003 mock: per-genus relative abundance (%)") +
  theme_minimal(base_size=11) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold",size=12),
        axis.text.y=element_text(face="italic"),
        axis.text.x=element_text(face="bold"))
ggsave(file.path(OUT,"figures","fig_atcc_deviation_heatmap.png"), p, width=7.5, height=7, dpi=200)
ggsave(file.path(OUT,"figures","fig_atcc_deviation_heatmap.pdf"), p, width=7.5, height=7)

cat("\n== ATCC genus relative abundance (%) ==\n")
comp %>% arrange(match(taxon,row_ord)) %>%
  mutate(across(c(Expected,BLCA,`BLCA@80`,Emu,NB),~round(100*.,2))) %>%
  as.data.frame() %>% print(row.names=FALSE)
cat("\nWrote figures/fig_atcc_deviation_heatmap.{png,pdf}\n")
