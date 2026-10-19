# =====================================================================
# HiFiTaxa benchmark: shared helpers, sourced by every script in benchmark/analysis/.
# Paths come from environment variables (see Paths below and README.md).
# =====================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr); library(readr)
  library(patchwork); library(scales)
})

# --- Paths -------------------------------------------------------------------
# HIFITAXA_DATA (required): directory holding the pipeline run directories listed in
#   README.md (results_atcc_16s, results_fungi_its, benchmark_holdout_16s_seed42,
#   benchmark_holdout_its_seed42). They are not part of the repository.
# HIFITAXA_OUT (default <this folder>/regenerated): where figures/ and tables/ are written.
# HIFITAXA_REPO: a HiFiTaxa checkout, used only for the fungal mock truth list in
#   benchmark/mock_exclude/. Found automatically when these scripts run from
#   benchmark/analysis/ of a checkout; otherwise it must be set.
#     HIFITAXA_DATA=/path/to/run_dirs HIFITAXA_OUT=/path/to/out Rscript 01_holdout_accuracy.R
BR   <- Sys.getenv("HIFITAXA_DATA", unset = "")
OUT  <- Sys.getenv("HIFITAXA_OUT",  unset = file.path(.SCRIPT_DIR, "regenerated"))
REPO <- Sys.getenv("HIFITAXA_REPO", unset = "")
if (REPO == "") {
  .cand <- normalizePath(file.path(.SCRIPT_DIR, "..", ".."), mustWork = FALSE)
  if (file.exists(file.path(.cand, "benchmark", "mock_exclude", "fungi_103_species.txt"))) REPO <- .cand
}
if (BR == "" || !dir.exists(BR)) stop(sprintf(
  "Benchmark data not found (HIFITAXA_DATA = '%s').
Set HIFITAXA_DATA to the directory holding the run directories listed in benchmark/analysis/README.md:
  HIFITAXA_DATA=/path/to/run_dirs Rscript <script>.R", BR))
if (REPO == "") stop(
  "HIFITAXA_REPO is not set and this script is not running from benchmark/analysis/ of a HiFiTaxa checkout.
Set HIFITAXA_REPO to a HiFiTaxa clone (https://github.com/alegarritano/HiFiTaxa); it provides the
fungal mock truth list benchmark/mock_exclude/fungi_103_species.txt.")

# every script ggsave()s into OUT/figures and write.csv()s into OUT/tables
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "tables"),  recursive = TRUE, showWarnings = FALSE)
RANKS <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
TRUTH_KEY <- c(Domain="superkingdom",Phylum="phylum",Class="class",Order="order",
               Family="family",Genus="genus",Species="species")
# False positive threshold, as a fraction of the reads classified at that rank (0.001 = 0.1%, the
# manuscript value). Supplementary Table S11 was produced with FP_THRESH=0 and FP_THRESH=0.0001
# (07_mock_retrieval.R names its outputs by the threshold when it is not the default).
FP_THRESH <- as.numeric(Sys.getenv("FP_THRESH", unset = "0.001"))
if (is.na(FP_THRESH) || FP_THRESH < 0) stop("FP_THRESH must be a non negative number, for example 0.001")

## ---------- normalisation / synonym folding (matches manuscript + repo) ------
strip_ani <- function(x) str_replace_all(x, "_[A-Z]+(?=(\\s|$))", "")  # Cereibacter_A -> Cereibacter
GENUS_ALIAS <- c(G047199095="Escherichia",ECMA0423="Escherichia",Pauljensenia="Schaalia",
                 Thaxterogaster="Cortinarius",Phlegmacium="Cortinarius",Calonarius="Cortinarius",
                 # REVISED 2026-10: decision 1. UNITE files the mock member Calocybe obscurissima under
                 # Rugosomyces obscurissimus, so the genus Rugosomyces is folded to Calocybe for this mock.
                 Rugosomyces="Calocybe")
# REVISED 2026-10: decision 3. On the 16S mock a species call is credited only when it is the
# NCBI name or the GTDB r232 name of that mock strain (evidence per strain in RESCORE_REPORT.md).
#   REMOVED "ECMA0423 sp047199055" -> E. coli : a separate GTDB r232 species; the mock strain
#           (K-12 MG1655, GCF_000005845.2) is not in it.
#   REMOVED "Escherichia boydii" -> E. coli   : a separate GTDB r232 species (46 genomes, type strain
#           ATCC 8700 = NCBI Shigella boydii); the mock strain is not in it.
#   REMOVED "Pauljensenia odontolytica" -> Schaalia odontolytica : GTDB r232 places the mock strain
#           ATCC 17982 in Pauljensenia dentiphila, so P. odontolytica is another species.
#   KEPT    "G047199095 sp047199095" (GTDB r232 name of MG1655), "Pauljensenia dentiphila" (GTDB r232
#           name of ATCC 17982), "Pseudomonas aeruginosa" (author decision; ATCC label of ATCC 9027).
#   ADDED   "Schaalia dentiphila" (current NCBI name of ATCC 17982) and "Bacillus cereus" (NCBI name
#           of ATCC 10987, which GTDB r232 calls Bacillus_A pacificus). Both follow from the same rule.
SPECIES_ALIAS <- c("G047199095 sp047199095"="Escherichia coli",
                   "Pseudomonas aeruginosa"="Pseudomonas paraeruginosa",
                   "Pauljensenia dentiphila"="Schaalia odontolytica",
                   "Schaalia dentiphila"="Schaalia odontolytica",
                   "Bacillus cereus"="Bacillus pacificus",
                   # REVISED 2026-10: decision 1. Two fungal mock members appear in UNITE under a
                   # spelling or synonym variant; each pair is folded to the mock member.
                   "Cortinarius ochraceobrunneum"="Cortinarius ochraceobrunneus",
                   "Phlegmacium ochraceobrunneum"="Cortinarius ochraceobrunneus",
                   "Rugosomyces obscurissimus"="Calocybe obscurissima")
                   # NOTE: Staphylococcus argenteus is NOT folded to aureus. They are distinct
                   # validly published species; the mock has no argenteus, so a BLCA argenteus
                   # call is scored honestly as a false positive, not merged away.
# REVISED 2026-10: decision 3. Species of the two placeholder genera that GTDB r232 split from
# Escherichia keep their GTDB name at species level (they are reported, not renamed to
# "Escherichia sp..."). Genus level folding (GENUS_ALIAS) is unchanged.
KEEP_GTDB_NAME_GENUS <- c("ECMA0423","G047199095")
clean <- function(x){x<-str_trim(x);x[x %in% c("","NA","Unclassified","unclassified","s__","g__","NA NA")]<-NA;x}
norm_genus <- function(g){g<-strip_ani(str_trim(g));ifelse(g %in% names(GENUS_ALIAS),GENUS_ALIAS[g],g)}
norm_species <- function(s){
  s<-str_trim(s); s<-ifelse(s %in% names(SPECIES_ALIAS),SPECIES_ALIAS[s],s); s<-strip_ani(s)
  s<-ifelse(s %in% names(SPECIES_ALIAS),SPECIES_ALIAS[s],s)
  # REVISED 2026-10: decision 3. no genus renaming for KEEP_GTDB_NAME_GENUS at species level
  t<-str_split_fixed(s," ",2); g<-ifelse(t[,1] %in% names(GENUS_ALIAS) & !(t[,1] %in% KEEP_GTDB_NAME_GENUS),GENUS_ALIAS[t[,1]],t[,1])
  o<-ifelse(t[,2]=="",g,paste(g,t[,2])); ifelse(o %in% names(SPECIES_ALIAS),SPECIES_ALIAS[o],o)}
# REVISED 2026-10: decisions 3 and 4. The placeholder word must now be a whole word. The old
# pattern "^(sp\\.?|...)" matched any epithet that only BEGINS with sp (spinosa, sphaerica) and the
# GTDB genome defined species (sp047199055), which are real species. Labels ending in " sp" and
# labels containing "Incertae sedis" (UNITE "GS01 phy Incertae sedis sp") are also placeholders.
is_placeholder_sp <- function(s){ep<-str_split_fixed(s," ",2)[,2]
  ep=="" | str_detect(ep,"^(sp|spp|unidentified|uncultured|indet|incertae)\\.?($|\\s)") |
    str_detect(s,regex("(\\ssp\\.?$|incertae[ _]sedis)",ignore_case=TRUE))}
# REVISED 2026-10: decision 4. UNITE genus labels of the form "Agaricales_gen_Incertae_sedis".
is_placeholder_genus <- function(g) str_detect(g,regex("incertae[ _]sedis",ignore_case=TRUE))

## ---------- truth lineage parser (test_10.taxonomy / reference_90.taxonomy) ---
# rows: id <tab> superkingdom:..;phylum:..;...;species:..;   (k__ never used here)
read_lineage_tsv <- function(path) {
  raw <- read_tsv(path, col_names=c("id","lin"), show_col_types=FALSE, progress=FALSE)
  out <- tibble(id=raw$id)
  for (r in RANKS) out[[r]] <- clean(str_match(raw$lin, paste0("(?:^|;)",TRUTH_KEY[[r]],":([^;]*)"))[,2])
  out
}

## ---------- classifier prediction parsers ----------
read_blca <- function(dir_or_tab, path_conf=NULL) {
  tabp <- if (is.null(path_conf)) file.path(dir_or_tab,"taxonomy_blca","blca_taxonomy_table.csv") else dir_or_tab
  cfp  <- if (is.null(path_conf)) file.path(dir_or_tab,"taxonomy_blca","blca_taxonomy_confidence.csv") else path_conf
  tab <- read_csv(tabp, show_col_types=FALSE, progress=FALSE); names(tab)[1]<-"id"
  for (r in RANKS) tab[[r]] <- clean(tab[[r]])
  conf <- read_csv(cfp, show_col_types=FALSE, progress=FALSE); names(conf)[1]<-"id"
  list(tab=tab, conf=conf)
}
# NB best_taxonomy.tsv: Feature ID / Taxon(d__|k__ ..;s__ ..) / Confidence
read_nb <- function(dir_or_path) {
  p <- if (dir.exists(dir_or_path)) file.path(dir_or_path,"taxonomy_nb","best_taxonomy.tsv") else dir_or_path
  nb <- read_tsv(p, show_col_types=FALSE, progress=FALSE); names(nb)[1:2]<-c("id","Taxon")
  prefix <- c(Domain="(?:d|k)__",Phylum="p__",Class="c__",Order="o__",Family="f__",Genus="g__",Species="s__")
  out <- tibble(id=nb$id)
  for (r in RANKS) out[[r]] <- clean(str_match(nb$Taxon, paste0(prefix[[r]],"([^;]*)"))[,2])
  out
}

## ---------- abundance readers ----------
# dada2 feature table: line 1 is a '# Constructed...' comment, line 2 header (#OTU ID + samples)
read_feature_table <- function(dir) {
  f <- file.path(dir,"dada2","feature_table.tsv")
  ft <- read_tsv(f, comment="# Constructed", show_col_types=FALSE, progress=FALSE)
  names(ft)[1] <- "id"
  samp <- setdiff(names(ft), "id")
  tibble(id=ft$id, count=rowSums(sapply(ft[samp], as.numeric), na.rm=TRUE))
}
# taxon -> mean-per-sample relative abundance from ASV (id,taxon) + counts
relab_from_counts <- function(ids, taxa, counts) {
  d <- tibble(id=ids, taxon=taxa) %>% inner_join(counts, by="id") %>%
       filter(!is.na(taxon)) %>% group_by(taxon) %>% summarise(count=sum(count),.groups="drop")
  d$relab <- d$count/sum(d$count); d %>% select(taxon,relab)
}
# Emu wide table (16S): emu_<rank>_table.tsv  (rank cols + samples), rows "rank:taxon"
read_emu <- function(dir, rank) {
  f <- file.path(dir,"taxonomy_emu",sprintf("emu_%s_table.tsv",rank))
  e <- read_tsv(f, show_col_types=FALSE, progress=FALSE)
  samp <- setdiff(names(e), c("species","genus","family","order","class","phylum","superkingdom"))
  tax <- clean(str_replace(e[[rank]], paste0("^",rank,":"), ""))
  # REVISED 2026-10: decision 2. Mean over ALL libraries: a library in which Emu does not report
  # the taxon (empty cell) counts as 0. The old rowMeans(..., na.rm=TRUE) skipped those libraries.
  # (read_emits below pools reads over the libraries, so it never skipped a library: unchanged.)
  em <- as.matrix(e[samp]); em[is.na(em)] <- 0
  m <- rowMeans(em); tibble(taxon=tax, relab=m/sum(m)) %>% filter(!is.na(taxon))
}
# EMITS full per-sample tables (ITS): taxon|acc|SH|type|k__..;s__.. <tab> relab <tab> reads
# rank='species' -> s__ ; rank='genus' -> g__ . Combine samples by summed reads.
read_emits <- function(dir, rank) {
  fs <- list.files(file.path(dir,"taxonomy_emits","per_sample"),
                   pattern=sprintf("_emits_%s\\.tsv$",rank), full.names=TRUE)
  pref <- if (rank=="species") "s__" else "g__"
  rows <- list()
  for (f in fs) {
    d <- suppressWarnings(read_tsv(f, col_names=c("taxstr","relab","reads"),
                                   show_col_types=FALSE, progress=FALSE))
    d <- d[!is.na(d$reads),]
    lin <- str_split_fixed(d$taxstr, "\\|", 5)[,5]
    tax <- str_replace_all(clean(str_match(lin, paste0(pref,"([^;]*)"))[,2]), "_", " ")
    rows[[f]] <- tibble(taxon=tax, reads=as.numeric(d$reads))
  }
  d <- bind_rows(rows) %>% filter(!is.na(taxon)) %>% group_by(taxon) %>%
       summarise(reads=sum(reads),.groups="drop")
  d$relab <- d$reads/sum(d$reads); d %>% select(taxon,relab)
}

## ---------- truth taxa sets (mocks) ----------
atcc_species <- unique(norm_species(c(
  "Acinetobacter baumannii","Bacillus pacificus","Phocaeicola vulgatus","Bifidobacterium adolescentis",
  "Clostridium beijerinckii","Cutibacterium acnes","Deinococcus radiodurans","Enterococcus faecalis",
  "Escherichia coli","Helicobacter pylori","Lactobacillus gasseri","Neisseria meningitidis",
  "Porphyromonas gingivalis","Pseudomonas paraeruginosa","Cereibacter sphaeroides","Schaalia odontolytica",
  "Staphylococcus aureus","Staphylococcus epidermidis","Streptococcus agalactiae","Streptococcus mutans")))
atcc_genus <- unique(norm_genus(str_split_fixed(atcc_species," ",2)[,1]))
# ATCC staggered expected relative abundances (compare_blca_emu_taxaplot.py)
ATCC_TRUTH <- c("Acinetobacter baumannii"=0.0018,"Bacillus pacificus"=0.0180,"Phocaeicola vulgatus"=0.0002,
  "Bifidobacterium adolescentis"=0.0002,"Clostridium beijerinckii"=0.0180,"Cutibacterium acnes"=0.0018,
  "Deinococcus radiodurans"=0.0002,"Enterococcus faecalis"=0.0002,"Escherichia coli"=0.1800,
  "Helicobacter pylori"=0.0018,"Lactobacillus gasseri"=0.0018,"Neisseria meningitidis"=0.0018,
  "Porphyromonas gingivalis"=0.1800,"Pseudomonas paraeruginosa"=0.0180,"Cereibacter sphaeroides"=0.1800,
  "Schaalia odontolytica"=0.0002,"Staphylococcus aureus"=0.0180,"Staphylococcus epidermidis"=0.1800,
  "Streptococcus agalactiae"=0.0180,"Streptococcus mutans"=0.1800)
.fung_truth <- file.path(REPO,"benchmark/mock_exclude/fungi_103_species.txt")
if (!file.exists(.fung_truth)) stop(sprintf(
  "Mock truth list not found at '%s'.
Set HIFITAXA_REPO to a HiFiTaxa checkout (git clone https://github.com/alegarritano/HiFiTaxa).
The list is read whenever 00_common.R is sourced, so it is required even by the holdout scripts.", .fung_truth))
.fr <- readLines(.fung_truth); .fr<-str_trim(.fr)
.fr <- .fr[.fr!="" & !str_starts(.fr,"#")]
fung_species <- unique(norm_species(.fr))
fung_genus   <- unique(norm_genus(str_split_fixed(fung_species," ",2)[,1]))

message(sprintf("[common] ATCC truth: %d species / %d genera | Fungal truth: %d species / %d genera",
                length(atcc_species), length(atcc_genus), length(fung_species), length(fung_genus)))

## ---------- manuscript figure grammar (colors, palette, theme) ----------
CLF_COL <- c("Expected"="#888888","Truth"="#888888",
             "BLCA"="#0072B2","BLCA@80"="#56B4E9","Emu"="#D55E00",
             "NB"="#CC79A7","EMITS"="#009E73")
## Unassigned lightened from #888888 to #BBBBBB: at #888888 the Incorrect/Unassigned
## pair collapsed under protanopia (dE 5.9). #BBBBBB separates them on lightness
## (dE 15.6 protan), which survives CVD. Verified with the dataviz palette validator.
CAT_COL <- c("Correct"="#44AA99","Incorrect"="#CC6677","Unassigned"="#BBBBBB",
             "Recovered"="#44AA99","Missed"="#bdbdbd","False positive"="#CC6677",
             "Detected"="#44AA99")
# Okabe-Ito extended qualitative palette for taxa (top-N), + Other/Unassigned
TAXA_PAL <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
              "#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#882255",
              "#661100","#6699CC","#AA4499","#44AA55","#DDDDDD","#000000")
OTHER_COL <- "#cccccc"; UNASSIGNED_COL <- "#888888"
theme_ms <- function(base=12) theme_classic(base_size=base) +
  theme(plot.title=element_text(face="bold",hjust=0),
        strip.background=element_rect(fill="grey95",colour="black"),
        strip.text=element_text(face="bold"))

## ---------- BLCA@80: BLCA table with genus/species masked below conf 80 -----
read_blca_at80 <- function(dir) {
  b <- read_blca(dir); tab <- b$tab; cf <- b$conf
  cf <- cf[match(tab$id, cf$id),]
  for (r in RANKS) { cc<-suppressWarnings(as.numeric(cf[[paste0(r,"_conf")]])); cc[is.na(cc)]<-0
                     tab[[r]][cc<80] <- NA }
  tab
}

## ---------- composition: read-weighted relab per taxon incl 'Unassigned' ----
# reads with NA taxon are folded into an 'Unassigned' category; result sums to 1.
comp_reads <- function(ids, taxa, counts) {
  d <- tibble(id=ids, taxon=taxa) %>% inner_join(counts, by="id")
  d$taxon[is.na(d$taxon)] <- "Unassigned"
  d %>% group_by(taxon) %>% summarise(count=sum(count),.groups="drop") %>%
    mutate(relab=count/sum(count)) %>% select(taxon,relab)
}
