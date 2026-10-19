#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Figure S5: classifier run time.
# Parsed from the Nextflow `pipeline_trace.txt` written by each benchmark run.
#
# Two measures, both in the same unit (time), so they share one axis:
#   Total task time  = sum of every task's `realtime` for that classifier
#                      (the COMPUTE cost; concurrency independent)
#   Wall clock       = max(submit + duration) - min(submit) across its tasks
#                      (the ELAPSED time a user waits, given the job's cores)
# BLCA parallelises across query chunks, so its wall clock is far below its
# total task time; NB is a single task; Emu/EMITS run one task per sample.
#
# Only the CLASSIFICATION processes are attributed to a classifier. The shared
# front end (QC, cutadapt/itsxrust, DADA2 denoise, import) is shown separately.
#
# Design: magnitude spanning ~6 orders of magnitude => log10 axis and a DOT plot.
# Bars on a log axis would misrepresent magnitude by length, so points are used,
# each directly labelled with a human readable duration.
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
## ---- parse Nextflow duration strings: "1d 2h 39m 39s", "1m 14s", "155ms", "5s" ----
parse_dur <- function(x) vapply(x, function(s) {
  if (is.na(s) || s %in% c("", "-")) return(NA_real_)
  tot <- 0
  m <- regmatches(s, regexpr("[0-9.]+ms", s))          # milliseconds first
  if (length(m) == 1) { tot <- tot + as.numeric(sub("ms","",m))/1000; s <- sub("[0-9.]+ms","",s) }
  for (nm in c("d","h","m","s")) {
    mult <- switch(nm, d=86400, h=3600, m=60, s=1)
    mm <- regmatches(s, regexpr(paste0("[0-9.]+", nm), s))
    if (length(mm) == 1) tot <- tot + as.numeric(sub(nm,"",mm)) * mult
  }
  tot
}, numeric(1), USE.NAMES=FALSE)

fmt_dur <- function(s) ifelse(s < 60,    sprintf("%.1f s", s),
                       ifelse(s < 3600,  sprintf("%.0f m", s/60),
                       ifelse(s < 86400, sprintf("%.1f h", s/3600),
                                         sprintf("%.1f d", s/86400))))

## ---- map a Nextflow process name to a classifier (or the shared front end) ----
classify_proc <- function(p) {
  p <- sub("^.*:", "", p)          # drop "taxonomy_only:" prefix
  p <- sub(" \\(.*$", "", p)       # drop " (3)" task index
  dplyr::case_when(
    p %in% c("blca_classify","merge_blca")                                  ~ "BLCA",
    grepl("^nb_classify", p)                                                ~ "NB",
    p %in% c("emu_classify","emu_collate")                                  ~ "Emu",
    p == "emits_classify"                                                   ~ "EMITS",
    TRUE                                                                    ~ "Pipeline front end")
}

read_trace <- function(path, dataset) {
  tr <- read_tsv(path, show_col_types=FALSE, progress=FALSE)
  tr %>% transmute(dataset,
                   proc     = name,
                   step     = classify_proc(name),
                   realtime = parse_dur(realtime),
                   submit   = as.POSIXct(submit, format="%Y-%m-%d %H:%M:%OS", tz="UTC"),
                   dur      = parse_dur(duration)) %>%
    filter(!is.na(realtime))
}

BRp <- BR
runs <- tribble(
  ~dataset,                                   ~path,
  "ATCC 16S mock\n(192 libraries, 58 ASVs)",  file.path(BRp,"results_atcc_16s/pipeline_trace.txt"),
  "Fungal ITS mock\n(2 libraries, 138 ASVs)", file.path(BRp,"results_fungi_its/pipeline_trace.txt"),
  "16S holdout\n(95,731 query sequences)",    file.path(BRp,"benchmark_holdout_16s_seed42/results/pipeline_trace.txt"),
  "ITS holdout\n(10,213 query sequences)",    file.path(BRp,"benchmark_holdout_its_seed42/results/pipeline_trace.txt"))

tr <- bind_rows(lapply(seq_len(nrow(runs)), function(i) read_trace(runs$path[i], runs$dataset[i])))

summ <- tr %>% group_by(dataset, step) %>%
  summarise(n_tasks        = n(),
            total_task_sec = sum(realtime, na.rm=TRUE),
            wall_sec       = as.numeric(difftime(max(submit + dur, na.rm=TRUE),
                                                 min(submit, na.rm=TRUE), units="secs")),
            .groups="drop")

dir.create(file.path(OUT,"tables"), showWarnings=FALSE, recursive=TRUE)
# flatten the facet-strip newlines so the supplementary table is a tidy CSV
write_csv(summ %>% mutate(dataset=gsub("\n", " ", dataset),
                          total_task_sec=round(total_task_sec,1), wall_sec=round(wall_sec,1),
                          total_task_time=fmt_dur(total_task_sec), wall_clock=fmt_dur(wall_sec)),
          file.path(OUT,"tables","classifier_runtimes.csv"))

## ---- figure: classifiers only, dot plot on a log10 time axis ----
d <- summ %>% filter(step != "Pipeline front end") %>%
  pivot_longer(c(total_task_sec, wall_sec), names_to="measure", values_to="sec") %>%
  mutate(measure = factor(measure, levels=c("total_task_sec","wall_sec"),
                          labels=c("Total task time (compute)","Wall clock (elapsed)")),
         step    = factor(step, levels=c("BLCA","NB","Emu","EMITS")),
         dataset = factor(dataset, levels=runs$dataset))

brks <- c(1, 10, 60, 600, 3600, 21600, 86400, 604800, 2592000)
labs <- c("1 s","10 s","1 min","10 min","1 h","6 h","1 d","7 d","30 d")
# upper limit must clear the largest value (BLCA 16S holdout total = 46.2 d ~ 4.0e6 s)
# plus room for its direct label, or the point is silently dropped.
p <- ggplot(d, aes(sec, step, colour=step)) +
  geom_point(size=3) +
  geom_text(aes(label=fmt_dur(sec)), hjust=-0.25, size=8/.pt, colour="grey20") +
  facet_grid(dataset ~ measure, scales="free_y", space="free_y") +
  scale_x_log10(breaks=brks, labels=labs,
                limits=c(0.5, 3e7), expand=expansion(mult=c(0.02,0.02))) +
  scale_colour_manual(values=CLF_COL, guide="none") +
  labs(x="Run time (log scale)", y=NULL,
       title="Classifier run time",
       subtitle="NCI Gadi, 48 cores per job. BLCA parallelises across query chunks, so its wall clock sits well below its total task time.") +
  theme_ms(11) +
  theme(panel.grid.major.y=element_line(colour="grey92", linetype="dotted"),
        panel.grid.major.x=element_line(colour="grey94"),
        panel.grid.minor=element_blank(),
        strip.text.y=element_text(angle=0, size=8),
        plot.subtitle=element_text(size=9, colour="grey30"))

ggsave(file.path(OUT,"figures","fig_classifier_runtimes.png"), p, width=11, height=8, dpi=200)
ggsave(file.path(OUT,"figures","fig_classifier_runtimes.pdf"), p, width=11, height=8)

cat("\n== Classifier run times ==\n")
summ %>% mutate(total_task_time=fmt_dur(total_task_sec), wall_clock=fmt_dur(wall_sec)) %>%
  select(dataset, step, n_tasks, total_task_time, wall_clock) %>%
  mutate(dataset=gsub("\n"," ",dataset)) %>% as.data.frame() %>% print(row.names=FALSE)
cat("\nWrote figures/fig_classifier_runtimes.{png,pdf} + tables/classifier_runtimes.csv\n")
