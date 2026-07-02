library(dada2)

# Read from command line argument the input FASTQ
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript learnError.R <input_fastq> <cpu>", call. = FALSE)
}
fnFs <- args[1]
if (!file.exists(fnFs)) {
  stop(sprintf("Input FASTQ not found: %s", fnFs), call. = FALSE)
}
cpu <- as.integer(args[2])
if (is.na(cpu) || cpu < 1) {
  stop("cpu (arg 2) must be a positive integer.", call. = FALSE)
}

# Cap CPU to 16. Higher doesn't really help
if (cpu > 16) {
  cpu <- 16
}
# Learn the error rates
errF <- learnErrors(fnFs, multithread=cpu,  errorEstimationFunction=dada2:::PacBioErrfun)
err_plot <- plotErrors(errF)
pdf("plot_error_model.pdf", width=12, height=8, useDingbats=FALSE)
print(err_plot)
dev.off()

# Save as RDS
saveRDS(errF, file="errorfun.rds")
