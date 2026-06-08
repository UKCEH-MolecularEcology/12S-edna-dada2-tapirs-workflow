# =========================
# 02_dada2_asv.R
# DADA2 pipeline and ASV export
# =========================

.args       <- commandArgs(trailingOnly = FALSE)
.script_dir <- dirname(normalizePath(sub("--file=", "", .args[grep("--file=", .args)])))
source(file.path(.script_dir, "00_config.R"))

library(dada2)
library(Biostrings)
library(tibble)
library(readr)

ts <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")

# -------------------------
# Locate cutadapt outputs
# -------------------------
fnFs <- sort(list.files(
  cutadapt_dir,
  pattern     = "_R1_001.*\\.(fastq|fq)(\\.gz)?$",
  full.names  = TRUE
))
fnRs <- sort(list.files(
  cutadapt_dir,
  pattern     = "_R2_001.*\\.(fastq|fq)(\\.gz)?$",
  full.names  = TRUE
))

if (length(fnFs) == 0 || length(fnRs) == 0) {
  stop("No cutadapt output files found in: ", cutadapt_dir)
}

sample.names <- sub("_R1_001.*$", "", basename(fnFs))

filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))

# -------------------------
# Quality filter and trim
# -------------------------
ts("STEP 1/6: Filtering and trimming reads")
out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen  = c(75, 75),
  maxN      = 0,
  maxEE     = c(as.numeric(Sys.getenv("SM_MAXEE_F", "2")), as.numeric(Sys.getenv("SM_MAXEE_R", "5"))),
  truncQ    = 2,
  rm.phix   = TRUE,
  compress  = TRUE,
  multithread = threads
)

write.csv(out, file.path(results_dir, "filtering_summary.csv"))

# Remove samples that produced no reads after filtering
passed  <- file.exists(filtFs) & file.exists(filtRs) &
           (file.info(filtFs)$size > 0) & (file.info(filtRs)$size > 0)

if (!any(passed)) stop("No reads survived quality filtering across any sample.")
if (sum(!passed) > 0) {
  cat("WARNING:", sum(!passed), "sample(s) had no reads after filtering and are excluded.\n")
}

fnFs_passed  <- fnFs[passed]
sample.names <- sample.names[passed]
filtFs       <- filtFs[passed]
filtRs       <- filtRs[passed]

# -------------------------
# Learn errors
# -------------------------
ts("STEP 2/6: Learning error rates (", length(filtFs), "samples,", threads, "threads)")
ts("  Learning forward error rates ...")
errF <- learnErrors(filtFs, nbases = 1e8, multithread = threads)
ts("  Learning reverse error rates ...")
errR <- learnErrors(filtRs, nbases = 1e8, multithread = threads)
ts("  Error rates learned")

# -------------------------
# Per-sample: dereplicate, denoise, merge
# Processing one sample at a time minimises peak memory use.
# -------------------------
names(filtFs) <- sample.names
names(filtRs) <- sample.names

ts("STEP 3/6: Per-sample dereplicate / denoise / merge (", length(sample.names), "samples)")

mergers          <- vector("list", length(sample.names))
names(mergers)   <- sample.names
denoisedF_counts <- setNames(integer(length(sample.names)), sample.names)
denoisedR_counts <- setNames(integer(length(sample.names)), sample.names)

for (i in seq_along(sample.names)) {
  sam <- sample.names[i]
  ts(sprintf("  [%d/%d] %s", i, length(sample.names), sam))
  derepF <- derepFastq(filtFs[[sam]], verbose = FALSE)
  ddF    <- dada(derepF, err = errF, multithread = threads, verbose = FALSE)
  derepR <- derepFastq(filtRs[[sam]], verbose = FALSE)
  ddR    <- dada(derepR, err = errR, multithread = threads, verbose = FALSE)
  mergers[[sam]]          <- mergePairs(ddF, derepF, ddR, derepR, verbose = TRUE)
  denoisedF_counts[[sam]] <- sum(getUniques(ddF))
  denoisedR_counts[[sam]] <- sum(getUniques(ddR))
}

# -------------------------
# Sequence table + chimera removal
# -------------------------
ts("STEP 4/6: Making sequence table")
seqtab <- makeSequenceTable(mergers)
write.csv(seqtab, file.path(results_dir, "seqtab_raw.csv"))
ts("  Sequence table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs")

ts("STEP 5/6: Removing chimeras")
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method      = "consensus",
  multithread = threads,
  verbose     = TRUE
)
write.csv(seqtab.nochim, file.path(results_dir, "seqtab_nochim.csv"))

# -------------------------
# Export ASVs
# -------------------------
ts("STEP 6/6: Exporting ASVs —", ncol(seqtab.nochim), "sequences after chimera removal")
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))

dna <- DNAStringSet(asv_seqs)
names(dna) <- asv_ids
writeXStringSet(dna, filepath = file.path(results_dir, "asvs.nochim.fasta"), format = "fasta")

asv_lookup <- tibble(ASV = asv_ids, sequence = asv_seqs, length = nchar(asv_seqs))
write_tsv(asv_lookup, file.path(results_dir, "asv_lookup.tsv"))

seqtab_asv             <- seqtab.nochim
colnames(seqtab_asv)   <- asv_ids

write.csv(
  data.frame(sample = rownames(seqtab_asv), seqtab_asv, check.names = FALSE),
  file.path(results_dir, "seqtab_asv.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(ASV = colnames(seqtab_asv), t(seqtab_asv), check.names = FALSE),
  file.path(results_dir, "seqtab_asv_transposed.csv"),
  row.names = FALSE
)
saveRDS(seqtab_asv,    file.path(results_dir, "seqtab_asv.rds"))
saveRDS(seqtab.nochim, file.path(results_dir, "seqtab_nochim.rds"))
saveRDS(errF,          file.path(results_dir, "errF.rds"))
saveRDS(errR,          file.path(results_dir, "errR.rds"))

# -------------------------
# Read tracking table
# -------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out[fnFs_passed, , drop = FALSE],
  denoisedF_counts,
  denoisedR_counts,
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file.path(results_dir, "dada2_read_tracking.csv"))

ts("DADA2 complete")
