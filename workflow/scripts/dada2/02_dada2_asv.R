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
# Checkpoint helpers
# Saved to a hidden dir not declared as Snakemake outputs,
# so they survive job failures and allow resuming mid-script.
# -------------------------
chk_dir  <- file.path(results_dir, ".checkpoints")
dir.create(chk_dir, showWarnings = FALSE, recursive = TRUE)
chk_path <- function(name) file.path(chk_dir, paste0(name, ".rds"))
has_chk  <- function(name) file.exists(chk_path(name))
save_chk <- function(obj, name) { saveRDS(obj, chk_path(name)); ts("  Checkpoint saved:", name) }
load_chk <- function(name) { ts("  Loading checkpoint:", name); readRDS(chk_path(name)) }

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
# STEP 1/6: Quality filter and trim
# -------------------------
if (has_chk("filter")) {
  ts("STEP 1/6: Filtering — loading checkpoint")
  chk          <- load_chk("filter")
  out          <- chk$out
  fnFs_passed  <- chk$fnFs_passed
  sample.names <- chk$sample.names
  filtFs       <- chk$filtFs
  filtRs       <- chk$filtRs
} else {
  ts("STEP 1/6: Filtering and trimming reads (", length(fnFs), "samples,", threads, "threads)")
  out <- filterAndTrim(
    fnFs, filtFs,
    fnRs, filtRs,
    truncLen    = c(75, 75),
    maxN        = 0,
    maxEE       = c(as.numeric(Sys.getenv("SM_MAXEE_F", "2")), as.numeric(Sys.getenv("SM_MAXEE_R", "2"))),
    truncQ      = 2,
    rm.phix     = TRUE,
    compress    = TRUE,
    multithread = threads
  )
  write.csv(out, file.path(results_dir, "filtering_summary.csv"))

  passed <- file.exists(filtFs) & file.exists(filtRs) &
            (file.info(filtFs)$size > 0) & (file.info(filtRs)$size > 0)

  if (!any(passed)) stop("No reads survived quality filtering across any sample.")
  if (sum(!passed) > 0)
    cat("WARNING:", sum(!passed), "sample(s) had no reads after filtering and are excluded.\n")

  fnFs_passed  <- fnFs[passed]
  sample.names <- sample.names[passed]
  filtFs       <- filtFs[passed]
  filtRs       <- filtRs[passed]

  save_chk(list(out=out, fnFs_passed=fnFs_passed,
                sample.names=sample.names, filtFs=filtFs, filtRs=filtRs), "filter")
}

# -------------------------
# STEP 2/6: Learn errors
# -------------------------
if (has_chk("errors")) {
  ts("STEP 2/6: Error rates — loading checkpoint")
  errs <- load_chk("errors")
  errF <- errs$errF
  errR <- errs$errR
} else {
  ts("STEP 2/6: Learning error rates (", length(filtFs), "samples,", threads, "threads)")
  ts("  Learning forward error rates ...")
  errF <- learnErrors(filtFs, nbases = 1e8, multithread = threads)
  ts("  Learning reverse error rates ...")
  errR <- learnErrors(filtRs, nbases = 1e8, multithread = threads)
  ts("  Error rates learned")
  save_chk(list(errF=errF, errR=errR), "errors")
}

# -------------------------
# STEP 3/6: Per-sample dereplicate / denoise / merge
# -------------------------
if (has_chk("denoised")) {
  ts("STEP 3/6: Denoising — loading checkpoint")
  dn               <- load_chk("denoised")
  mergers          <- dn$mergers
  denoisedF_counts <- dn$denoisedF_counts
  denoisedR_counts <- dn$denoisedR_counts
} else {
  ts("STEP 3/6: Per-sample dereplicate / denoise / merge (", length(sample.names), "samples)")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

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
  save_chk(list(mergers=mergers,
                denoisedF_counts=denoisedF_counts,
                denoisedR_counts=denoisedR_counts), "denoised")
}

# -------------------------
# STEP 4/6: Sequence table
# -------------------------
if (has_chk("seqtab")) {
  ts("STEP 4/6: Sequence table — loading checkpoint")
  seqtab <- load_chk("seqtab")
} else {
  ts("STEP 4/6: Making sequence table")
  seqtab <- makeSequenceTable(mergers)
  write.csv(seqtab, file.path(results_dir, "seqtab_raw.csv"))
  ts("  Sequence table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs")
  save_chk(seqtab, "seqtab")
}

# -------------------------
# STEP 5/6: Chimera removal
# -------------------------
if (has_chk("nochim")) {
  ts("STEP 5/6: Chimera removal — loading checkpoint")
  seqtab.nochim <- load_chk("nochim")
} else {
  ts("STEP 5/6: Removing chimeras")
  seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method      = "consensus",
    multithread = threads,
    verbose     = TRUE
  )
  write.csv(seqtab.nochim, file.path(results_dir, "seqtab_nochim.csv"))
  save_chk(seqtab.nochim, "nochim")
}

# -------------------------
# STEP 6/6: Export ASVs
# -------------------------
ts("STEP 6/6: Exporting ASVs —", ncol(seqtab.nochim), "sequences after chimera removal")

asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))

dna <- DNAStringSet(asv_seqs)
names(dna) <- asv_ids
writeXStringSet(dna, filepath = file.path(results_dir, "asvs.nochim.fasta"), format = "fasta")

asv_lookup <- tibble(ASV = asv_ids, sequence = asv_seqs, length = nchar(asv_seqs))
write_tsv(asv_lookup, file.path(results_dir, "asv_lookup.tsv"))

seqtab_asv           <- seqtab.nochim
colnames(seqtab_asv) <- asv_ids

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
# out rownames are basenames of input paths
# -------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out[basename(fnFs_passed), , drop = FALSE],
  denoisedF_counts,
  denoisedR_counts,
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(track, file.path(results_dir, "dada2_read_tracking.csv"))

ts("DADA2 complete")
