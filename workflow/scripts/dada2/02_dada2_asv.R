# =========================
# 02_dada2_asv.R
# DADA2 pipeline and ASV export
# Matches Lauren's script structure: list-based dada(), default learnErrors nbases.
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
  pattern    = "_R1_001.*\\.(fastq|fq)(\\.gz)?$",
  full.names = TRUE
))
fnRs <- sort(list.files(
  cutadapt_dir,
  pattern    = "_R2_001.*\\.(fastq|fq)(\\.gz)?$",
  full.names = TRUE
))

if (length(fnFs) == 0 || length(fnRs) == 0) {
  stop("No cutadapt output files found in: ", cutadapt_dir)
}

sample.names <- sub("_R1_001.*$", "", basename(fnFs))

filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_dir, paste0(sample.names, "_R_filt.fastq.gz"))

# -------------------------
# STEP 1/5: Quality filter and trim
# -------------------------
if (has_chk("filter")) {
  ts("STEP 1/5: Filtering — loading checkpoint")
  chk          <- load_chk("filter")
  out          <- chk$out
  fnFs_passed  <- chk$fnFs_passed
  sample.names <- chk$sample.names
  filtFs       <- chk$filtFs
  filtRs       <- chk$filtRs
} else {
  ts("STEP 1/5: Filtering and trimming reads (", length(fnFs), "samples,", threads, "threads)")
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
# STEP 2/5: Learn errors
# -------------------------
if (has_chk("errors")) {
  ts("STEP 2/5: Error rates — loading checkpoint")
  errs <- load_chk("errors")
  errF <- errs$errF
  errR <- errs$errR
} else {
  ts("STEP 2/5: Learning error rates (", length(filtFs), "samples,", threads, "threads)")
  errF <- learnErrors(filtFs, nbases = 1e8, multithread = threads)
  errR <- learnErrors(filtRs, nbases = 1e8, multithread = threads)
  ts("  Error rates learned")
  save_chk(list(errF=errF, errR=errR), "errors")
}

# -------------------------
# STEP 3/5: Dereplicate, denoise, merge (all samples at once, matching Lauren's approach)
# -------------------------
if (has_chk("denoised")) {
  ts("STEP 3/5: Denoising — loading checkpoint")
  dn      <- load_chk("denoised")
  mergers <- dn$mergers
} else {
  ts("STEP 3/5: Dereplicating", length(sample.names), "samples")
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

  derepFs <- derepFastq(filtFs, verbose = TRUE)
  derepRs <- derepFastq(filtRs, verbose = TRUE)
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names

  ts("STEP 3/5: Running DADA2 (forward)")
  dadaFs <- dada(derepFs, err = errF, multithread = threads)
  ts("STEP 3/5: Running DADA2 (reverse)")
  dadaRs <- dada(derepRs, err = errR, multithread = threads)

  ts("STEP 3/5: Merging pairs")
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

  save_chk(list(mergers=mergers), "denoised")
}

# -------------------------
# STEP 4/5: Sequence table and chimera removal
# -------------------------
if (has_chk("nochim")) {
  ts("STEP 4/5: Sequence table + chimera removal — loading checkpoint")
  chk2          <- load_chk("nochim")
  seqtab        <- chk2$seqtab
  seqtab.nochim <- chk2$seqtab.nochim
} else {
  ts("STEP 4/5: Making sequence table")
  seqtab <- makeSequenceTable(mergers)
  write.csv(seqtab, file.path(results_dir, "seqtab_raw.csv"))
  ts("  Sequence table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs")

  ts("STEP 4/5: Removing chimeras")
  seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method      = "consensus",
    multithread = threads,
    verbose     = TRUE
  )
  write.csv(seqtab.nochim, file.path(results_dir, "seqtab_nochim.csv"))

  save_chk(list(seqtab=seqtab, seqtab.nochim=seqtab.nochim), "nochim")
}

# -------------------------
# STEP 5/5: Export ASVs
# -------------------------
ts("STEP 5/5: Exporting ASVs —", ncol(seqtab.nochim), "sequences after chimera removal")

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

ts("DADA2 complete")
