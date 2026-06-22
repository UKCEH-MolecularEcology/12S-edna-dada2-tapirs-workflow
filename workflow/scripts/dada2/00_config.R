# =========================
# 00_config.R
# Central configuration — reads from environment variables when run via
# Snakemake shell rules, or falls back to standalone defaults.
#
# Env vars (set by Snakemake rules):
#   SM_INPUT_DIR      path to raw FASTQ directory
#   SM_RESULTS_DIR    path to results output directory  (results/dada2)
#   SM_REF_DB_DIR     path to SINTAX reference database directory
#   SM_THREADS        number of threads
#   SM_SINTAX_CUTOFF  SINTAX confidence cutoff (default 0.7)
# =========================

# Primers (12S Riaz V5 region)
FWD_base <- "ACTGGGATTAGATACCCC"
REV <- "TAGAACAGGCTCCTCTAG"

# If fusion tag was used, first base of FWD primer is A or G (degenerate R)
fusion_tag <- Sys.getenv("SM_FUSION_TAG", unset = "no")
FWD <- if (tolower(fusion_tag) == "yes") {
  paste0("R", substring(FWD_base, 2))
} else {
  FWD_base
}

# Threads
threads <- as.integer(
  Sys.getenv("SM_THREADS", unset = as.character(max(1L, min(10L, parallel::detectCores()))))
)

# Directories
raw_dir       <- Sys.getenv("SM_INPUT_DIR",   unset = "RawSeqs")
results_dir   <- Sys.getenv("SM_RESULTS_DIR", unset = "results")
reference_dir <- Sys.getenv("SM_REF_DB_DIR",  unset = "reference_db")

cutadapt_dir     <- file.path(results_dir, "cutadapt")
filtered_dir     <- file.path(results_dir, "filtered")
unzipped_ref_dir <- file.path(results_dir, "reference_db_unzipped")

dir.create(results_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(cutadapt_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(filtered_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(unzipped_ref_dir, showWarnings = FALSE, recursive = TRUE)

# External tools
cutadapt_bin <- "cutadapt"

# usearch: prefer SM_USEARCH_BIN env var, then bin/usearch in repo root,
# then /hdd0/susbus/tools/usearch, then fall back to PATH
.usearch_candidates <- c(
  Sys.getenv("SM_USEARCH_BIN", unset = ""),
  file.path(.script_dir, "..", "..", "..", "bin", "usearch"),
  "/hdd0/susbus/tools/usearch",
  "usearch"
)
usearch_bin <- Filter(function(p) p != "" && (file.exists(p) || p == "usearch"),
                      .usearch_candidates)[[1]]

# SINTAX confidence cutoff
sintax_cutoff <- as.numeric(Sys.getenv("SM_SINTAX_CUTOFF", unset = "0.7"))

# Reference databases (SINTAX-formatted FASTA)
reference_dbs <- data.frame(
  db_name = c("INBO", "MIDORI", "CLARE"),
  db_file = c(
    "all_seqs_INBO_riaz_amplified.sintax.fasta",
    "MIDORI2_LONGEST_NUC_GB269_srRNA_SINTAX.fasta.gz",
    "12s_verts.trimmed_RiazV5.sintax.fasta"
  ),
  stringsAsFactors = FALSE
)
