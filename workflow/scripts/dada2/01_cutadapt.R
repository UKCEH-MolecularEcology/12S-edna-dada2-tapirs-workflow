# =========================
# 01_cutadapt.R
# Primer trimming with cutadapt
# =========================

# Locate and source config from same directory as this script (works with Rscript /abs/path)
.args       <- commandArgs(trailingOnly = FALSE)
.script_dir <- dirname(normalizePath(sub("--file=", "", .args[grep("--file=", .args)])))
source(file.path(.script_dir, "00_config.R"))

revcomp <- function(x) {
  comp <- c(
    A = "T", C = "G", G = "C", T = "A",
    R = "Y", Y = "R", S = "S", W = "W",
    K = "M", M = "K", B = "V", V = "B",
    D = "H", H = "D", N = "N"
  )
  x <- toupper(gsub("\\s+", "", x))
  paste(rev(comp[strsplit(x, "")[[1]]]), collapse = "")
}

FWD_rc <- revcomp(FWD)
REV_rc <- revcomp(REV)

r1 <- list.files(
  raw_dir,
  pattern    = "_R1_001.*\\.(fastq|fq)(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

r2   <- sub("_R1_001", "_R2_001", r1, fixed = TRUE)
keep <- file.exists(r2)

if (!any(keep)) stop("No matching R1/R2 pairs found in: ", raw_dir)

r1 <- r1[keep]
r2 <- r2[keep]

trim_left <- as.integer(Sys.getenv("SM_TRIM_LEFT", unset = "0"))

for (i in seq_along(r1)) {
  cat("Running cutadapt on:", basename(r1[i]), "\n")

  out_r1 <- file.path(cutadapt_dir, basename(r1[i]))
  out_r2 <- file.path(cutadapt_dir, basename(r2[i]))

  if (trim_left > 0) {
    # Fixed 5' trimming — remove primer by position, no primer detection
    args <- c(
      "--cores",          as.character(threads),
      "--minimum-length", "75",
      "-u",  as.character(trim_left),
      "-U",  as.character(trim_left),
      "-o",  out_r1,
      "-p",  out_r2,
      r1[i], r2[i]
    )
  } else {
    # Primer detection via sequence matching
    args <- c(
      "--cores",           as.character(threads),
      "--discard-untrimmed",
      "--pair-filter=any",
      "--minimum-length",  "75",
      "-g",  FWD,
      "-G",  REV,
      "-a",  REV_rc,
      "-A",  FWD_rc,
      "-o",  out_r1,
      "-p",  out_r2,
      r1[i], r2[i]
    )
  }

  res <- system2(cutadapt_bin, args = args, stdout = TRUE, stderr = TRUE)
  cat(paste(res, collapse = "\n"), "\n")
}

cat("cutadapt complete\n")
