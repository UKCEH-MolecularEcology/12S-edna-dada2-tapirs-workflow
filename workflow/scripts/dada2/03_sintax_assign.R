# =========================
# 03_sintax_assign.R
# vsearch SINTAX classification against three reference databases
# (replaces usearch -sintax with vsearch --sintax; flags use -- prefix)
# =========================

.args       <- commandArgs(trailingOnly = FALSE)
.script_dir <- dirname(normalizePath(sub("--file=", "", .args[grep("--file=", .args)])))
source(file.path(.script_dir, "00_config.R"))

library(tibble)
library(dplyr)
library(readr)
library(stringr)

asv_fasta       <- file.path(results_dir, "asvs.nochim.fasta")
asv_lookup_file <- file.path(results_dir, "asv_lookup.tsv")
seqtab_asv_file <- file.path(results_dir, "seqtab_asv.rds")

if (!file.exists(asv_fasta))       stop("Missing ASV fasta: ",  asv_fasta)
if (!file.exists(asv_lookup_file)) stop("Missing ASV lookup: ", asv_lookup_file)
if (!file.exists(seqtab_asv_file)) stop("Missing ASV table: ",  seqtab_asv_file)

asv_lookup <- read_tsv(asv_lookup_file, show_col_types = FALSE)
seqtab_asv <- readRDS(seqtab_asv_file)

# =========================
# Helper functions
# =========================

resolve_reference <- function(ref_path) {
  if (!file.exists(ref_path)) stop("Reference database not found: ", ref_path)

  if (grepl("\\.gz$", ref_path, ignore.case = TRUE)) {
    out_file <- file.path(
      unzipped_ref_dir,
      sub("\\.gz$", "", basename(ref_path), ignore.case = TRUE)
    )
    if (!file.exists(out_file)) {
      cat("Unzipping reference database:", basename(ref_path), "\n")
      con_in  <- gzfile(ref_path, open = "rb")
      con_out <- file(out_file,   open = "wb")
      repeat {
        bytes <- readBin(con_in, what = raw(), n = 1e6)
        if (length(bytes) == 0) break
        writeBin(bytes, con_out)
      }
      close(con_in)
      close(con_out)
    }
    return(out_file)
  }
  ref_path
}

# vsearch uses --sintax (double dash), unlike usearch which used -sintax
run_sintax <- function(vsearch, asv_fa, db_fa, out_tsv, cutoff = 0.7, threads = 8) {
  args <- c(
    "--sintax",         asv_fa,
    "--db",             db_fa,
    "--tabbedout",      out_tsv,
    "--strand",         "both",
    "--sintax_cutoff",  as.character(cutoff),
    "--threads",        as.character(threads)
  )
  res <- system2(vsearch, args = args, stdout = TRUE, stderr = TRUE)
  cat(paste(res, collapse = "\n"), "\n")
}

read_sintax_simple <- function(path, db_name) {
  x <- read.delim(path, header = FALSE, sep = "\t",
                  stringsAsFactors = FALSE, fill = TRUE)
  if (ncol(x) < 2) stop("Unexpected SINTAX output format in: ", path)
  tibble(ASV = x[[1]], taxonomy = x[[2]], database = db_name)
}

extract_rank_name <- function(x, prefix) {
  m <- str_match(x, paste0(prefix, ":([^,(]+)\\(([^\\)]+)\\)"))
  m[, 2]
}

extract_rank_conf <- function(x, prefix) {
  m <- str_match(x, paste0(prefix, ":([^,(]+)\\(([^\\)]+)\\)"))
  suppressWarnings(as.numeric(m[, 3]))
}

parse_sintax <- function(df) {
  df %>% mutate(
    kingdom      = extract_rank_name(taxonomy, "k"),
    kingdom_conf = extract_rank_conf(taxonomy, "k"),
    phylum       = extract_rank_name(taxonomy, "p"),
    phylum_conf  = extract_rank_conf(taxonomy, "p"),
    class        = extract_rank_name(taxonomy, "c"),
    class_conf   = extract_rank_conf(taxonomy, "c"),
    order        = extract_rank_name(taxonomy, "o"),
    order_conf   = extract_rank_conf(taxonomy, "o"),
    family       = extract_rank_name(taxonomy, "f"),
    family_conf  = extract_rank_conf(taxonomy, "f"),
    genus        = extract_rank_name(taxonomy, "g"),
    genus_conf   = extract_rank_conf(taxonomy, "g"),
    species      = extract_rank_name(taxonomy, "s"),
    species_conf = extract_rank_conf(taxonomy, "s")
  )
}

summarise_assignments <- function(df, db_name) {
  tibble(
    database           = db_name,
    n_asv              = nrow(df),
    assigned_family    = sum(!is.na(df$family)),
    assigned_genus     = sum(!is.na(df$genus)),
    assigned_species   = sum(!is.na(df$species)),
    mean_genus_conf    = mean(df$genus_conf,   na.rm = TRUE),
    mean_species_conf  = mean(df$species_conf, na.rm = TRUE)
  )
}

make_species_abundance <- function(parsed_df, seqtab_asv, out_file) {
  tax_map  <- parsed_df %>% select(ASV, species)
  keep_asv <- tax_map$ASV[!is.na(tax_map$species)]
  if (length(keep_asv) == 0) {
    warning("No species assignments found for ", out_file)
    return(NULL)
  }
  mat  <- seqtab_asv[, keep_asv, drop = FALSE]
  sp   <- tax_map$species[match(colnames(mat), tax_map$ASV)]
  sp_split <- split(seq_len(ncol(mat)), sp)
  sp_mat   <- sapply(sp_split, function(idx) {
    if (length(idx) == 1) mat[, idx] else rowSums(mat[, idx, drop = FALSE])
  })
  if (is.vector(sp_mat)) {
    sp_mat <- matrix(sp_mat, ncol = 1)
    colnames(sp_mat) <- names(sp_split)
    rownames(sp_mat) <- rownames(mat)
  }
  sp_df <- data.frame(sample = rownames(sp_mat), sp_mat, check.names = FALSE)
  write.csv(sp_df, out_file, row.names = FALSE)
  invisible(sp_df)
}

make_asv_taxonomy_abundance <- function(parsed_df, seqtab_asv, out_file) {
  seqtab_df <- data.frame(ASV = colnames(seqtab_asv), t(seqtab_asv), check.names = FALSE)
  merged <- parsed_df %>%
    select(ASV, taxonomy,
           kingdom, kingdom_conf, phylum, phylum_conf,
           class, class_conf, order, order_conf,
           family, family_conf, genus, genus_conf,
           species, species_conf) %>%
    left_join(seqtab_df, by = "ASV")
  write.csv(merged, out_file, row.names = FALSE)
  invisible(merged)
}

# =========================
# Run SINTAX across databases
# =========================

parsed_results <- list()
summary_list   <- list()

for (i in seq_len(nrow(reference_dbs))) {
  db_name  <- reference_dbs$db_name[i]
  db_file  <- file.path(reference_dir, reference_dbs$db_file[i])

  cat("Processing database:", db_name, "\n")
  db_resolved <- resolve_reference(db_file)

  out_tsv <- file.path(results_dir, paste0("asv_sintax_", db_name, ".tsv"))

  run_sintax(
    vsearch = vsearch_bin,
    asv_fa  = asv_fasta,
    db_fa   = db_resolved,
    out_tsv = out_tsv,
    cutoff  = sintax_cutoff,
    threads = threads
  )

  sx   <- read_sintax_simple(out_tsv, db_name)
  sx_p <- parse_sintax(sx)

  write_tsv(sx_p, file.path(results_dir, paste0("asv_sintax_", db_name, "_parsed.tsv")))

  asv_tax <- asv_lookup %>% left_join(sx_p, by = "ASV")
  write_tsv(asv_tax, file.path(results_dir, paste0("asv_taxonomy_", db_name, ".tsv")))

  make_species_abundance(
    parsed_df  = sx_p,
    seqtab_asv = seqtab_asv,
    out_file   = file.path(results_dir, paste0("species_abundance_", db_name, ".csv"))
  )

  make_asv_taxonomy_abundance(
    parsed_df  = sx_p,
    seqtab_asv = seqtab_asv,
    out_file   = file.path(results_dir, paste0("asv_taxonomy_abundance_", db_name, ".csv"))
  )

  parsed_results[[db_name]] <- sx_p
  summary_list[[db_name]]   <- summarise_assignments(sx_p, db_name)
}

# =========================
# Summary and comparison tables
# =========================

summary_tbl <- bind_rows(summary_list)
write_tsv(summary_tbl, file.path(results_dir, "sintax_database_summary.tsv"))
print(summary_tbl)

compare_tbl <- asv_lookup
for (db_name in names(parsed_results)) {
  df <- parsed_results[[db_name]] %>%
    select(
      ASV,
      !!paste0("taxonomy_",    db_name) := taxonomy,
      !!paste0("family_",      db_name) := family,
      !!paste0("genus_",       db_name) := genus,
      !!paste0("species_",     db_name) := species,
      !!paste0("genus_conf_",  db_name) := genus_conf,
      !!paste0("species_conf_",db_name) := species_conf
    )
  compare_tbl <- compare_tbl %>% left_join(df, by = "ASV")
}
write_tsv(compare_tbl, file.path(results_dir, "asv_taxonomy_compare.tsv"))

cat("SINTAX complete\n")
