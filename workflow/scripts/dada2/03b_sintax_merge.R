# =========================
# 03b_sintax_merge.R
# Aggregate per-database SINTAX results into comparison and summary tables.
# Runs after all 03a_sintax_run.R jobs complete.
# =========================

.args       <- commandArgs(trailingOnly = FALSE)
.script_dir <- dirname(normalizePath(sub("--file=", "", .args[grep("--file=", .args)])))
source(file.path(.script_dir, "00_config.R"))

library(tibble)
library(dplyr)
library(readr)

ts <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")

db_names_str <- Sys.getenv("SM_SINTAX_DBS")
if (db_names_str == "") stop("SM_SINTAX_DBS not set")
db_names <- trimws(strsplit(db_names_str, ",")[[1]])

asv_lookup_file <- file.path(results_dir, "asv_lookup.tsv")
if (!file.exists(asv_lookup_file)) stop("Missing ASV lookup: ", asv_lookup_file)
asv_lookup <- read_tsv(asv_lookup_file, show_col_types = FALSE)

ts("Merging SINTAX results — databases:", paste(db_names, collapse = ", "))

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

parsed_results <- list()
summary_list   <- list()

for (db_name in db_names) {
  parsed_file <- file.path(results_dir, paste0("asv_sintax_", db_name, "_parsed.tsv"))
  if (!file.exists(parsed_file)) stop("Missing parsed SINTAX file: ", parsed_file)
  sx_p <- read_tsv(parsed_file, show_col_types = FALSE)
  parsed_results[[db_name]] <- sx_p
  summary_list[[db_name]]   <- summarise_assignments(sx_p, db_name)
  ts(sprintf("  Loaded %s: %d ASVs, %d assigned to species",
             db_name, nrow(sx_p), summary_list[[db_name]]$assigned_species))
}

summary_tbl <- bind_rows(summary_list)
write_tsv(summary_tbl, file.path(results_dir, "sintax_database_summary.tsv"))
print(summary_tbl)

compare_tbl <- asv_lookup
for (db_name in names(parsed_results)) {
  df <- parsed_results[[db_name]] %>%
    select(
      ASV,
      !!paste0("taxonomy_",     db_name) := taxonomy,
      !!paste0("family_",       db_name) := family,
      !!paste0("genus_",        db_name) := genus,
      !!paste0("species_",      db_name) := species,
      !!paste0("genus_conf_",   db_name) := genus_conf,
      !!paste0("species_conf_", db_name) := species_conf
    )
  compare_tbl <- compare_tbl %>% left_join(df, by = "ASV")
}
write_tsv(compare_tbl, file.path(results_dir, "asv_taxonomy_compare.tsv"))

ts("SINTAX merge complete")
