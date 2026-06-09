# =========================
# 06_blank_cleanup_from_workflow.R
# Blank cleanup: LOD thresholds from PCR, extraction, and site blanks.
# Fish taxa are intentionally kept.
# =========================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)

results_dir <- Sys.getenv("SM_RESULTS_DIR", unset = "results")

clare_file <- file.path(results_dir, "asv_taxonomy_abundance_CLARE.csv")
ncl_file   <- file.path(results_dir, "ncl_matrix_raw.csv")

if (!file.exists(clare_file)) stop("Missing CLARE output: ", clare_file)
if (!file.exists(ncl_file))   stop("Missing cleanup matrix: ", ncl_file,
                                   ". Run 05_make_cleanup_input_from_CLARE.R first.")

# ============================================================
# 1. Read files
# ============================================================

clare <- read_csv(clare_file, show_col_types = FALSE)
ncl   <- read_csv(ncl_file,   show_col_types = FALSE)

sample_cols <- setdiff(names(ncl), c("species", "taxonomy"))

# ============================================================
# 2. Classify columns by blank type
# ============================================================

classify_sample_type <- function(nm) {
  nm_u <- toupper(nm)
  if (grepl("PCR_BLANK|PCRBLANK|NEG|POS",                   nm_u)) return("pcr_blank")
  if (grepl("EXTRACTION_BLANK|EXTRACTIONBLANK|_EB_|_EB[0-9]", nm_u)) return("extraction_blank")
  if (grepl("BLANK",                                          nm_u)) return("site_blank")
  "sample"
}

get_site_code <- function(nm) {
  # Site blanks (SPGBLANK, Sylphium_blank) are field/spring blanks with no
  # per-site match. All samples share a single "global" pool.
  "global"
}

sample_info <- tibble(
  sample      = sample_cols,
  sample_type = map_chr(sample_cols, classify_sample_type),
  site        = map_chr(sample_cols, get_site_code)
)

write_csv(sample_info, file.path(results_dir, "cleanup_sample_classification.csv"))

# ============================================================
# 3. Lowest accepted taxonomy + confidence (from CLARE ASV output)
# ============================================================

lowest_tax <- clare %>%
  mutate(across(any_of(c("species","genus","family","order","class","phylum","kingdom")),
                ~ na_if(.x, ""))) %>%
  transmute(
    taxonomy,
    lowest_rank = case_when(
      !is.na(species) ~ "species",
      !is.na(genus)   ~ "genus",
      !is.na(family)  ~ "family",
      !is.na(order)   ~ "order",
      !is.na(class)   ~ "class",
      !is.na(phylum)  ~ "phylum",
      !is.na(kingdom) ~ "kingdom",
      TRUE            ~ NA_character_
    ),
    lowest_name = case_when(
      !is.na(species) ~ species,
      !is.na(genus)   ~ genus,
      !is.na(family)  ~ family,
      !is.na(order)   ~ order,
      !is.na(class)   ~ class,
      !is.na(phylum)  ~ phylum,
      !is.na(kingdom) ~ kingdom,
      TRUE            ~ NA_character_
    ),
    lowest_conf = case_when(
      !is.na(species) ~ species_conf,
      !is.na(genus)   ~ genus_conf,
      !is.na(family)  ~ family_conf,
      !is.na(order)   ~ order_conf,
      !is.na(class)   ~ class_conf,
      !is.na(phylum)  ~ phylum_conf,
      !is.na(kingdom) ~ kingdom_conf,
      TRUE            ~ NA_real_
    )
  ) %>%
  distinct()

ncl_tax <- ncl %>%
  select(species, taxonomy) %>%
  left_join(lowest_tax, by = "taxonomy")

write_csv(ncl_tax, file.path(results_dir, "cleanup_taxonomy_lookup.csv"))

# ============================================================
# 4. Long format
# ============================================================

ncl_long <- ncl %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "reads") %>%
  mutate(reads = as.numeric(reads)) %>%
  left_join(sample_info, by = "sample") %>%
  left_join(ncl_tax,     by = c("species", "taxonomy"))

# ============================================================
# 5. Blank statistics (lab + site)
# ============================================================

lab_blank_stats <- ncl_long %>%
  filter(sample_type %in% c("pcr_blank", "extraction_blank")) %>%
  group_by(species) %>%
  summarise(
    mean_lab_blank_reads = mean(reads, na.rm = TRUE),
    sd_lab_blank_reads   = sd(reads,   na.rm = TRUE),
    max_lab_blank_reads  = max(reads,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_lab_blank_reads = ifelse(is.na(sd_lab_blank_reads), 0, sd_lab_blank_reads),
    lab_LOD = mean_lab_blank_reads + 3  * sd_lab_blank_reads,
    lab_LOQ = mean_lab_blank_reads + 10 * sd_lab_blank_reads
  )

site_blank_stats <- ncl_long %>%
  filter(sample_type == "site_blank", !is.na(site)) %>%
  group_by(species, site) %>%
  summarise(
    mean_site_blank_reads = mean(reads, na.rm = TRUE),
    sd_site_blank_reads   = sd(reads,   na.rm = TRUE),
    max_site_blank_reads  = max(reads,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_site_blank_reads = ifelse(is.na(sd_site_blank_reads), 0, sd_site_blank_reads),
    site_LOD = mean_site_blank_reads + 3  * sd_site_blank_reads,
    site_LOQ = mean_site_blank_reads + 10 * sd_site_blank_reads
  )

write_csv(lab_blank_stats,  file.path(results_dir, "cleanup_lab_blank_stats.csv"))
write_csv(site_blank_stats, file.path(results_dir, "cleanup_site_blank_stats.csv"))

# ============================================================
# 6. Diagnostic table
# ============================================================

diagnostic_tbl <- ncl_long %>%
  filter(sample_type == "sample") %>%
  group_by(species, taxonomy, lowest_rank, lowest_name, lowest_conf, site) %>%
  summarise(
    max_sample_reads   = max(reads,  na.rm = TRUE),
    mean_sample_reads  = mean(reads, na.rm = TRUE),
    n_positive_samples = sum(reads > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(site_blank_stats, by = c("species", "site")) %>%
  left_join(lab_blank_stats,  by = "species") %>%
  mutate(
    site_blank_to_sample_ratio = ifelse(max_sample_reads > 0,
                                        max_site_blank_reads / max_sample_reads, NA_real_),
    lab_blank_to_sample_ratio  = ifelse(max_sample_reads > 0,
                                        max_lab_blank_reads  / max_sample_reads, NA_real_)
  ) %>%
  arrange(desc(site_blank_to_sample_ratio), desc(lab_blank_to_sample_ratio))

write_csv(diagnostic_tbl, file.path(results_dir, "cleanup_blank_diagnostic_table.csv"))

# ============================================================
# 7. Apply cleanup to true samples
# ============================================================

sample_long <- ncl_long %>%
  filter(sample_type == "sample") %>%
  left_join(lab_blank_stats  %>% select(species, lab_LOD),        by = "species") %>%
  left_join(site_blank_stats %>% select(species, site, site_LOD), by = c("species", "site")) %>%
  mutate(
    lab_LOD          = ifelse(is.na(lab_LOD),  0, lab_LOD),
    site_LOD         = ifelse(is.na(site_LOD), 0, site_LOD),
    reads_lab_clean  = ifelse(reads < lab_LOD,                   0, reads),
    reads_site_clean = ifelse(reads < site_LOD,                  0, reads),
    reads_both_clean = ifelse(reads < pmax(lab_LOD, site_LOD),   0, reads)
  )

# ============================================================
# 8. Remove human and domestic livestock; fish are intentionally kept.
# ============================================================
exclude_taxa    <- c("Homo_sapiens", "Sus_scrofa", "Canis_lupus",
                     "Ovis_aries", "Bos_taurus", "Bovidae")
exclude_pattern <- paste(exclude_taxa, collapse = "|")
sample_long     <- sample_long %>%
  filter(!grepl(exclude_pattern, taxonomy, ignore.case = TRUE))

# ============================================================
# 9. Write cleaned matrices
# ============================================================

make_matrix <- function(long_df, value_col, out_file) {
  out <- long_df %>%
    select(species, taxonomy, lowest_rank, lowest_name, lowest_conf,
           sample, !!sym(value_col)) %>%
    pivot_wider(names_from = sample, values_from = !!sym(value_col), values_fill = 0)
  write_csv(out, out_file)
  invisible(out)
}

make_matrix(sample_long, "reads_lab_clean",  file.path(results_dir, "ncl_cleaned_labLOD.csv"))
make_matrix(sample_long, "reads_site_clean", file.path(results_dir, "ncl_cleaned_siteLOD.csv"))
make_matrix(sample_long, "reads_both_clean", file.path(results_dir, "ncl_cleaned_bothLOD.csv"))

# ============================================================
# 10. Presence / absence versions
# ============================================================

make_pa <- function(infile, outfile) {
  x           <- read_csv(infile, show_col_types = FALSE)
  meta        <- c("species","taxonomy","lowest_rank","lowest_name","lowest_conf")
  samp_cols   <- setdiff(names(x), meta)
  x[samp_cols] <- lapply(x[samp_cols], function(z) as.integer(z > 0))
  write_csv(x, outfile)
}

make_pa(file.path(results_dir, "ncl_cleaned_labLOD.csv"),
        file.path(results_dir, "ncl_cleaned_labLOD_pa.csv"))
make_pa(file.path(results_dir, "ncl_cleaned_siteLOD.csv"),
        file.path(results_dir, "ncl_cleaned_siteLOD_pa.csv"))
make_pa(file.path(results_dir, "ncl_cleaned_bothLOD.csv"),
        file.path(results_dir, "ncl_cleaned_bothLOD_pa.csv"))

# ============================================================
# 11. Long format output
# ============================================================

ncl_long_out <- sample_long %>%
  select(species, taxonomy, lowest_rank, lowest_name, lowest_conf,
         sample, site, reads, reads_lab_clean, reads_site_clean, reads_both_clean)

write_csv(ncl_long_out, file.path(results_dir, "ncl_cleaned_long.csv"))

# ============================================================
# 12. Summary
# ============================================================

summary_tbl <- sample_long %>%
  summarise(
    total_raw_reads          = sum(reads,             na.rm = TRUE),
    total_lab_clean_reads    = sum(reads_lab_clean,   na.rm = TRUE),
    total_site_clean_reads   = sum(reads_site_clean,  na.rm = TRUE),
    total_both_clean_reads   = sum(reads_both_clean,  na.rm = TRUE),
    n_nonzero_raw            = sum(reads > 0,             na.rm = TRUE),
    n_nonzero_lab_clean      = sum(reads_lab_clean > 0,   na.rm = TRUE),
    n_nonzero_site_clean     = sum(reads_site_clean > 0,  na.rm = TRUE),
    n_nonzero_both_clean     = sum(reads_both_clean > 0,  na.rm = TRUE)
  )

write_csv(summary_tbl, file.path(results_dir, "cleanup_summary.csv"))

cat("Cleanup complete\n")
