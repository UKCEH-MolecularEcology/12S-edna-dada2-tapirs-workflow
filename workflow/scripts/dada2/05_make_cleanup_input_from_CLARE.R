# =========================
# 05_make_cleanup_input_from_CLARE.R
# Build taxon-level matrix from CLARE SINTAX output
# =========================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

results_dir <- Sys.getenv("SM_RESULTS_DIR", unset = "results")

infile  <- file.path(results_dir, "asv_taxonomy_abundance_CLARE.csv")
outfile <- file.path(results_dir, "ncl_matrix_raw.csv")

if (!file.exists(infile)) stop("Missing input file: ", infile)

x <- read_csv(infile, show_col_types = FALSE)

meta_cols <- c(
  "ASV", "taxonomy",
  "kingdom", "kingdom_conf",
  "phylum",  "phylum_conf",
  "class",   "class_conf",
  "order",   "order_conf",
  "family",  "family_conf",
  "genus",   "genus_conf",
  "species", "species_conf"
)

sample_cols <- setdiff(names(x), meta_cols)

# Assign lowest available taxonomic label, preferring species
x2 <- x %>%
  mutate(
    species = na_if(species, ""),
    genus   = na_if(genus,   ""),
    family  = na_if(family,  ""),
    order   = na_if(order,   ""),
    class   = na_if(class,   ""),
    phylum  = na_if(phylum,  ""),
    kingdom = na_if(kingdom, ""),
    taxon_label = case_when(
      !is.na(species) ~ species,
      !is.na(genus)   ~ genus,
      !is.na(family)  ~ family,
      !is.na(order)   ~ order,
      !is.na(class)   ~ class,
      !is.na(phylum)  ~ phylum,
      !is.na(kingdom) ~ kingdom,
      TRUE ~ ASV
    )
  )

# Collapse ASVs sharing the same label
ncl_raw <- x2 %>%
  group_by(taxon_label) %>%
  summarise(
    taxonomy = dplyr::first(na.omit(taxonomy)),
    across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(species = taxon_label)

write_csv(ncl_raw, outfile)
cat("Wrote: ", outfile, "\n", sep = "")
