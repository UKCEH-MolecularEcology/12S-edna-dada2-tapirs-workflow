# ==============================================================================
# DADA2 — BLANK CLEANUP (scripts 05 + 06)
# Outputs go to results/dada2/
# ==============================================================================

rule dada2_cleanup_prep:
    input:
        "results/dada2/asv_taxonomy_abundance_CLARE.csv"
    output:
        "results/dada2/ncl_matrix_raw.csv"
    params:
        results_dir = "results/dada2"
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_RESULTS_DIR="{params.results_dir}" \
        Rscript {workflow.basedir}/scripts/dada2/05_make_cleanup_input_from_CLARE.R
        """


rule dada2_blank_cleanup:
    input:
        ncl   = "results/dada2/ncl_matrix_raw.csv",
        clare = "results/dada2/asv_taxonomy_abundance_CLARE.csv"
    output:
        lab  = "results/dada2/ncl_cleaned_labLOD.csv",
        site = "results/dada2/ncl_cleaned_siteLOD.csv",
        both = "results/dada2/ncl_cleaned_bothLOD.csv"
    params:
        results_dir = "results/dada2"
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_RESULTS_DIR="{params.results_dir}" \
        Rscript {workflow.basedir}/scripts/dada2/06_blank_cleanup_from_workflow.R
        """
