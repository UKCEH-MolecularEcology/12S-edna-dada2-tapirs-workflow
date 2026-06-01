# ==============================================================================
# DADA2 — SINTAX TAXONOMY (vsearch --sintax, three databases)
# Outputs go to results/dada2/
# ==============================================================================

rule dada2_sintax:
    input:
        fasta  = "results/dada2/asvs.nochim.fasta",
        rds    = "results/dada2/seqtab_asv.rds",
        lookup = "results/dada2/asv_lookup.tsv"
    output:
        compare  = "results/dada2/asv_taxonomy_compare.tsv",
        summary  = "results/dada2/sintax_database_summary.tsv",
        clare    = "results/dada2/asv_taxonomy_abundance_CLARE.csv"
    params:
        results_dir = "results/dada2",
        ref_db_dir  = config["dada2_ref_db"],
        cutoff      = config.get("sintax_cutoff", 0.7)
    threads:
        config.get("dada2_threads", 10)
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_RESULTS_DIR="{params.results_dir}" \
        SM_REF_DB_DIR="{params.ref_db_dir}" \
        SM_SINTAX_CUTOFF="{params.cutoff}" \
        SM_THREADS="{threads}" \
        Rscript {workflow.basedir}/scripts/dada2/03_sintax_assign.R
        """
