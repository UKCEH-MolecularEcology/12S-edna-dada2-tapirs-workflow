# ==============================================================================
# DADA2 — ASV INFERENCE
# Outputs go to results/dada2/
# ==============================================================================

rule dada2_asv:
    input:
        "results/dada2/.cutadapt.done"
    output:
        fasta      = "results/dada2/asvs.nochim.fasta",
        rds        = "results/dada2/seqtab_asv.rds",
        lookup     = "results/dada2/asv_lookup.tsv",
        seqtab_csv = "results/dada2/seqtab_asv.csv",
        filtering  = "results/dada2/filtering_summary.csv"
    params:
        results_dir = "results/dada2",
        ref_db_dir  = config["dada2_ref_db"],
        maxEE_F     = config.get("dada2_maxEE_F", 2),
        maxEE_R     = config.get("dada2_maxEE_R", 5)
    threads:
        config.get("dada2_threads", 10)
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_RESULTS_DIR="{params.results_dir}" \
        SM_REF_DB_DIR="{params.ref_db_dir}" \
        SM_THREADS="{threads}" \
        SM_MAXEE_F="{params.maxEE_F}" \
        SM_MAXEE_R="{params.maxEE_R}" \
        Rscript {workflow.basedir}/scripts/dada2/02_dada2_asv.R
        """
