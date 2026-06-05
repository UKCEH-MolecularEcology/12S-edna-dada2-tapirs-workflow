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
        ref_db_dir  = config["dada2_ref_db"]
    threads:
        config.get("dada2_threads", 10)
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_RESULTS_DIR="{params.results_dir}" \
        SM_REF_DB_DIR="{params.ref_db_dir}" \
        SM_THREADS="{threads}" \
        Rscript {workflow.basedir}/scripts/dada2/02_dada2_asv.R
        """
