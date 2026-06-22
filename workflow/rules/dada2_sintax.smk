# ==============================================================================
# DADA2 — SINTAX TAXONOMY (usearch -sintax, three databases in parallel)
# Outputs go to results/dada2/
# Install usearch first: bash bin/install_usearch.sh
# ==============================================================================

SINTAX_DBS = config.get("sintax_databases", ["INBO", "MIDORI", "CLARE"])

rule dada2_sintax_db:
    """Run usearch SINTAX for one reference database (parallelised across databases)."""
    input:
        fasta  = "results/dada2/asvs.nochim.fasta",
        rds    = "results/dada2/seqtab_asv.rds",
        lookup = "results/dada2/asv_lookup.tsv"
    output:
        raw      = "results/dada2/asv_sintax_{db}.tsv",
        parsed   = "results/dada2/asv_sintax_{db}_parsed.tsv",
        taxonomy = "results/dada2/asv_taxonomy_{db}.tsv",
        species  = "results/dada2/species_abundance_{db}.csv",
        abund    = "results/dada2/asv_taxonomy_abundance_{db}.csv"
    wildcard_constraints:
        db = "|".join(SINTAX_DBS)
    params:
        results_dir = "results/dada2",
        ref_db_dir  = config["dada2_ref_db"],
        cutoff      = config.get("sintax_cutoff", 0.7),
        usearch_bin = config.get("usearch_bin", "bin/usearch")
    threads:
        config.get("sintax_threads", config.get("dada2_threads", 10))
    log:
        "logs/dada2_sintax_{db}.log"
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_RESULTS_DIR="{params.results_dir}" \
        SM_REF_DB_DIR="{params.ref_db_dir}" \
        SM_SINTAX_CUTOFF="{params.cutoff}" \
        SM_THREADS="{threads}" \
        SM_DB_NAME="{wildcards.db}" \
        SM_USEARCH_BIN="{params.usearch_bin}" \
        Rscript {workflow.basedir}/scripts/dada2/03a_sintax_run.R > {log} 2>&1
        """


rule dada2_sintax_merge:
    """Aggregate per-database SINTAX results into comparison and summary tables."""
    input:
        parsed = expand("results/dada2/asv_sintax_{db}_parsed.tsv", db=SINTAX_DBS),
        lookup = "results/dada2/asv_lookup.tsv"
    output:
        compare = "results/dada2/asv_taxonomy_compare.tsv",
        summary = "results/dada2/sintax_database_summary.tsv"
    params:
        results_dir = "results/dada2",
        sintax_dbs  = ",".join(SINTAX_DBS)
    log:
        "logs/dada2_sintax_merge.log"
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_RESULTS_DIR="{params.results_dir}" \
        SM_SINTAX_DBS="{params.sintax_dbs}" \
        Rscript {workflow.basedir}/scripts/dada2/03b_sintax_merge.R > {log} 2>&1
        """
