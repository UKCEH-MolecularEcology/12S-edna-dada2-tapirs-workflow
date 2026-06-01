# ==============================================================================
# DADA2 — PRIMER TRIMMING (cutadapt)
# Outputs go to results/dada2/cutadapt/
# ==============================================================================

_d2_results = "results/dada2"

rule dada2_cutadapt:
    input:
        r1 = expand(config["input_dir"] + "/{sample}_R1_001.fastq.gz", sample=SAMPLES)
    output:
        done = touch("results/dada2/.cutadapt.done")
    params:
        input_dir   = config["input_dir"],
        results_dir = _d2_results,
        ref_db_dir  = config["dada2_ref_db"]
    threads:
        config.get("dada2_threads", 10)
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_INPUT_DIR="{params.input_dir}" \
        SM_RESULTS_DIR="{params.results_dir}" \
        SM_REF_DB_DIR="{params.ref_db_dir}" \
        SM_THREADS="{threads}" \
        Rscript {workflow.basedir}/scripts/dada2/01_cutadapt.R
        """
