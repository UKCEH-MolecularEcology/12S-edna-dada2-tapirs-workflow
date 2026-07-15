# ==============================================================================
# DADA2 — PRIMER TRIMMING (cutadapt)
# Outputs go to results/dada2/cutadapt/
# ==============================================================================

_d2_results = "results/dada2"

rule dada2_cutadapt:
    input:
        r1 = expand(config["input_dir"] + "/{sample}" + R1_SUFFIX, sample=SAMPLES)
    output:
        done = touch("results/dada2/.cutadapt.done")
    params:
        input_dir   = config["input_dir"],
        results_dir = _d2_results,
        ref_db_dir  = config["dada2_ref_db"],
        fusion_tag = config.get("fusion_tag", "no"),
        trim_left  = config.get("trim_left", 0),
        primers_already_removed = config.get("primers_already_removed", "no")
    threads:
        config.get("dada2_threads", 10)
    log:
        "logs/dada2_cutadapt.log"
    conda:
        "../envs/dada2.yaml"
    shell:
        """
        SM_INPUT_DIR="{params.input_dir}" \
        SM_RESULTS_DIR="{params.results_dir}" \
        SM_REF_DB_DIR="{params.ref_db_dir}" \
        SM_THREADS="{threads}" \
        SM_FUSION_TAG="{params.fusion_tag}" \
        SM_TRIM_LEFT="{params.trim_left}" \
        SM_PRIMERS_ALREADY_REMOVED="{params.primers_already_removed}" \
        Rscript {workflow.basedir}/scripts/dada2/01_cutadapt.R > {log} 2>&1
        """
