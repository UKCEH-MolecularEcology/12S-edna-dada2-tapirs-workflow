# ==============================================================================
# TAPIRS — VSEARCH (dereplication, denoising/clustering, chimera detection)
# Outputs go to results/tapirs/
# ==============================================================================

rule vsearch_dereplicate:
    input:
        fa = "results/tapirs/05_fasta/{LIBRARIES}/{SAMPLES}.fasta"
    output:
        derep        = "results/tapirs/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta",
        cluster_file = "results/tapirs/08_clusters/derep/{LIBRARIES}/{SAMPLES}.cluster.tsv"
    log:
        "logs/vsearch_dereplicate/{LIBRARIES}/{SAMPLES}.log"
    conda:
        "../envs/tapirs.yaml"
    shell:
        "vsearch \
        --derep_fulllength {input.fa} \
        --sizeout \
        --minuniquesize {config[VSEARCH_minuniqsize]} \
        --output {output.derep} \
        --fasta_width 0 \
        --uc {output.cluster_file} \
        > {log} 2>&1"


rule vsearch_cluster:
    input:
        derep = "results/tapirs/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta"
              if config['cluster_method'] == "cluster" else []
    output:
        cluster      = "results/tapirs/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta",
        cluster_file = "results/tapirs/08_clusters/cluster/{LIBRARIES}/{SAMPLES}.cluster.tsv"
    log:
        "logs/vsearch_cluster/{LIBRARIES}/{SAMPLES}.log"
    conda:
        "../envs/tapirs.yaml"
    shell:
        "vsearch \
        --cluster_fast {input.derep} \
        --sizein --sizeout \
        --query_cov {config[VSEARCH_query_cov]} \
        --id {config[VSEARCH_cluster_id]} \
        --strand both \
        --centroids {output.cluster} \
        --fasta_width 0 \
        --uc {output.cluster_file} \
        > {log} 2>&1"


rule vsearch_denoise:
    input:
        "results/tapirs/06_dereplicated/{LIBRARIES}/{SAMPLES}.derep.fasta"
        if config['cluster_method'] == "denoise" else []
    output:
        seqs            = "results/tapirs/07_denoised/{LIBRARIES}/{SAMPLES}.denoise.fasta",
        denoise_results = "results/tapirs/08_clusters/denoise/{LIBRARIES}/{SAMPLES}.denoise.tsv"
    log:
        "logs/vsearch_denoise/{LIBRARIES}/{SAMPLES}.log"
    conda:
        "../envs/tapirs.yaml"
    shell:
        "vsearch \
        --cluster_unoise {input} \
        --sizein --sizeout \
        --minsize {config[VSEARCH_minsize]} \
        --unoise_alpha {config[VSEARCH_unoise_alpha]} \
        --id {config[VSEARCH_unoise_id]} \
        --centroids {output.seqs} \
        --fasta_width 0 \
        --uc {output.denoise_results} \
        > {log} 2>&1"


if config['chimera_detection'] == "ref":

    rule vsearch_uchime_ref:
        input:
            cluster = "results/tapirs/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta"
                    if config['cluster_method'] == "cluster"
                    else "results/tapirs/07_denoised/{LIBRARIES}/{SAMPLES}.denoise.fasta"
        output:
            nonchimeras = "results/tapirs/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta",
            chimeras    = "results/tapirs/08_dechimera/{LIBRARIES}/{SAMPLES}.chimera.fasta"
        log:
            "logs/vsearch_uchime_ref/{LIBRARIES}/{SAMPLES}.log"
        conda:
            "../envs/tapirs.yaml"
        shell:
            "vsearch \
            --uchime_ref {input.cluster} \
            --db {config[dechim_blast_db]} \
            --chimeras {output.chimeras} \
            --borderline {output.chimeras} \
            --mindiffs {config[VSEARCH_mindiffs]} \
            --mindiv {config[VSEARCH_mindiv]} \
            --fasta_width 0 \
            --nonchimeras {output.nonchimeras} \
            > {log} 2>&1"

elif config['chimera_detection'] == "denovo":

    rule vsearch_uchime3_denovo:
        input:
            cluster = "results/tapirs/07_clustered/{LIBRARIES}/{SAMPLES}.cluster.fasta"
                    if config['cluster_method'] == "cluster"
                    else "results/tapirs/07_denoised/{LIBRARIES}/{SAMPLES}.denoise.fasta"
        output:
            nonchimeras = "results/tapirs/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta",
            chimeras    = "results/tapirs/08_dechimera/{LIBRARIES}/{SAMPLES}.chimera.fasta"
        log:
            "logs/vsearch_uchime3_denovo/{LIBRARIES}/{SAMPLES}.log"
        conda:
            "../envs/tapirs.yaml"
        shell:
            "vsearch \
            --uchime3_denovo {input.cluster} \
            --abskew {config[VSEARCH_abskew]} \
            --chimeras {output.chimeras} \
            --borderline {output.chimeras} \
            --fasta_width 0 \
            --nonchimeras {output.nonchimeras} \
            > {log} 2>&1"


rule vsearch_rereplicate:
    input:
        nonchimeras = "results/tapirs/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta"
    output:
        rerep = "results/tapirs/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta"
    log:
        "logs/vsearch_rereplicate/{LIBRARIES}/{SAMPLES}.log"
    conda:
        "../envs/tapirs.yaml"
    shell:
        "vsearch \
        --rereplicate {input.nonchimeras} \
        --fasta_width 0 \
        --output {output.rerep} \
        > {log} 2>&1"
