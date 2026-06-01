# ==============================================================================
# TAPIRS — KRAKEN2 + TAXONOMY
# Outputs go to results/tapirs/
# ==============================================================================

rule kraken2:
    input:
        reads = "results/tapirs/09_rereplicated/{LIBRARIES}/{SAMPLES}.rerep.fasta"
    output:
        reports = "results/tapirs/kraken2/reports/{LIBRARIES}/{SAMPLES}.txt",
        outputs = "results/tapirs/kraken2/outputs/{LIBRARIES}/{SAMPLES}.krk"
    conda:
        "../envs/tapirs.yaml"
    shell:
        """
        if [ -s {input.reads} ]; then
            kraken2 --db {config[kraken2_db]} {input.reads} \
                --threads {config[kraken2_threads]} \
                --confidence {config[kraken2_confidence]} \
                --report {output.reports} \
                --output {output.outputs}
        else
            touch {output.reports}
            touch {output.outputs}
        fi
        """


rule taxonomy_to_kraken2:
    input:
        config['taxdump'] + '/names.dmp',
        kraken2 = "results/tapirs/kraken2/outputs/{LIBRARIES}/{SAMPLES}.krk"
    params:
        taxdump = config['taxdump']
    output:
        kraken2_tax = "results/tapirs/kraken2_tax/{LIBRARIES}/{SAMPLES}.krk.tax.tsv"
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/tapirs/taxonomy_to_kraken2.py"


rule kraken2_to_tsv:
    input:
        lca   = expand("results/tapirs/kraken2_tax/{combo}.krk.tax.tsv",  combo=real_combos),
        rerep = expand("results/tapirs/09_rereplicated/{combo}.rerep.fasta", combo=real_combos)
    params:
        lowest_rank  = config['lowest_taxonomic_rank'],
        highest_rank = config['highest_taxonomic_rank'],
        rerep_dir    = "results/tapirs/09_rereplicated"
    output:
        tsv = "results/tapirs/" + config['my_experiment'] + "_kraken2_conf" + str(config['kraken2_confidence']).split('.')[1] + "_" + config['cluster_method'] + ".tsv"
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/tapirs/mlca_to_tsv.py"


rule kraken2_to_tsv_full_lineage:
    input:
        lca   = expand("results/tapirs/kraken2_tax/{combo}.krk.tax.tsv",  combo=real_combos),
        rerep = expand("results/tapirs/09_rereplicated/{combo}.rerep.fasta", combo=real_combos)
    params:
        lowest_rank  = config['lowest_taxonomic_rank'],
        highest_rank = config['highest_taxonomic_rank'],
        rerep_dir    = "results/tapirs/09_rereplicated"
    output:
        tsv = "results/tapirs/" + config['my_experiment'] + "_kraken2_conf" + str(config['kraken2_confidence']).split('.')[1] + "_" + config['cluster_method'] + "_full_lineage.tsv"
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/tapirs/mlca_to_tsv_full_lineage.py"
