# ==============================================================================
# TAPIRS — MLCA / LCA
# Outputs go to results/tapirs/
# ==============================================================================

rule mlca:
    input:
        blast = "results/tapirs/blast_tax/{LIBRARIES}/{SAMPLES}.blast.tax.tsv"
    output:
        lca = "results/tapirs/mlca/{LIBRARIES}/{SAMPLES}.lca.tsv"
    params:
        bitscore = config['MLCA_bitscore'],
        identity = config['MLCA_identity'],
        coverage = config['MCLA_coverage'],
        majority = config['MLCA_majority'],
        min_hits = config['MLCA_hits']
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/tapirs/mlca.py"


rule mlca_to_tsv:
    input:
        lca   = expand("results/tapirs/mlca/{combo}.lca.tsv",              combo=real_combos),
        rerep = expand("results/tapirs/09_rereplicated/{combo}.rerep.fasta", combo=real_combos)
    params:
        lowest_rank  = config['lowest_taxonomic_rank'],
        highest_rank = config['highest_taxonomic_rank'],
        rerep_dir    = "results/tapirs/09_rereplicated"
    output:
        tsv = "results/tapirs/" + config['my_experiment'] + "_blast" + str(config['MLCA_identity']) + "_" + config['cluster_method'] + ".tsv"
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/tapirs/mlca_to_tsv.py"


rule mlca_to_tsv_full_lineage:
    input:
        lca   = expand("results/tapirs/mlca/{combo}.lca.tsv",              combo=real_combos),
        rerep = expand("results/tapirs/09_rereplicated/{combo}.rerep.fasta", combo=real_combos)
    params:
        lowest_rank  = config['lowest_taxonomic_rank'],
        highest_rank = config['highest_taxonomic_rank'],
        rerep_dir    = "results/tapirs/09_rereplicated"
    output:
        tsv = "results/tapirs/" + config['my_experiment'] + "_blast" + str(config['MLCA_identity']) + "_" + config['cluster_method'] + "_full_lineage.tsv"
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/tapirs/mlca_to_tsv_full_lineage.py"
