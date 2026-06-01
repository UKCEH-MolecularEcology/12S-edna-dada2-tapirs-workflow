# ==============================================================================
# TAPIRS — BLAST + TAXONOMY
# Outputs go to results/tapirs/
# ==============================================================================

rule blast:
    input:
        query = "results/tapirs/08_dechimera/{LIBRARIES}/{SAMPLES}.nc.fasta"
    output:
        blast = "results/tapirs/blast/{LIBRARIES}/{SAMPLES}.blast.tsv"
    params:
        outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'",
        evalue    = float(config["BLAST_min_evalue"])
    conda:
        "../envs/tapirs.yaml"
    shell:
        "blastn \
        -query {input.query} \
        -db {config[blast_db]} \
        -outfmt {params.outformat} \
        -perc_identity {config[BLAST_min_perc_ident]} \
        -evalue {params.evalue} \
        -max_target_seqs {config[BLAST_max_target_seqs]} \
        -num_threads {config[BLAST_threads]} \
        -out {output.blast}"


rule taxonomy_to_blast:
    input:
        config['taxdump'] + '/names.dmp',
        blast = "results/tapirs/blast/{LIBRARIES}/{SAMPLES}.blast.tsv"
    params:
        taxdump = config['taxdump']
    output:
        blast_tax = "results/tapirs/blast_tax/{LIBRARIES}/{SAMPLES}.blast.tax.tsv"
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/tapirs/taxonomy_to_blast.py"
