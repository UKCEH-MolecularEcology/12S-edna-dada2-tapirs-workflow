# ==============================================================================
# QC REPORTS — FastQC on raw reads, MultiQC aggregation, read count summary
# ==============================================================================

rule fastqc:
    input:
        config["input_dir"] + "/{sample}_R{read}_001.fastq.gz"
    output:
        html = "results/qc/fastqc/{sample}_R{read}_001_fastqc.html",
        zip  = "results/qc/fastqc/{sample}_R{read}_001_fastqc.zip"
    wildcard_constraints:
        read = "[12]"
    conda:
        "../envs/qc.yaml"
    shell:
        "mkdir -p results/qc/fastqc && fastqc {input} --outdir results/qc/fastqc/ --quiet"


rule multiqc:
    input:
        fastqc      = expand("results/qc/fastqc/{sample}_R{read}_001_fastqc.zip",
                             sample=SAMPLES, read=["1", "2"]),
        fastp_trim  = expand("results/tapirs/02_trimmed/fastp_trimmed_reports/{combo}.fastp.json",
                             combo=real_combos),
        fastp_merge = expand("results/tapirs/03_merged/fastp_merged_reports/{combo}.merged.fastp.json",
                             combo=real_combos)
    output:
        "results/qc/multiqc_report.html"
    params:
        outdir = "results/qc"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc \
          results/qc/fastqc/ \
          results/tapirs/02_trimmed/fastp_trimmed_reports/ \
          results/tapirs/03_merged/fastp_merged_reports/ \
          -o {params.outdir} \
          --filename multiqc_report.html \
          --force --quiet
        """


rule read_summary:
    input:
        fastp_trim_jsons = expand(
            "results/tapirs/02_trimmed/fastp_trimmed_reports/{combo}.fastp.json",
            combo=real_combos
        ),
        rerep_fastas = expand(
            "results/tapirs/09_rereplicated/{combo}.rerep.fasta",
            combo=real_combos
        ),
        dada2_filtering = "results/dada2/filtering_summary.csv",
        dada2_seqtab    = "results/dada2/seqtab_asv.csv"
    output:
        summary = "results/read_counts_summary.tsv"
    conda:
        "../envs/tapirs.yaml"
    script:
        "../scripts/read_summary.py"
