# ==============================================================================
# TAPIRS — QUALITY CONTROL
# Outputs go to results/tapirs/
# ==============================================================================

rule fastp_trim_reads:
    input:
        read1 = config['input_data'] + "/{LIBRARIES}/{SAMPLES}" + R1_SUFFIX,
        read2 = config['input_data'] + "/{LIBRARIES}/{SAMPLES}" + R2_SUFFIX
    output:
        R1trimmed  = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.trimmed.fastq",
        R2trimmed  = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.R2.trimmed.fastq",
        R1unpaired = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.unpaired.fastq",
        R2unpaired = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.R2.unpaired.fastq",
        failed     = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.trimmed.failed.fastq",
        json       = "results/tapirs/02_trimmed/fastp_trimmed_reports/{LIBRARIES}/{SAMPLES}.fastp.json",
        html       = "results/tapirs/02_trimmed/fastp_trimmed_reports/{LIBRARIES}/{SAMPLES}.fastp.html"
    log:
        "logs/fastp_trim/{LIBRARIES}/{SAMPLES}.log"
    conda:
        "../envs/tapirs.yaml"
    shell:
        "fastp \
        --disable_adapter_trimming \
        --in1 {input.read1} \
        --in2 {input.read2} \
        --out1 {output.R1trimmed} \
        --out2 {output.R2trimmed} \
        --unpaired1 {output.R1unpaired} \
        --unpaired2 {output.R2unpaired} \
        --failed_out {output.failed} \
        -j {output.json} \
        -h {output.html} \
        --qualified_quality_phred {config[FASTP_qual_phred]} \
        --unqualified_percent_limit {config[Fastp_unqualified_percent_limit]} \
        --average_qual {config[FASTP_qual_phred]} \
        --cut_tail \
        --cut_window_size {config[FASTP_window_size]} \
        --cut_mean_quality {config[FASTP_qual_phred]} \
        --trim_poly_g \
        --trim_poly_x \
        --poly_g_min_len {config[FASTP_poly_g_min]} \
        --poly_x_min_len {config[FASTP_poly_x_min]} \
        --length_required {config[FASTP_len_required]} \
        --trim_front1 {config[FASTP_trim_front1]} \
        --trim_front2 {config[FASTP_trim_front2]} \
        --overlap_diff_percent_limit {config[FASTP_diff_percent_limit]} \
        --max_len1 {config[FASTP_max_len1]} \
        --max_len2 {config[FASTP_max_len2]} \
        > {log} 2>&1"


rule fastp_merge_reads:
    input:
        R1trimmed = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.trimmed.fastq",
        R2trimmed = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.R2.trimmed.fastq"
    output:
        merged    = "results/tapirs/03_merged/{LIBRARIES}/{SAMPLES}.merged.fastq",
        R1unmerged = "results/tapirs/03_merged/{LIBRARIES}/{SAMPLES}.R1.unmerged.fastq",
        R2unmerged = "results/tapirs/03_merged/{LIBRARIES}/{SAMPLES}.R2.unmerged.fastq",
        json      = "results/tapirs/03_merged/fastp_merged_reports/{LIBRARIES}/{SAMPLES}.merged.fastp.json",
        html      = "results/tapirs/03_merged/fastp_merged_reports/{LIBRARIES}/{SAMPLES}.merged.fastp.html"
    log:
        "logs/fastp_merge/{LIBRARIES}/{SAMPLES}.log"
    conda:
        "../envs/tapirs.yaml"
    shell:
        "fastp \
        --disable_quality_filtering \
        --disable_adapter_trimming \
        --in1 {input.R1trimmed} \
        --in2 {input.R2trimmed} \
        --out1 {output.R1unmerged} \
        --out2 {output.R2unmerged} \
        --merge \
        --merged_out {output.merged} \
        --overlap_len_require {config[FASTP_min_overlap]} \
        --overlap_diff_limit {config[FASTP_diff_limit]} \
        --overlap_diff_percent_limit {config[FASTP_diff_percent_limit]} \
        --length_limit {config[FASTP_length_limit]} \
        --length_required {config[FASTP_len_required]} \
        -j {output.json} \
        -h {output.html} \
        --correction \
        > {log} 2>&1"


rule merge_forward_reads:
    input:
        merged     = "results/tapirs/03_merged/{LIBRARIES}/{SAMPLES}.merged.fastq",
        R1unpaired = "results/tapirs/02_trimmed/{LIBRARIES}/{SAMPLES}.R1.unpaired.fastq",
        R1unmerged = "results/tapirs/03_merged/{LIBRARIES}/{SAMPLES}.R1.unmerged.fastq"
    output:
        fq = "results/tapirs/04_forward_merged/{LIBRARIES}/{SAMPLES}.forward.merged.fastq"
    run:
        filenames = [input.merged, input.R1unpaired, input.R1unmerged]
        with open(str(output.fq), 'w') as outfile:
            for fname in filenames:
                with open(str(fname)) as infile:
                    outfile.write(infile.read())


rule seqkit_fq2fa:
    input:
        fq = "results/tapirs/04_forward_merged/{LIBRARIES}/{SAMPLES}.forward.merged.fastq"
    output:
        fa = "results/tapirs/05_fasta/{LIBRARIES}/{SAMPLES}.fasta"
    log:
        "logs/seqkit_fq2fa/{LIBRARIES}/{SAMPLES}.log"
    conda:
        "../envs/tapirs.yaml"
    shell:
        "seqkit fq2fa {input.fq} -o {output.fa} > {log} 2>&1"
