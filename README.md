# 12S-edna-dada2-tapirs-workflow

Merged Snakemake workflow that runs the **Tapirs** metabarcoding pipeline and the **DADA2/SINTAX** pipeline simultaneously from a single folder of paired-end 12S fastq.gz files.

## Outputs

| Directory | Pipeline | Contents |
|-----------|----------|----------|
| `results/tapirs/` | Tapirs | OTU tables (BLAST and/or Kraken2), all intermediates |
| `results/dada2/`  | DADA2  | ASV FASTA, abundance tables, SINTAX taxonomy, blank-cleaned matrices |

## Quick start

### 0. Install Snakemake

The workflow requires Snakemake 7.32.4. Install it into a dedicated conda environment:

```bash
conda create -p /hdd0/susbus/tools/conda_envs/snakemake -c conda-forge -c bioconda \
  snakemake=7.32.4 -y
conda activate /hdd0/susbus/tools/conda_envs/snakemake
```

### 1. Edit the config

```yaml
# config/config.yaml
input_dir: "/path/to/your/fastq_folder"   # flat folder with *_R1_001.fastq.gz pairs
my_experiment: "my_run"
analysis_method: "both"                   # "blast", "kraken2", or "both"
```

Update the four database paths (blast_db, kraken2_db, taxdump, dada2_ref_db).

### 2. Run

```bash
cd /prj/DECODE/12S_pipeline/12S-edna-dada2-tapirs-workflow

# Dry run
snakemake -s workflow/Snakefile --use-conda --cores 64 \
  --conda-prefix /prj/DECODE/conda_envs \
  --conda-frontend conda \
  -rpk -n

# Full run
snakemake -s workflow/Snakefile --use-conda --cores 64 \
  --conda-prefix /prj/DECODE/conda_envs \
  --conda-frontend conda \
  -rpk

# HPC (SLURM example)
snakemake -s workflow/Snakefile --use-conda --cores 64 \
  --conda-prefix /prj/DECODE/conda_envs \
  --conda-frontend conda \
  --cluster "sbatch -c {threads} --mem=32G -t 12:00:00" \
  --jobs 50 -rpk
```

> **Note:** `--conda-frontend conda` is required to avoid a mamba/conda version conflict.
> Conda environments are built on first run into `/prj/DECODE/conda_envs/` and reused thereafter.

Conda environments are created automatically on first run inside `.snakemake/conda/`.

## Input file naming

Files must follow the standard Illumina naming convention:

```
{SAMPLE}_R1_001.fastq.gz
{SAMPLE}_R2_001.fastq.gz
```

Blank controls are auto-detected by name pattern:

| Type | Pattern (case-insensitive) |
|------|---------------------------|
| PCR blank | `PCR_BLANK`, `PCRBLANK`, `NEG`, `POS` |
| Extraction blank | `EXTRACTION_BLANK`, `EXTRACTIONBLANK`, `_EB_`, `_EB[0-9]` |
| Site blank | any remaining name containing `BLANK` |

## Tapirs pipeline (results/tapirs/)

Steps: fastp trim → fastp merge → vsearch dereplicate → vsearch UNOISE3 denoise → chimera detection → BLAST/Kraken2 → MLCA taxonomy → OTU tables

Key output files:
- `results/tapirs/{experiment}_blast{id}_{method}.tsv` — OTU abundance table (BLAST)
- `results/tapirs/{experiment}_blast{id}_{method}_full_lineage.tsv` — full lineage variant
- `results/tapirs/{experiment}_kraken2_conf{conf}_{method}.tsv` — OTU abundance table (Kraken2)
- `results/tapirs/{experiment}_kraken2_conf{conf}_{method}_full_lineage.tsv`

## DADA2 pipeline (results/dada2/)

Steps: cutadapt primer trim → DADA2 filter/denoise/merge → chimera removal → vsearch SINTAX (INBO, MIDORI2, CLARE databases) → taxon-level collapse → LOD-based blank cleanup

Key output files (file names unchanged from original pipeline):

| File | Description |
|------|-------------|
| `asvs.nochim.fasta` | Final ASV sequences |
| `asv_lookup.tsv` | ASV ID ↔ sequence mapping |
| `seqtab_asv.csv` | ASV × sample abundance matrix |
| `asv_taxonomy_abundance_CLARE.csv` | ASV + CLARE taxonomy + abundance |
| `asv_taxonomy_compare.tsv` | Three-database comparison |
| `ncl_cleaned_labLOD.csv` | Lab-blank LOD cleaned matrix |
| `ncl_cleaned_siteLOD.csv` | Site-blank LOD cleaned matrix |
| `ncl_cleaned_bothLOD.csv` | Combined LOD cleaned matrix (**main dataset**) |
| `ncl_cleaned_long.csv` | Long-format with all cleaned versions |

## Tools

| Tool | Version | Purpose |
|------|---------|---------|
| fastp | ≥0.23 | Tapirs QC and merging |
| vsearch | ≥2.22 | Tapirs clustering; SINTAX for DADA2 |
| BLAST | ≥2.9 | Tapirs sequence search |
| kraken2 | ≥2.1 | Tapirs k-mer classification |
| seqkit | any | FASTQ→FASTA conversion |
| dada2 | ≥1.26 | ASV inference |
| cutadapt | ≥4 | DADA2 primer trimming |

`usearch` is **not required** — vsearch provides a free, conda-installable SINTAX implementation with compatible output format.

## Notes

- `results/tapirs/` and `results/dada2/` are written simultaneously; both sub-workflows are fully independent and run in parallel under Snakemake.
- `config/samples.tsv` is auto-generated from the input directory; do not edit manually.
- Human and domestic livestock taxa are removed during DADA2 blank cleanup; fish are retained.
