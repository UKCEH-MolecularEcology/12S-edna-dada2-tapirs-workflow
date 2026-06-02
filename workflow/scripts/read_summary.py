# ------------------------------------------------------------------------------
# read_summary.py
# Build a per-sample read count summary across both pipelines.
#
# Columns:
#   sample                  sample name
#   raw_reads               raw read pairs (from fastp trim JSON)
#   tapirs_after_trim       read pairs surviving fastp quality trim
#   tapirs_rereplicated     reads in rereplicated FASTA (= reads entering BLAST/Kraken2)
#   dada2_after_cutadapt    read pairs surviving primer trimming (from filtering_summary.csv)
#   dada2_after_filter      read pairs surviving DADA2 quality filter
#   dada2_final             total reads assigned to ASVs in sequence table
# ------------------------------------------------------------------------------

import json
import os
import re
import pandas as pd

trim_jsons   = sorted(snakemake.input.fastp_trim_jsons)
rerep_fastas = sorted(snakemake.input.rerep_fastas)
filter_csv   = snakemake.input.dada2_filtering
seqtab_csv   = snakemake.input.dada2_seqtab
outfile      = snakemake.output.summary


def parse_trim_json(path):
    with open(path) as f:
        d = json.load(f)
    # fastp counts R1 + R2 separately; divide by 2 for read pairs
    raw     = d['summary']['before_filtering']['total_reads'] // 2
    trimmed = d['summary']['after_filtering']['total_reads']  // 2
    return raw, trimmed


def count_fasta(path):
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        return sum(1 for line in f if line.startswith('>'))


def sample_from_path(path, suffix):
    base = os.path.basename(path)
    return base[:-len(suffix)] if base.endswith(suffix) else os.path.splitext(base)[0]


records = {}

# ── Raw reads + Tapirs trimmed ─────────────────────────────────────────────
for path in trim_jsons:
    sample = sample_from_path(path, '.fastp.json')
    raw, trimmed = parse_trim_json(path)
    records[sample] = {
        'sample':             sample,
        'raw_reads':          raw,
        'tapirs_after_trim':  trimmed,
    }

# ── Tapirs rereplicated reads ──────────────────────────────────────────────
for path in rerep_fastas:
    sample = sample_from_path(path, '.rerep.fasta')
    records.setdefault(sample, {'sample': sample})
    records[sample]['tapirs_rereplicated'] = count_fasta(path)

# ── DADA2 filtering summary ────────────────────────────────────────────────
# filterAndTrim rownames are the input (cutadapt) file paths
filt = pd.read_csv(filter_csv, index_col=0)
for idx, row in filt.iterrows():
    sample = re.sub(r'_R1_001.*$', '', os.path.basename(str(idx)))
    records.setdefault(sample, {'sample': sample})
    records[sample]['dada2_after_cutadapt'] = int(row['reads.in'])
    records[sample]['dada2_after_filter']   = int(row['reads.out'])

# ── DADA2 final read counts from seqtab ───────────────────────────────────
seqtab   = pd.read_csv(seqtab_csv)
asv_cols = [c for c in seqtab.columns if c != 'sample']
for _, row in seqtab.iterrows():
    sample = str(row['sample'])
    records.setdefault(sample, {'sample': sample})
    records[sample]['dada2_final'] = int(row[asv_cols].sum())

# ── Write output ───────────────────────────────────────────────────────────
cols = [
    'sample',
    'raw_reads',
    'tapirs_after_trim',
    'tapirs_rereplicated',
    'dada2_after_cutadapt',
    'dada2_after_filter',
    'dada2_final',
]

df = pd.DataFrame(list(records.values()))
for c in cols:
    if c not in df.columns:
        df[c] = pd.NA

df = df[cols].sort_values('sample').reset_index(drop=True)
df.to_csv(outfile, sep='\t', index=False)
