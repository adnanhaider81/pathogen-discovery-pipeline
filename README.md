# Pathogen Discovery Pipeline

A reproducible metagenomic confirmatory mapping pipeline for pathogen discovery.

This repo is panel-driven: it can be used for bacteria, DNA viruses, fungi, and parasites by editing the taxon lists in `configs/`.
For RNA viruses, use an RNA library prep (or an RNA workflow) before confirmatory mapping.

## Workflow

1. Download host reference genomes (for host subtraction) from NCBI using the NCBI Datasets CLI
2. Map reads to host references and extract read pairs where both mates are unmapped (non-host)
3. Optional low-complexity read filtering (BBduk entropy filter)
4. Download curated pathogen reference genomes (plus optional virus genomes)
5. Build a combined competitive reference
6. Competitive mapping of non-host reads to the combined pathogen panel
7. Summarize breadth and depth per reference per sample

## Installation

```bash
conda env create -f environment.yml
conda activate pathogen_discovery_pipeline
```

## Configure panels

- `configs/hosts/host_taxa.txt`
- `configs/panels/pathogen_taxa.txt`
- `configs/panels/virus_taxa.txt` (optional)

## Run

Default FASTQ naming for integer sample IDs matches Illumina bcl2fastq:
`1_S1_L001_R1_001.fastq.gz` and `1_S1_L001_R2_001.fastq.gz`.

```bash
FASTQ_DIR=/path/to/fastqs \
FIRST_SAMPLE=1 LAST_SAMPLE=22 \
THREADS=14 MIN_MAPQ=20 \
OUTDIR=results \
bin/pathogen_discovery_pipeline.sh
```

Custom sample IDs:
```bash
printf "3\n4\n8\n9\n14\n" > samples.txt
FASTQ_DIR=/path/to/fastqs SAMPLE_LIST_FILE=samples.txt OUTDIR=results bin/pathogen_discovery_pipeline.sh
```

Custom FASTQ naming (Python format string using `{i}`):
```bash
R1_TEMPLATE="S{i}_R1.fastq.gz" R2_TEMPLATE="S{i}_R2.fastq.gz" \
FASTQ_DIR=/path/to/fastqs OUTDIR=results bin/pathogen_discovery_pipeline.sh
```

## Output

Main summary: `OUTDIR/confirm_summary_competitive.tsv`

See `docs/interpretation.md` and `docs/panel_design.md`.
