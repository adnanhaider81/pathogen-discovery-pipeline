# Pathogen Discovery Pipeline

A panel driven metagenomic confirmatory mapping workflow with host subtraction and competitive mapping.
Use it for bacteria, DNA viruses, fungi, and parasites by editing the taxon lists in `configs/`.
For RNA viruses, you need an RNA compatible library prep (or a separate RNA workflow) and an appropriate reference panel.

## What this pipeline is for

This workflow is designed for the step after an initial screen (for example Kraken2 Bracken, Kaiju, or Centrifuge).
It answers: do my reads preferentially map to a curated set of candidate organisms, and what is the coverage pattern?

It is not intended to replace full outbreak investigation, culture, targeted PCR, or epidemiology.
It produces evidence to prioritize confirmatory assays.

## Features

- Host subtraction (map to host genomes, keep read pairs where both mates are unmapped)
- Optional low complexity read filtering (BBduk entropy filter)
- Automatic download of reference genomes from NCBI using NCBI Datasets CLI
- Competitive mapping to a combined reference (reduces false positives vs mapping one genome at a time)
- Per sample, per reference summary metrics:
  - mapped primary reads
  - mapped primary reads with MAPQ >= MIN_MAPQ
  - length weighted mean depth
  - breadth percent coverage
- Simple, reproducible conda environment

## Repository layout

- `bin/pathogen_discovery_pipeline.sh` : main pipeline script
- `environment.yml` : conda environment
- `configs/hosts/host_taxa.txt` : host taxa to subtract
- `configs/panels/pathogen_taxa.txt` : bacteria fungi parasites panel
- `configs/panels/virus_taxa.txt` : optional virus panel
- `docs/interpretation.md` : how to interpret breadth depth and reduce false positives
- `docs/panel_design.md` : panel design guidance
- `docs/troubleshooting.md` : common issues and fixes

## Installation

Requirements: Linux, conda or mamba, and enough disk for reference genomes.

Create the environment:
```bash
conda env create -f environment.yml
conda activate pathogen_discovery_pipeline
```

Tip: `mamba env create -f environment.yml` is often faster.

## Configure panels

Edit these files:
- `configs/hosts/host_taxa.txt`
- `configs/panels/pathogen_taxa.txt`
- `configs/panels/virus_taxa.txt` (optional)

Each file is one taxon per line. Lines starting with `#` are comments.
Use NCBI recognized scientific names.

Recommendation: for every likely pathogen, also include several close relatives.
This helps prevent mis assignment when your sample contains a near neighbor not in the panel.

## Input FASTQs

The default mode assumes integer sample IDs and Illumina bcl2fastq style names:
- `1_S1_L001_R1_001.fastq.gz`
- `1_S1_L001_R2_001.fastq.gz`

If your naming differs, set templates using a Python format string with `{i}`:
```bash
R1_TEMPLATE="S{i}_R1.fastq.gz" R2_TEMPLATE="S{i}_R2.fastq.gz"
```

## Run

Basic run (sample IDs 1..22):
```bash
FASTQ_DIR=/path/to/fastqs \
  FIRST_SAMPLE=1 LAST_SAMPLE=22 \
  THREADS=14 MIN_MAPQ=20 \
  OUTDIR=results \
  bin/pathogen_discovery_pipeline.sh
```

Run a custom subset:
```bash
printf "3\n4\n8\n9\n14\n" > samples.txt
FASTQ_DIR=/path/to/fastqs SAMPLE_LIST_FILE=samples.txt OUTDIR=results bin/pathogen_discovery_pipeline.sh
```

Useful toggles:
- `MASK_REF=1` masks low complexity in the combined pathogen reference (dustmasker)
- `FILTER_LOW_COMPLEX_READS=1` filters low complexity reads after host subtraction (BBduk)
- `MIN_MAPQ=20` controls mapping quality cutoff used for coverage metrics
- `STRONG_MIN_MAPQ20`, `STRONG_MIN_BREADTH_PCT`, `STRONG_MIN_MEAN_DEPTH` control STRONG_SIGNAL labeling
- `WEAK_MIN_MAPQ20`, `WEAK_MIN_BREADTH_PCT` control WEAK_SIGNAL labeling

Example stricter thresholds:
```bash
STRONG_MIN_MAPQ20=1000 STRONG_MIN_BREADTH_PCT=0.2 STRONG_MIN_MEAN_DEPTH=0.01 \
  FASTQ_DIR=/path/to/fastqs OUTDIR=results bin/pathogen_discovery_pipeline.sh
```

## Outputs

Main summary:
- `OUTDIR/confirm_summary_competitive.tsv`

Intermediate outputs:
- `OUTDIR/bams/<sample>/host.unsorted.bam` : host alignment BAM
- `OUTDIR/nonhost_fastq/<sample>.nonhost_R1.fastq.gz` and `_R2.fastq.gz` : non host read pairs
- `OUTDIR/bams/<sample>/combined.bam` : competitive mapping BAM
- `OUTDIR/metrics/<sample>/coverage_by_contig.tsv` : `samtools coverage` per combined contig
- `OUTDIR/metrics/<sample>/summary_rows.tsv` : per reference rows appended to the main summary

## Interpretation quick guide

- Competitive mapping is stronger evidence than separate per genome mapping.
- Breadth of coverage across the genome is often more informative than mean depth at low abundance.
- Many reads with tiny breadth can indicate repeats, low complexity, or contamination.
- If a high consequence organism appears, do not report at species level unless close relatives are included and marker confirmation is done.

See `docs/interpretation.md` for detailed guidance and recommended follow ups.

## How to cite

If you use this pipeline in a report or publication, please cite:

Haider, S. A. (2025). Pathogen Discovery Pipeline (Version 0.1.0) [Software]. National Institute of Health, Islamabad. https://github.com/adnanhaider81/pathogen-discovery-pipeline

A machine readable citation file is provided as `CITATION.cff`.

## License

MIT license is included.

## Acknowledgements

This pipeline relies on open source bioinformatics tools including minimap2, samtools, BBMap suite, BLAST, and the NCBI Datasets CLI.
