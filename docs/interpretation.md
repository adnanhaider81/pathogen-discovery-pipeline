# Interpretation guide

This workflow performs confirmatory competitive mapping. It is designed to reduce false positives relative to mapping each candidate genome independently.
Even so, metagenomic confirmation needs careful interpretation.

## What the summary columns mean

In `confirm_summary_competitive.tsv` each row is a sample by reference.

- `mapped_primary`
  Number of primary alignments to that reference (MAPQ not applied).
- `mapped_primary_mapq20` (or MAPQ cutoff you set)
  Number of primary alignments with MAPQ >= MIN_MAPQ.
  MAPQ is a proxy for uniqueness and helps suppress ambiguous placements.
- `mean_depth`
  Length weighted mean depth across the reference, aggregated from `samtools coverage` per contig.
  In low abundance metagenomes, mean depth is usually very small.
- `breadth_covered_pct`
  Percent of the reference covered by at least 1 read at MAPQ >= MIN_MAPQ.
  Breadth is often the most useful metric for distinguishing real signals from repeat driven artifacts.
- `signal_label`
  A convenience label computed from thresholds you can tune.

## Interpreting breadth and depth

Common patterns:

1. Real signal (good)
   - Breadth is measurable and distributed across many loci
   - Reads align with reasonable MAPQ
   - Breadth increases with sequencing depth or with additional samples

2. Repeat or low complexity driven signal (often false positive)
   - Many reads but breadth remains extremely small
   - Reads cluster in a few regions or contigs
   - Signal changes substantially when you mask references or filter low complexity reads

3. Cross mapping due to missing close relatives
   - A species in the panel absorbs reads from a near neighbor not present
   - Adding close relatives redistributes the reads and reduces incorrect species calls
   - In this situation, report at genus or complex level unless marker confirmation is done

## How to raise confidence

1. Add close relatives
   For each likely organism, include multiple close relatives. Examples:
   - Pasteurellaceae: Mannheimia, Pasteurella, Bibersteinia, Histophilus
   - Bacillus: multiple Bacillus cereus group genomes rather than a single Bacillus anthracis
   - Poxviruses: multiple orthopox viruses if one is suspected

2. Use masking and low complexity filtering
   - `MASK_REF=1` runs dustmasker on the combined pathogen reference
   - `FILTER_LOW_COMPLEX_READS=1` runs BBduk entropy filtering after host subtraction
   If a call disappears after these steps, treat it as low confidence.

3. Use contigs for confirmation
   If you assemble contigs, map contigs to the combined pathogen reference.
   Long high identity alignments spanning non repetitive regions are strong evidence.
   Contig based confirmation is often more convincing than short read mapping alone.

4. Use negative controls
   Include extraction blanks and library blanks when possible.
   Low level signals that also appear in blanks are often contamination.

5. Validate with targeted assays
   For top candidates, use organism specific PCR, qPCR, and or culture.
   For bacteria, consider 16S or marker gene targets.
   For viruses, consider amplicon based confirmation.

## Host subtraction and limitations

- Host subtraction can remove the majority of reads in respiratory swabs.
  This is expected and not a failure.
- Very high host fraction can reduce sensitivity.
  Re extraction, deeper sequencing, or improved sampling can help.
- The host panel should match your sample source. Add relevant hosts if needed.
- This pipeline maps and subtracts host first, which often speeds up and improves specificity relative to using host decoys only.

## Reporting guidance

- Prefer genus level reporting unless close relatives are included and markers confirm species.
- Avoid high consequence claims from low breadth evidence.
- Provide the exact panel used, the MAPQ cutoff, and the breadth and depth values.
- State that results are metagenomic evidence requiring orthogonal confirmation.
