# Troubleshooting

## NCBI Datasets download fails

- Check internet connectivity.
- Update the conda environment: `conda env update -f environment.yml`.
- Some taxa do not have a suitable assembly at the selected level. The pipeline will warn and skip.

## Pipeline runs but summary looks empty

Common causes:
- Reads are almost entirely host and non host FASTQs are tiny.
- MIN_MAPQ is too strict for your panel.
- The panel FASTAs failed to download, leaving an empty combined reference.
  Check `OUTDIR/pathogen_refs/pathogen_manifest.tsv`.

## Memory killed during indexing

If your combined reference becomes very large:
- Reduce the panel size.
- Prefer minimap2 mapping without heavy indexing options.
- Run on a machine with more RAM.

## Unexpected pathogen labels

- Add close relatives for that group.
- Enable `MASK_REF=1` and `FILTER_LOW_COMPLEX_READS=1`.
- Confirm with contigs and targeted PCR.

## RNA viruses are missing

If your library prep is DNA focused, RNA viruses will be absent or very low.
Use RNA extraction and RNA library prep for RNA virus discovery.
