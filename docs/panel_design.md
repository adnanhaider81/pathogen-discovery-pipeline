# Panel design

Your panel is a hypothesis list. Start broad, then tighten based on results.

## Principles

1. Include plausible pathogens for the syndrome, host species, geography, and season.
2. Include close relatives for each plausible organism.
   This is critical to reduce mis assignment.
3. Keep the panel small enough to run efficiently, but large enough to be informative.
4. Use NCBI recognized taxon names. The pipeline uses the NCBI Datasets CLI.

## Practical workflow

1. Build a broad panel for the syndrome.
2. Run on a few representative samples.
3. Inspect which groups appear as candidates.
4. Add more close relatives to the candidate groups.
5. Rerun confirmatory mapping.
6. Validate top candidates with targeted assays.

## Tips by organism type

Bacteria:
- Add multiple species within the same genus if one appears.
- For high consequence species, include near neighbors and report conservatively.

DNA viruses:
- Include multiple strains or representative genomes when available.
- For poxviruses and herpesviruses, include several species in the same genus.

Fungi:
- Fungal genomes can be large and repetitive.
  Expect lower specificity unless you include close relatives and use masking.

Parasites:
- Many protozoan genomes contain repeats. Masking and filtering are especially helpful.
