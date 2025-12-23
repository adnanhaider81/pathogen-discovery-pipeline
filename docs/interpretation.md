# Interpretation

Competitive mapping improves specificity compared to mapping each pathogen separately, but confirmation still needs care.

## Practical signals

- Breadth of coverage across a genome is often more informative than mean depth at low abundance.
- A very small breadth with many reads can come from repeats, low complexity, or contamination.

## Increase confidence

1. Include close relatives in the panel
2. Mask low complexity in the reference (dustmasker)
3. Filter low-complexity reads (BBduk entropy filter)
4. Confirm with contigs when possible
5. Use negative controls (extraction and library blanks)
