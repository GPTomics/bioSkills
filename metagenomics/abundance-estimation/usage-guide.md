# Abundance Estimation - Usage Guide

## Overview
Bracken re-estimates species abundance from a Kraken2 report by redistributing reads stranded at higher taxonomic ranks back down to species. But Bracken is step one: the resulting table is a composition (the sequencer fixed the total), so making a defensible statement requires choosing an estimand, handling the compositional structure (CLR, zero replacement), deciding normalization, using a reference-frame test, and - for any "increased/decreased" claim - anchoring to absolute load. Bracken read fractions and MetaPhlAn percentages are different physical quantities and must not be merged.

## Prerequisites
```bash
conda install -c bioconda bracken kraken2
pip install scikit-bio pandas
# R compositional tooling (optional, for sparse-table zero replacement and DA):
# install.packages('zCompositions'); BiocManager::install(c('ALDEx2','ANCOMBC'))
```

Conceptual prerequisites:
- A relative-abundance table is compositional. "Taxon X increased" is undefined without a reference frame or an external load anchor.
- Bracken `-r` must match the bracken-build `-l` and the actual post-trim read length; it is not auto-detected.
- No library-size normalization fixes the genome-size confound - only a coverage estimand or a genome-normalized profiler does.
- Decide rarefaction per analysis: defensible for diversity, contested for differential abundance.

## Quick Start
Tell your AI agent what you want to do:
- "Run Bracken on my Kraken2 reports at read length 150 and build a species abundance matrix"
- "CLR-transform my abundance table and handle the zeros properly before stats"
- "Help me decide whether to rarefy for this analysis"
- "I want to claim a taxon increased - what do I need beyond relative abundances?"

## Example Prompts

### Bracken to matrix
> "I have Kraken2 reports for 24 gut samples at 150 bp. Run Bracken at species level, combine into one abundance matrix, and flag any species whose abundance came mostly from redistributed rather than directly assigned reads."

### Compositional handling
> "Transform my species-by-sample count matrix with CLR after Bayesian-multiplicative zero replacement, and compute an Aitchison distance matrix for ordination."

### Relative vs absolute
> "My relative abundances say Bacteroides dropped after treatment. We also have flow-cytometry cell counts. Convert to absolute load and tell me whether Bacteroides actually decreased or just lost share to a bloom."

### Cross-tool comparison
> "I have both Bracken and MetaPhlAn profiles for the same samples. Should I merge them into one table?"

## What the Agent Will Do
1. Run Bracken with the correct read length and level, then combine per-sample outputs into a matrix.
2. Flag species dominated by redistributed reads (likely fabricated relatives of database-absent taxa).
3. Replace zeros (multiplicative/Bayesian) and CLR-transform before any multivariate or correlation analysis.
4. Recommend a normalization, and decide rarefaction per the downstream analysis.
5. Use a reference-frame DA method (ALDEx2/ANCOM-BC) rather than naive tests on proportions.
6. Convert to absolute load only when an external anchor exists, and state the estimand throughout.

## Tips
- Distrust a species with large `added_reads` but tiny `kraken_assigned_reads` - it may be a redistribution artifact.
- `fraction_total_reads` is a fraction of classified-and-retained reads; differing host fractions make denominators non-comparable.
- Use CLR + Aitchison distance, never Pearson or Bray-Curtis on raw proportions for co-occurrence.
- GMPR or CSS survive sparse tables that break DESeq/edgeR median-of-ratios.
- Keep Bracken (read fraction) and MetaPhlAn (cell fraction) in separate tables.

## Output Columns

| Column | Description |
|--------|-------------|
| kraken_assigned_reads | Direct Kraken2 assignments |
| added_reads | Reads redistributed from higher levels |
| new_est_reads | Total estimated reads (assigned + added) |
| fraction_total_reads | Proportion of classified-and-retained reads |

## Related Skills

- kraken-classification - Generates the Kraken2 report; owns the read-count-is-not-abundance reframe
- metaphlan-profiling - Cell-fraction abundance; a different estimand
- metagenome-visualization - Diversity, ordination, and DA tool mechanics
- contamination-controls - Host-fraction handling that affects the denominator
- workflows/metagenomics-pipeline - End-to-end profiling and abundance

## Resources
- [Bracken GitHub](https://github.com/jenniferlu717/Bracken)
- [Bracken Paper](https://peerj.com/articles/cs-104/)
- [scikit-bio composition module](https://scikit.bio)
