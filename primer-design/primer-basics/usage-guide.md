# PCR Primer Design - Usage Guide

## Overview

This skill designs and ranks PCR primer pairs for a target template with primer3-py. The central idea it enforces: primer3 is a LOCAL optimizer over the single template it is handed, so its top-ranked pair is the lowest-penalty candidate under the supplied bounds, not a guarantee that the primers amplify the target uniquely in the genome. It covers the seq_args/global_args constraint model, why Tm is a salt/concentration-dependent nearest-neighbor prediction rather than a fixed sequence property, 3'-end and GC-clamp mechanism, 5'-tail handling, masking SNPs under primer 3' ends, and diagnosing runs that return no pairs. Use it for standard PCR, cloning, genotyping, and sequencing primers; route genome-wide specificity to primer-specificity, dimer/hairpin checks to primer-validation, and qPCR assays to qpcr-primers.

## Prerequisites

```bash
pip install primer3-py
```

- Input is a target template sequence (a string or a FASTA loaded with biopython).
- Tm is computed under specific salt, Mg2+, dNTP, and oligo concentrations -- supply the real reaction values so the reported Tm is meaningful; the same sequence reports a different Tm under different conditions.
- primer3 does NOT check genome-wide specificity; a separate in-silico PCR / BLAST step is mandatory before ordering (primer-specificity).
- Coordinates are 0-based and every interval is [start, length], not [start, end].

## Quick Start

Tell your AI agent what you want to do:
- "Design PCR primers to amplify the central 600 bp of this template, Tm 58-62 C, GC 40-60%"
- "Find primers that flank exon 3 and force the amplicon to cover it"
- "Design primers but keep them off the SNP at position 220"
- "Design allele-specific (ARMS) primers to genotype this SNP"
- "What annealing temperature should I run these primers at?"
- "Why did primer3 return no pairs for this region?"

## Example Prompts

### Standard PCR design
> "I have a 2 kb template FASTA. Design primer pairs that amplify the central 600 bp region, Tm 58-62 C, GC 40-60%, and give me the top three ranked pairs with their Tm, product size, and pair penalty."

### Flank or target a feature
> "Design primers that flank positions 800-1100 so the amplicon fully covers that exon, and keep both primers within the clean region 600-1300."

### Avoid variants under the primers
> "Design primers for this region but exclude the common SNPs at positions 540 and 612 so I do not get allele dropout at the 3' end."

### Cloning with tails
> "Design the template-binding cores for Gibson assembly primers with 20 bp overlaps, and tell me to check the full tailed oligos for dimers separately."

### Allele-specific genotyping (ARMS)
> "Design an allele-specific PCR to genotype the SNP at position 300: one forward primer per allele with the discriminating base at the 3' end and a second destabilizing mismatch near it, sharing a common reverse primer."

### Annealing temperature
> "My primers have predicted Tm 60 and 61 C. What annealing temperature should I start at, and when should I use touchdown PCR instead of a fixed Ta?"

### Diagnose a failed run
> "primer3 returned zero pairs for my GC-rich template. Turn on the explain flag and tell me which single constraint to loosen first."

## What the Agent Will Do

1. Load the template sequence and confirm coordinates are 0-based with [start, length] intervals.
2. Place per-template data in seq_args (SEQUENCE_TEMPLATE plus any TARGET/INCLUDED/EXCLUDED/OVERLAP_JUNCTION constraint) and run-wide settings in global_args (Tm/GC/size bounds, salt, weights).
3. Call design_primers and read the ranked pairs from the flat result dict, reporting Tm, GC, product size, and pair penalty.
4. Flag a Tm mismatch between the two primers, an over-stable 3' end, or a SNP under a primer.
5. Diagnose a zero-pair run from the explain tallies and loosen one constraint at a time.
6. Hand the chosen pair to primer-validation (dimers/hairpins) and primer-specificity (genome off-target) before ordering.

## Tips

- Match the two primer Tms to within about 2 C (set PRIMER_PAIR_MAX_DIFF_TM); a mismatched pair lets one strand dominate.
- Tighten the default GC window (20-80) to 40-60%; extremes prime poorly and predict Tm badly.
- Set a GC clamp of 1, but cap PRIMER_MAX_END_STABILITY (library default 100 is effectively off) so the 3' end is not so stable it misprimes.
- Use a BOUND (MIN/MAX) to forbid something and a WEIGHT to merely discourage it; over-tightening bounds is what returns zero pairs.
- For a primer with a 5' tail, design the binding core in primer3 so the Tm reflects only the annealing region, then append the tail and re-check dimers on the full oligo.
- Never order the top pair without a genome specificity pass; primer3 has not looked outside the supplied template.
- Predicted Tm is not the annealing temperature; start Ta about 3-5 C below the lower primer Tm, optimize with a gradient, and use touchdown PCR when specificity is hard.
- For allele-specific (ARMS) genotyping, put the discriminating base at the 3' terminus and add a second mismatch at the -2 or -3 position to widen the allele discrimination.
- For a GC-rich or structured template that will not amplify, DMSO, betaine, or 7-deaza-dGTP lower the effective Tm and disrupt secondary structure -- a reagent lever to try before redesigning.

## Related Skills

- primer-validation - Check chosen primers for dimers, hairpins, and 3'-end stability
- primer-specificity - Confirm the pair amplifies only the target genome-wide
- qpcr-primers - Design qPCR primers and hydrolysis/molecular-beacon probes
- database-access/entrez-fetch - Fetch the template sequence to design against
- sequence-manipulation/seq-objects - Reverse-complement and extract subsequences
- sequence-io/read-sequences - Read the target FASTA
