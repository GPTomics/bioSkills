# primer-design

## Overview

Design, validate, and specificity-check PCR and qPCR primers with primer3-py. The category enforces three traps practitioners miss: Tm is a salt/concentration-dependent nearest-neighbor prediction, not a fixed sequence property; primer3 scores thermodynamics on the supplied template only and does NOT check genome specificity (a separate mandatory step); and a qPCR assay is a quantitative measurement device whose validity rests on amplification efficiency, not just amplification.

**Tool type:** python | **Primary tools:** primer3-py, MFEprimer-3.0, Primer-BLAST, UCSC isPcr

## Skills

| Skill | Description |
|-------|-------------|
| primer-basics | Design and rank PCR primer pairs under Tm/GC/size/complementarity constraints on a template |
| primer-validation | Validate chosen oligos for hairpins, self-/cross-dimers, and 3'-end stability at reaction conditions |
| primer-specificity | Confirm a primer pair amplifies only the target genome-wide via pair-aware in-silico PCR |
| qpcr-primers | Co-design qPCR primers and TaqMan/molecular-beacon probes for a valid quantitative assay |

## Method Selection

| Scenario | Method | Skill |
|----------|--------|-------|
| Standard PCR / cloning / sequencing amplicon | primer3 design_primers under Tm/GC/size constraints | primer-basics |
| Check chosen primers for dimers/hairpins | calc_hairpin / calc_homodimer / calc_heterodimer / calc_end_stability | primer-validation |
| Troubleshoot a primer-dimer band or smear | 3'-end dimer check at the real reaction conditions | primer-validation |
| Confirm primers amplify only the target genome-wide | pair-aware in-silico PCR (MFEprimer / Primer-BLAST / isPcr) | primer-specificity |
| Real-time / quantitative assay (TaqMan or SYBR) | co-designed primers + probe, short amplicon, efficiency | qpcr-primers |
| Expression assay, avoid genomic DNA | exon-junction-spanning primers + genome pseudogene check | qpcr-primers + primer-specificity |

## Example Prompts

- "Design PCR primers to amplify the central region of this template, Tm 58-62 C"
- "Validate this primer pair for dimers at my reaction conditions"
- "Check whether these primers amplify anything besides my target in the human genome"
- "Design a TaqMan qPCR assay with a probe Tm 8-10 C above the primers"
- "Make exon-spanning primers that will not amplify genomic DNA"

## Requirements

```bash
pip install primer3-py
# for genome specificity (primer-specificity): MFEprimer / BLAST / UCSC isPcr, or NCBI Primer-BLAST (web)
conda install -c bioconda mfeprimer blast
```

## Related Skills

- **sequence-manipulation** - Reverse-complement, subsequence extraction, reading frames
- **database-access** - Fetch template sequences and build/query BLAST databases
- **read-alignment** - Align candidate amplicons / reads to a genome
- **differential-expression** - Downstream analysis that qPCR assays validate against
