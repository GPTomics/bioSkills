# Motif Search - Usage Guide

## Overview

This skill enables AI agents to help you find patterns, motifs, and subsequences in biological sequences using Biopython. It covers exact matches, overlapping-match counting, regex and IUPAC degenerate patterns, position weight matrix (PWM/PSSM) scoring with reproducible thresholds, and reading motif matrices from JASPAR, MEME, and TRANSFAC files.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Find all occurrences of GAATTC, including overlapping ones"
- "Find all matches to GATNNTC where N is any base"
- "Score this promoter for a JASPAR binding motif on both strands"
- "Pick a match threshold at a 1% false-positive rate"

## Example Prompts

### Exact Pattern Search
> "Find all positions where GAATTC occurs in my sequence, counting overlapping hits"

### Ambiguous / Degenerate Pattern
> "Find all matches to the pattern GATNNTC where N is any base"

### Regulatory Elements
> "Search for TATA box variants in the first 500 bp"

### Binding-Site Scoring
> "Build a PWM from these aligned binding sites and scan my sequence for matches above a 1% false-positive threshold"

### Reading Motif Matrices
> "Load this JASPAR PFM and report its consensus, then scan my promoter on both strands"

### Both Strands
> "Score the forward and reverse strands separately and report which strand each hit is on"

## What the Agent Will Do

1. Choose the search method from the pattern type (exact, degenerate IUPAC, regex, or PWM)
2. Use an overlap-aware method when overlapping hits matter (lookahead or nt_search)
3. For a PWM, set pseudocounts and background before reading the PSSM, then derive a threshold from the score distribution
4. Scan both strands and label each hit's strand
5. Return positions, matched subsequences, and scores

## Pattern Types

- **Exact**: Fixed sequence like `GAATTC`
- **IUPAC**: With ambiguity codes like `GATNNTC` (N = any base)
- **Regex**: Flexible patterns like `TATA[AT]A[AT]`
- **PWM/PSSM**: Probabilistic scoring matrices for graded binding-site matches

## Tips

- `str.count`/`re.findall` miss overlapping hits; use `re.finditer(r'(?=(motif))', seq)` or `Bio.SeqUtils.nt_search`
- `nt_search` expands IUPAC ambiguity in the query and reports overlaps; it returns `[pattern, pos1, pos2, ...]` (just `[pattern]` when there are no hits)
- Set a motif's pseudocounts (0.5 or sqrt(N)) before reading its PSSM, or count-0 cells produce -inf scores; the PSSM is recomputed on each access, so set pseudocounts and background first
- A PSSM score is log2-odds in bits, not a probability; derive a threshold from `pssm.distribution(...).threshold_fpr(fpr)` rather than picking one by eye
- Always search both strands; `pssm.search` defaults to scanning both, returning reverse-strand hits at negative positions
- Match the file format string exactly: CIS-BP/HOMER/HOCOMOCO use `pfm-four-columns`, not `cisbp`/`homer`/`hocomoco`; four-columns vs four-rows is the common mix-up

## Related Skills

- seq-objects - Create Seq objects for searching
- reverse-complement - Reverse-complement the target to search the opposite strand
- transcription-translation - ORF and codon-context motifs in coding sequences
- sequence-properties - GC content and per-sequence properties around hits
- restriction-analysis/restriction-sites - Restriction enzyme recognition sites
- chip-seq/motif-analysis - De novo motif discovery and enrichment in peak sets
- database-access/entrez-fetch - Download motif matrices from JASPAR/NCBI
