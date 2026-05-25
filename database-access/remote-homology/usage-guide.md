# Remote Homology Usage Guide

## Overview

Detect distant homologs using profile and structure-aware methods. Encodes the modern landscape: PSI-BLAST and jackhmmer for iterative profile search; HHblits/HHsearch for profile-profile (PDB70, Pfam); MMseqs2 and DIAMOND as the modern blastp replacements (100-10,000x faster); Foldseek (van Kempen 2024) for structure-aware homology via the 3Di alphabet, with ProstT5 for sequence-only access to structural search.

## Prerequisites

```bash
conda install -c bioconda hmmer mmseqs2 diamond hhsuite foldseek blast
hmmsearch -h | head -3       # HMMER 3.4+
mmseqs version               # MMseqs2 15+
diamond --version            # DIAMOND 2.1+
hhblits -h | head -3         # HH-suite3 3.3+
foldseek --version           # Foldseek 9+
```

For Foldseek databases:
```bash
mkdir -p foldseek_dbs tmp
foldseek databases --help    # list available DBs
```

## Quick Start

- "Find structural homologs of this protein via Foldseek against AlphaFoldDB Swiss-Prot subset"
- "Run jackhmmer 3 iterations against UniRef90; save the final HMM for downstream hmmsearch"
- "Search a protein query with MMseqs2 at -s 7.5 (HMMER-equivalent sensitivity) instead of BLAST"
- "Annotate domains with hmmscan against Pfam-A using --cut_ga calibrated thresholds"
- "DIAMOND --ultra-sensitive against UniRef90 for a metagenomic protein set"

## Example Prompts

### Foldseek for divergent structural homologs

> "I have a protein where blastp finds nothing significant. Predict the structure (or use ProstT5 for sequence-only) and search AlphaFoldDB with Foldseek easy-search. Report TM-scores and probability cutoff -- hits with prob > 0.9 are structurally confident."

### PSI-BLAST with drift protection

> "Run psiblast for 3 iterations against UniRef90 with -inclusion_ethresh 0.002 (stricter than the 0.005 default). Save the PSSM with -out_pssm and the included sequence set so I can audit for paralog contamination. Don't iterate to convergence -- 4+ iterations drift."

### MMseqs2 as PSI-BLAST replacement

> "Same iterative profile search as PSI-BLAST but use MMseqs2 with --num-iterations 3 and -s 7.5 (HMMER-equivalent sensitivity). It's 100x faster."

### Pfam domain annotation

> "Annotate domains of these 5,000 protein sequences using hmmscan against Pfam-A with --cut_ga. Use the calibrated gathering thresholds, not arbitrary E-value cutoffs."

### HHsearch against PDB70

> "Build an HHblits profile of this protein against UniRef30 (3 iterations), then hhsearch against PDB70 for the deepest possible structural homology to known PDB entries."

### DIAMOND for metagenomic scale

> "I have 1 million predicted ORFs from a metagenome. Use DIAMOND blastp --ultra-sensitive against UniRef90 with -p 32. blastp would take days; DIAMOND will finish in an hour."

## What the Agent Will Do

1. Identify whether the problem is single-query distant homology (PSI-BLAST, jackhmmer, Foldseek) or large-batch (MMseqs2, DIAMOND).
2. Default to MMseqs2 (-s 7.5) or DIAMOND --ultra-sensitive over BLAST for any batch >50 sequences.
3. For very deep homology to known structures, run Foldseek (with ProstT5 if no structure available) against AlphaFoldDB or PDB100.
4. For domain annotation, use hmmscan vs Pfam-A with --cut_ga.
5. For PSI-BLAST, cap iterations at 3 and use stricter inclusion threshold; save PSSM for reproducibility.
6. For HHsearch profile-profile, build query MSA via HHblits first.
7. Recommend mixing sequence + structure evidence for the toughest cases.
8. Note when fold-level similarity alone (Foldseek to a TIM barrel) is not evidence of homology.

## Tips

- The twilight zone of sequence homology is 20-35% pairwise identity (Rost 1999). Below that, structure beats sequence -- reach for Foldseek.
- ProstT5 (Heinzinger 2024) is the bridge: sequence -> 3Di alphabet directly, no AF2 step needed for Foldseek.
- For metagenomic-scale work, MMseqs2 `easy-cluster` is the modern de-facto for clustering hundreds of millions of sequences.
- PSI-BLAST is **non-deterministic** in detail (input order matters). For reproducibility, save the PSSM and re-use with `-in_pssm`.
- HHsearch vs PDB70 is still the gold standard for deepest possible homology to PDB; for AlphaFoldDB coverage, prefer Foldseek.
- DIAMOND v2 adds frameshift-aware mode (`--frameshift 15`) for nanopore / PacBio long-read protein search.
- Common folds (TIM barrel, Rossmann, ABC) appear across many superfamilies. A high Foldseek score on a common fold is necessary but not sufficient for homology.
- The Foldseek-multimer mode (Foldseek 9+) searches protein complexes against multimer databases -- relevant for interactome work.

## Related Skills

- blast-searches - Remote BLAST baseline
- local-blast - Local BLAST+ for moderate-scale work
- ortholog-inference - Distinct topic: orthology vs homology (RBH, OrthoFinder, OMA)
- alignment/multiple-alignment - Build MSAs for HMM profiles
- structural-biology/alphafold-predictions - Predict structures for Foldseek queries
