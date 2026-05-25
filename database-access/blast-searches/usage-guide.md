# BLAST Searches Usage Guide

## Overview

Run remote BLAST searches against NCBI's servers using `Bio.Blast.NCBIWWW`. Encodes the program-vs-database decision matrix, Karlin-Altschul E-value math and why bit-scores beat E-values across databases, the `max_target_seqs` misinterpretation that has biased many published workflows (Shah 2019), composition-based-statistics modes, word-size choice for cross-species vs intra-species DNA, and reproducibility traps with `nt`/`nr`.

## Prerequisites

```bash
pip install biopython
```

No API key required for remote BLAST; one search at a time is the polite cap.

## Quick Start

- "Identify this unknown DNA sequence -- BLASTN against refseq_select for stability, not nt"
- "Find mammalian homologs of this protein in Swiss-Prot using blastp with composition-based statistics mode 2"
- "Run blastp for this short peptide (12 aa) using PAM30 matrix and word size 2"
- "Use dc-megablast for this cross-species mRNA query -- megablast won't find diverged homologs"
- "Search a sequence with hitlist_size=500 and post-filter top 10 by bit-score, avoiding the max_target_seqs trap"

## Example Prompts

### Picking the right program

> "I have an unknown DNA sequence from an environmental sample. Run blastn against refseq_select_rna with word_size=11 and E-value cutoff 1e-10. Don't use megablast -- it requires 28-nt exact match seeds and will miss divergent homologs."

### Reproducibility-safe database choice

> "I'm BLASTing for a publication. Use refseq_select instead of nr/nt because those change daily and aren't reproducible. If I need the broader search, record today's date and archive a frozen copy of the database snapshot."

### Avoiding the max_target_seqs trap

> "Set hitlist_size=500 not 10. The Bio.Blast hitlist_size parameter maps to BLAST+'s max_target_seqs, which is an early-termination heuristic, not 'top N by E-value' (Shah et al. 2019). Then post-filter to the top 10 by bit-score in Python."

### Short-peptide search

> "BLAST this 12-aa proteomics peptide against swissprot. Use matrix='PAM30', word_size=2, expect=1000, composition_based_statistics=3. Default BLOSUM62/word=3 will miss everything for queries this short."

### Cross-database analysis caveat

> "I have results from a blastp vs nr and a blastp vs swissprot for the same query. Don't compare E-values across databases -- they scale with database size. Sort by bit-score for cross-DB comparison."

## What the Agent Will Do

1. Choose program from the (query molecule, target molecule) pair and the species-distance question.
2. Choose database from the reproducibility-vs-coverage tradeoff (refseq_select for stability; nt/nr only with snapshot date).
3. Set word_size, matrix, and composition-based-statistics appropriate for query length and protein composition.
4. Set hitlist_size large (500+) to dodge the max_target_seqs trap; post-filter top N.
5. Use entrez_query for pre-filtering by organism (faster and more meaningful E-values than post-filter).
6. Parse XML with NCBIXML.read(); compute identity and coverage from HSP attributes.
7. Sort by bit-score, not E-value, when comparing across databases.
8. Recommend defection to local-blast / DIAMOND / MMseqs2 for >50 sequences.

## Tips

- Bit-score is database-size normalized and is the correct cross-database metric. E-value is not.
- For cross-species DNA homology: use `dc-megablast` (discontiguous), never default megablast (word=28 misses divergent hits).
- For protein remote homology (E in 10^-3 to 10^-1 range), switch to PSI-BLAST / jackhmmer / Foldseek -- see `remote-homology`.
- `entrez_query='Mammalia[Organism]'` filters BLAST's search space before the search runs. Faster than post-filtering and gives correct (database-size-adjusted) E-values.
- For publication: never use `nt`/`nr` without recording snapshot date or archiving the database. `refseq_select` is the safe default.
- Remote BLAST submits to the NCBI queue. For >5 sequences/min, switch to local BLAST. For >1000 sequences, switch to DIAMOND or MMseqs2.
- The "twilight zone" of homology (Rost 1999 *Protein Eng* 12:85) is 20-35% identity; below 20% sequence-only search is unreliable -- need structure (Foldseek) or HMM.
- Compositional bias inflates E-values for low-complexity regions. CBS=2 (default) handles most; consider hard-masking with `filter='S'` for severe bias.

## Related Skills

- local-blast - Local BLAST+ for >50 sequences or custom databases
- remote-homology - PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek for distant homology
- ortholog-inference - RBH, OrthoFinder, OMA for orthology calls
- sequence-io/read-sequences - Load query FASTA files
- entrez-fetch - Fetch full GenBank records for top hits
