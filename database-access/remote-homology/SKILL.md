---
name: bio-remote-homology
description: Detect distant homologs using profile and structure-aware methods that go beyond standard BLAST. Use when sequence identity falls into the twilight zone (<35% pairwise), when BLAST fails to find homologs that should exist, when working at metagenomic scale (DIAMOND, MMseqs2), or when structure beats sequence (Foldseek). Covers PSI-BLAST (iterative PSSM), jackhmmer (iterative HMM), HHblits/HHsearch (profile-profile), DIAMOND, MMseqs2, and Foldseek (3Di structural alphabet, van Kempen 2024).
tool_type: mixed
primary_tool: HMMER
---

## Version Compatibility

Reference examples tested with: NCBI BLAST+ 2.15+, HMMER 3.4+, MMseqs2 15+, DIAMOND 2.1+, HH-suite3 3.3+, Foldseek 9+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then introspect signatures

If a flag is unrecognized or behavior changes, introspect with `--help` and adapt the example to match the installed version rather than retrying.

# Remote Homology

**"Find homologs my BLAST missed"** -> Standard BLAST detects similarity reliably down to ~35% pairwise identity (the "twilight zone", Rost 1999 *Protein Eng* 12:85). Below that, profile methods (PSSMs, HMMs) and structure-aware methods (Foldseek) recover homologs that pairwise alignment misses.

This skill covers the decision: which method, when, against what database. The competition has shifted substantially since 2015: PSI-BLAST is no longer the de-facto standard; MMseqs2 and DIAMOND have replaced BLAST in most large-scale workflows; Foldseek (van Kempen et al. 2024 *Nat Biotechnol* 42:243) detects homologs no sequence method can reach by searching with a 3Di structural alphabet derived from AlphaFold/ESMFold predictions.

- CLI: `psiblast`, `jackhmmer`, `hmmsearch`, `hhblits`, `mmseqs`, `diamond`, `foldseek`
- Python: `Bio.SearchIO` for output parsing; tool-specific clients exist but subprocess is preferred
- Web: HHpred (HHblits webserver), Foldseek webserver, ColabFold for paired structure search

## Required Setup

```bash
# Install via conda
conda install -c bioconda hmmer mmseqs2 diamond hhsuite foldseek
# BLAST+ (separate)
conda install -c bioconda blast

# Verify
hmmsearch -h | head -3       # HMMER 3.4+
mmseqs version               # MMseqs2 15+
diamond --version            # DIAMOND 2.1+
hhblits -h | head -3         # HH-suite3 3.3+
foldseek --version           # Foldseek 9+
```

## Decision matrix: which method when

| Question | Best tool | Why | Sensitivity / Speed |
|---|---|---|---|
| Quick all-vs-all proteome | MMseqs2 or DIAMOND | 100-10,000x faster than BLAST at comparable sensitivity | Highest throughput, near-BLAST sensitivity |
| Identify distant protein homolog (single query) | jackhmmer | Iterative HMM; usually beats PSI-BLAST | Higher sensitivity than PSI-BLAST |
| Distant homology where structure available | Foldseek | 3Di alphabet finds homologs sequence misses | Finds hits PSI-BLAST/HMMER cannot |
| Profile-profile comparison (PDB70 / Pfam) | HHblits + HHsearch | Profile vs profile is most sensitive when target also has profile | Best sensitivity for very-deep homology |
| Domain assignment | hmmscan against Pfam-A | Curated, calibrated thresholds | Standard practice |
| Metagenomic protein clustering | MMseqs2 `easy-cluster` | Scales to >1B sequences | Production-grade |
| ORF search vs metagenome | DIAMOND `blastx --frameshift` | Frameshift-aware; long reads | Best for noisy long reads |
| Structure-aware homology (no AF2 prediction available) | Foldseek + ProstT5 | Predicts 3Di alphabet from sequence via PLM | Skip the AF2 step |

## Foldseek: the 2024 revolution

Foldseek (van Kempen, Kim, Tumescheit et al. 2024 *Nat Biotechnol* 42:243) searches protein structures by representing each residue's local geometry as a 21-letter "3Di" alphabet, then running BLAST-style alignment in this alphabet. Two consequences:
- **4-5 orders of magnitude faster than DALI** (the previous gold-standard structure aligner).
- **Finds homologs that sequence methods cannot**: when sequence has diverged past detection but structure is preserved, Foldseek recovers the homology. Reported sensitivity vs traditional structure search is comparable; sensitivity vs sequence methods is dramatically higher at low identity.

Two access modes:
1. **Have a structure** (PDB or AlphaFold): `foldseek easy-search query.pdb db_dir result.m8 tmp_dir`
2. **Sequence only, no structure**: use **ProstT5** (Heinzinger et al. 2024) to embed sequence to 3Di alphabet directly, skipping AF2 entirely: `foldseek databases ProstT5 prostt5_db tmp` then `foldseek easy-search seq.fa db result.m8 tmp --prostt5-model prostt5_db`

The major prebuilt Foldseek databases (AlphaFoldDB, PDB100, ESMAtlas) are downloadable via `foldseek databases`.

## PSI-BLAST: still useful, but watch the drift

PSI-BLAST (Altschul et al. 1997 *Nucleic Acids Res* 25:3389) builds a PSSM iteratively: each iteration includes hits below `-inclusion_ethresh` (default 0.005) in the next PSSM. Convergence is when no new hits cross the threshold. **Stopping at convergence is often the wrong call** -- iterations 2-3 are usually optimal; iterations 4+ frequently drift into paralog inclusion, contaminating the PSSM.

| Parameter | Default | Postdoc tuning |
|---|---|---|
| `-num_iterations` | 1 | 3 for most workflows; >3 risks drift |
| `-inclusion_ethresh` | 0.005 | 0.002 if specificity matters (Altschul 1997 recommendation) |
| `-evalue` | 10 | 0.01 for reporting cutoff |
| `-num_threads` | 1 | 8 for large DBs |

PSI-BLAST is also **non-deterministic** in detail: different input order or DB version can produce different PSSMs. For reproducibility, save the PSSM (`-out_pssm pssm.asn`) and re-use with `-in_pssm`.

## HMMER 3 (hmmsearch, jackhmmer)

HMMER 3 (Eddy 2011 *PLoS Comput Biol* 7:e1002195) is profile HMM search. Two main workflows:
- **`hmmsearch profile.hmm seqdb`**: search a database with a known HMM (Pfam, custom).
- **`jackhmmer query.fa seqdb`**: iterative search like PSI-BLAST but with full HMM math. Typically higher sensitivity than PSI-BLAST at the same number of iterations.

For domain assignment, `hmmscan query.fa Pfam-A.hmm` is the canonical pipeline. Pfam HMMs come with calibrated gathering thresholds (`-gathering`) -- use them instead of arbitrary E-value cutoffs.

```bash
# Build domain database once
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

# Annotate query against Pfam-A with gathering threshold (calibrated cutoff)
hmmscan --cut_ga --domtblout query.domtbl Pfam-A.hmm query.fa
```

## HHblits / HHsearch (HH-suite3)

HHblits (Remmert et al. 2012 *Nat Methods* 9:173; HH-suite3: Steinegger et al. 2019 *BMC Bioinformatics* 20:473) is profile-profile alignment. The most sensitive method when both query and target have an HMM representation. Standard pipeline:
1. Build query MSA with `hhblits` against UniRef30 (or HHblits' default DB).
2. Convert MSA to query HMM.
3. Search against a profile DB (PDB70 for structure, Pfam-A for domains) with `hhsearch`.

Output is in HHM format. For very deep homology (the structural twilight zone), HHsearch vs PDB70 is still the gold standard.

## MMseqs2 (the modern protein search workhorse)

MMseqs2 (Steinegger & Soding 2017 *Nat Biotechnol* 35:1026) replaces BLAST in nearly all large-scale workflows.

Key advantages:
- **Speed**: 400-10,000x faster than blastp at similar sensitivity.
- **Sensitivity**: At `-s 7.5`, matches HMMER sensitivity (but on raw sequence, not profiles).
- **Iterative profile search**: `mmseqs search --num-iterations 3` matches PSI-BLAST behavior at much higher speed.

```bash
# Build target DB
mmseqs createdb target.fasta targetDB
mmseqs createindex targetDB tmp

# Sensitive search
mmseqs easy-search query.fa targetDB results.m8 tmp -s 7.5 --num-iterations 3

# All-vs-all clustering at 50% sequence identity
mmseqs easy-cluster all_proteins.fa cluster tmp --min-seq-id 0.5 -c 0.8
```

The `-s` parameter trades sensitivity for speed: 1.0 (fast), 4.0 (default), 7.5 (HMMER-like sensitivity).

## DIAMOND (the modern blastp replacement)

DIAMOND (Buchfink et al. 2015 *Nat Methods* 12:59; v2: Buchfink et al. 2021 *Nat Methods* 18:366) is the de-facto replacement for blastp on large-scale workflows.

| Feature | DIAMOND v2 | blastp |
|---|---|---|
| Speed | 100-10,000x faster | baseline |
| Sensitivity (default) | ~95% of blastp | baseline |
| `--ultra-sensitive` | >99% of blastp | baseline |
| Frameshift-aware (long reads) | Yes (`--frameshift 15`) | No |
| GPU support | No (CPU-only) | No |

```bash
# Build DIAMOND DB
diamond makedb --in nr.fa -d nr

# Sensitive search
diamond blastp -d nr -q query.fa -o results.tsv \
        --more-sensitive -e 1e-10 -p 16 \
        --outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore stitle

# Long-read frameshift-aware (for nanopore/PacBio metagenomics)
diamond blastx -d nr -q longreads.fa --frameshift 15 -o reads.tsv --outfmt 6
```

For protein remote homology in 2026, DIAMOND `--ultra-sensitive` or MMseqs2 `-s 7.5` should be the default before reaching for BLAST.

## Iterative HMMER (jackhmmer)

```bash
# 3 iterations, save checkpoint HMM at each iteration
jackhmmer -N 3 --chkhmm iter.hmm --tblout hits.tbl query.fa uniref90.fa

# After convergence (no new hits below threshold) use the final HMM for downstream searches
hmmsearch iter-3.hmm target.fa > hits.txt
```

## Code patterns

### Foldseek search against AlphaFoldDB

**Goal:** Find structural homologs of a protein structure (or sequence via ProstT5) in AlphaFold's predicted structure database.

**Approach:** Download AlphaFoldDB structure DB (or ProstT5 for sequence-only); search with foldseek `easy-search`; parse m8 tabular output.

**Reference (Foldseek 9+):**
```bash
#!/bin/bash
# Reference: foldseek 9+ | Verify API if version differs

mkdir -p foldseek_dbs tmp
# Download the AlphaFoldDB Swiss-Prot subset (~few GB; full AFDB is much larger)
foldseek databases Alphafold/Swiss-Prot afdb_sp foldseek_dbs/tmp

# Structure-vs-structure search
foldseek easy-search query.pdb foldseek_dbs/afdb_sp results.m8 tmp \
         --format-output query,target,fident,alnlen,evalue,bits,prob,qtmscore,ttmscore

head -5 results.m8
# qtmscore/ttmscore are TM-score equivalents from local Foldseek alignment.
# Hits with prob > 0.9 are confidently structurally homologous.
```

### Sequence-only Foldseek via ProstT5

```bash
foldseek databases ProstT5 prostt5_model tmp
foldseek easy-search query.fa foldseek_dbs/afdb_sp seq_results.m8 tmp \
         --prostt5-model prostt5_model --threads 8
```

This is the path to take when only a protein sequence is available -- ProstT5 (a protein language model) predicts the 3Di alphabet directly from sequence.

### PSI-BLAST with saved PSSM

**Goal:** Build a position-specific scoring matrix iteratively, then re-use it for downstream searches.

**Approach:** 3 iterations against UniRef90 (or nr); save ASN.1 + ASCII PSSM; subsequent searches use `-in_pssm`.

**Reference (NCBI BLAST+ 2.15+):**
```bash
psiblast -query distant_protein.fa -db uniref90 \
         -num_iterations 3 \
         -inclusion_ethresh 0.002 \
         -evalue 0.01 \
         -num_threads 8 \
         -out_pssm distant.pssm.asn \
         -out_ascii_pssm distant.pssm.txt \
         -out psiblast_results.txt

# Reuse saved PSSM in subsequent searches against a different DB
psiblast -in_pssm distant.pssm.asn -db swissprot \
         -out swissprot_via_pssm.txt
```

### MMseqs2 sensitive iterative search

**Goal:** PSI-BLAST-equivalent iterative profile search, but 100x faster.

**Approach:** `mmseqs search --num-iterations 3 -s 7.5`.

**Reference (MMseqs2 15+):**
```bash
mmseqs createdb query.fa queryDB
mmseqs createdb uniref90.fa uniref90DB
mmseqs createindex uniref90DB tmp

mmseqs search queryDB uniref90DB resultDB tmp \
       --num-iterations 3 \
       -s 7.5 \
       -e 1e-5 \
       --threads 16

mmseqs convertalis queryDB uniref90DB resultDB results.m8 \
       --format-output query,target,fident,alnlen,evalue,bits
```

### Pfam domain annotation (canonical)

```bash
# One-time prep
hmmpress Pfam-A.hmm

# Annotate
hmmscan --cut_ga --domtblout query.domtbl --cpu 8 Pfam-A.hmm query.fa

# Filter: gathering threshold passes are already significance-validated
awk '!/^#/ {print $1, $2, $4, $5, $7, $8, $13}' query.domtbl | head
# columns: target_name, accession, query_name, accession, full_evalue, full_score, i_evalue
```

### HHsearch against PDB70 (deepest homology to PDB)

```bash
# Build query MSA via HHblits vs UniRef30
hhblits -i query.fa -d uniref30 -oa3m query.a3m -n 3 -cpu 8

# Search PDB70 with the query profile
hhsearch -i query.a3m -d pdb70 -o query.hhr -cpu 8

head -30 query.hhr   # Top hits with probability + alignment statistics
```

### DIAMOND ultra-sensitive on a metagenome

```bash
diamond makedb --in uniref90.fa -d uniref90
diamond blastp -d uniref90 -q metagenome_proteins.fa -o hits.tsv \
        --ultra-sensitive -e 1e-5 -p 32 \
        --outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore stitle
```

## Failure modes

### PSI-BLAST profile drift
- **Trigger:** Iterating to convergence (5+ iterations).
- **Mechanism:** Each iteration includes hits below threshold; eventually paralogs and divergent family members contaminate the PSSM.
- **Symptom:** Later iterations return many implausible hits; functional inference goes wrong.
- **Fix:** Cap at 3 iterations; inspect the saved PSSM and the included sequence set; use stricter `-inclusion_ethresh 0.001`.

### Foldseek "structure but no homology" hits
- **Trigger:** Searching small fragments or highly conserved folds (TIM barrels, Rossmann folds).
- **Mechanism:** Structural fold is preserved across deep divergence; superfamily hits exist without true homology.
- **Symptom:** High structural similarity to functionally unrelated proteins.
- **Fix:** Combine Foldseek hits with sequence-based evidence; check shared catalytic residues; consider that fold-level similarity is necessary but not sufficient for homology.

### MMseqs2 default sensitivity
- **Trigger:** `mmseqs easy-search` without `-s`.
- **Mechanism:** Default `-s 4.0` is fast but misses remote homologs.
- **Symptom:** Equivalent to a fast BLAST; misses what HMMER would find.
- **Fix:** Set `-s 7.5` for distant homology; `-s 5.7` is a middle ground.

### DIAMOND default mode lossy
- **Trigger:** `diamond blastp` without `--more-sensitive` or `--ultra-sensitive`.
- **Mechanism:** Default mode trades ~5% sensitivity for speed vs blastp.
- **Symptom:** Hits BLAST would find are missing.
- **Fix:** Use `--more-sensitive` for general work, `--ultra-sensitive` for remote homology.

### Profile method on a low-complexity query
- **Trigger:** Query has signal peptide, coiled-coil, or repeat region.
- **Mechanism:** Profile is dominated by low-complexity columns; false hits to other low-complexity proteins.
- **Symptom:** Many high-scoring hits to unrelated low-complexity proteins.
- **Fix:** Mask the low-complexity region (SEG: `segmasker -infmt fasta -in query.fa`) before building the profile.

### HHblits database version drift
- **Trigger:** Using a UniRef30 database from a different release than the PDB70 search DB.
- **Mechanism:** Profile statistics depend on the DB's amino acid distribution.
- **Symptom:** Hit probabilities are miscalibrated.
- **Fix:** Use UniRef30 and PDB70 from the same MMseqs2 / HH-suite release.

### Foldseek without ProstT5 for sequence-only query
- **Trigger:** No structure available; tried to predict with AF2 first (slow).
- **Mechanism:** ProstT5 predicts 3Di alphabet directly from sequence, skipping AF2.
- **Symptom:** Days-long AF2 prediction step for a query that could be Foldseek'd in seconds.
- **Fix:** Use `foldseek databases ProstT5 ...` and `--prostt5-model`.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| PSI-BLAST returns implausible hits | Profile drift (too many iterations) | Cap at 3 iterations; tighter `-inclusion_ethresh` |
| MMseqs2 hits all unrelated | Default sensitivity too low | `-s 7.5` |
| DIAMOND misses BLAST hits | Default mode lossy | `--more-sensitive` |
| Foldseek hits structurally unrelated proteins | Common fold, no homology | Cross-check with sequence and functional residues |
| HHblits prefilter no hits | Query MSA too sparse | Add `-n 4` iterations; check input |
| jackhmmer ConvergenceError | Loop bug pre-v3.4 | Upgrade HMMER |

## References

- Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. (1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Res* 25:3389-3402.
- Rost B. (1999) Twilight zone of protein sequence alignments. *Protein Eng* 12:85-94.
- Eddy SR. (2011) Accelerated profile HMM searches. *PLoS Comput Biol* 7:e1002195.
- Remmert M, Biegert A, Hauser A, Soding J. (2012) HHblits: lightning-fast iterative protein sequence searching by HMM-HMM alignment. *Nat Methods* 9:173-175.
- Buchfink B, Xie C, Huson DH. (2015) Fast and sensitive protein alignment using DIAMOND. *Nat Methods* 12:59-60.
- Steinegger M, Soding J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol* 35:1026-1028.
- Steinegger M, Meier M, Mirdita M, Vohringer H, Haunsberger SJ, Soding J. (2019) HH-suite3 for fast remote homology detection and deep protein annotation. *BMC Bioinformatics* 20:473.
- Buchfink B, Reuter K, Drost HG. (2021) Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods* 18:366-368.
- van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Soding J, Steinegger M. (2024) Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* 42:243-246.
- Heinzinger M, Weissenow K, Sanchez JG, Henkel A, Mirdita M, Steinegger M, Rost B. (2024) Bilingual language model for protein sequence and structure. *NAR Genom Bioinform* 6:lqae150.

## Related Skills

- blast-searches - Remote BLAST against NCBI; baseline for closer homologs
- local-blast - Local BLAST+ for moderate-scale workflows
- ortholog-inference - Orthology calls (RBH, OrthoFinder, OMA, Compara)
- alignment/multiple-alignment - Build MSAs for HMM profiles
- structural-biology/alphafold-predictions - Predict structures for Foldseek input
- structural-biology/modern-structure-prediction - ESMFold, ColabFold pipelines
