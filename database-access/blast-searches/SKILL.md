---
name: bio-blast-searches
description: Run remote BLAST searches against NCBI servers using Biopython Bio.Blast.NCBIWWW. Use when identifying unknown sequences, finding homologs, picking the correct BLAST program (blastn/blastp/blastx/tblastn/tblastx/psiblast/megablast/dc-megablast), interpreting Karlin-Altschul E-values, avoiding the max_target_seqs trap (Shah 2019), choosing composition-based statistics, or limiting searches by organism. Covers RID lifecycle, database choice (nt/nr/refseq_select/swissprot), word-size and CBS taxonomy.
tool_type: python
primary_tool: Bio.Blast.NCBIWWW
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, NCBI BLAST+ 2.15+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Blast.NCBIWWW.qblast)` to check signatures
- CLI: `blastn -version` then `blastn -help`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# BLAST Searches (Remote)

**"Find similar sequences in NCBI's database"** -> Submit a query to NCBI's remote BLAST servers; receive a Request ID (RID); poll for completion; parse the XML hit table. Best for one-off identification of a few sequences. For >50 sequences, switch to `local-blast` or DIAMOND/MMseqs2 in `remote-homology`.

The two most consequential decisions: **which program** (defines query+target molecule types and word-size defaults) and **which database** (defines the search space and therefore E-value baselines). The third most important: do NOT misuse `max_target_seqs` -- it is an early-termination heuristic, not a "give me the top N hits" filter (Shah et al. 2019).

- Python: `NCBIWWW.qblast(program, db, sequence)` + `NCBIXML.read(handle)` (BioPython)
- CLI: `blastn -remote -db nt -query seq.fa -out hits.xml -outfmt 5` (BLAST+)
- Web: https://blast.ncbi.nlm.nih.gov/Blast.cgi (RID lookup)

## Required Setup

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
```

No API key needed for remote BLAST itself, but NCBI's general rate-limit ethic still applies -- one search at a time, polite waiting, no parallelism.

## Program decision (query vs database molecule)

| Program | Query | Target | Word size default | Use case |
|---|---|---|---|---|
| `blastn` | DNA | DNA | 11 | General DNA similarity |
| `megablast` | DNA | DNA | 28 | High-identity DNA (>=95%) -- PCR primer hits, contamination |
| `dc-megablast` | DNA | DNA | 11 (discontiguous) | Cross-species mRNA (sensitive, gapped) |
| `blastp` | Protein | Protein | 3 (6 also valid) | General protein homology |
| `blastx` | DNA | Protein | 3 | Translated DNA query vs protein DB; ORF discovery |
| `tblastn` | Protein | DNA | 3 | Protein query vs translated DB; find unannotated CDS |
| `tblastx` | DNA | DNA | 3 (both translated) | Most expensive; deep cross-species coding similarity |
| `psiblast` | Protein | Protein | 3 | Iterative PSSM-based remote homology -- see `remote-homology` |

**The misuse to avoid:** using default `blastn` (word=11) for cross-species DNA where `dc-megablast` is the right tool. Or using `megablast` (word=28) for cross-species homology where it will miss every divergent hit. The most-misused BLAST parameter according to literature.

## Database decision (search space)

| Database (`db=`) | Content | Size (2026 approx) | Stable for reproducibility? |
|---|---|---|---|
| `nt` | Non-redundant nucleotide (all GenBank+EMBL+DDBJ) | ~250 GB | NO -- changes daily |
| `nr` | Non-redundant protein | ~300 GB | NO -- changes daily |
| `refseq_select` | One curated rep per species (RNA + protein) | small | YES -- versioned releases |
| `refseq_rna` | RefSeq mRNA | ~10 GB | YES |
| `refseq_protein` | RefSeq protein | small | YES |
| `swissprot` | UniProt Swiss-Prot (reviewed) | small | YES -- monthly releases |
| `pdb` | Protein structures | small | YES |
| `refseq_genomic` | RefSeq genomic | huge | YES |
| `env_nr` / `env_nt` | Environmental (metagenomic) | huge | YES |

**For publication reproducibility, never search `nt` or `nr`** without recording the snapshot date and ideally archiving a frozen copy. Default to `refseq_select` for any cross-species homology question; switch to `nt`/`nr` only when curated coverage is insufficient.

## E-value interpretation (Karlin-Altschul)

E-value = K * m * n * exp(-lambda * S), where m = effective query length, n = effective database size, lambda and K are scoring-matrix-dependent constants (Karlin & Altschul 1990 PNAS 87:2264).

| E-value | Bit-score (BLOSUM62, protein) | Interpretation |
|---|---|---|
| < 1e-50 | > 200 | Strong; almost certainly homologous |
| 1e-50 to 1e-10 | 100-200 | Significant; likely homolog |
| 1e-10 to 1e-3 | 50-100 | Marginal; check identity + coverage |
| 0.01 to 10 | 30-50 | Possible remote homolog; needs profile method |
| > 10 | < 30 | Random; not meaningful |

**Key implication of E = K * m * n * exp(-lambda * S):** the same alignment against a 100x larger database has a 100x larger E-value. Cross-database E-value comparison is meaningless. Bit-score is database-size normalized and is the right cross-database metric.

For protein remote homology where E is marginal (10^-3 to 10^-1), reach for profile methods: PSI-BLAST, jackhmmer, HHblits, or Foldseek -- see `remote-homology` skill.

## Composition-Based Statistics (CBS)

Compositional bias inflates significance for low-complexity proteins. The CBS modes (Yu et al. 2006 *Nucleic Acids Res* 34:5966):

| `composition_based_statistics` | Mode | Use when |
|---|---|---|
| 0 | Off | Almost never |
| 1 | F&S 2002 score adjustment | Legacy compatibility |
| 2 | Yu&Altschul 2005 conditional score adjustment | **Default since BLAST+ 2.2.17** -- correct for most cases |
| 3 | Universal statistics | Short queries (< 30 aa) where mode 2 over-corrects |

For protein queries under 30 aa, switch to CBS=3. For protein with known compositional bias (e.g. coiled-coil regions, signal peptides), CBS=2 is appropriate but consider hard-masking with SEG.

## The `max_target_seqs` trap

**The misuse**: `max_target_seqs=10` is interpreted as "return the 10 most significant hits". It is not. The flag is an **early termination** parameter that affects which hits the search ever considers, not which it ultimately reports (Shah N, Nute MG, Warnow T, Pop M. (2019) Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows. *Bioinformatics* 35:1613-1614).

**Consequences:**
- Setting `max_target_seqs=10` can return entirely different hits than `max_target_seqs=500` then filtering to top 10 by E-value.
- The "top 10" by E-value as reported may not be the actual top 10.

**Correct pattern:** set `hitlist_size` (Bio.Blast parameter name) large (1000+), then post-filter to the top N by E-value or bit-score in Python.

## Word size, gap costs, and matrix

| Search | Word size | Matrix (protein) | Gap (open, extend) |
|---|---|---|---|
| megablast (high identity DNA) | 28 | n/a | 0, 0 (linear) |
| blastn (sensitive DNA) | 11 | n/a | 5, 2 |
| blastp default | 3 | BLOSUM62 | 11, 1 |
| blastp distant | 2 | BLOSUM45 | 14, 2 |
| Short peptides (<30 aa) | 2 | PAM30 or BLOSUM45 | 9, 1 |

For very short query proteins (e.g. proteomics-identified peptides), BLOSUM45 + word=2 + PAM30 substitution matrix is more sensitive than the default. Use `matrix='PAM30'` for searches against `swissprot`.

## RID lifecycle

| Phase | Server state | Client action |
|---|---|---|
| Submit | RID created, queued | NCBIWWW.qblast() returns handle |
| Running | Queue + compute | Poll status |
| Done | RID + results retained | Fetch XML |
| Expired | RID purged | 24-36h after completion |

`NCBIWWW.qblast()` handles polling internally with a fixed retry interval. For long-running searches (>5 min) or batches, submit and capture the RID, then poll independently to avoid blocking. The RID is visible at `https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=...` for 24-36 hours.

## Code patterns

### Standard remote BLASTN with reproducible parameters

**Goal:** Run BLASTN with explicit, paper-quality parameters.

**Approach:** Specify program, database (refseq_select for stability), word size, expect, and a large hitlist_size to dodge the max_target_seqs trap.

**Reference (BioPython 1.83+):**
```python
from Bio.Blast import NCBIWWW, NCBIXML

handle = NCBIWWW.qblast(
    program='blastn',
    database='refseq_select_rna',
    sequence=query_seq,
    expect=1e-10,
    word_size=11,
    hitlist_size=500,  # large; filter top-N downstream
    format_type='XML',
)
record = NCBIXML.read(handle); handle.close()
top10 = sorted(record.alignments, key=lambda a: a.hsps[0].expect)[:10]
```

### Protein search with organism restriction

**Goal:** Find mammalian homologs of a query protein in Swiss-Prot.

**Approach:** `entrez_query` filters the BLAST search space pre-execution; faster and more meaningful E-values than post-filtering.

**Reference (BioPython 1.83+):**
```python
handle = NCBIWWW.qblast(
    program='blastp',
    database='swissprot',
    sequence=protein_seq,
    entrez_query='Mammalia[Organism]',
    expect=1e-5,
    composition_based_statistics=2,
    hitlist_size=200,
)
record = NCBIXML.read(handle); handle.close()
```

### Short peptide search

```python
handle = NCBIWWW.qblast(
    program='blastp',
    database='swissprot',
    sequence=peptide_seq,  # < 30 aa
    matrix_name='PAM30',
    word_size=2,
    expect=1000,  # short queries need permissive cutoff
    composition_based_statistics=3,
    hitlist_size=100,
)
```

### Save XML for re-parsing

```python
handle = NCBIWWW.qblast('blastn', 'refseq_select_rna', query)
with open('blast.xml', 'w') as f:
    f.write(handle.read())
handle.close()

with open('blast.xml') as f:
    record = NCBIXML.read(f)
```

### Hit extraction with identity + coverage filtering

**Goal:** Return structured top hits with biological metrics, not just E-values.

**Approach:** Walk alignments + first HSP; compute identity and query coverage as fractions; sort by bit-score (database-size invariant) not E-value.

**Reference (BioPython 1.83+):**
```python
def top_hits(record, min_identity=0.5, min_coverage=0.7, top_n=10):
    qlen = record.query_length
    hits = []
    for aln in record.alignments:
        hsp = aln.hsps[0]
        ident = hsp.identities / hsp.align_length
        cov = hsp.align_length / qlen
        if ident >= min_identity and cov >= min_coverage:
            hits.append({
                'accession': aln.accession,
                'title': aln.title,
                'evalue': hsp.expect,
                'bits': hsp.bits,
                'identity': ident,
                'coverage': cov,
            })
    return sorted(hits, key=lambda h: -h['bits'])[:top_n]
```

### Programmatic RID polling for long jobs

```python
import time

handle = NCBIWWW.qblast('tblastn', 'nr', query, hitlist_size=500, format_type='XML')
# Bio.Blast handles polling internally; for explicit control use the REST API directly
# or save and re-parse the RID URL
```

## Failure modes

### `max_target_seqs` misinterpretation
- **Trigger:** Setting `hitlist_size=10` and assuming top 10 by E-value.
- **Mechanism:** It's an early-termination param; can miss legitimate top hits.
- **Symptom:** Different "top 10" between hitlist=10 and hitlist=500 filtered.
- **Fix:** Always set `hitlist_size=500+` and post-filter; cite Shah 2019.

### Cross-database E-value comparison
- **Trigger:** Comparing E from a `nt` search against E from a `swissprot` search.
- **Mechanism:** E scales linearly with database size; comparison is meaningless.
- **Symptom:** Misleading rankings between two analyses.
- **Fix:** Compare bit-scores instead, or set the same database for both.

### Megablast for cross-species
- **Trigger:** Default `megablast` (word=28) on a cross-species DNA query.
- **Mechanism:** Word size 28 requires 28-nt exact match to seed; cross-species mRNA has too much divergence.
- **Symptom:** Zero hits or only hits to the same species.
- **Fix:** Use `dc-megablast` (discontiguous) or `blastn` with word=11.

### Reproducibility loss against `nt`/`nr`
- **Trigger:** Manuscript says "BLASTed against nt"; reviewer re-runs 3 weeks later.
- **Mechanism:** Databases change daily; new genomes deposited.
- **Symptom:** Different hit set, different paper conclusions.
- **Fix:** Use `refseq_select` for reproducibility, or record snapshot date + archive subset.

### Server timeout on large queries
- **Trigger:** Multi-megabase query or batch submission.
- **Mechanism:** Remote BLAST has a per-query compute budget.
- **Symptom:** Job stuck in queue, eventually fails.
- **Fix:** Split into smaller queries; or switch to `local-blast` / DIAMOND / MMseqs2.

### Compositional bias inflates E
- **Trigger:** Protein query with low-complexity region (coiled-coil, signal peptide).
- **Mechanism:** Default CBS=2 handles most cases, but extreme bias still inflates scores.
- **Symptom:** Many "significant" hits to unrelated low-complexity proteins.
- **Fix:** Confirm CBS=2 is on; consider hard-masking with `filter='S'` (SEG).

### Empty FASTA defline submitted
- **Trigger:** Sending `sequence` as a raw string without `>id\n`.
- **Mechanism:** BLAST treats as anonymous query; some downstream parsers misbehave.
- **Symptom:** Hits returned but `record.query` is None.
- **Fix:** Always pass FASTA with a defline; or pass a `SeqRecord`.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Stuck > 5 min | Large query or busy queue | Submit RID, poll separately; or use local |
| URLError / timeout | Network or NCBI maintenance | Retry with backoff; status at status.ncbi.nlm.nih.gov |
| No hits | Wrong program / database type | Verify query and DB molecule types match |
| Empty XML | RID expired | Re-submit; RIDs purge after 24-36h |
| 1000s of low-complexity hits | CBS disabled or extreme bias | CBS=2; consider SEG filter |
| Cross-DB E mismatch | Comparing E across DBs | Use bit-score instead |

## References

- Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. (1990) Basic local alignment search tool. *J Mol Biol* 215:403-410.
- Karlin S, Altschul SF. (1990) Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. *Proc Natl Acad Sci USA* 87:2264-2268.
- Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. (1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Res* 25:3389-3402.
- Yu YK, Gertz EM, Agarwala R, Schaffer AA, Altschul SF. (2006) Retrieval accuracy, statistical significance and compositional similarity in protein sequence database searches. *Nucleic Acids Res* 34:5966-5973.
- Shah N, Nute MG, Warnow T, Pop M. (2019) Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows. *Bioinformatics* 35:1613-1614.
- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. (2009) BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.

## Related Skills

- local-blast - Faster, unlimited local BLAST+ pipelines and database build
- remote-homology - PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek for distant homology
- ortholog-inference - Reciprocal best hit, OrthoFinder, OMA for orthology
- sequence-io/read-sequences - Load query sequences from FASTA
- entrez-fetch - Fetch full records for BLAST hits
