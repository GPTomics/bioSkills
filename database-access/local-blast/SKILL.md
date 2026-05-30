---
name: bio-local-blast
description: Build local BLAST databases and run searches using NCBI BLAST+ command-line tools. Use when running >50 queries, building custom databases with -parse_seqids and -taxid, downloading prebuilt NCBI databases via update_blastdb.pl, choosing -task variants (megablast/dc-megablast/blastn/blastn-short), tuning soft/hard masking, scaling threads, or extracting hits with blastdbcmd. Encodes BLAST v5 vs v4 database format, taxonomy filtering, makeblastdb pitfalls.
tool_type: cli
primary_tool: BLAST+
---

## Version Compatibility

Reference examples tested with: NCBI BLAST+ 2.15+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `blastn -version` then `blastn -help` to confirm flags
- CLI: `makeblastdb -help` to confirm database build options

If a flag is unrecognized or behavior changes, introspect with `-help` and adapt the example to match the installed version rather than retrying.

# Local BLAST

**"Run BLAST locally for speed and control"** -> Build or download a BLAST+ database, run the appropriate program with carefully chosen `-task`, masking, and thread settings, parse tabular output. Local BLAST is the right tool when remote is rate-limited or when the database must be reproducible (frozen).

The biggest mistakes are (a) using `nt`/`nr` without realizing they're >250 GB and grow weekly, (b) not building with `-parse_seqids` and then being unable to extract hit sequences with `blastdbcmd`, (c) using default `blastn` for cross-species when `dc-megablast` is correct, and (d) thinking `-num_threads 32` will scale -- past ~16 threads BLAST is I/O bound.

- CLI: `makeblastdb`, `blastn`/`blastp`, `blastdbcmd`, `update_blastdb.pl` (NCBI BLAST+)
- Python: `subprocess` wrapper (preferred); `Bio.Blast.Applications` was deprecated and removed -- do not use

## Installation

```bash
# conda (preferred)
conda install -c bioconda blast

# macOS
brew install blast

# Ubuntu
sudo apt install ncbi-blast+

# Verify
blastn -version    # NCBI BLAST+ 2.15+ expected
update_blastdb.pl --showall pretty | head
```

## Database format: v5 vs v4

NCBI introduced BLAST database v5 in BLAST+ 2.10 (2020). v5 includes taxonomy indexing directly in the database files, enabling `-taxids` and `-taxidlist` filtering without a companion file. v4 databases require `taxonomy4blast.sqlite3` to be present and discoverable.

| Feature | v4 | v5 |
|---|---|---|
| Default for prebuilt NCBI dbs | No (legacy) | Yes (since 2020) |
| `-taxids`, `-taxidlist` support | No | Yes |
| `blastdbcmd -taxids` | No | Yes |
| New `-info` output fields | No | Yes |

`update_blastdb.pl` downloads v5 by default. When building a database manually with `makeblastdb`, v5 format requires `-blastdb_version 5`. **Always pass `-blastdb_version 5` and `-parse_seqids` when building from scratch.**

## `makeblastdb` flag taxonomy

| Flag | Effect | When |
|---|---|---|
| `-dbtype nucl` or `-dbtype prot` | Required | Always |
| `-parse_seqids` | Indexes accessions so `blastdbcmd -entry <acc>` works | Almost always (downstream extraction) |
| `-hash_index` | Speeds up extraction by accession | Large dbs |
| `-blastdb_version 5` | Use v5 format | Always |
| `-taxid 9606` | Single taxid for all seqs | Single-species DB |
| `-taxid_map file.tsv` | Per-sequence taxid mapping (seqid<TAB>taxid) | Multi-species DB |
| `-mask_data masking.asnb` | Apply precomputed soft-masking | Production pipelines |
| `-title "..."` | Free-text label | Cosmetic |
| `-out path/prefix` | DB file path prefix | Always |

```bash
makeblastdb -in reference.fasta -dbtype nucl \
            -blastdb_version 5 \
            -parse_seqids \
            -hash_index \
            -title "Custom reference 2026-05" \
            -out custom_db
```

## `-task` taxonomy (the most-misused BLAST setting)

For `blastn`, the `-task` flag picks among heuristics with different word sizes and gap parameters.

| `-task` | Word | Gapped | Use case | Mistake to avoid |
|---|---|---|---|---|
| `megablast` (default) | 28 | linear | >=95% identity, intra-species, primer hits, contamination check | Used for cross-species and misses everything |
| `dc-megablast` | 11 (discontiguous) | yes | Cross-species mRNA homology | Underused -- this is what `blastn` "should" be for cross-species |
| `blastn` | 11 | yes | General sensitive DNA | Slower than dc-megablast for same job |
| `blastn-short` | 7 | yes | Queries <50 nt (primers, small RNAs) | Default megablast can't seed at length 7 |
| `rmblastn` | 11 | yes | Repeat masking; bundled with RepeatModeler | Specialized |

For `blastp`:

| `-task` | Word | Use case |
|---|---|---|
| `blastp` (default) | 3 | General protein similarity |
| `blastp-fast` | 6 | Faster, less sensitive |
| `blastp-short` | 2 | Peptides <30 aa, with PAM30 + word_size=2 typical |

## Soft vs hard masking

| Setting | Effect on seed | Effect on extension | Effect on score |
|---|---|---|---|
| `-soft_masking true` (default for several tasks) | Skip masked positions when seeding | Allow extension through masked | Score includes masked positions |
| `-soft_masking false` + `-dust yes` / `-seg yes` | Skip masked positions when seeding | Skip masked positions in extension | Score excludes masked positions |
| Hard-mask in input FASTA (N or X) | Hard exclusion everywhere | Hard exclusion | Treated as mismatches |

Soft masking is correct for almost all cases. Hard masking creates artificial mismatches at masked boundaries and can split true alignments. The exception: searching against a database of repeats explicitly, where hard masking on the query is the right choice.

## Thread scaling

BLAST+ parallelizes per-query (with `-num_threads`) but is I/O bound past ~16 threads on most hardware. For >100,000 query batches the better answer is splitting the input FASTA into N chunks and running N parallel `blastn` invocations -- this saturates CPUs better than `-num_threads 64`.

| Threads | Typical speedup vs single | Notes |
|---|---|---|
| 1-8 | Near-linear | Default sweet spot |
| 8-16 | Sub-linear (1.5-2x over 8) | Useful on big SMP boxes |
| 16-32 | Diminishing returns | I/O bound for most DBs |
| 32+ | Often slower | Cache thrash + I/O contention |

For massive workflows, prefer **DIAMOND** (Buchfink et al. 2021 *Nat Methods* 18:366) or **MMseqs2** (Steinegger & Soding 2017 *Nat Biotechnol* 35:1026) -- 100-10,000x faster than BLASTP at comparable sensitivity. See `remote-homology` skill.

## Output format reference (`-outfmt`)

| `-outfmt` | Description | Use |
|---|---|---|
| 0 | Pairwise (default; human-readable) | Debugging, inspection |
| 5 | XML | Programmatic parsing (Bio.SearchIO) |
| 6 | Tabular (no header) | Most pipelines |
| 7 | Tabular with comment headers | Self-documenting |
| 11 | ASN.1 binary | Re-parse with later versions |

Custom tabular fields:
```bash
blastn -query q.fa -db db -outfmt "6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore staxids sscinames stitle"
```

Field key fields for analysis:
- `pident` = percent identity over the HSP (NOT the query); for query-level, use `qcovhsp`
- `qcovs` = total query coverage by all HSPs of this subject (the "coverage" most users want)
- `qcovhsp` = query coverage by best HSP alone (use when there's only one HSP per hit)
- `staxids` = taxonomy IDs (v5 only); critical for any "what species" workflow

## Prebuilt NCBI databases via `update_blastdb.pl`

```bash
# List available
update_blastdb.pl --showall pretty | grep -E 'refseq|swissprot|nt|nr'

# Download (with decompress)
update_blastdb.pl --decompress refseq_select_rna

# Download specific volume of split database
update_blastdb.pl --decompress refseq_protein

# Download with parallelism
update_blastdb.pl --decompress --num_threads 4 refseq_select_rna
```

Sizes (approximate, 2026):
- `refseq_select_rna`: ~5 GB
- `refseq_protein`: ~30 GB
- `swissprot`: <1 GB
- `nt`: ~250 GB
- `nr`: ~300 GB

For most use cases, `refseq_select_*` is the right starting point. `nt`/`nr` are storage-heavy and reproducibility-hostile.

## Code patterns

### Build and search a custom protein database

**Goal:** Build a BLAST+ protein database from a custom FASTA and search against it.

**Approach:** `makeblastdb` with v5 + parse_seqids + hash_index; `blastp` with explicit outfmt.

**Reference (NCBI BLAST+ 2.15+):**
```bash
#!/bin/bash
# Reference: NCBI BLAST+ 2.15+ | Verify API if version differs

REF=reference_proteins.fasta
DB=ref_prot_db
QUERY=query.fasta
OUT=hits.tsv

makeblastdb -in "$REF" -dbtype prot \
            -blastdb_version 5 -parse_seqids -hash_index \
            -title "$REF $(date +%Y-%m-%d)" \
            -out "$DB"

blastp -query "$QUERY" -db "$DB" \
       -evalue 1e-10 \
       -num_threads 8 \
       -max_target_seqs 500 \
       -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle" \
       -out "$OUT"

# Top hit per query by bit-score (column 7)
sort -k1,1 -k7,7gr "$OUT" | awk '!seen[$1]++' > top_hit_per_query.tsv
```

### Cross-species DNA with dc-megablast

```bash
blastn -query mouse_cdna.fa -db human_refseq_rna \
       -task dc-megablast \
       -word_size 11 \
       -evalue 1e-10 \
       -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore" \
       -num_threads 8 \
       -out cross_species.tsv
```

### Short primer search

```bash
blastn -query primers.fa -db genome_db \
       -task blastn-short \
       -word_size 7 \
       -evalue 1000 \
       -outfmt 6 \
       -out primer_hits.tsv
```

### Taxonomy-filtered search (BLAST v5 only)

```bash
# Restrict to specific taxids
blastp -query query.fa -db nr \
       -taxids 9606,10090,10116 \
       -outfmt "6 qseqid sseqid staxids sscinames evalue bitscore" \
       -out mammalian_hits.tsv

# Or to a taxid subtree (NCBI BLAST+ 2.13+)
echo 9606 > human_only.txt
blastp -query query.fa -db nr -taxidlist human_only.txt -outfmt 6 -out human_hits.tsv
```

### Extract subject sequences for top hits

```bash
# Requires database built with -parse_seqids
cut -f2 top_hit_per_query.tsv | sort -u > hit_accessions.txt
blastdbcmd -db ref_prot_db -entry_batch hit_accessions.txt -out hits.fasta

# Pull a range of a sequence
blastdbcmd -db genome_db -entry NC_000001.11 -range 1000000-1001000 -out region.fa
```

### Reciprocal best hit (RBH) for ortholog candidates

See `ortholog-inference` skill for the principled treatment. Quick version:

```bash
blastp -query A.fa -db B_db -outfmt 6 -evalue 1e-5 -num_threads 8 \
       -max_target_seqs 5 -out A_vs_B.tsv
blastp -query B.fa -db A_db -outfmt 6 -evalue 1e-5 -num_threads 8 \
       -max_target_seqs 5 -out B_vs_A.tsv

# Best forward + reverse, intersect
awk '!seen[$1]++ {print $1"\t"$2}' A_vs_B.tsv | sort > A_best
awk '!seen[$1]++ {print $1"\t"$2}' B_vs_A.tsv | sort > B_best
awk 'NR==FNR{a[$1]=$2; next} a[$2]==$1' A_best B_best > rbh.tsv
```

This works but does NOT handle paralog mis-pairs from gene duplication; for that use OrthoFinder or OMA (in `ortholog-inference`).

### Python wrapper with version pinning

```python
import subprocess
import shutil


def require_tool(name, min_version=None):
    if not shutil.which(name):
        raise RuntimeError(f'{name} not on PATH')
    out = subprocess.run([name, '-version'], capture_output=True, text=True)
    print(f'  {out.stdout.strip().splitlines()[0]}')


def run_blast(query, db, out, program='blastp', evalue=1e-10, threads=8, hitlist=500):
    require_tool(program)
    cmd = [program, '-query', query, '-db', db, '-out', out,
           '-evalue', str(evalue),
           '-num_threads', str(threads),
           '-max_target_seqs', str(hitlist),
           '-outfmt', '6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore stitle']
    subprocess.run(cmd, check=True)


def parse_tabular(path):
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'qcovs', 'qcovhsp', 'evalue', 'bitscore', 'stitle']
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.rstrip('\n').split('\t')
            d = dict(zip(cols, vals))
            for k in ('pident', 'qcovs', 'qcovhsp', 'evalue', 'bitscore'):
                d[k] = float(d[k])
            d['length'] = int(d['length'])
            rows.append(d)
    return rows
```

## Failure modes

### `nt`/`nr` size shock
- **Trigger:** `update_blastdb.pl --decompress nt` without realizing the size.
- **Mechanism:** `nt` is ~250 GB compressed, ~1 TB indexed.
- **Symptom:** Disk fills mid-download; partial DB unusable.
- **Fix:** Use `refseq_select` for most workflows; only pull `nt`/`nr` with intent and >1 TB free.

### Missing `-parse_seqids`
- **Trigger:** Built DB without `-parse_seqids`; later try `blastdbcmd -entry`.
- **Mechanism:** Without the parsed index, `blastdbcmd` can't look up by accession.
- **Symptom:** `Error: ... not found in database`.
- **Fix:** Rebuild with `-parse_seqids` (cheap if FASTA still on disk).

### Wrong `-task` for the question
- **Trigger:** Default `blastn` for cross-species mRNA (word=11 but ungapped seeding).
- **Mechanism:** Discontiguous seed (`dc-megablast`) is much more sensitive across species.
- **Symptom:** Far fewer hits than the question warrants.
- **Fix:** Use `-task dc-megablast` for cross-species; `-task megablast` only for >=95% identity.

### Thread saturation
- **Trigger:** `-num_threads 64` on a 32-core box.
- **Mechanism:** I/O bound past ~16 threads; cache thrash hurts past CPU count.
- **Symptom:** No speedup or slowdown.
- **Fix:** Cap at 8-16; split FASTA and run parallel processes instead for very large batches.

### v4 database, expecting v5 features
- **Trigger:** Old prebuilt DB; `-taxids` flag returns "Taxonomy database not available".
- **Mechanism:** v4 needs `taxonomy4blast.sqlite3` companion; v5 has taxonomy indexed in DB.
- **Symptom:** Taxonomy filtering silently no-ops or errors.
- **Fix:** Re-download with `update_blastdb.pl --decompress` (gets v5); or use v5 explicitly when building.

### Soft-masking confusion
- **Trigger:** Hard-masking input (replacing repeats with N or X) instead of using `-dust`/`-seg`.
- **Mechanism:** Hard-mask creates artificial mismatches at boundaries.
- **Symptom:** True alignments split into multiple short HSPs.
- **Fix:** Pass unmasked FASTA + soft-mask via `-soft_masking true` + `-dust yes`/`-seg yes`.

### `max_target_seqs` truncation
- **Trigger:** `-max_target_seqs 10` (Shah et al. 2019 *Bioinformatics* 35:1613).
- **Mechanism:** Early termination, not top-N filter.
- **Symptom:** Different top-10 than `-max_target_seqs 500` + post-filter.
- **Fix:** Set `-max_target_seqs` large (500+); filter top N in awk/Python.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `BLAST Database error` | DB path wrong, or alias missing | `blastdbcmd -db <db> -info` to confirm |
| `Error: entry not found` | Built without `-parse_seqids` | Rebuild |
| Taxonomy filter no-op | v4 DB | Upgrade to v5 |
| Threads >16 not faster | I/O bound | Split input + parallel invocations |
| `nt` download fills disk | Database is huge | Use refseq_select |
| `Sequence too short` | Query < word_size | Use `-task blastn-short` (word=7) |
| Out of memory | Single large query | Reduce `-num_threads`, split query |

## References

- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. (2009) BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
- Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. (1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Res* 25:3389-3402.
- Shah N, Nute MG, Warnow T, Pop M. (2019) Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows. *Bioinformatics* 35:1613-1614.
- Boratyn GM, Camacho C, Cooper PS, et al. (2013) BLAST: a more efficient report with usability improvements. *Nucleic Acids Res* 41:W29-W33.

## Related Skills

- blast-searches - Remote BLAST against NCBI servers
- remote-homology - PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek for distant homology
- ortholog-inference - Reciprocal best hit, OrthoFinder, OMA for ortholog calls
- sequence-io/read-sequences - Load query/reference FASTA files
- batch-downloads - Download large reference FASTA sets before makeblastdb
