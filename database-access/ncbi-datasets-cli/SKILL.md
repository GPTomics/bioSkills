---
name: bio-ncbi-datasets-cli
description: Download genome assemblies, gene records, and ortholog data from NCBI using the modern Datasets v2 CLI (replaces assembly_summary.txt scraping and many EFetch workflows). Use when bulk-pulling genome assemblies, gene metadata across species, ortholog sets, or BLAST databases; when E-utilities are too slow for genome-scale work; or when automatic checksum verification, parallel download, and clean accession-driven retrieval are required. Encodes the JSON-lines output format, dataformat conversion, --dehydrated for cloud workflows, and when Datasets is/isn't the right tool.
tool_type: cli
primary_tool: NCBI Datasets CLI
---

## Version Compatibility

Reference examples tested with: NCBI Datasets CLI 16.0+ (2024), dataformat 16.0+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `datasets --version`, `dataformat --version`
- Subcommand help: `datasets <subcommand> --help`

If a subcommand or flag is unrecognized, run `datasets --help` and adapt. The CLI is under active development; major releases (v15 -> v16) added subcommands and renamed flags.

# NCBI Datasets CLI

**"Pull genome / gene / ortholog data from NCBI in 2026"** -> The Datasets v2 CLI (launched 2023) is the official, supported bulk endpoint for genome and gene-centric data. It replaces the prior best-practice of scraping `assembly_summary.txt` + parallel FTP + manual checksum verification. For genome-scale data, it is strictly better than E-utilities (EFetch).

The CLI is not the right answer for everything. PubMed, SRA reads, and custom Entrez queries still belong to E-utilities. The defection rule: **if the question is about genome assemblies, gene records, or pre-computed orthologs, use Datasets; otherwise stay with E-utilities**.

- CLI: `datasets download genome accession GCF_...`
- CLI: `datasets summary gene symbol BRCA1 --taxon human`
- Python: `subprocess` wrapper; Python client `ncbi-datasets-pylib` (experimental as of 2024)

## Installation

```bash
# conda
conda install -c conda-forge ncbi-datasets-cli

# Or direct download (Linux, macOS, Windows binaries)
curl -O https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets

datasets --version    # 16.0+ expected
dataformat --version  # bundled companion tool
```

## What's in scope (use Datasets) vs out of scope (use E-utilities or other tools)

| Question | Datasets | Use instead |
|---|---|---|
| Genome assembly download | yes | — |
| All reference genomes for a taxon | yes | — |
| Gene record metadata (multi-species) | yes | — |
| Ortholog data for a gene | yes (`datasets summary gene ... --ortholog`) | OrthoDB / Compara for tree-aware orthology |
| Virus data (assemblies, metadata) | yes (`datasets download virus`) | — |
| Annotation files (GFF3, GTF) for a genome | yes | — |
| Protein records (curated, with cross-refs) | partial | UniProt REST for richer annotation |
| PubMed | no | `entrez-search` / `entrez-fetch` |
| SRA reads | no | `sra-data` |
| BLAST | no | `blast-searches` / `local-blast` |
| Custom Entrez queries | no | `entrez-search` |
| Pre-computed alignments (Compara) | no | `ensembl-rest` |

## Subcommand taxonomy

| Subcommand | Purpose | Example |
|---|---|---|
| `datasets summary genome` | Metadata only; JSON output | `datasets summary genome accession GCF_000001405.40` |
| `datasets download genome` | Download data files | `datasets download genome accession GCF_...` |
| `datasets summary gene` | Gene record metadata | `datasets summary gene symbol BRCA1 --taxon human` |
| `datasets download gene` | Download gene products | `datasets download gene symbol BRCA1 --taxon human` |
| `datasets summary taxonomy` | Taxonomy info | `datasets summary taxonomy taxon human` |
| `datasets download virus` | Virus assemblies/proteins | `datasets download virus genome taxon SARS-CoV-2` |
| `dataformat tsv` / `dataformat excel` | Convert JSON-lines to tabular | `dataformat tsv gene-summary` |

`datasets summary` always returns JSON-lines on stdout (one object per record). `datasets download` produces a `.zip` (default) or a "dehydrated" stub for cloud workflows.

## Key parameters (download)

| Flag | Effect |
|---|---|
| `--filename out.zip` | Where to write the archive |
| `--include genome,gff3,gtf,protein,cds,rna,seq-report` | Which file types to include |
| `--reference` | Restrict to reference assemblies only (one per species) |
| `--annotated` | Restrict to annotated assemblies |
| `--assembly-source RefSeq` / `GenBank` / `all` | Database source |
| `--assembly-level chromosome,complete` | Assembly quality level |
| `--released-after 2024-01-01` | Date filter |
| `--dehydrated` | Skip data; download just stubs + URL list (for parallel pull) |
| `--api-key XXX` | Optional API key (raises rate limit) |
| `--no-progressbar` | For non-interactive use |

For very large pulls (1000+ genomes), `--dehydrated` is the right choice: download the metadata stubs first, then run `datasets rehydrate` later or pull URLs in parallel from the manifest.

## JSON-lines output + dataformat

`datasets summary` returns JSON-lines (one JSON object per line) on stdout. Pipe through `dataformat tsv` for tabular:

```bash
datasets summary genome taxon "Escherichia coli" --reference --as-json-lines \
  | dataformat tsv genome --fields accession,organism-name,assembly-level,scaffold-n50 \
  > ecoli_refs.tsv
```

`dataformat` subcommands match summary types: `genome`, `gene`, `virus-genome`, etc. The `--fields` list is documented per type via `dataformat tsv <type> --help`.

## When to use --dehydrated for cloud workflows

The "dehydrated" mode separates data discovery from data transfer:

1. **Discover**: `datasets download genome taxon human --reference --dehydrated --filename human.zip` (fast; ~MB).
2. **Inspect**: `unzip -p human.zip ncbi_dataset/fetch.txt` -- a TSV of all URLs to pull.
3. **Pull**: either `datasets rehydrate --directory ./human/` or use `aria2c --input-file=fetch.txt` for parallel pull.

This is essential for HPC / cloud pipelines where inspection of the pending transfer is needed before committing the I/O.

## Checksum verification (automatic)

`datasets` verifies MD5 checksums for every downloaded file automatically. Rehydrate workflows also verify. If a file fails checksum, Datasets retries up to 3 times then errors. This replaces the `md5sum -c` step that was required with assembly_summary.txt-based scraping.

## Code patterns

### Download a single reference genome

**Goal:** Get human reference assembly with genome + GTF + protein + CDS.

**Approach:** `datasets download genome accession ... --include ...`.

**Reference (NCBI Datasets CLI 16.0+):**
```bash
#!/bin/bash
# Reference: NCBI Datasets CLI 16.0+ | Verify API if version differs

datasets download genome accession GCF_000001405.40 \
    --include genome,gff3,gtf,protein,cds,seq-report \
    --filename human_grch38.zip

unzip -q human_grch38.zip -d human_grch38/
ls -lh human_grch38/ncbi_dataset/data/GCF_000001405.40/
```

### Bulk download all reference bacterial genomes

**Goal:** Pull every RefSeq reference bacterial assembly with annotation.

**Approach:** `--dehydrated` first for inspection; rehydrate with parallel pull.

**Reference (NCBI Datasets CLI 16.0+):**
```bash
#!/bin/bash
# Step 1: dehydrated discovery
datasets download genome taxon Bacteria \
    --reference --annotated --assembly-source RefSeq \
    --include genome,gff3,protein \
    --dehydrated --filename bact_refs.zip

unzip -q bact_refs.zip -d bact_refs/
wc -l bact_refs/ncbi_dataset/fetch.txt   # how many files will be pulled

# Step 2: parallel pull via aria2 (or datasets rehydrate)
aria2c --input-file=bact_refs/ncbi_dataset/fetch.txt \
       --dir=bact_refs/ncbi_dataset/data/ \
       --max-concurrent-downloads=8 \
       --retry-wait=5
```

### Gene metadata across species

```bash
datasets summary gene symbol BRCA1 \
    --taxon Mammalia \
    --as-json-lines \
  | dataformat tsv gene --fields gene-id,symbol,taxname,description,nomenclature-authority,chromosomes \
  > brca1_mammals.tsv

head brca1_mammals.tsv
```

### Find orthologs for a gene

```bash
datasets summary gene symbol BRCA1 --taxon human --ortholog --as-json-lines \
  | dataformat tsv gene --fields gene-id,symbol,taxname,description \
  > brca1_orthologs.tsv
```

`--ortholog` returns NCBI's ortholog set (a single representative per species; tree-aware orthology with multiple co-orthologs is in `ortholog-inference` / Compara / OMA).

### Filter assemblies by quality and date

```bash
datasets summary genome taxon "Salmonella enterica" \
    --assembly-level chromosome,complete \
    --released-after 2024-01-01 \
    --as-json-lines \
  | dataformat tsv genome --fields accession,organism-name,assembly-level,scaffold-n50,submission-date \
  > sal_2024.tsv
```

### Python wrapper with checksum + retry awareness

**Reference (NCBI Datasets CLI 16.0+):**
```python
import subprocess
import json
from pathlib import Path


def datasets_summary(subcommand, *args):
    '''Run `datasets summary` and parse JSON-lines stdout.'''
    cmd = ['datasets', 'summary', subcommand, *args, '--as-json-lines']
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return [json.loads(line) for line in out.stdout.strip().split('\n') if line]


def datasets_download(subcommand, *args, out='dataset.zip', include=None):
    cmd = ['datasets', 'download', subcommand, *args, '--filename', out]
    if include:
        cmd += ['--include', ','.join(include)]
    subprocess.run(cmd, check=True)
    return Path(out)


genomes = datasets_summary('genome', 'taxon', 'Escherichia coli', '--reference')
print(f'{len(genomes)} reference E. coli assemblies')
for g in genomes[:3]:
    acc = g.get('accession')
    n50 = g.get('assemblyStats', {}).get('contigN50')
    print(f'  {acc}  N50={n50}')

datasets_download('genome', 'accession', 'GCF_000005845.2',
                  out='ecoli_k12.zip',
                  include=['genome', 'gff3', 'protein'])
```

### Comparison vs E-utilities

```python
# E-utilities path: ESearch in assembly db -> ESummary -> manual FTP pull
#   ~30 API calls + manual md5 + serial download
# Datasets path:
#   datasets download genome accession GCF_...  # one command, automatic md5, parallel inside
```

For genome workflows, Datasets is 5-50x faster than the equivalent E-utilities pipeline and far more reliable.

## Failure modes

### Choosing Datasets for the wrong question
- **Trigger:** Trying to pull raw SRA reads via Datasets.
- **Mechanism:** Datasets covers genome/gene/ortholog, not raw reads.
- **Symptom:** Subcommand not found or empty result.
- **Fix:** Use `sra-data` skill (prefetch/fasterq-dump) for raw reads.

### `--reference` filter loses too much
- **Trigger:** Bulk pull of "all assemblies for a species"; `--reference` returns one per species.
- **Mechanism:** Reference subset is the canonical single representative.
- **Symptom:** Far fewer assemblies than expected for a species with hundreds of submissions.
- **Fix:** Drop `--reference` for full set; add `--assembly-level chromosome,complete` for quality filter instead.

### Dehydrated workflow forgotten
- **Trigger:** 1000-genome pull without `--dehydrated`.
- **Mechanism:** Datasets downloads serially within one ZIP; can take hours.
- **Symptom:** Slow; no parallelism; one giant ZIP.
- **Fix:** Use `--dehydrated` + aria2c with `--max-concurrent-downloads`.

### dataformat field name guessing
- **Trigger:** `dataformat tsv genome --fields foo,bar` with invented field names.
- **Mechanism:** Field names are constrained per summary type.
- **Symptom:** "Unknown field" error.
- **Fix:** `dataformat tsv genome --help` lists valid field names; pull JSON-lines and inspect with `jq` to discover fields.

### Old assembly_summary.txt-based scripts still in use
- **Trigger:** Legacy pipeline scraping `https://ftp.ncbi.nlm.nih.gov/genomes/all/refseq/...`.
- **Mechanism:** Pre-2023 best practice; FTP listing parsing is fragile.
- **Symptom:** Slow; brittle; no checksums; broken when NCBI restructures FTP.
- **Fix:** Switch to Datasets CLI; the FTP path still works but Datasets is the supported modern path.

### API key not used for high-volume
- **Trigger:** 1000+ summary calls in a loop without `--api-key`.
- **Mechanism:** NCBI rate-limits unauthenticated bulk traffic.
- **Symptom:** Throttling; slow downloads.
- **Fix:** Pass `--api-key YOUR_KEY` to bulk commands; obtain from `https://www.ncbi.nlm.nih.gov/account/settings/`.

### CLI version drift
- **Trigger:** Using Datasets v14 with v16 docs.
- **Mechanism:** Subcommands and flags renamed between major versions.
- **Symptom:** "Unknown flag" or different output structure.
- **Fix:** Pin to v16+; `conda update ncbi-datasets-cli`.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| "command not found: datasets" | Not installed | `conda install -c conda-forge ncbi-datasets-cli` |
| Subcommand not found | Old version | Upgrade to v16+ |
| Slow 1000-genome pull | Serial download | Use `--dehydrated` + aria2c |
| "Unknown field" in dataformat | Wrong field name | Check `dataformat <type> --help` |
| Throttled bulk pull | No API key | Pass `--api-key` |
| `--reference` returns 1 per species | By design | Drop the flag or use `--assembly-level` |
| MD5 mismatch retried | Network issue | Datasets retries automatically; persistent failure -> investigate network |

## References

- O'Leary NA, Cox E, Holmes JB, et al. (2024) Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. *Sci Data* 11:732.
- NCBI Datasets documentation: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
- NCBI. Datasets CLI usage. https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/

## Related Skills

- entrez-search - For PubMed, custom queries, and non-genome data
- entrez-fetch - For single-record fetches outside genome/gene scope
- batch-downloads - Bulk E-utilities (when not genome-scale)
- sra-data - Raw sequencing reads (NOT covered by Datasets)
- ensembl-rest - Ensembl REST as alternative for Ensembl-native species
- ortholog-inference - Compara/OMA/OrthoDB for tree-aware orthology
