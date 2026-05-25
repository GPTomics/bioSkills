---
name: bio-sra-data
description: Download raw sequencing reads from NCBI SRA using sra-tools (prefetch, fasterq-dump, vdb-validate) or the ENA mirror. Use when pulling FASTQ for SRR/ERR/DRR accessions, deciding between SRA-direct, ENA mirror, or AWS/GCP cloud mirror (STRIDES), handling --include-technical for 10x and other single-cell records, validating with MD5/vdb-validate, navigating SRR/SRX/SRS/SRP/PRJNA hierarchy, or finding accessions via pysradb. Encodes SRA cloud-egress economics, the fasterq-dump uncompressed-scratch trap, and the --max-size default that silently truncates large prefetches.
tool_type: cli
primary_tool: sra-tools
---

## Version Compatibility

Reference examples tested with: sra-tools 3.0+ (fasterq-dump, prefetch, vdb-validate, vdb-config), pysradb 2.2+, ENA portal API 2.0+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `fasterq-dump --version`, `prefetch --version`
- Python: `pip show pysradb`

If a flag is unrecognized or behavior changes, run `<tool> --help` and adapt.

# SRA Data

**"Download FASTQ from this SRA accession"** -> Two paths exist in 2026: the **SRA toolkit** (NCBI's official, with prefetch + fasterq-dump) and the **ENA mirror** (EMBL-EBI's mirror with direct FASTQ download, often faster). For >1 TB workflows, a third path: **AWS Open Data** (STRIDES program) where same-region EC2 pulls SRA data with zero egress cost.

The single most impactful decision is **where to pull from**. SRA-direct is the default but ENA is faster more often than not, and AWS Open Data is the right answer for cloud-native analysis pipelines.

- CLI: `prefetch SRR...`, `fasterq-dump SRR...`, `vdb-validate SRR...` (sra-tools)
- CLI: `curl https://ftp.sra.ebi.ac.uk/...` (ENA mirror; direct FASTQ)
- CLI: `aws s3 cp s3://sra-pub-run-odp/sra/SRR.../SRR... ./SRR....sra ...` (STRIDES; object is unsuffixed; same-region free)
- Python: `pysradb` for metadata; `subprocess` for download

## Required Setup

```bash
# sra-tools (toolkit)
conda install -c bioconda sra-tools           # 3.0+
fasterq-dump --version                        # confirm

# Configure cache location (default ~/ncbi/ -- often too small)
vdb-config --cfg                              # show current config
vdb-config --set /repository/user/main/public/root=/data/sra_cache

# Optional: pysradb for metadata
pip install pysradb
```

For STRIDES cloud:
```bash
# AWS CLI (no NCBI auth needed for public buckets)
aws s3 ls s3://sra-pub-run-odp/sra/SRR12345678/ --no-sign-request
```

## Decision matrix: where to pull from

| Source | When best | Speed | Cost |
|---|---|---|---|
| **ENA mirror** (FTP/Aspera) | Default for most workflows | Often fastest; direct FASTQ (no SRA->FASTQ conversion needed) | Free; no rate limit observed |
| **SRA toolkit + AWS STRIDES** | Same-region EC2/EKS | Fastest within AWS us-east-1 | Free egress within region; small storage cost |
| **SRA toolkit + GCP STRIDES** | Same-region GCP Compute Engine | Fastest within GCP us-central1 | Free egress within region |
| **SRA-direct (prefetch + fasterq-dump)** | On-prem; small downloads; need SRA-format access | Variable; can be slow off-peak fails | Free; NCBI throttles by IP |
| **Aspera (`ascp`)** | Institutional accounts only | Faster than HTTPS on long links | NCBI public Aspera retired 2019; ENA public Aspera retired ~2023; institutional use still possible |

**Default recommendation**: **ENA mirror** for off-cloud, **STRIDES (AWS/GCP)** for in-cloud analysis. SRA-direct only when neither is available or when SRA format itself is needed (e.g. for re-extraction of technical reads).

## SRA accession hierarchy

| Prefix | Type | Granularity |
|---|---|---|
| SRR / ERR / DRR | Run | One sequencing run (file-level) |
| SRX / ERX / DRX | Experiment | Library prep + sequencing strategy |
| SRS / ERS / DRS | Sample | Biological sample |
| SRP / ERP / DRP | Study | Project (deprecated; superseded by BioProject) |
| PRJNA / PRJEB / PRJDB | BioProject | Top-level project ID |
| SAMN / SAMEA / SAMD | BioSample | Biological sample (cross-archive) |

Conversion is via SRA metadata: `pysradb metadata <ID>` or `efetch -db sra -id <UID> -rettype runinfo`.

The actual download unit is SRR/ERR/DRR (runs). The BioProject (PRJNA...) is the convenient top-level handle for "pull all data for paper X".

## fasterq-dump vs fastq-dump

`fasterq-dump` (sra-tools 2.10+) is the multi-threaded successor. **Always prefer it**, with two exceptions noted below.

| Aspect | fasterq-dump | fastq-dump |
|---|---|---|
| Threads | Multi (`-e N`) | Single |
| Speed | ~5-10x faster | Baseline |
| Disk overhead | Writes uncompressed FASTQ to scratch (~3x final size) | In-place; lower scratch |
| Compression | NOT built-in (post-process with pigz) | `--gzip` flag built-in |
| Single-cell technical reads | `--include-technical` works | Some 10x records need fastq-dump for full extraction |
| 10x split semantics | Sometimes incomplete | Sometimes the only way to get all reads |

The **uncompressed-scratch trap**: `fasterq-dump` writes uncompressed FASTQ first, then leaves it uncompressed. A 100 GB compressed FASTQ needs ~300 GB of scratch space + 300 GB of final output. Either compress post-hoc with `pigz` or use `--mem` to control RAM/disk tradeoff.

## prefetch and the `--max-size` trap

`prefetch` downloads `.sra` files to the configured cache before extraction. Default `--max-size 20G` silently skips runs larger than 20 GB.

```bash
# Wrong: silently skips runs >20 GB
prefetch SRR12345678

# Right: set max-size explicitly to your largest expected size
prefetch SRR12345678 --max-size 100G -p
```

For unknown-size queues, set max-size to a generous upper bound (e.g. `--max-size 200G`) or query metadata first with `pysradb metadata`.

## ENA mirror: direct FASTQ URLs

ENA stores FASTQ files directly (no SRA-format intermediate). Discover URLs via the ENA portal API:

```bash
curl 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR12345678&result=read_run&fields=fastq_ftp,fastq_md5,read_count&format=tsv'
```

Returns TSV with semicolon-separated paired-end URLs and md5 checksums.

Direct download:
```bash
curl -O 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/078/SRR12345678/SRR12345678_1.fastq.gz'
```

ENA's mirror is typically faster than SRA's because (a) it's hosted on Aspera-aware servers, (b) the FASTQ is pre-compressed (no SRA->FASTQ conversion needed), (c) EMBL-EBI's bandwidth is generous. For most downloads in 2026, ENA is the right default.

## Single-cell / 10x quirks

10x Genomics records include "technical reads" (cell barcodes, UMIs) interleaved with biological reads. Default `fasterq-dump` (or `fastq-dump`) skips them. To get all reads:

```bash
# fasterq-dump with technical reads
fasterq-dump SRR12345678 --include-technical --split-files -p -O ./fastq/

# Some 10x records require fastq-dump -- check sra-stat first
sra-stat --xml SRR12345678 | grep -E '(spotCount|baseCount|tag)'
```

For 10x v3, expect 3 files per run: R1 (barcode+UMI), R2 (cDNA), I1 (index). For 10x v2: R1 (barcode), R2 (UMI+cDNA), I1.

## MD5 / vdb-validate

Always verify downloads.

```bash
# vdb-validate for SRA-format files (toolkit path)
vdb-validate SRR12345678

# md5sum for ENA FASTQ files
md5sum -c <(echo "<expected_md5>  SRR12345678_1.fastq.gz")
```

ENA provides md5 in the portal API response. SRA-toolkit's `vdb-validate` is the equivalent for `.sra` files (different file format).

## Cloud (STRIDES) access

NCBI's STRIDES initiative mirrored SRA data to AWS Open Data (us-east-1) and GCP (us-central1). Same-region pulls have zero egress cost.

```bash
# List SRA cloud-hosted files (no NCBI auth needed)
aws s3 ls s3://sra-pub-run-odp/sra/SRR12345678/ --no-sign-request

# Direct copy to EC2 in us-east-1. The STRIDES object is named without a `.sra`
# suffix (just SRR12345678); rename on copy to keep fasterq-dump happy.
aws s3 cp s3://sra-pub-run-odp/sra/SRR12345678/SRR12345678 ./SRR12345678.sra --no-sign-request

# Then fasterq-dump locally
fasterq-dump ./SRR12345678.sra -p -e 8
```

For cloud-native analysis pipelines (Nextflow on AWS Batch, Cromwell, etc.), STRIDES is the right path.

## Code patterns

### Single SRR via ENA mirror (preferred default)

**Goal:** Download paired-end FASTQ for one SRR; verify md5; minimal dependencies.

**Approach:** Query ENA portal API for FASTQ URLs and md5; download with curl; verify with md5sum.

**Reference (ENA portal API 2.0+, curl):**
```bash
#!/bin/bash
SRR="${1:-SRR12345678}"
OUT="${2:-./fastq}"
mkdir -p "${OUT}"

# Get FASTQ URLs + md5 from ENA portal API
META=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -1)
URLS=$(echo "${META}" | cut -f1 | tr ';' '\n')
MD5S=$(echo "${META}" | cut -f2 | tr ';' '\n')

i=0
while read url; do
    fname="${OUT}/$(basename ${url})"
    expected_md5=$(echo "${MD5S}" | sed -n "$((i+1))p")
    echo "Downloading ${fname}"
    curl -sL -o "${fname}" "https://${url}"
    actual_md5=$(md5sum "${fname}" | awk '{print $1}')
    if [ "${actual_md5}" != "${expected_md5}" ]; then
        echo "MD5 MISMATCH ${fname}: expected ${expected_md5}, got ${actual_md5}"
        exit 1
    fi
    echo "  md5 OK"
    i=$((i+1))
done <<< "${URLS}"
```

### prefetch + fasterq-dump (SRA toolkit, classic)

```bash
#!/bin/bash
SRR="${1:-SRR12345678}"
OUT="${2:-./fastq}"
THREADS="${3:-8}"
mkdir -p "${OUT}"

# prefetch with explicit max-size (default 20G silently skips larger)
prefetch "${SRR}" --max-size 100G -p

# Validate SRA file
vdb-validate "${SRR}" || { echo "Validation FAILED"; exit 1; }

# Extract FASTQ (multi-threaded; uncompressed scratch ~3x final size)
fasterq-dump "${SRR}" -O "${OUT}" -e "${THREADS}" -p --split-files

# Compress post-hoc (fasterq-dump does NOT compress)
pigz -p "${THREADS}" "${OUT}/${SRR}"_*.fastq

# Cleanup SRA cache if you don't need it
# rm -rf ~/ncbi/sra/${SRR}.sra
```

### Batch via pysradb metadata

**Goal:** Convert a list of GSE / BioProject / SRX IDs to SRR run accessions.

**Approach:** pysradb metadata returns a full hierarchy table; pull SRR column.

**Reference (pysradb 2.2+):**
```python
from pysradb import SRAweb
import pandas as pd


def gse_to_srr(gse):
    db = SRAweb()
    df = db.gse_to_srp(gse)
    if df.empty:
        return []
    srp = df['study_accession'].iloc[0]
    runs = db.srp_to_srr(srp)
    return runs['run_accession'].tolist()


def bioproject_to_runs(prjna):
    db = SRAweb()
    return db.sra_metadata(prjna, detailed=True)


def batch_resolve(ids):
    db = SRAweb()
    rows = []
    for id in ids:
        try:
            meta = db.sra_metadata(id, detailed=True)
            rows.append(meta)
        except Exception as e:
            print(f'{id}: {e}')
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# Resolve a GSE to all its SRRs
srrs = gse_to_srr('GSE123456')
print(f'GSE123456 -> {len(srrs)} SRRs')
```

### Cloud (STRIDES) via AWS

```bash
#!/bin/bash
# Run from EC2 in us-east-1 for zero egress
SRR="${1:-SRR12345678}"

# Check if available on AWS Open Data
aws s3 ls "s3://sra-pub-run-odp/sra/${SRR}/" --no-sign-request

# Download .sra (then extract locally)
aws s3 cp "s3://sra-pub-run-odp/sra/${SRR}/${SRR}" "./${SRR}.sra" --no-sign-request

fasterq-dump "./${SRR}.sra" -p -e 8 --split-files
pigz -p 8 "${SRR}"_*.fastq
```

### 10x single-cell with technical reads

```bash
#!/bin/bash
SRR="${1:-SRR_10x_run}"
OUT="${2:-./fastq_10x}"
mkdir -p "${OUT}"

# Get all reads including technical (barcode/UMI/index)
fasterq-dump "${SRR}" --include-technical --split-files -p -O "${OUT}" -e 8

# 10x v3 expects: R1 (28-bp barcode+UMI), R2 (cDNA), I1 (sample index)
ls -la "${OUT}/${SRR}"_*.fastq
pigz -p 8 "${OUT}/${SRR}"_*.fastq
```

## Failure modes

### prefetch --max-size silent skip
- **Trigger:** Default 20 GB limit; run is 50 GB.
- **Mechanism:** prefetch returns success but downloads nothing.
- **Symptom:** vdb-validate or fasterq-dump fails because no file exists.
- **Fix:** Always set `--max-size` explicitly to a generous upper bound (e.g. 200G).

### fasterq-dump scratch space exhaustion
- **Trigger:** Run is 100 GB compressed; scratch dir has 200 GB free.
- **Mechanism:** fasterq-dump writes ~300 GB uncompressed, fills disk.
- **Symptom:** "out of disk space" mid-extraction.
- **Fix:** Use a scratch dir with 4-5x the compressed size; or use `--mem` to trade memory for disk; or stick with `fastq-dump --gzip` (slower but lower scratch).

### 10x technical reads missing
- **Trigger:** Default `fasterq-dump` on a 10x record.
- **Mechanism:** Technical reads (barcodes, UMIs) are skipped by default.
- **Symptom:** Only the cDNA file (R2) appears; CellRanger / STARsolo errors.
- **Fix:** Add `--include-technical`; verify with `sra-stat --xml` first.

### SRA-direct slowness during US business hours
- **Trigger:** Downloading from NCBI 9 AM-5 PM ET weekdays.
- **Mechanism:** NCBI bandwidth contention; institutional users have priority.
- **Symptom:** kbps-level download speeds.
- **Fix:** Switch to ENA mirror or AWS STRIDES; run outside US business hours.

### Aspera deprecation
- **Trigger:** Old script using `ascp` against `anonftp@ftp.ncbi.nlm.nih.gov`.
- **Mechanism:** NCBI retired public Aspera in 2019; ENA followed ~2023; only institutional accounts retain support.
- **Symptom:** Connection refused or auth fails.
- **Fix:** Switch to HTTPS (slower but works); for fastest cloud transfer use STRIDES (AWS/GCP).

### Cloud egress costs surprise
- **Trigger:** STRIDES pull from EC2 in us-west-2 against bucket in us-east-1.
- **Mechanism:** Cross-region egress is charged.
- **Symptom:** Unexpected AWS bill.
- **Fix:** Match compute region to bucket region (us-east-1 for AWS, us-central1 for GCP).

### vdb-config not persisted across containers
- **Trigger:** Docker container without persisted `~/.ncbi/user-settings.mkfg`.
- **Mechanism:** Cache config is per-user, per-home; container rebuild loses it.
- **Symptom:** Cache fills container's small layer; download fails.
- **Fix:** Mount a host volume at `~/.ncbi/` and persist user-settings.mkfg; or set `--temp` and `-O` explicitly in commands.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| "item not found" | Invalid accession or not in current SRA | Verify; check ENA mirror |
| Scratch disk full mid-extraction | fasterq-dump uncompressed write | Use larger scratch or fastq-dump --gzip |
| Slow SRA-direct download | Business-hours contention | ENA or STRIDES |
| 10x reads missing | --include-technical not set | Add the flag |
| Container loses cache config | vdb-config not persisted | Mount ~/.ncbi as volume |
| prefetch returns "success" but no file | --max-size silent skip | Set --max-size explicitly |
| AWS bill on STRIDES | Cross-region pull | Match compute region |

## References

- NCBI. SRA Toolkit documentation. https://github.com/ncbi/sra-tools/wiki
- NCBI. STRIDES program. https://datascience.nih.gov/strides
- Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. (2011) The sequence read archive. *Nucleic Acids Res* 39:D19-D21.
- Cochrane G, Karsch-Mizrachi I, Takagi T; International Nucleotide Sequence Database Collaboration. (2016) The International Nucleotide Sequence Database Collaboration. *Nucleic Acids Res* 44:D48-D50.
- Choudhary S. (2019) pysradb: A Python package to query next-generation sequencing metadata and data from NCBI Sequence Read Archive. *F1000Research* 8:532.

## Related Skills

- entrez-search - Search the SRA db for accessions before downloading
- geo-data - GEO Series often link to SRA; gds -> sra ELink
- read-qc/quality-reports - QC the downloaded FASTQ
- read-qc/fastp-workflow - Adapter trim downloaded FASTQ
- ncbi-datasets-cli - Modern bulk path for genome data (NOT for SRA reads)
