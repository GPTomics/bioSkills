---
name: bio-geo-data
description: Query and download from NCBI Gene Expression Omnibus (GEO) and EMBL-EBI's BioStudies/ArrayExpress mirror. Use when finding expression datasets, navigating SuperSeries vs SubSeries, choosing between series-matrix (submitter-normalized) and raw supplementary files, downloading via GEOparse (Python) or GEOquery (R/Bioconductor), linking GEO to SRA for raw reads, or distinguishing GSE/GSM/GPL/GDS record types. Encodes the SuperSeries trap, the series-matrix normalization-trust caveat, GEOmetadb deprecation, ArrayExpress migration to BioStudies, and processed-vs-raw decision matrix.
tool_type: mixed
primary_tool: Bio.Entrez
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, GEOparse 2.0+, R Bioconductor GEOquery 2.70+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython geoparse` then introspect signatures
- R: `packageVersion('GEOquery')`

If the GSE structure doesn't match expectations (missing fields, malformed series matrix), re-fetch from FTP directly and inspect the SOFT or MINiML file as source of truth.

# GEO Data

**"Pull expression data from GEO accession GSE..."** -> GEO stores Series (GSE), Samples (GSM), Platforms (GPL), and curated DataSets (GDS, frozen 2018). The single most consequential decision is **processed (series matrix) vs raw (supplementary files / linked SRA)** — the answer turns on how much trust the submitter's normalization deserves.

The single most-missed gotcha: **SuperSeries**. A GSE may be a meta-container (`!Series_relation = SuperSeries of: GSExxxxx`) holding multiple sub-studies on different platforms. Naively pulling samples from a SuperSeries gives mixed Affymetrix + Illumina + RNA-seq, mis-batched.

- Python: `Entrez.esearch(db='gds')`, GEOparse for full series download
- R: `GEOquery::getGEO()` (Bioconductor; more mature than GEOparse)
- CLI: `wget` from `ftp.ncbi.nlm.nih.gov/geo/series/...`

## Required Setup

```bash
pip install biopython GEOparse pandas
# OR for R-side:
# R: BiocManager::install('GEOquery')
```

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional'
```

## GEO record taxonomy

| Prefix | Type | Granularity | What's in it |
|---|---|---|---|
| GSE | Series | One study | Title, summary, design, links to GSMs, supplementary files |
| GSM | Sample | One biological/technical sample | Submitter metadata, per-sample processed data, link to raw SRA |
| GPL | Platform | One array / sequencer | Probe annotations or sequencer model |
| GDS | DataSet | Curated, normalized subset of one GSE | Re-normalized expression matrix (frozen 2018; new GDS no longer created) |
| GSEXXX SuperSeries | Series meta-container | Wraps multiple SubSeries | `!Series_relation = SuperSeries of: ...` |

**`GDS` is dead-as-format**: NCBI stopped creating new GDS records in 2018. Existing GDS still queryable but use GSE for anything current.

## The SuperSeries trap

A SuperSeries (GSE) wraps multiple SubSeries, often with different platforms. Detection:

```python
# Read the !Series_relation field from SOFT format
from Bio import Entrez
h = Entrez.esummary(db='gds', id='200122288')   # example
r = Entrez.read(h)[0]; h.close()
print(r.get('summary'))   # may or may not flag SuperSeries
# Definitive check: download SOFT and grep:
#   curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122288/soft/GSE122288_family.soft.gz | zgrep Series_relation
```

A `SuperSeries of: GSE12345` line means the SuperSeries' samples are the union of all SubSeries — almost certainly mixed-platform / mixed-batch. Process each SubSeries independently.

Symmetric trap: a paper may cite a SubSeries (`SubSeries of: GSEsuper`) where the wider context is essential — check both directions.

## Decision matrix: processed vs raw vs SRA

| Question | Source | Trust level |
|---|---|---|
| "I want expression values; submitter normalization is fine" | Series matrix (`GSE_series_matrix.txt.gz`) | Trust submitter's normalization |
| "I want raw Affymetrix CEL files and to do my own RMA" | Supplementary files (`suppl/`) | Re-normalize locally |
| "I want raw RNA-seq FASTQ" | pysradb `gse_to_srp -> srp_to_srr` (Entrez gds->sra ELink unreliable) | Always raw; processed at submitter is rarely re-usable |
| "I want submitter-provided counts (RNA-seq)" | Supplementary files (usually a `*_counts.txt.gz`) | Trust at risk; submitter pipelines vary |
| "I want a curated subset across many studies" | Use ArchS4 (https://archs4.org) or recount3 | Curated re-processing |

**Default to raw whenever possible.** For Affymetrix: CEL + locally-run RMA is far more reliable than the submitter's "normalized" matrix. For RNA-seq: SRA FASTQ + locally-run alignment/quantification is the only reproducible path; submitter counts often use a private pipeline.

## Series matrix files

A series matrix (`GSE12345_series_matrix.txt.gz`) is a header (sample metadata as `!Sample_*` lines) plus a sample-by-feature expression table. The format is fragile and the values' provenance is whatever the submitter chose. Critical caveats:

- For Affymetrix: the matrix is usually RMA-normalized but submitters sometimes apply additional transforms (log2, scaling, batch correction).
- For RNA-seq: the matrix is sometimes log-CPM, sometimes raw counts, sometimes VST/rlog — read `!Series_overall_design` and `!Sample_data_processing` to know.
- The header has `!Sample_characteristics_ch1` rows that hold the metadata of interest — these are submitter-formatted strings, often inconsistent within one series.

## SOFT vs MINiML

| Format | Content | Parser support |
|---|---|---|
| **SOFT** (`*_family.soft.gz`) | Plain-text, key=value style | GEOparse (Python), GEOquery (R), Entrez Direct |
| **MINiML** (`*_family.xml.tgz`) | XML-structured | GEOparse, GEOquery, custom XML |

Both contain the same content. SOFT is the legacy, MINiML the XML successor. GEOparse handles SOFT well; for very large series (1000+ samples) MINiML's XML structure is slower to parse.

## GEOparse vs GEOquery

| Aspect | GEOparse (Python) | GEOquery (R/Bioconductor) |
|---|---|---|
| Maturity | OK; some known supplementary-file fetch issues since ~2022 | Mature; Bioconductor-supported |
| Output | `GEOparse.GSE` object with `gsms`, `gpls`, `metadata` dicts | `ExpressionSet` or list per platform |
| Supplementary files | `gse.download_supplementary_files()` (sometimes flakey) | `getGEOSuppFiles(gse)` (more reliable) |
| Integration | Pandas DataFrames | Bioconductor ecosystem |
| When | Python-first pipelines | R-first / use ExpressionSet downstream |

For production GEO workflows in R, GEOquery is the stable choice. For Python, GEOparse is the only option but verify file counts after download.

## GEOmetadb status

GEOmetadb (Zhu 2008) was a SQLite mirror of GEO metadata enabling fast SQL queries. **Unmaintained since 2020**; downloads still work but data is stale. Modern replacement: pysradb (`pysradb gse_to_srp`, `pysradb metadata`) covers most of the GEO->SRA mapping; for full GEO queries fall back to Entrez gds.

## ArrayExpress -> BioStudies migration (2020)

ArrayExpress (EMBL-EBI's microarray archive, mirroring GEO) was migrated into BioStudies in 2020. Old `E-MTAB-####` accessions still resolve but the API moved:

| Old (pre-2020) | New (BioStudies) |
|---|---|
| `https://www.ebi.ac.uk/arrayexpress/...` | `https://www.ebi.ac.uk/biostudies/...` |
| ArrayExpress REST | BioStudies REST: `https://www.ebi.ac.uk/biostudies/api/v1/...` |

For new workflows, use BioStudies. For legacy ArrayExpress URLs in old papers, redirect via BioStudies.

## Code patterns

### Search GEO for studies matching a query

**Goal:** Find GSE accessions matching keywords + organism + study type.

**Approach:** ESearch on `gds` db with field-qualified terms; filter to `gse[Entry Type]`; summarize with ESummary.

**Reference (BioPython 1.83+):**
```python
from Bio import Entrez
import time

Entrez.email = 'researcher@institution.edu'


def search_geo(term, study_type='gse', organism=None, max_results=50):
    full_term = f'{term} AND {study_type}[Entry Type]'
    if organism:
        full_term += f' AND {organism}[Organism]'
    h = Entrez.esearch(db='gds', term=full_term, retmax=max_results)
    s = Entrez.read(h); h.close()
    if not s['IdList']:
        return []
    h = Entrez.esummary(db='gds', id=','.join(s['IdList']))
    summaries = Entrez.read(h); h.close()
    return summaries


for s in search_geo('breast cancer RNA-seq', organism='Homo sapiens', max_results=10):
    # Surface SuperSeries
    relation = s.get('summary', '')
    is_super = 'SuperSeries' in str(relation)
    print(f"  {s['Accession']:12} {s['n_samples']:>4} samples  {'[SuperSeries]' if is_super else '':12}  {s['title'][:60]}")
```

### Detect SuperSeries before pulling data

**Goal:** Avoid mixing platforms by detecting SuperSeries structure first.

**Approach:** Download SOFT family file and read `!Series_relation` keys.

```python
import gzip
import urllib.request


def check_super_or_sub_series(gse):
    prefix = gse[:-3] + 'nnn'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/soft/{gse}_family.soft.gz'
    urllib.request.urlretrieve(url, f'{gse}.soft.gz')
    super_of = []
    sub_of = None
    with gzip.open(f'{gse}.soft.gz', 'rt') as f:
        for line in f:
            if line.startswith('!Series_relation'):
                if 'SuperSeries of' in line:
                    super_of.append(line.split('SuperSeries of: ')[1].strip())
                elif 'SubSeries of' in line:
                    sub_of = line.split('SubSeries of: ')[1].strip()
            if line.startswith('^SAMPLE'):
                break   # Speed: don't read past header
    return {'super_of': super_of, 'sub_of': sub_of}


print(check_super_or_sub_series('GSE122288'))
# {'super_of': ['GSExxxxx', 'GSEyyyyy'], 'sub_of': None}  -> SuperSeries; process subseries separately
```

### Download series matrix with submitter caveat

```python
import gzip
import pandas as pd


def download_series_matrix(gse):
    prefix = gse[:-3] + 'nnn'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz'
    urllib.request.urlretrieve(url, f'{gse}_matrix.txt.gz')
    return f'{gse}_matrix.txt.gz'


def parse_series_matrix(path):
    metadata = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!'):
                key, *vals = line.rstrip('\n').split('\t')
                metadata[key] = [v.strip('"') for v in vals]
        expr = pd.read_csv(f, sep='\t', index_col=0, comment='!')
    # Series matrix values are whatever submitter chose -- check metadata['!Sample_data_processing']
    return metadata, expr


meta, expr = parse_series_matrix(download_series_matrix('GSE123456'))
print('Sample-level data processing notes:')
for note in set(meta.get('!Sample_data_processing', [])):
    print(f'  - {note}')
```

### Link GEO Series to SRA runs (preferred path: pysradb)

```python
from pysradb import SRAweb


def gse_to_srr(gse):
    db = SRAweb()
    srp_df = db.gse_to_srp(gse)
    if srp_df.empty:
        return []
    srp = srp_df['study_accession'].iloc[0]
    srr_df = db.srp_to_srr(srp)
    return srr_df['run_accession'].tolist()


srrs = gse_to_srr('GSE123456')
print(f'GSE123456 -> {len(srrs)} SRR runs')
```

### GEOparse: full Series download

```python
import GEOparse


def get_gse(gse_id, dest='./geo_cache'):
    gse = GEOparse.get_GEO(geo=gse_id, destdir=dest)
    print(f'{gse_id}: {len(gse.gsms)} samples, {len(gse.gpls)} platforms')
    for gsm_name, gsm in list(gse.gsms.items())[:3]:
        print(f'  {gsm_name}: {gsm.metadata.get("title", ["?"])[0]}')
    return gse


# Supplementary files (raw data) -- verify file count manually after
gse = get_gse('GSE123456')
gse.download_supplementary_files(directory='./geo_cache')
```

### R: GEOquery (more reliable supplementary download)

```r
# Reference: Bioconductor GEOquery 2.70+ | Verify API if version differs
library(GEOquery)

gse <- getGEO('GSE123456', GSEMatrix = TRUE)
length(gse)             # one ExpressionSet per platform
head(pData(gse[[1]]))   # sample metadata
head(exprs(gse[[1]]))   # expression matrix (submitter-normalized -- verify processing notes)

# Raw / supplementary files
supp_dir <- getGEOSuppFiles('GSE123456', baseDir = './geo_cache')
list.files(rownames(supp_dir))
```

### Find datasets by PubMed citation

```python
def geo_from_pubmed(pmid):
    h = Entrez.elink(dbfrom='pubmed', db='gds', id=pmid)
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    gds_ids = [l['Id'] for l in r[0]['LinkSetDb'][0]['Link']]
    h = Entrez.esummary(db='gds', id=','.join(gds_ids))
    summaries = Entrez.read(h); h.close()
    return summaries
```

## Failure modes

### SuperSeries pulled as one experiment
- **Trigger:** GSE accession from a paper; turns out to be a SuperSeries wrapping multiple platforms.
- **Mechanism:** Default download merges all samples without flagging the structure.
- **Symptom:** Downstream batch correction can't recover the mixed-platform structure; spurious "batch" effects.
- **Fix:** Always check `!Series_relation` in SOFT before pulling; process SubSeries independently.

### Series matrix is not what it appears to be
- **Trigger:** Series matrix downloaded; treated as RMA-normalized when submitter applied additional transforms.
- **Mechanism:** Series matrix contents are at submitter's discretion.
- **Symptom:** Re-analysis gives different answers than the published paper.
- **Fix:** Read `!Sample_data_processing` to know what's in the matrix; re-normalize from raw if in doubt.

### Submitter-provided RNA-seq counts mis-trusted
- **Trigger:** Using a `*_counts.txt.gz` supplementary file as the count matrix.
- **Mechanism:** Submitter's pipeline (aligner, GTF version, counting strategy) is rarely documented.
- **Symptom:** Counts don't agree with re-quantification from SRA FASTQ.
- **Fix:** Pull SRA FASTQ + re-quantify with a known pipeline (Salmon, kallisto, STAR + featureCounts).

### Platform GPL mismatch
- **Trigger:** One GSE with multiple platforms; series matrix split across multiple files.
- **Mechanism:** `GSE_series_matrix.txt.gz` is the merged one; per-platform are `GSE-GPLxxx_series_matrix.txt.gz`.
- **Symptom:** "Missing samples" or NaN-heavy expression matrix.
- **Fix:** Download per-platform matrix files; check `!Series_platform_id` count.

### GEOparse supplementary files flakey
- **Trigger:** `gse.download_supplementary_files()` silently misses files.
- **Mechanism:** Known issue with the GEOparse FTP enumeration since ~2022.
- **Symptom:** Local cache missing CEL or counts files.
- **Fix:** Use R GEOquery or direct FTP `wget -r` on the suppl/ subdirectory.

### ArrayExpress URL rot
- **Trigger:** Old paper links `https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1234/`.
- **Mechanism:** ArrayExpress migrated to BioStudies in 2020.
- **Symptom:** 404 or redirect.
- **Fix:** Use `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1234`.

### GEOmetadb stale
- **Trigger:** Old pipeline downloads `GEOmetadb.sqlite` for fast queries.
- **Mechanism:** GEOmetadb unmaintained since 2020.
- **Symptom:** Missing recent series; outdated annotations.
- **Fix:** Switch to pysradb for SRA-linked queries; Entrez gds for full GEO.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Empty IdList for `gse[entry_type]` | Wrong field name | Use `gse[Entry Type]` (case-sensitive) |
| Matrix file has no expression data | SuperSeries with no aggregate matrix | Pull per-SubSeries matrices |
| Submitter "normalized" matrix gives different result than paper | Hidden submitter transforms | Re-process from raw |
| 404 on ArrayExpress URL | Migrated to BioStudies | Use new BioStudies URL |
| GEOparse missing CEL files | Known flake | Use R GEOquery or direct FTP |
| GEOmetadb-based pipeline missing recent series | DB unmaintained | Switch to pysradb / Entrez |

## References

- Edgar R, Domrachev M, Lash AE. (2002) Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Res* 30:207-210.
- Barrett T, Wilhite SE, Ledoux P, et al. (2013) NCBI GEO: archive for functional genomics data sets - update. *Nucleic Acids Res* 41:D991-D995.
- Davis S, Meltzer PS. (2007) GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics* 23:1846-1847.
- Gumienny R. GEOparse: Python library to parse GEO databases. https://github.com/guma44/GEOparse (no journal publication).
- Sarkans U, Gostev M, Athar A, et al. (2018) The BioStudies database--one stop shop for all data supporting a life sciences study. *Nucleic Acids Res* 46:D1266-D1270.
- Lachmann A, Torre D, Keenan AB, et al. (2018) Massive mining of publicly available RNA-seq data from human and mouse. *Nat Commun* 9:1366. (ARCHS4)
- Wilks C, Zheng SC, Chen FY, et al. (2021) recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol* 22:323.

## Related Skills

- entrez-search - General gds search
- entrez-link - gds <-> pubmed, bioproject links (gds->sra ELink is unreliable; use pysradb)
- sra-data - Download raw FASTQ from GEO-linked SRA runs
- expression-matrix/normalization - Re-normalize raw expression data
- rna-quantification/alignment-free-quant - Salmon/kallisto re-quantification of GEO/SRA data
- ensembl-rest - Cross-reference Ensembl IDs in series-matrix files
