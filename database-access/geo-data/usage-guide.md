# GEO Data Usage Guide

## Overview

Query NCBI GEO (and EMBL-EBI's BioStudies/ArrayExpress mirror) for expression datasets. Encodes the SuperSeries trap (a GSE may wrap multiple SubSeries on different platforms; default download mixes them), the series-matrix normalization-trust caveat (submitter-normalized; re-process from raw for reproducibility), processed-vs-raw decision (Affymetrix CEL + locally-run RMA; RNA-seq via SRA FASTQ), GEOparse vs GEOquery trade-off, GEOmetadb deprecation (2020), and ArrayExpress migration to BioStudies (2020).

## Prerequisites

```bash
pip install biopython GEOparse pandas pysradb
# OR for R:
# R: BiocManager::install('GEOquery')
```

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional'
```

## Quick Start

- "Search GEO for human single-cell RNA-seq from 2024; flag any SuperSeries"
- "Detect whether GSE122288 is a SuperSeries before downloading; if so, list its SubSeries"
- "Resolve GSE123456 to SRA run accessions via pysradb; hand off to sra-data for FASTQ download"
- "Download the series matrix for GSE123456 and tell me what 'normalization' the submitter applied"
- "Pull supplementary CEL files for GSE12345 via R's GEOquery -- GEOparse has been flakey on suppl downloads"

## Example Prompts

### Detecting SuperSeries before pulling

> "Before I download GSE122288 as a single experiment, check its SOFT family file for !Series_relation. If it's a SuperSeries, list each SubSeries and recommend processing them independently to avoid mixing platforms."

### Processed-vs-raw decision

> "I want expression values for GSE123456. Check whether the technology is Affymetrix (-> download CEL files and re-do RMA) or RNA-seq (-> resolve to SRA and re-quantify with Salmon). Don't trust the submitter's series-matrix values for downstream stats."

### GEO -> SRA linkage

> "Find SRA accessions for GSE123456 using pysradb (more reliable than the gds -> sra ELink). Return a list of SRR run IDs to hand off to the sra-data skill."

### Cross-reference from publication

> "Find all GEO datasets cited in PMID 35412348 via pubmed -> gds ELink. Summarize each with title, sample count, platform, and SuperSeries status."

### Submitter-data-processing audit

> "Download the series matrix for GSE123456 and dump every unique value of !Sample_data_processing. Tell me whether the matrix is raw counts, log-CPM, VST, or RMA-normalized -- the answer determines whether I can use it as-is."

## What the Agent Will Do

1. Search gds db with field-qualified terms (`gse[Entry Type]`, `Homo sapiens[Organism]`, `expression profiling by high throughput sequencing[GDS Type]`).
2. For any GSE returned, check `!Series_relation` in SOFT to detect SuperSeries before pulling.
3. Pick the right download path: series matrix for fast-and-trusting; supplementary files for raw Affymetrix / submitter counts; SRA-link for RNA-seq raw FASTQ.
4. Read `!Sample_data_processing` to surface what's actually in the series matrix.
5. For R-side analyses, recommend GEOquery (Bioconductor) over GEOparse for supplementary file reliability.
6. For SRA hand-off, use pysradb to resolve GSE -> SRP -> SRR; pass run list to sra-data skill.
7. Warn on stale GEOmetadb usage; recommend pysradb / Entrez gds.
8. For ArrayExpress accessions (E-MTAB-*), use the new BioStudies URL.

## Tips

- The SuperSeries trap is the single biggest GEO mistake. ALWAYS check `!Series_relation` for `SuperSeries of: ...` before treating a GSE as one experiment.
- Submitter-normalized series matrices vary widely. For Affymetrix, default to CEL + locally-run RMA. For RNA-seq, default to SRA FASTQ + locally-run quantification.
- GEOparse (Python) is OK for SOFT parsing but flakey on supplementary-file download since ~2022. For Python pipelines that need suppl files, use direct FTP (`wget -r` on `suppl/`).
- For R workflows, GEOquery is the mature choice; getGEOSuppFiles is more reliable than GEOparse's equivalent.
- `[Entry Type]` is case-sensitive: `gse[Entry Type]` works; `gse[entry_type]` returns empty.
- GEOmetadb is unmaintained since 2020. Use pysradb for GSE<->SRA mapping; Entrez gds for full GEO search.
- ArrayExpress was migrated into BioStudies in 2020; old `E-MTAB-*` accessions still resolve via `https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-XXXX`.
- For curated re-processed RNA-seq across thousands of GEO studies, ARCHS4 (https://archs4.org) and recount3 are alternatives with consistent pipelines.

## Related Skills

- entrez-search - General gds search
- entrez-link - gds <-> sra / pubmed / bioproject links
- sra-data - Pull raw FASTQ from GEO-linked SRA
- expression-matrix/normalization - Re-normalize raw expression data
- rna-quantification/alignment-free-quant - Salmon/kallisto re-quantification
- ensembl-rest - Cross-reference Ensembl IDs in series matrices
