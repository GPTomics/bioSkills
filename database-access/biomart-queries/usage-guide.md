# BioMart Queries Usage Guide

## Overview

Bulk-query Ensembl BioMart for cross-database ID mapping, coordinate tables, and ortholog wide tables. Encodes: BioMart vs Ensembl REST decision (BioMart for >5K rows, REST for per-record), the mart-vs-dataset hierarchy, attribute/filter discovery, many-to-many row-multiplication when joining cross-references, R biomaRt vs Python pybiomart trade-off, and version pinning via `useEnsembl(version=110)`.

## Prerequisites

```bash
pip install pybiomart pandas
# R alternative (more mature):
# BiocManager::install('biomaRt')
```

For published work, pin the Ensembl release in code.

## Quick Start

- "Convert 5,000 Ensembl Gene IDs to HGNC + RefSeq + UniProt in one query"
- "Pull all protein-coding genes on chr17 with coordinates and biotype"
- "Get a wide ortholog table: human Ensembl Gene ID, mouse ortholog, zebrafish ortholog"
- "Fetch GO term annotations for a list of genes (long format)"
- "Pin BioMart to release 110 for reproducibility"

## Example Prompts

### Bulk ID mapping

> "I have 8,000 Ensembl Gene IDs. Convert them to HGNC symbol, NCBI Entrez Gene ID, RefSeq mRNA accessions, and Swiss-Prot UniProt accessions. Use one pybiomart query against hsapiens_gene_ensembl with filters={'ensembl_gene_id': [...]}. Don't loop Ensembl REST -- that's 8,000 sequential calls."

### Coordinate table

> "Pull all protein-coding genes on chromosome 17 with: ensembl_gene_id, external_gene_name, start, end, strand, biotype. One BioMart query."

### Ortholog wide table

> "Build a wide table: for every human protein-coding gene, the mouse ortholog Ensembl ID and the zebrafish ortholog Ensembl ID, both filtered to ortholog_one2one only. Use mmusculus_homolog_ensembl_gene + drerio_homolog_ensembl_gene attributes."

### GO annotations long format

> "Get GO term annotations for [TP53, BRCA1, MYC, EGFR]. Returns long format -- one row per (gene, GO term)."

### Version pinning

> "I'm writing a paper. Use R biomaRt with useEnsembl(version=110) so the query reproduces in 2030. Don't use the default which follows the current release."

## What the Agent Will Do

1. Switch from Ensembl REST to BioMart for any query expected to return >5,000 rows.
2. Discover marts and datasets via `server.marts` and `mart.datasets` before guessing names.
3. List available `attributes` and `filters` programmatically; don't hallucinate field names.
4. Filter by Ensembl Gene ID (stable) over gene symbol (HGNC renames break workflows).
5. For ortholog queries, prefer `ortholog_one2one` for high-confidence joins.
6. Pin Ensembl release with `useEnsembl(version=N)` for reproducibility.
7. Surface row-count surprises caused by many-to-many cross-ref joins (e.g. 5K genes -> 50K rows).
8. Chunk very large queries by chromosome to avoid timeouts.

## Tips

- BioMart query returns one flat TSV. No pagination, no per-record loop, no rate-limit cascade.
- Discover attributes and filters with `ds.attributes` / `ds.filters` -- don't guess.
- Cross-reference joins multiply rows. A gene with 10 RefSeq mRNAs produces 10 rows in a query that includes refseq_mrna. Filter to canonical isoforms downstream when one-per-gene is needed.
- HGNC symbols are unstable (MARCH1 -> MARCHF1 in 2020). Filter by `ensembl_gene_id` or `hgnc_id` for stability.
- For Ensembl release pinning: R `useEnsembl(version=110)`, or use archive host URL (`https://nov2020.archive.ensembl.org`).
- pybiomart (Python) is light; R biomaRt is more mature and Bioconductor-supported. For R-based pipelines, biomaRt is the canonical choice.
- For non-vertebrate species, switch to the Ensembl Genomes BioMart at `http://plants.ensembl.org/biomart/martservice`.
- For variant data, use `ENSEMBL_MART_SNP` not `ENSEMBL_MART_ENSEMBL`.
- For real-time per-record queries, Ensembl REST is the better tool. BioMart is batch-oriented.

## Related Skills

- ensembl-rest - Per-record Ensembl queries (complement to BioMart)
- ortholog-inference - Compara orthologs with confidence semantics
- uniprot-access - Alternative ID-mapping via UniProt
- ncbi-datasets-cli - NCBI-side bulk path
- entrez-search - NCBI alternative for non-Ensembl queries
