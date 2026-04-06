---
name: bio-expression-matrix-gene-id-mapping
description: Convert between gene identifier systems including Ensembl, Entrez, HGNC symbols, and UniProt. Use when mapping IDs for pathway analysis or matching different data sources.
tool_type: mixed
primary_tool: biomaRt
---

## Version Compatibility

Reference examples tested with: pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Gene ID Mapping

## Python: mygene

**Goal:** Convert between gene identifier systems (Ensembl, Entrez, Symbol, UniProt) using the MyGene.info API.

**Approach:** Query mygene with source IDs, specifying scopes and target fields, to build an ID mapping dictionary.

**"Convert my Ensembl gene IDs to gene symbols"** → Query a gene annotation service to map between identifier systems, handling one-to-many mappings.

```python
import mygene
import pandas as pd

mg = mygene.MyGeneInfo()

# Ensembl to Symbol
ensembl_ids = ['ENSG00000141510', 'ENSG00000012048', 'ENSG00000141736']
results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
mapping = {r['query']: r.get('symbol', None) for r in results}
# {'ENSG00000141510': 'TP53', 'ENSG00000012048': 'BRCA1', 'ENSG00000141736': 'ERBB2'}

# Symbol to Entrez
symbols = ['TP53', 'BRCA1', 'ERBB2']
results = mg.querymany(symbols, scopes='symbol', fields='entrezgene', species='human')
mapping = {r['query']: r.get('entrezgene', None) for r in results}

# Ensembl to multiple fields
results = mg.querymany(ensembl_ids, scopes='ensembl.gene',
    fields=['symbol', 'entrezgene', 'uniprot'], species='human')
```

## Python: pyensembl

**Goal:** Map gene identifiers using a local Ensembl database for offline, fast lookups.

**Approach:** Load a specific Ensembl release and query gene objects by ID or name.

```python
from pyensembl import EnsemblRelease

# Load Ensembl release (downloads automatically first time)
ensembl = EnsemblRelease(110, species='human')  # or 'mouse'

# Gene ID to symbol
gene = ensembl.gene_by_id('ENSG00000141510')
print(gene.gene_name)  # TP53

# Symbol to gene ID
gene = ensembl.genes_by_name('TP53')[0]
print(gene.gene_id)  # ENSG00000141510

# Batch conversion
def ensembl_to_symbol(ensembl_ids, release=110):
    ens = EnsemblRelease(release, species='human')
    mapping = {}
    for eid in ensembl_ids:
        try:
            gene = ens.gene_by_id(eid.split('.')[0])  # Remove version
            mapping[eid] = gene.gene_name
        except ValueError:
            mapping[eid] = None
    return mapping
```

## Python: gseapy

```python
import gseapy as gp

# Ensembl to Symbol using Enrichr
gene_list = ['ENSG00000141510', 'ENSG00000012048']
converted = gp.biomart.ensembl2name(gene_list, organism='hsapiens')
```

## R: biomaRt

**Goal:** Map gene identifiers using the Ensembl BioMart web service in R.

**Approach:** Connect to the Ensembl BioMart and retrieve attribute mappings for a list of gene IDs.

```r
library(biomaRt)

# Connect to Ensembl
ensembl <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')

# Ensembl to Symbol
ensembl_ids <- c('ENSG00000141510', 'ENSG00000012048', 'ENSG00000141736')
results <- getBM(
    attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'),
    filters='ensembl_gene_id',
    values=ensembl_ids,
    mart=ensembl
)

# Symbol to Ensembl
symbols <- c('TP53', 'BRCA1', 'ERBB2')
results <- getBM(
    attributes=c('hgnc_symbol', 'ensembl_gene_id'),
    filters='hgnc_symbol',
    values=symbols,
    mart=ensembl
)

# All available attributes
listAttributes(ensembl)
```

## R: org.db Packages

**Goal:** Map gene identifiers using Bioconductor organism annotation packages for fast local lookups.

**Approach:** Use mapIds from AnnotationDbi with organism-specific org.db packages.

```r
library(org.Hs.eg.db)  # Human
library(AnnotationDbi)

# Ensembl to Symbol
ensembl_ids <- c('ENSG00000141510', 'ENSG00000012048')
symbols <- mapIds(org.Hs.eg.db, keys=ensembl_ids, keytype='ENSEMBL', column='SYMBOL')

# Symbol to Entrez
symbols <- c('TP53', 'BRCA1')
entrez <- mapIds(org.Hs.eg.db, keys=symbols, keytype='SYMBOL', column='ENTREZID')

# Available keytypes
keytypes(org.Hs.eg.db)
# ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, SYMBOL, UNIPROT, etc.
```

## Apply Mapping to Count Matrix

**Goal:** Replace gene IDs in a count matrix index with a different identifier type.

**Approach:** Map IDs via mygene, update the DataFrame index, and aggregate duplicates by summing.

**"Convert the gene IDs in my count matrix from Ensembl to symbols"** → Map the row index to a new ID type, handling version suffixes and duplicate mappings by summation.

```python
import pandas as pd
import mygene

def map_count_matrix_ids(counts, from_type='ensembl.gene', to_type='symbol', species='human'):
    '''Map gene IDs in count matrix index.'''
    mg = mygene.MyGeneInfo()

    # Remove version numbers from Ensembl IDs
    clean_ids = [g.split('.')[0] for g in counts.index]

    # Query mygene
    results = mg.querymany(clean_ids, scopes=from_type, fields=to_type, species=species)

    # Build mapping
    mapping = {}
    for r in results:
        if to_type in r:
            mapping[r['query']] = r[to_type]

    # Apply mapping
    new_index = [mapping.get(g.split('.')[0], g) for g in counts.index]
    counts_mapped = counts.copy()
    counts_mapped.index = new_index

    # Handle duplicates (sum)
    counts_mapped = counts_mapped.groupby(counts_mapped.index).sum()

    return counts_mapped

# Usage
counts_symbols = map_count_matrix_ids(counts, 'ensembl.gene', 'symbol')
```

## R Equivalent

**Goal:** Replace gene IDs in an R count matrix using biomaRt with duplicate aggregation.

**Approach:** Query BioMart for the mapping, merge with the count matrix, and sum duplicate rows.

```r
library(biomaRt)

map_count_matrix_ids <- function(counts, from_type='ensembl_gene_id', to_type='hgnc_symbol') {
    ensembl <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')

    # Remove version numbers
    clean_ids <- gsub('\\..*', '', rownames(counts))

    # Get mapping
    mapping <- getBM(
        attributes=c(from_type, to_type),
        filters=from_type,
        values=clean_ids,
        mart=ensembl
    )

    # Merge and aggregate duplicates
    counts$gene_id <- clean_ids
    merged <- merge(counts, mapping, by.x='gene_id', by.y=from_type, all.x=TRUE)
    merged$gene_id <- NULL

    # Use symbol as rowname, sum duplicates
    rownames(merged) <- merged[[to_type]]
    merged[[to_type]] <- NULL
    counts_mapped <- aggregate(. ~ rownames(merged), data=merged, FUN=sum)
    rownames(counts_mapped) <- counts_mapped[,1]
    counts_mapped <- counts_mapped[,-1]

    return(counts_mapped)
}
```

## Handle Unmapped IDs

**Goal:** Track and gracefully handle gene IDs that fail to map to the target identifier system.

**Approach:** Keep original IDs for unmapped genes and report mapping success rate.

```python
def robust_id_mapping(gene_ids, from_type, to_type, species='human'):
    '''Map IDs with fallback for unmapped genes.'''
    import mygene
    mg = mygene.MyGeneInfo()

    clean_ids = [g.split('.')[0] for g in gene_ids]
    results = mg.querymany(clean_ids, scopes=from_type, fields=to_type, species=species)

    mapping = {}
    unmapped = []
    for r in results:
        original = gene_ids[clean_ids.index(r['query'])]
        if to_type in r:
            mapping[original] = r[to_type]
        else:
            mapping[original] = original  # Keep original if unmapped
            unmapped.append(original)

    print(f'Mapped: {len(gene_ids) - len(unmapped)}/{len(gene_ids)}')
    print(f'Unmapped: {len(unmapped)}')

    return mapping, unmapped
```

## Common ID Types and Database Selection

| Type | Example | Use Case | Stability |
|------|---------|----------|-----------|
| Ensembl Gene | ENSG00000141510 | RNA-seq, GTF files | Stable across releases (versioned) |
| Ensembl Transcript | ENST00000269305 | Transcript-level analysis | Stable (versioned) |
| Entrez Gene | 7157 | NCBI databases, KEGG | Stable (never reused) |
| HGNC Symbol | TP53 | Human readable display | Changes frequently |
| UniProt | P04637 | Protein databases | Stable (versioned releases) |
| RefSeq | NM_000546 | NCBI RefSeq | Stable (versioned) |

### Database Selection Guide

| Scenario | Recommended ID | Why |
|----------|---------------|-----|
| Computational key / primary index | Ensembl Gene ID | Stable, versioned, consistent with GTF |
| Pathway analysis (KEGG, Reactome) | Entrez Gene ID | Required by most pathway databases |
| GO enrichment | Entrez or Ensembl | Both supported by clusterProfiler |
| Display labels (plots, tables) | HGNC Symbol | Human-readable |
| Cross-database integration | Ensembl | Best-connected hub across databases |
| Protein-level analysis | UniProt | Primary protein database |

Best practice: use stable IDs (Ensembl or Entrez) as computational keys. Use symbols only as display labels. Always pin mappings to a specific database release and archive the cross-reference table for reproducibility.

### Gene Symbol Instability

Gene symbols change regularly as nomenclature committees update names. NCBI updates daily; Bioconductor org.db packages update every 6 months. For the most current mappings, query mygene.info or download gene_info from NCBI FTP directly. Never use symbols as the primary key in a pipeline -- always join on stable IDs and add symbols as a display column.

## PAR Gene Complications

Pseudo-autosomal region (PAR) genes exist on both X and Y chromosomes with identical sequences. In Ensembl GTF files, PAR genes have coordinates on both chromosomes, potentially creating duplicate entries in count matrices. Reads from PAR regions cannot be unambiguously assigned to X or Y.

```python
# Check for PAR gene duplicates in a count matrix
par_genes_human = ['SHOX', 'IL3RA', 'SLC25A6', 'P2RY8', 'AKAP17A', 'ASMT', 'DHRSX']
duplicated_ids = counts.index[counts.index.duplicated()].unique()
if len(duplicated_ids) > 0:
    print(f'Duplicate gene entries found: {len(duplicated_ids)}')
    # Sum duplicates (standard approach for PAR genes)
    counts = counts.groupby(counts.index).sum()
```

Some reference genomes mask the Y-chromosome PAR to avoid double-counting. Check whether the GTF includes PAR genes on both chromosomes before building count matrices.

## Cross-Species Ortholog Mapping

**Goal:** Map gene IDs between species for cross-species comparisons or integration.

**Approach:** Use Ensembl Compara (via biomaRt) to find orthologs, selecting the appropriate stringency level.

```r
library(biomaRt)

human <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
mouse <- useEnsembl(biomart='genes', dataset='mmusculus_gene_ensembl')

# Human to mouse one-to-one orthologs
orthologs <- getLDS(
    attributes=c('hgnc_symbol', 'ensembl_gene_id'),
    filters='ensembl_gene_id',
    values=human_gene_ids,
    mart=human,
    attributesL=c('mgi_symbol', 'ensembl_gene_id'),
    martL=mouse
)
```

| Strategy | When to use | Trade-off |
|----------|------------|-----------|
| One-to-one orthologs only | Cross-species scRNA-seq integration | Most conservative; loses genes without clear orthologs |
| Include one-to-many | Broader gene coverage needed | Must select: highest homology confidence or highest expression |
| Include many-to-many | Maximum inclusivity | Introduces ambiguity; use with caution |

For cross-species scRNA-seq integration, use only one-to-one orthologs (standard practice).

## Build tx2gene for tximport

**Goal:** Create the transcript-to-gene mapping table required by tximport for gene-level summarization from Salmon/kallisto output.

**Approach:** Extract transcript-gene relationships from a GTF file or Ensembl BioMart.

```r
# From GTF (recommended for consistency with quantification index)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('annotation.gtf')
k <- keys(txdb, keytype='TXNAME')
tx2gene <- AnnotationDbi::select(txdb, k, 'GENEID', 'TXNAME')

# From BioMart
library(biomaRt)
mart <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
tx2gene <- getBM(
    attributes=c('ensembl_transcript_id_version', 'ensembl_gene_id_version'),
    mart=mart
)
colnames(tx2gene) <- c('TXNAME', 'GENEID')
```

```python
# Python: extract tx2gene from GTF
import pandas as pd

def tx2gene_from_gtf(gtf_path):
    '''Extract transcript-to-gene mapping from GTF.'''
    records = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#') or '\ttranscript\t' not in line:
                continue
            attrs = line.strip().split('\t')[8]
            gene_id = [a.split('"')[1] for a in attrs.split(';') if 'gene_id' in a][0]
            tx_id = [a.split('"')[1] for a in attrs.split(';') if 'transcript_id' in a][0]
            records.append({'TXNAME': tx_id, 'GENEID': gene_id})
    return pd.DataFrame(records).drop_duplicates()
```

## Related Skills

- expression-matrix/counts-ingest - Load count data
- expression-matrix/metadata-joins - Add annotations
- rna-quantification/tximport-workflow - Uses tx2gene mapping
- pathway-analysis/go-enrichment - Requires Entrez IDs
- pathway-analysis/kegg-pathways - Requires Entrez IDs
