# Gene ID Mapping - Usage Guide

## Overview
Convert between different gene identifier systems (Ensembl, Symbol, Entrez, UniProt) using mygene, biomaRt, and org.db packages for cross-database integration.

## Prerequisites
```bash
pip install mygene pandas pyensembl
```

```r
install.packages("BiocManager")
BiocManager::install(c("biomaRt", "org.Hs.eg.db", "AnnotationDbi"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Convert my Ensembl gene IDs to gene symbols for visualization"
- "Map gene IDs to Entrez for KEGG pathway analysis"
- "Convert my count matrix index from Ensembl to symbols"

## Example Prompts
### Basic Conversion
> "Convert these Ensembl IDs to gene symbols: ENSG00000141510, ENSG00000133703"

> "Map my count matrix gene IDs from Ensembl to Entrez for pathway analysis"

### Handling Edge Cases
> "Convert my Ensembl IDs to symbols but keep the original ID if no mapping is found"

> "Handle one-to-many mappings when converting to UniProt IDs"

### Species-Specific
> "Map mouse Ensembl IDs (ENSMUSG) to gene symbols"

> "Convert my zebrafish gene IDs using the appropriate database"

## What the Agent Will Do
1. Identify the source ID type and target ID type
2. Select appropriate mapping tool (mygene for Python, biomaRt/org.db for R)
3. Clean IDs (remove version suffixes like .15 from ENSG00000141510.15)
4. Perform batch query with caching for efficiency
5. Handle unmapped IDs and one-to-many mappings appropriately

## Common Scenarios

| From | To | When |
|------|-----|------|
| Ensembl | Symbol | Display/visualization |
| Ensembl | Entrez | KEGG/GO enrichment |
| Symbol | Ensembl | Match to GTF |
| Entrez | UniProt | Protein analysis |

## Python Workflow

```python
import pandas as pd
import mygene

class GeneMapper:
    def __init__(self, species='human'):
        self.mg = mygene.MyGeneInfo()
        self.species = species
        self.cache = {}

    def map_ids(self, ids, from_type, to_type):
        cache_key = (tuple(ids), from_type, to_type)
        if cache_key in self.cache:
            return self.cache[cache_key]

        clean_ids = [str(g).split('.')[0] for g in ids]
        results = self.mg.querymany(clean_ids, scopes=from_type, fields=to_type, species=self.species, verbose=False)

        mapping = {}
        for r in results:
            if to_type in r:
                val = r[to_type]
                if isinstance(val, list):
                    val = val[0]
                mapping[r['query']] = val

        self.cache[cache_key] = mapping
        return mapping

    def convert_counts(self, counts, from_type, to_type):
        mapping = self.map_ids(counts.index, from_type, to_type)
        new_index = [mapping.get(str(g).split('.')[0], g) for g in counts.index]
        result = counts.copy()
        result.index = new_index
        result = result[~result.index.duplicated(keep='first')]
        return result

mapper = GeneMapper('human')
counts_symbol = mapper.convert_counts(counts, 'ensembl.gene', 'symbol')
counts_entrez = mapper.convert_counts(counts, 'ensembl.gene', 'entrezgene')
```

## R Workflow

```r
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)

convert_ids_biomart <- function(ids, from_attr, to_attr, dataset='hsapiens_gene_ensembl') {
    ensembl <- useEnsembl(biomart='genes', dataset=dataset)
    results <- getBM(attributes=c(from_attr, to_attr), filters=from_attr, values=ids, mart=ensembl)
    mapping <- setNames(results[[to_attr]], results[[from_attr]])
    return(mapping)
}

convert_ids_orgdb <- function(ids, from_keytype, to_column, orgdb=org.Hs.eg.db) {
    mapping <- mapIds(orgdb, keys=ids, keytype=from_keytype, column=to_column, multiVals='first')
    return(mapping)
}

convert_counts <- function(counts, from_keytype, to_column) {
    clean_ids <- gsub('\\..*', '', rownames(counts))
    mapping <- convert_ids_orgdb(clean_ids, from_keytype, to_column)
    new_names <- ifelse(is.na(mapping[clean_ids]), clean_ids, mapping[clean_ids])
    rownames(counts) <- new_names
    counts <- aggregate(. ~ rownames(counts), data=counts, FUN=sum)
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    return(counts)
}
```

## Handling Edge Cases

### Unmapped IDs
```python
def safe_map(counts, mapper, from_type, to_type):
    mapping = mapper.map_ids(counts.index, from_type, to_type)
    new_index = []
    for g in counts.index:
        clean_g = str(g).split('.')[0]
        new_index.append(mapping.get(clean_g, g))
    counts.index = new_index
    return counts
```

### One-to-Many Mappings
```python
results = mg.querymany(['ENSG00000141510'], scopes='ensembl.gene', fields='uniprot.Swiss-Prot', species='human')

for r in results:
    uniprots = r.get('uniprot', {}).get('Swiss-Prot', [])
    if isinstance(uniprots, str):
        uniprots = [uniprots]
    print(f"{r['query']} -> {uniprots}")
```

### Deprecated/Retired IDs
```python
from pyensembl import EnsemblRelease

for release in [110, 100, 90, 75]:
    try:
        ens = EnsemblRelease(release, species='human')
        gene = ens.gene_by_id(ensembl_id)
        print(f'Found in release {release}: {gene.gene_name}')
        break
    except:
        continue
```

## Species-Specific Databases

| Species | org.db Package | Ensembl Dataset |
|---------|----------------|-----------------|
| Human | org.Hs.eg.db | hsapiens_gene_ensembl |
| Mouse | org.Mm.eg.db | mmusculus_gene_ensembl |
| Rat | org.Rn.eg.db | rnorvegicus_gene_ensembl |
| Zebrafish | org.Dr.eg.db | drerio_gene_ensembl |
| Fly | org.Dm.eg.db | dmelanogaster_gene_ensembl |
| Worm | org.Ce.eg.db | celegans_gene_ensembl |

## Validation

```python
def validate_mapping(original_ids, mapping, expected_mapped_pct=0.8):
    mapped = sum(1 for k, v in mapping.items() if v is not None)
    pct = mapped / len(original_ids)
    print(f'Mapped: {mapped}/{len(original_ids)} ({pct:.1%})')
    if pct < expected_mapped_pct:
        print(f'Warning: mapping rate below {expected_mapped_pct:.0%}')
        print('Check: correct species? correct ID type?')
    return pct >= expected_mapped_pct
```

## Tips
- Always batch queries - query multiple IDs at once rather than one at a time
- Cache results for reuse across analyses
- Use local databases (org.db packages) for faster lookups than API calls
- Remove version numbers from Ensembl IDs before mapping (ENSG00000141510.15 -> ENSG00000141510)
- Validate mapping rates - low rates often indicate wrong species or ID type
