# Gene ID Mapping Usage Guide

## Overview

Different databases and tools use different gene identifiers. Mapping between these systems is essential for integrating data from multiple sources and for using pathway analysis tools that require specific ID types.

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
        '''Map gene IDs with caching.'''
        cache_key = (tuple(ids), from_type, to_type)
        if cache_key in self.cache:
            return self.cache[cache_key]

        clean_ids = [str(g).split('.')[0] for g in ids]
        results = self.mg.querymany(clean_ids, scopes=from_type,
            fields=to_type, species=self.species, verbose=False)

        mapping = {}
        for r in results:
            if to_type in r:
                val = r[to_type]
                if isinstance(val, list):
                    val = val[0]  # Take first if multiple
                mapping[r['query']] = val

        self.cache[cache_key] = mapping
        return mapping

    def convert_counts(self, counts, from_type, to_type):
        '''Convert count matrix index.'''
        mapping = self.map_ids(counts.index, from_type, to_type)
        new_index = [mapping.get(str(g).split('.')[0], g) for g in counts.index]
        result = counts.copy()
        result.index = new_index
        result = result[~result.index.duplicated(keep='first')]  # or groupby().sum()
        return result

# Usage
mapper = GeneMapper('human')

# Ensembl to Symbol
counts_symbol = mapper.convert_counts(counts, 'ensembl.gene', 'symbol')

# Ensembl to Entrez
counts_entrez = mapper.convert_counts(counts, 'ensembl.gene', 'entrezgene')
```

## R Workflow

```r
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Method 1: biomaRt (more complete but slower)
convert_ids_biomart <- function(ids, from_attr, to_attr, dataset='hsapiens_gene_ensembl') {
    ensembl <- useEnsembl(biomart='genes', dataset=dataset)
    results <- getBM(
        attributes=c(from_attr, to_attr),
        filters=from_attr,
        values=ids,
        mart=ensembl
    )
    mapping <- setNames(results[[to_attr]], results[[from_attr]])
    return(mapping)
}

# Method 2: org.db (faster, local)
convert_ids_orgdb <- function(ids, from_keytype, to_column, orgdb=org.Hs.eg.db) {
    mapping <- mapIds(orgdb, keys=ids, keytype=from_keytype, column=to_column,
        multiVals='first')
    return(mapping)
}

# Convert count matrix
convert_counts <- function(counts, from_keytype, to_column) {
    clean_ids <- gsub('\\..*', '', rownames(counts))
    mapping <- convert_ids_orgdb(clean_ids, from_keytype, to_column)
    new_names <- ifelse(is.na(mapping[clean_ids]), clean_ids, mapping[clean_ids])
    rownames(counts) <- new_names
    # Sum duplicates
    counts <- aggregate(. ~ rownames(counts), data=counts, FUN=sum)
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    return(counts)
}
```

## Handling Edge Cases

### Unmapped IDs
```python
# Keep original ID if no mapping found
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
# Some IDs map to multiple targets
results = mg.querymany(['ENSG00000141510'], scopes='ensembl.gene',
    fields='uniprot.Swiss-Prot', species='human')

# Handle multiple results
for r in results:
    uniprots = r.get('uniprot', {}).get('Swiss-Prot', [])
    if isinstance(uniprots, str):
        uniprots = [uniprots]
    print(f"{r['query']} -> {uniprots}")
```

### Deprecated/Retired IDs
```python
# Use archived Ensembl for old IDs
from pyensembl import EnsemblRelease

# Try current release first, fall back to older
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

## Performance Tips

1. **Batch queries**: Always query multiple IDs at once
2. **Cache results**: Store mappings for reuse
3. **Use local databases**: org.db packages faster than API calls
4. **Remove version numbers**: Clean IDs before mapping

```python
# Efficient batch query
ids = counts.index.tolist()
results = mg.querymany(ids, scopes='ensembl.gene', fields='symbol,entrezgene',
    species='human', as_dataframe=True)
```

## Validation

```python
def validate_mapping(original_ids, mapping, expected_mapped_pct=0.8):
    '''Check mapping quality.'''
    mapped = sum(1 for k, v in mapping.items() if v is not None)
    pct = mapped / len(original_ids)
    print(f'Mapped: {mapped}/{len(original_ids)} ({pct:.1%})')
    if pct < expected_mapped_pct:
        print(f'Warning: mapping rate below {expected_mapped_pct:.0%}')
        print('Check: correct species? correct ID type?')
    return pct >= expected_mapped_pct
```
