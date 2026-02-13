---
name: bio-small-rna-seq-target-prediction
description: Predict miRNA target genes using sequence-based algorithms and database lookups. Use when identifying potential mRNA targets of differentially expressed or functionally important miRNAs.
tool_type: mixed
primary_tool: miRanda
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# miRNA Target Prediction

**"Predict target genes for my miRNAs"** → Identify potential mRNA targets of miRNAs using sequence-based thermodynamic alignment and database lookups from validated interaction repositories.
- CLI: `miranda miRNA.fa UTR.fa -sc 140 -en -20` for de novo prediction
- Python: API queries to miRTarBase, TargetScan for validated targets

## miRanda Algorithm

**Goal:** Predict miRNA-mRNA target interactions using thermodynamic alignment scoring.

**Approach:** Run miRanda to align miRNA sequences against 3' UTR sequences with minimum score and energy thresholds.

```bash
# Run miRanda for target prediction
miranda miRNA.fa UTRs.fa \
    -sc 140 \
    -en -20 \
    -out predictions.txt

# Options:
# -sc 140: Minimum alignment score (default 140)
# -en -20: Maximum free energy threshold (kcal/mol)
# Higher score and lower energy = stronger prediction
```

## Parse miRanda Output

**Goal:** Extract miRNA-target interaction records from miRanda output into a structured DataFrame.

**Approach:** Parse the tab-delimited output lines starting with '>' to extract miRNA, target, score, energy, and position fields.

```python
import pandas as pd

def parse_miranda(output_file):
    '''Parse miRanda output file'''
    results = []
    with open(output_file) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    results.append({
                        'mirna': parts[0].lstrip('>'),
                        'target': parts[1],
                        'score': float(parts[2]),
                        'energy': float(parts[3]),
                        'position': parts[4]
                    })
    return pd.DataFrame(results)
```

## TargetScan Database Lookup

**Goal:** Retrieve conserved miRNA target predictions from the TargetScan database.

**Approach:** Query the downloadable TargetScan context++ score file by miRNA family name and rank by prediction score.

```python
import requests
import pandas as pd

def query_targetscan(mirna_family):
    '''Query TargetScan for predicted targets

    Note: TargetScan uses miRNA family names (e.g., miR-21-5p)
    '''
    # TargetScan provides downloadable files
    # For human: https://www.targetscan.org/vert_80/vert_80_data_download/
    targetscan_file = 'Predicted_Targets_Context_Scores.txt'

    df = pd.read_csv(targetscan_file, sep='\t')
    targets = df[df['miRNA family'] == mirna_family]
    return targets.sort_values('context++ score')
```

## miRDB Database Lookup

**Goal:** Retrieve machine-learning-based miRNA target predictions from miRDB.

**Approach:** Query the miRDB prediction file by miRNA ID and filter for high-confidence targets (score >= 80).

```python
def query_mirdb(mirna_id):
    '''Query miRDB for target predictions

    miRDB uses machine learning for target prediction
    Score > 80 indicates high confidence
    '''
    # Download from http://mirdb.org/download.html
    mirdb_file = 'miRDB_v6.0_prediction_result.txt'

    df = pd.read_csv(mirdb_file, sep='\t', header=None,
                     names=['mirna', 'target', 'score'])
    targets = df[df['mirna'] == mirna_id]
    return targets[targets['score'] >= 80].sort_values('score', ascending=False)
```

## Combine Multiple Databases

**Goal:** Identify high-confidence miRNA targets predicted by multiple independent algorithms.

**Approach:** Compute the intersection of predictions across miRanda, TargetScan, and miRDB, keeping targets found in at least N databases.

```python
def consensus_targets(mirna, min_databases=2):
    '''Find targets predicted by multiple databases

    More reliable targets are predicted by multiple algorithms
    '''
    miranda_targets = set(query_miranda_targets(mirna))
    targetscan_targets = set(query_targetscan_targets(mirna))
    mirdb_targets = set(query_mirdb_targets(mirna))

    # Count predictions per target
    all_targets = miranda_targets | targetscan_targets | mirdb_targets
    consensus = []

    for target in all_targets:
        count = sum([
            target in miranda_targets,
            target in targetscan_targets,
            target in mirdb_targets
        ])
        if count >= min_databases:
            consensus.append({
                'target': target,
                'n_databases': count,
                'miranda': target in miranda_targets,
                'targetscan': target in targetscan_targets,
                'mirdb': target in mirdb_targets
            })

    return pd.DataFrame(consensus).sort_values('n_databases', ascending=False)
```

## Python miRNA Target Prediction

**Goal:** Retrieve experimentally validated miRNA-target interactions from miRTarBase.

**Approach:** Load the miRTarBase Excel download and filter by miRNA name to get validated targets with experimental evidence types.

```python
# Using mirtarbase package for validated targets
def get_validated_targets(mirna):
    '''Get experimentally validated targets from miRTarBase'''
    # Download from https://mirtarbase.cuhk.edu.cn/
    mirtarbase_file = 'miRTarBase_MTI.xlsx'

    df = pd.read_excel(mirtarbase_file)
    validated = df[df['miRNA'] == mirna]
    return validated[['Target Gene', 'Experiments', 'Support Type']]
```

## Seed Match Analysis

**Goal:** Find miRNA seed region complementary matches within a 3' UTR sequence.

**Approach:** Extract the 7-mer seed (positions 2-8), compute its reverse complement, and scan the UTR for all occurrences.

```python
from Bio.Seq import Seq

def find_seed_matches(mirna_seq, utr_seq):
    '''Find seed matches in UTR sequence

    Seed region: positions 2-8 of miRNA (7-mer)
    '''
    mirna = Seq(mirna_seq)
    utr = Seq(utr_seq)

    # Get seed (positions 2-8, 0-indexed: 1-7)
    seed = str(mirna[1:8])
    seed_rc = str(Seq(seed).reverse_complement())

    matches = []
    start = 0
    while True:
        pos = str(utr).find(seed_rc, start)
        if pos == -1:
            break
        matches.append(pos)
        start = pos + 1

    return matches
```

## Functional Enrichment of Targets

**Goal:** Identify biological functions enriched among predicted miRNA target genes.

**Approach:** Run GO and KEGG enrichment analysis on the target gene list using Enrichr via gseapy.

```python
def enrich_target_genes(targets, background=None):
    '''Run GO enrichment on predicted target genes'''
    import gseapy as gp

    enr = gp.enrichr(
        gene_list=targets,
        gene_sets=['GO_Biological_Process_2021', 'KEGG_2021_Human'],
        organism='Human'
    )
    return enr.results
```

## Related Skills

- differential-mirna - Get DE miRNAs for target prediction
- pathway-analysis/go-enrichment - Enrich target gene functions
- database-access/entrez-fetch - Query biological databases
