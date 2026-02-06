---
name: bio-rna-structure-ncrna-search
description: Searches for non-coding RNA homologs and classifies RNA families using Infernal covariance model searches against the Rfam database. Identifies structured RNAs by sequence and secondary structure conservation. Use when querying sequences against Rfam, building custom covariance models for novel RNA families, or classifying non-coding transcripts by family.
tool_type: cli
primary_tool: Infernal
---

# ncRNA Search

Search for non-coding RNA homologs and classify RNA families using covariance models (CMs). Infernal scores both sequence and secondary structure conservation, making it more sensitive than sequence-only methods for structured RNAs.

## Rfam Database Setup

```bash
# Download current Rfam covariance models (~500 MB compressed)
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz

# Press the CM database (required before searching)
cmpress Rfam.cm

# Download clan information (for resolving overlapping hits)
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
```

## cmscan: Query Sequences Against Rfam

Search one or more query sequences against the full Rfam database to classify ncRNAs by family.

```bash
# Basic cmscan against Rfam
cmscan --cpu 8 --tblout results.tbl --fmt 2 Rfam.cm query.fa > results.out

# With clan overlap resolution (removes redundant hits from same clan)
cmscan --cpu 8 --tblout results.tbl --fmt 2 --clanin Rfam.clanin --oclan Rfam.cm query.fa > results.out

# Strict E-value threshold for high-confidence hits
cmscan --cpu 8 --tblout results.tbl --fmt 2 --clanin Rfam.clanin --oclan \
    -E 1e-5 --incE 1e-5 Rfam.cm query.fa > results.out
```

### cmscan Key Options

| Option | Description |
|--------|-------------|
| `--tblout` | Tabular output (easier to parse) |
| `--fmt 2` | Format 2 tabular output (includes model coordinates) |
| `-E` | Report E-value threshold (default: 10.0) |
| `--incE` | Inclusion E-value for significant hits (default: 0.01) |
| `--clanin` | Clan file for resolving overlapping family hits |
| `--oclan` | Enable clan overlap resolution |
| `--cut_ga` | Use Rfam gathering threshold (curated per-family cutoff) |
| `--cpu` | Number of threads |
| `--noali` | Skip alignment output (faster for large searches) |

### Gathering Threshold vs E-value

```bash
# Gathering threshold: curated per-family cutoff, recommended for Rfam
# Provides consistent sensitivity per family as calibrated by Rfam curators
cmscan --cut_ga --tblout results.tbl --fmt 2 --clanin Rfam.clanin --oclan \
    Rfam.cm query.fa > results.out

# E-value based: unified threshold, useful for custom databases
cmscan -E 1e-3 --tblout results.tbl Rfam.cm query.fa > results.out
```

## cmsearch: Search Specific CM Against Sequence Database

Search a specific covariance model against a sequence database (inverse of cmscan).

```bash
# Extract a single family CM from Rfam
cmfetch Rfam.cm RF00005 > tRNA.cm
cmpress tRNA.cm

# Search for tRNAs in a genome
cmsearch --cpu 8 --tblout trna_hits.tbl tRNA.cm genome.fa > trna_hits.out

# Search with bit score threshold instead of E-value
cmsearch --cpu 8 -T 30.0 --tblout hits.tbl tRNA.cm genome.fa > hits.out
```

## Building Custom Covariance Models

For novel RNA families not in Rfam, build a custom CM from a structure-annotated alignment.

### Step 1: Prepare Stockholm Alignment

```
# STOCKHOLM 1.0
#=GF AC   MYFAM00001
#=GF DE   My novel RNA family
seq1   GGGCUAUUAGCUCAGUUGGUUAGAGC
seq2   GGGCUAUAAGCUCAGUUGGAUAGAGC
seq3   GGGCUAUUAGCUCAGUUGGUUAGAGC
#=GC SS_cons  ((((....((((......))))))))
//
```

### Step 2: Build and Calibrate

```bash
# Build CM from alignment
cmbuild my_family.cm alignment.sto

# Calibrate E-value statistics (compute-intensive but essential for accurate E-values)
cmcalibrate --cpu 8 my_family.cm

# Press for searching
cmpress my_family.cm
```

### Step 3: Search

```bash
# Search genome with custom CM
cmsearch --cpu 8 --tblout custom_hits.tbl my_family.cm target_sequences.fa > custom_hits.out

# Iterative search: use hits to refine alignment and rebuild CM
cmsearch -A new_hits.sto my_family.cm target_sequences.fa
# Then manually curate new_hits.sto and rebuild
```

### cmbuild Options

| Option | Description |
|--------|-------------|
| `--hand` | Use reference annotation for consensus (trust SS_cons exactly) |
| `--enone` | Turn off entropy weighting |
| `-n` | Name the CM |
| `--ere` | Target mean match state relative entropy |

## Parsing Infernal Output

### Tabular Output (--tblout --fmt 2)

```python
import pandas as pd

def parse_cmscan_tblout(tblout_file):
    '''Parse Infernal cmscan --fmt 2 tabular output.'''
    rows = []
    with open(tblout_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 18:
                continue
            rows.append({
                'target_name': fields[0],
                'target_accession': fields[1],
                'query_name': fields[2],
                'query_accession': fields[3],
                'mdl_type': fields[4],
                'mdl_from': int(fields[5]),
                'mdl_to': int(fields[6]),
                'seq_from': int(fields[7]),
                'seq_to': int(fields[8]),
                'strand': fields[9],
                'trunc': fields[10],
                'pass': fields[11],
                'gc': float(fields[12]),
                'bias': float(fields[13]),
                'score': float(fields[14]),
                'evalue': float(fields[15]),
                'inc': fields[16],
                'description': ' '.join(fields[17:])
            })
    df = pd.DataFrame(rows)
    return df


def filter_significant_hits(df, evalue_threshold=1e-5):
    '''Filter for significant hits and sort by score.'''
    significant = df[df['evalue'] <= evalue_threshold].copy()
    significant = significant.sort_values('score', ascending=False)
    return significant


def summarize_families(df):
    '''Summarize ncRNA family assignments.'''
    summary = df.groupby('target_name').agg(
        count=('query_name', 'count'),
        best_score=('score', 'max'),
        best_evalue=('evalue', 'min')
    ).sort_values('count', ascending=False)
    return summary
```

### Extract Hit Sequences

```python
from Bio import SeqIO

def extract_hit_sequences(fasta_file, hits_df, output_file):
    '''Extract sequences for cmscan/cmsearch hits.'''
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    records = []
    for _, hit in hits_df.iterrows():
        seq_record = seqs[hit['query_name']]
        start, end = sorted([hit['seq_from'], hit['seq_to']])
        subseq = seq_record[start-1:end]
        if hit['strand'] == '-':
            subseq = subseq.reverse_complement()
        subseq.id = f'{hit["query_name"]}_{start}_{end}_{hit["target_name"]}'
        subseq.description = f'family={hit["target_name"]} score={hit["score"]:.1f} E={hit["evalue"]:.1e}'
        records.append(subseq)
    SeqIO.write(records, output_file, 'fasta')
    print(f'Extracted {len(records)} hit sequences to {output_file}')
```

## Quality Thresholds

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| E-value (Rfam scan) | < 1e-5 | High-confidence family assignment |
| Gathering threshold | --cut_ga | Rfam-curated per-family cutoffs, recommended default |
| Bit score | > 20 | Minimum for reportable hits in custom searches |
| Truncation | != 5'/3' | Hits at sequence edges may be truncated; check completeness |

## Related Skills

- secondary-structure-prediction - Predict structures for novel ncRNA candidates
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation pipelines
- alignment/msa-statistics - Evaluate alignment quality for CM building
