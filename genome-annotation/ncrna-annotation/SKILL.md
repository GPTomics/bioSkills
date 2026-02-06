---
name: bio-genome-annotation-ncrna-annotation
description: Identify non-coding RNAs including tRNAs, rRNAs, snoRNAs, and regulatory RNAs using Infernal covariance model searches against Rfam and tRNAscan-SE for tRNA prediction. Use when performing genome-wide ncRNA annotation with assembly input producing GFF output.
tool_type: cli
primary_tool: Infernal
---

# Non-Coding RNA Annotation

Identify and annotate non-coding RNAs in genome assemblies using Infernal (general ncRNAs via Rfam covariance models) and tRNAscan-SE (specialized tRNA detection).

## Infernal / cmscan

Infernal uses covariance models (CMs) from Rfam to identify ncRNA families by both sequence and secondary structure similarity.

### Rfam Database Setup

```bash
# Download Rfam covariance models (~600 MB compressed)
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz

# Press the CM database (required for cmscan)
cmpress Rfam.cm

# Download clan information (for overlap resolution)
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
```

### Basic cmscan

```bash
# Search genome against all Rfam families
# --cut_ga: Use Rfam gathering threshold (recommended, family-specific)
# --rfam: Run in Rfam mode (skip slow alignment for large models)
# --nohmmonly: Use CM scoring only (more accurate, slower)
# --fmt 2: Output in table format
cmscan \
    --cut_ga \
    --rfam \
    --nohmmonly \
    --tblout ncrna_hits.tbl \
    --fmt 2 \
    --cpu 16 \
    --clanin Rfam.clanin \
    Rfam.cm \
    genome.fasta > ncrna_hits.out
```

### Key Options

| Option | Description |
|--------|-------------|
| `--cut_ga` | Use Rfam gathering threshold (recommended) |
| `--rfam` | Speed optimization for Rfam-scale searches |
| `--nohmmonly` | Force CM-only scoring (more sensitive for structured RNAs) |
| `--clanin` | Rfam clan file for resolving overlapping hits |
| `--tblout` | Tabular output file |
| `--fmt 2` | Extended table format with accessions |
| `--cpu` | CPU threads |
| `-E` | E-value threshold (default: 10; --cut_ga is preferred) |

### Convert cmscan Output to GFF3

```bash
# Using Infernal's built-in conversion
grep -v '^#' ncrna_hits.tbl | \
awk 'BEGIN{OFS="\t"} {
    if ($10 == "+") strand = "+"; else strand = "-";
    if ($10 == "+") {start=$8; end=$9} else {start=$9; end=$8};
    print $4, "Infernal", "ncRNA", start, end, $17, strand, ".", "ID="$2";Name="$3";rfam_acc="$2";evalue="$17
}' > ncrna_infernal.gff3
```

### E-value Guidelines

| Category | Typical E-value | Notes |
|----------|----------------|-------|
| High confidence | < 1e-10 | Unambiguous family assignment |
| Moderate | 1e-10 to 1e-5 | Likely real, verify context |
| Rfam GA threshold | family-specific | Recommended cutoff with --cut_ga |
| Marginal | > 1e-3 | May be pseudogenes or degraded copies |

## tRNAscan-SE

Specialized tRNA detector using covariance models. More sensitive and specific for tRNAs than Infernal alone.

### Basic Usage

```bash
# Bacterial/archaeal mode
tRNAscan-SE -B -o trna_results.txt --gff trna.gff3 genome.fasta

# Eukaryotic mode
tRNAscan-SE -E -o trna_results.txt --gff trna.gff3 genome.fasta

# General mode (auto-detect domain)
tRNAscan-SE -G -o trna_results.txt --gff trna.gff3 genome.fasta
```

### Key Options

| Option | Description |
|--------|-------------|
| `-B` | Bacterial mode |
| `-A` | Archaeal mode |
| `-E` | Eukaryotic mode |
| `-G` | General (mixed/unknown) mode |
| `-o` | Tabular output |
| `--gff` | GFF3 output |
| `--detail` | Detailed output with isotype info |
| `--thread` | CPU threads |
| `-Q` | Covariance model scoring only (slower, more accurate) |

### Domain-Specific Models

| Domain | Flag | Notes |
|--------|------|-------|
| Bacteria | `-B` | Shorter introns, smaller genomes |
| Archaea | `-A` | Split tRNAs, bulge-helix-bulge introns |
| Eukaryota | `-E` | Longer introns, pseudogenes common |
| Organelle | `-O` | Mitochondrial/chloroplast tRNAs |
| Mitochondrial | `-M` | Mammalian mitochondrial-specific |

### Expected tRNA Counts

| Organism Type | Expected tRNAs | Notes |
|---------------|----------------|-------|
| Bacteria | 30-90 | Fewer isotypes, some shared |
| Archaea | 30-60 | Similar to bacteria |
| Yeast | 275-400 | Many isodecoders |
| Nematode | 600-900 | Gene family expansion |
| Mammal | 400-600 | Many pseudogenes |
| Plant | 500-1,000+ | Large gene families |

## barrnap (rRNA Detection)

barrnap predicts rRNA genes using HMM models. Note: barrnap has been unmaintained since 2018. Bakta handles rRNA detection internally for prokaryotes.

```bash
# Detect rRNAs
barrnap --kingdom bac genome.fasta > rrna.gff3

# Kingdoms: bac (bacteria), arc (archaea), euk (eukaryote), mito (mitochondria)
barrnap --kingdom euk --threads 4 genome.fasta > rrna.gff3
```

## Combining ncRNA Annotations

```python
import pandas as pd
from collections import defaultdict

def parse_infernal_tbl(tbl_file):
    '''Parse Infernal cmscan tabular output.'''
    hits = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 18:
                continue
            strand = '+' if parts[9] == '+' else '-'
            start, end = (int(parts[7]), int(parts[8])) if strand == '+' else (int(parts[8]), int(parts[7]))
            hits.append({
                'target': parts[0],
                'rfam_acc': parts[1],
                'rfam_name': parts[2],
                'seqid': parts[3],
                'start': start,
                'end': end,
                'strand': strand,
                'evalue': float(parts[16]),
                'score': float(parts[15]),
                'rna_type': classify_rfam(parts[2]),
            })
    return pd.DataFrame(hits)

def classify_rfam(rfam_name):
    '''Classify Rfam family into broad ncRNA categories.'''
    name_lower = rfam_name.lower()
    if 'rrna' in name_lower or 'ssu' in name_lower or 'lsu' in name_lower:
        return 'rRNA'
    if 'trna' in name_lower:
        return 'tRNA'
    if 'snorna' in name_lower or 'snord' in name_lower or 'snora' in name_lower:
        return 'snoRNA'
    if 'mir' in name_lower:
        return 'miRNA'
    if 'riboswitch' in name_lower or 'thermoregulator' in name_lower:
        return 'riboswitch'
    if 'ires' in name_lower or 'leader' in name_lower:
        return 'cis-reg'
    return 'other_ncRNA'

def parse_trnascan_gff(gff_file):
    '''Parse tRNAscan-SE GFF3 output.'''
    trnas = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'tRNA':
                continue
            attrs = dict(item.split('=') for item in parts[8].split(';') if '=' in item)
            trnas.append({
                'seqid': parts[0],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'strand': parts[6],
                'isotype': attrs.get('isotype', 'unknown'),
                'anticodon': attrs.get('anticodon', 'unknown'),
                'score': float(parts[5]) if parts[5] != '.' else 0,
                'rna_type': 'tRNA',
            })
    return pd.DataFrame(trnas)

def combine_ncrna_annotations(infernal_tbl, trnascan_gff, output_gff):
    '''Combine Infernal and tRNAscan-SE results, preferring tRNAscan for tRNAs.'''
    infernal_df = parse_infernal_tbl(infernal_tbl)
    trna_df = parse_trnascan_gff(trnascan_gff)

    # Remove Infernal tRNA hits (tRNAscan-SE is more accurate for tRNAs)
    infernal_no_trna = infernal_df[infernal_df['rna_type'] != 'tRNA']

    summary = defaultdict(int)
    for rna_type in infernal_no_trna['rna_type'].unique():
        summary[rna_type] = len(infernal_no_trna[infernal_no_trna['rna_type'] == rna_type])
    summary['tRNA'] = len(trna_df)

    print('=== ncRNA Summary ===')
    for rna_type, count in sorted(summary.items()):
        print(f'  {rna_type}: {count}')
    print(f'  Total: {sum(summary.values())}')

    return infernal_no_trna, trna_df, summary
```

## Troubleshooting

### cmscan Is Slow
- Use `--rfam` flag for Rfam-mode optimizations
- Use `--FZ 2` to further speed up (less sensitive)
- Split genome into chunks and run in parallel

### Missing Expected ncRNAs
- Verify Rfam.cm is current version
- Check that cmpress was run successfully (look for .i1m, .i1p, .i1f, .i1i files)
- For mitochondrial/chloroplast ncRNAs, search organelle sequence separately

### tRNAscan-SE Finds Too Many Pseudogenes
- Normal for eukaryotic genomes (especially mammals)
- Filter by score: high-confidence tRNAs typically score > 50 bits
- Use `--detail` flag to identify pseudogenes in output

## Related Skills

- prokaryotic-annotation - Bakta includes ncRNA annotation
- eukaryotic-gene-prediction - Gene prediction excludes ncRNAs
- rna-structure/ncrna-search - Targeted ncRNA homology searches
