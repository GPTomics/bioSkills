---
name: bio-genome-annotation-functional-annotation
description: Assign GO terms, KEGG orthologs, Pfam domains, and EC numbers to predicted proteins using eggNOG-mapper and InterProScan. Produces functional summaries for downstream pathway and enrichment analysis. Use when adding functional annotation to predicted genes or characterizing protein functions in a new genome.
tool_type: cli
primary_tool: eggNOG-mapper
---

# Functional Annotation

Assign functional annotations (GO terms, KEGG orthologs, Pfam domains, EC numbers) to predicted protein sequences using eggNOG-mapper and InterProScan.

## eggNOG-mapper

### Database Setup

```bash
# Download eggNOG v5.0 database (~44 GB)
# Required for local searches; use --data_dir to specify location
download_eggnog_data.py --data_dir /path/to/eggnog_db -y

# Download DIAMOND database only (~9 GB, faster setup)
download_eggnog_data.py --data_dir /path/to/eggnog_db -y -D

# Download taxon-specific databases (optional, smaller)
download_eggnog_data.py --data_dir /path/to/eggnog_db -y -t 2 # Bacteria
download_eggnog_data.py --data_dir /path/to/eggnog_db -y -t 2759 # Eukaryota
```

### Basic Usage

```bash
emapper.py \
    -i predicted_proteins.faa \
    --output functional_annot \
    --output_dir eggnog_out \
    --data_dir /path/to/eggnog_db \
    --cpu 16 \
    -m diamond
```

### Key Options

| Option | Description |
|--------|-------------|
| `-i` | Input protein FASTA |
| `--output` | Output file prefix |
| `--data_dir` | Path to eggNOG database |
| `-m` | Search mode: diamond (fast), mmseqs (sensitive), hmmer |
| `--cpu` | CPU threads |
| `--tax_scope` | Taxonomic scope (auto, Bacteria, Eukaryota, etc.) |
| `--go_evidence` | GO evidence filter (experimental, non-electronic, all) |
| `--target_orthologs` | Ortholog type (one2one, all) |
| `--seed_ortholog_evalue` | E-value cutoff (default: 0.001) |
| `--seed_ortholog_score` | Min bit score (default: 60) |
| `--override` | Overwrite existing output |

### With Taxonomic Scope

```bash
# Restrict to bacterial orthologs for a prokaryotic genome
emapper.py \
    -i proteins.faa \
    --output annot \
    --output_dir eggnog_out \
    --data_dir /path/to/eggnog_db \
    --cpu 16 \
    -m diamond \
    --tax_scope Bacteria \
    --go_evidence non-electronic
```

### Output Files

```
eggnog_out/
├── annot.emapper.annotations    # Main annotation table
├── annot.emapper.hits           # DIAMOND/mmseqs hits
├── annot.emapper.seed_orthologs # Best orthologs
└── annot.emapper.pfam           # Pfam domain annotations
```

### Key Output Columns

| Column | Content |
|--------|---------|
| seed_ortholog | Best matching ortholog |
| evalue | E-value of best hit |
| GOs | GO term annotations |
| EC | Enzyme Commission numbers |
| KEGG_ko | KEGG ortholog IDs |
| KEGG_Pathway | KEGG pathway mappings |
| COG_category | COG functional category |
| PFAMs | Pfam domain annotations |
| Description | Functional description |

## InterProScan

InterProScan searches multiple protein signature databases simultaneously.

### Basic Usage

```bash
interproscan.sh \
    -i predicted_proteins.faa \
    -o interpro_results.tsv \
    -f tsv,gff3 \
    -cpu 16 \
    -goterms \
    -pa
```

### Key Options

| Option | Description |
|--------|-------------|
| `-i` | Input protein FASTA |
| `-o` | Output file |
| `-f` | Output formats: tsv, gff3, xml, json |
| `-cpu` | CPU threads |
| `-goterms` | Include GO term mappings |
| `-pa` | Include pathway annotations |
| `-appl` | Specific applications to run (comma-separated) |
| `-dp` | Disable precalculated match lookup |

### Select Specific Databases

```bash
# Run only Pfam, TIGRFAM, and CDD
interproscan.sh \
    -i proteins.faa \
    -o interpro_results.tsv \
    -f tsv,gff3 \
    -cpu 16 \
    -goterms -pa \
    -appl Pfam,TIGRFAM,CDD
```

### Available Applications

| Application | Description |
|-------------|-------------|
| Pfam | Protein families |
| TIGRFAM | Functionally equivalent protein families |
| SUPERFAMILY | Structural domain assignments |
| CDD | Conserved Domain Database |
| PANTHER | Protein classification |
| Gene3D | Structural domain predictions |
| Coils | Coiled-coil predictions |
| MobiDBLite | Disordered regions |
| SignalP | Signal peptides |
| TMHMM | Transmembrane helices |

## Merging eggNOG and InterProScan Results

```python
import pandas as pd

def parse_eggnog(annotations_file):
    '''Parse eggNOG-mapper annotations output.'''
    df = pd.read_csv(annotations_file, sep='\t', comment='#',
                     header=None, skiprows=5)
    col_names = [
        'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs',
        'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name',
        'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module',
        'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy',
        'BiGG_Reaction', 'PFAMs'
    ]
    df.columns = col_names[:len(df.columns)]
    return df

def parse_interproscan_tsv(tsv_file):
    '''Parse InterProScan TSV output.'''
    col_names = [
        'protein_id', 'md5', 'length', 'analysis', 'signature_acc',
        'signature_desc', 'start', 'stop', 'score', 'status', 'date',
        'interpro_acc', 'interpro_desc', 'go_terms', 'pathways'
    ]
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=col_names)
    return df

def merge_annotations(eggnog_file, interpro_file):
    '''Merge eggNOG and InterProScan annotations per protein.'''
    eggnog_df = parse_eggnog(eggnog_file)
    interpro_df = parse_interproscan_tsv(interpro_file)

    interpro_summary = interpro_df.groupby('protein_id').agg({
        'signature_acc': lambda x: ','.join(x.dropna().unique()),
        'interpro_acc': lambda x: ','.join(x.dropna().unique()),
        'go_terms': lambda x: '|'.join(x.dropna().unique()),
    }).reset_index()
    interpro_summary.columns = ['query', 'interpro_signatures', 'interpro_ids', 'interpro_go']

    merged = eggnog_df.merge(interpro_summary, on='query', how='outer')

    merged['all_go'] = merged.apply(
        lambda row: combine_go_terms(row.get('GOs', ''), row.get('interpro_go', '')), axis=1
    )
    return merged

def combine_go_terms(eggnog_go, interpro_go):
    '''Combine GO terms from both sources, removing duplicates.'''
    terms = set()
    for go_str in [eggnog_go, interpro_go]:
        if pd.notna(go_str) and go_str != '-':
            terms.update(t.strip() for t in str(go_str).replace('|', ',').split(',') if t.strip().startswith('GO:'))
    return ','.join(sorted(terms)) if terms else '-'
```

## Annotation Statistics

```python
def annotation_summary(merged_df):
    '''Summarize functional annotation coverage.'''
    total = len(merged_df)
    has_go = (merged_df['all_go'] != '-').sum()
    has_kegg = merged_df['KEGG_ko'].notna().sum() if 'KEGG_ko' in merged_df else 0
    has_pfam = merged_df['PFAMs'].notna().sum() if 'PFAMs' in merged_df else 0
    has_ec = merged_df['EC'].notna().sum() if 'EC' in merged_df else 0
    has_desc = (merged_df['Description'] != '-').sum() if 'Description' in merged_df else 0

    print(f'Total proteins: {total}')
    print(f'With GO terms: {has_go} ({has_go/total:.1%})')
    print(f'With KEGG orthologs: {has_kegg} ({has_kegg/total:.1%})')
    print(f'With Pfam domains: {has_pfam} ({has_pfam/total:.1%})')
    print(f'With EC numbers: {has_ec} ({has_ec/total:.1%})')
    print(f'With description: {has_desc} ({has_desc/total:.1%})')

    # Annotation coverage target: >60% with at least one functional term
    has_any = ((merged_df['all_go'] != '-') | merged_df['PFAMs'].notna() | merged_df['KEGG_ko'].notna()).sum()
    print(f'With any annotation: {has_any} ({has_any/total:.1%})')
```

## Troubleshooting

### Low Annotation Rate
- Check protein sequence quality (no fragmented ORFs)
- Try broader taxonomic scope (--tax_scope auto)
- Run both eggNOG-mapper and InterProScan and merge results

### eggNOG Database Errors
- Verify database version matches emapper version
- Re-download with `download_eggnog_data.py --data_dir /path -y`

### InterProScan Memory Issues
- Reduce batch size with `-b` option
- Split input FASTA into smaller chunks

## Related Skills

- prokaryotic-annotation - Bakta includes basic functional annotation
- eukaryotic-gene-prediction - Produces protein sequences for functional annotation
- pathway-analysis/go-enrichment - Enrichment analysis using GO annotations
- pathway-analysis/kegg-pathways - Pathway mapping with KEGG orthologs
