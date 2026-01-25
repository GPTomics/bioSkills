---
name: bio-crispr-library-design
description: CRISPR library design for genetic screens. Covers sgRNA selection, library composition, control design, and oligo ordering. Use when designing custom sgRNA libraries for knockout, activation, or interference screens.
tool_type: python
primary_tool: crispor
---

# Library Design

## sgRNA Selection Criteria

```python
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

def score_sgrna(sequence, pam='NGG'):
    '''Score sgRNA based on multiple criteria.'''
    scores = {}

    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    scores['gc_content'] = 1 - abs(gc_content - 0.5) * 2

    if len(sequence) >= 4:
        has_poly_t = 'TTTT' in sequence
        scores['poly_t'] = 0 if has_poly_t else 1

    starts_with_g = sequence.startswith('G')
    scores['start_g'] = 1 if starts_with_g else 0.5

    scores['length'] = 1 if len(sequence) == 20 else 0.8

    overall = np.mean(list(scores.values()))
    return overall, scores

def design_sgrnas_for_gene(gene_sequence, n_guides=4, pam='NGG'):
    '''Design sgRNAs targeting a gene.'''
    candidates = []

    pam_pattern = pam.replace('N', '[ACGT]')
    import re

    for strand in ['+', '-']:
        seq = gene_sequence if strand == '+' else str(Seq(gene_sequence).reverse_complement())

        for match in re.finditer(f'([ACGT]{{20}})({pam_pattern})', seq):
            sgrna = match.group(1)
            position = match.start()

            if strand == '-':
                position = len(seq) - position - 23

            score, details = score_sgrna(sgrna)

            candidates.append({
                'sequence': sgrna,
                'pam': match.group(2),
                'strand': strand,
                'position': position,
                'score': score,
                'gc_content': (sgrna.count('G') + sgrna.count('C')) / 20,
                **details
            })

    candidates_df = pd.DataFrame(candidates)
    candidates_df = candidates_df.sort_values('score', ascending=False)

    return candidates_df.head(n_guides)

gene_seq = 'ATGCGATCGATCGATCGATCGAATCGATCGATCGAGGCGATCGATCGATCGATCGAATCGATCGATCGAGGCGATCGATCGATCGATCGAATCGATCGATCGAGG'
guides = design_sgrnas_for_gene(gene_seq, n_guides=5)
print(guides[['sequence', 'position', 'strand', 'score', 'gc_content']])
```

## Library Composition

```python
def design_library(gene_list, guides_per_gene=4, include_controls=True):
    '''Design complete library for gene list.'''
    library = []

    for gene in gene_list:
        gene_data = get_gene_sequence(gene)

        guides = design_sgrnas_for_gene(gene_data['sequence'], n_guides=guides_per_gene)

        for idx, guide in guides.iterrows():
            library.append({
                'gene': gene,
                'gene_id': gene_data.get('ensembl_id', ''),
                'guide_number': idx + 1,
                'sequence': guide['sequence'],
                'pam': guide['pam'],
                'position': guide['position'],
                'strand': guide['strand'],
                'score': guide['score'],
                'type': 'targeting'
            })

    if include_controls:
        controls = design_control_guides()
        library.extend(controls)

    return pd.DataFrame(library)

def get_gene_sequence(gene_name):
    '''Fetch gene sequence (placeholder - use Ensembl API or local files).'''
    return {
        'sequence': 'ATGC' * 250,
        'ensembl_id': f'ENSG_{hash(gene_name) % 100000:05d}'
    }

genes = ['TP53', 'BRCA1', 'KRAS', 'MYC', 'CDK4']
library = design_library(genes, guides_per_gene=4)
print(f'Library size: {len(library)} guides')
print(f'Genes: {library["gene"].nunique()}')
```

## Control Guide Design

```python
def design_control_guides(n_nontargeting=100, n_essential=20, n_nonessential=20):
    '''Design control guides for library.'''
    controls = []

    for i in range(n_nontargeting):
        sequence = generate_nontargeting_sequence()
        controls.append({
            'gene': f'NonTargeting_{i+1}',
            'gene_id': '',
            'guide_number': 1,
            'sequence': sequence,
            'pam': 'NGG',
            'position': -1,
            'strand': '',
            'score': 0,
            'type': 'non-targeting'
        })

    essential_genes = ['RPS3', 'RPL11', 'EIF3A', 'POLR2A', 'CDK1']
    for gene in essential_genes[:n_essential]:
        controls.append({
            'gene': gene,
            'gene_id': '',
            'guide_number': 1,
            'sequence': get_validated_guide(gene),
            'pam': 'NGG',
            'position': 0,
            'strand': '+',
            'score': 1,
            'type': 'essential-control'
        })

    nonessential_genes = ['AAVS1', 'ROSA26']
    for gene in nonessential_genes[:n_nonessential]:
        controls.append({
            'gene': gene,
            'gene_id': '',
            'guide_number': 1,
            'sequence': get_validated_guide(gene),
            'pam': 'NGG',
            'position': 0,
            'strand': '+',
            'score': 1,
            'type': 'safe-harbor-control'
        })

    return controls

def generate_nontargeting_sequence(length=20):
    '''Generate random non-targeting sequence.'''
    while True:
        seq = ''.join(np.random.choice(['A', 'C', 'G', 'T'], length))
        gc = (seq.count('G') + seq.count('C')) / length
        if 0.4 <= gc <= 0.6 and 'TTTT' not in seq:
            return seq

def get_validated_guide(gene):
    '''Get validated guide sequence for control gene.'''
    validated = {
        'RPS3': 'GAGCTTCTTCAGCAGCATGG',
        'RPL11': 'GAAACAGGGCATCATCTACG',
        'EIF3A': 'GTGCAAGAGGATGATGACAA',
        'AAVS1': 'GGGGCCACTAGGGACAGGAT',
        'ROSA26': 'GAAGATGGGCGGGAGTCTTC'
    }
    return validated.get(gene, generate_nontargeting_sequence())
```

## Off-Target Analysis

```python
def check_offtargets(guide_sequence, genome_index, max_mismatches=3):
    '''Check for potential off-target sites.'''
    from subprocess import run
    import tempfile

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f'>guide\n{guide_sequence}\n')
        query_file = f.name

    result = run(
        ['bowtie', '-a', '-n', str(max_mismatches), '-l', '20', genome_index, '-f', query_file],
        capture_output=True, text=True
    )

    offtargets = []
    for line in result.stdout.strip().split('\n'):
        if line:
            fields = line.split('\t')
            offtargets.append({
                'chromosome': fields[2],
                'position': int(fields[3]),
                'strand': fields[1],
                'mismatches': int(fields[7]) if len(fields) > 7 else 0
            })

    return offtargets

def filter_by_offtargets(library_df, genome_index, max_offtargets=10):
    '''Filter library to remove guides with too many off-targets.'''
    filtered = []

    for _, guide in library_df.iterrows():
        offtargets = check_offtargets(guide['sequence'], genome_index)
        n_offtargets = len([ot for ot in offtargets if ot['mismatches'] <= 2])

        if n_offtargets <= max_offtargets:
            guide_dict = guide.to_dict()
            guide_dict['n_offtargets'] = n_offtargets
            filtered.append(guide_dict)

    return pd.DataFrame(filtered)
```

## Oligo Design for Cloning

```python
def design_oligos(library_df, vector='lentiGuide-Puro'):
    '''Design oligos for library cloning.'''
    vector_specs = {
        'lentiGuide-Puro': {
            'forward_prefix': 'CACCG',
            'forward_suffix': '',
            'reverse_prefix': 'AAAC',
            'reverse_suffix': 'C'
        },
        'pLKO': {
            'forward_prefix': 'CCGG',
            'forward_suffix': 'CTCGAG',
            'reverse_prefix': 'AATTCTCGAG',
            'reverse_suffix': ''
        }
    }

    spec = vector_specs.get(vector, vector_specs['lentiGuide-Puro'])

    oligos = []
    for _, guide in library_df.iterrows():
        seq = guide['sequence']

        forward = spec['forward_prefix'] + seq + spec['forward_suffix']
        reverse = spec['reverse_prefix'] + str(Seq(seq).reverse_complement()) + spec['reverse_suffix']

        oligos.append({
            'guide_id': f"{guide['gene']}_{guide['guide_number']}",
            'gene': guide['gene'],
            'guide_sequence': seq,
            'forward_oligo': forward,
            'reverse_oligo': reverse,
            'type': guide.get('type', 'targeting')
        })

    return pd.DataFrame(oligos)

oligos = design_oligos(library)
oligos.to_csv('library_oligos.csv', index=False)
print(f'Designed {len(oligos)} oligo pairs')
```

## Pool Design for Synthesis

```python
def design_array_oligos(library_df, array_format='12K'):
    '''Design array oligos for pooled synthesis.'''
    formats = {
        '12K': {'capacity': 12000, 'length': 200},
        '92K': {'capacity': 92000, 'length': 150},
        '244K': {'capacity': 244000, 'length': 60}
    }

    spec = formats[array_format]

    primer_5 = 'AGGCTTGGATTTCTATAACTTCGTATAGCATACATTATACGAAGTTAT'
    primer_3 = 'ATAACTTCGTATAATGTATGCTATACGAAGTTATCTTGGATTTCTAGA'
    scaffold = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'

    array_oligos = []
    for _, guide in library_df.iterrows():
        full_oligo = primer_5 + guide['sequence'] + scaffold + primer_3

        if len(full_oligo) > spec['length']:
            print(f"Warning: {guide['gene']} oligo too long for {array_format}")
            continue

        array_oligos.append({
            'id': f"{guide['gene']}_{guide['guide_number']}",
            'sequence': full_oligo,
            'length': len(full_oligo)
        })

    if len(array_oligos) > spec['capacity']:
        print(f"Warning: Library ({len(array_oligos)}) exceeds {array_format} capacity ({spec['capacity']})")

    return pd.DataFrame(array_oligos)

array_oligos = design_array_oligos(library, '92K')
array_oligos.to_csv('array_synthesis.csv', index=False)
```

## Library QC

```python
def qc_library(library_df):
    '''Quality control checks for library design.'''
    qc = {}

    qc['total_guides'] = len(library_df)
    qc['unique_genes'] = library_df[library_df['type'] == 'targeting']['gene'].nunique()
    qc['guides_per_gene'] = library_df[library_df['type'] == 'targeting'].groupby('gene').size().describe()

    gc_contents = library_df['sequence'].apply(lambda x: (x.count('G') + x.count('C')) / len(x))
    qc['gc_mean'] = gc_contents.mean()
    qc['gc_std'] = gc_contents.std()
    qc['gc_range'] = (gc_contents.min(), gc_contents.max())

    has_poly_t = library_df['sequence'].apply(lambda x: 'TTTT' in x)
    qc['poly_t_count'] = has_poly_t.sum()

    type_counts = library_df['type'].value_counts()
    qc['control_ratio'] = type_counts.get('non-targeting', 0) / len(library_df)

    return qc

qc = qc_library(library)
print('Library QC:')
for key, value in qc.items():
    print(f'  {key}: {value}')
```

## Alternative PAM Systems

The examples above use SpCas9 with NGG PAM. Alternative systems expand targeting range:

| System | PAM | Use Case |
|--------|-----|----------|
| SpCas9 | NGG | Standard, most validated |
| SpCas9-NG | NG | Relaxed PAM requirement |
| SpRY | NRN/NYN | Near-PAMless, broadest targeting |
| Cas12a (Cpf1) | TTTV | AT-rich regions, staggered cuts |
| SaCas9 | NNGRRT | AAV delivery (smaller gene) |

For alternative PAMs, modify the `design_sgrnas_for_gene()` function:

```python
# Cas12a example (TTTV PAM, 23nt guide)
def design_cas12a_guides(gene_sequence, n_guides=4):
    pam_pattern = 'TTT[ACG]'  # TTTV
    guide_length = 23

    for match in re.finditer(f'({pam_pattern})([ACGT]{{{guide_length}}})', gene_sequence):
        pam = match.group(1)
        guide = match.group(2)
        # Cas12a cuts downstream of guide
        # ...
```

## Related Skills

- mageck-analysis - Analyze screen results
- crispresso-editing - Validate editing efficiency
- screen-qc - QC sequencing data
- hit-calling - Identify screen hits
