import numpy as np
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq

# === CONFIGURATION ===
output_dir = Path('library_design/')
output_dir.mkdir(exist_ok=True)

np.random.seed(42)

# === 1. DEFINE TARGET GENES ===
print('Defining target gene list...')

target_genes = ['TP53', 'BRCA1', 'BRCA2', 'KRAS', 'NRAS', 'BRAF', 'MYC', 'MYCN', 'CDK4', 'CDK6',
                'RB1', 'PTEN', 'PIK3CA', 'AKT1', 'MTOR', 'EGFR', 'ERBB2', 'MET', 'ALK', 'ROS1']

print(f'Target genes: {len(target_genes)}')

# === 2. GENERATE MOCK sgRNA CANDIDATES ===
print('\nDesigning sgRNAs for each gene...')

def generate_sgrna_sequence(length=20):
    bases = ['A', 'C', 'G', 'T']
    seq = ''.join(np.random.choice(bases, length))
    return seq

def score_sgrna(sequence):
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    gc_score = 1 - abs(gc_content - 0.5) * 2

    poly_t_penalty = 0 if 'TTTT' in sequence else 1
    start_g_bonus = 1 if sequence.startswith('G') else 0.8

    return gc_score * poly_t_penalty * start_g_bonus, gc_content

guides_per_gene = 4
all_guides = []

for gene in target_genes:
    candidates = []
    for _ in range(50):
        seq = generate_sgrna_sequence()
        score, gc = score_sgrna(seq)
        candidates.append({'sequence': seq, 'score': score, 'gc_content': gc})

    candidates_df = pd.DataFrame(candidates)
    candidates_df = candidates_df.sort_values('score', ascending=False)
    top_guides = candidates_df.head(guides_per_gene)

    for i, (_, guide) in enumerate(top_guides.iterrows()):
        all_guides.append({
            'gene': gene,
            'guide_number': i + 1,
            'sequence': guide['sequence'],
            'score': guide['score'],
            'gc_content': guide['gc_content'],
            'type': 'targeting'
        })

print(f'Targeting guides: {len(all_guides)}')

# === 3. ADD CONTROLS ===
print('\nAdding control guides...')

n_nontargeting = 50
for i in range(n_nontargeting):
    while True:
        seq = generate_sgrna_sequence()
        score, gc = score_sgrna(seq)
        if score > 0.5:
            break

    all_guides.append({
        'gene': f'NonTargeting_{i+1:03d}',
        'guide_number': 1,
        'sequence': seq,
        'score': 0,
        'gc_content': gc,
        'type': 'non-targeting'
    })

essential_genes = ['RPL11', 'RPS3', 'EIF3A', 'POLR2A', 'CDK1', 'SF3B1', 'U2AF1', 'PRPF8', 'SNRPD1', 'SNRPE']
for gene in essential_genes:
    seq = generate_sgrna_sequence()
    _, gc = score_sgrna(seq)
    all_guides.append({
        'gene': gene,
        'guide_number': 1,
        'sequence': seq,
        'score': 1.0,
        'gc_content': gc,
        'type': 'essential-control'
    })

for gene in ['AAVS1', 'ROSA26']:
    seq = generate_sgrna_sequence()
    _, gc = score_sgrna(seq)
    all_guides.append({
        'gene': gene,
        'guide_number': 1,
        'sequence': seq,
        'score': 1.0,
        'gc_content': gc,
        'type': 'safe-harbor'
    })

library_df = pd.DataFrame(all_guides)
print(f'Total library size: {len(library_df)}')

# === 4. DESIGN OLIGOS ===
print('\nDesigning cloning oligos...')

forward_prefix = 'CACCG'
reverse_prefix = 'AAAC'
reverse_suffix = 'C'

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

library_df['forward_oligo'] = forward_prefix + library_df['sequence']
library_df['reverse_oligo'] = reverse_prefix + library_df['sequence'].apply(reverse_complement) + reverse_suffix

# === 5. LIBRARY QC ===
print('\n=== LIBRARY QC ===')

print(f"\nGuide type distribution:")
print(library_df['type'].value_counts())

targeting = library_df[library_df['type'] == 'targeting']
print(f"\nTargeting guides per gene:")
guides_per = targeting.groupby('gene').size()
print(f"  Mean: {guides_per.mean():.1f}")
print(f"  Min: {guides_per.min()}")
print(f"  Max: {guides_per.max()}")

print(f"\nGC content distribution:")
print(f"  Mean: {library_df['gc_content'].mean():.1%}")
print(f"  Std: {library_df['gc_content'].std():.1%}")
print(f"  Range: {library_df['gc_content'].min():.1%} - {library_df['gc_content'].max():.1%}")

poly_t_count = library_df['sequence'].apply(lambda x: 'TTTT' in x).sum()
print(f"\nPoly-T sequences: {poly_t_count} ({poly_t_count/len(library_df):.1%})")

control_pct = (library_df['type'] != 'targeting').sum() / len(library_df)
print(f"Control fraction: {control_pct:.1%}")

# === 6. EXPORT ===
print('\n=== EXPORTING LIBRARY ===')

library_df.to_csv(output_dir / 'library_design.csv', index=False)

oligo_order = library_df[['gene', 'guide_number', 'forward_oligo', 'reverse_oligo', 'type']].copy()
oligo_order['oligo_id'] = library_df['gene'] + '_g' + library_df['guide_number'].astype(str)
oligo_order.to_csv(output_dir / 'oligo_order.csv', index=False)

summary = {
    'total_guides': len(library_df),
    'targeting_genes': library_df[library_df['type'] == 'targeting']['gene'].nunique(),
    'targeting_guides': (library_df['type'] == 'targeting').sum(),
    'nontargeting_controls': (library_df['type'] == 'non-targeting').sum(),
    'essential_controls': (library_df['type'] == 'essential-control').sum(),
    'safeharbor_controls': (library_df['type'] == 'safe-harbor').sum(),
    'mean_gc': library_df['gc_content'].mean(),
    'guides_per_gene': guides_per_gene
}

pd.DataFrame([summary]).to_csv(output_dir / 'library_summary.csv', index=False)

print(f"\nLibrary design saved to {output_dir}/")
print(f"  - library_design.csv: Full library with scores")
print(f"  - oligo_order.csv: Oligo ordering format")
print(f"  - library_summary.csv: Summary statistics")
