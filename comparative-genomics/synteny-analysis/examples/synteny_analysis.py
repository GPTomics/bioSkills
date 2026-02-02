'''Synteny analysis with MCScanX'''

import subprocess
import pandas as pd
from collections import defaultdict


def prepare_gff_for_mcscanx(gff_file, output_file, species_prefix):
    '''Convert GFF3 to MCScanX format

    MCScanX format: species  gene_id  chromosome  start  end
    Tab-separated, sorted by chromosome and position
    '''
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] == 'gene':
                chrom = parts[0]
                start, end = int(parts[3]), int(parts[4])
                attrs = parts[8]
                gene_id = None
                for attr in attrs.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1]
                        break
                if gene_id:
                    genes.append({
                        'species': species_prefix,
                        'gene': gene_id,
                        'chrom': f'{species_prefix}{chrom}',
                        'start': start,
                        'end': end
                    })

    df = pd.DataFrame(genes)
    df = df.sort_values(['chrom', 'start'])
    df.to_csv(output_file, sep='\t', index=False, header=False)
    print(f'Wrote {len(df)} genes to {output_file}')

    return output_file


def run_blastp(fasta1, fasta2, output_file, evalue=1e-10, num_threads=4):
    '''Run all-vs-all BLASTP for homology detection

    Parameters:
    - evalue: 1e-10 for closely related species, 1e-5 for distant
    - outfmt 6: Tabular format required by MCScanX
    '''
    # Make combined database
    db_name = 'combined_db'
    subprocess.run(f'cat {fasta1} {fasta2} > combined.faa', shell=True)
    subprocess.run(f'makeblastdb -in combined.faa -dbtype prot -out {db_name}', shell=True)

    # Run BLASTP
    # -max_target_seqs 5: Top 5 hits per query (balance sensitivity vs speed)
    # -outfmt 6: Tabular output for MCScanX
    cmd = f'blastp -query combined.faa -db {db_name} -out {output_file} -evalue {evalue} -num_threads {num_threads} -max_target_seqs 5 -outfmt 6'
    subprocess.run(cmd, shell=True)

    return output_file


def run_mcscanx(prefix, min_genes=5, max_gaps=25):
    '''Run MCScanX synteny detection

    Parameters:
    - min_genes (default 5): Minimum genes per syntenic block
      Lower values (3) find more blocks but include noise
      Higher values (10) give high-confidence blocks only
    - max_gaps (default 25): Maximum intervening genes allowed
      Larger values accommodate gene loss/rearrangement
    '''
    cmd = f'MCScanX -s {min_genes} -m {max_gaps} {prefix}'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f'MCScanX error: {result.stderr}')
        return None

    return f'{prefix}.collinearity'


def parse_collinearity_file(collinearity_file):
    '''Parse MCScanX output into structured data'''
    blocks = []
    current_block = None

    with open(collinearity_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('## Alignment'):
                if current_block:
                    blocks.append(current_block)
                parts = line.split()
                block_num = int(parts[2].rstrip(':'))
                score = int(parts[3].split('=')[1])
                e_value = float(parts[4].split('=')[1])
                n_genes = int(parts[5])
                plus_minus = parts[6] if len(parts) > 6 else 'plus'
                current_block = {
                    'block_id': block_num,
                    'score': score,
                    'e_value': e_value,
                    'n_genes': n_genes,
                    'orientation': plus_minus,
                    'gene_pairs': []
                }
            elif current_block and line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 3:
                    gene1 = parts[1]
                    gene2 = parts[2]
                    current_block['gene_pairs'].append((gene1, gene2))

    if current_block:
        blocks.append(current_block)

    return blocks


def summarize_synteny(blocks):
    '''Summarize synteny analysis results'''
    print('Synteny Analysis Summary')
    print('=' * 50)

    total_blocks = len(blocks)
    total_genes = sum(b['n_genes'] for b in blocks)
    avg_genes = total_genes / total_blocks if total_blocks > 0 else 0

    print(f'\nSyntenic blocks found: {total_blocks}')
    print(f'Total gene pairs: {total_genes}')
    print(f'Average genes per block: {avg_genes:.1f}')

    # Size distribution
    sizes = [b['n_genes'] for b in blocks]
    print(f'\nBlock size range: {min(sizes)} - {max(sizes)} genes')

    # Orientation
    plus = sum(1 for b in blocks if b['orientation'] == 'plus')
    minus = total_blocks - plus
    print(f'\nOrientation: {plus} plus, {minus} minus (inverted)')

    # Top blocks by size
    print('\nLargest syntenic blocks:')
    sorted_blocks = sorted(blocks, key=lambda x: x['n_genes'], reverse=True)
    for i, block in enumerate(sorted_blocks[:5], 1):
        print(f"  {i}. Block {block['block_id']}: {block['n_genes']} genes "
              f"(score={block['score']}, {block['orientation']})")

    return {
        'total_blocks': total_blocks,
        'total_gene_pairs': total_genes,
        'avg_block_size': avg_genes
    }


def detect_wgd_signature(blocks, ks_values=None):
    '''Detect whole-genome duplication signatures

    WGD indicators:
    - Multiple large syntenic blocks with same Ks
    - 2:1 or 4:1 chromosome mapping ratios
    - Parallel diagonal lines in dot plot
    '''
    # Large blocks (>20 genes) suggest major conserved regions
    # Multiple such blocks between same chromosome pairs indicate WGD
    large_blocks = [b for b in blocks if b['n_genes'] >= 20]
    print(f'\nLarge blocks (>20 genes): {len(large_blocks)}')

    if ks_values:
        import numpy as np
        from scipy.stats import gaussian_kde

        # Ks peak detection for WGD dating
        # Peak at Ks ~0.8-1.2 common for ancient WGD
        ks_clean = [k for k in ks_values if 0.01 < k < 2.0]
        if len(ks_clean) > 50:
            kde = gaussian_kde(ks_clean)
            x = np.linspace(0.01, 2.0, 200)
            y = kde(x)
            peaks = x[np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1]
            print(f'Ks peaks detected at: {peaks}')


if __name__ == '__main__':
    print('Synteny Analysis with MCScanX')
    print('=' * 50)

    # Example: Simulated blocks for demonstration
    example_blocks = [
        {'block_id': 1, 'score': 2500, 'e_value': 0, 'n_genes': 50,
         'orientation': 'plus', 'gene_pairs': [('At1g01010', 'Os1g01010')]},
        {'block_id': 2, 'score': 1800, 'e_value': 0, 'n_genes': 36,
         'orientation': 'minus', 'gene_pairs': [('At2g01010', 'Os3g01010')]},
        {'block_id': 3, 'score': 1200, 'e_value': 0, 'n_genes': 24,
         'orientation': 'plus', 'gene_pairs': [('At3g01010', 'Os2g01010')]},
        {'block_id': 4, 'score': 800, 'e_value': 0, 'n_genes': 16,
         'orientation': 'plus', 'gene_pairs': [('At4g01010', 'Os4g01010')]},
        {'block_id': 5, 'score': 400, 'e_value': 1e-100, 'n_genes': 8,
         'orientation': 'minus', 'gene_pairs': [('At5g01010', 'Os5g01010')]},
    ]

    summarize_synteny(example_blocks)
    detect_wgd_signature(example_blocks)

    print('\n\nTo run on real data:')
    print('1. prepare_gff_for_mcscanx(gff1, "species1.gff", "sp1")')
    print('2. prepare_gff_for_mcscanx(gff2, "species2.gff", "sp2")')
    print('3. run_blastp(pep1.fa, pep2.fa, "combined.blast")')
    print('4. Combine: cat species1.gff species2.gff > combined.gff')
    print('5. run_mcscanx("combined")')
    print('6. blocks = parse_collinearity_file("combined.collinearity")')
