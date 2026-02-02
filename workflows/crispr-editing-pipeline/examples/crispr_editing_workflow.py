#!/usr/bin/env python3
'''
CRISPR editing pipeline: target to delivery-ready constructs.
Supports knockout, base editing, and HDR knockin designs.
Requires: biopython, pandas, numpy, matplotlib, primer3-py
'''
import re
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from pathlib import Path

# Thresholds
ACTIVITY_THRESHOLD = 0.6    # Minimum guide activity score (standard cutoff)
SPECIFICITY_THRESHOLD = 0.7 # Minimum specificity score (off-target penalty)
GC_MIN, GC_MAX = 0.4, 0.7   # Optimal GC content range
HOMOLOGY_ARM_LENGTH = 800   # Standard for plasmid HDR donors
SSODN_ARM_LENGTH = 50       # For single-strand oligo donors

OUTPUT_DIR = Path('crispr_design_results')


def find_all_guides(sequence, pam='NGG'):
    '''Find all potential gRNA target sites with specified PAM.'''
    guides = []
    seq_str = str(sequence).upper()

    # Forward strand: 20bp + NGG
    for match in re.finditer(r'(?=([ATCG]{20}[ATCG]GG))', seq_str):
        pos = match.start()
        target = match.group(1)[:20]
        pam_seq = match.group(1)[20:23]
        guides.append({
            'sequence': target,
            'pam': pam_seq,
            'position': pos,
            'strand': '+',
            'cut_site': pos + 17  # Cas9 cuts 3bp upstream of PAM
        })

    # Reverse strand: CCN + 20bp (complement of NGG)
    for match in re.finditer(r'(?=(CC[ATCG][ATCG]{20}))', seq_str):
        pos = match.start()
        full = match.group(1)
        target = str(Seq(full[3:23]).reverse_complement())
        pam_seq = str(Seq(full[0:3]).reverse_complement())
        guides.append({
            'sequence': target,
            'pam': pam_seq,
            'position': pos,
            'strand': '-',
            'cut_site': pos + 6
        })

    return pd.DataFrame(guides)


def calculate_activity_score(guide_seq):
    '''Calculate guide activity score using Rule Set 2-like heuristics.'''
    score = 0.5

    # GC content (optimal: 40-70%)
    gc = (guide_seq.count('G') + guide_seq.count('C')) / len(guide_seq)
    if GC_MIN <= gc <= GC_MAX:
        score += 0.2
    elif gc < 0.3 or gc > 0.8:
        score -= 0.2

    # Poly-T terminator (>4 T's terminates Pol III)
    if 'TTTT' in guide_seq:
        score -= 0.3

    # G at position 20 (PAM-proximal) preferred
    if guide_seq[-1] == 'G':
        score += 0.1

    # Avoid GG at positions 19-20 (can cause issues)
    if guide_seq[-2:] == 'GG':
        score -= 0.1

    # Seed region GC (positions 12-20)
    seed = guide_seq[11:20]
    seed_gc = (seed.count('G') + seed.count('C')) / len(seed)
    if 0.4 <= seed_gc <= 0.7:
        score += 0.1

    # Purine at position 20 slightly preferred
    if guide_seq[-1] in 'AG':
        score += 0.05

    return min(1.0, max(0.0, score))


def calculate_specificity_score(guide_seq, offtargets=None):
    '''Calculate specificity score (simplified without genome search).'''
    # In production, this would use Cas-OFFinder results
    # Here we use sequence-based heuristics
    score = 1.0

    # Penalize low-complexity sequences
    if len(set(guide_seq)) < 4:
        score -= 0.3

    # Penalize homopolymers
    for base in 'ATCG':
        if base * 4 in guide_seq:
            score -= 0.1

    # Penalize repetitive dinucleotides
    for di in ['AT', 'TA', 'GC', 'CG']:
        if (di * 3) in guide_seq:
            score -= 0.1

    return max(0.0, score)


def design_knockout(guides_df, target_sequence, output_dir):
    '''Design knockout experiment with validation primers.'''
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    best_guide = guides_df.iloc[0]

    design = {
        'edit_type': 'knockout',
        'guide_sequence': best_guide['sequence'],
        'pam': best_guide['pam'],
        'position': int(best_guide['position']),
        'strand': best_guide['strand'],
        'cut_site': int(best_guide['cut_site']),
        'activity_score': float(best_guide['activity_score']),
        'specificity_score': float(best_guide['specificity_score']),
        'expected_outcome': 'Frameshift indel leading to nonsense-mediated decay',
        'validation_strategy': 'PCR + gel/Sanger sequencing of ~300bp amplicon'
    }

    # Validation amplicon coordinates (~300bp centered on cut site)
    cut = design['cut_site']
    design['validation_amplicon'] = {
        'start': max(0, cut - 150),
        'end': min(len(target_sequence), cut + 150),
        'expected_wt_size': 300,
        'note': 'Indels will cause size shift on gel or mixed peaks in Sanger'
    }

    with open(output_dir / 'knockout_design.json', 'w') as f:
        json.dump(design, f, indent=2)

    return design


def design_base_edit(guides_df, target_sequence, target_position, edit_type='CBE', output_dir='.'):
    '''Design base editing experiment.

    CBE editing window: positions 4-8 (counting from PAM-distal end)
    ABE editing window: positions 4-7 (slightly narrower)
    '''
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    target_base = target_sequence[target_position].upper()
    compatible_guides = []

    for _, guide in guides_df.iterrows():
        # Calculate position of target base within protospacer
        if guide['strand'] == '+':
            pos_in_guide = target_position - guide['position']
        else:
            pos_in_guide = guide['position'] + 20 - target_position

        # Check if in editing window (4-8 for CBE, 4-7 for ABE)
        window_start, window_end = (4, 8) if edit_type == 'CBE' else (4, 7)
        if window_start <= pos_in_guide <= window_end:
            # Check base compatibility
            if edit_type == 'CBE' and target_base in ['C', 'G']:
                guide_copy = guide.copy()
                guide_copy['target_position_in_guide'] = pos_in_guide
                compatible_guides.append(guide_copy)
            elif edit_type == 'ABE' and target_base in ['A', 'T']:
                guide_copy = guide.copy()
                guide_copy['target_position_in_guide'] = pos_in_guide
                compatible_guides.append(guide_copy)

    if not compatible_guides:
        return {'error': 'No compatible guides found for target position'}

    compatible_df = pd.DataFrame(compatible_guides)
    best_guide = compatible_df.sort_values('activity_score', ascending=False).iloc[0]

    design = {
        'edit_type': edit_type,
        'guide_sequence': best_guide['sequence'],
        'pam': best_guide['pam'],
        'target_position': int(target_position),
        'target_base': target_base,
        'position_in_guide': int(best_guide['target_position_in_guide']),
        'editing_window': f'positions {window_start}-{window_end}',
        'expected_edit': f"{target_base}>{('T' if edit_type == 'CBE' else 'G')}",
        'activity_score': float(best_guide['activity_score']),
        'bystander_warning': 'Check for other C (CBE) or A (ABE) bases in editing window'
    }

    # Check for bystander edits
    guide_seq = best_guide['sequence']
    target_bases = 'C' if edit_type == 'CBE' else 'A'
    bystanders = []
    for i in range(window_start, window_end + 1):
        if i != best_guide['target_position_in_guide']:
            if guide_seq[i-1] == target_bases:
                bystanders.append(i)

    design['potential_bystander_positions'] = bystanders

    with open(output_dir / f'{edit_type.lower()}_design.json', 'w') as f:
        json.dump(design, f, indent=2)

    return design


def design_hdr_knockin(guides_df, target_sequence, insert_sequence, output_dir='.', arm_length=800):
    '''Design HDR knockin with donor template.'''
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    best_guide = guides_df.iloc[0]
    cut_site = int(best_guide['cut_site'])

    # Extract homology arms
    left_start = max(0, cut_site - arm_length)
    left_arm = target_sequence[left_start:cut_site]
    right_arm = target_sequence[cut_site:cut_site + arm_length]

    # Mutate PAM in donor to prevent re-cutting
    # For right arm, PAM is ~3bp after cut site
    pam_pos_in_right = 3
    right_arm_mutated = right_arm[:pam_pos_in_right] + 'A' + right_arm[pam_pos_in_right+1:]  # NGG -> NAG

    # Assemble donor
    donor_template = left_arm + insert_sequence + right_arm_mutated

    design = {
        'edit_type': 'HDR_knockin',
        'guide_sequence': best_guide['sequence'],
        'pam': best_guide['pam'],
        'cut_site': cut_site,
        'activity_score': float(best_guide['activity_score']),
        'left_arm_length': len(left_arm),
        'right_arm_length': len(right_arm_mutated),
        'insert_length': len(insert_sequence),
        'total_donor_length': len(donor_template),
        'pam_mutation': 'NGG -> NAG (prevents re-cutting)',
        'delivery_recommendation': 'plasmid' if len(donor_template) > 200 else 'ssODN'
    }

    # Save donor template
    with open(output_dir / 'hdr_donor_template.fasta', 'w') as f:
        f.write(f">HDR_donor_template length={len(donor_template)}\n")
        for i in range(0, len(donor_template), 80):
            f.write(donor_template[i:i+80] + '\n')

    with open(output_dir / 'hdr_knockin_design.json', 'w') as f:
        json.dump(design, f, indent=2)

    return design


def visualize_guides(guides_df, target_length, exon_coords=None, output_path='guide_landscape.pdf'):
    '''Visualize guide positions and scores along target.'''
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), gridspec_kw={'height_ratios': [1, 2]})

    # Top panel: Target structure
    ax1 = axes[0]
    ax1.axhline(y=0.5, color='lightgray', linewidth=15, solid_capstyle='butt')

    if exon_coords:
        for start, end in exon_coords:
            ax1.axhline(y=0.5, xmin=start/target_length, xmax=end/target_length,
                       color='steelblue', linewidth=25, solid_capstyle='butt')

    ax1.set_xlim(0, target_length)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel('Target')
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)

    # Bottom panel: Guide scores
    ax2 = axes[1]
    colors = ['#2ecc71' if s >= ACTIVITY_THRESHOLD else '#f39c12' if s >= 0.4 else '#e74c3c'
              for s in guides_df['activity_score']]

    ax2.scatter(guides_df['position'], guides_df['activity_score'], c=colors, s=60, alpha=0.7, edgecolor='white')
    ax2.axhline(y=ACTIVITY_THRESHOLD, color='#2ecc71', linestyle='--', alpha=0.7, label=f'Threshold ({ACTIVITY_THRESHOLD})')
    ax2.set_xlim(0, target_length)
    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel('Position (bp)', fontsize=12)
    ax2.set_ylabel('Activity Score', fontsize=12)
    ax2.legend(loc='lower right')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f'Saved guide landscape to {output_path}')
    return fig


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Example target sequence (BRCA1 exon 2 region)
    target_sequence = '''ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTG'''.replace('\n', '')

    print('=== CRISPR Editing Pipeline ===')
    print(f'Target length: {len(target_sequence)} bp')

    # Step 1: Find all guides
    print('\n=== Finding guide RNAs ===')
    guides = find_all_guides(target_sequence)
    print(f'Found {len(guides)} potential guide sites')

    # Step 2: Score guides
    print('\n=== Scoring guides ===')
    guides['activity_score'] = guides['sequence'].apply(calculate_activity_score)
    guides['specificity_score'] = guides['sequence'].apply(calculate_specificity_score)
    guides['combined_score'] = guides['activity_score'] * 0.6 + guides['specificity_score'] * 0.4

    # Filter and rank
    good_guides = guides[guides['activity_score'] >= ACTIVITY_THRESHOLD].copy()
    good_guides = good_guides.sort_values('combined_score', ascending=False)

    print(f'High-scoring guides (activity >= {ACTIVITY_THRESHOLD}): {len(good_guides)}')

    # Save all guides
    guides.to_csv(OUTPUT_DIR / 'all_guides.tsv', sep='\t', index=False)
    good_guides.to_csv(OUTPUT_DIR / 'filtered_guides.tsv', sep='\t', index=False)

    # Step 3: Design experiments
    print('\n=== Designing experiments ===')

    # 3a: Knockout
    print('\n--- Knockout Design ---')
    ko_design = design_knockout(good_guides, target_sequence, OUTPUT_DIR)
    print(f"Guide: {ko_design['guide_sequence']}")
    print(f"Activity score: {ko_design['activity_score']:.2f}")

    # 3b: Base editing (find a C for CBE)
    print('\n--- Base Editing Design (CBE) ---')
    c_positions = [i for i, base in enumerate(target_sequence) if base == 'C']
    if c_positions:
        target_c = c_positions[50]  # Pick a C in middle of sequence
        cbe_design = design_base_edit(good_guides, target_sequence, target_c, 'CBE', OUTPUT_DIR)
        if 'error' not in cbe_design:
            print(f"Guide: {cbe_design['guide_sequence']}")
            print(f"Target position: {cbe_design['target_position']}")
            print(f"Position in guide: {cbe_design['position_in_guide']}")
        else:
            print(cbe_design['error'])

    # 3c: HDR knockin
    print('\n--- HDR Knockin Design ---')
    flag_tag = 'GATTACAAGGATGACGACGATAAG'  # FLAG tag sequence
    hdr_design = design_hdr_knockin(good_guides, target_sequence, flag_tag, OUTPUT_DIR, arm_length=100)
    print(f"Guide: {hdr_design['guide_sequence']}")
    print(f"Donor length: {hdr_design['total_donor_length']} bp")
    print(f"Delivery: {hdr_design['delivery_recommendation']}")

    # Visualization
    print('\n=== Generating visualizations ===')
    visualize_guides(guides, len(target_sequence),
                    exon_coords=[(0, 150), (200, 350), (400, len(target_sequence))],
                    output_path=OUTPUT_DIR / 'guide_landscape.pdf')

    print('\n=== Pipeline Complete ===')
    print(f'Results saved to: {OUTPUT_DIR}/')


if __name__ == '__main__':
    main()
