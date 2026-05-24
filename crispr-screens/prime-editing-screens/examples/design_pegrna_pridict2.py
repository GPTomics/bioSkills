# Reference: PRIDICT2 (uzh-dqbm-cmi/PRIDICT2), pandas 2.2+, biopython 1.83+ | Verify API if version differs
#
# pegRNA library design with PRIDICT2 efficiency prediction.
# PRIDICT2 is invoked via CLI (pridict2_pegRNA_design.py batch).
# Filters to high-efficiency pegRNAs before library synthesis.

import pandas as pd
import subprocess
from pathlib import Path
from Bio.Seq import Seq
import re

# === INPUTS ===
# Intended variants: variant_id, chrom, pos, ref, alt, context
variants_df = pd.read_csv('intended_variants.csv')

# === PRIDICT2 CLI invocation ===
# Build the batch input CSV expected by PRIDICT2: sequence_name, sequence (with (REF/ALT) notation)
# Then call: python pridict2_pegRNA_design.py batch --input-fname variants.csv --output-dir out/ --cores 8 --summarize
# This example uses a placeholder predict_pridict2() that simulates the CLI output.

# === STEP 1: GENERATE pegRNA CANDIDATES ===
def find_pegrna_candidates(chrom, pos, ref, alt, context_seq):
    '''Find pegRNAs that can install this variant.
    Returns list of (spacer, pbs, rtt) candidates.'''
    # PE requires NGG PAM within ~30 nt of edit
    candidates = []
    edit_position_in_context = 30  # variant at position 30 of 60-nt context
    # Scan PAMs in 30-nt upstream + downstream
    for strand in ['+', '-']:
        seq = context_seq if strand == '+' else str(Seq(context_seq).reverse_complement())
        for pam_match in re.finditer(r'(?=([ACGT]GG))', seq):
            pam_pos = pam_match.start()
            cut_pos = pam_pos - 3
            # Distance from cut to edit
            edit_dist = abs(edit_position_in_context - cut_pos)
            if edit_dist > 30:
                continue
            # Spacer = 20 nt upstream of PAM
            if pam_pos < 20:
                continue
            spacer = seq[pam_pos-20:pam_pos]
            # PBS = 11-13 nt downstream of cut (toward edit)
            pbs_length = 12
            pbs_start = cut_pos
            pbs = seq[pbs_start:pbs_start+pbs_length]
            # RTT = encodes intended edit; 10-20 nt
            rtt_length = 15
            # Replace ref with alt in RTT
            rtt_template = list(seq[pbs_start+pbs_length:pbs_start+pbs_length+rtt_length])
            # Position of edit in RTT
            edit_in_rtt = edit_position_in_context - (pbs_start + pbs_length)
            if 0 <= edit_in_rtt < rtt_length:
                rtt_template[edit_in_rtt] = alt
                rtt = ''.join(rtt_template)
            else:
                continue
            candidates.append({
                'spacer': spacer,
                'pbs': pbs,
                'rtt': rtt,
                'pam_strand': strand,
                'edit_dist': edit_dist,
            })
    return candidates

# === STEP 2: PRIDICT2 PREDICTION ===
# In production: write candidates to a CSV with PRIDICT2's batch format
# (sequence_name, sequence with (REF/ALT) notation) and call:
#   python pridict2_pegRNA_design.py batch --input-fname candidates.csv \
#       --output-dir predictions/ --cores 8 --summarize
# Then read predictions/pridict2_summary.csv. The placeholder below simulates.
def predict_pridict2(spacer, pbs, rtt, context):
    '''Placeholder for PRIDICT2 batch CLI output parsing.
    In real use: parse predictions/pridict2_summary.csv after running PRIDICT2 CLI.'''
    pbs_gc = (pbs.count('G') + pbs.count('C')) / len(pbs)
    sim_efficiency = max(0.05, min(0.95, 0.4 + 0.3 * (0.5 - abs(pbs_gc - 0.45))))
    return {'efficiency': sim_efficiency, 'indel_rate': 0.02, 'scaffold_incorp': 0.03}

# === STEP 3: SCREEN ALL CANDIDATES ===
all_pegrnas = []
for _, var in variants_df.iterrows():
    candidates = find_pegrna_candidates(var['chrom'], var['pos'],
                                          var['ref'], var['alt'], var['context'])
    for c in candidates:
        pred = predict_pridict2(c['spacer'], c['pbs'], c['rtt'], var['context'])
        all_pegrnas.append({
            'variant_id': var['variant_id'],
            'spacer': c['spacer'],
            'pbs': c['pbs'],
            'rtt': c['rtt'],
            'pam_strand': c['pam_strand'],
            'edit_dist': c['edit_dist'],
            'pridict2_efficiency': pred['efficiency'],
            'pridict2_indel': pred['indel_rate'],
            'pridict2_scaffold': pred['scaffold_incorp'],
        })

pegrnas_df = pd.DataFrame(all_pegrnas)

# === STEP 4: FILTER AND RANK ===
# Mathis 2024 PRIDICT2 threshold: efficiency >50%
efficient = pegrnas_df[pegrnas_df['pridict2_efficiency'] > 0.50]

# Pick top 3 per variant for library
top3 = (efficient.sort_values(['variant_id', 'pridict2_efficiency'],
                                ascending=[True, False])
                  .groupby('variant_id')
                  .head(3))

# === STEP 5: REPORT ===
n_variants_total = variants_df['variant_id'].nunique()
n_variants_covered = top3['variant_id'].nunique()
print(f'Variants attempted: {n_variants_total}')
print(f'Variants with >=1 pegRNA passing PRIDICT2 >50%: {n_variants_covered}')
print(f'Total pegRNAs (top 3 per variant): {len(top3)}')

# Flag variants without coverage
uncovered = set(variants_df['variant_id']) - set(top3['variant_id'])
print(f'Uncovered variants requiring alternative (BE, SpRY-PE, or accept exclusion): {len(uncovered)}')

# === STEP 6: EXPORT ===
top3.to_csv('peg_library_filtered.csv', index=False)
if uncovered:
    pd.DataFrame({'variant_id': list(uncovered)}).to_csv('uncovered_variants.csv', index=False)

print('\nNext steps:')
print('  1. Pilot top pegRNAs at chromatin-accessible loci first')
print('  2. Synthesize library via oligo pool')
print('  3. Run pooled PE screen at MOI 0.3 in PE-validated cell line')
print('  4. Endpoint amplicon sequencing + CRISPResso2 PE mode')
print('  5. Filter to >5% intended-edit pegRNAs')
print('  6. MAGeCK MLE per-variant hit calling')
