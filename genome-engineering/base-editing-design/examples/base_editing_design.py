'''Design guides for cytosine and adenine base editing'''

from Bio.Seq import Seq
import re

# Editing window positions (1-indexed from PAM-distal end)
# For 20nt spacer: position 1 is 5' end, position 20 is adjacent to PAM
# Editing window for BE4max (CBE): positions 4-8
# Editing window for ABE8e: positions 4-7 (narrower than CBE)
CBE_WINDOW = (4, 8)
ABE_WINDOW = (4, 7)

# Position-dependent editing efficiency
# Based on published characterization data
# 1.0 = maximum efficiency at that position
CBE_POSITION_EFFICIENCY = {
    1: 0.05, 2: 0.10, 3: 0.20,
    4: 0.70, 5: 0.90, 6: 1.00,  # Peak at position 6
    7: 0.85, 8: 0.50,
    9: 0.20, 10: 0.10
}

ABE_POSITION_EFFICIENCY = {
    1: 0.02, 2: 0.05, 3: 0.15,
    4: 0.60, 5: 0.95, 6: 1.00,  # Peak at positions 5-6
    7: 0.70,
    8: 0.20, 9: 0.05
}

# Sequence context preferences for CBE
# 5' neighbor of target C affects efficiency
CBE_CONTEXT_SCORES = {
    'T': 1.0,   # TC: most preferred
    'C': 0.8,   # CC: good
    'A': 0.6,   # AC: moderate
    'G': 0.4,   # GC: less efficient
}


def find_cbe_targets(sequence, target_c_position):
    '''Find guides that place a C in the CBE editing window

    Args:
        sequence: DNA sequence containing the target C
        target_c_position: 0-indexed position of C to edit

    Returns:
        List of guide options sorted by bystander count
    '''
    sequence = sequence.upper()
    guides = []

    for pam_match in re.finditer(r'(?=(.GG))', sequence):
        pam_pos = pam_match.start()
        spacer_start = pam_pos - 20
        if spacer_start < 0:
            continue

        # Position in spacer (1-indexed from 5' end)
        c_pos_in_spacer = target_c_position - spacer_start + 1

        if CBE_WINDOW[0] <= c_pos_in_spacer <= CBE_WINDOW[1]:
            spacer = sequence[spacer_start:pam_pos]

            # Find bystander Cs
            bystanders = []
            for i in range(CBE_WINDOW[0] - 1, CBE_WINDOW[1]):
                if i < len(spacer) and spacer[i] == 'C':
                    if (spacer_start + i) != target_c_position:
                        bystanders.append(i + 1)

            # Score sequence context
            context_score = 0.5
            if c_pos_in_spacer > 1:
                neighbor = spacer[c_pos_in_spacer - 2]  # 5' neighbor
                context_score = CBE_CONTEXT_SCORES.get(neighbor, 0.5)

            guides.append({
                'spacer': spacer,
                'pam': sequence[pam_pos:pam_pos + 3],
                'target_position': c_pos_in_spacer,
                'efficiency': CBE_POSITION_EFFICIENCY.get(c_pos_in_spacer, 0.1),
                'context_score': context_score,
                'bystander_positions': bystanders,
                'bystander_count': len(bystanders),
                'strand': '+'
            })

    return sorted(guides, key=lambda x: (x['bystander_count'], -x['efficiency']))


def find_abe_targets(sequence, target_a_position):
    '''Find guides that place an A in the ABE editing window'''
    sequence = sequence.upper()
    guides = []

    for pam_match in re.finditer(r'(?=(.GG))', sequence):
        pam_pos = pam_match.start()
        spacer_start = pam_pos - 20
        if spacer_start < 0:
            continue

        a_pos_in_spacer = target_a_position - spacer_start + 1

        if ABE_WINDOW[0] <= a_pos_in_spacer <= ABE_WINDOW[1]:
            spacer = sequence[spacer_start:pam_pos]

            bystanders = []
            for i in range(ABE_WINDOW[0] - 1, ABE_WINDOW[1]):
                if i < len(spacer) and spacer[i] == 'A':
                    if (spacer_start + i) != target_a_position:
                        bystanders.append(i + 1)

            guides.append({
                'spacer': spacer,
                'pam': sequence[pam_pos:pam_pos + 3],
                'target_position': a_pos_in_spacer,
                'efficiency': ABE_POSITION_EFFICIENCY.get(a_pos_in_spacer, 0.1),
                'bystander_positions': bystanders,
                'bystander_count': len(bystanders),
                'strand': '+'
            })

    return sorted(guides, key=lambda x: (x['bystander_count'], -x['efficiency']))


def predict_editing_outcome(spacer, editor='CBE'):
    '''Predict all bases that will be edited in the window

    Returns list of predicted edits with efficiency scores
    '''
    edits = []
    window = CBE_WINDOW if editor == 'CBE' else ABE_WINDOW
    target_base = 'C' if editor == 'CBE' else 'A'
    result_base = 'T' if editor == 'CBE' else 'G'
    efficiency_map = CBE_POSITION_EFFICIENCY if editor == 'CBE' else ABE_POSITION_EFFICIENCY

    for i in range(window[0] - 1, window[1]):
        if i < len(spacer) and spacer[i] == target_base:
            pos = i + 1
            edits.append({
                'position': pos,
                'change': f'{target_base}>{result_base}',
                'efficiency': efficiency_map.get(pos, 0.1)
            })

    return edits


if __name__ == '__main__':
    # Example: Design CBE guide for a C>T conversion
    target_seq = 'ATCGATCGATCGATCGATCGCATCGATCGATCGATCGATCGATCGATCGAGG'
    target_c_pos = 20  # Position of C to edit

    print(f'Target sequence: {target_seq}')
    print(f'Target C position: {target_c_pos} (base: {target_seq[target_c_pos]})')
    print()

    # Find CBE guides
    cbe_guides = find_cbe_targets(target_seq, target_c_pos)

    print(f'Found {len(cbe_guides)} CBE guide options:')
    for i, g in enumerate(cbe_guides[:3], 1):
        print(f"\n{i}. Spacer: {g['spacer']}")
        print(f"   Target at position {g['target_position']} (efficiency: {g['efficiency']:.0%})")
        print(f"   Context score: {g['context_score']:.1f}")
        print(f"   Bystander Cs: {g['bystander_positions'] if g['bystander_positions'] else 'None'}")

        # Predict editing outcome
        edits = predict_editing_outcome(g['spacer'], 'CBE')
        print(f'   Predicted edits:')
        for e in edits:
            is_target = '(TARGET)' if e['position'] == g['target_position'] else '(bystander)'
            print(f"      Position {e['position']}: {e['change']} ({e['efficiency']:.0%}) {is_target}")
