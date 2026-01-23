'''miRNA target prediction using multiple approaches'''

import pandas as pd
from Bio.Seq import Seq

def find_seed_matches(mirna_seq, utr_seq, seed_type='7mer-m8'):
    '''Find miRNA seed matches in UTR sequence

    Seed types:
    - 8mer: Perfect 8nt match (positions 1-8)
    - 7mer-m8: Match at positions 2-8 plus position 1 A
    - 7mer-A1: Match at positions 2-8 plus A at position 1
    - 6mer: Match at positions 2-7

    Position 1 = 5' end of miRNA
    '''
    mirna = Seq(mirna_seq)
    utr = Seq(utr_seq.upper())

    # Seed is positions 2-8 (0-indexed: 1-8)
    seed_7mer = str(mirna[1:8])
    seed_7mer_rc = str(Seq(seed_7mer).reverse_complement())

    matches = []
    start = 0
    while True:
        pos = str(utr).find(seed_7mer_rc, start)
        if pos == -1:
            break

        # Check for 8mer (A at target position 1)
        if pos > 0 and str(utr)[pos - 1] == 'A':
            match_type = '8mer' if str(utr)[pos + 7:pos + 8] == str(Seq(str(mirna[0])).complement()) else '7mer-A1'
        else:
            match_type = '7mer-m8'

        matches.append({
            'position': pos,
            'match_type': match_type,
            'seed_match': seed_7mer_rc
        })
        start = pos + 1

    return matches

def parse_miranda_output(filepath):
    '''Parse miRanda prediction output

    miRanda scores:
    - Score >= 140: Default threshold
    - Energy <= -20 kcal/mol: Stable duplex
    Higher score + lower energy = stronger prediction
    '''
    results = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    results.append({
                        'mirna': parts[0].lstrip('>'),
                        'target': parts[1],
                        'score': float(parts[2]),
                        'energy': float(parts[3])
                    })
    df = pd.DataFrame(results)

    # Filter by default thresholds
    SCORE_THRESHOLD = 140
    ENERGY_THRESHOLD = -20
    filtered = df[(df['score'] >= SCORE_THRESHOLD) & (df['energy'] <= ENERGY_THRESHOLD)]
    return filtered

def load_targetscan(filepath, mirna_family):
    '''Load TargetScan predictions for a miRNA family

    TargetScan context++ score:
    - < -0.4: Very high confidence
    - < -0.2: High confidence
    - < 0: Predicted target
    '''
    df = pd.read_csv(filepath, sep='\t')
    targets = df[df['miRNA family'] == mirna_family].copy()
    targets = targets.sort_values('context++ score')
    return targets

def load_mirdb(filepath, mirna_id, score_threshold=80):
    '''Load miRDB predictions

    miRDB score interpretation:
    - >= 80: High confidence (recommended)
    - 60-80: Medium confidence
    - < 60: Low confidence
    '''
    df = pd.read_csv(filepath, sep='\t', header=None,
                     names=['mirna', 'target_gene', 'score'])
    targets = df[(df['mirna'] == mirna_id) & (df['score'] >= score_threshold)]
    return targets.sort_values('score', ascending=False)

def find_consensus_targets(miranda_targets, targetscan_targets, mirdb_targets, min_databases=2):
    '''Find targets predicted by multiple databases

    Consensus targets (2+ databases) are more reliable
    3/3 consensus = high confidence validated computationally
    '''
    all_targets = set(miranda_targets) | set(targetscan_targets) | set(mirdb_targets)

    consensus = []
    for target in all_targets:
        in_miranda = target in miranda_targets
        in_targetscan = target in targetscan_targets
        in_mirdb = target in mirdb_targets
        n_databases = sum([in_miranda, in_targetscan, in_mirdb])

        if n_databases >= min_databases:
            consensus.append({
                'target': target,
                'n_databases': n_databases,
                'miranda': in_miranda,
                'targetscan': in_targetscan,
                'mirdb': in_mirdb
            })

    return pd.DataFrame(consensus).sort_values('n_databases', ascending=False)

# Example usage
if __name__ == '__main__':
    # Example miRNA sequence (hsa-miR-21-5p)
    mirna = 'UAGCUUAUCAGACUGAUGUUGA'

    # Example 3' UTR with seed match
    utr = 'ACGTACGTATAAGCTATCGTACGT'

    matches = find_seed_matches(mirna, utr)
    print(f'Seed matches found: {len(matches)}')
    for m in matches:
        print(f"  Position {m['position']}: {m['match_type']}")
