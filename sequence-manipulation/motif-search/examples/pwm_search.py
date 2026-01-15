'''Probabilistic motif searching with PWM/PSSM'''
from Bio import motifs
from Bio.Seq import Seq

# Create a motif from known instances
print('=== Create Motif from Instances ===')
instances = [Seq('TACAA'), Seq('TACGA'), Seq('TACTA'), Seq('TGCAA')]
m = motifs.create(instances)

# Consensus sequences
print(f'Consensus: {m.consensus}')
print(f'Degenerate consensus: {m.degenerate_consensus}')
print(f'Anticonsensus: {m.anticonsensus}')

# Position frequency matrix
print('\n=== Position Frequency Matrix ===')
print(m.counts)

# Position weight matrix (normalized frequencies)
print('\n=== Position Weight Matrix ===')
pwm = m.counts.normalize(pseudocounts=0.5)
print(pwm)

# Position-specific scoring matrix (log-odds)
print('\n=== PSSM (Log-Odds) ===')
pssm = pwm.log_odds()
print(pssm)

# Information content
print('\n=== Information Content ===')
print(f'Mean IC: {pssm.mean():.3f} bits')
print(f'Max possible score: {pssm.max:.3f}')
print(f'Min possible score: {pssm.min:.3f}')

# Search a sequence with the PSSM
print('\n=== PSSM Search ===')
test_seq = Seq('ATGCTACAAGCTACGATGCAACTACTA')
print(f'Searching: {test_seq}')
print(f'Threshold: 0.0')

for position, score in pssm.search(test_seq, threshold=0.0):
    subseq = test_seq[position:position + len(m.consensus)]
    print(f'Position {position}: {subseq} (score: {score:.2f})')

# Search both strands
print('\n=== Search Both Strands ===')
for position, score in pssm.search(test_seq, threshold=0.0, both=True):
    if position >= 0:
        strand = '+'
        subseq = test_seq[position:position + len(m.consensus)]
    else:
        strand = '-'
        pos = -position - len(m.consensus)
        subseq = test_seq.reverse_complement()[pos:pos + len(m.consensus)]
    print(f'{strand} strand position {abs(position)}: {subseq} (score: {score:.2f})')

# Calculate threshold from score distribution
print('\n=== Threshold Calculation ===')
sd = pssm.distribution()
print(f'Threshold for 1% FPR: {sd.threshold_fpr(0.01):.2f}')
print(f'Threshold for 10% FNR: {sd.threshold_fnr(0.1):.2f}')
