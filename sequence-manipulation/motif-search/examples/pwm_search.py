'''Overlapping motif counting and PWM/PSSM binding-site scoring on both strands'''
# Reference: biopython 1.83+ | Verify API if version differs
import re
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

print('=== Overlapping-Match Trap ===')
target = 'AAGCGCGCGAA'
print(f'target: {target}, motif: GCGC')
print(f"str.count (wrong, drops overlaps): {target.count('GCGC')}")
print(f"lookahead starts (correct): {[m.start(1) for m in re.finditer(r'(?=(GCGC))', target)]}")

result = nt_search(target, 'GCGC')
positions = result[1:] if len(result) > 1 else []
print(f'nt_search (expands IUPAC, reports overlaps): {result} -> starts {positions}')

degenerate = nt_search('AAGATCCTC', 'GATNNTC')
print(f'nt_search degenerate GATNNTC: {degenerate}')

print('\n=== Build PWM / PSSM with Pseudocounts ===')
instances = [Seq('TACAA'), Seq('TACGA'), Seq('TACTA'), Seq('TGCAA')]
m = motifs.create(instances)
m.pseudocounts = 0.5  # avoids -inf log-odds from count-0 cells; set BEFORE reading m.pssm
m.background = None    # uniform 0.25 background
pssm = m.pssm          # log2-odds, recomputed from pseudocounts/background on access

print(f'consensus: {m.consensus}')
print(f'degenerate consensus: {m.degenerate_consensus}')
print(f'mean information content: {pssm.mean():.3f} bits')

print('\n=== Reproducible Threshold from Score Distribution ===')
dist = pssm.distribution(background=m.background, precision=10**4)
fpr_threshold = dist.threshold_fpr(0.01)  # 1% false-positive rate, not a hand-picked cutoff
print(f'threshold at 1% FPR: {fpr_threshold:.3f} bits')

print('\n=== Scan Both Strands ===')
seq = Seq('ATGCTACAAGCTTGTAGCTACGA')
rc_pssm = pssm.reverse_complement()
hits = []
for position, score in pssm.search(seq, threshold=fpr_threshold, both=False):
    hits.append(('+', position, seq[position:position + len(m)], score))
for position, score in rc_pssm.search(seq, threshold=fpr_threshold, both=False):
    window = seq[position:position + len(m)]
    hits.append(('-', position, window.reverse_complement(), score))

for strand, position, site, score in sorted(hits, key=lambda h: h[1]):
    print(f'{strand} strand pos {position}: {site} (score {score:.2f})')
