'''Find enzyme pairs for directional cloning'''

from Bio import SeqIO
from Bio.Restriction import Analysis, CommOnly

vector = SeqIO.read('vector.fasta', 'fasta')
insert = SeqIO.read('insert.fasta', 'fasta')

print(f'Vector: {vector.id} ({len(vector.seq)} bp)')
print(f'Insert: {insert.id} ({len(insert.seq)} bp)')
print('=' * 60)

vec_analysis = Analysis(CommOnly, vector.seq, linear=False)
ins_analysis = Analysis(CommOnly, insert.seq)

vec_once = set(vec_analysis.once_cutters().keys())
ins_non = set(ins_analysis.only_dont_cut())

candidates = vec_once & ins_non

print(f'\nEnzymes that cut vector once AND do not cut insert: {len(candidates)}')

five_prime = [e for e in candidates if e.is_5overhang()]
three_prime = [e for e in candidates if e.is_3overhang()]
blunt = [e for e in candidates if e.is_blunt()]

print(f"\n5' overhang enzymes ({len(five_prime)}):")
for e in sorted(five_prime, key=lambda x: str(x))[:10]:
    sites = vec_analysis.full()[e]
    print(f'  {e}: site={e.site}, ovhg={e.ovhgseq}, vec_pos={sites[0]}')

print(f"\n3' overhang enzymes ({len(three_prime)}):")
for e in sorted(three_prime, key=lambda x: str(x))[:10]:
    sites = vec_analysis.full()[e]
    print(f'  {e}: site={e.site}, ovhg={e.ovhgseq}, vec_pos={sites[0]}')

print(f'\nBlunt enzymes ({len(blunt)}):')
for e in sorted(blunt, key=lambda x: str(x))[:10]:
    sites = vec_analysis.full()[e]
    print(f'  {e}: site={e.site}, vec_pos={sites[0]}')

if five_prime and three_prime:
    e1, e2 = five_prime[0], three_prime[0]
    print(f"\n\nSuggested pair for directional cloning:")
    print(f"  5' end: {e1} ({e1.site})")
    print(f"  3' end: {e2} ({e2.site})")
