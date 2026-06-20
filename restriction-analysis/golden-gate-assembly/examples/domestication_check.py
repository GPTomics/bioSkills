'''Scan parts for internal Type IIS sites (Golden Gate domestication)'''
# Reference: biopython 1.83+ (API verified on 1.86) | Verify API if version differs
# A part with an internal copy of the assembly enzyme's site is fragmented during the
# one-pot reaction. enzyme.search() finds sites on BOTH strands (the site is asymmetric).

from Bio import SeqIO
from Bio.Restriction import BsaI, BsmBI, BbsI, SapI

ENZYMES = [
    (BsaI, 'Golden Gate / MoClo Level 1'),
    (BsmBI, 'MoClo Level 0 / Level 2'),
    (BbsI, 'Golden Gate alternative'),
    (SapI, '3-nt overhang, in-frame fusions'),
]

records = list(SeqIO.parse('parts.fasta', 'fasta'))
print(f'Screening {len(records)} part(s) for internal Type IIS sites')
print('=' * 60)

for record in records:
    print(f'\n{record.id} ({len(record.seq)} bp):')
    for enzyme, role in ENZYMES:
        hits = enzyme.search(record.seq)
        if hits:
            print(f'  {enzyme} ({enzyme.site}) [{role}]: {len(hits)} internal site(s) at {hits} -> DOMESTICATE')
        else:
            print(f'  {enzyme} ({enzyme.site}) [{role}]: clean')
