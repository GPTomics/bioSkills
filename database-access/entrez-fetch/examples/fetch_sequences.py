'''Fetch sequences with version-pinned accessions, gbwithparts for WGS, and fasta_cds_aa for one-shot CDS translation.'''
# Reference: biopython 1.83+, entrez direct 21.0+ | Verify API if version differs
from Bio import Entrez, SeqIO
import time

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def sniff_then_parse(text, expected_start, fmt):
    '''Guard against HTML error pages returned with HTTP 200.'''
    if not text.lstrip().startswith(expected_start):
        raise RuntimeError(f'Unexpected response: starts {text[:80]!r}')
    from io import StringIO
    return SeqIO.read(StringIO(text), fmt)


print('=== Version-pinned GenBank for reproducibility ===')
h = Entrez.efetch(db='nucleotide', id='NM_007294.4', rettype='gb', retmode='text')
gb = sniff_then_parse(h.read(), 'LOCUS', 'genbank'); h.close()
print(f'Accession: {gb.id}')
print(f'Definition: {gb.description}')
print(f'Length: {len(gb.seq)} nt')
print(f'CDS features: {sum(1 for f in gb.features if f.type == "CDS")}')
time.sleep(DELAY)

print('\n=== fasta_cds_aa: one-shot CDS translation ===')
h = Entrez.efetch(db='nucleotide', id='NC_000913.3', rettype='fasta_cds_aa', retmode='text')
proteins = list(SeqIO.parse(h, 'fasta')); h.close()
print(f'NC_000913.3 (E. coli K-12) CDS-translated proteins: {len(proteins)}')
print(f'First: {proteins[0].id}  length={len(proteins[0].seq)}')
time.sleep(DELAY)

print('\n=== Multi-accession batch with mixed-type guard ===')
ids = ['NM_007294.4', 'NM_000059.4', 'NM_000546.6']
assert all(id.count('.') == 1 for id in ids), 'Mix of versioned/unversioned IDs in batch'
h = Entrez.efetch(db='nucleotide', id=','.join(ids), rettype='fasta', retmode='text')
records = list(SeqIO.parse(h, 'fasta')); h.close()
for r in records:
    print(f'  {r.id}: {len(r.seq)} nt')
assert len(records) == len(ids), f'Truncated batch: requested {len(ids)} got {len(records)}'
time.sleep(DELAY)

print('\n=== WGS gotcha: gb vs gbwithparts ===')
wgs_accession = 'NZ_CP063325.1'  # a complete genome accession; gbwithparts ensures sequence is inlined
h = Entrez.efetch(db='nucleotide', id=wgs_accession, rettype='gbwithparts', retmode='text')
text = h.read(); h.close()
print(f'gbwithparts byte length: {len(text):,} (expect large -- sequence inlined)')
time.sleep(DELAY)
