'''Reproducible BLASTN against refseq_select with large hitlist_size + post-filter, demonstrating bit-score sorting.'''
# Reference: biopython 1.83+, ncbi blast+ 2.15+ | Verify API if version differs
from Bio.Blast import NCBIWWW, NCBIXML

# Human beta-globin partial CDS
QUERY = '''>HBB_partial
ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG'''


def top_n_by_bitscore(record, n=10, min_identity=0.7, min_coverage=0.5):
    qlen = record.query_length
    hits = []
    for aln in record.alignments:
        hsp = aln.hsps[0]
        ident = hsp.identities / hsp.align_length
        cov = hsp.align_length / qlen
        if ident >= min_identity and cov >= min_coverage:
            hits.append({
                'accession': aln.accession,
                'title': aln.title,
                'evalue': hsp.expect,
                'bits': hsp.bits,
                'identity': ident,
                'coverage': cov,
            })
    return sorted(hits, key=lambda h: -h['bits'])[:n]


print('Submitting BLASTN to NCBI (30-60s typical)...')
handle = NCBIWWW.qblast(
    program='blastn',
    database='refseq_select_rna',  # reproducible; nt would not be
    sequence=QUERY,
    expect=1e-10,
    word_size=11,                  # blastn default; use 28 only for >=95% identity targets
    hitlist_size=500,              # large; avoids max_target_seqs trap (Shah 2019)
    format_type='XML',
)
record = NCBIXML.read(handle); handle.close()

print(f'\nQuery: {record.query[:60]}')
print(f'Query length: {record.query_length}')
print(f'Database: {record.database}')
print(f'Total alignments returned: {len(record.alignments)}')

print('\nTop 10 by bit-score (database-size invariant):')
for i, hit in enumerate(top_n_by_bitscore(record, n=10), 1):
    print(f'  {i:>2}. {hit["accession"]:<14}  bits={hit["bits"]:>6.1f}  E={hit["evalue"]:.1e}  '
          f'id={hit["identity"]:.2f}  cov={hit["coverage"]:.2f}')
    print(f'      {hit["title"][:90]}')
