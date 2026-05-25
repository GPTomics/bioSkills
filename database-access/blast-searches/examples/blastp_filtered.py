'''BLASTP with entrez_query pre-filter, composition-based statistics, and Swiss-Prot for reproducibility.'''
# Reference: biopython 1.83+, ncbi blast+ 2.15+ | Verify API if version differs
from Bio.Blast import NCBIWWW, NCBIXML

# Human hemoglobin alpha chain
HBA_HUMAN = ('MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH'
             'GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR')


def filter_top(record, top=10, min_cov=0.8):
    qlen = record.query_length
    out = []
    for aln in record.alignments:
        hsp = aln.hsps[0]
        cov = hsp.align_length / qlen
        if cov >= min_cov:
            out.append((aln, hsp))
    return sorted(out, key=lambda ah: -ah[1].bits)[:top]


print('=== BLASTP against Swiss-Prot, Mammalia pre-filter, CBS mode 2 ===')
handle = NCBIWWW.qblast(
    program='blastp',
    database='swissprot',
    sequence=HBA_HUMAN,
    entrez_query='Mammalia[Organism]',
    expect=1e-10,
    composition_based_statistics=2,  # Yu&Altschul 2005; default since BLAST+ 2.2.17
    matrix_name='BLOSUM62',
    hitlist_size=200,
    format_type='XML',
)
record = NCBIXML.read(handle); handle.close()

print(f'Total hits returned: {len(record.alignments)}')
print(f'\nTop 10 by bit-score, filtered to coverage >= 0.8:')
for aln, hsp in filter_top(record, top=10):
    pct = 100 * hsp.identities / hsp.align_length
    print(f'  {aln.accession:<14} bits={hsp.bits:>6.1f}  E={hsp.expect:.1e}  '
          f'{pct:>5.1f}% id  cov={hsp.align_length / record.query_length:.2f}')
    print(f'      {aln.title[:90]}')
