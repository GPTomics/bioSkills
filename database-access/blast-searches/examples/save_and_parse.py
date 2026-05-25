'''Persist BLAST XML for later re-parsing; demonstrates bit-score vs E-value sorting and identity/coverage filtering.'''
# Reference: biopython 1.83+, ncbi blast+ 2.15+ | Verify API if version differs
from Bio.Blast import NCBIWWW, NCBIXML

QUERY = '''>HBB_partial
ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG'''


def run_and_save(sequence, out_xml, program='blastn', database='refseq_select_rna', hitlist=500):
    print(f'Running {program} against {database} (hitlist_size={hitlist}; avoids max_target_seqs trap)...')
    handle = NCBIWWW.qblast(program, database, sequence,
                            hitlist_size=hitlist, format_type='XML',
                            expect=1e-10)
    with open(out_xml, 'w') as f:
        f.write(handle.read())
    handle.close()
    print(f'XML saved to {out_xml}')


def parse_hits(xml_path):
    with open(xml_path) as f:
        record = NCBIXML.read(f)
    qlen = record.query_length
    hits = []
    for aln in record.alignments:
        hsp = aln.hsps[0]
        hits.append({
            'accession': aln.accession,
            'title': aln.title,
            'evalue': hsp.expect,
            'bits': hsp.bits,
            'identity': hsp.identities / hsp.align_length,
            'coverage': hsp.align_length / qlen,
            'q_range': (hsp.query_start, hsp.query_end),
            's_range': (hsp.sbjct_start, hsp.sbjct_end),
        })
    return hits, record


run_and_save(QUERY, 'blast_results.xml')

hits, record = parse_hits('blast_results.xml')
print(f'\nParsed {len(hits)} alignments from saved XML')
print(f'Query length: {record.query_length}')

print('\nTop 10 by bit-score (correct for cross-DB comparison):')
for h in sorted(hits, key=lambda x: -x['bits'])[:10]:
    print(f'  {h["accession"]:<14}  bits={h["bits"]:>6.1f}  E={h["evalue"]:.1e}  '
          f'id={h["identity"]:.2f}  cov={h["coverage"]:.2f}')

print('\nTop 10 by E-value (same query+DB -- ranking should agree with bit-score):')
for h in sorted(hits, key=lambda x: x['evalue'])[:10]:
    print(f'  {h["accession"]:<14}  E={h["evalue"]:.1e}  bits={h["bits"]:>6.1f}')
