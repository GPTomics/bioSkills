'''Read A2M / A3M alignments (HMMER, HHsuite, ColabFold) and extract match-only columns.

A2M: each sequence has identical match-column count; lowercase = insert state.
A3M: insert columns are not padded; lowercase characters appear inline per sequence.
'''
# Reference: biopython 1.83+ | Verify API if version differs

from Bio import AlignIO

def match_only_columns(a2m_alignment):
    return [
        ''.join(c for c in str(record.seq) if c.isupper() or c == '-')
        for record in a2m_alignment
    ]

if __name__ == '__main__':
    alignment = AlignIO.read('hhsearch_output.a2m', 'fasta')
    print(f'A2M alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

    match_columns = match_only_columns(alignment)
    print(f'After dropping insert states: {len(match_columns[0])} match columns')
    for record, match_only in list(zip(alignment, match_columns))[:5]:
        inserts = sum(1 for c in str(record.seq) if c.islower())
        print(f'{record.id}: {len(match_only)} match, {inserts} inserts')
