'''Clean alignment by removing gappy columns and sequences'''

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'Original: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

def remove_gappy_columns(alignment, threshold=0.5):
    num_seqs = len(alignment)
    keep_columns = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        gap_fraction = column.count('-') / num_seqs
        if gap_fraction < threshold:
            keep_columns.append(col_idx)

    new_records = []
    for record in alignment:
        new_seq = ''.join(str(record.seq)[i] for i in keep_columns)
        new_records.append(SeqRecord(Seq(new_seq), id=record.id, description=record.description))
    return MultipleSeqAlignment(new_records)

def filter_by_gap_content(alignment, max_gap_fraction=0.2):
    filtered = []
    for record in alignment:
        gap_fraction = str(record.seq).count('-') / len(record.seq)
        if gap_fraction <= max_gap_fraction:
            filtered.append(record)
    return MultipleSeqAlignment(filtered)

# threshold=0.5: Remove columns with >=50% gaps. Standard cutoff for phylogenetics.
# Use 0.3 for stringent filtering; 0.7 if expecting many indels.
cleaned = remove_gappy_columns(alignment, threshold=0.5)
print(f'After column cleaning: {len(cleaned)} sequences, {cleaned.get_alignment_length()} columns')

# max_gap_fraction=0.2: Remove sequences with >20% gaps. Excludes fragmentary sequences.
# Use 0.1 for high-quality datasets; 0.3-0.4 for diverse or ancient sequences.
cleaned = filter_by_gap_content(cleaned, max_gap_fraction=0.2)
print(f'After sequence filtering: {len(cleaned)} sequences, {cleaned.get_alignment_length()} columns')

AlignIO.write(cleaned, 'cleaned_alignment.fasta', 'fasta')
print('Saved to cleaned_alignment.fasta')
