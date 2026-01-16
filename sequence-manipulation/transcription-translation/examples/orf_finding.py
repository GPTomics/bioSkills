'''Finding open reading frames in all six frames'''
from Bio.Seq import Seq

def six_frame_translation(seq):
    '''Translate a sequence in all six reading frames'''
    frames = []
    for strand, s in [('+', seq), ('-', seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(s) - frame) // 3)
            fragment = s[frame:frame + length]
            frames.append((strand, frame, fragment.translate()))
    return frames

def find_orfs(seq, min_protein_length=10):
    '''Find all ORFs (start codon to stop codon) in all six frames'''
    orfs = []
    for strand, s in [('+', seq), ('-', seq.reverse_complement())]:
        for frame in range(3):
            trans = str(s[frame:].translate())
            aa_start = 0
            while True:
                start = trans.find('M', aa_start)
                if start == -1:
                    break
                stop = trans.find('*', start)
                if stop == -1:
                    stop = len(trans)
                orf = trans[start:stop]
                if len(orf) >= min_protein_length:
                    nt_start = start * 3 + frame
                    nt_end = stop * 3 + frame
                    orfs.append({'strand': strand, 'frame': frame, 'nt_start': nt_start, 'nt_end': nt_end, 'protein': orf})
                aa_start = start + 1
    return orfs

# Example sequence
seq = Seq('ATGCGATCGATCGATGTTTTGGCATTAAGATCGATCGATCGATCG')

print('=== Six Frame Translation ===')
for strand, frame, protein in six_frame_translation(seq):
    print(f'{strand}{frame}: {protein}')

print('\n=== ORFs (min 3 aa) ===')
orfs = find_orfs(seq, min_protein_length=3)
for orf in orfs:
    print(f"{orf['strand']} frame {orf['frame']} [{orf['nt_start']}:{orf['nt_end']}]: {orf['protein']}")
