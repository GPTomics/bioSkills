#!/usr/bin/env python3
'''Strain identification using sourmash'''

import sourmash
from pathlib import Path
import sys

def compute_signature(fasta_path, ksize=31, scaled=1000):
    '''Compute sourmash signature for a genome'''
    mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
    for record in sourmash.load_file_as_signatures(fasta_path):
        return record
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                mh.add_sequence(line.strip().upper())
    return sourmash.SourmashSignature(mh, name=Path(fasta_path).stem)

def compare_genomes(fasta_files):
    '''Compare multiple genomes using sourmash'''
    sigs = []
    for fasta in fasta_files:
        print(f'Computing signature for {fasta}...')
        sig = sourmash.load_one_signature(str(fasta), ksize=31)
        sigs.append(sig)

    print('\n=== Pairwise Similarities ===')
    print('Genome1\tGenome2\tContainment\tJaccard')
    for i, sig1 in enumerate(sigs):
        for j, sig2 in enumerate(sigs):
            if i < j:
                containment = sig1.contained_by(sig2)
                jaccard = sig1.jaccard(sig2)
                print(f'{sig1.name}\t{sig2.name}\t{containment:.4f}\t{jaccard:.4f}')

def main():
    input_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('.')
    fasta_files = list(input_dir.glob('*.sig'))

    if not fasta_files:
        print('No .sig files found. Creating signatures from FASTA files...')
        fasta_files = list(input_dir.glob('*.fasta')) + list(input_dir.glob('*.fa'))
        for fasta in fasta_files:
            cmd = f'sourmash sketch dna -p k=31,scaled=1000 {fasta}'
            import subprocess
            subprocess.run(cmd.split(), check=True)
        fasta_files = list(input_dir.glob('*.sig'))

    if fasta_files:
        compare_genomes(fasta_files)
    else:
        print('No FASTA or signature files found')

if __name__ == '__main__':
    main()
