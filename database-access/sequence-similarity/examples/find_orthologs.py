import subprocess
import argparse
import os

def run_blast(query, db, output, evalue=1e-10, threads=4):
    cmd = [
        'blastp', '-query', query, '-db', db,
        '-outfmt', '6', '-evalue', str(evalue),
        '-max_target_seqs', '1', '-num_threads', str(threads),
        '-out', output
    ]
    subprocess.run(cmd, check=True)

def parse_blast(filename):
    best_hits = {}
    with open(filename) as f:
        for line in f:
            parts = line.strip().split('\t')
            query, subject = parts[0], parts[1]
            if query not in best_hits:
                best_hits[query] = subject
    return best_hits

def find_rbh(forward, reverse):
    orthologs = []
    for a, b in forward.items():
        if b in reverse and reverse[b] == a:
            orthologs.append((a, b))
    return orthologs

def main():
    parser = argparse.ArgumentParser(description='Find orthologs via reciprocal best hits')
    parser.add_argument('species_a', help='FASTA file for species A')
    parser.add_argument('species_b', help='FASTA file for species B')
    parser.add_argument('-e', '--evalue', type=float, default=1e-10, help='E-value threshold')
    parser.add_argument('-t', '--threads', type=int, default=4, help='CPU threads')
    parser.add_argument('-o', '--output', default='orthologs.txt', help='Output file')

    args = parser.parse_args()

    print('Building BLAST databases...')
    subprocess.run(['makeblastdb', '-in', args.species_a, '-dbtype', 'prot', '-out', 'db_a'], check=True)
    subprocess.run(['makeblastdb', '-in', args.species_b, '-dbtype', 'prot', '-out', 'db_b'], check=True)

    print('Running forward BLAST...')
    run_blast(args.species_a, 'db_b', 'forward.txt', args.evalue, args.threads)

    print('Running reverse BLAST...')
    run_blast(args.species_b, 'db_a', 'reverse.txt', args.evalue, args.threads)

    print('Finding reciprocal best hits...')
    forward = parse_blast('forward.txt')
    reverse = parse_blast('reverse.txt')
    orthologs = find_rbh(forward, reverse)

    with open(args.output, 'w') as f:
        for a, b in orthologs:
            f.write(f'{a}\t{b}\n')

    print(f'Found {len(orthologs)} ortholog pairs')
    print(f'Output: {args.output}')

    for f in ['db_a.phr', 'db_a.pin', 'db_a.psq', 'db_b.phr', 'db_b.pin', 'db_b.psq', 'forward.txt', 'reverse.txt']:
        if os.path.exists(f):
            os.remove(f)

if __name__ == '__main__':
    main()
