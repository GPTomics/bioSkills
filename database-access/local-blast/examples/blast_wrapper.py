'''Python wrapper for local BLAST operations'''
import subprocess
import os
from pathlib import Path

def make_blast_db(fasta_file, db_name, db_type='nucl', parse_seqids=True):
    '''Create a BLAST database from a FASTA file'''
    cmd = ['makeblastdb', '-in', fasta_file, '-dbtype', db_type, '-out', db_name]
    if parse_seqids:
        cmd.append('-parse_seqids')

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {result.stderr}")
    print(f"Created {db_type} database: {db_name}")

def run_blast(query, db, output, program='blastn', evalue=1e-5, threads=4,
              max_targets=10, outfmt='6'):
    '''Run a BLAST search'''
    cmd = [
        program,
        '-query', query,
        '-db', db,
        '-out', output,
        '-outfmt', str(outfmt),
        '-evalue', str(evalue),
        '-num_threads', str(threads),
        '-max_target_seqs', str(max_targets)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{program} failed: {result.stderr}")
    print(f"Results saved to: {output}")

def parse_blast_tabular(filename):
    '''Parse BLAST tabular output (outfmt 6)'''
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    hits = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            values = line.strip().split('\t')
            if len(values) >= 12:
                hit = dict(zip(columns, values[:12]))
                hit['pident'] = float(hit['pident'])
                hit['evalue'] = float(hit['evalue'])
                hit['length'] = int(hit['length'])
                hit['bitscore'] = float(hit['bitscore'])
                hits.append(hit)
    return hits

def get_best_hits(hits):
    '''Get best hit per query (lowest E-value)'''
    best = {}
    for hit in hits:
        qid = hit['qseqid']
        if qid not in best or hit['evalue'] < best[qid]['evalue']:
            best[qid] = hit
    return list(best.values())

def filter_hits(hits, min_identity=90, min_length=100, max_evalue=1e-10):
    '''Filter hits by criteria'''
    return [h for h in hits
            if h['pident'] >= min_identity
            and h['length'] >= min_length
            and h['evalue'] <= max_evalue]

# Example usage
if __name__ == '__main__':
    # Create database
    make_blast_db('reference.fasta', 'ref_db', 'nucl')

    # Run BLAST
    run_blast('query.fasta', 'ref_db', 'results.tsv', threads=8)

    # Parse and analyze
    hits = parse_blast_tabular('results.tsv')
    print(f"\nTotal hits: {len(hits)}")

    best = get_best_hits(hits)
    print(f"Unique queries with hits: {len(best)}")

    good_hits = filter_hits(hits, min_identity=95)
    print(f"Hits with >=95% identity: {len(good_hits)}")

    print("\nTop 5 hits:")
    for hit in sorted(hits, key=lambda x: x['evalue'])[:5]:
        print(f"  {hit['qseqid']} -> {hit['sseqid']}: {hit['pident']:.1f}%, E={hit['evalue']:.2e}")
