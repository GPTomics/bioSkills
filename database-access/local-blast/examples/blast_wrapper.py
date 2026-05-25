'''Python wrapper for BLAST+ via subprocess: builds v5 DB, runs with -task, parses tabular with qcovs/qcovhsp distinction.'''
# Reference: ncbi blast+ 2.15+ | Verify API if version differs
import subprocess
import shutil


def require(tool):
    if not shutil.which(tool):
        raise RuntimeError(f'{tool} not on PATH; conda install -c bioconda blast')


def make_blast_db(fasta, name, dbtype='nucl', parse_seqids=True, version=5):
    require('makeblastdb')
    cmd = ['makeblastdb', '-in', fasta, '-dbtype', dbtype, '-out', name,
           '-blastdb_version', str(version), '-hash_index',
           '-title', f'{fasta} ({dbtype})']
    if parse_seqids:
        cmd.append('-parse_seqids')
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'makeblastdb failed: {r.stderr}')


def run_blast(query, db, out, program='blastn', task=None, evalue=1e-10,
              threads=8, hitlist=500, soft_masking=True):
    '''hitlist=500 with post-filter avoids the max_target_seqs early-termination trap (Shah 2019).'''
    require(program)
    cmd = [program, '-query', query, '-db', db, '-out', out,
           '-evalue', str(evalue),
           '-num_threads', str(threads),
           '-max_target_seqs', str(hitlist),
           '-soft_masking', 'true' if soft_masking else 'false',
           '-outfmt', '6 qseqid sseqid pident length qcovs qcovhsp evalue bitscore staxids sscinames stitle']
    if task:
        cmd += ['-task', task]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f'{program} failed: {r.stderr}')


def parse_tabular(path):
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'qcovs', 'qcovhsp',
            'evalue', 'bitscore', 'staxids', 'sscinames', 'stitle']
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.rstrip('\n').split('\t')
            d = dict(zip(cols, vals + [''] * (len(cols) - len(vals))))
            for k in ('pident', 'qcovs', 'qcovhsp', 'evalue', 'bitscore'):
                d[k] = float(d[k]) if d[k] else 0.0
            d['length'] = int(d['length']) if d['length'] else 0
            rows.append(d)
    return rows


def top_by_bitscore_per_query(rows, n=1):
    by_q = {}
    for r in rows:
        by_q.setdefault(r['qseqid'], []).append(r)
    for q in by_q:
        by_q[q] = sorted(by_q[q], key=lambda x: -x['bitscore'])[:n]
    return by_q


def filter_hits(rows, min_pident=30.0, min_qcovs=50.0, max_evalue=1e-5):
    return [r for r in rows
            if r['pident'] >= min_pident
            and r['qcovs'] >= min_qcovs
            and r['evalue'] <= max_evalue]


if __name__ == '__main__':
    make_blast_db('reference.fasta', 'ref_db', dbtype='nucl')
    run_blast('query.fasta', 'ref_db', 'results.tsv', program='blastn',
              task='dc-megablast', threads=8)
    rows = parse_tabular('results.tsv')
    print(f'Total HSPs: {len(rows)}')
    print(f'Unique queries with any hit: {len(set(r["qseqid"] for r in rows))}')

    good = filter_hits(rows, min_pident=70.0, min_qcovs=80.0)
    print(f'After identity>=70 + qcovs>=80 + E<=1e-5: {len(good)}')

    print('\nTop hit per query by bit-score:')
    for q, top in list(top_by_bitscore_per_query(rows, n=1).items())[:5]:
        h = top[0]
        print(f'  {q} -> {h["sseqid"]}  bits={h["bitscore"]:.1f}  qcovs={h["qcovs"]:.0f}%  E={h["evalue"]:.1e}')
