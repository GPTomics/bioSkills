'''Permutation test for interval-overlap significance with pybedtools.'''
# Reference: pybedtools 0.10+, bedtools 2.31+ | Verify API if version differs
import argparse
import pybedtools

N_PERMUTATIONS = 1000   # >=1000 gives a stable empirical p down to ~0.001; fewer cannot resolve small p


def permutation_pvalue(query_bed, annotation_bed, genome_file, workspace_bed, blacklist_bed, n=N_PERMUTATIONS):
    '''Empirical p that query and annotation overlap more than a size-preserving, workspace-restricted null.'''
    query = pybedtools.BedTool(query_bed).sort()
    annotation = pybedtools.BedTool(annotation_bed).sort()
    observed = query.jaccard(annotation)['jaccard']
    null = [query.shuffle(g=genome_file, incl=workspace_bed, excl=blacklist_bed, chrom=True).sort().jaccard(annotation)['jaccard'] for _ in range(n)]
    hits = sum(x >= observed for x in null)
    p = (hits + 1) / (n + 1)   # +1 correction so a permutation p is never exactly 0 (Phipson & Smyth 2010)
    return {'observed': observed, 'null_mean': sum(null) / len(null), 'hits': hits, 'n': n, 'pvalue': p}


def main():
    ap = argparse.ArgumentParser(description='Permutation overlap-significance test (matched null).')
    ap.add_argument('--query', required=True)
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--genome', required=True, help='chrom<TAB>size file')
    ap.add_argument('--workspace', required=True, help='accessible regions = the universe (NOT the whole genome)')
    ap.add_argument('--blacklist', required=True, help='ENCODE blacklist + assembly gaps to exclude')
    ap.add_argument('--n', type=int, default=N_PERMUTATIONS)
    args = ap.parse_args()
    result = permutation_pvalue(args.query, args.annotation, args.genome, args.workspace, args.blacklist, args.n)
    print(result)


if __name__ == '__main__':
    main()
