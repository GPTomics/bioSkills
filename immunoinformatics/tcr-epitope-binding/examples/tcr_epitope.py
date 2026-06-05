'''TCR specificity by the honest route: cluster + database lookup, not de-novo prediction.

General TCR-epitope prediction for unseen epitopes collapses to near-random because
labeled data is dominated by a few immunodominant epitopes and there is no true
negative set. The defensible task is unsupervised clustering ("these TCRs likely
share a specificity") plus database lookup, propagating labels by guilt-by-association
to clusters containing a known member. This script clusters CDR3b by edit distance,
annotates by exact VDJdb match (confidence-filtered), and keeps HLA as a covariate.
For truly de-novo epitopes, rank-and-validate; never report a confident per-pair call.
'''
# Reference: pandas 2.2+, scipy 1.12+ | Verify API if version differs

import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def levenshtein(s1, s2):
    if len(s1) < len(s2):
        s1, s2 = s2, s1
    prev = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        cur = [i + 1]
        for j, c2 in enumerate(s2):
            cur.append(min(prev[j + 1] + 1, cur[j] + 1, prev[j] + (c1 != c2)))
        prev = cur
    return prev[-1]


def cluster_cdr3(cdr3_list, max_edits=2):
    '''Edit-distance clustering as a dependency-free stand-in for tcrdist3/GLIPH2.
    Cluster within one HLA-coherent cohort; the same CDR3 on a different allele is a
    different specificity.'''
    n = len(cdr3_list)
    dist = [[levenshtein(cdr3_list[i], cdr3_list[j]) for j in range(n)] for i in range(n)]
    z = linkage(squareform(dist, checks=False), method='average')
    return dict(zip(cdr3_list, fcluster(z, t=max_edits, criterion='distance')))


def lookup_vdjdb(cdr3_list, vdjdb, min_confidence=1):
    '''Exact CDR3b match against confidence-filtered VDJdb; report the hit and HLA, not
    a binding probability. Near-matches are handled by clustering, not claimed as binders.'''
    db = vdjdb[vdjdb['vdjdb.score'] >= min_confidence]
    hits = db[db['cdr3'].isin(set(cdr3_list))]
    return hits[['cdr3', 'antigen.epitope', 'antigen.species', 'mhc.a']]


def annotate_clusters(clusters, hits):
    '''Propagate a known epitope label to every member of a cluster containing a hit.'''
    epi_by_cdr3 = dict(zip(hits['cdr3'], hits['antigen.epitope']))
    cluster_epi = {}
    for cdr3, cl in clusters.items():
        if cdr3 in epi_by_cdr3:
            cluster_epi[cl] = epi_by_cdr3[cdr3]
    return {cdr3: cluster_epi.get(cl, 'unknown') for cdr3, cl in clusters.items()}


if __name__ == '__main__':
    repertoire = ['CASSIRSSYEQYF', 'CASSIRSAYEQYF', 'CASSLAPGATNEKLFF', 'CASSPGTGGYEQYF']
    vdjdb = pd.DataFrame({
        'cdr3': ['CASSIRSSYEQYF'],
        'antigen.epitope': ['GILGFVFTL'],
        'antigen.species': ['InfluenzaA'],
        'mhc.a': ['HLA-A*02:01'],
        'vdjdb.score': [2],
    })
    clusters = cluster_cdr3(repertoire)
    hits = lookup_vdjdb(repertoire, vdjdb)
    annotate_clusters(clusters, hits)
