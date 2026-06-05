'''Build correct HLA class II allele strings and apply class II %Rank cutoffs.

The decision-relevant nuance for class II is allele nomenclature and the DQ/DP
heterodimer pairing trap: a heterozygous donor expresses DR (single-chain, low
combinatorial burden) plus DQ/DP heterodimers, and only documented alpha/beta
pairings should be scored - never the full DQA1 x DQB1 cross product. This helper
formats alleles for NetMHCIIpan and classifies output by the looser class II
thresholds (strong <= 1%, weak <= 5%), distinct from class I (0.5%/2.0%).
'''
# Reference: NetMHCIIpan 4.3+, pandas 2.2+ | Verify API if version differs

import pandas as pd


def dr_allele(beta):
    '''DR is single-chain; the beta allele names the molecule. beta like "DRB1*01:01".'''
    return 'DRB1_' + beta.split('*')[1].replace(':', '')


def dq_dp_allele(alpha, beta):
    '''NetMHCIIpan heterodimer string, e.g. HLA-DQA10501-DQB10201.
    alpha/beta like "DQA1*05:01" / "DQB1*02:01".'''
    a = alpha.split('*')[0] + alpha.split('*')[1].replace(':', '')
    b = beta.split('*')[0] + beta.split('*')[1].replace(':', '')
    return f'HLA-{a}-{b}'


def class_ii_alleles(genotype):
    '''genotype: dict with keys DRB1 (list), and optional documented DQ/DP pairs as
    lists of (alpha, beta) tuples. Returning explicit pairs avoids inventing
    non-existent trans heterodimers.'''
    alleles = [dr_allele(b) for b in genotype.get('DRB1', [])]
    alleles += [dq_dp_allele(a, b) for a, b in genotype.get('DQ_pairs', [])]
    alleles += [dq_dp_allele(a, b) for a, b in genotype.get('DP_pairs', [])]
    return alleles


def classify_class_ii(rank):
    '''Class II %Rank cutoffs (NetMHCIIpan convention). LOWER = stronger.'''
    if rank <= 1.0:
        return 'strong'
    if rank <= 5.0:
        return 'weak'
    return 'non-binder'


def parse_netmhciipan_xls(path, rank_col='Rank'):
    '''Parse a NetMHCIIpan -xls table and flag DQ/DP calls as lower-confidence
    than DR (isotype confidence DR > DP > DQ).'''
    df = pd.read_csv(path, sep='\t')
    df['call'] = df[rank_col].apply(classify_class_ii)
    df['isotype'] = df['Allele'].str.contains('DQ').map({True: 'DQ', False: ''}).replace('', None)
    df['lower_confidence'] = df['Allele'].str.contains('DQ|DP')
    return df


if __name__ == '__main__':
    genotype = {
        'DRB1': ['DRB1*01:01', 'DRB1*15:01'],
        'DQ_pairs': [('DQA1*05:01', 'DQB1*02:01')],   # documented cis pair only
        'DP_pairs': [('DPA1*01:03', 'DPB1*04:01')],
    }
    class_ii_alleles(genotype)
    [classify_class_ii(r) for r in (0.4, 1.0, 3.0, 8.0)]
