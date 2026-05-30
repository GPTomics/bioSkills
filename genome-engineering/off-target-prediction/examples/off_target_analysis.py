'''Prepare a Cas-OFFinder search, parse its output, and CFD-rank candidate off-targets.

The CFD matrix is NOT hand-typed here: it is the published Doench 2016 table that ships
with CRISPOR / the Doench code as mismatch_score.pkl + pam_scores.pkl. cfd_score loads
those tables and takes the product of per-position mismatch penalties x the PAM penalty.
CFD is a SpCas9/NGG relative ranker for comparing guides, NOT a calibrated cutting
probability. The candidate list from Cas-OFFinder is a hypothesis set (predicted), not a
verdict -- validate real risk with an empirical assay and amplicon deep-seq.
'''
# Reference: cas-offinder 3.0+, pandas 2.2+ | Verify API if version differs

import pickle
import pandas as pd

CFD_MISMATCH_PKL = 'mismatch_score.pkl'   # ships with CRISPOR / Doench 2016 code; do not fabricate
CFD_PAM_PKL = 'pam_scores.pkl'


def prepare_cas_offinder_input(guides, genome_dir, pam='NRG', guide_len=20,
                               max_mismatches=4, dna_bulge=0, rna_bulge=0):
    '''Build a Cas-OFFinder input file string.

    pam 'NRG' also catches NAG/NGG off-targets. The query line is the guide bases plus N's
    for the PAM positions, so its length equals the pattern length. A bulge line is emitted
    only when a bulge size is requested (native support requires Cas-OFFinder >= 3.0.0).
    max_mismatches: 4 balances speed/sensitivity; bulges/variants can rescue more-distant sites.
    '''
    pattern = 'N' * guide_len + pam
    lines = [genome_dir]
    if dna_bulge or rna_bulge:
        lines.append(f'{dna_bulge} {rna_bulge}')
    lines.append(pattern)
    pam_pad = 'N' * len(pam)
    for g in guides:
        lines.append(f'{g.upper()}{pam_pad} {max_mismatches}')
    return '\n'.join(lines) + '\n'


def parse_cas_offinder_output(filepath):
    '''Parse Cas-OFFinder output (no header line; positions are 0-based, Bowtie convention).

    Column count depends on mode. No-bulge output has 6: pattern, chrom, position, site,
    strand, mismatches. Bulge mode (Cas-OFFinder >= 3.0.0) has 8: bulge_type, pattern, site,
    chrom, position, strand, mismatches, bulge_size. Dispatch on the observed width.
    '''
    df = pd.read_csv(filepath, sep='\t', header=None, comment='#')
    if df.shape[1] == 8:
        df.columns = ['bulge_type', 'pattern', 'site', 'chrom', 'position', 'strand', 'mismatches', 'bulge_size']
    else:
        df.columns = ['pattern', 'chrom', 'position', 'site', 'strand', 'mismatches']
    return df


def load_cfd_tables(mismatch_pkl=CFD_MISMATCH_PKL, pam_pkl=CFD_PAM_PKL):
    with open(mismatch_pkl, 'rb') as f:
        mismatch = pickle.load(f)
    with open(pam_pkl, 'rb') as f:
        pam = pickle.load(f)
    return mismatch, pam


def cfd_score(guide, off_protospacer, off_pam, mismatch_table, pam_table):
    '''CFD = product of per-position mismatch penalties x PAM penalty (Doench 2016 tables).

    Keys follow the published convention: mismatch 'rX:dY,pos' (RNA base X vs DNA-complement
    base Y at 1-based position), PAM the last two PAM bases. A matching position contributes 1.
    '''
    guide, off = guide.upper().replace('T', 'U'), off_protospacer.upper()
    score = 1.0
    for i, (g, o) in enumerate(zip(guide, off.replace('T', 'U')), 1):
        if g != o:
            comp = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A'}[o]   # DNA base complementary to the off-target base
            score *= mismatch_table.get(f'r{g}:d{comp},{i}', 1.0)
    return score * pam_table.get(off_pam[-2:].upper(), 1.0)


def aggregate_specificity(cfd_scores):
    '''CRISPOR CFD specificity on a 0-100 scale: 100/(1 + sum(off-target CFDs)), CFDs on 0-1.

    Equivalently 10000/(100 + 100*sum) -- the CRISPOR/Hsu aggregate folded for 0-1 inputs.
    Higher = more specific (100 = no off-targets). Use to compare candidate GUIDES; >~80 is a
    research heuristic, NOT a clinical safety pass.
    '''
    return 100.0 / (1.0 + sum(cfd_scores))


if __name__ == '__main__':
    guides = ['GGCCGACCTGTCGCTGACGC', 'CGCCAGCGTCAGCGACAGGT']
    print('Cas-OFFinder input (NRG PAM, 2 nt DNA + RNA bulges, <=4 mismatches):')
    print(prepare_cas_offinder_input(guides, '/path/to/genome_dir',
                                     max_mismatches=4, dna_bulge=2, rna_bulge=2))

    # In practice: off_targets = parse_cas_offinder_output('output.txt'); load tables; score each site.
    # Here, illustrate the aggregation math with per-site CFDs as the published table would return
    # them (these are example values, not a fabricated matrix):
    example_site_cfds = [0.71, 0.18, 0.04, 0.02]
    spec = aggregate_specificity(example_site_cfds)
    print(f'Aggregate specificity for the example candidate set: {spec:.1f} '
          f'(compare among guides; >~80 = good for guide selection, not a safety verdict)')
    print('Highest-CFD candidate (0.71) is a design-time concern -> nominate for empirical '
          'validation (predicted -> detected -> validated); report the LoD with any "no off-target" claim.')
