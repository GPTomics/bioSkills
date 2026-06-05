'''Rank neoantigen candidates within a patient with feature transparency.

Immunogenicity is the least-solved layer: dedicated scores reach only ~AUROC 0.6-0.7
and correlated poorly with validated immunogenicity in TESLA, where presentation
strength, abundance, agretopicity, and foreignness carried the signal. So this script
does NOT emit a single composite verdict. It filters on the non-negotiable
expression/clonality gates, computes agretopicity defensively (flagging anchor
inflation and WT-denominator instability), and orders candidates within one patient
while keeping every feature visible. Scores are within-patient only; never compare
across patients or alleles.
'''
# Reference: pandas 2.2+ | Verify API if version differs

import pandas as pd


def filter_expression_clonality(df, tpm_min=1.0, vaf_min=0.25):
    '''Filters, not scores: an unexpressed or low-VAF peptide is never displayed,
    so it is gated out before ranking matters.'''
    return df[(df['gene_tpm'] >= tpm_min) & (df['rna_vaf'] >= vaf_min)].copy()


def defensive_dai(df, wt='wt_ic50', mt='mt_ic50', anchor='mutation_at_anchor', wt_cap=5000):
    '''DAI = IC50_WT / IC50_MT, with both traps flagged: anchor mutations inflate DAI
    without changing the TCR-facing surface, and a barely-presented WT makes the ratio
    numerical noise.'''
    out = df.copy()
    out['dai'] = out[wt] / out[mt]
    out['dai_anchor_artifact'] = out[anchor]
    out['dai_unstable'] = out[wt] > wt_cap
    out['dai_trustworthy'] = ~out['dai_anchor_artifact'] & ~out['dai_unstable']
    return out


def rank_within_patient(df):
    '''Order by presentation, abundance, then quality - the axes TESLA found carry the
    signal - keeping features side by side for human curation rather than collapsing
    them into one over-trusted number.'''
    return df.sort_values(['presentation_rank', 'gene_tpm', 'foreignness'],
                          ascending=[True, False, False])


if __name__ == '__main__':
    demo = pd.DataFrame({
        'gene': ['TP53', 'KRAS', 'BRAF', 'EGFR', 'PIK3CA'],
        'mt_peptide': ['HMTEVVRHC', 'VVVGADGVGK', 'LATEKSRWSG', 'KITDFGLAKL', 'STRDPLSEIT'],
        'presentation_rank': [0.3, 0.8, 1.5, 0.6, 2.2],
        'mt_ic50': [45, 120, 350, 90, 600],
        'wt_ic50': [1500, 800, 360, 8000, 700],
        'mutation_at_anchor': [False, False, True, False, False],
        'foreignness': [0.7, 0.4, 0.2, 0.5, 0.3],
        'gene_tpm': [25, 150, 80, 200, 12],
        'rna_vaf': [0.45, 0.35, 0.30, 0.40, 0.15],
    })

    kept = filter_expression_clonality(demo)
    scored = defensive_dai(kept)
    ranked = rank_within_patient(scored)
    ranked[['gene', 'presentation_rank', 'gene_tpm', 'dai', 'dai_trustworthy', 'foreignness']]
