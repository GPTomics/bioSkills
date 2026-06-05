'''Prioritize pVACseq neoantigen candidates with the downstream filters that matter.

Binding is the entry gate, not the answer. This script encodes the load-bearing
downstream steps: drop candidates on HLA-LOH-lost alleles (silent invalidator),
rank by cancer cell fraction (CCF, not raw VAF), and weight agretopicity/foreignness
quality. The output is a tier-1 hypothesis list for MS + functional validation, not
a final answer. Input columns follow the pVACseq aggregate report.
'''
# Reference: pVACtools 4.1+, pandas 2.2+ | Verify API if version differs

import pandas as pd


def drop_lost_alleles(df, lost_alleles, allele_col='HLA Allele'):
    '''A peptide assigned to an HLA allele the tumor deleted (LOHHLA) is not weakly
    presented - it is not presented at all. Run LOHHLA first; this step errors
    silently if skipped.'''
    return df[~df[allele_col].isin(set(lost_alleles))].copy()


def add_quality(df, wt='Median WT IC50 Score', mt='Median MT IC50 Score'):
    '''Agretopicity / DAI = IC50_WT / IC50_MT (>1 = mutant binds better, surface not
    tolerized). Anchor-position mutations inflate DAI without changing the TCR-facing
    surface, so DAI is paired with, not substituted for, anchor evaluation.'''
    out = df.copy()
    out['agretopicity'] = out[wt] / out[mt].clip(lower=0.1)
    return out


def add_ccf(df, purity, vaf_col='Tumor DNA VAF', cn_col='Local Copy Number'):
    '''Clonality needs cancer cell fraction, not raw VAF: CCF = VAF * (CN_local*purity
    + 2*(1-purity)) / purity, capped at 1. Clonal (CCF ~1) targets beat subclonal.'''
    out = df.copy()
    cn = out[cn_col] if cn_col in out.columns else 2
    out['ccf'] = (out[vaf_col] * (cn * purity + 2 * (1 - purity)) / purity).clip(upper=1.0)
    out['clonal'] = out['ccf'] >= 0.8
    return out


def prioritize(df):
    '''Rank within a patient by presentation strength, clonality, expression, and
    quality. Composite scores are for ordering a candidate list, never for absolute
    go/no-go across patients.'''
    out = df.copy()
    out['binding_score'] = 1 - (out['Median MT IC50 Score'].clip(upper=5000) / 5000)
    out['agretopicity_score'] = out['agretopicity'].clip(upper=10) / 10
    out['expr_score'] = (out['Gene Expression'] / out['Gene Expression'].max()).fillna(0)
    out['priority'] = (0.30 * out['binding_score'] + 0.25 * out['ccf']
                       + 0.25 * out['agretopicity_score'] + 0.20 * out['expr_score'])
    return out.sort_values('priority', ascending=False)


if __name__ == '__main__':
    demo = pd.DataFrame({
        'Gene Name': ['TP53', 'KRAS', 'BRAF', 'PIK3CA', 'EGFR', 'NRAS'],
        'Mutation': ['R175H', 'G12D', 'V600E', 'E545K', 'L858R', 'Q61K'],
        'MT Epitope Seq': ['HMTEVVRHC', 'VVVGADGVGK', 'LATEKSRWSG', 'STRDPLSEIT', 'KITDFGLAKL', 'ILDTAGKEEY'],
        'HLA Allele': ['HLA-A*02:01', 'HLA-A*02:01', 'HLA-B*07:02', 'HLA-A*02:01', 'HLA-A*24:02', 'HLA-B*07:02'],
        'Median MT IC50 Score': [45, 120, 350, 85, 420, 180],
        'Median WT IC50 Score': [1500, 800, 350, 450, 500, 600],
        'Tumor DNA VAF': [0.45, 0.35, 0.25, 0.15, 0.40, 0.20],
        'Gene Expression': [25, 150, 80, 45, 200, 30],
    })

    lost = ['HLA-A*24:02']  # from LOHHLA: this allele was deleted in the tumor
    candidates = drop_lost_alleles(demo, lost)
    candidates = add_quality(candidates)
    candidates = add_ccf(candidates, purity=0.6)
    ranked = prioritize(candidates)
    ranked[['Gene Name', 'Mutation', 'HLA Allele', 'Median MT IC50 Score', 'agretopicity', 'ccf', 'priority']]
