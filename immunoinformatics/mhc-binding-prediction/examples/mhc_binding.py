'''MHC class I binding/presentation prediction with MHCflurry.

Scores peptides against a patient genotype, classifies by %Rank (not raw nM,
which is allele-biased), and tiles a protein into candidate epitopes. The
Class1PresentationPredictor returns columns: peptide, sample_name, affinity (nM),
best_allele, processing_score, presentation_score; the affinity_percentile (%Rank)
column is only added when predict(..., include_affinity_percentile=True). Lower
affinity nM = stronger; lower affinity_percentile = stronger; higher
presentation_score = more likely naturally presented.
'''
# Reference: mhcflurry 2.1+, pandas 2.2+ | Verify API if version differs

import pandas as pd
from mhcflurry import Class1PresentationPredictor


def classify_by_percentile(affinity_percentile):
    '''Class I %Rank cutoffs (NetMHCpan convention). %Rank is comparable across
    alleles; raw nM is not. Strong <= 0.5%, weak <= 2.0%.'''
    if affinity_percentile <= 0.5:
        return 'strong'
    if affinity_percentile <= 2.0:
        return 'weak'
    return 'non-binder'


def predict_presentation(peptides, genotype, predictor=None):
    '''Batched presentation prediction for one patient genotype.

    genotype: list of HLA class I alleles, e.g. ['HLA-A*02:01', 'HLA-B*07:02'].
    Reports the best-presenting allele per peptide.
    '''
    predictor = predictor or Class1PresentationPredictor.load()
    df = predictor.predict(peptides=list(peptides), alleles={'patient': list(genotype)},
                           include_affinity_percentile=True, verbose=0)
    df['call'] = df['affinity_percentile'].apply(classify_by_percentile)
    return df


def scan_protein(protein_seq, genotype, lengths=(8, 9, 10, 11), percentile_cutoff=2.0):
    '''Tile a protein into 8-11mers (9mers dominate real ligands), score all
    windows in one batched call, return windows at/under the weak-binder cutoff.'''
    predictor = Class1PresentationPredictor.load()
    windows = [(protein_seq[i:i + k], i + 1, k) for k in lengths for i in range(len(protein_seq) - k + 1)]
    peptides = [w[0] for w in windows]
    pos = {w[0]: (w[1], w[2]) for w in windows}
    df = predictor.predict(peptides=peptides, alleles={'patient': list(genotype)},
                           include_affinity_percentile=True, verbose=0)
    df['position'] = df['peptide'].map(lambda p: pos[p][0])
    df['length'] = df['peptide'].map(lambda p: pos[p][1])
    df['call'] = df['affinity_percentile'].apply(classify_by_percentile)
    return df[df['affinity_percentile'] <= percentile_cutoff].sort_values('affinity_percentile')


def rank_neoantigen_candidates(df, expression_tpm):
    '''Guard against EL/MS abundance bias: presentation models over-rank peptides
    from abundant proteins, so a lowly expressed neoantigen can be real yet
    under-ranked. Join measured expression and judge within-target rather than
    trusting the presentation score against the proteome. expression_tpm maps
    peptide -> source-gene TPM.'''
    df = df.copy()
    df['expression_tpm'] = df['peptide'].map(expression_tpm).fillna(0.0)
    df['expressed'] = df['expression_tpm'] >= 1.0
    return df.sort_values(['expressed', 'affinity_percentile'], ascending=[False, True])


if __name__ == '__main__':
    known_epitopes = ['GILGFVFTL', 'NLVPMVATV', 'SIINFEKL']  # flu M1, CMV pp65, ovalbumin
    genotype = ['HLA-A*02:01', 'HLA-A*03:01', 'HLA-B*07:02']

    predictor = Class1PresentationPredictor.load()
    calls = predict_presentation(known_epitopes, genotype, predictor)
    calls[['peptide', 'best_allele', 'affinity', 'affinity_percentile', 'presentation_score', 'call']]

    test_protein = 'MSIINFEKLAAAGILGFVFTLVSSAYNLVPMVATVQTLNF'
    hits = scan_protein(test_protein, ['HLA-A*02:01'], lengths=(9,))
    hits[['peptide', 'position', 'best_allele', 'affinity_percentile', 'call']].head(10)
