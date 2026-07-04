# Reference: tcrdist3 0.2+ olga 1.2+ | Verify API if version differs
'''Responsible TCR-beta specificity annotation, clustering, and a generation-probability null.

Runs end to end with a tiny synthetic repertoire and a synthetic VDJdb-style table so it needs
no downloads. Database annotation is pure pandas and always runs. tcrdist3 clustering and OLGA
Pgen are import-guarded: if the package is absent the step is skipped with a message, so the
core lesson (a match is a hypothesis, tested against V/HLA concordance and a Pgen null) still runs.
'''
import numpy as np
import pandas as pd


def synthetic_repertoire():
    return pd.DataFrame({
        'cdr3_b_aa': ['CASSIRSSYEQYF', 'CASSLGQAYEQYF', 'CASSPGTGGYEQYF', 'CASSFGTEAFF', 'CASSWDRGNTGELFF'],
        'v_b_gene': ['TRBV19*01', 'TRBV7-9*01', 'TRBV19*01', 'TRBV28*01', 'TRBV6-5*01'],
        'j_b_gene': ['TRBJ2-7*01', 'TRBJ2-7*01', 'TRBJ2-7*01', 'TRBJ1-1*01', 'TRBJ2-2*01'],
        'count': [812, 47, 33, 5, 19]})


def synthetic_vdjdb():
    # A curated-style export: CDR3 + V + restricting HLA + epitope + confidence score 0-3.
    return pd.DataFrame({
        'cdr3_b_aa': ['CASSIRSSYEQYF', 'CASSIRSSYEQYF', 'CASSPGTGGYEQYF', 'CASSKTGGSNEQFF'],
        'v_b_gene': ['TRBV19*01', 'TRBV5-1*01', 'TRBV19*01', 'TRBV4-1*01'],
        'mhc_a': ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-A*02:01', 'HLA-A*11:01'],
        'antigen_epitope': ['GILGFVFTL', 'GILGFVFTL', 'NLVPMVATV', 'AVFDRKSDAK'],
        'vdjdb_score': [3, 1, 0, 2]})


def annotate_by_db(repertoire, vdjdb, donor_hla, min_confidence=1):
    # vdjdb_score 0-3: 0 = critical info missing, 3 = independently validated; >=1 drops single-observation noise.
    db = vdjdb[vdjdb['vdjdb_score'] >= min_confidence]
    hits = repertoire.merge(db, on=['cdr3_b_aa', 'v_b_gene'], how='inner', suffixes=('', '_db'))  # V concordance, not CDR3 alone
    carries_hla = hits['mhc_a'].apply(lambda a: any(a.startswith(h) for h in donor_hla))  # restricting HLA must be in the donor
    hits = hits[carries_hla].copy()
    hits['annotation_confidence'] = 'hypothesis'  # a curated match, never a specificity call
    return hits


def beta_neighborhoods(clone_df, radius=50):
    # radius in TCRdist units; ~50 is a common meta-clonotype inclusion radius (Mayer-Blackwell 2021 eLife 10:e68605).
    from tcrdist.repertoire import TCRrep
    tr = TCRrep(cell_df=clone_df, organism='human', chains=['beta'], db_file='alphabeta_gammadelta_db.tsv')
    neighbors = [set(np.where(row <= radius)[0]) for row in tr.pw_beta]
    return tr, neighbors


def load_beta_pgen_model(model_dir):
    import os
    import olga.load_model as load_model
    import olga.generation_probability as generation_probability
    gen_data = load_model.GenomicDataVDJ()
    gen_data.load_igor_genomic_data(os.path.join(model_dir, 'model_params.txt'),
                                    os.path.join(model_dir, 'V_gene_CDR3_anchors.csv'),
                                    os.path.join(model_dir, 'J_gene_CDR3_anchors.csv'))
    gen_model = load_model.GenerativeModelVDJ()
    gen_model.load_and_process_igor_model(os.path.join(model_dir, 'model_marginals.txt'))
    return generation_probability.GenerationProbabilityVDJ(gen_model, gen_data)


def main():
    repertoire = synthetic_repertoire()
    vdjdb = synthetic_vdjdb()
    donor_hla = ['HLA-A*02']  # donor carries A*02 only, so B*07 and A*11 records must be dropped

    hits = annotate_by_db(repertoire, vdjdb, donor_hla, min_confidence=1)
    print('Database annotation (hypotheses, V + HLA concordant, score >= 1):')
    print(hits[['cdr3_b_aa', 'v_b_gene', 'mhc_a', 'antigen_epitope', 'vdjdb_score', 'annotation_confidence']].to_string(index=False))
    print('\nNote: the CASSPGTGGYEQYF NLVPMVATV record (score 0) and the non-A*02 records were correctly excluded.')

    clone_df = repertoire.drop_duplicates(subset=['cdr3_b_aa', 'v_b_gene', 'j_b_gene']).reset_index(drop=True)
    try:
        tr, neighbors = beta_neighborhoods(clone_df, radius=50)
        sizes = [len(n) for n in neighbors]
        print('\nTCRdist beta neighborhood sizes (radius 50 units, hypotheses not antigen labels):', sizes)
    except ImportError:
        print('\n[skip] tcrdist3 not installed; clustering step omitted. Install with: pip install tcrdist3')

    try:
        import os
        model_dir = os.path.join(os.path.dirname(__import__('olga').__file__), 'default_models', 'human_T_beta')
        model = load_beta_pgen_model(model_dir)
        for _, r in hits.drop_duplicates(subset=['cdr3_b_aa']).iterrows():
            p = model.compute_aa_CDR3_pgen(r['cdr3_b_aa'], r['v_b_gene'], r['j_b_gene'])
            print('Pgen', r['cdr3_b_aa'], '=', '%.2e' % p, '(high Pgen => recurs and matches by chance, down-weight)')
    except ImportError:
        print('\n[skip] olga not installed; Pgen null omitted. Install with: pip install olga')


if __name__ == '__main__':
    main()
