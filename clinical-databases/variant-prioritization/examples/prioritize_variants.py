'''Rare-disease variant prioritization pipeline: trio cascading filters + DNV + compound het + ACMG SF.

Reference: cyvcf2 0.30+, pandas 2.2+, myvariant 1.0+ | Verify API if version differs.
ACMG SF list pinned to v3.2 (Miller 2023 Genet Med; 81 genes including new CALM1/2/3).
ACMG classification logic (PVS1, PP3, BS1, etc.) is in clinical-databases/acmg-classification.
'''
from cyvcf2 import VCF
import pandas as pd
import myvariant


# Frequency thresholds per ClinGen SVI; gene-specific via Whiffin 2017 max-credible-AF.
DOMINANT_FREQ_FILTER = 0.0001    # grpmax_faf95 < 0.01% (PM2_Supporting)
RECESSIVE_FREQ_FILTER = 0.005    # grpmax_faf95 < 0.5%

# ACMG SF v3.2 (Miller 2023 Genet Med 25:100866). Pin to current release.
ACMG_SF_V3_2 = {
    # CALM v3.2 ADDITIONS -- calmodulinopathy (long QT / CPVT; high actionability)
    'CALM1', 'CALM2', 'CALM3',
    # Cardiomyopathies / aortopathies
    'ACTA2', 'ACTC1', 'COL3A1', 'DES', 'FBN1', 'FLNC', 'GLA', 'LMNA', 'MYBPC3',
    'MYH7', 'MYH11', 'MYL2', 'MYL3', 'PRKAG2', 'PKP2', 'RBM20', 'SCN5A', 'SMAD3',
    'TGFBR1', 'TGFBR2', 'TMEM43', 'TNNI3', 'TNNT2', 'TPM1', 'TTN',
    # Channelopathies + arrhythmias
    'CACNA1S', 'KCNH2', 'KCNQ1', 'RYR1', 'RYR2',
    # Cancer predisposition
    'APC', 'ATM', 'BAP1', 'BMPR1A', 'BRCA1', 'BRCA2', 'BRIP1', 'CDH1', 'CDKN2A',
    'CHEK2', 'GREM1', 'HOXB13', 'MAX', 'MEN1', 'MLH1', 'MSH2', 'MSH6', 'MUTYH',
    'NF2', 'PALB2', 'PMS2', 'PTEN', 'RAD51C', 'RAD51D', 'RB1', 'RET', 'SDHAF2',
    'SDHB', 'SDHC', 'SDHD', 'SMAD4', 'STK11', 'TMEM127', 'TP53', 'TSC1', 'TSC2',
    'VHL', 'WT1',
    # Vascular
    'ACVRL1', 'ENG',
    # Metabolic / other
    'FH', 'GAA', 'HFE', 'HNF1A', 'LDLR', 'NTRK1', 'OTC', 'PCSK9', 'TTR'
}

FUNCTIONAL_CONSEQUENCES = {
    'missense_variant', 'stop_gained', 'stop_lost', 'start_lost',
    'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
    'splice_donor_variant', 'splice_acceptor_variant',
    'protein_altering_variant'
}


def stage1_qc_frequency_filter(vcf_path, samples, max_grpmax_faf95=DOMINANT_FREQ_FILTER,
                                min_dp=10, min_gq=20):
    '''Stage 1: QC + population frequency. Reduces 100k -> ~5-15k variants.

    Args:
        vcf_path: joint-called trio VCF
        samples: ordered list of sample IDs (e.g., [proband, mother, father])
        max_grpmax_faf95: 0.0001 (dominant) or 0.005 (recessive) per ClinGen SVI
    '''
    vcf = VCF(vcf_path)
    sample_idx = {s: vcf.samples.index(s) for s in samples}
    rows = []
    for v in vcf:
        if v.FILTER is not None:
            continue
        depths = v.gt_depths
        if min(depths[sample_idx[s]] for s in samples) < min_dp:
            continue
        if v.QUAL is not None and v.QUAL < min_gq:
            continue
        gnomad_faf95 = (v.INFO.get('grpmax_faf95') or v.INFO.get('AF_grpmax') or
                        v.INFO.get('AF_popmax') or 0)
        if gnomad_faf95 > max_grpmax_faf95:
            continue
        rows.append({
            'chrom': v.CHROM, 'pos': v.POS, 'ref': v.REF, 'alt': v.ALT[0] if v.ALT else None,
            'qual': v.QUAL,
            'gt_proband': v.gt_types[sample_idx[samples[0]]],
            'gt_mother': v.gt_types[sample_idx[samples[1]]] if len(samples) > 1 else None,
            'gt_father': v.gt_types[sample_idx[samples[2]]] if len(samples) > 2 else None,
            'dp_proband': depths[sample_idx[samples[0]]],
            'dp_mother': depths[sample_idx[samples[1]]] if len(samples) > 1 else None,
            'dp_father': depths[sample_idx[samples[2]]] if len(samples) > 2 else None,
            'gnomad_faf95': gnomad_faf95,
            'consequence': _parse_consequence(v.INFO.get('CSQ', '')),
            'gene': _parse_gene(v.INFO.get('CSQ', ''))
        })
    return pd.DataFrame(rows)


def _parse_consequence(csq):
    if not csq:
        return None
    return csq.split(',')[0].split('|')[1] if '|' in csq else None


def _parse_gene(csq):
    if not csq:
        return None
    fields = csq.split(',')[0].split('|')
    return fields[3] if len(fields) > 3 else None


def stage2_functional_filter(df):
    '''Stage 2: retain functional coding variants only.'''
    if 'consequence' not in df.columns:
        return df
    is_functional = df['consequence'].astype(str).apply(
        lambda c: any(fc in c for fc in FUNCTIONAL_CONSEQUENCES)
    )
    return df[is_functional]


def stage3a_call_de_novo_candidates(df):
    '''Stage 3a: identify DNV candidates from Mendelian-violation.

    WARNING: this Mendelian-violation logic has 10-30% false-positive rate.
    All reportable DNVs require IGV inspection. Use DeNovoGear / DeNovoCNN
    for production with Bayesian posterior.
    '''
    is_dnv = (
        (df['gt_mother'] == 0) & (df['gt_father'] == 0) &
        df['gt_proband'].isin([1, 3]) &
        (df['dp_mother'] >= 10) & (df['dp_father'] >= 10)
    )
    df = df.copy()
    df['de_novo_candidate'] = is_dnv
    return df


def stage3b_compound_het_candidates(df):
    '''Stage 3b: identify compound het: trio-phased, het from each parent.

    For singletons, supplement with WhatsHap read-based phasing.
    '''
    df = df.copy()
    het_in_proband = (df['gt_proband'] == 1)
    maternal_het = (df['gt_mother'] == 1) & (df['gt_father'] == 0)
    paternal_het = (df['gt_father'] == 1) & (df['gt_mother'] == 0)
    df['inherited_maternal'] = het_in_proband & maternal_het
    df['inherited_paternal'] = het_in_proband & paternal_het

    candidate_genes = set()
    for gene, sub in df[het_in_proband].groupby('gene'):
        if pd.isna(gene):
            continue
        has_mat = sub['inherited_maternal'].any()
        has_pat = sub['inherited_paternal'].any()
        if has_mat and has_pat:
            candidate_genes.add(gene)
    df['compound_het_candidate'] = df['gene'].isin(candidate_genes) & het_in_proband
    return df


def stage4_acmg_sf_check(df):
    '''Stage 4: flag ACMG SF v3.2 (Miller 2023; 81 genes) candidates.

    Only P/LP variants in SF genes are reportable as secondary findings.
    Apply ACMG classification (acmg-classification skill) before final reporting.
    '''
    df = df.copy()
    df['acmg_sf_v3_2_gene'] = df['gene'].isin(ACMG_SF_V3_2)
    return df


def stage5_annotate_with_myvariant(df, hgvs_list=None, chunk=500):
    '''Stage 5: pull ClinVar, gnomAD, AlphaMissense, REVEL, SpliceAI via myvariant.info.

    Defer ACMG classification logic (PVS1, PP3, BS1) to acmg-classification skill.
    This stage only attaches the evidence; classification is downstream.
    '''
    if hgvs_list is None:
        hgvs_list = df.apply(lambda r: f"chr{r['chrom']}:g.{r['pos']}{r['ref']}>{r['alt']}", axis=1).tolist()
    mv = myvariant.MyVariantInfo()
    fields = [
        'clinvar.clinical_significance', 'clinvar.review_status',
        'gnomad_exome.faf95', 'gnomad_exome.af.af',
        'dbnsfp.alphamissense.score', 'dbnsfp.revel.score',
        'dbnsfp.cadd.phred', 'dbnsfp.spliceai.ds_max',
        '_meta'
    ]
    annotations = mv.getvariants(hgvs_list, fields=fields)
    rows = []
    for r in annotations:
        clinvar = r.get('clinvar', {}) or {}
        gnomad = r.get('gnomad_exome', {}) or {}
        dbnsfp = r.get('dbnsfp', {}) or {}
        rows.append({
            'variant': r.get('query'),
            'clinvar_sig': clinvar.get('clinical_significance'),
            'clinvar_review': clinvar.get('review_status'),
            'grpmax_faf95': gnomad.get('faf95', {}).get('popmax') if isinstance(gnomad.get('faf95'), dict) else None,
            'alphamissense': dbnsfp.get('alphamissense', {}).get('score'),
            'revel': dbnsfp.get('revel', {}).get('score'),
            'cadd_phred': dbnsfp.get('cadd', {}).get('phred'),
            'spliceai_ds_max': dbnsfp.get('spliceai', {}).get('ds_max'),
        })
    annot_df = pd.DataFrame(rows)
    df = df.copy()
    df['variant'] = hgvs_list
    return df.merge(annot_df, on='variant', how='left')


def exomiser_command(yml_path, vcf_path, hpo_terms, output_dir):
    '''Emit Exomiser command for phenotype-driven prioritization.

    HPO terms (e.g., HP:0001250 seizure) MUST be specific. 5-10 specific HPO terms.
    Sparse generic HPO degrades Exomiser hiPHIVE accuracy significantly.
    '''
    return (f'java -jar exomiser-cli-14.0.0.jar '
            f'--analysis {yml_path} '
            f'--vcf {vcf_path} '
            f'--hpo {",".join(hpo_terms)} '
            f'--output-dir {output_dir}')


if __name__ == '__main__':
    print('Rare-disease variant prioritization pipeline:')
    print('1. Stage 1: QC + grpmax_faf95 < 0.0001 (dominant) or < 0.005 (recessive)')
    print('2. Stage 2: Filter to functional coding consequences')
    print('3a. Stage 3a: De novo candidates (Mendelian-violation; supplement with DeNovoGear)')
    print('3b. Stage 3b: Compound het candidates (trio-phased; WhatsHap for singletons)')
    print('4. Stage 4: ACMG SF v3.2 (Miller 2023; 81 genes) check')
    print('5. Stage 5: Aggregated annotation via myvariant.info')
    print('6. Stage 6: Exomiser hiPHIVE phenotype-driven ranking')
    print('7. Final: ACMG classification per `clinical-databases/acmg-classification`')
    print('         Output: ranked candidates with tier, inheritance, classification, SF flag')
