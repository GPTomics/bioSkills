'''ClinVar query patterns: REST API, local VCF, ClinGen Allele Registry.

Reference: requests 2.31+, cyvcf2 0.30+ | Verify API if version differs.
Handles 2024 v2 XML schema (VariationArchive anchor; tripartite Germline/Somatic/Oncogenicity).
'''
import requests
import time
from cyvcf2 import VCF
import pandas as pd

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
ALLELE_REG = 'https://reg.clinicalgenome.org'

REVIEW_STATUS_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, single submitter': 1,
    'criteria provided, conflicting interpretations': 1,
    'no assertion criteria provided': 0,
    'no classification provided': 0
}


def clinvar_summary(variation_id):
    '''Fetch VCV-level summary by ClinVar VariationID via E-utilities esummary.'''
    r = requests.get(f'{EUTILS}/esummary.fcgi',
                     params={'db': 'clinvar', 'id': variation_id, 'retmode': 'json'},
                     timeout=30)
    r.raise_for_status()
    payload = r.json()['result'].get(str(variation_id), {})
    germline = payload.get('germline_classification', {})
    return {
        'vcv': payload.get('accession'),
        'name': payload.get('title'),
        'germline_class': germline.get('description'),
        'germline_review_status': germline.get('review_status'),
        'germline_last_evaluated': germline.get('last_evaluated'),
        'somatic_clinical_impact': payload.get('clinical_impact_classification', {}).get('description'),
        'oncogenicity': payload.get('oncogenicity_classification', {}).get('description'),
        'star_rating': REVIEW_STATUS_STARS.get(germline.get('review_status', '').lower(), 0)
    }


def car_id(hgvs_g):
    '''Resolve HGVS-g to ClinGen Allele Registry CA ID. Build/transcript-agnostic.'''
    r = requests.put(f'{ALLELE_REG}/allele',
                     headers={'Content-Type': 'text/plain'},
                     data=hgvs_g, timeout=30)
    if not r.ok:
        return None
    payload = r.json()
    at_id = payload.get('@id', '')
    return at_id.rsplit('/', 1)[-1] if at_id else None


def lookup_local_vcf(clinvar_vcf, chrom, pos, ref, alt):
    '''Query a local clinvar.vcf.gz by GRCh38 coords. Returns VCV-level INFO -- NOT condition-specific.

    For condition-stratified analysis use RCV-level XML parsing.
    '''
    vcf = VCF(clinvar_vcf)
    for v in vcf(f'{chrom}:{pos}-{pos}'):
        if v.REF == ref and alt in v.ALT:
            info = v.INFO
            return {
                'allele_id': info.get('ALLELEID'),
                'clnsig': info.get('CLNSIG'),
                'clnsig_conf': info.get('CLNSIGCONF'),
                'clnrevstat': info.get('CLNREVSTAT'),
                'clndn': info.get('CLNDN'),
                'clnvc': info.get('CLNVC'),
                'clnhgvs': info.get('CLNHGVS'),
                'oncdn': info.get('ONCDN'),
                'scidn': info.get('SCIDN')
            }
    return None


def batch_resolve_to_car_then_clinvar(hgvs_list, sleep=0.34):
    '''Resolve HGVS-g list to CA IDs, then query ClinVar via linked VariationID.

    sleep=0.34s -> ~3 req/s to stay under default NCBI rate limit without API key.
    '''
    rows = []
    for hgvs in hgvs_list:
        ca = car_id(hgvs)
        if ca is None:
            rows.append({'hgvs': hgvs, 'ca_id': None, 'vcv': None})
            time.sleep(sleep)
            continue
        time.sleep(sleep)
        car_payload = requests.get(f'{ALLELE_REG}/allele/{ca}', timeout=30).json()
        externals = car_payload.get('externalRecords', {})
        clinvar_records = externals.get('ClinVarVariations', [])
        variation_id = clinvar_records[0]['variationId'] if clinvar_records else None
        record = {'hgvs': hgvs, 'ca_id': ca, 'variation_id': variation_id}
        if variation_id is not None:
            record.update(clinvar_summary(variation_id))
            time.sleep(sleep)
        rows.append(record)
    return pd.DataFrame(rows)


def filter_by_star_and_freshness(df, min_star=2, max_age_months=36):
    '''Keep variants with star >= min_star AND last_evaluated within max_age_months.

    Star rationale: ClinGen SVI operational guidance treats star>=2 as acceptable for
    clinical action without further review. Star=3 (VCEP) supersedes all lower-star.
    Freshness rationale: Yauy 2022 documents ~1247 classification changes per monthly release.
    '''
    df = df.copy()
    df['last_evaluated_date'] = pd.to_datetime(df['germline_last_evaluated'], errors='coerce')
    today = pd.Timestamp.today()
    df['age_months'] = (today - df['last_evaluated_date']).dt.days / 30.44
    keep = (df['star_rating'] >= min_star) & (df['age_months'] <= max_age_months)
    return df[keep]


def parse_clnsig_conflict(clnsig_conf):
    '''Split CLNSIGCONF into individual conflicting calls.

    The clinical meaning of "Conflicting interpretations" depends on which calls conflict:
    - P vs LP : usually immaterial (both clinically actionable)
    - P vs VUS : clinically meaningful
    - P vs B/LB : major conflict; requires resolution
    '''
    if clnsig_conf is None:
        return []
    calls = []
    for entry in str(clnsig_conf).split('|'):
        if '(' in entry:
            label, count = entry.rsplit('(', 1)
            calls.append({'label': label.strip('_'), 'submitter_count': int(count.rstrip(')'))})
    pathogenic_codes = {'Pathogenic', 'Likely_pathogenic'}
    benign_codes = {'Benign', 'Likely_benign'}
    has_path = any(c['label'] in pathogenic_codes for c in calls)
    has_benign = any(c['label'] in benign_codes for c in calls)
    has_vus = any(c['label'] == 'Uncertain_significance' for c in calls)
    if has_path and has_benign:
        conflict_severity = 'severe'
    elif has_path and has_vus:
        conflict_severity = 'meaningful'
    elif has_path:
        conflict_severity = 'minor'
    else:
        conflict_severity = 'non_pathogenic_only'
    return {'calls': calls, 'severity': conflict_severity}


def annotate_vcf_with_bcftools(input_vcf, clinvar_vcf, output_vcf):
    '''Emit the bcftools command to annotate a user VCF with 2024-schema INFO fields.'''
    return (
        f'bcftools annotate '
        f'-a {clinvar_vcf} '
        f'-c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVC,INFO/CLNHGVS,'
        f'INFO/CLNSIGCONF,INFO/ONCDN,INFO/SCIDN '
        f'{input_vcf} -O z -o {output_vcf}'
    )


if __name__ == '__main__':
    example_hgvs = 'NC_000017.11:g.43094464G>A'
    ca = car_id(example_hgvs)
    print(f'CA ID for {example_hgvs}: {ca}')
    summary = clinvar_summary(17661)
    print(f'Germline classification: {summary["germline_class"]} ({summary["star_rating"]} stars)')
    print(f'Last evaluated: {summary["germline_last_evaluated"]}')
