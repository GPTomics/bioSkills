#!/usr/bin/env python3
'''Research-triage prioritization of variants for clinical review.

This ranks candidates for a curator to review; it is NOT an ACMG classification.
ClinVar significance is treated as a LEAD (review status carried separately), and a
SINGLE calibrated predictor (REVEL) supplies computational evidence at PP3/BP4
supporting strength -- stacking correlated predictors would fake independence and
over-call pathogenic (Pejaver 2022, AJHG 109:2163).
'''
# Reference: cyvcf2 0.30+, bcftools 1.19+ | Verify API if version differs

from cyvcf2 import VCF
import csv
import sys

def prioritize_variant(v):
    '''Triage score from ClinVar leads, grpmax filtering AF, and one calibrated predictor.'''
    score = 0
    reasons = []

    # ClinVar assertion is a LEAD, not evidence: Pathogenic (+10) > Likely pathogenic (+8).
    clnsig = str(v.INFO.get('CLNSIG', ''))
    if 'Pathogenic' in clnsig and 'Likely' not in clnsig:
        score += 10
        reasons.append('ClinVar_Pathogenic_lead')
    elif 'Likely_pathogenic' in clnsig:
        score += 8
        reasons.append('ClinVar_Likely_Pathogenic_lead')

    # Frequency is per-disease; grpmax filtering AF (lower 95% CI), not global AF (Whiffin 2017).
    # <0.0001 dominant-plausible (+5), <0.01 recessive-plausible (+3); absent treated as rare.
    faf = v.INFO.get('fafmax_faf95_max', v.INFO.get('gnomAD_AF', 0.0)) or 0.0
    if faf < 0.0001:
        score += 5
        reasons.append('Rare_FAF<0.0001')
    elif faf < 0.01:
        score += 3
        reasons.append('Uncommon_FAF<0.01')

    # ONE calibrated missense predictor. REVEL PP3_Supporting >=0.644, BP4_Supporting <=0.290
    # (Pejaver 2022 well-reproduced supporting thresholds; higher strengths from the supplement).
    revel = v.INFO.get('REVEL', None)
    if revel is not None and revel >= 0.644:
        score += 2
        reasons.append('REVEL_PP3supporting')
    elif revel is not None and revel <= 0.290:
        score -= 2
        reasons.append('REVEL_BP4supporting')

    # LoF is a PVS1 candidate only after the Abou Tayoun tree (mechanism + NMD + MANE);
    # here it just raises triage priority for curator review.
    consequence = str(v.INFO.get('Consequence', v.INFO.get('ANN', '')))
    if 'stop_gained' in consequence or 'frameshift' in consequence:
        score += 5
        reasons.append('LoF_candidate_PVS1_needs_review')
    elif 'missense' in consequence:
        score += 2
        reasons.append('Missense')

    return score, reasons

def filter_clinical_variants(vcf_path, output_path, min_score=5):
    '''Write triage-ranked variants scoring at or above min_score, carrying review status.'''
    vcf = VCF(vcf_path)
    results = []
    for v in vcf:
        score, reasons = prioritize_variant(v)
        if score >= min_score:
            results.append({
                'chrom': v.CHROM, 'pos': v.POS, 'ref': v.REF, 'alt': ','.join(v.ALT),
                'score': score, 'reasons': ';'.join(reasons),
                'clnsig': v.INFO.get('CLNSIG', '.'),
                'clnrevstat': v.INFO.get('CLNREVSTAT', '.'),
                'gene': v.INFO.get('SYMBOL', v.INFO.get('Gene', '.'))
            })
    vcf.close()

    results.sort(key=lambda x: x['score'], reverse=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys() if results else [], delimiter='\t')
        writer.writeheader()
        writer.writerows(results)

    print(f'Found {len(results)} variants with triage score >= {min_score}')
    print(f'Results written to {output_path}')
    return results

if __name__ == '__main__':
    vcf_path = sys.argv[1] if len(sys.argv) > 1 else 'annotated.vcf.gz'
    output_path = sys.argv[2] if len(sys.argv) > 2 else 'clinical_variants.tsv'
    min_score = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    filter_clinical_variants(vcf_path, output_path, min_score)
