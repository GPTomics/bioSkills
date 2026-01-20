#!/usr/bin/env python3
'''Filter and prioritize variants for clinical interpretation'''

from cyvcf2 import VCF
import csv

def prioritize_variant(v):
    '''Score variant based on clinical criteria. Weights reflect ACMG evidence strength.'''
    score = 0
    reasons = []

    # ClinVar: Pathogenic (+10) > Likely pathogenic (+8) per ACMG strong evidence
    clnsig = v.INFO.get('CLNSIG', '')
    if 'Pathogenic' in str(clnsig):
        score += 10
        reasons.append('ClinVar_Pathogenic')
    elif 'Likely_pathogenic' in str(clnsig):
        score += 8
        reasons.append('ClinVar_Likely_Pathogenic')

    # AF: <0.001 rare (+5), <0.01 uncommon (+3); thresholds per gnomAD guidance
    af = v.INFO.get('gnomAD_AF', v.INFO.get('AF', 1.0))
    if af is not None and af < 0.001:
        score += 5
        reasons.append('Rare_AF<0.001')
    elif af is not None and af < 0.01:
        score += 3
        reasons.append('Uncommon_AF<0.01')

    # CADD: >25 high deleterious (~0.3% most deleterious), >20 moderate (~1%)
    cadd = v.INFO.get('CADD_PHRED', 0)
    if cadd and cadd > 25:
        score += 4
        reasons.append('HighCADD>25')
    elif cadd and cadd > 20:
        score += 2
        reasons.append('ModerateCADD>20')

    # Consequence: LoF (+5) > missense (+2); LoF generally more disruptive
    consequence = v.INFO.get('Consequence', v.INFO.get('ANN', ''))
    if 'stop_gained' in str(consequence) or 'frameshift' in str(consequence):
        score += 5
        reasons.append('LoF')
    elif 'missense' in str(consequence):
        score += 2
        reasons.append('Missense')

    return score, reasons

def filter_clinical_variants(vcf_path, output_path, min_score=5):
    '''Filter VCF for clinically relevant variants'''
    vcf = VCF(vcf_path)

    results = []
    for v in vcf:
        score, reasons = prioritize_variant(v)

        if score >= min_score:
            results.append({
                'chrom': v.CHROM,
                'pos': v.POS,
                'ref': v.REF,
                'alt': ','.join(v.ALT),
                'score': score,
                'reasons': ';'.join(reasons),
                'clnsig': v.INFO.get('CLNSIG', '.'),
                'gene': v.INFO.get('SYMBOL', v.INFO.get('Gene', '.'))
            })

    vcf.close()

    results.sort(key=lambda x: x['score'], reverse=True)

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys() if results else [], delimiter='\t')
        writer.writeheader()
        writer.writerows(results)

    print(f'Found {len(results)} variants with score >= {min_score}')
    print(f'Results written to {output_path}')

    return results

if __name__ == '__main__':
    import sys
    vcf_path = sys.argv[1] if len(sys.argv) > 1 else 'annotated.vcf.gz'
    output_path = sys.argv[2] if len(sys.argv) > 2 else 'clinical_variants.tsv'
    min_score = int(sys.argv[3]) if len(sys.argv) > 3 else 5

    filter_clinical_variants(vcf_path, output_path, min_score)
