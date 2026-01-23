'''Query gnomAD for population allele frequencies'''

import myvariant
import pandas as pd

mv = myvariant.MyVariantInfo()

def get_gnomad_frequencies(variant_hgvs):
    '''Get gnomAD exome and genome frequencies via myvariant.info

    Returns dict with global and population-specific frequencies.
    None values indicate variant not observed in that dataset.
    '''
    result = mv.getvariant(variant_hgvs, fields=['gnomad_exome.af', 'gnomad_genome.af'])
    if not result:
        return None

    exome = result.get('gnomad_exome', {}).get('af', {})
    genome = result.get('gnomad_genome', {}).get('af', {})

    return {
        'exome_af': exome.get('af') if isinstance(exome, dict) else None,
        'genome_af': genome.get('af') if isinstance(genome, dict) else None,
        'exome_af_nfe': exome.get('af_nfe') if isinstance(exome, dict) else None,
        'exome_af_afr': exome.get('af_afr') if isinstance(exome, dict) else None,
        'exome_af_eas': exome.get('af_eas') if isinstance(exome, dict) else None,
    }

# Example variants - mix of common and rare
variants = [
    'chr7:g.140453136A>T',   # BRAF V600E - rare somatic, very rare germline
    'rs1800566',             # NQO1 P187S - common variant
    'rs104894155',           # CFTR F508del - rare, pathogenic
]

print('gnomAD frequency lookup:\n')
records = []
for v in variants:
    freqs = get_gnomad_frequencies(v)
    if freqs:
        records.append({'variant': v, **freqs})
        exome_af = freqs['exome_af']
        genome_af = freqs['genome_af']
        # Use max of exome/genome for conservative rare variant filtering
        max_af = max(filter(None, [exome_af, genome_af]), default=None)
        print(f'{v}:')
        print(f'  Exome AF: {exome_af if exome_af else "absent"}')
        print(f'  Genome AF: {genome_af if genome_af else "absent"}')
        print(f'  Max AF: {max_af if max_af else "absent"}')
    else:
        print(f'{v}: Not found in gnomAD')
    print()

# Filtering thresholds for rare disease analysis
# ACMG PM2: Absent or AF < 0.01 (1%) supports pathogenicity
RARE_THRESHOLD = 0.01       # 1% - standard rare disease cutoff
VERY_RARE_THRESHOLD = 0.001 # 0.1% - stringent filtering

def is_rare(af, threshold=RARE_THRESHOLD):
    '''Check if variant is rare (supports ACMG PM2)'''
    if af is None:
        return True  # Absent from gnomAD = rare
    return af < threshold

print('Rarity assessment:')
for r in records:
    max_af = max(filter(None, [r['exome_af'], r['genome_af']]), default=None)
    rare = is_rare(max_af)
    very_rare = is_rare(max_af, VERY_RARE_THRESHOLD)
    print(f"  {r['variant']}: {'RARE' if rare else 'COMMON'} (AF={max_af})")
