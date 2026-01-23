'''Query dbSNP for rsID information and coordinate mapping'''

import myvariant
import pandas as pd

mv = myvariant.MyVariantInfo()

def rsid_lookup(rsid):
    '''Look up rsID and return key annotations'''
    result = mv.getvariant(rsid, fields=['dbsnp', 'clinvar.clinical_significance', 'gnomad_exome.af.af'])
    if not result:
        return None

    dbsnp = result.get('dbsnp', {})
    return {
        'rsid': dbsnp.get('rsid'),
        'chrom': dbsnp.get('chrom'),
        'pos_hg38': dbsnp.get('hg38', {}).get('start') if isinstance(dbsnp.get('hg38'), dict) else None,
        'pos_hg19': dbsnp.get('hg19', {}).get('start') if isinstance(dbsnp.get('hg19'), dict) else None,
        'ref': dbsnp.get('ref'),
        'alt': dbsnp.get('alt'),
        'gene': dbsnp.get('gene', {}).get('symbol') if isinstance(dbsnp.get('gene'), dict) else None,
        'class': dbsnp.get('class'),
        'clinvar': result.get('clinvar', {}).get('clinical_significance'),
        'gnomad_af': result.get('gnomad_exome', {}).get('af', {}).get('af')
    }

# Single rsID lookup
rsid = 'rs121913527'
info = rsid_lookup(rsid)
print(f'rsID lookup: {rsid}')
if info:
    print(f"  Location: chr{info['chrom']}:{info['pos_hg38']} (GRCh38)")
    print(f"  Change: {info['ref']}>{info['alt']}")
    print(f"  Gene: {info['gene']}")
    print(f"  Class: {info['class']}")
    print(f"  ClinVar: {info['clinvar']}")
    print(f"  gnomAD AF: {info['gnomad_af']}")

# Batch rsID lookup
rsids = ['rs121913527', 'rs1800566', 'rs104894155', 'rs28897696']
print(f'\nBatch lookup of {len(rsids)} rsIDs:')

results = mv.getvariants(rsids, fields=['dbsnp.rsid', 'dbsnp.chrom', 'dbsnp.hg38', 'dbsnp.gene.symbol'])

records = []
for r in results:
    dbsnp = r.get('dbsnp', {})
    hg38 = dbsnp.get('hg38', {})
    gene = dbsnp.get('gene', {})
    records.append({
        'query': r.get('query'),
        'rsid': dbsnp.get('rsid'),
        'chrom': dbsnp.get('chrom'),
        'pos': hg38.get('start') if isinstance(hg38, dict) else None,
        'gene': gene.get('symbol') if isinstance(gene, dict) else None
    })

df = pd.DataFrame(records)
print(df.to_string(index=False))

# Map coordinates to rsID
def coords_to_rsid(chrom, pos, ref, alt):
    '''Find rsID for a genomic coordinate'''
    hgvs = f'chr{chrom}:g.{pos}{ref}>{alt}'
    result = mv.getvariant(hgvs, fields=['dbsnp.rsid'])
    if result:
        return result.get('dbsnp', {}).get('rsid')
    return None

print('\nCoordinate to rsID mapping:')
test_coords = [
    ('7', 140453136, 'A', 'T'),   # BRAF V600E
    ('17', 7577120, 'C', 'T'),    # TP53
]
for chrom, pos, ref, alt in test_coords:
    rsid = coords_to_rsid(chrom, pos, ref, alt)
    print(f'  chr{chrom}:{pos}{ref}>{alt} -> {rsid or "not found"}')
