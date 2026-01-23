'''Query myvariant.info for variant annotations from multiple databases'''

import myvariant
import pandas as pd

mv = myvariant.MyVariantInfo()

# Single variant query - returns all available annotations
variant = 'chr7:g.140453136A>T'  # BRAF V600E
result = mv.getvariant(variant)
print(f'Query: {variant}')
print(f'ClinVar: {result.get("clinvar", {}).get("clinical_significance", "N/A")}')
print(f'gnomAD AF: {result.get("gnomad_exome", {}).get("af", {}).get("af", "N/A")}')
print(f'CADD phred: {result.get("cadd", {}).get("phred", "N/A")}')

# Batch query with specific fields
variants = [
    'rs121913527',      # BRAF V600E
    'rs1800566',        # NQO1 P187S
    'rs104894155',      # CFTR F508del
    'rs28897696',       # BRCA1
]

# Request only the fields we need - faster and smaller response
results = mv.getvariants(
    variants,
    fields=['clinvar.clinical_significance', 'gnomad_exome.af.af', 'cadd.phred', 'dbsnp.rsid']
)

# Parse results into DataFrame
records = []
for r in results:
    # Handle nested dictionary access safely
    clinvar = r.get('clinvar', {})
    gnomad = r.get('gnomad_exome', {})

    records.append({
        'query': r.get('query'),
        'rsid': r.get('dbsnp', {}).get('rsid'),
        'clinvar_sig': clinvar.get('clinical_significance') if isinstance(clinvar, dict) else None,
        'gnomad_af': gnomad.get('af', {}).get('af') if isinstance(gnomad, dict) else None,
        'cadd_phred': r.get('cadd', {}).get('phred')
    })

df = pd.DataFrame(records)
print('\nBatch query results:')
print(df.to_string(index=False))

# Search for pathogenic variants in a gene
# Useful for finding known disease variants
search_results = mv.query(
    'clinvar.gene.symbol:BRCA1 AND clinvar.clinical_significance:Pathogenic',
    size=10,  # Limit results; increase for comprehensive search
    fields=['clinvar.clinical_significance', 'clinvar.hgvs.coding']
)
print(f'\nFound {search_results["total"]} pathogenic BRCA1 variants in ClinVar')
for hit in search_results['hits'][:3]:
    hgvs = hit.get('clinvar', {}).get('hgvs', {}).get('coding', 'N/A')
    print(f'  {hgvs}')
