'''dbSNP Build 156 query patterns: rsID resolution, merge-chain following, SPDI conversion.

Reference: requests 2.31+, myvariant 1.0+ | Verify Variation Services API if Build differs.
Build 156 (Sept 2022) is the current schema; pre-Build-156 E-utilities returns thin legacy summary.
'''
import requests
import time
import myvariant
import pandas as pd

VARSVC = 'https://api.ncbi.nlm.nih.gov/variation/v0'

REFSEQ_GRCH38 = {
    '1': 'NC_000001.11', '2': 'NC_000002.12', '3': 'NC_000003.12', '4': 'NC_000004.12',
    '5': 'NC_000005.10', '6': 'NC_000006.12', '7': 'NC_000007.14', '8': 'NC_000008.11',
    '9': 'NC_000009.12', '10': 'NC_000010.11', '11': 'NC_000011.10', '12': 'NC_000012.12',
    '13': 'NC_000013.11', '14': 'NC_000014.9', '15': 'NC_000015.10', '16': 'NC_000016.10',
    '17': 'NC_000017.11', '18': 'NC_000018.10', '19': 'NC_000019.10', '20': 'NC_000020.11',
    '21': 'NC_000021.9', '22': 'NC_000022.11', 'X': 'NC_000023.11', 'Y': 'NC_000024.10'
}


def refsnp(rsid, sleep=0.34):
    '''Fetch full Build 156 RefSNP JSON. sleep=0.34s -> ~3 req/s without API key.'''
    rs_int = str(rsid).lstrip('rs')
    r = requests.get(f'{VARSVC}/refsnp/{rs_int}', timeout=30)
    time.sleep(sleep)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def resolve_merge_chain(rsid, max_hops=10):
    '''Follow multi-hop merge chain. Cycle-safe with max_hops cap.

    Multi-hop merges (rs3 -> rs2 -> rs1 -> rs0) appear as multiple RsMergeArch rows;
    Variation Services handles this via merged_snapshot_data, but rare edge cases
    require explicit traversal.
    '''
    seen = set()
    current = str(rsid).lstrip('rs')
    for _ in range(max_hops):
        if current in seen:
            return {'error': 'merge_cycle', 'chain': list(seen)}
        seen.add(current)
        payload = refsnp(current)
        if payload is None:
            return {'status': 'not_found', 'final_rsid': current, 'chain': list(seen)}
        if payload.get('is_withdrawn'):
            return {'status': 'withdrawn', 'final_rsid': current, 'chain': list(seen),
                    'reason': payload.get('withdrawn_release', {})}
        if payload.get('primary_snapshot_data') is not None:
            return {'status': 'resolved', 'final_rsid': payload.get('refsnp_id'),
                    'chain': list(seen)}
        merged = payload.get('merged_snapshot_data', [])
        if not merged:
            return {'status': 'orphan', 'final_rsid': current, 'chain': list(seen)}
        current = str(merged[0].get('merged_into', ''))
    return {'error': 'hop_limit_exceeded', 'chain': list(seen)}


def alleles_grch38(payload):
    '''Extract alleles from RefSNP JSON for GRCh38 only. Returns list of {ref, alt, spdi}.

    A single rsID with len(alleles) > 2 is multi-allelic -- 6-8% of dbSNP rsIDs.
    Naive rsID -> variant mappings break for these sites.
    '''
    if payload is None or payload.get('is_withdrawn'):
        return []
    primary = payload.get('primary_snapshot_data', {})
    out = []
    for placement in primary.get('placements_with_allele', []):
        seq_id_traits = placement.get('seq_id_traits_by_assembly', [{}])[0]
        if 'GRCh38' not in seq_id_traits.get('assembly_name', ''):
            continue
        for allele in placement.get('alleles', []):
            spdi = allele.get('allele', {}).get('spdi', {})
            if spdi.get('inserted_sequence') == spdi.get('deleted_sequence'):
                continue
            out.append({
                'ref': spdi.get('deleted_sequence'),
                'alt': spdi.get('inserted_sequence'),
                'spdi': f"{spdi.get('seq_id')}:{spdi.get('position')}:{spdi.get('deleted_sequence')}:{spdi.get('inserted_sequence')}",
                'pos_0based': spdi.get('position')
            })
    return out


def vcf_to_canonical_spdi(chrom, pos, ref, alt):
    '''VCF (1-based) -> canonical right-aligned SPDI (0-based).

    The Variant Overprecision Correction Algorithm normalizes left/right-aligned
    representations to a single canonical form -- required for cross-database joins.
    '''
    refseq = REFSEQ_GRCH38.get(str(chrom).lstrip('chr'))
    if refseq is None:
        return None
    raw_spdi = f'{refseq}:{pos - 1}:{ref}:{alt}'
    r = requests.get(f'{VARSVC}/spdi/{raw_spdi}/canonical_representative', timeout=30)
    time.sleep(0.34)
    if not r.ok:
        return None
    return r.json().get('data', {}).get('spdi')


def spdi_to_rsid(spdi_str):
    '''SPDI -> rsID if a cluster exists at that position with matching allele.'''
    r = requests.get(f'{VARSVC}/spdi/{spdi_str}/rsids', timeout=30)
    time.sleep(0.34)
    if not r.ok:
        return None
    rsids = r.json().get('data', {}).get('rsids', [])
    return rsids[0] if rsids else None


def hgvs_to_canonical_spdi(hgvs):
    '''HGVS-g / HGVS-c -> canonical SPDI via Variation Services contextuals.'''
    r = requests.get(f'{VARSVC}/hgvs/{hgvs}/contextuals', timeout=30)
    time.sleep(0.34)
    if not r.ok:
        return None
    contextuals = r.json().get('data', {}).get('spdis', [])
    return contextuals[0] if contextuals else None


def batch_normalize_rsids(rsids):
    '''Normalize a list of rsIDs: resolve merges, flag multi-allelic, return DataFrame.'''
    mv = myvariant.MyVariantInfo()
    fields = ['dbsnp', 'gnomad_exome.af.af', 'gnomad_genome.af.af', 'clinvar.clinical_significance']
    aggregated = mv.getvariants(rsids, fields=fields)
    rows = []
    for entry in aggregated:
        rsid = entry.get('query')
        merge_result = resolve_merge_chain(rsid)
        alleles = alleles_grch38(refsnp(merge_result.get('final_rsid', rsid)))
        rows.append({
            'input_rsid': rsid,
            'canonical_rsid': merge_result.get('final_rsid'),
            'status': merge_result.get('status'),
            'chain_length': len(merge_result.get('chain', [])),
            'is_multiallelic': len(alleles) > 1 and len({a['ref'] for a in alleles}) > 1,
            'n_alleles': len(alleles),
            'gnomad_exome_af': entry.get('gnomad_exome', {}).get('af', {}).get('af'),
            'gnomad_genome_af': entry.get('gnomad_genome', {}).get('af', {}).get('af'),
            'clinvar_sig': entry.get('clinvar', {}).get('clinical_significance')
        })
    return pd.DataFrame(rows)


def alfa_population_frequencies(rsid):
    '''Extract ALFA per-population AFs from RefSNP JSON.

    ALFA aggregates dbGaP studies across 12 ancestry groups. Use when:
    - Variant is array-genotyped (gnomAD may miss it)
    - Consent-respecting frequency lookup needed
    Do NOT use for rare-variant FAF95 -- use gnomAD instead.
    '''
    payload = refsnp(rsid)
    if payload is None:
        return None
    freq_records = (payload.get('primary_snapshot_data', {})
                    .get('allele_annotations', [{}])[0]
                    .get('frequency', []))
    alfa = [f for f in freq_records if 'ALFA' in f.get('study_name', '')]
    out = {}
    for record in alfa:
        ancestry = record.get('common_name', 'Unknown')
        total = record.get('total_count')
        count = record.get('allele_count')
        out[ancestry] = {
            'allele': record.get('observation', {}).get('inserted_sequence'),
            'count': count,
            'total': total,
            'af': count / total if total else None
        }
    return out


if __name__ == '__main__':
    apoe_e4 = 'rs429358'
    chain = resolve_merge_chain(apoe_e4)
    print(f'rs429358 (APOE e4): status={chain["status"]}, chain_length={len(chain.get("chain", []))}')
    payload = refsnp(chain['final_rsid'])
    alleles = alleles_grch38(payload)
    print(f'GRCh38 alleles: {alleles}')

    spdi = vcf_to_canonical_spdi('17', 43094464, 'G', 'A')
    print(f'BRCA1 chr17:43094464:G>A canonical SPDI: {spdi}')
    rsid = spdi_to_rsid(spdi)
    print(f'  -> rsID: {rsid}')
