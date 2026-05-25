'''KEGG Orthology (KO) lookup: gene -> KO -> all orthologous genes across species.'''
# Reference: requests 2.31+ | Verify API if version differs
import requests
import time

BASE = 'https://rest.kegg.jp'
SLEEP = 0.3  # KEGG: be polite; no strict published limit


def ko_for_gene(species, gene):
    '''species: KEGG code (hsa=human, mmu=mouse, dme=fly). gene: NCBI Gene ID or KEGG locus.'''
    r = requests.get(f'{BASE}/link/ko/{species}:{gene}')
    r.raise_for_status()
    return [line.split('\t')[1].replace('ko:', '') for line in r.text.strip().split('\n') if line]


def genes_for_ko(ko_id):
    '''All KEGG genes annotated with this KO.'''
    r = requests.get(f'{BASE}/link/genes/{ko_id}')
    r.raise_for_status()
    return [line.split('\t')[1] for line in r.text.strip().split('\n') if line]


def ko_info(ko_id):
    '''Description + pathway list for a KO.'''
    r = requests.get(f'{BASE}/get/{ko_id}')
    r.raise_for_status()
    return r.text


print('=== Human TP53 (Gene ID 7157) -> KEGG Orthology ===')
ko_ids = ko_for_gene('hsa', '7157')
print(f'KOs: {ko_ids}')
time.sleep(SLEEP)

if ko_ids:
    print(f'\n=== KO {ko_ids[0]} info (truncated) ===')
    info = ko_info(ko_ids[0])
    for line in info.split('\n')[:15]:
        print(f'  {line}')
    time.sleep(SLEEP)

    print(f'\n=== Members of {ko_ids[0]} across all species ===')
    members = genes_for_ko(ko_ids[0])
    print(f'Total members: {len(members)}')
    species_counts = {}
    for m in members:
        sp = m.split(':')[0]
        species_counts[sp] = species_counts.get(sp, 0) + 1
    print(f'Species represented: {len(species_counts)}')
    print('Top species by member count:')
    for sp, count in sorted(species_counts.items(), key=lambda x: -x[1])[:10]:
        print(f'  {sp:>6}: {count}')

print('\n=== License note ===')
print('KEGG REST API: free for academic use; commercial use requires paid license.')
print('For commercial workflows, prefer eggNOG (permissive) or OrthoDB.')
