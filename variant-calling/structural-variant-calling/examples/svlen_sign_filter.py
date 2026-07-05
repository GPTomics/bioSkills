#!/usr/bin/env python3
# Reference: Python 3.10+ stdlib | Verify API if version differs
'''Demonstrate the SVLEN sign-convention pitfall in SV VCF filtering.

Deletions carry a NEGATIVE SVLEN by the historical VCF convention (Manta, DELLY),
so a naive `SVLEN >= min_size` filter silently drops every deletion. The correct
filter uses abs(SVLEN). This script shows both outcomes on a small in-memory VCF
so the difference is concrete, and needs no external tools to run.'''

import re

MIN_SV_SIZE = 50  # bp; VCF/Manta/DELLY convention: SVs are >=50 bp, smaller events are indels

VCF_RECORDS = [
    ('chr1', 1000, 'DEL_big',   'DEL', -8200),   # real deletion, negative SVLEN
    ('chr1', 20000, 'DEL_small', 'DEL', -30),     # 30 bp deletion, below the 50 bp SV floor
    ('chr2', 5000, 'INS_big',   'INS', 4100),     # insertion, positive SVLEN
    ('chr2', 9000, 'DUP_big',   'DUP', 15000),    # duplication, positive SVLEN
    ('chr3', 300, 'DEL_huge',   'DEL', -120000),  # large deletion, negative SVLEN
]


def parse_svlen(info):
    match = re.search(r'SVLEN=(-?\d+)', info)
    return int(match.group(1)) if match else 0


def build_vcf_lines():
    lines = []
    for chrom, pos, sv_id, svtype, svlen in VCF_RECORDS:
        end = pos + abs(svlen) if svtype != 'INS' else pos
        info = f'SVTYPE={svtype};SVLEN={svlen};END={end}'
        lines.append(f'{chrom}\t{pos}\t{sv_id}\tN\t<{svtype}>\t.\tPASS\t{info}')
    return lines


def filter_naive(lines):
    return [l for l in lines if parse_svlen(l.split('\t')[7]) >= MIN_SV_SIZE]


def filter_correct(lines):
    return [l for l in lines if abs(parse_svlen(l.split('\t')[7])) >= MIN_SV_SIZE]


lines = build_vcf_lines()
naive = filter_naive(lines)
correct = filter_correct(lines)

naive_ids = {l.split('\t')[2] for l in naive}
correct_ids = {l.split('\t')[2] for l in correct}

print(f'Total records: {len(lines)}')
print(f"Naive  'SVLEN >= {MIN_SV_SIZE}' kept: {sorted(naive_ids)}")
print(f"Correct 'abs(SVLEN) >= {MIN_SV_SIZE}' kept: {sorted(correct_ids)}")
print(f'Deletions silently dropped by the naive filter: {sorted(correct_ids - naive_ids)}')

# The correct filter drops only DEL_small (30 bp, genuinely sub-SV); the naive filter also
# loses every real deletion (DEL_big, DEL_huge) purely from the negative-sign convention.
assert naive_ids == {'INS_big', 'DUP_big'}, 'naive filter should keep only positive-SVLEN records'
assert correct_ids == {'DEL_big', 'INS_big', 'DUP_big', 'DEL_huge'}, 'abs() filter should keep all >=50 bp events'
print('\nEquivalent bcftools filter: bcftools view -i \'ABS(SVLEN) >= 50\' svs.vcf')
