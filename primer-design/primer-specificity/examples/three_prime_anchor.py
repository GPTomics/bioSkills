'''Why BLAST is the wrong tool for primer specificity: the 3' terminus, not overall similarity,
decides whether a candidate site primes. This offline demo contrasts an overall-similarity score
(calc_heterodimer dG, what a similarity search tracks) against the 3'-anchor score
(calc_end_stability dG) for candidate sites. The real genome-wide step is pair-aware in-silico PCR
(mfeprimer / isPcr / Primer-BLAST) against the correct database -- see SKILL.md. No file IO, no network.'''
# Reference: primer3-py 2.3+ | Verify API if version differs

import primer3

COMP = str.maketrans('ACGT', 'TGCA')

def revcomp(s):
    return s.translate(COMP)[::-1]

def mutate(s, i):
    return s[:i] + ('A' if s[i] != 'A' else 'C') + s[i + 1:]

primer = 'GTCTCCTCTGACTTCAACAGCG'
site = revcomp(primer)             # the strand the primer anneals to; primer 3' base pairs site[0]

sites = {
    'on-target (perfect)': site,
    'internal mismatch (3prime intact)': mutate(site, len(site) // 2),
    '3prime-terminal mismatch': mutate(site, 0),
}

print('primer:', primer)
print(f'{"site":36s} {"overall dG":>11s} {"3prime dG":>11s}')
for label, s in sites.items():
    overall = primer3.calc_heterodimer(primer, s).dg / 1000     # what overall-similarity tracks
    anchor = primer3.calc_end_stability(primer, s).dg / 1000    # the 3'-anchor BLAST ignores
    print(f'{label:36s} {overall:8.2f}    {anchor:8.2f}  kcal/mol')

print('\nThe 3-prime-mismatch site keeps a strong OVERALL dG (a similarity search would surface it as')
print('a fine hit) but its 3-prime ANCHOR dG collapses -- so it will not prime. Overall similarity is')
print('blind to the 3-prime end; confirm specificity with pair-aware in-silico PCR (genome + transcriptome for RT-qPCR).')
