'''Validate a primer pair for PCR'''

import primer3

forward = 'GTCTCCTCTGACTTCAACAGCG'
reverse = 'ACCACCCTGTTGCTGTAGCCAA'

print('=== Primer Pair Validation ===\n')
print(f'Forward: {forward}')
print(f'Reverse: {reverse}\n')

fwd_tm = primer3.calc_tm(forward)
rev_tm = primer3.calc_tm(reverse)
print(f'Forward Tm: {fwd_tm:.1f}C')
print(f'Reverse Tm: {rev_tm:.1f}C')
print(f'Tm Difference: {abs(fwd_tm - rev_tm):.1f}C')

fwd_hairpin = primer3.calc_hairpin(forward)
rev_hairpin = primer3.calc_hairpin(reverse)
print(f'\nForward hairpin Tm: {fwd_hairpin.tm:.1f}C (dG: {fwd_hairpin.dg:.0f})')
print(f'Reverse hairpin Tm: {rev_hairpin.tm:.1f}C (dG: {rev_hairpin.dg:.0f})')

fwd_homodimer = primer3.calc_homodimer(forward)
rev_homodimer = primer3.calc_homodimer(reverse)
print(f'\nForward homodimer Tm: {fwd_homodimer.tm:.1f}C')
print(f'Reverse homodimer Tm: {rev_homodimer.tm:.1f}C')

heterodimer = primer3.calc_heterodimer(forward, reverse)
print(f'\nHeterodimer Tm: {heterodimer.tm:.1f}C (dG: {heterodimer.dg:.0f})')

print('\n=== Summary ===')
issues = []
if abs(fwd_tm - rev_tm) > 2:
    issues.append('Tm difference > 2C')
if fwd_hairpin.tm > 45 or rev_hairpin.tm > 45:
    issues.append('Hairpin Tm > 45C')
if fwd_homodimer.tm > 45 or rev_homodimer.tm > 45:
    issues.append('Homodimer Tm > 45C')
if heterodimer.tm > 45:
    issues.append('Heterodimer Tm > 45C')

if issues:
    print('Issues found:')
    for issue in issues:
        print(f'  - {issue}')
else:
    print('No significant issues detected')
