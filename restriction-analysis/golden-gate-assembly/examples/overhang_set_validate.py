'''Validate a Golden Gate fusion-overhang set'''
# Reference: biopython 1.83+ (API verified on 1.86) | Verify API if version differs
# Design rules: every junction overhang distinct; none palindromic (self-ligates); no overhang
# equal to the reverse complement of another (cross-ligates); avoid homopolymers (low fidelity,
# Potapov 2018). For many fragments, use curated high-fidelity overhang sets.

from Bio.Seq import Seq


def validate_overhang_set(overhangs):
    issues = []
    if len(set(overhangs)) != len(overhangs):
        issues.append('duplicate overhangs -> parts assemble in more than one order')

    for o in overhangs:
        if o == str(Seq(o).reverse_complement()):
            issues.append(f'{o}: palindromic -> self-ligates')
        if len(set(o)) == 1:
            issues.append(f'{o}: homopolymer -> low ligation fidelity')

    rc = {o: str(Seq(o).reverse_complement()) for o in overhangs}
    for a in overhangs:
        for b in overhangs:
            if a != b and rc[a] == b:
                issues.append(f'{a} is the reverse complement of {b} -> cross-ligates')

    return issues or ['overhang set OK']


# A clean four-junction set vs a deliberately broken one
for label, overhangs in [
    ('clean set', ['AATG', 'GCTT', 'TACT', 'AGGA']),
    ('broken set', ['AATG', 'AATG', 'GATC', 'GGGG']),   # duplicate, palindrome, homopolymer
]:
    print(f'{label}: {overhangs}')
    for line in validate_overhang_set(overhangs):
        print(f'  {line}')
    print()
