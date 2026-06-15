'''Validate a primer pair for hairpins and dimers at reaction conditions, weighting the 3' end.
A dimer that ties up the recessed 3' end is polymerase-extendable into primer-dimer; rank by dG at
the annealing temperature and by 3'-end involvement, not by global dG/Tm. Self-contained, no file IO.'''
# Reference: primer3-py 2.3+ | Verify API if version differs

import primer3

fwd = 'GTCTCCTCTGACTTCAACAGCG'
rev = 'ACCACCCTGTTGCTGTAGCCAA'   # this demo pair is deliberately ~3.7 C Tm-mismatched, so the pair-Tm flag fires

COND = dict(mv_conc=50.0, dv_conc=3.0, dntp_conc=0.8, dna_conc=250.0, temp_c=60.0)  # reaction salt/Mg/dNTP/oligo + Ta
TM_COND = {k: COND[k] for k in ('mv_conc', 'dv_conc', 'dntp_conc', 'dna_conc')}      # calc_tm takes no temp_c
HAIRPIN_MARGIN = 10      # flag a hairpin whose Tm is within this many C of Ta; it is not denatured at anneal
DIMER_DG_KCAL = -9.0     # flag whole-molecule dimers below this (kcal/mol); heuristic, condition-dependent
END_DG_KCAL = -5.0       # stricter line for 3'-end-anchored (extendable) dimers
TA = COND['temp_c']

def report(label, res):
    if not res.structure_found:
        print(f'{label}: no structure')
        return None
    print(f'{label}: Tm={res.tm:.1f}C dG={res.dg / 1000:.2f} kcal/mol')   # .dg is cal/mol
    return res

print('=== per-primer self-structure (at reaction conditions, temp_c = Ta) ===')
warnings = []
for name, seq in [('fwd', fwd), ('rev', rev)]:
    hp = report(f'{name} hairpin', primer3.calc_hairpin(seq, **COND))
    if hp and hp.tm > TA - HAIRPIN_MARGIN:
        warnings.append(f'{name} hairpin Tm within {HAIRPIN_MARGIN}C of Ta')
    hd = report(f'{name} homodimer', primer3.calc_homodimer(seq, **COND))
    if hd and hd.dg / 1000 < DIMER_DG_KCAL:
        warnings.append(f'{name} homodimer dG < {DIMER_DG_KCAL} kcal/mol')

print('\n=== cross-dimer and 3'"'"'-end (the extendable, lethal class) ===')
het = report('heterodimer', primer3.calc_heterodimer(fwd, rev, **COND))
if het and het.dg / 1000 < DIMER_DG_KCAL:
    warnings.append('heterodimer dG below threshold')
for a, b, who in [(fwd, rev, 'fwd'), (rev, fwd, 'rev')]:   # check BOTH 3' ends: calc_end_stability scores arg1's 3' end
    end = primer3.calc_end_stability(a, b, **COND)
    print(f"{who} 3'-end stability dG={end.dg / 1000:.2f} kcal/mol (extendable-dimer risk)")
    if end.dg / 1000 < END_DG_KCAL:
        warnings.append(f"{who} 3'-end dimer is stable -> extendable into primer-dimer")

dtm = abs(primer3.calc_tm(fwd, **TM_COND) - primer3.calc_tm(rev, **TM_COND))
print(f'\npair Tm difference={dtm:.1f}C')
if dtm > 2:
    warnings.append('pair Tm difference > 2C (one strand will dominate)')

print('\n=== verdict ===')
print('\n'.join(f'  - {w}' for w in warnings) if warnings else '  no flags at these conditions')
