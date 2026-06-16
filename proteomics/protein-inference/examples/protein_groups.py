'''Parsimony grouping + PICKED protein-group FDR on a labeled peptide->protein map.

Self-contained demo (no input files): a small INLINE peptide-to-protein mapping
including decoy proteins, so the picked-FDR logic runs end to end. With a real
search, replace SAMPLE_MAP with the peptide->protein lists parsed from an idXML
(pyopenms.IdXMLFile) or a MaxQuant proteinGroups.txt, and replace the toy scores
with the inference engine's group probabilities.
'''
# Reference: pandas 2.2+ | Verify API if version differs
import pandas as pd

DECOY_PREFIX = 'DECOY_'  # marks reversed/shuffled decoy proteins
FDR_CUTOFF = 0.01  # 1% protein-group FDR, the community standard for the LIST

# peptide -> proteins it maps to. 'unique' peptides map to one protein, 'shared'
# (degenerate) peptides map to several. Decoy proteins carry the DECOY_ prefix.
SAMPLE_MAP = {
    'AAAEUNIQUER': ['P_A'],
    'BBBAUNIQUEK': ['P_A'],
    'SHAREDABK': ['P_A', 'P_B'],
    'CCCBUNIQUEK': ['P_B'],
    'DDDSUBSETK': ['P_C', 'P_A'],
    'EEEONLYK': ['P_D'],
    'ZZZDECOYK': [DECOY_PREFIX + 'P_A'],   # decoy counterpart of target P_A; picking pairs them
    'YYYDECOYR': [DECOY_PREFIX + 'P_A'],
    'XXXDECOYK': [DECOY_PREFIX + 'P_B'],   # decoy counterpart of target P_B
}


def apply_parsimony(peptide_protein_map):
    protein_peptides = {}
    for pep, prots in peptide_protein_map.items():
        for p in prots:
            protein_peptides.setdefault(p, set()).add(pep)
    all_peptides = set(peptide_protein_map)
    covered = set()
    selected = []
    while covered != all_peptides:
        best = max(protein_peptides, key=lambda p: len(protein_peptides[p] - covered))
        new_coverage = protein_peptides[best] - covered
        if not new_coverage:
            break
        selected.append(best)
        covered |= new_coverage
    return selected, protein_peptides


def build_groups(selected, protein_peptides, peptide_protein_map):
    # collapse proteins with identical observed-peptide evidence into one group
    by_evidence = {}
    for prot in selected:
        by_evidence.setdefault(frozenset(protein_peptides[prot]), []).append(prot)
    groups = []
    for peptides, prots in by_evidence.items():
        n_unique = sum(1 for pep in peptides if len(peptide_protein_map[pep]) == 1)
        # leading protein: most unique peptides first, then accession (canonical proxy)
        prots_sorted = sorted(prots, key=lambda p: (-sum(1 for pep in protein_peptides[p] if len(peptide_protein_map[pep]) == 1), p))
        groups.append({
            'leading_protein': prots_sorted[0],
            'accessions': prots_sorted,
            'n_peptides': len(peptides),
            'n_unique_peptides': n_unique,
            'score': len(peptides) + n_unique,  # toy score; use engine probability with real data
            'is_decoy': all(a.startswith(DECOY_PREFIX) for a in prots_sorted),
        })
    return groups


def picked_group_fdr(groups, cutoff=FDR_CUTOFF, decoy_prefix=DECOY_PREFIX):
    # PICKED FDR: pair each target group with its decoy counterpart (same accessions
    # minus the prefix) and keep only the higher-scoring of the pair before counting.
    by_base = {}
    for g in groups:
        base = frozenset(a.replace(decoy_prefix, '') for a in g['accessions'])
        if base not in by_base or g['score'] > by_base[base]['score']:
            by_base[base] = g
    picked = sorted(by_base.values(), key=lambda g: g['score'], reverse=True)
    targets = decoys = 0
    for g in picked:
        if g['is_decoy']:
            decoys += 1
        else:
            targets += 1
        g['fdr'] = decoys / targets if targets else 1.0
    running_min = 1.0
    for g in reversed(picked):  # enforce monotone q-values from the bottom up
        running_min = min(running_min, g['fdr'])
        g['qvalue'] = running_min
    return picked


selected, protein_peptides = apply_parsimony(SAMPLE_MAP)
print(f'Parsimony kept {len(selected)} proteins out of {len(protein_peptides)} candidates')
print(f'Dropped (subsumable): {sorted(set(protein_peptides) - set(selected))}')

groups = build_groups(selected, protein_peptides, SAMPLE_MAP)
picked = picked_group_fdr(groups)

report = pd.DataFrame(picked).sort_values('score', ascending=False)
report

passing = report[(~report['is_decoy']) & (report['qvalue'] <= FDR_CUTOFF)]
print(f'\nProtein groups passing {FDR_CUTOFF:.0%} picked-group FDR: {len(passing)}')
print(report[['leading_protein', 'n_peptides', 'n_unique_peptides', 'is_decoy', 'qvalue']].to_string(index=False))
