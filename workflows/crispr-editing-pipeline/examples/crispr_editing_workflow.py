#!/usr/bin/env python3
'''CRISPR editing pipeline: orchestrate guide -> off-target -> modality -> donor -> validation.

This orchestrator does the deterministic, local parts honestly (PAM enumeration, hard filters,
a Bae-style out-of-frame estimate, a codon-checked HDR block) and ROUTES the parts that need
trained models or a genome to the real tools: on-target ranking to CRISPOR (context-valid
model), off-target nomination to Cas-OFFinder/CRISPRme, base-editor outcomes to BE-Hive,
prime-editing ranking to PRIDICT/DeepPrime. It does not fabricate activity or specificity scores.
'''
# Reference: biopython 1.83+, pandas 2.2+, matplotlib 3.8+ | Verify API if version differs
import re
from math import exp
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Data import CodonTable

GC_MIN, GC_MAX = 0.40, 0.70    # outside this on-target activity falls off (Doench 2014); hard-dropped here
LENGTH_WEIGHT = 20.0           # MMEJ deletion-length decay (Bae 2014 framework)
EDIT_TO_CUT_MAX = 10           # HDR incorporation falls sharply beyond ~10 bp (Paquet 2016)
FWD = CodonTable.standard_dna_table.forward_table
OUTPUT_DIR = Path('crispr_design_results')


def find_guides(seq):
    '''Enumerate SpCas9 (NGG) spacers on both strands; cut ~3 bp 5' of the PAM (forward coords).'''
    seq, n, rows = seq.upper(), len(seq), []
    for m in re.finditer(r'(?=([ACGT]GG))', seq):
        p = m.start()
        if p >= 20:
            rows.append({'spacer': seq[p - 20:p], 'pam': seq[p:p + 3], 'cut': p - 3, 'strand': '+'})
    rc = str(Seq(seq).reverse_complement())
    for m in re.finditer(r'(?=([ACGT]GG))', rc):
        p = m.start()
        if p >= 20:
            rows.append({'spacer': rc[p - 20:p], 'pam': rc[p:p + 3], 'cut': n - (p - 3), 'strand': '-'})
    return pd.DataFrame(rows)


def u6_valid(spacer):
    gc = sum(c in 'GC' for c in spacer) / len(spacer)
    return 'TTTT' not in spacer and GC_MIN <= gc <= GC_MAX   # TTTT terminates Pol III


def out_of_frame(local, cut):
    '''Bae-style microhomology out-of-frame estimate (frameshift fraction; >66 preferred for KO).'''
    left, right, wo, wt, seen = local[:cut], local[cut:], 0.0, 0.0, set()
    for k in range(min(len(left), len(right)), 2, -1):
        for i in range(len(left) - k + 1):
            mh = left[i:i + k]
            for j in range(len(right) - k + 1):
                if right[j:j + k] == mh:
                    dl, key = (cut + j) - i, (i, cut + j)
                    if dl > 0 and key not in seen:
                        seen.add(key)
                        w = exp(-dl / LENGTH_WEIGHT) * ((k - sum(c in 'GC' for c in mh)) + 2 * sum(c in 'GC' for c in mh))
                        wt += w
                        wo += w if dl % 3 else 0
    return 100.0 * wo / wt if wt else 0.0


def codon_checked_pam_block(donor, pam_start, cds_offset):
    '''Disrupt the NGG PAM with a synonymous change so the edited allele cannot be re-cut.'''
    for g in (pam_start + 2, pam_start + 1):
        if donor[g] != 'G':
            continue
        ci = cds_offset + 3 * ((g - cds_offset) // 3)
        codon = donor[ci:ci + 3]
        if len(codon) < 3:
            continue
        for base in 'ACT':
            cand = codon[:g - ci] + base + codon[g - ci + 1:]
            if FWD.get(cand) == FWD.get(codon) and FWD.get(codon):
                return donor[:ci] + cand + donor[ci + 3:], f'silent PAM block {codon}->{cand}'
    return donor, 'no synonymous PAM block -- use silent SEED mutations instead'


def plot_outcome_landscape(guides, target_len, path):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.scatter(guides['cut'], guides['out_of_frame'], s=40, alpha=0.7)
    ax.axhline(66, ls='--', alpha=0.6, label='out-of-frame 66 (frameshift-reliable)')
    ax.set_xlim(0, target_len); ax.set_ylim(0, 100)
    ax.set_xlabel('cut position (bp)'); ax.set_ylabel('predicted out-of-frame %')
    ax.legend(loc='lower right')
    fig.tight_layout(); fig.savefig(path, dpi=150); print(f'saved {path}')


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    target = ('ATGGCCGAGACCAAGTTCGACGAGCTGAAGGCCTACATCGAGAAGCTGGGCGAGAAGGGCTTCGAC'
              'GAGGTGAAGCTGTACGGCGACGAGCTGAAGTTCTACGGCGAGAAGCTGTTCGAGAAGCTGTTCTAA')

    # Stage 1 -- guide design (rank by frameshift fraction; route on-target to CRISPOR)
    guides = find_guides(target)
    guides = guides[guides['spacer'].map(u6_valid)].copy()
    guides['out_of_frame'] = [out_of_frame(target[max(0, c - 30):c + 30], c - max(0, c - 30)) for c in guides['cut']]
    guides = guides.sort_values('out_of_frame', ascending=False)
    guides.to_csv(OUTPUT_DIR / 'guides_ranked.tsv', sep='\t', index=False)
    print(f'Stage 1: {len(guides)} filtered guides; rank on-target with CRISPOR, carry 3-6.')

    # Stage 2 -- off-target (route to Cas-OFFinder / CRISPRme; in-silico nomination != verdict)
    print('Stage 2: run Cas-OFFinder (+bulges) / CRISPRme; reject low-mm high-CFD off-targets in genes.')

    # Stage 3a -- knockout: take the most frameshift-prone guide
    ko = guides.iloc[0]
    print(f"Stage 3a (KO): {ko['spacer']} PAM={ko['pam']} cut={ko['cut']} out-of-frame~{ko['out_of_frame']:.0f}")

    # Stage 3d -- HDR knock-in with a mandatory codon-checked blocking mutation
    pam_start = next(i for i in range(20, len(target) - 3) if target[i + 1:i + 3] == 'GG'
                     and codon_checked_pam_block(target, i, 0)[0] != target)
    blocked, note = codon_checked_pam_block(target, pam_start, 0)
    print(f'Stage 3d (HDR): {note}; edit must sit within {EDIT_TO_CUT_MAX} bp of the cut; '
          'pick format by cell type (ssODN/lssDNA/AAV; HITI for post-mitotic).')

    # Stage 4 -- validation
    print('Stage 4: quantify intended edit + indels by amplicon deep-seq (CRISPResso2); report purity + LoD.')

    plot_outcome_landscape(guides, len(target), OUTPUT_DIR / 'outcome_landscape.pdf')
    print(f'Done. Results in {OUTPUT_DIR}/')


if __name__ == '__main__':
    main()
