# Primer Validation - Usage Guide

## Overview

This skill validates chosen PCR/qPCR oligos for intramolecular thermodynamic liabilities with primer3-py: hairpins, self-dimers, cross-dimers, and 3'-end stability. The central idea it enforces: a "dimer-free" verdict is a thermodynamic prediction at the supplied salt, Mg2+, dNTP, and oligo conditions and the temperature of evaluation, and the 3' end is the lethal locus because a 3'-anchored dimer is polymerase-extendable into primer-dimer. So structures are ranked by their dG at the annealing temperature and by 3'-end involvement, not by a single global score. Use it to check primer pairs before ordering, troubleshoot dimers or smears, or screen many oligos; route genome-wide off-target checking to primer-specificity and primer design to primer-basics.

## Prerequisites

```bash
pip install primer3-py
```

- Input is one or more oligo sequences (primers and/or a probe), as strings.
- Supply the real reaction conditions: monovalent salt, Mg2+, dNTP, and oligo concentration, plus temp_c set to the annealing temperature. The same primer reads "fine" or "dimer-prone" depending on these.
- ThermoResult dG is in cal/mol (divide by 1000 for kcal/mol), and `.structure_found` must be checked before trusting any number.
- This skill checks the oligos against themselves and each other; it does NOT check whether they bind off-target in the genome (that is primer-specificity).

## Quick Start

Tell your AI agent what you want to do:
- "Check this primer pair for hairpins and dimers at my qPCR conditions (50 mM monovalent, 3 mM Mg, anneal 60 C)"
- "Is the 3' end of my forward and reverse primer forming a stable dimer?"
- "Screen these 12 oligos and flag any with strong secondary structure"
- "My PCR gives a low-molecular-weight band; could it be a primer-dimer?"

## Example Prompts

### Pre-order pair check
> "Validate this forward/reverse pair before I order them. Use my reaction conditions (50 mM K, 3 mM Mg, 0.8 mM dNTP, 250 nM primer) and evaluate at the 60 C anneal step. Tell me the hairpin, homodimer, and heterodimer dG and whether the 3' ends pair."

### Troubleshoot a primer-dimer
> "I see a strong low-Tm peak in my SYBR melt curve and a small band on the gel. Check whether my primers form a 3'-end heterodimer that could be extending into primer-dimer, and show the structure."

### Batch screening
> "Here are 20 candidate primers. Compute hairpin and homodimer Tm for each with the fast Tm-only functions and give me the ones above 45 C to inspect in detail."

### Tailed primers
> "My primers carry EcoRI and BamHI 5' tails. Validate the full tailed oligos for dimers, since the restriction sites are palindromic."

## What the Agent Will Do

1. Read the oligo sequence(s) and the reaction conditions (salt, Mg2+, dNTP, oligo concentration, annealing temperature).
2. Run calc_hairpin and calc_homodimer on each oligo and calc_heterodimer on the pair at those conditions and temp_c.
3. Compute calc_end_stability to expose 3'-end-anchored, extendable dimers and print the ASCII structure so the 3'-end pairing is visible.
4. Gate every result on `.structure_found` and convert dG from cal/mol to kcal/mol before comparing to thresholds.
5. Compare the two primer Tms and flag a mismatch greater than ~2 C.
6. Report flags as conditions-dependent (not absolute pass/fail) and route a redesign to primer-basics or an off-target check to primer-specificity.

## Tips

- Evaluate at the annealing temperature (set temp_c), not the 37 C default; a hairpin that melts well below Ta is harmless.
- Weight 3'-end structures most heavily; a weak dimer that pairs the 3' ends beats a strong one with free 3' ends for causing artifact.
- Convert dG to kcal/mol (divide cal/mol by 1000) before applying any threshold.
- Check the heterodimer explicitly; two individually clean primers can still cross-dimer.
- For a tailed primer, validate the full oligo, not just the binding core; palindromic sites and Gibson arms dimerize.
- Treat thresholds as flags for inspection, not verdicts; read the ASCII structure and judge at the real conditions.

## Related Skills

- primer-basics - Design Tm-matched primer pairs (redesign if validation fails)
- primer-specificity - Genome-wide off-target / in-silico PCR (a different question)
- qpcr-primers - Co-design qPCR primers and probes, including probe self-structure
- sequence-manipulation/seq-objects - Reverse-complement and assemble tailed oligos to validate
