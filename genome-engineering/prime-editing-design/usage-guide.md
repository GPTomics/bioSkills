# Prime Editing Design - Usage Guide

## Overview

Designs pegRNAs and nicking guides for prime editing (PE) -- scarless point mutations, small insertions/deletions, and any of the 12 base conversions without a double-strand break or donor. The skill chooses the nick/strand, designs the PBS and RTT as a per-locus *panel* (there is no universal optimum), selects the PE system (PE2/PE3/PE3b/PE4/PE5/ PEmax/PE7), adds the free wins (PAM-disrupting and MMR-evading silent edits, a tevopreQ1 3' motif), ranks with PRIDICT/DeepPrime, and routes large insertions to twinPE/PASTE. Its discipline is that PE efficiency is a cellular-genetics problem (flap resolution + mismatch repair), so the *system* choice and a tested panel matter more than any single oligo.

## Prerequisites

```bash
pip install biopython
# PrimeDesign: Docker-only (no pip) -- docker run ... pinellolab/primedesign primedesign_cli ...
#   Edit is encoded inline: (ref/edit) substitution, (+ins) insertion, (-del) deletion.
# PRIDICT2.0 / DeepPrime: web servers / authors' code for ranking pegRNA panels.
```

## Quick Start

Tell your AI agent what you want to do:
- "Design a pegRNA to correct the G551D mutation in CFTR"
- "Give me a PBS x RTT panel to test for this edit, not just one pegRNA"
- "My pegRNA edits at 2% -- what system or design changes would help?"
- "Should I use PE3 or PE3b here, and how do I avoid indels?"
- "Insert a 2 kb reporter -- is prime editing the right tool?"

## Example Prompts

### Point mutations and small edits

> "Design pegRNAs to install the protective APOE Christchurch variant. Give me a PBS x RTT panel, add a PAM-disrupting silent edit if possible, and tell me which PE system to use in MMR-proficient iPSCs."

> "Correct the sickle-cell mutation (HBB E6V) with prime editing -- and explain whether a base editor would be a better choice for this conversion."

### System and efficiency

> "My PE3 design gives good editing but lots of indels. How do I fix that?"

> "This locus is GC-rich -- how should that change my PBS length?"

> "We're in HCT116 -- do I still benefit from MLH1dn (PE4/PE5)?"

### Large edits

> "I need to knock in a 1.5 kb cassette -- route me to the right method."

## What the Agent Will Do

1. Establish the edit, cell type, MMR status, and delivery (expressed vs synthetic pegRNA).
2. Decide whether a base editor (transition) or integrase method (large insert) is the better tool before committing to PE.
3. Choose the nick/strand and design a PBS x RTT panel, enforcing the 5'-G-prepend and don't-start-on-C rules.
4. Add the free wins: a PAM-disrupting silent edit and 1-2 MMR-evading silent edits.
5. Select the PE system by the axis the problem needs (MMR -> PE4/PE5; indels -> PE3b; stability -> epegRNA/PE7).
6. Rank the panel with PRIDICT2.0/DeepPrime (a prior, not a verdict) and recommend testing.
7. For PE3/PE3b, design the second nicking guide at ~40-100 bp on the appropriate strand.

## Tips

- There is no universal PBS/RTT -- design and test a panel; the optimum is set by local GC/Tm and the locus. emitting a single pegRNA understates the work; deliver a ranked panel.
- The PE system carries the order of magnitude: MLH1dn (PE4/PE5) for MMR-active cells, epegRNA/PE7 for pegRNA stability, PEmax for the protein. PE5max+epegRNA is the modern default.
- PE3 raises indels (its second nick is a transient near-DSB); use PE3b when the edit makes/breaks a protospacer, or stay on PE2 when indels are unacceptable.
- Take the free wins: a PAM-disrupting silent edit stops re-nicking; MMR-evading silent edits raise yield -- both trivial in coding sequence.
- Report edit:indel purity, not efficiency alone; in MMR-deficient lines MLH1dn does nothing.
- Prediction models are priors trained on HEK293T + small edits; weight them less off-distribution and still test. Chromatin can sink a perfect design.
- Use PrimeDesign's exact inline notation; verify it against the repo -- it is the most-hallucinated PE detail.

## Related Skills

- base-editing-design - Preferred for the single transition a base editor can make
- grna-design - Generic spacer scoring; plain-nuclease knockout when precision is unneeded
- off-target-prediction - pegRNA spacer and PE3 nicking-guide off-target considerations
- hdr-template-design - Large-insertion alternative (HDR/HITI)
- crispr-screens/prime-editing-screens - Pooled prime-editing screen analysis
- crispr-screens/crispresso-editing - Quantify intended-edit vs indel rates from amplicons
- variant-calling/variant-annotation - Identify the pathogenic variant to correct or install
