# CLIP Motif Analysis - Usage Guide

## Overview

Discover RBP binding motifs from CLIP-seq peaks or single-nucleotide crosslink sites. The fundamental confound is the uracil bias of UV254 crosslinking - U is the most-crosslinked base, so naive logos centered on CL positions are U-enriched even for non-U-binding RBPs. Modern CL-registered tools (mCross, PEKA) correct for this; legacy peak-based tools (HOMER, MEME, STREME) need GC-matched backgrounds. In vitro RBNS Kd from Dominguez 2018 is the orthogonal validation standard.

## Prerequisites

```bash
conda install -c bioconda homer meme bedtools
# CL-registered tools
git clone https://github.com/chaolinzhanglab/mCross
pip install peka
# RBPamp
pip install rbpamp
```

## Quick Start

Tell your AI agent what you want to do:
- "Run HOMER de novo motif on my CLIPper peaks with a GC-matched 3' UTR background"
- "Register motifs to crosslink positions with mCross"
- "Confirm mCross result with PEKA orthogonal method"
- "Compare in vivo CLIP motif to in vitro RBNS Kd (Dominguez 2018)"
- "Scan known RBP motifs from CISBP-RNA / ATtRACT with FIMO"
- "STREME because MEME is too slow on 5000 peaks"
- "Generate position-specific logo with kpLogo"

## Example Prompts

### De Novo Discovery

> "HOMER findMotifs.pl on my stringent peak set with -rna -len 5,6,7,8 and GC-matched 3' UTR background"

> "STREME on 10000 peak sequences - MEME is too slow"

### CL-Registered Motif

> "Run mCross on PureCLIP single-nt sites; show me the CL-position-offset histogram and the registered PWM"

> "PEKA for positional k-mer enrichment as orthogonal confirmation of mCross"

### Known Motif Validation

> "FIMO scan against CISBP-RNA PWM for PTBP1 at threshold 1e-5"

> "Show me the top RBNS Kd-ranked motif for HNRNPK from Dominguez 2018 and compare to my CLIP top hit"

### Diagnostics

> "Why is my HOMER top motif all AU? Check for background mismatch"

> "Inspect motif information content per position; if mean IC < 0.5 something is wrong"

> "Validate the U-bias correction in my mCross output"

### Reconciliation

> "HOMER says one motif, RBNS says another - reconcile in vivo vs in vitro"

> "Multiple modes in mCross - is this multimodal RBP (SRSF1) or noise?"

## What the Agent Will Do

1. Extract peak sequences with `bedtools getfasta -s` preserving strand
2. Generate GC-matched background from expressed 3' UTRs / introns / CDS (not naive shuffled)
3. Run HOMER `-rna -len 5,6,7,8` and STREME for de novo
4. For mCross / PEKA: provide single-nt CL sites from PureCLIP or CTK CITS
5. FIMO scan against CISBP-RNA / ATtRACT for known motifs (if RBP is published)
6. Compare to RBNS Kd (Dominguez 2018) if available
7. Report QC: information content per position, replicate concordance, CL-offset histogram for mCross
8. Reconcile: 3 independent tools agreeing on core motif = high confidence

## Tips

- **U bias is real and universal.** ~60-80% of CL events at U. Cross-check central column of PWM.
- **HOMER background matters.** Auto-shuffled is GC-blind; provide GC-matched.
- **mCross needs single-nt CL sites.** Peak BED breaks the model.
- **PEKA does not need input control.** Internal normalization vs low-count crosslinks.
- **In vitro RBNS Kd is the gold reference.** Dominguez 2018 = 78 RBPs.
- **Skipper window motifs and CLIPper peak motifs may differ.** Window dilutes; peak sharpens.
- **STREME > MEME for large peak sets.** MEME is O(N^2).
- **Information content > 1.0 bits per position = real motif.** < 0.5 = no signal.
- **FIMO threshold 1e-5 for in vivo data.** Default 1e-4 too lenient.
- **Multimode is biology.** SRSF1 GGA half-sites, PTBP1 multimer; report both.

## Related Skills

- clip-seq/clip-peak-calling - Upstream peaks
- clip-seq/crosslink-site-detection - Single-nt CL sites for mCross / PEKA
- clip-seq/binding-site-annotation - Region context
- clip-seq/clip-deep-learning - RBPNet / DeepRiPe sequence-to-binding
- chip-seq/motif-analysis - DNA-protein analogue
