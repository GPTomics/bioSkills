# tRF and piRNA Profiling - Usage Guide

## Overview

Profile the non-miRNA small RNAs in a small RNA-seq library - tRNA-derived fragments (tRFs/tsRNAs), piRNAs, and rRNA/snoRNA-derived species - using MINTmap, unitas, SPORTS, and proTRAC. The governing reality is that small RNA-seq is a size cut over more than eight biochemically distinct classes, and the library kit decides what is captured: standard ligation needs a 5'-monophosphate and a 3'-OH, so 5'-OH and cyclic-phosphate species (tRNA halves, many tRFs) are silently absent unless a T4 PNK prep was used. Detection is not function: any abundant structured RNA sheds degradation fragments into the size window, so functionality must be earned through end precision, strand bias, phasing, the piRNA ping-pong signature, or AGO/PIWI loading. Multimapping over redundant tRNA and piRNA loci is endemic, which is why tRF tools report exclusive versus ambiguous counts and piRNAs are quantified at cluster level.

## Prerequisites

```bash
conda install -c bioconda mintmap unitas protrac
# SPORTS1.0 and ShortStack (plants) are installed separately from their repositories
pip install numpy pandas
```

## Quick Start

Tell your AI agent:
- "Annotate every small-RNA class in my library and give me the composition"
- "Quantify tRFs with MINTmap and use the exclusive counts"
- "Detect piRNA clusters and test the ping-pong signature"
- "Is this abundant species a real tRF or just tRNA degradation?"
- "My library is TruSeq - can it even see tRNA halves?"

## Example Prompts

### Class Annotation

> "Run unitas to annotate miRNA, piRNA, tRF, rRF, and snoRNA fractions in my library"

> "My TruSeq library shows almost no tRNA halves - is that real or an assay artifact?"

### tRF Quantification

> "Profile tRFs with MINTmap and report exclusive versus ambiguous counts separately"

> "Classify my tRFs as tRF-5, tRF-3, tRF-1, i-tRF, or tRNA-half"

### piRNA Analysis

> "Detect piRNA clusters with proTRAC and check for strand asymmetry and 1U bias"

> "Compute the ping-pong z-score to confirm an active piRNA pathway"

> "I see abundant piRNAs in a plasma sample - is that plausible?"

## What the Agent Will Do

1. Check whether the prep could capture the target class (TruSeq vs T4 PNK / PANDORA-seq) before interpreting absence
2. Annotate all small-RNA classes (unitas/SPORTS) and read the class composition
3. Quantify tRFs at locus resolution with MINTmap, trusting exclusive counts
4. Detect piRNA clusters (proTRAC) and test the ping-pong signature for an active pathway
5. Separate processed functional species from degradation using end precision, strand bias, and phasing

## Tips

- Small RNA-seq is not a miRNA assay - it pools 8+ classes, and the kit sets the composition
- Absence of a 5'-OH/cyclic-phosphate class (tRNA halves, many tRFs) in a TruSeq library is an assay artifact, not low expression - rescue needs T4 PNK prep
- High read count is not evidence of function; require end precision, strand bias, phasing, ping-pong, or AGO/PIWI loading
- tRNA and piRNA loci are redundant - use MINTmap exclusive counts and quantify piRNAs at cluster/family level
- The ping-pong signature is a sharp 10-nt 5'-5' overlap with 1U (primary) and 10A (secondary) - it evidences the active SECONDARY (transposon) pathway specifically
- A flat ping-pong score does not mean no piRNAs: primary piRNAs are phased (1U, Zucchini-dependent), and adult testis is >95% pachytene piRNAs that are non-ping-pong - test phasing too
- Operationalize phasing as a peak in the same-strand 5'-to-5' distance histogram at the modal piRNA length (~26-28 nt); proTRAC reports it natively (it is the primary-pathway complement to the ping-pong test)
- piRNAs carry a 3' 2'-O-methyl that suppresses standard adapter ligation, so low piRNA yield can be an end-chemistry artifact, not low abundance
- Abundant piRNAs in somatic or plasma data are a red flag for misannotation (often tRFs or Y-RNA fragments matching piRBase by chance)
- rRNA-derived fragments are the hardest class - even tiny rRNA decay yields huge counts, so demand reproducibility
- For plants, use ShortStack (24-nt siRNA, phasiRNA biology), not animal miRNA tools

## Related Skills

- smrna-preprocessing - Size windows and end chemistry that determine class capture
- mirdeep2-analysis - tRF/rRF stacks are miRDeep2 false positives
- mirge3-analysis - Known miRNAs and a basic tRF module
- differential-mirna - Count-based DE for tRF/piRNA matrices
- genome-annotation/ncrna-annotation - tRNA/rRNA/snoRNA locus annotation
