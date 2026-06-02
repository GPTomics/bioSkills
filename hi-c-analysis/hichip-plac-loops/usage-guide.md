# Protein-Directed and Targeted 3C Loop Calling - Usage Guide

## Overview

This skill calls statistically significant chromatin loops from protein-directed and targeted 3C assays - HiChIP, PLAC-seq, Capture Hi-C (PCHi-C), and ChIA-PET - where the contact background is anchored at the assayed protein's binding sites (CTCF, cohesin, H3K27ac, YY1) or at captured baits (promoters), making 1D coverage wildly non-uniform. Generic Hi-C loop callers (cooltools dots, Juicer HiCCUPS) use a local uniform/donut background that is the wrong null for this data and will call a "loop" at every coverage spike. The dedicated callers - FitHiChIP, MAPS, hichipper for HiChIP/PLAC-seq, and CHiCAGO for Capture Hi-C - each regress out the per-anchor coverage bias jointly with the genomic-distance decay before testing significance. The skill also covers the with-vs-without separate ChIP-seq anchor decision, the loose-vs-stringent foreground/background choice, and differential looping between conditions via diffloop.

## Prerequisites

```bash
# FitHiChIP (default; config-driven CLI) - install by git clone or Docker (no bioconda package)
git clone https://github.com/ay-lab/FitHiChIP
# upstream valid pairs come from HiC-Pro
conda install -c bioconda hic-pro
# MAPS (PLAC-seq/HiChIP regression model)
git clone https://github.com/ijuric/MAPS
# hichipper (restriction-aware QC + loops)
pip install hichipper
# CHiCAGO + diffloop (Bioconductor)
# R: BiocManager::install(c('Chicago', 'diffloop'))
```

ChIP-seq peaks for anchors are called separately (see chip-seq/peak-calling). FitHiChIP needs valid pairs in HiC-Pro format, a peak file, and a chrom-sizes file. CHiCAGO needs capture-design files (baitmap, rmap, NPB/NBaitsPB/proxOE) from `makeDesignFiles.py`.

## Quick Start

Tell your AI agent what you want to do:
- "Call loops from my H3K27ac HiChIP with FitHiChIP using my ChIP-seq peaks as anchors"
- "I have PLAC-seq data - run MAPS at 5kb"
- "QC my HiChIP library and call loops with hichipper, then compare conditions with diffloop"
- "This is Promoter-Capture Hi-C - call PIRs with CHiCAGO at score 5"
- "Should I use peak-to-peak or peak-to-all background for CTCF HiChIP?"
- "Don't run cooltools dots on my HiChIP - pick the right caller"

## Example Prompts

### FitHiChIP loop calling
> "I have a cohesin HiChIP library as HiC-Pro valid pairs and an independent CTCF ChIP-seq peak set. Call loops with FitHiChIP at 5kb using the stringent peak-to-peak background and coverage-bias regression, FDR 0.01, 20kb-2Mb distance range, and merge nearby contacts."

### MAPS for PLAC-seq
> "Run MAPS on my H3K4me3 PLAC-seq data at 5kb with a 1Mb binning range, using my MACS2 peaks as anchors and the per-bin genomic-feature/bias file, FDR<=0.01."

### Choosing the background
> "My HiChIP is H3K27ac (broad) versus another that is CTCF (sharp). For each, tell me whether to use FitHiChIP loose/peak-to-all or stringent/peak-to-peak and why, and whether to anchor on ChIP-seq peaks or HiChIP-derived peaks."

### Capture Hi-C with CHiCAGO
> "I have Promoter-Capture Hi-C processed with HiCUP. Build the CHiCAGO experiment from my capture design files, run the pipeline to fit the Brownian+technical background, and export significant promoter-interacting regions at score >= 5 in WashU format."

### Differential loops
> "I have hichipper loop calls for two WT and two KO HiChIP replicates. Build a union loop set, drop sub-20kb loops, and test for differential looping between WT and KO with diffloop."

### Avoiding the wrong tool
> "Someone ran cooltools dots on this HiChIP and got a loop at every peak. Explain why that null is wrong for protein-directed data and re-call the loops correctly."

## What the Agent Will Do

1. Identify the assay (HiChIP/PLAC-seq vs Capture Hi-C/PCHi-C) and route to the correct caller family - never generic Hi-C dots/HiCCUPS.
2. Establish anchors: prefer an independent ChIP-seq peak set; fall back to HiChIP-derived peaks (hichipper/HiChIP-Peaks) only when no ChIP exists, flagging the circularity.
3. For HiChIP/PLAC-seq, write a FitHiChIP config (resolution, distance range, foreground IntType, background UseP2PBackgrnd, bias type, FDR) or configure MAPS; choose loose/peak-to-all for broad marks and stringent/peak-to-peak for sharp factors.
4. For Capture Hi-C, build the CHiCAGO experiment from the capture design files and run the two-component-background pipeline, exporting PIRs at score >= 5.
5. Set a lower distance floor (~20kb) to exclude the self-ligation/diagonal regime.
6. QC the library (hichipper restriction-aware metrics; replicate reproducibility of called loops).
7. For condition comparisons, build a union loop set and run an edgeR-style count test (diffloop / FitHiChIP DiffAnalysis), not a per-pixel matrix subtraction.

## Tips

- Generic loop callers (cooltools dots, Juicer HiCCUPS) are the WRONG null for HiChIP/PLAC-seq/Capture data - their uniform background reads coverage spikes as loops. Use FitHiChIP/MAPS/CHiCAGO.
- The hard, unique work is the coverage+distance significance model and lives here; peak calling lives in chip-seq/peak-calling. This skill consumes peaks and produces FDR-controlled loops.
- Prefer an independent ChIP-seq peak set as anchors. Using peaks derived from the same HiChIP reads couples the peak and loop error structures.
- Match the foreground/background to the mark: loose/peak-to-all for broad marks (H3K27ac), stringent/peak-to-peak for sharp factors (CTCF/cohesin).
- FitHiChIP is config-driven: the caller, background, bias model, and FDR are set in the key=value file passed with `-C`, not on the command line.
- Capture Hi-C is asymmetric (bait x other-end) - a symmetric caller cannot apply; CHiCAGO's per-bait two-component background is mandatory.
- Set a lower distance threshold (~20kb) or the diagonal/self-ligation regime floods the call set.
- Harmonize chromosome naming (`chr1` vs `1`) across valid pairs, peak file, and chrom sizes; mismatches silently zero out the anchors.
- Differential looping has no DESeq2 analog: union the loops, count per replicate, test the delta (diffloop).

## Related Skills

- loop-calling - The bulk in-situ Hi-C counterpart (cooltools dots / chromosight); correct null there, wrong null here
- hic-differential - Matrix-level condition comparison framing behind differential loops
- contact-pairs - Produces the valid pairs (HiC-Pro / pairtools) these callers consume
- hic-data-io - Cooler handling of the contact maps upstream of anchored loop calling
- chip-seq/peak-calling - Calls the ChIP-seq peaks used as independent loop anchors
- chip-seq/peak-annotation - Annotate loop anchors with TFs/genes
- atac-seq/enhancer-gene-linking - Enhancer-promoter contacts complement HiChIP/PCHi-C loops
- genome-intervals/overlap-significance - Test loop-anchor enrichment at features against a structured null
