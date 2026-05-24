# Library Design - Usage Guide

## Overview

Design pooled sgRNA libraries for CRISPR knockout, interference (CRISPRi), activation (CRISPRa), Cas12a multiplex, base-editor, or prime-editor screens. Covers chemistry selection, on-target scoring (Rule Set 2 / Azimuth / DeepSpCas9 / CRISPRon), off-target scoring (CFD / MIT / Elevation), TSS-relative positioning, control composition, oligo synthesis layout, and post-cloning QC.

## Prerequisites

```bash
pip install crispor biopython pandas numpy
# Optional modern predictors
pip install azimuth-2  # Microsoft Research Rule Set 2 implementation
# CLI tools
conda install -c bioconda crispor flashfry
# For Cas12a / multiplex
pip install crisprDesignData  # Bioconductor data
```

Required inputs: gene list (HGNC symbols or Ensembl IDs), target genome assembly (GRCh38 / GRCm39 / project-specific), screen chemistry (Cas9 / CRISPRi / CRISPRa / Cas12a / BE / PE), and either FANTOM5 CAGE peaks (for CRISPRi/a) or coding-exon coordinates (for Cas9 KO).

## Quick Start

Tell the AI agent what to design:
- "Design a custom Brunello-style Cas9 knockout library for 250 kinases with 4 guides per gene plus 500 non-targeting controls"
- "Design a CRISPRi library against 1,200 lncRNAs using Dolcetto positioning rules against the FANTOM5 TSS"
- "Build a CRISPRa Calabrese-style library targeting -150 to -75 of TSS for 800 transcription factors"
- "Design an enAsCas12a paralog library covering 600 paralog pairs as 4-guide in4mer arrays"
- "Design a base-editor library tiling editing windows across all exons of BRCA1 and BRCA2 for variant scanning"
- "Diagnose why my freshly cloned plasmid pool has Gini 0.28 and ~3% zero-count guides"

## Example Prompts

### Cas9 Knockout Libraries

> "Design a focused Cas9 KO library targeting 350 DNA damage response genes. Use Rule Set 2 / Azimuth on-target scoring and CFD off-target filtering. Target the first 5-65% of each protein, exclude guides with GC outside 30-70% or poly-T runs, 4 guides per gene, 500 non-targeting controls, 50 AAVS1 controls, 50 CEGv2 reference essentials, 50 NEGv1 non-essentials."

> "I'm screening in HCT116 cells. Should I use DeepSpCas9 instead of Azimuth Rule Set 2 for guide ranking? Compare the two and recommend."

> "I have a custom gene list of 80 paralog pairs. Should I use Cas9 single-KO or Cas12a multiplex, and why?"

### CRISPRi / CRISPRa Libraries

> "Design a Dolcetto-style CRISPRi library for 600 transcription factors. Resolve TSS from FANTOM5 highest-rank CAGE peak per gene; if FANTOM5 lacks a peak, fall back to Ensembl canonical TSS and flag those genes for review. Target the Dolcetto window (-50 to +300 from TSS). 6 guides per gene plus 1,000 NTCs."

> "Compare Calabrese vs Horlbeck CRISPRa TSS targeting windows for my activation screen in iPSC-derived neurons. We need to detect modest fold-change activation."

> "My CRISPRi screen has weak dropout signal even on RPL/RPS genes. Audit guide positioning against current FANTOM5 / matched neuron CAGE data and re-design any guides outside ±100 from the empirical TSS."

### Cas12a Multiplex / Paralog Libraries

> "Build an Inzolia-style enAsCas12a 4-guide array library covering 400 paralog pairs in the receptor tyrosine kinase family. Include singleton controls (gene A alone, gene B alone, double-NTC) so we can score genetic interaction = double_KO_LFC - sum(single_KO_LFC)."

> "Design an in4mer triple-KO library testing all triplets within a 25-gene synthetic-lethality hypothesis."

### Base / Prime Editor Libraries

> "Design a CBE saturation library across BRCA1 exons 1-23 tiling every NGG-adjacent spacer that places at least one C in editing positions 4-8. Flag bystander Cs and annotate predicted amino acid change."

> "Build a prime-editor library to install 320 specific ClinVar variants in MLH1, MSH2, MSH6, PMS2. Use PRIDICT2 to score pegRNA candidates and select the top-3 scored pegRNAs per intended edit."

### Library Diagnostics

> "My plasmid pool shows Gini 0.28, skew ratio 6, and 2.5% zero-count guides on Brunello-style sequencing. Diagnose the cause (PCR bias, synthesis defect, cloning bottleneck, library age) and recommend remediation."

> "We see ERBB2 dropping out as 'essential' in HER2-amplified SK-BR-3 cells. Is this a library-design issue or a copy-number bias problem?"

## What the Agent Will Do

1. Confirm chemistry from goal: KO vs knockdown vs activation vs paralog vs variant-function
2. Resolve target coordinates: coding exons for Cas9; FANTOM5 highest-rank CAGE TSS for CRISPRi/a; full transcript for tiling
3. Enumerate PAM-adjacent protospacers in the relevant window per chemistry
4. Score on-target with Rule Set 2 / Azimuth (Cas9), Horlbeck rules (CRISPRi/a), or PRIDICT2 (PE)
5. Score off-target with CFD (genome-wide for Cas9; CRISPOR aggregate for completeness)
6. Filter on GC, poly-T, length, position in CDS / TSS window
7. Select top N guides per gene (4 for Brunello/TKOv3; 6 for Avana/Dolcetto/Calabrese)
8. Append controls: NTCs (~1%), safe-harbor (50-100), CEGv2 reference essentials, NEGv1 non-essentials, olfactory pseudo-targeting
9. Build oligos: subpool primers + BsmBI overhangs + scaffold; check Twist 92K/244K length limits
10. Output library table + oligo synthesis order + expected coverage / cell number / sequencing depth

## Tips

- For new screens in well-characterized cancer lines, use the off-the-shelf Brunello / Dolcetto / Calabrese / Inzolia pool (Addgene) rather than re-designing from scratch -- you inherit decades of community validation including the calibration of MAGeCK / BAGEL2 / Chronos against these specific libraries.
- For tissue-specific lines (neurons, hepatocytes, primary T cells), re-derive TSS from matched CAGE/GRO-seq before CRISPRi/a design; the single most common silent failure of CRISPRi screens is mis-positioned guides against an alternative TSS.
- Plasmid-pool sequencing (200-500 reads/sgRNA, Gini <0.1, ≥99% guide detection) is non-negotiable; everything downstream is normalized against this baseline.
- Subpool synthesis with 10-20k oligos per subpool stretches a single Twist 92K order across multiple focused screens.
- For paralog screens, the difficulty is not the library but the GI scoring: double_LFC - (single_A_LFC + single_B_LFC). Include all three singleton conditions in the library; without them, GI scoring is impossible post hoc.
- Avoid CRISPR-Cas12a libraries when running a screen with Cas9 cell lines or vice versa; the enzymes and PAM are different and orthogonal libraries cannot be mixed.
- For variant-function screens (BE/PE), restrict guides a priori to those predicted >50% editing efficiency; the analytic framework needs editing-efficient guides for variant calls to be confident.

## Library Specification Reference Table

| Library | Modality | Guides/Gene | Coverage Need | Typical Use |
|---------|----------|-------------|---------------|-------------|
| Brunello | Cas9 KO | 4 | 500x | Modern KO standard |
| TKOv3 | Cas9 KO | 4 | 500x | BAGEL2-calibrated |
| Avana | Cas9 KO | 6 | 500x | DepMap legacy (pre-2021) |
| Dolcetto | CRISPRi | 6 | 500x | Knockdown of essential/cuttable-toxicity genes |
| Calabrese | CRISPRa | 6 | 500x | GoF activation |
| Inzolia | Cas12a array | 4-array | 500x array | Paralog buffering / GI screens |

## Control Composition

| Control type | Count in 70k library | Function |
|--------------|----------------------|----------|
| Non-targeting (NTC) | 500-1,000 (~1%) | Primary null distribution |
| Safe-harbor (AAVS1) | 50-100 | Cas9 cut-toxicity baseline |
| Olfactory receptors | 50-100 | Orthogonal non-expressed null |
| CEGv2 reference essentials | 50-100 | Positive control (BAGEL2 calibration) |
| NEGv1 reference non-essentials | 50-100 | Negative control (BAGEL2 calibration) |

## Related Skills

- crispr-screens/screen-qc - Library skew, Gini, replicate correlation, essentialome PR-AUC
- crispr-screens/mageck-analysis - Standard analysis pipeline for the designed library
- crispr-screens/combinatorial-screens - Cas12a multiplex / paralog-pair library design
- crispr-screens/base-editing-analysis - GRACE-style base-editor library design and editing-window analysis
- crispr-screens/prime-editing-screens - PRIDICT2-optimized pegRNA libraries
- crispr-screens/copy-number-correction - Filter amplicon-driven artifacts in cancer-cell-line screens
- crispr-screens/in-vivo-screens - Bottleneck math for animal-model focused libraries
