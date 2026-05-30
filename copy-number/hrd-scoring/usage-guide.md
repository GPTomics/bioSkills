# HRD Scoring Usage Guide

## Overview

Homologous recombination deficiency leaves characteristic copy-number scars. Three are quantified and summed into an HRD score: loss of heterozygosity (LOH segments > 15 Mb), large-scale state transitions (LST), and telomeric allelic imbalance (TAI). A high score predicts response to platinum chemotherapy and PARP inhibitors and underlies approved companion diagnostics. The decisive subtlety is that the score is a *scar*, a record of past HR deficiency, so it can remain high in a tumor that has since restored HR function (for example by a BRCA reversion mutation). This skill covers scarHRD scar computation, the whole-genome HRDetect and CHORD models, and the ploidy/purity caveats.

## Prerequisites

```bash
R -e "remotes::install_github('sztup/scarHRD')"   # GitHub-only
conda install -c bioconda sequenza-utils          # allele-specific input for scarHRD
# HRDetect / CHORD: install from their respective repositories (whole-genome data only)
```

Inputs: allele-specific copy number (a Sequenza `.seqz` file or an ASCAT-style allele-specific segment table) from an adequately pure tumor, with a known ploidy. HRDetect and CHORD additionally need whole-genome somatic mutation calls.

## Quick Start

Tell the AI agent what to do:
- "Compute the HRD score (LOH, LST, TAI) from this Sequenza output with scarHRD"
- "Decide between scar-based HRD and HRDetect for my whole-genome tumor data"
- "Explain why this BRCA-wild-type tumor scores HRD-high"
- "Correct my LST count for whole-genome doubling"
- "Interpret an HRD-high result in a tumor that failed PARP-inhibitor therapy"

## Example Prompts

### Computing the score

> "Run scarHRD on this Sequenza .seqz file in GRCh38, report LOH, LST, TAI and the sum, and confirm the score was computed with the correct ploidy."

> "I have whole-genome sequencing for this tumor. Recommend HRDetect over the scar score and explain what additional evidence it integrates."

### Interpretation

> "This tumor scores HRD-high but is BRCA wild-type. Work through whether this is a true HR lesion elsewhere in the pathway or a false-high from whole-genome doubling."

> "An HRD-high patient progressed on a PARP inhibitor. Explain how a scar score can be correct yet the tumor no longer HR-deficient."

### Method and assay

> "My HRD score came from a targeted panel. Explain why I cannot compare it directly to the Myriad myChoice GIS >= 42 cutoff."

## What the Agent Will Do

1. Confirm the input is allele-specific copy number from an adequately pure sample
2. Supply the tumor ploidy so LST is correctly normalized
3. Compute LOH, LST, and TAI and their sum with scarHRD (or HRDetect/CHORD for WGS)
4. Cross-check an LST-driven high score against whole-genome-doubling status
5. Interpret the score as evidence of past HR deficiency, not current state
6. Use the assay-validated cutoff, not a cross-assay threshold

## Tips

- The HRD score needs allele-specific copy number; total CN or relative log2 cannot give LOH or TAI.
- LST rises with ploidy; always compute the score with the correct ploidy so whole-genome doubling does not inflate it.
- A high score reflects *past* HR deficiency; BRCA reversion restores HR function but the scars persist, so HRD-high is not a guarantee of drug response.
- For whole-genome data, HRDetect and CHORD are more accurate because they add SNV and rearrangement signature evidence.
- Companion-diagnostic cutoffs (GIS >= 42) are validated per assay and are not portable to a different panel or to WGS.
- Allele-specific calling fails below ~30-40% purity; do not score low-purity samples.

## Related Skills

- copy-number/allele-specific-copy-number - Allele-specific copy number input
- copy-number/subclonal-copy-number - Whole-genome-doubling detection for LST correction
- copy-number/recurrent-cnv - Copy-number signatures including the HRD signature
- copy-number/cnv-annotation - Annotating HR-pathway gene copy-number status
- clinical-databases/somatic-signatures - SNV mutational signature 3
- clinical-databases/variant-prioritization - BRCA1/2 and HR-pathway variant interpretation
