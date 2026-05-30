# Focal Amplification and ecDNA Usage Guide

## Overview

A depth caller reports that an oncogene is amplified and at what amplitude, but not the structure. Focal amplifications come in distinct architectures: extrachromosomal DNA (ecDNA, circular and episomal), breakage-fusion-bridge cycles, homogeneously staining regions, and simple linear amplification. The architecture determines the biology: ecDNA lacks a centromere, segregates unequally, and can surge in copy number under selection, making it a structural basis for therapy resistance. Resolving architecture requires the breakpoint graph, which is what AmpliconArchitect and the AmpliconSuite pipeline reconstruct. This skill covers seed selection, graph reconstruction, and ecDNA classification.

## Prerequisites

```bash
# AmpliconSuite-pipeline bundles AmpliconArchitect and AmpliconClassifier
conda install -c bioconda ampliconsuite          # or install from GitHub
conda install -c bioconda cnvkit samtools        # CNV seeding
# Configure $AA_DATA_REPO (reference download) and a free academic Mosek license
```

Inputs: a whole-genome sequencing BAM (AmpliconArchitect is designed for WGS, not panels/WES); copy-number seeds (high-CN focal regions, ideally from a vetted CNVkit callset); a genome build matching the BAM.

## Quick Start

Tell the AI agent what to do:
- "Reconstruct the architecture of this focal MYC amplification from WGS"
- "Determine whether this high-copy EGFR amplification is ecDNA or a chromosomal HSR"
- "Run the AmpliconSuite pipeline end-to-end on this tumor BAM"
- "Select copy-number seeds for AmpliconArchitect from my CNVkit calls"
- "Explain why a depth caller cannot tell me whether an amplification is ecDNA"

## Example Prompts

### Architecture reconstruction

> "Run AmpliconSuite end-to-end on this tumor WGS BAM in GRCh38 and report which amplicons AmpliconClassifier labels as ecDNA versus BFB versus HSR."

> "Seed AmpliconArchitect from my vetted CNVkit focal calls and reconstruct the breakpoint graph for the chr8 amplicon."

### Interpretation

> "This EGFR amplicon has copy number around 40 and spans three non-contiguous segments. Assess whether this pattern is consistent with ecDNA and what would confirm it."

> "AmpliconClassifier returned 'unknown' for a clearly amplified locus. Explain the short-read resolution limit and what orthogonal assay to use."

### Scope

> "My depth caller says this region is amplified. Explain why that does not establish ecDNA and what breakpoint-graph evidence is required."

## What the Agent Will Do

1. Confirm the input is adequate-coverage whole-genome sequencing
2. Build or vet copy-number seeds: focal, high-CN regions only
3. Run AmpliconArchitect to reconstruct the breakpoint graph and optimize balanced flow
4. Classify each amplicon with AmpliconClassifier (ecDNA / BFB / HSR / linear)
5. Cross-check ecDNA calls against the closed-cycle graph signature
6. Recommend orthogonal confirmation (FISH, single-cell, optical mapping) for high-stakes calls

## Tips

- Depth establishes that a region is amplified and how much; it never establishes the architecture, since ecDNA and a chromosomal HSR look identical in a depth profile.
- AmpliconArchitect reconstructs the regions it is seeded with; garbage seeds (noisy or arm-level) produce garbage amplicons.
- An ecDNA call requires a closed-cycle breakpoint graph plus an AmpliconClassifier label, ideally with orthogonal confirmation.
- AmpliconArchitect needs whole-genome sequencing and adequate coverage; it is not for panels or WES.
- Set the genome build explicitly; AmpliconArchitect was historically hg19-centric.
- Complex amplicons exceed short-read resolution; escalate to optical mapping or long reads when the graph is fragmented.

## Related Skills

- copy-number/cnvkit-analysis - Generates copy-number seeds for amplicon reconstruction
- copy-number/recurrent-cnv - Cohort-level recurrent focal amplification
- copy-number/allele-specific-copy-number - Absolute copy number of amplified loci
- copy-number/cnv-annotation - Oncogene annotation of amplified regions
- copy-number/subclonal-copy-number - Subclonal dynamics of ecDNA copy number
- long-read-sequencing/structural-variants - Long-read resolution of complex amplicons
