# VDJtools Analysis - Usage Guide

## Overview

VDJtools performs standardized post-analysis of TCR/BCR clonotype tables (from MiXCR and other RepSeq pipelines): diversity, clonality, overlap, segment usage, spectratype, and longitudinal tracking. The dominant analytical fact is that a repertoire is a depth-limited sample of an enormous clonal population, so nearly every metric -- richness, Shannon, clonality, Jaccard, shared counts -- is confounded by sequencing depth and by library chemistry. Samples must be brought to a common depth (DownSample, or the resampled CalcDiversityStats table / rarefaction curves) before any cross-sample claim, diversity should be read as a Hill profile rather than one index, overlap should use depth-robust metrics (Morisita-Horn, F2) with a fixed clonotype match key, and public clonotypes reflect generation probability more than antigen selection. immunarch is the actively maintained R alternative that consumes the same tables.

## Prerequisites

```bash
# VDJtools (fat JAR, needs Java 8+)
wget https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.zip
unzip vdjtools-1.2.1.zip
java -version

# One-time: install the R packages the Plot* modules call
java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar RInstall

# Optional modern R alternative
Rscript -e "install.packages('immunarch')"
```

## Quick Start

Tell your AI agent what you want to do:
- "Downsample my repertoires to equal depth, then compare diversity"
- "Report a Hill diversity profile (q=0, q=1, q=2) per sample"
- "Compute Morisita-Horn overlap between conditions on nucleotide clonotypes"
- "Make a spectratype and V-J usage plot for each sample"
- "Track clonotypes across timepoints and flag expansions"

## Example Prompts

### Diversity

> "Bring all samples to a common read depth and compare observed diversity, exp-Shannon and inverse Simpson between responders and non-responders."

> "Draw rarefaction curves so I can see whether the diversity difference survives at equal depth."

> "Summarize clonality per sample, but pair it with inverse Simpson so it is not misleading."

### Overlap and sharing

> "Compute pairwise Morisita-Horn overlap on amino-acid clonotypes and cluster the samples."

> "Find clonotypes shared between pre- and post-treatment, using nucleotide identity with V and J matching."

> "Which public clonotypes are shared across patients, and how should I interpret them given generation probability?"

### Structure and tracking

> "Plot the CDR3-length spectratype and tell me whether it looks polyclonal or skewed."

> "Compare V-gene usage across samples, accounting for possible primer bias."

> "Track the top clonotypes across the vaccination time course and identify expansions."

## What the Agent Will Do

1. Convert upstream clonotype tables to VDJtools format and build the metadata file.
2. Filter to functional clonotypes and, if needed, decontaminate cross-sample chimeras.
3. Downsample all samples to a common depth (or plan to read the resampled diversity table).
4. Run CalcDiversityStats and report the Hill profile from the resampled table; add rarefaction curves.
5. Run CalcPairwiseDistances with a fixed match key, report Morisita-Horn/F2, and cluster.
6. Generate spectratype, segment-usage, and clonal-tracking outputs and summarize.

## Tips

- Compare diversity only on depth-normalized data: use the `.resampled.txt` table or `DownSample` first. Raw comparison across unequal depth is the most common invalidating error.
- Report q=0, q=1 and q=2 together (a Hill profile), never a single index. inverseSimpson (q=2) is the most depth-robust.
- Do not feed chao1/efronThisted frequency-only or rare-clone-filtered data; without genuine singleton/doubleton counts they are unstable, and PCR error inflates them.
- Prefer Morisita-Horn or F2 for overlap; Jaccard and public counts are dominated by the shallower sample.
- Hold the clonotype match key (`-i`) constant across a study and state it in every figure; aa keys inflate sharing via convergent recombination.
- Public clonotypes are mostly high-Pgen convergent sequences, not antigen-specific; condition on Pgen or annotate against VDJdb before biological claims.
- Compare V/J usage only within one library protocol; multiplex-PCR primer bias masquerades as biology.
- Increase Java heap for large datasets: `java -Xmx8g -jar vdjtools.jar ...`.

## Related Skills

- mixcr-analysis - Generate input clonotype tables
- repertoire-visualization - Rarefaction, spectratype and overlap figures
- immcantation-analysis - BCR-aware diversity and clonal analysis
- specificity-annotation - Pgen-aware interpretation of public clonotypes
- experimental-design/sample-size - Sequencing depth and power planning
- workflows/tcr-pipeline - End-to-end orchestration
