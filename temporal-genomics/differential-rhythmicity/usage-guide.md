# Differential Rhythmicity

## Overview

Differential rhythmicity asks how an oscillation CHANGES between conditions, genotypes, treatments, tissues, or ages - not whether a feature is rhythmic in one condition (that is single-condition detection, the sibling temporal-genomics/circadian-rhythms skill). It is a statistically distinct question, and its dominant failure mode is the detect-then-Venn anti-pattern: subtracting two independently thresholded rhythmicity lists systematically overestimates reprogramming because near p~0.05 the lists differ mostly from threshold noise, not biology. The correct approach fits a single model spanning both conditions and tests the condition x time interaction directly, borrowing strength across conditions. A second distinction that matters is differential EXPRESSION (a change in mean level, the condition main effect) versus differential RHYTHMICITY (a change in the oscillation, the condition x time interaction) - conflating them mislabels mean shifts as rhythm reprogramming. Because a difference test estimates an interaction, differential rhythmicity needs a matched, evenly-sampled, replicated grid in every condition and more power than detection.

## Prerequisites

- R 4.2+
- Install: `install.packages(c('limorhyde', 'circacompare'))`
- DODR is archived on CRAN (removed 2020); install from the archive with `remotes::install_version('DODR')`, or rely on it being pulled in as a compareRhythms dependency
- Bioconductor: `BiocManager::install(c('limma', 'edgeR', 'DESeq2', 'rain'))`
- GitHub tools: `remotes::install_github(c('naef-lab/dryR', 'bharathananth/compareRhythms', 'diffCircadian/diffCircadian'))`
- Single-condition rhythm detection should already be run (see temporal-genomics/circadian-rhythms) so the amplitude filter and reference-group rhythmicity are known.

## Quick Start

Tell your AI agent what you want to do:
- "Test which genes change their 24h rhythm between wild-type and knockout"
- "Classify each gene as gain, loss, phase change, or amplitude change relative to control"
- "Separate genes whose mean level shifts from genes whose oscillation changes"
- "Compare rhythms between young and aged tissue without intersecting two rhythmicity lists"
- "Estimate the amplitude and phase difference for these specific clock genes with p-values"

## Example Prompts

### Genome-scale differential rhythmicity
> "I have RNA-seq time-courses for wild-type and clock-mutant liver, sampled every 4h for 48h with 2 replicates. Fit a model that tests the genotype x time interaction so I get differential rhythmicity separately from differential expression, and rank the genes that changed their rhythm."

### Direct outcome-class classification
> "Classify every gene as gain-of-rhythm, loss-of-rhythm, phase change, amplitude change, or unchanged relative to the normal-chow reference, using a method built to replace the Venn-diagram approach."

### More than two conditions
> "I have three feeding regimes sampled across the day. Assign each gene a parsimonious model of which conditions share amplitude and phase, rather than running pairwise tests."

### Targeted confirmation
> "For these eight candidate genes, give me explicit estimates and p-values for the difference in mesor, amplitude, and phase between treatment and control."

### Guarding a claim
> "My knockout shows lower bulk amplitude - help me decide whether that is loss of the cell-autonomous rhythm or loss of synchrony between cells before I call it arrhythmic."

## What the Agent Will Do

1. Confirm the design: matched, evenly-spaced timepoints across conditions, >=2 cycles, >=6 samples/cycle, and replicates in each condition; flag unmatched grids that break the interaction model.
2. Choose a method: a genome-scale screen (LimoRhyde for two conditions integrated with a DE pipeline; dryR for >=2 conditions; compareRhythms for direct gain/loss/change/same) plus optional targeted confirmation (CircaCompare, diffCircadian) on hits.
3. Fit one model across conditions and test the condition x time interaction (differential rhythmicity), keeping the condition main effect (differential expression) as a separate result.
4. Amplitude-filter, then classify significant features into the four outcome classes by comparing per-condition amplitude and phase.
5. Guard the interpretation: report reduced ensemble amplitude rather than "arrhythmic" when synchrony is a candidate cause, and distinguish driven (LD/masking) changes from endogenous clock changes.

## Tips

- Never define "lost rhythm" as one rhythmicity list minus another; use an interaction or model-selection test that spans both conditions.
- The interaction F-test indicates THAT a rhythm changed, not which class; read per-condition amplitude and phase to assign gain/loss/phase/amplitude.
- A condition main-effect hit is differential expression (mean shift), not differential rhythmicity; report the two separately.
- Match the sampling grid across conditions; unmatched or unevenly-spaced timepoints make the condition x time terms non-estimable or confounded.
- Differential rhythmicity needs more power than detection; if a design was barely adequate to detect a rhythm, it is too thin to confidently call a rhythm change.
- Reduced bulk amplitude can be desynchrony rather than loss of per-cell rhythm; report "reduced ensemble amplitude" and separate causes with single-cell or imaging data.
- Under a light-dark cycle a rhythm change can be a driven/masking change; endogenous rewiring claims need free-running (constant-darkness) data.
- compareRhythms `deseq2`/`edger`/`voom` expect raw counts; use `mod_sel`/`limma`/`cosinor` for normalized or microarray data.

## Related Skills

temporal-genomics/circadian-rhythms - Single-condition rhythm detection and estimation; run first, then compare across conditions here
differential-expression/timeseries-de - Temporal differential expression (a trend over time), which is differential expression, not differential rhythmicity
temporal-genomics/temporal-clustering - Group differentially-rhythmic genes by the shape of their change
single-cell/preprocessing - Entry to cell-level analysis, to separate reduced ensemble amplitude (desynchrony) from true loss of per-cell rhythm
