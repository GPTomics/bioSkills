# Differential m6A Analysis - Usage Guide

## Overview

Identify differential m6A methylation between conditions from MeRIP-seq using exomePeak2 differential mode, QNB beta-binomial, RADAR, MeTDiff, and (in the defensible paired-symmetric case) edgeR / DESeq2 on featureCounts-on-peaks matrices. Covers paired / interaction / batch designs, the stoichiometry-vs-expression-vs-IP-efficiency confound that all MeRIP differential methods inherit, normalisation choice, the McIntyre 2020 reproducibility caveat, effect-size filtering as a guardrail, and orthogonal-validation routes for absolute stoichiometry claims.

## Prerequisites

```r
BiocManager::install(c('exomePeak2', 'GenomicFeatures', 'BSgenome.Hsapiens.UCSC.hg38',
                       'DESeq2', 'edgeR', 'Rsubread', 'rtracklayer', 'GenomicAlignments',
                       'GenomicRanges'))
devtools::install_github('lzcyzm/QNB')
devtools::install_github('scottzijiezhang/RADAR')
devtools::install_github('compgenomics/MeTPeak')
```

Reference inputs:

- Coordinate-sorted indexed GENOME BAM files for IP and input per replicate per condition (from merip-preprocessing)
- Matching GENCODE / Ensembl GTF
- Sample sheet: sample, condition, biological replicate, antibody clone, antibody lot, batch
- (Optional) Pre-called peak BED for QNB / RADAR (often generated in m6a-peak-calling)

## Quick Start

- "Run exomePeak2 differential mode with 3 control and 3 treatment paired IP/input BAMs"
- "Run QNB beta-binomial differential on pre-called peaks for a small-N (2 vs 2) design"
- "Build a volcano plot of differential m6A peaks with effect-size and FDR thresholds"
- "Include antibody lot as a fixed effect in the differential design matrix"
- "Cross-validate differential calls against published GLORI sites"
- "Distinguish differential m6A from expression change at peak-bearing transcripts"

## Example Prompts

### Standard Differential

> "Run exomePeak2 differential mode on 3 control IP/input pairs and 3 treatment IP/input pairs against GRCh38 + GENCODE v44 GTF; use BSgenome for GC correction; export differential peaks with log2FC and FDR."

> "Run QNB beta-binomial differential on pre-called exomePeak2 peaks with featureCounts-derived count matrices; condition factor is ctrl vs treat."

### Complex Designs

> "Build a paired design with patient as blocking factor; 4 patients each contributing pre-treatment and post-treatment IP/input pairs."

> "Build an interaction design: genotype × treatment with 3 replicates per cell."

> "Include `antibody_lot` as a fixed effect in the design because the study used two lots."

### Validation and Effect-Size Filtering

> "Apply |log2FC| >= 0.5 and FDR < 0.05 filters; report effect-size distribution; cross-check direction concordance across replicates."

> "Cross-validate the top 20 differential peaks against published GLORI sites; flag any whose direction disagrees."

### Visualisation

> "Build a volcano plot of differential m6A peaks with effect-size and FDR thresholds highlighted; annotate top 10 by gene name."

> "Build an MA plot showing the relationship between mean IP/input ratio and log2 fold-change to diagnose normalisation."

> "Per-peak boxplot of IP/input ratio across replicates for the top 5 differential peaks."

### Distinguishing Stoichiometry from Expression

> "For each differential m6A peak, compute the input log2FC alongside the IP/input log2FC; flag peaks where input changes track IP changes (likely expression-driven)."

### RADAR / MeTDiff Sensitivity

> "Run RADAR diffIP() on the same data as a sensitivity analysis; report concordance with exomePeak2 differential calls."

> "Run MeTDiff (paired with MeTPeak peak calling); report 2-tool consensus differential set."

## What the Agent Will Do

1. Verify BAM and GTF chromosome naming consistency (chr1 vs 1)
2. Build TxDb from matched GTF
3. Inspect sample-sheet metadata for batch / antibody-lot confounding with condition
4. Build `experiment_design` data.frame including condition + any blocking factors (batch, lot, patient)
5. Run exomePeak2 differential mode as primary analysis
6. Run QNB (for small-N) or RADAR (for replicate-variance-aware) as sensitivity analysis
7. Compute count matrices via featureCounts if QNB / DESeq2 path is requested
8. Apply effect-size threshold (|log2FC| >= 0.5 minimum) AND FDR < 0.05 AND replicate-direction concordance
9. Cross-check differential signal against input log2FC to flag expression-driven false positives
10. Render volcano + MA plots; per-peak boxplot for top hits
11. Cross-reference top differential peaks against m6A-Atlas / REPIC / GLORI if available
12. Flag the stoichiometry-vs-enrichment caveat in the final report
13. Recommend orthogonal validation route (GLORI / SAC-seq / m6Anet) for any absolute-stoichiometry claim

## Tips

- exomePeak2 differential mode is the field default for standard 3-vs-3; QNB for small-N (2-vs-2); RADAR when replicate variance is a concern.
- N=2 is under-powered. N=3 is the practical minimum (McIntyre 2020). N=4-5 preferred.
- Always include antibody lot / batch / IP day in the design matrix when they vary across conditions.
- MeRIP differential reports CHANGES IN ENRICHMENT RATIO, NOT changes in absolute stoichiometry. Use "increased / decreased enrichment" terminology unless orthogonally validated.
- For absolute stoichiometry claims at named loci, run GLORI / SAC-seq / m6Anet per-read.
- Effect-size filter is non-negotiable. |log2FC| >= 0.5 minimum, >= 1 for high-stakes claims.
- Cross-method concordance is the reproducibility anchor. Run >=2 methods and report concordant set.
- McIntyre 2020 documented ~45% median between-lab peak overlap. Differential calls within this noise envelope routinely do not replicate.
- exomePeak2 has NO `mode=` argument; differential is triggered by populating `bam_treated_ip` + `bam_treated_input` alongside the control-arm `bam_ip` + `bam_input`. `peak_calling_mode` is a separate argument controlling locus scope (`'exon' | 'full_transcript' | 'whole_genome'`).
- featureCounts strand-specificity matters; verify protocol and pass `strandSpecific=` correctly.

## Related Skills

- merip-preprocessing - Upstream IP/input BAM preparation; antibody-lot metadata originates here
- m6a-peak-calling - Peak calling that produces input to differential
- m6anet-analysis - ONT direct-RNA orthogonal validation
- modification-visualization - Volcano / MA / per-peak boxplot rendering
- differential-expression/deseq2-basics - Canonical DE design philosophy
- differential-expression/de-results - Post-DE ranking and interpretation
- chip-seq/differential-binding - General IP-vs-input differential framework
- rna-quantification/featurecounts-counting - Count matrix construction
- data-visualization/volcano-and-ma-plots - Volcano / MA plot recipes
- pathway-analysis/go-enrichment - GO enrichment on differential gene lists
