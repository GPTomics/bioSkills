# Overlap Significance - Usage Guide

## Overview

This skill turns a genomic-interval overlap count into a defensible statistic. A raw overlap ("847 of my 1,000 peaks overlap enhancers") is meaningless without a null model, and because the genome is structured - genes, regulatory elements, and callable space all cluster in the same high-GC, high-mappability territory - almost any two feature sets look "significantly co-localized" under a uniform-random null. The skill covers the full taxonomy of region-set significance methods: bedtools fisher (a fast analytic screen), bedtools shuffle + jaccard permutation, GAT (isochore/GC-conditioned simulation with FDR), regioneR (flexible permutation with clustering-preserving randomization and localZScore), LOLA (universe-relative Fisher against region databases), and GREAT/rGREAT (regulatory-domain ontology enrichment from regions). The dominant lever throughout is the universe/background choice - getting it right matters more than which test you pick.

## Prerequisites

CLI and Python:

```bash
conda install -c bioconda bedtools gat
pip install pybedtools
```

R/Bioconductor:

```r
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
BiocManager::install(c('regioneR', 'LOLA', 'rGREAT'))
```

Notes: pybedtools needs a `bedtools` binary on PATH. GAT requires a workspace (accessible regions) file and benefits from an isochore/GC segmentation. The ENCODE blacklist for your assembly (e.g. hg38) and a chrom-sizes/genome file are needed before any shuffle-based test. Bioconductor packages are pinned to the Bioconductor release - record it with results.

## Quick Start

Tell your AI agent what you want to do:
- "Is the overlap between my peaks and enhancers more than chance?"
- "Run a 1000-permutation colocalization test with regioneR against hg38, excluding the blacklist"
- "Enrich my ATAC peaks against ENCODE TFBS region sets with LOLA"
- "Get GO terms for my peak set with GREAT using a proper background"
- "Test enrichment against many tracks with GC control and FDR using GAT"
- "Are these two CNV call sets concordant at 50% reciprocal overlap?"

## Example Prompts

### Permutation colocalization
> "I have ChIP peaks and a CTCF binding-site BED for hg38. Test whether they co-localize more than expected with regioneR: mask the ENCODE blacklist, use circularRandomizeRegions because my peaks are clustered, 1000 permutations, and report the p-value, z-score, and a localZScore profile."

### Universe-relative enrichment
> "Enrich my ATAC peak set against the LOLA Core hg38 region database, using the union of all my called peak sets as the userUniverse rather than the whole genome, and give me the top enrichments by odds ratio with FDR-adjusted q-values."

### Composition-controlled, many tracks
> "Run GAT to test my peaks against a directory of annotation tracks, with an accessible-genome workspace and GC isochores, 10000 samples, nucleotide-overlap counter, and report fold change, empirical p, and qvalue per track."

### Ontology from regions
> "Run rGREAT on my distal peak set against GO biological process using an accessible-region background, and list only terms significant by both the binomial and hypergeometric tests."

### Triage then validate
> "Quickly screen whether my peaks and enhancers overlap with bedtools fisher, and if it's low, follow up with a proper permutation test - I know fisher's analytic null is weak."

## What the Agent Will Do

1. Establish the **question**: overlap-enrichment of two region sets, region-database enrichment, or ontology-from-regions - and route GWAS/eQTL statistical colocalization to causal-genomics instead.
2. Define the **universe/background** explicitly - the callable/accessible pool the query was drawn from, not the whole genome - because this dominates the result.
3. Prepare inputs: sort all BEDs consistently, harmonize chromosome naming, exclude the ENCODE blacklist and assembly gaps.
4. Pick the method from the decision tree: fisher to triage; regioneR/GAT for permutation; LOLA for region-database enrichment; rGREAT for ontology.
5. Choose a null that preserves region size, the workspace, and (for clustered queries) clustering via circular randomization or isochore conditioning.
6. Run with an adequate permutation count (>=1000) or sample count (>=10000) and apply multiple-testing control across tracks/terms.
7. Report the permutation p-value/z-score (or odds ratio + q), not the raw count, and run the shuffle-your-own-query sanity check.

## Tips

- Spend your effort on the universe/background, not on which library computes the overlap - the tools largely agree on the arithmetic.
- The honest answer is usually *less* significant than the naive test; deflation under a proper null is the method working, not failing.
- bedtools fisher is a screen, never the reported result - validate any low p by simulation.
- Use circularRandomizeRegions (regioneR) when the query is autocorrelated/clustered; plain randomization inflates significance.
- Always exclude the ENCODE blacklist and assembly gaps before shuffling; place regions only within the accessible workspace.
- For "same event" CNV/SV concordance use 50% reciprocal overlap (`-f 0.5 -r`); one-sided fractions let a giant call swallow a tiny one.
- For GREAT, trust a term only if both the binomial and hypergeometric tests are significant, and always supply a real background.

## Related Skills

- interval-arithmetic - The intersect/shuffle/jaccard/fisher mechanics this skill turns into a test
- bed-file-basics - BED format, coordinate systems, and sorting the inputs every test requires
- proximity-operations - Nearest-feature assignment when the question is distance, not overlap enrichment
- chip-seq/peak-calling - Source of the peak query sets tested for enrichment
- chip-seq/peak-annotation - Assign enriched peaks to genes/features
- atac-seq/atac-peak-calling - Source of ATAC peak sets and the accessible-region universe
- pathway-analysis/go-enrichment - Gene-list ontology enrichment; GREAT is the region-based analog
- causal-genomics/colocalization-analysis - GWAS/eQTL statistical colocalization, a distinct problem from interval overlap
- data-visualization/genome-tracks - Render the query and annotation tracks behind an enrichment claim
