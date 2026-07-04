# Gene Family Evolution - Usage Guide

## Overview

Gene-family birth-death models on phylogenies (Hahn 2005; Csurös 2010) detect lineage-specific gene-family expansions and contractions. **CAFE5** (Mendes 2020 Bioinformatics 36:5516) is the modern standard with gamma-distributed rate categories for biologically realistic modeling. **Annotation heterogeneity is the dominant confounder** -- different annotation pipelines predict different gene counts, producing apparent "expansions" that are pipeline artifacts. Consistent annotation + BUSCO/Compleasm completeness filtering are mandatory.

For HGT-affected gene families (prokaryotes), CAFE5's birth-death-only model misses transfer; use ALE/GeneRax/AleRax instead ([[gene-tree-species-tree-reconciliation]]).

## Prerequisites

```bash
# Modern CAFE5
conda install -c bioconda cafe

# Alternative methods
# Count (Csurös 2010; Java)
wget http://www.iro.umontreal.ca/~csuros/gene_content/count.tar.gz && tar xf count.tar.gz

# BadiRate (Librado 2012)
git clone https://github.com/PauloRoldan/badirate

# Tree time-calibration
conda install -c bioconda treepl  # or ape::chronos in R
# Or LSD2 (faster)

# Time-source database
# TimeTree (web): http://www.timetree.org

# Functional enrichment
Rscript -e "
install.packages('BiocManager')
BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'topGO'))
"

# OrthoFinder input
conda install -c bioconda orthofinder
```

## Quick Start

Tell the AI agent what to do:
- "Run CAFE5 on OrthoFinder HOG counts to identify gene families with lineage-specific expansions"
- "Identify ancestral gene-family counts at internal nodes using Count"
- "Run CAFE5-error mode with explicit annotation error rates per species"
- "Test for trait-correlated rate shifts in gene families using RERconverge"

## Example Prompts

### Standard CAFE5 Analysis

> "I have OrthoFinder HOG counts for 30 mammalian species. Build a time-calibrated species tree with treePL using divergence dates from TimeTree. Run CAFE5 with 4 gamma rate categories and identify families with lineage-specific rate shifts (p < 0.05, FDR-corrected). Annotate top expanded / contracted families with functional categories."

### Annotation-Heterogeneity Correction

> "The 30 mammalian genomes have heterogeneous annotation pipelines (some annotated with NCBI RefSeq, others with Ensembl, others with BRAKER3). Run CAFE5-error mode with empirical per-species annotation error rates. Compare results with strict CAFE5 to assess heterogeneity impact."

### Convergent Rate Shifts

> "Identify gene families with rate shifts in echolocating species (bats and dolphins) using RERconverge on the CAFE5 family rates. Cross-validate with CSUBST for convergent substitution patterns. Report top convergent gene families with functional annotation."

## What the Agent Will Do

1. **Validate inputs**: OrthoFinder HOG matrix; species tree; per-genome BUSCO completeness
2. **Filter genomes** with BUSCO < 90% (or document)
3. **Normalize annotation pipeline** if possible; report error rates if not
4. **Time-calibrate species tree** with treePL / LSD2 / TimeTree
5. **Run CAFE5** with gamma rate categories (`-k 4`)
6. **CAFE5-error mode** if annotation heterogeneity is present
7. **FDR-correct** across families
8. **Identify lineage-specific expansions / contractions** at significant rate shifts
9. **Functional enrichment** with clusterProfiler / topGO
10. **Cross-validate** with Count for ancestral state reconstruction; ALE for HGT-affected families
11. **Report**: per-family lambda, significant families, functional categories, ancestral state estimates
12. **Caveats**: annotation heterogeneity, assembly fragmentation, HGT, polyploidy

## Tips

- CAFE5 (Mendes 2020) is the modern standard; CAFE4 deprecated
- Filter input to multi-copy orthogroups (max count >= 2); single-copy contribute no information
- Tree MUST be ultrametric (time-calibrated); use treePL / LSD2 / ape::chronos
- Negative branch lengths cause failures; calibrate properly
- Require >= 100 families for CAFE5 (preferably > 1000)
- Annotation heterogeneity is the single largest confounder; re-annotate consistently if possible
- BUSCO completeness >= 90% mandatory; exclude or correct fragmented genomes
- Multi-testing correction: FDR (Benjamini-Hochberg) across families
- For HGT-affected bacterial families, switch to ALE / GeneRax / AleRax
- For polyploid species, assign subgenomes first; analyze each subgenome separately
- For convergent rate shifts, RERconverge complements CAFE5; CSUBST for amino-acid convergence
- Whale.jl native WGD modeling for WGD-affected lineages
- Report robust lambda with and without outlier families; manual annotation of outliers
- Time calibration from TimeTree is convenient; or compute from molecular clock

## Related Skills

comparative-genomics/ortholog-inference - OrthoFinder HOG matrix is CAFE5 input
comparative-genomics/gene-tree-species-tree-reconciliation - ALE per-family DTL, complement
comparative-genomics/whole-genome-duplication - Post-WGD retention bias
comparative-genomics/positive-selection - Selection within expanded families
comparative-genomics/ancestral-reconstruction - Ancestral count reconstruction
phylogenetics/divergence-dating - Time-calibrated tree for CAFE5
pathway-analysis/go-enrichment - Functional enrichment of expanded families
pathway-analysis/gsea - GSEA on family expansions
