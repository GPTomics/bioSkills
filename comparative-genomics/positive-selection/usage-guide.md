# Positive Selection - Usage Guide

## Overview

dN/dS-based detection of positive selection uses codon models to identify gene-wide, branch-specific, or per-site adaptive evolution. **PAML codeml** (Yang 1997) provides site, branch, and branch-site model A (Zhang 2005); **HyPhy** provides BUSTED, BUSTED-S, **BUSTED-MH** (Lucaci 2023 multi-hit aware), MEME, FEL, FUBAR, aBSREL, RELAX, and GARD. **GARD recombination pre-screen is mandatory** before any selection test; recombination breakpoints inflate false positives 5-50x. Alignment errors (PRANK codon-aware MSA + PREQUAL / HmmCleaner segment filter) and gBGC are the next-largest confounders.

For McDonald-Kreitman framework, **asymptotic alpha** (Messer-Petrov 2013) corrects for slightly deleterious bias; **polyDFE / GRAPES** provide full DFE inference; **impMKT** improves gene-level evidence (Murga-Moreno 2022).

## Prerequisites

```bash
# Codon model engines
conda install -c bioconda paml hyphy gard
# RDP5 (viral recombination): http://web.cbio.uct.ac.za/~darren/rdp.html

# Codon-aware alignment
conda install -c bioconda prank

# Alignment filtering (REQUIRED)
conda install -c bioconda prequal hmmcleaner

# McDonald-Kreitman framework
Rscript -e "install.packages('asymptoticMK')"
# polyDFE: https://github.com/paula-tataru/polyDFE
# GRAPES: https://github.com/BioPP/grapes

# Convergent evolution
Rscript -e "remotes::install_github('nclark-lab/RERconverge')"
pip install csubst
# PhyloAcc: https://phyloacc.github.io/

# Python
pip install ete4 pyhyphy
```

## Quick Start

Tell the AI agent what to do:
- "Test for positive selection in this gene across mammals using PAML site models + HyPhy MEME"
- "Run GARD recombination screen before any selection test"
- "Apply BUSTED-MH (multi-hit aware) to handle viral / Plasmodium / Trypanosoma genes"
- "Compute asymptotic alpha for adaptive substitution rate from polymorphism + divergence data"

## Example Prompts

### Standard Gene Selection Test

> "Test gene X for positive selection across 20 mammalian orthologs. Pipeline: PRANK codon-aware MSA -> PREQUAL filter -> GARD recombination screen -> BUSTED-MH (gene-wide) + MEME (sites) + aBSREL (branches). FDR-correct across multiple genes. Cross-validate with PAML M8 vs M8a. Report selected sites and lineages."

### Genome-Wide Selection Scan

> "Run a genome-wide selection scan on 8000 single-copy orthologs across 30 species. Pipeline: codeml M0 baseline + BUSTED-MH per gene + MEME for sites. FDR-correct across genes. Report top adaptive genes with multi-method support."

### McDonald-Kreitman Framework

> "I have polymorphism (SFS) and divergence data for 5000 Drosophila genes. Compute asymptotic alpha (corrected for slightly deleterious bias) and the full DFE via polyDFE. Report per-gene alpha and 95% CI. Cross-validate with GRAPES."

## What the Agent Will Do

1. **Pre-process alignment**: PRANK codon-aware MSA; PREQUAL segment filter (mandatory)
2. **Pre-screen recombination**: GARD (mandatory before any selection test)
3. **Gene-wide screen**: BUSTED-MH (multi-hit aware)
4. **Site-level test**: MEME (episodic) + FEL/FUBAR (pervasive)
5. **Branch-level test**: aBSREL for unspecified branches; codeml branch-site for pre-specified
6. **Cross-validate**: PAML vs HyPhy concordance
7. **FDR-correct** across genes
8. **For population genetics**: asymptotic alpha + polyDFE + GRAPES
9. **For convergent evolution**: CSUBST (amino acid) + RERconverge (categorical trait)
10. **Report**: significant sites + branches + genes; per-method support; multi-testing correction; saturation check
11. **Caveats**: recombination, alignment error, gBGC, multi-hit, saturation, foreground specification

## Tips

- GARD pre-screen is MANDATORY before any selection test (Anisimova 2003)
- Use PRANK (Loytynoja 2014) for codon-aware MSA; MACSE V2 for frameshift-tolerant
- Apply PREQUAL or HmmCleaner segment-level filter; block filters (Gblocks) lose information
- BUSTED-MH (Lucaci 2023) corrects for multi-hit substitutions; preferred over basic BUSTED for viral / Plasmodium / Trypanosoma
- Branch-site test LRT uses 50:50 chi-square mixture (critical value 2.71 at p=0.05, NOT 3.84)
- For genome scans, FDR (Benjamini-Hochberg) across genes; site-level p<0.1 default
- aBSREL for branch-specific selection without foreground pre-specification (adaptive multi-test)
- For specified foreground, codeml branch-site mod A vs A1
- Check dS < 1.5 per branch; saturated dS makes codon analysis unreliable
- gBGC inflates dN/dS at high-recombination regions; W->S substitution ratio diagnostic
- BEB posterior > 0.95 for significant sites; > 0.99 highly significant
- HyPhy MEME is the only method for per-site episodic selection
- For non-coding accelerated evolution, use PhyloAcc / phyloP-acc (not codon)
- For convergent substitution, CSUBST / RERconverge
- McDonald-Kreitman asymptotic alpha (Messer-Petrov 2013) handles slightly deleterious bias
- polyDFE / GRAPES for full DFE + demographic correction

## Related Skills

comparative-genomics/ortholog-inference - Single-copy ortholog alignments
comparative-genomics/ancestral-reconstruction - Ancestral sequence at selected branches
comparative-genomics/gene-tree-species-tree-reconciliation - Reconciled gene trees as input
alignment/multiple-alignment - PRANK / MACSE codon-aware MSA
alignment/alignment-trimming - PREQUAL / HmmCleaner filtering
phylogenetics/modern-tree-inference - Tree inference for codeml
population-genetics/selection-statistics - SFS-based alpha + DFE
causal-genomics/heritability-partitioning - LDSC partition with selection annotations
variant-calling/variant-annotation - Functional annotation of selected sites
