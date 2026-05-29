# experimental-design

## Overview

Skills for designing genomics experiments so that statistical inference is valid by construction: choosing the experimental unit, randomizing and blocking, balancing technical variation, powering the study, sizing the sample, and controlling error rates across thousands of tests. The category covers the discovery/high-dimensional-omics regime (FDR, genes/peaks/cells); for the confirmatory/regulated-trial regime (FWER, closed testing, regulated sample size) see clinical-biostatistics.

**Tool type:** r | **Primary tools:** designit, RNASeqPower, ssizeRNA, qvalue, sva

## Skills

| Skill | Description |
|-------|-------------|
| randomization-blocking | Experimental unit, pseudoreplication, randomization, blocking, factorial/split-plot/nested designs |
| batch-design | Balancing technical variation against biology; no-rescue theorem; correction-method choice; SVA detection |
| power-analysis | Per-gene negative-binomial power for sequencing assays; simulation; depth-vs-replicate; post-hoc-power fallacy |
| sample-size | Minimum biological replicates at a target power and FDR; pilot dispersions; scRNA-seq donors-vs-cells |
| multiple-testing | FDR vs FWER; BH/BY dependence; q-value/pi0; local FDR; IHW; independent filtering; GWAS threshold |

## Example Prompts

- "What is my experimental unit, and is my n the number of mice or the number of cells?"
- "Help me assign 24 samples to 3 sequencing batches so batch is orthogonal to condition."
- "All my tumors were sequenced in one run and normals in another; can ComBat fix this?"
- "How many biological replicates for RNA-seq to detect 1.5-fold changes at 80% power, FDR 0.05?"
- "Should I sequence deeper or add samples on a fixed budget?"
- "How many donors versus cells for a scRNA-seq disease-versus-control comparison?"
- "Should I use BH, BY, or q-value when my test statistics are correlated?"
- "Use IHW with mean expression to gain power over plain Benjamini-Hochberg."
- "A reviewer says my study is underpowered because observed power was 0.3; how do I respond?"
- "Is it valid to filter out low-count genes before testing to boost power?"

## Requirements

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install(c('RNASeqPower', 'ssizeRNA', 'PROPER', 'qvalue', 'IHW', 'sva', 'RUVSeq', 'OSAT', 'limma', 'edgeR', 'DESeq2'))

# CRAN
install.packages(c('designit', 'lme4', 'lmerTest', 'pwr'))

# powsimR for scRNA-seq power (GitHub; pin a commit)
# remotes::install_github('bvieth/powsimR')
```

```bash
# Python equivalents
pip install statsmodels scipy numpy pandas
```

## Related Skills

- **clinical-biostatistics** - Confirmatory/regulated-trial power, sample size, and multiplicity (FWER, closed testing)
- **differential-expression** - Runs the DE test and executes batch correction (ComBat-seq, RUVSeq, SVA)
- **single-cell** - scRNA-seq preprocessing, pseudobulk aggregation, and batch integration
- **machine-learning** - Model validation and data leakage, of which batch confounding is one source
- **read-qc** - Quality control that verifies an experiment's design succeeded
