---
name: bio-ribo-seq-translation-efficiency
description: Quantify translation efficiency (TE) as ribosome occupancy relative to mRNA abundance and test for differential TE between conditions. Use when separating translational from transcriptional regulation, distinguishing genuine translational control from buffering, or choosing between riborex, Xtail, anota2seq, and DESeq2 interaction models.
tool_type: mixed
primary_tool: riborex
---

## Version Compatibility

Reference examples tested with: riborex 2.4+, xtail 1.1+, anota2seq 1.24+, DESeq2 1.42+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Translation Efficiency

**"Calculate translation efficiency from my Ribo-seq and RNA-seq"** -> Compute footprint density relative to mRNA density per gene and test which genes change translation independently of transcription, distinguishing real translational control from buffering.
- R: `riborex` (DESeq2/edgeR backend), `Xtail`, or `anota2seq` for differential TE
- Python: per-gene TE ratio for ranking/visualization only

TE = RPF density / mRNA density over the SAME region. Both assays must come from matched samples and be counted over the CDS. TE isolates translational regulation and is a relative translation-rate proxy at steady state; occupancy is not protein output.

## The central trap: a ratio is for ranking, not testing

The naive per-gene ratio (TPM_ribo/TPM_rna) is fine for ranking and plots but WRONG for differential testing: it ignores count heteroskedasticity, treating a gene with 5 reads like one with 5000. Differential TE is NOT "compute TE per condition then test the difference" -- it is a CONDITION x ASSAY INTERACTION on raw counts with proper negative-binomial dispersion modeling, where log2FC(TE) = log2FC(RPF) - log2FC(mRNA). The whole differential-TE field exists to do this interaction correctly.

## Mode of regulation: control vs buffering

When both assays move, there are distinct biological modes that a single TE fold-change cannot separate:

| Mode | RPF | mRNA | TE | Meaning |
|------|-----|------|-----|---------|
| mRNA abundance | up | up | ~flat | transcriptional, not translational |
| translation (forwarded) | up | flat | up | genuine translational control -> protein changes |
| buffering | flat | up | down | translation absorbs the mRNA change, protein held constant |

Buffering (a homeostatic mechanism) and genuine translational activation can produce the SAME |log2FC(TE)|. Calling a buffered gene "translationally activated" is a wrong conclusion. Only anota2seq formally names the mode, by regressing translated mRNA on total mRNA (analysis of partial variance).

## Differential-TE tool selection

| Tool | Statistic | Names buffering | Best when |
|------|-----------|-----------------|-----------|
| riborex | wraps DESeq2/edgeR/Voom on a merged interaction design | no | fast drop-in for DESeq2 users |
| Xtail | two pipelines (FC-vs-FC, ratio-vs-ratio), reports the more conservative | partial (won't call a buffered gene a hit) | conservative differential-TE calls + clean plots |
| anota2seq | per-mRNA APV + random variance model | YES | mode-of-regulation biology; the postdoc-grade choice |
| RiboDiff | NB GLM, shared dispersion by default | no | few replicates; CLI pipeline |
| DESeq2 interaction | ~assay+condition+assay:condition, Wald or LRT | no (post-hoc) | full control, custom contrasts, batch terms |

## Quick per-gene TE (ranking screen only)

**Goal:** Rank genes by TE for a quick look, not for inference.

**Approach:** Normalize both assays, take the log2 ratio over the CDS with a pseudocount.

```python
import numpy as np

def log2_te(ribo_cds_tpm, rna_cds_tpm, pseudocount=0.1):
    '''Per-gene log2 TE for ranking/plots. Both inputs counted over the CDS.

    Pseudocount 0.1 TPM avoids log(0) and dampens low-count noise. Not for testing.
    '''
    return np.log2((ribo_cds_tpm + pseudocount) / (rna_cds_tpm + pseudocount))
```

Count BOTH assays over the CDS. Using full-transcript RNA against CDS-only RPF introduces a UTR-length confound (long-UTR genes look low-TE). Exclude the first ~15 and last ~5 codons of the CDS so initiation and termination peaks do not dominate the RPF count.

## Differential TE with riborex

**Goal:** Test differential TE reusing a familiar DE engine.

**Approach:** Pass CDS count matrices and condition vectors; riborex builds the interaction design internally and returns DESeq2-format results.

```r
library(riborex)

# rna_counts / ribo_counts: genes x samples integer CDS counts
res <- riborex(rnaCntTable = rna_counts, riboCntTable = ribo_counts,
               rnaCond = c("ctrl", "ctrl", "treat", "treat"),
               riboCond = c("ctrl", "ctrl", "treat", "treat"),
               engine = "DESeq2")
sig <- res[which(res$padj < 0.05), ]   # log2FoldChange is the TE change
```

Engines are `"DESeq2"` (default), `"edgeR"`, `"edgeRD"`, `"Voom"` (Voom is single-factor only).

## Differential TE with anota2seq (names the mode)

**Goal:** Separate translation, buffering, and mRNA-abundance regulation.

**Approach:** Provide translated (RPF) and total (RNA) matrices, run the pipeline, then classify each gene's mode.

```r
library(anota2seq)

ads <- anota2seqDataSetFromMatrix(dataP = ribo_counts, dataT = rna_counts,
                                  phenoVec = c("ctrl", "ctrl", "treat", "treat"),
                                  dataType = "RNAseq", normalize = TRUE)
ads <- anota2seqRun(ads, useRVM = TRUE)
ads <- anota2seqRegModes(ads)   # one mode per gene: translation > abundance > buffering
translation_hits <- anota2seqGetOutput(ads, analysis = "translation",
                                       output = "selected", selContrast = 1)
```

## Differential TE with a DESeq2 interaction

**Goal:** Full control over the interaction model.

**Approach:** Merge RPF and RNA counts, fit the interaction, and select the interaction coefficient by name from `resultsNames` (never hardcode it).

```r
library(DESeq2)
counts <- cbind(ribo_counts, rna_counts)
coldata <- data.frame(
    condition = factor(rep(c("ctrl", "ctrl", "treat", "treat"), 2)),
    assay = factor(rep(c("ribo", "rna"), each = 4)))
dds <- DESeqDataSetFromMatrix(counts, coldata, ~ assay + condition + assay:condition)
dds <- DESeq(dds)

# The interaction name is auto-generated from factor levels; pick it programmatically.
# DESeq2 renders interaction coefficients with a DOT (e.g. assayrna.conditiontreat),
# while main effects use underscores -- so match the dot, not the formula's colon.
nm <- grep("\\.", resultsNames(dds), value = TRUE)
res_te <- results(dds, name = nm)
```

Size factors are estimated PER ASSAY (ribo among ribo, RNA among RNA); the implicit assumption is that the median gene's TE is unchanged. If a global translational shift is expected (e.g. mTOR inhibition), median normalization is violated and spike-ins are needed to anchor absolute scale.

## Confounders to check

mRNA isoform switching changes the CDS/UTR counting region between conditions; UTR changes that alter uORF usage can make a main-ORF TE change SECONDARY to uORF regulation rather than direct translational control. Cross-check called ORFs and uORFs (see orf-detection) before attributing a TE shift to the main ORF.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Low-count genes dominate the hit list | t-test/ratio on log-TE | Use count-based GLM (riborex/Xtail/anota2seq/DESeq2) |
| Long-UTR genes systematically low TE | RNA counted over full transcript, RPF over CDS | Count both over the CDS |
| `results(dds, name='conditiontreat.assayribo')` errors | Hardcoded interaction name | Select from `resultsNames(dds)` by the "." term (interaction coefficients render with a dot, not the formula's colon) |
| Unstable dispersion or anota2seq RVM warnings | Too few replicates (n=2 as in the examples) | Use >=3 replicates per condition per assay; n=2 is illustrative only |
| Buffered gene reported as translationally activated | Single TE fold-change cannot separate modes | Use anota2seq mode-of-regulation |
| TE shifts vanish or invert globally | Global translational change breaks median normalization | Add spike-ins; do not assume median TE unchanged |
| Initiation peak inflates RPF counts | Whole-CDS counting including start/stop peaks | Trim first ~15 / last ~5 codons |

## Related Skills

- ribosome-periodicity - Calibrate P-site offsets for CDS footprint counts
- orf-detection - Rule out uORF-driven (secondary) TE changes
- rna-quantification/featurecounts-counting - Generate matched RNA-seq CDS counts
- differential-expression/deseq2-basics - Count-based DE foundations

## References

- Li W, Wang W, Uren PJ, Penalva LOF, Smith AD. 2017. Riborex: fast and flexible identification of differential translation from Ribo-seq data. Bioinformatics 33(11):1735-1737. doi:10.1093/bioinformatics/btx047
- Xiao Z, Zou Q, Liu Y, Yang X. 2016. Genome-wide assessment of differential translations with ribosome profiling data. Nat Commun 7:11194. doi:10.1038/ncomms11194
- Oertlin C, Lorent J, Murie C, Furic L, Topisirovic I, Larsson O. 2019. Generally applicable transcriptome-wide analysis of translation using anota2seq. Nucleic Acids Res 47(12):e70. doi:10.1093/nar/gkz223
- Zhong Y, Karaletsos T, Drewe P, et al. 2017. RiboDiff: detecting changes of mRNA translation efficiency from ribosome footprints. Bioinformatics 33(1):139-141. doi:10.1093/bioinformatics/btw585
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15(12):550. doi:10.1186/s13059-014-0550-8
