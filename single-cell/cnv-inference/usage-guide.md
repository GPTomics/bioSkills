# Copy-Number Inference from Single-Cell RNA-seq - Usage Guide

## Overview

This skill infers large-scale copy-number alterations from tumor single-cell or single-nucleus RNA-seq by treating averaged expression over genomic windows as a proxy for DNA copy number. It separates malignant from normal cells and calls subclones at chromosome-arm / large-segment (~5 Mb) resolution. It covers reference-based expression smoothing (inferCNV), reference-free segmentation (copyKAT, SCEVAN), and haplotype-aware allele-plus-expression inference (Numbat). It is distinct from DNA/WES-based copy-number (copy-number/cnvkit-analysis), which measures DNA read depth and resolves focal events.

## Prerequisites

```r
# inferCNV (Bioconductor)
BiocManager::install('infercnv')

# copyKAT and SCEVAN (GitHub)
remotes::install_github('navinlabcode/copykat')
remotes::install_github('AntonioDeFalco/SCEVAN')

# Numbat (CRAN); preprocessing also needs cellsnp-lite and Eagle2 on the PATH
install.packages('numbat')
```

## Quick Start

Tell your AI agent what you want to do:
- "Which cells in my tumor scRNA-seq are malignant?"
- "Infer chromosome-arm CNVs from my tumor expression matrix"
- "Run inferCNV using my T cells and myeloid cells as the normal reference"
- "Call tumor subclones with allele-aware Numbat"
- "Pick a CNV-inference method when I have no annotated normal cells"

## Example Prompts

### Malignant vs normal
> "Separate malignant from normal cells in this tumor sample using a CNV-inference method"
> "I annotated immune and stromal cells; use them as the inferCNV reference and call which cells are aneuploid"

### Method choice
> "I have no normal cells annotated - which reference-free CNV caller should I use?"
> "Should I use inferCNV, copyKAT, or Numbat for this dataset, and why?"

### Subclones and alleles
> "Run Numbat to resolve subclones and copy-neutral LOH from my BAM and counts"
> "My copyKAT subclones look unstable - how do I confirm them?"

### Interpretation and pitfalls
> "My immune cells are being called aneuploid - what went wrong with my reference?"
> "Run CNV inference per patient instead of on my integrated cross-patient object"

## What the Agent Will Do

1. Confirm the analysis is per sample / per patient, never on a cross-patient integrated object (tumor karyotypes are patient-private)
2. Decide reference-based vs reference-free by whether confident non-malignant cells are annotated in the sample
3. Decide expression-only vs allele-aware by whether subclone resolution or copy-neutral LOH matters and whether a BAM is available
4. For reference-based inferCNV, assemble the counts matrix, annotation file, and gene-ordering file and name the normal groups
5. Run the chosen method with droplet-appropriate parameters, denoising, and (inferCNV) the HMM
6. Classify cells as malignant vs normal and validate the call with lineage markers, mutations, or allele evidence, not the CNV heatmap alone
7. Treat subclone calls as hypotheses and reconcile expression-only calls against allele-aware (Numbat) or DNA evidence

## Tips

- **The reference is everything** - a tumor-contaminated or mismatched normal reference fabricates CNVs in every other cell; pick confident in-sample non-malignant lineages.
- **Run per patient** - integrating across patients before inference erases the patient-private CNV signal the analysis depends on.
- **Expression is a proxy, not DNA** - resolution is chromosome-arm (~5 Mb); for focal amplifications and deletions use DNA-based copy-number (copy-number/cnvkit-analysis).
- **Absence of CNV is not normality** - CNV-quiet and low-grade tumors look flat; confirm malignancy with allele evidence or markers.
- **Malignant calling is a clustering decision** - support it with orthogonal evidence, since the CNV heatmap alone is a noisy hypothesis.
- **Use cutoff 0.1 for droplet data** - 10x is sparse; reserve cutoff 1 for full-length Smart-seq in inferCNV.
- **Anchor copyKAT when mostly aneuploid** - pass known-normal barcodes via norm.cell.names so it can find a diploid baseline.
- **Watch cell-cycle stripes** - proliferating cells create banding that mimics CNV; account for cycle and denoise.
- **Use enough reference cells** - tens to hundreds, not a handful, or the noisy baseline mean fabricates CNVs everywhere.
- **Match reference and tumor sex** - X-inactivation compensates most chrX expression, but an opposite-sex (external/shipped/pooled) reference still fabricates a uniform sex-chromosome CNV via chrY presence/absence, XIST, and escape genes; match sexes or drop the sex chromosomes.
- **Distrust HLA and Ig/TCR segments** - high, variable expression at MHC (6p) and immunoglobulin/TCR loci mimics CNV segments that track lineage, not copy number.
- **Balanced WGD looks copy-neutral** - per-cell normalization hides a uniform doubling from every expression method; confirm ploidy with allele or DNA evidence.

## Related Skills

- single-cell/preprocessing - QC and normalization that precede CNV inference
- single-cell/clustering - Provides the cell groups and the malignant-vs-normal clustering the CNV call refines
- single-cell/cell-annotation - Identifies the non-malignant lineages used as the normal reference
- single-cell/batch-integration - Why CNV inference must run per patient before any cross-patient integration
- copy-number/cnvkit-analysis - DNA/WES-based copy-number that resolves focal events, the orthogonal contrast to this expression proxy
