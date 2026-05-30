# CNV Annotation Usage Guide

## Overview

CNV annotation adds biological and clinical context to copy number calls: which genes a segment overlaps, whether those genes are dosage-sensitive, whether the CNV is common in the population, and whether it carries known cancer drivers. The central principle is that overlap is not consequence: a CNV spanning a gene matters only if the gene is dosage-sensitive, the affected portion is functionally important, and (for focal events) the gene is a driver rather than a passenger. Annotated CNVs feed clinical classification (germline-cnv-interpretation) and driver analysis (recurrent-cnv).

## Prerequisites

```bash
conda install -c bioconda bedtools annotsv
pip install pybedtools pandas pysam
# R route: BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
```

Inputs: CNV segments (CNVkit `.cns`, a BED, or a VCF); a gene model (GENCODE/RefSeq); optionally the ClinGen dosage map, COSMIC Cancer Gene Census, gnomAD-SV/DGV catalogs, and a ClinVar VCF. AnnotSV bundles most reference databases on first download.

## Quick Start

Tell the AI agent what to do:
- "Annotate my CNV segments with overlapping genes and dosage-sensitivity scores"
- "Run AnnotSV on my CNV VCF for comprehensive annotation and ranking"
- "Filter out common population CNVs using gnomAD-SV with reciprocal overlap"
- "Identify the likely driver gene of this focal amplification, not all 30 it spans"
- "Flag whether each gene is wholly or only partially deleted"

## Example Prompts

### Gene and dosage annotation

> "Annotate these CNV segments with overlapping genes, record what fraction of each gene is covered, and tag genes with ClinGen haploinsufficiency and triplosensitivity scores."

> "Run AnnotSV on this CNV VCF in GRCh38 and explain how to read the ACMG rank."

### Driver vs passenger

> "This focal amplification spans 30 genes. Identify the likely driver using the COSMIC Cancer Gene Census and explain why the passengers should not all be reported."

### Filtering and context

> "Filter my CNV callset against gnomAD-SV using 50% reciprocal overlap and explain why reciprocal overlap matters for benign filtering."

> "My CNVs are on GRCh37 and my gene model is GRCh38. Diagnose the build mismatch and recommend a safe liftOver approach."

## What the Agent Will Do

1. Verify CNV and annotation genome builds match
2. Intersect CNV segments with gene models, recording overlap extent and exon hits
3. Tag genes with ClinGen dosage scores; flag whole-gene vs partial overlap
4. Annotate cancer driver role (CGC/OncoKB) and flag direction-consistent driver hits
5. Filter against population SV catalogs with reciprocal overlap
6. Run AnnotSV for comprehensive one-pass annotation and triage ranking
7. Hand off to germline-cnv-interpretation or recurrent-cnv as appropriate

## Tips

- Overlap is not consequence: require dosage evidence (ClinGen HI/TS) and record whole-gene-vs-partial status before calling a gene affected.
- For focal amplifications, report the driver (recurrence peak + known oncogene), not every gene the amplicon spans.
- Always confirm CNVs and the gene model are on the same genome build; a mismatch silently returns wrong genes.
- Use reciprocal overlap (`-f 0.5 -r`) for population-frequency filtering so a tiny CNV is not matched to a huge population CNV.
- Parse ClinVar `CLNSIG` against its controlled vocabulary; SNV/indel pathogenicity does not transfer to a CNV.
- CNV-gene pathway enrichment is biased by gene-dense CNV-prone loci; treat it as hypothesis-generating.

## Related Skills

- copy-number/germline-cnv-interpretation - ACMG/ClinGen points-based CNV classification
- copy-number/recurrent-cnv - GISTIC2 recurrence peaks and driver identification
- copy-number/cnvkit-analysis - Generates the CNV segments to annotate
- copy-number/cnv-visualization - Visualizing annotated CNVs
- pathway-analysis/go-enrichment - GO/KEGG enrichment methodology
- genome-intervals/bed-file-basics - BED interval operations
- clinical-databases/clinvar-lookup - Querying ClinVar for variant pathogenicity
