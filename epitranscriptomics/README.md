# epitranscriptomics

## Overview

Detection, quantification, differential analysis, and visualisation of RNA modifications (primarily m6A / N6-methyladenosine, with secondary coverage of m6Am and m5C) from MeRIP-seq / m6A-seq antibody-IP data and Oxford Nanopore direct-RNA-sequencing (ONT DRS) signal data. Covers the full pipeline from raw FASTQ / POD5 to publication figures with explicit handling of the antibody-vs-chemistry tradeoff and the m6A-vs-m6Am cross-reactivity that defines the field's main analytical pitfall.

**Tool type:** mixed | **Primary tools:** exomePeak2, MeTPeak, MACS3, m6Anet, xPore, Guitar, deepTools

## Skills

| Skill | Description |
|-------|-------------|
| merip-preprocessing | Trim, splice-aware-align (STAR / HISAT2), and QC paired IP/input MeRIP-seq libraries; saturation curves; explicit do-NOT-dedup convention for non-UMI MeRIP |
| m6a-peak-calling | Call m6A peaks with exomePeak2 / MeTPeak / MACS3; confirm DRACH motif as a sanity check (NOT a filter); flag 5'UTR m6A-vs-m6Am ambiguity |
| m6a-differential | Differential m6A between conditions via exomePeak2 / QNB / RADAR / MeTDiff; handle stoichiometry-vs-expression-vs-IP-efficiency confound; effect-size filtering |
| m6anet-analysis | ONT direct-RNA m6A calling with m6Anet, xPore, Nanocompore, ELIGOS, and Dorado native modification calling; per-site vs per-read interpretation; DRACH-only constraint |
| modification-visualization | Guitar metagene (stop-codon enrichment QC anchor), peak-centred heatmaps, IP/Input paired browser tracks (pyGenomeTracks / ggcoverage), DRACH sequence logos, 5'UTR / CDS / 3'UTR stacked bars |

## Example Prompts

- "Align my paired MeRIP-seq IP and input FASTQ libraries with STAR splice-aware; do NOT deduplicate (non-UMI MeRIP); generate replicate-correlation and plotFingerprint IP-enrichment QC; build saturation curves so peak counts are comparable across libraries"
- "Call m6A peaks with exomePeak2 against GRCh38 + GENCODE; cross-check with MeTPeak; confirm DRACH enrichment via HOMER as a sanity check; flag 5'UTR peaks as m6A-vs-m6Am ambiguous"
- "Run differential m6A analysis between 3 control and 3 treatment IP/input pairs with exomePeak2 differential mode; include antibody lot as a fixed effect; apply |log2FC| >= 0.5 and FDR < 0.05; cross-validate top hits against published GLORI sites"
- "Detect m6A from RNA004 ONT direct-RNA-sequencing with the full m6Anet pipeline (Dorado basecalling -> minimap2 transcriptome alignment -> nanopolish eventalign -> m6anet dataprep -> m6anet inference); filter for n_reads >= 20 AND probability_modified >= 0.9"
- "Compare m6A between WT and METTL3-KO with xPore Bayesian diffmod on ONT direct-RNA data"
- "Render the canonical Guitar metagene plot confirming stop-codon enrichment; pair with the 5'UTR / CDS / 3'UTR stacked bar; render a DRACH sequence logo from peak centres"
- "Build a pyGenomeTracks browser figure of paired IP / Input / log2 IP-over-Input tracks at the METTL3 locus"

## Requirements

```bash
# Short-read MeRIP-seq pipeline
conda install -c bioconda star hisat2 samtools deeptools preseq fastp multiqc macs3 macs2 homer bedtools picard

# ONT direct-RNA pipeline
conda install -c bioconda minimap2 nanopolish
pip install m6anet ont-fast5-api xpore nanocompore
# Dorado: download binary from https://github.com/nanoporetech/dorado

# R / Bioconductor
BiocManager::install(c(
    'exomePeak2', 'GenomicFeatures', 'BSgenome.Hsapiens.UCSC.hg38',
    'TxDb.Hsapiens.UCSC.hg38.knownGene', 'rtracklayer', 'GenomicRanges',
    'GenomicAlignments', 'ChIPseeker', 'Guitar', 'ComplexHeatmap',
    'DESeq2', 'edgeR', 'Rsubread'
))

# GitHub-only R packages
devtools::install_github('compgenomics/MeTPeak')
devtools::install_github('lzcyzm/QNB')
devtools::install_github('scottzijiezhang/RADAR')

# Visualisation R packages
install.packages(c('ggcoverage', 'ggseqlogo', 'ggplot2', 'dplyr', 'circlize'))
```

## Related Skills

- **chip-seq** - IP-vs-input peak-calling concepts (MACS3 lineage); FRiP / fingerprint QC; browser-track + peak-centred heatmap patterns transfer directly
- **clip-seq** - miCLIP / m6A-CLIP antibody-CLIP single-nucleotide methods for orthogonal MeRIP validation
- **differential-expression** - Design-matrix philosophy for paired / interaction / batch designs (m6a-differential defers here)
- **read-alignment** - General STAR / HISAT2 / BWA splice-aware alignment
- **long-read-sequencing** - Dorado / Guppy basecalling, signal QC, general ONT mechanics upstream of m6anet-analysis
- **rna-quantification** - featureCounts / count-matrix construction for downstream differential
- **data-visualization** - General ggplot2 / heatmap / browser-track / sequence-logo primitives
- **variant-calling** - A-to-I editing (REDItools / JACUSA) lives there; can confound m6A calls at overlapping loci
- **pathway-analysis** - GO / GSEA / KEGG enrichment on peak-bearing gene lists downstream
