# chip-seq

## Overview

Postdoc-grade ChIP-seq analysis covering peak calling, QC, differential binding, motif discovery, peak annotation, super-enhancers, visualization, CUT&RUN/CUT&Tag, spike-in normalization, chromatin state segmentation, deep-learning models, and allele-specific binding. Embeds the antibody-validation framework, hyper-ChIPable artifact awareness, three-distinct-normalization-problems framing, ENCODE Nself/Nt consistency rules, and current SOTA (chromBPNet, EnFormer, SEACR, ChromHMM v1.27, JASPAR 2026).

**Tool type:** mixed | **Primary tools:** MACS3 (CLI), HOMER (CLI), SEACR (CLI), ChIPseeker (R), DiffBind (R), csaw (R), ChromHMM (Java), chromBPNet (Python), BaalChIP (R), WASP (Python), deepTools (CLI), pyGenomeTracks (CLI), MEME suite (CLI)

## Skills

| Skill | Description |
|-------|-------------|
| peak-calling | MACS3/MACS2/HOMER/SPP narrow & broad peak calling with ENCODE IDR vs naive overlap; per-tool failure modes |
| chipseq-qc | FRiP, NSC/RSC, NRF/PBC1/PBC2, deepTools fingerprint (JS/AUC), ChIPQC, fragment-size diagnostic, hyper-ChIPable detection, antibody validation framework |
| differential-binding | DiffBind, DESeq2, edgeR, csaw (sliding windows), NormR, MAnorm2 with three-distinct-normalization-problems framing (composition / trended / global shift) |
| motif-analysis | HOMER, MEME-ChIP (STREME / CentriMo / TOMTOM / FIMO), monaLisa, AME with background-selection theory; TF-MoDISco cross-reference |
| peak-annotation | ChIPseeker, HOMER, pyranges with nearest-TSS vs host-gene, ENCODE cCRE classification (PLS/pELS/dELS/CTCF-only/DNase-H3K4me3), GREAT, ChIP-Enrich |
| super-enhancers | ROSE/ROSE2/LILY/HOMER -style super with H3K27ac vs MED1 vs BRD4 marker choice; CRC reconstruction (Saint-Andre 2016); spike-in for cross-condition |
| chipseq-visualization | deepTools, pyGenomeTracks, Gviz, EnrichedHeatmap, ChIPseeker, IGV batch with bigWig normalization decisions (CPM/BPM/RPGC/spike-in scaled) |
| cut-and-run-tag | CUT&RUN (Skene Henikoff 2017) and CUT&Tag (Kaya-Okur 2019) with SEACR vs MACS2 consensus per btaf375 2025 benchmark; E. coli spike-in carryover |
| spike-in-normalization | ChIP-Rx (Drosophila); RRPM vs Rx-Input; integration with DiffBind/DESeq2/edgeR/csaw; Hammond Norris 2024 failure modes; ChIPseqSpikeInFree post-hoc |
| chromatin-state-segmentation | ChromHMM v1.27, Segway, EpiSegMix, EpiLogos, IDEAS, full-stack ChromHMM (Vu Ernst 2022); state-count selection logic |
| chip-deep-learning | BPNet, chromBPNet, EnFormer, DeepSEA, JASPAR 2026 Deep Learning collection; in silico mutagenesis; TF-MoDISco motif syntax discovery |
| allele-specific-binding | WASP filter (mandatory), RASQUAL, BaalChIP, AlleleSeq with imprinted-loci + X-inactivation + copy-number-imbalance pitfalls |

## Example Prompts

- "Call peaks from chip.bam with input.bam using ENCODE-compliant MACS2 + IDR pipeline"
- "Compute FRiP, NSC, RSC, NRF/PBC1/PBC2 and verify fragment-size diagnostic before peak calling"
- "Detect hyper-ChIPable artifacts by intersecting peaks against top-1% input signal regions"
- "Run differential binding for HDAC-inhibitor experiment with spike-in normalization (global shift)"
- "Use csaw sliding windows for global-shift-robust differential binding on H3K27me3"
- "Find motifs centrally enriched at CTCF peak summits with MEME-ChIP CentriMo (test direct vs tethered binding)"
- "Annotate peaks to ENCODE cCREs (PLS / pELS / dELS / CTCF-only)"
- "Call super-enhancers with ROSE2; spike-in normalize before cross-condition BET-inhibitor comparison"
- "Generate spike-in-scaled bigWigs with bamCoverage --scaleFactor (mutually exclusive with --normalizeUsing)"
- "Run CUT&RUN pipeline: bowtie2 with Henikoff parameters, SEACR + MACS2 consensus, E. coli spike-in fraction"
- "Run CUT&Tag with --keep-dup all (PCR duplicates contain biology at low cycle counts)"
- "Learn 15-state ChromHMM model from H3K4me3 / H3K27ac / H3K4me1 / H3K36me3 / H3K27me3 + input"
- "Predict variant effects on FOXA1 binding using JASPAR 2026 BPNet model + ensemble"
- "Apply WASP filter then BaalChIP for allele-specific binding in cancer (CN-aware)"
- "Validate chromBPNet variant predictions against measured BaalChIP ASB"
- "Apply the ENCODE Nself/Nt consistency rule to decide whether the library passes"

## Requirements

```bash
# CLI tools
conda install -c bioconda macs3 macs2 homer samtools bedtools idr phantompeakqualtools \
    deeptools meme bowtie2 cutadapt featurecounts pygenometracks

# CUT&RUN/CUT&Tag
git clone https://github.com/FredHutch/SEACR.git    # SEACR_1.4+

# Spike-in
# ChIPseqSpikeInFree:  BiocManager::install('ChIPseqSpikeInFree')
# SpikeFlow (Snakemake wrapper): https://github.com/sebastian-gregoricchio/SpikeFlow

# Chromatin state segmentation
wget http://compbio.mit.edu/ChromHMM/ChromHMM.zip   # Java 8+
unzip ChromHMM.zip

# Deep learning (Python; GPU recommended)
pip install chrombpnet bpnet enformer-pytorch tfmodisco-lite shap

# Allele-specific binding
git clone https://github.com/bmvdgeijn/WASP.git    # WASP (mandatory upstream)
# RASQUAL: https://github.com/dg13/rasqual

# ROSE2 for super-enhancers
git clone https://github.com/linlabbcm/rose2.git
```

```r
# Bioconductor
BiocManager::install(c('DiffBind', 'DESeq2', 'edgeR', 'csaw', 'apeglm', 'rtracklayer',
                       'ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db',
                       'ChIPQC', 'rGREAT', 'chipenrich', 'monaLisa', 'BSgenome.Hsapiens.UCSC.hg38',
                       'JASPAR2024', 'BaalChIP', 'AllelicImbalance',
                       'GenomicFeatures', 'GenomicRanges', 'normr',
                       'Gviz', 'EnrichedHeatmap', 'ComplexHeatmap',
                       'ChIPseqSpikeInFree'))
```

```bash
# Python (general)
pip install pandas pyranges pysam pybedtools numpy scipy matplotlib
pip install logomaker pydeseq2
```

External resources:
- [ENCODE ChIP-seq pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline) - canonical pipeline v2.1.6
- [ENCODE blacklist v2](https://github.com/Boyle-Lab/Blacklist) - Amemiya 2019
- [SCREEN: Search Candidate cis-Regulatory Elements](https://screen.encodeproject.org/) - ENCODE cCRE registry
- [JASPAR 2026](https://jaspar.genereg.net/) - PWM + deep-learning collection (1259 BPNet ChIP models)
- [Henikoff Lab CUT&RUN/CUT&Tag protocols](https://yezhengstat.github.io/CUTTag_tutorial/)
- [chromBPNet](https://github.com/kundajelab/chrombpnet)
- [Hammond Norris 2024 spike-in review](https://pmc.ncbi.nlm.nih.gov/articles/PMC12266361/)

## Related Skills

- **atac-seq** - Parallel ATAC-seq workflows (no input control; Tn5 shift; same DiffBind/csaw/IDR/spike-in patterns)
- **alignment-files** - BAM filtering, deduplication, indexing upstream of ChIP-seq
- **pathway-analysis** - Functional enrichment of peak-associated genes
- **causal-genomics** - ASB and DL variant predictions feed fine-mapping
- **genome-intervals** - BED / GTF / bigWig manipulation
- **methylation-analysis** - Combine ChIP-seq with DNA methylation
- **read-alignment** - Bowtie2 / BWA / chromap upstream
- **read-qc** - Adapter trimming (especially aggressive for CUT&Tag short fragments)
- **data-visualization** - Heatmaps, genome tracks, browser views
