# bioSkills

A collection of skills that guide AI coding agents (Claude Code, Codex, Gemini) through common bioinformatics tasks.

## Project Goal

This repository provides AI agents with expert knowledge for bioinformatics workflows. Each skill contains code patterns, best practices, and examples that help agents generate correct, idiomatic code for common tasks.

Target users range from undergrads learning computational biology to PhD researchers processing large-scale data. The skills cover the full spectrum from basic sequence manipulation to advanced analyses like single-cell RNA-seq and population genetics.

## Requirements

### Python
- Python 3.9+
- biopython, pysam, cyvcf2, pybedtools, pyBigWig, scikit-allel, anndata

```bash
pip install biopython pysam cyvcf2 pybedtools pyBigWig scikit-allel anndata mygene
```

### R/Bioconductor
Required for differential expression, single-cell, pathway analysis, and methylation skills.

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'edgeR', 'Seurat', 'clusterProfiler', 'methylKit'))
```

### CLI Tools
```bash
# macOS
brew install samtools bcftools blast minimap2 bedtools

# Ubuntu/Debian
sudo apt install samtools bcftools ncbi-blast+ minimap2 bedtools

# conda
conda install -c bioconda samtools bcftools blast minimap2 bedtools \
    fastp kraken2 metaphlan sra-tools bwa-mem2 bowtie2 star hisat2 \
    manta delly cnvkit macs3 tobias
```

## Installation

### Claude Code

```bash
git clone https://github.com/your-username/bioSkills.git
cd bioSkills
./install-claude.sh                              # Install globally
./install-claude.sh --project /path/to/project   # Or install to specific project
./install-claude.sh --list                       # List available skills
./install-claude.sh --validate                   # Validate all skills
./install-claude.sh --update                     # Only update changed skills
./install-claude.sh --uninstall                  # Remove all bio-* skills
```

### Codex CLI

```bash
./install-codex.sh                               # Install globally
./install-codex.sh --project /path/to/project    # Or install to specific project
./install-codex.sh --validate                    # Validate all skills
./install-codex.sh --update                      # Only update changed skills
./install-codex.sh --uninstall                   # Remove all bio-* skills
```

### Gemini CLI

```bash
./install-gemini.sh                              # Install globally
./install-gemini.sh --project /path/to/project   # Or install to specific project
./install-gemini.sh --validate                   # Validate all skills
./install-gemini.sh --update                     # Only update changed skills
./install-gemini.sh --uninstall                  # Remove all bio-* skills
```

Codex and Gemini installers convert to the Agent Skills standard (`examples/` -> `scripts/`, `usage-guide.md` -> `references/`).

## Skill Categories

| Category | Skills | Primary Tools | Description |
|----------|--------|---------------|-------------|
| **sequence-io** | 9 | Bio.SeqIO | Read, write, convert FASTA/FASTQ/GenBank and 40+ formats |
| **sequence-manipulation** | 7 | Bio.Seq, Bio.SeqUtils | Transcription, translation, motif search, sequence properties |
| **database-access** | 10 | Bio.Entrez, BLAST+, SRA toolkit, UniProt API | NCBI/UniProt queries, SRA downloads, BLAST, homology searches |
| **alignment-files** | 9 | samtools, pysam | SAM/BAM/CRAM viewing, sorting, filtering, statistics, validation |
| **variant-calling** | 13 | bcftools, cyvcf2, Manta, Delly, VEP, SnpEff | VCF/BCF calling, SVs, filtering, annotation, clinical interpretation |
| **alignment** | 4 | Bio.Align, Bio.AlignIO | Pairwise and multiple sequence alignment, MSA statistics, alignment I/O |
| **phylogenetics** | 5 | Bio.Phylo, IQ-TREE2, RAxML-ng | Tree I/O, visualization, ML inference with model selection, ultrafast bootstrap |
| **differential-expression** | 6 | DESeq2, edgeR, ggplot2, pheatmap | RNA-seq differential expression, visualization, batch correction |
| **structural-biology** | 6 | Bio.PDB, ESMFold, Chai-1 | PDB/mmCIF parsing, SMCRA navigation, geometric analysis, ML structure prediction |
| **single-cell** | 13 | Seurat, Scanpy, Pertpy, Cassiopeia | scRNA-seq QC, clustering, trajectory, communication, annotation, perturb-seq, lineage tracing |
| **pathway-analysis** | 6 | clusterProfiler, ReactomePA, rWikiPathways, enrichplot | GO, KEGG, Reactome, WikiPathways enrichment |
| **restriction-analysis** | 4 | Bio.Restriction | Restriction sites, mapping, enzyme selection |
| **methylation-analysis** | 4 | Bismark, methylKit, bsseq | Bisulfite alignment, methylation calling, DMRs |
| **chip-seq** | 7 | MACS3, ChIPseeker, DiffBind | Peak calling, annotation, differential binding, motifs, QC, super-enhancers |
| **metagenomics** | 7 | Kraken2, MetaPhlAn, Bracken, HUMAnN | Taxonomic classification, abundance estimation, functional profiling, AMR detection |
| **long-read-sequencing** | 8 | Dorado, minimap2, Clair3, modkit, IsoSeq3 | Basecalling, alignment, polishing, variant calling, SV calling, methylation, Iso-Seq |
| **read-qc** | 7 | FastQC, MultiQC, fastp, Trimmomatic, Cutadapt | Quality reports, adapter trimming, filtering, UMIs |
| **genome-intervals** | 7 | BEDTools, pybedtools, pyBigWig | BED/GTF operations, interval arithmetic, bedGraph, bigWig |
| **population-genetics** | 6 | PLINK, FlashPCA2, ADMIXTURE, scikit-allel | GWAS, biobank-scale PCA, admixture, selection statistics |
| **rna-quantification** | 4 | featureCounts, Salmon, kallisto, tximport | Gene/transcript quantification, count matrix QC |
| **read-alignment** | 4 | bwa-mem2, bowtie2, STAR, HISAT2 | Short-read alignment for DNA and RNA-seq |
| **expression-matrix** | 4 | pandas, anndata, scanpy, biomaRt | Count matrix handling, gene ID mapping |
| **copy-number** | 4 | CNVkit, GATK | CNV detection, visualization, annotation |
| **phasing-imputation** | 4 | Beagle, SHAPEIT5, bcftools | Haplotype phasing, genotype imputation |
| **atac-seq** | 6 | MACS3, DiffBind, chromVAR, TOBIAS | ATAC-seq peaks, differential accessibility, footprinting, TF motif deviation |
| **genome-assembly** | 8 | SPAdes, Flye, hifiasm, QUAST, BUSCO | Assembly, polishing, scaffolding, quality assessment |
| **primer-design** | 3 | primer3-py | PCR primer design, qPCR probes, validation |
| **spatial-transcriptomics** | 11 | Squidpy, SpatialData, Scanpy, scimap | Visium, Xenium, Slide-seq, spatial stats, domain detection, deconvolution, spatial proteomics |
| **hi-c-analysis** | 8 | cooler, cooltools, pairtools, HiCExplorer | Contact matrices, compartments, TADs, loops, differential |
| **workflows** | 28 | Various (workflow-specific) | End-to-end pipelines: RNA-seq, variants, ChIP-seq, scRNA-seq, spatial, Hi-C, proteomics, microbiome, CRISPR, metabolomics, multi-omics |
| **proteomics** | 9 | pyOpenMS, MSstats, limma, QFeatures | Mass spec data import, QC, quantification, differential abundance, PTM, DIA |
| **microbiome** | 6 | DADA2, phyloseq, ALDEx2, QIIME2 | 16S/ITS amplicon processing, taxonomy, diversity, differential abundance |
| **multi-omics-integration** | 4 | MOFA2, mixOmics, SNF | Cross-modality integration, factor analysis, network fusion |
| **crispr-screens** | 8 | MAGeCK, JACKS, CRISPResso2, BAGEL2 | Pooled screen analysis, sgRNA efficacy modeling, hit calling, base/prime editing |
| **metabolomics** | 8 | XCMS, MetaboAnalystR, lipidr, MS-DIAL | Peak detection, annotation, normalization, pathway mapping, lipidomics, targeted |
| **imaging-mass-cytometry** | 6 | steinbock, squidpy, napari | IMC preprocessing, segmentation, spatial analysis, annotation, QC |
| **flow-cytometry** | 8 | flowCore, CATALYST, CytoML | FCS handling, compensation, gating, clustering, differential, QC |
| **reporting** | 5 | RMarkdown, Quarto, Jupyter, MultiQC, matplotlib | Reproducible reports, QC aggregation, publication figures |
| **experimental-design** | 4 | RNASeqPower, ssizeRNA, qvalue, sva | Power analysis, sample size, multiple testing, batch design |
| **workflow-management** | 4 | Snakemake, Nextflow, cwltool, Cromwell | Scalable pipeline frameworks with containers |
| **data-visualization** | 11 | ggplot2, matplotlib, plotly, ComplexHeatmap | Publication-quality figures, heatmaps, interactive plots, genome tracks, circos, UpSet, volcano |
| **tcr-bcr-analysis** | 5 | MiXCR, VDJtools, Immcantation, scirpy | TCR/BCR repertoire analysis, clonotype assembly, diversity metrics |
| **small-rna-seq** | 5 | miRDeep2, miRge3, cutadapt, DESeq2 | miRNA/piRNA analysis, differential expression, target prediction |
| **ribo-seq** | 5 | Plastid, RiboCode, ORFik, riborex | Ribosome profiling, translation efficiency, ORF detection |
| **epitranscriptomics** | 5 | exomePeak2, MACS3, m6Anet, Guitar | RNA modifications (m6A), MeRIP-seq, ONT direct RNA |
| **clip-seq** | 5 | CLIPper, PureCLIP, umi_tools, HOMER | Protein-RNA interactions, crosslink detection, binding site motifs |
| **clinical-databases** | 10 | myvariant, requests, pandas, SigProfiler | Clinical variant queries, ClinVar/gnomAD, pharmacogenomics, TMB, HLA, PRS, signatures |

**Total: 330 skills across 47 categories**

## Example Usage

Once skills are deployed, ask your agent naturally:

```
"QC my paired-end FASTQ files and trim adapters"
"Align my RNA-seq reads to the human genome"
"Call variants from my whole-exome BAM files"
"Find structural variants in my tumor-normal pair"
"Annotate my VCF with clinical significance scores"
"Run differential expression on my count matrix, treated vs control"
"Analyze my time-series RNA-seq experiment"
"What pathways are enriched in my upregulated genes?"
"Load my 10X data, remove doublets, and cluster"
"Integrate these three scRNA-seq batches and find markers"
"Infer developmental trajectories from my single-cell data"
"What ligand-receptor pairs are active between my cell types?"
"Annotate cell types in my PBMC dataset"
"Call peaks from my ATAC-seq and check TSS enrichment"
"Find transcription factor footprints in my accessibility data"
"Run chromVAR to find variable TF motifs across my ATAC-seq samples"
"Compare ChIP-seq binding between treatment and control"
"Identify super-enhancers from my H3K27ac ChIP-seq"
"Find A/B compartments and TADs in my Hi-C data"
"Detect chromatin loops from my contact matrix"
"What species are in my metagenomic sample?"
"Process my 16S amplicon data and assign taxonomy"
"Find differentially abundant taxa between my sample groups"
"Screen my metagenome for antibiotic resistance genes"
"Run a GWAS on my case-control genotypes"
"Calculate population differentiation statistics"
"Call variants from my Oxford Nanopore data"
"Find structural variants in my long-read data"
"Assemble my PacBio HiFi reads and check completeness"
"Scaffold my draft assembly using Hi-C"
"Quantify proteins from my DIA mass spec data"
"Find differentially abundant proteins between conditions"
"Process my LC-MS metabolomics data and normalize"
"Annotate my metabolite features with pathway information"
"Analyze my CRISPR dropout screen and find hits"
"Jointly analyze multiple CRISPR screens with JACKS for sgRNA efficacy"
"Find spatially variable genes in my Visium data"
"Deconvolve cell types in my spatial transcriptomics spots"
"Cluster my CyTOF data and compare populations between groups"
"Segment cells from my imaging mass cytometry and analyze neighborhoods"
"Integrate my transcriptomics and proteomics datasets"
"Predict my protein structure with ESMFold"
"Build an ML phylogenetic tree with IQ-TREE2 and ultrafast bootstrap"
"Create a clustered heatmap of my top DE genes"
"Analyze my TCR repertoire and calculate diversity"
"Find differentially expressed miRNAs and their targets"
"Calculate translation efficiency from Ribo-seq data"
"Detect m6A modifications from MeRIP-seq"
"Find RBP binding sites from my CLIP-seq data"
"Look up clinical significance of my variants in ClinVar"
"Track clonal lineages using CRISPR barcodes"
"Analyze my Perturb-seq CRISPR screen"
"How many samples do I need for my RNA-seq experiment?"
"Calculate power for detecting 2-fold changes"
"Analyze my base editing experiment for C-to-T conversion"
"Call methylation from my nanopore BAM file"
"Generate a MultiQC report from my pipeline outputs"
"Export my figure at 300 DPI for journal submission"
```

The agent will select appropriate tools based on context.

## Contributing

See `CLAUDE.md` for development guidelines, file structure requirements, and quality standards.

Key requirements:
- SKILL.md must include "Use when..." in description
- `primary_tool` must be a single value (not comma-separated)
- Quick Start uses bullets; Example Prompts use blockquotes
- Examples must document magic numbers with rationale

## License

MIT License - see LICENSE file for details.
