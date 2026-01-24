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
| **sequence-manipulation** | 7 | Bio.Seq | Transcription, translation, motif search, sequence properties |
| **database-access** | 10 | Bio.Entrez, BLAST+, HMMER, SRA toolkit | NCBI/UniProt queries, SRA downloads, BLAST, homology searches |
| **alignment-files** | 9 | samtools, pysam | SAM/BAM/CRAM viewing, sorting, filtering, validation |
| **variant-calling** | 13 | bcftools, GATK, DeepVariant, Manta, Delly, VEP | VCF/BCF calling, SVs, filtering, annotation, clinical interpretation, joint calling |
| **alignment** | 4 | Bio.Align | Pairwise and multiple sequence alignment |
| **phylogenetics** | 5 | Bio.Phylo, IQ-TREE2 | Tree I/O, visualization, manipulation, distance matrices, ML inference |
| **differential-expression** | 6 | DESeq2, edgeR, sva, limma | RNA-seq differential expression, batch correction, time-series |
| **structural-biology** | 6 | Bio.PDB, ESMFold | PDB/mmCIF parsing, geometric analysis, AlphaFold, modern structure prediction |
| **single-cell** | 13 | Seurat, Scanpy, Signac, CellChat, Pertpy, Cassiopeia | scRNA-seq QC, clustering, trajectory, communication, annotation, integration, perturb-seq, lineage tracing |
| **pathway-analysis** | 6 | clusterProfiler | GO, KEGG, Reactome, WikiPathways enrichment |
| **restriction-analysis** | 4 | Bio.Restriction | Restriction sites, mapping, enzyme selection |
| **methylation-analysis** | 4 | Bismark, methylKit | Bisulfite alignment, methylation calling, DMRs |
| **chip-seq** | 7 | MACS3, ChIPseeker, HOMER, IDR, ROSE | Peak calling, annotation, differential binding, motifs, QC, super-enhancers |
| **metagenomics** | 7 | Kraken2, MetaPhlAn, AMRFinderPlus, MASH | Taxonomic classification, abundance estimation, AMR detection, strain tracking |
| **long-read-sequencing** | 7 | Dorado, minimap2, medaka, Clair3, IsoSeq3 | Basecalling, alignment, polishing, variant calling, SV calling, Iso-Seq |
| **read-qc** | 7 | FastQC, fastp, umi_tools, RSeQC | Quality reports, adapter trimming, filtering, UMIs, RNA-seq QC |
| **genome-intervals** | 7 | BEDTools, pybedtools, pyBigWig | BED/GTF operations, interval arithmetic, bedGraph, bigWig |
| **population-genetics** | 6 | PLINK, FlashPCA2, scikit-allel | GWAS, population structure, selection statistics |
| **rna-quantification** | 4 | Salmon, featureCounts | Gene/transcript quantification, count matrix QC |
| **read-alignment** | 4 | bwa-mem2, STAR | Short-read alignment for DNA and RNA-seq |
| **expression-matrix** | 4 | pandas, anndata | Count matrix handling, gene ID mapping |
| **copy-number** | 4 | CNVkit, GATK | CNV detection, visualization, annotation |
| **phasing-imputation** | 4 | Beagle, SHAPEIT5 | Haplotype phasing, genotype imputation |
| **atac-seq** | 6 | MACS3, TOBIAS, chromVAR | ATAC-seq peaks, QC, footprinting, nucleosome positioning, motif deviation |
| **genome-assembly** | 8 | SPAdes, Flye, hifiasm, YaHS, CheckM2 | Assembly, polishing, scaffolding, contamination detection |
| **primer-design** | 3 | primer3-py | PCR primer design, qPCR probes, validation |
| **spatial-transcriptomics** | 10 | Squidpy, SpatialData | Visium, Xenium, Slide-seq, spatial stats, domain detection, deconvolution, high-resolution |
| **hi-c-analysis** | 8 | cooler, cooltools, pairtools | Contact matrices, compartments, TADs, loops, differential |
| **workflows** | 28 | mixed | End-to-end pipelines: RNA-seq, variants, somatic, ChIP-seq, scRNA-seq, spatial, Hi-C, proteomics, microbiome, CRISPR, metabolomics, IMC, cytometry, multi-omics, TCR, small-RNA, Ribo-seq, MeRIP, CLIP |
| **proteomics** | 9 | pyOpenMS, MSstats, DIA-NN, limma | Mass spec data import, QC, quantification, differential abundance, PTM, DIA, spectral libraries |
| **microbiome** | 6 | DADA2, phyloseq, ALDEx2, QIIME2 | 16S/ITS amplicon processing, taxonomy, diversity, differential abundance |
| **multi-omics-integration** | 4 | MOFA2, mixOmics, SNF | Cross-modality integration, factor analysis, network fusion |
| **crispr-screens** | 7 | MAGeCK, JACKS, CRISPResso2 | Pooled screen analysis, guide counting, hit calling, QC, sgRNA efficacy |
| **metabolomics** | 8 | XCMS, MetaboAnalystR, lipidr, MS-DIAL | Peak detection, annotation, normalization, pathway mapping, lipidomics, targeted |
| **imaging-mass-cytometry** | 6 | steinbock, Cellpose, squidpy, napari | IMC preprocessing, segmentation, spatial analysis, annotation, QC |
| **flow-cytometry** | 8 | flowCore, CATALYST, diffcyt, flowAI | FCS handling, compensation, gating, clustering, differential, QC, bead normalization |
| **reporting** | 2 | RMarkdown, Quarto | Reproducible analysis reports in HTML, PDF, Word |
| **workflow-management** | 4 | Snakemake, Nextflow, CWL, WDL | Scalable pipeline frameworks with containers |
| **data-visualization** | 8 | ggplot2, ComplexHeatmap, plotly, pyGenomeTracks, Circos | Publication-quality figures, heatmaps, interactive plots, genome tracks, circos |
| **tcr-bcr-analysis** | 5 | MiXCR, VDJtools, Immcantation, scirpy | TCR/BCR repertoire analysis, clonotype assembly, diversity metrics |
| **small-rna-seq** | 5 | miRDeep2, miRge3, DESeq2, miRanda | miRNA/piRNA analysis, differential expression, target prediction |
| **ribo-seq** | 5 | Plastid, RiboCode, ORFik | Ribosome profiling, translation efficiency, ORF detection |
| **epitranscriptomics** | 5 | exomePeak2, m6Anet, Guitar | RNA modifications (m6A), MeRIP-seq, ONT direct RNA |
| **clip-seq** | 5 | CLIPper, PureCLIP, umi_tools | Protein-RNA interactions, binding site analysis, CLIP protocols |
| **clinical-databases** | 5 | myvariant, ClinVar API, gnomAD | Clinical variant queries, frequency databases, variant prioritization |

**Total: 312 skills across 46 categories**

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
