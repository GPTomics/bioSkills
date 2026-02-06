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
| **database-access** | 11 | Bio.Entrez, BLAST+, SRA toolkit, UniProt API, STRINGdb | NCBI/UniProt queries, SRA downloads, BLAST, homology searches, interaction databases |
| **alignment-files** | 9 | samtools, pysam | SAM/BAM/CRAM viewing, sorting, filtering, statistics, validation |
| **variant-calling** | 13 | bcftools, cyvcf2, Manta, Delly, VEP, SnpEff | VCF/BCF calling, SVs, filtering, annotation, clinical interpretation |
| **alignment** | 4 | Bio.Align, Bio.AlignIO | Pairwise and multiple sequence alignment, MSA statistics, alignment I/O |
| **phylogenetics** | 5 | Bio.Phylo, IQ-TREE2, RAxML-ng | Tree I/O, visualization, ML inference with model selection, ultrafast bootstrap |
| **differential-expression** | 6 | DESeq2, edgeR, ggplot2, pheatmap | RNA-seq differential expression, visualization, batch correction |
| **structural-biology** | 6 | Bio.PDB, ESMFold, Chai-1 | PDB/mmCIF parsing, SMCRA navigation, geometric analysis, ML structure prediction |
| **single-cell** | 14 | Seurat, Scanpy, Pertpy, Cassiopeia, MeboCost | scRNA-seq QC, clustering, trajectory, communication, annotation, perturb-seq, lineage tracing, metabolite communication |
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
| **alternative-splicing** | 6 | rMATS-turbo, SUPPA2, IsoformSwitchAnalyzeR | Splicing quantification, differential splicing, isoform switching, sashimi visualization |
| **chemoinformatics** | 7 | RDKit, DeepChem, AutoDock Vina | Molecular I/O, descriptors, similarity, ADMET, virtual screening, reaction enumeration |
| **liquid-biopsy** | 6 | ichorCNA, fgbio, VarDict, FinaleToolkit | cfDNA preprocessing, fragmentomics, tumor fraction, ctDNA mutations, longitudinal monitoring |
| **workflows** | 38 | Various (workflow-specific) | End-to-end pipelines: RNA-seq, variants, ChIP-seq, scRNA-seq, spatial, Hi-C, proteomics, microbiome, CRISPR, metabolomics, multi-omics, immunotherapy, outbreak, metabolic modeling, splicing, liquid biopsy, genome annotation, GRN, causal genomics |
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
| **data-visualization** | 12 | ggplot2, matplotlib, plotly, ComplexHeatmap, NetworkX | Publication-quality figures, heatmaps, interactive plots, genome tracks, circos, UpSet, volcano, networks |
| **tcr-bcr-analysis** | 5 | MiXCR, VDJtools, Immcantation, scirpy | TCR/BCR repertoire analysis, clonotype assembly, diversity metrics |
| **small-rna-seq** | 5 | miRDeep2, miRge3, cutadapt, DESeq2 | miRNA/piRNA analysis, differential expression, target prediction |
| **ribo-seq** | 5 | Plastid, RiboCode, ORFik, riborex | Ribosome profiling, translation efficiency, ORF detection |
| **epitranscriptomics** | 5 | exomePeak2, MACS3, m6Anet, Guitar | RNA modifications (m6A), MeRIP-seq, ONT direct RNA |
| **clip-seq** | 5 | CLIPper, PureCLIP, umi_tools, HOMER | Protein-RNA interactions, crosslink detection, binding site motifs |
| **clinical-databases** | 10 | myvariant, requests, pandas, SigProfiler | Clinical variant queries, ClinVar/gnomAD, pharmacogenomics, TMB, HLA, PRS, signatures |
| **genome-engineering** | 5 | crisprscan, Cas-OFFinder, PrimeDesign | CRISPR guide design, off-target prediction, prime/base editing, HDR templates |
| **systems-biology** | 5 | cobrapy, CarveMe, memote | Flux balance analysis, metabolic reconstruction, model curation, gene essentiality |
| **epidemiological-genomics** | 5 | mlst, TreeTime, TransPhylo, AMRFinderPlus | Pathogen typing, phylodynamics, transmission networks, AMR surveillance |
| **immunoinformatics** | 5 | mhcflurry, pVACtools, BepiPred | MHC binding prediction, neoantigen identification, epitope prediction |
| **comparative-genomics** | 5 | MCScanX, PAML, OrthoFinder | Synteny analysis, positive selection, ancestral reconstruction, ortholog inference |
| **genome-annotation** | 6 | Bakta, BRAKER3, eggNOG-mapper, RepeatMasker | Prokaryotic/eukaryotic annotation, functional assignment, repeats, ncRNA, annotation transfer |
| **gene-regulatory-networks** | 5 | pySCENIC, SCENIC+, WGCNA, CellOracle | Co-expression networks, regulon inference, multiomics GRN, perturbation simulation |
| **causal-genomics** | 5 | TwoSampleMR, coloc, susieR, MR-PRESSO | Mendelian randomization, colocalization, fine-mapping, mediation, pleiotropy |
| **rna-structure** | 3 | ViennaRNA, Infernal, ShapeMapper2 | RNA secondary structure prediction, ncRNA search, structure probing |
| **machine-learning** | 6 | sklearn, shap, lifelines, scvi-tools | Biomarker discovery, model interpretation, survival analysis, atlas mapping |

**Total: 412 skills across 60 categories**

## Example Usage

Once skills are deployed, ask your agent naturally. Here are examples across common workflows—the full collection covers 412 skills across 60 categories:

```
# RNA-seq & Differential Expression
"I have RNA-seq counts from treated vs control samples - find the differentially expressed genes"
"Run the complete RNA-seq pipeline from my FASTQ files to a list of DE genes"
"What biological pathways are enriched in my upregulated genes?"
"Run GSEA to see if whole pathways are up or down in my treatment"
"Align my paired-end RNA-seq reads to the human genome with STAR"
"Count reads per gene from my aligned BAM files"

# Single-Cell Analysis
"I just got my 10X scRNA-seq data - filter out low-quality cells and normalize"
"Cluster my single-cell data and help me figure out what cell types they are"
"Find marker genes for each cluster so I can annotate cell types"
"Reconstruct the differentiation trajectory and find branch points in my data"
"Which ligand-receptor pairs show active communication between my cell types?"

# Variant Calling & Clinical Genomics
"Call somatic variants from my tumor-normal BAM files"
"I found a BRCA1 variant in my patient - is it pathogenic according to ACMG guidelines?"
"Which of my variants are already known to be disease-causing in ClinVar?"
"What's the population frequency of this variant in gnomAD?"
"Annotate my VCF with gene names, functional effects, and clinical databases"
"Find structural variants like deletions and duplications in my WGS data"
"My patient has CYP2D6 variants - what's their metabolizer phenotype?"

# Epigenomics & Chromatin
"Call peaks from my ChIP-seq data for this transcription factor"
"Identify open chromatin regions from my ATAC-seq samples"
"Find differentially methylated regions between tumor and normal"
"Identify TADs and chromatin loops from my Hi-C contact matrix"

# Experimental Design & QC
"How many replicates do I need to detect 2-fold changes with 80% power?"
"Check the quality of my sequencing data before I start analysis"
"What's the alignment rate and coverage quality in my BAM files?"
"Generate a MultiQC report summarizing all my pipeline outputs"

# Protein & Structure
"Predict the 3D structure of my protein sequence using AlphaFold"
"Find differentially abundant proteins between my treatment conditions"

# CRISPR & Genome Engineering
"Design guide RNAs to knock out TP53 with minimal off-targets"
"Analyze my CRISPR dropout screen to find essential genes"

# Pipelines & Reproducibility
"Set up a Snakemake workflow so I can rerun this analysis on 50 samples"
"Run a complete biomarker discovery pipeline with proper cross-validation"
"Annotate my new genome assembly end-to-end from repeats to functional annotation"
"Run post-GWAS causal inference on my summary statistics"

# Sequences & Databases
"Download gene sequences and annotations from NCBI for my gene list"
"Design PCR primers to amplify a 500bp region of this gene"
"Read my FASTA files and extract specific sequences"

# Spatial & Tissue Analysis
"Identify spatially distinct tissue regions in my Visium data"
"What species are present in my metagenomic sample?"

# Long-Read Sequencing
"Basecall my Oxford Nanopore fast5 files with high accuracy"
"Assess if my genome assembly is complete and high quality"

# Specialized Analysis
"Analyze differential exon usage to detect alternative splicing changes"
"Extract and analyze TCR sequences from my T cell RNA-seq"
"Build a survival model to find genes associated with patient outcomes"
"Use machine learning to discover biomarkers that predict treatment response"
"Validate my predictive model with proper cross-validation to avoid overfitting"

# Immunotherapy & Cancer
"Predict neoantigens from tumor mutations for immunotherapy"
"Determine HLA type from RNA-seq for neoantigen prediction"
"Detect ultra-low frequency mutations in my liquid biopsy cfDNA"

# Genome Annotation
"Annotate my newly assembled bacterial genome with Bakta"
"Run BRAKER3 gene prediction on my eukaryotic assembly"
"Assign functional annotations with eggNOG-mapper and InterProScan"

# Gene Regulatory Networks
"Infer transcription factor regulons from my single-cell data with pySCENIC"
"Build a co-expression network and find hub genes with WGCNA"
"Simulate what happens if I knock out this transcription factor"

# Causal Genomics
"Run Mendelian randomization to test if BMI causes heart disease"
"Test whether my GWAS hit and eQTL share the same causal variant"
"Fine-map my GWAS locus to identify the most likely causal variants"

# RNA Structure
"Predict the secondary structure and folding energy of my RNA sequence"
"Search for ncRNA homologs in my transcript using Rfam"

# Phylogenetics & Evolution
"Build a phylogenetic tree and visualize evolutionary relationships"
"Find orthologs of my human gene across vertebrate species"
```

The agent will select appropriate tools based on context. See the Skill Categories table above for the complete list of available skills.

## Contributing

See `CLAUDE.md` for development guidelines, file structure requirements, and quality standards.

Key requirements:
- SKILL.md must include "Use when..." in description
- `primary_tool` must be a single value (not comma-separated)
- Quick Start uses bullets; Example Prompts use blockquotes
- Examples must document magic numbers with rationale

## License

MIT License - see LICENSE file for details.
