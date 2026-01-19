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

# conda (recommended for full bioinformatics stack)
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
```

### Codex CLI

```bash
./install-codex.sh                               # Install globally
./install-codex.sh --project /path/to/project    # Or install to specific project
```

### Gemini CLI

```bash
./install-gemini.sh                              # Install globally
./install-gemini.sh --project /path/to/project   # Or install to specific project
```

Codex and Gemini installers convert to the Agent Skills standard (`examples/` -> `scripts/`, `usage-guide.md` -> `references/`).

## Skill Categories

| Category | Skills | Primary Tools | Description |
|----------|--------|---------------|-------------|
| **sequence-io** | 9 | Bio.SeqIO | Read, write, convert FASTA/FASTQ/GenBank and 40+ formats |
| **sequence-manipulation** | 7 | Bio.Seq | Transcription, translation, motif search, sequence properties |
| **database-access** | 10 | Bio.Entrez, BLAST+, HMMER, SRA toolkit | NCBI/UniProt queries, SRA downloads, BLAST, homology searches |
| **alignment-files** | 9 | samtools, pysam | SAM/BAM/CRAM viewing, sorting, filtering, validation |
| **variant-calling** | 11 | bcftools, GATK, Manta, Delly, VEP | VCF/BCF variant calling, SVs, filtering, annotation |
| **alignment** | 4 | Bio.Align | Pairwise and multiple sequence alignment |
| **phylogenetics** | 4 | Bio.Phylo | Tree I/O, visualization, manipulation, distance matrices |
| **differential-expression** | 4 | DESeq2, edgeR | RNA-seq differential expression analysis |
| **structural-biology** | 5 | Bio.PDB | PDB/mmCIF parsing, geometric analysis, AlphaFold |
| **single-cell** | 7 | Seurat, Scanpy, Signac | scRNA-seq QC, doublets, clustering, markers, multimodal, scATAC |
| **pathway-analysis** | 6 | clusterProfiler | GO, KEGG, Reactome, WikiPathways enrichment |
| **restriction-analysis** | 4 | Bio.Restriction | Restriction sites, mapping, enzyme selection |
| **methylation-analysis** | 4 | Bismark, methylKit | Bisulfite alignment, methylation calling, DMRs |
| **chip-seq** | 6 | MACS3, ChIPseeker, HOMER, IDR | Peak calling, annotation, differential binding, motifs, QC |
| **metagenomics** | 5 | Kraken2, MetaPhlAn | Taxonomic classification, abundance estimation |
| **long-read-sequencing** | 5 | Dorado, minimap2, medaka | Basecalling, alignment, polishing, SV calling |
| **read-qc** | 7 | FastQC, fastp, umi_tools, RSeQC | Quality reports, adapter trimming, filtering, UMIs, RNA-seq QC |
| **genome-intervals** | 7 | BEDTools, pybedtools, pyBigWig | BED/GTF operations, interval arithmetic, bedGraph, bigWig |
| **population-genetics** | 6 | PLINK, scikit-allel | GWAS, population structure, selection statistics |
| **rna-quantification** | 4 | Salmon, featureCounts | Gene/transcript quantification, count matrix QC |
| **read-alignment** | 4 | bwa-mem2, STAR | Short-read alignment for DNA and RNA-seq |
| **expression-matrix** | 4 | pandas, anndata | Count matrix handling, gene ID mapping |
| **copy-number** | 4 | CNVkit, GATK | CNV detection, visualization, annotation |
| **phasing-imputation** | 4 | Beagle, SHAPEIT5 | Haplotype phasing, genotype imputation |
| **atac-seq** | 4 | MACS3, TOBIAS | ATAC-seq peaks, QC, footprinting |
| **genome-assembly** | 4 | SPAdes, Flye, QUAST | De novo assembly, polishing, QC |

**Total: 148 skills across 26 categories**

## Usage

Once skills are deployed, ask your agent naturally:

```
"Parse my FASTA file and show sequence lengths"
"Filter FASTQ reads with mean quality below 25"
"Translate this coding sequence to protein"
"Sort my BAM file by coordinate and index it"
"Call variants from this BAM against the reference"
"Search NCBI for all human BRCA1 sequences"
"Download the SRA run SRR1234567"
"Align these two protein sequences using BLOSUM62"
"Calculate the pairwise identity matrix for this alignment"
"Convert my Clustal alignment to PHYLIP format"
"Build a neighbor joining tree from this alignment"
"Root the tree using Mouse as outgroup"
"Draw this phylogenetic tree and save as PDF"
"Run DESeq2 on my count matrix comparing treated vs control"
"Create a volcano plot of my differential expression results"
"Extract genes with padj < 0.05 and |log2FC| > 1"
"Download PDB structure 4HHB and list its chains"
"Calculate the distance between residue 50 and 100 CA atoms"
"Superimpose these two structures and calculate RMSD"
"Remove all water molecules from this PDB file"
"Load my 10X data and run QC filtering"
"Cluster my single-cell data and generate a UMAP"
"Find marker genes for each cluster"
"Annotate cell types based on canonical markers"
"Run GO enrichment on my differentially expressed genes"
"Find enriched KEGG pathways for my gene list"
"Perform GSEA using MSigDB hallmark gene sets"
"Create a dot plot of my enrichment results"
"Find enzymes that cut my plasmid exactly once"
"What fragments will EcoRI produce from this sequence?"
"Find enzymes that don't cut my insert but cut the vector"
"Create a restriction map for my sequence"
"Align my bisulfite sequencing reads with Bismark"
"Find differentially methylated regions between treatment and control"
"Call peaks from my ChIP-seq BAM with MACS3"
"Annotate my peaks with nearby genes"
"Find differential binding between conditions"
"Classify my metagenomic reads with Kraken2"
"Profile my metagenome with MetaPhlAn"
"Create a stacked bar chart of species abundances"
"Align my Nanopore reads with minimap2"
"Polish my assembly with medaka"
"Find structural variants from my long reads"
"Generate QC report for my long read data"
"Run FastQC on my FASTQ files and summarize with MultiQC"
"Trim adapters from my paired-end reads with Cutadapt"
"Run fastp on my RNA-seq data"
"Find peaks that overlap with promoters"
"Convert bedGraph to bigWig"
"Extract gene coordinates from GTF"
"Convert VCF to PLINK binary format"
"Run GWAS on my case-control data"
"Perform PCA for population structure"
"Calculate Fst between populations"
"Count reads per gene from my BAM files"
"Quantify transcripts with Salmon"
"Import Salmon results into R for DESeq2"
"Check for sample outliers before DE analysis"
"Align my DNA reads to the human genome with bwa-mem2"
"Run STAR on my RNA-seq data"
"Build a STAR index with my GTF file"
"Load my featureCounts output into a dataframe"
"Map Ensembl IDs to gene symbols"
"Join sample metadata with my count matrix"
"Detect CNVs from my exome BAM files"
"Run CNVkit on tumor-normal pairs"
"Plot copy number profile for a sample"
"Phase my VCF file with Beagle"
"Impute missing genotypes using 1000 Genomes"
"Call peaks from my ATAC-seq BAM files"
"Calculate TSS enrichment score"
"Run TF footprinting with TOBIAS"
"Call structural variants with Manta"
"Detect deletions and inversions with Delly"
"Extract UMIs and deduplicate my reads"
"Annotate variants with VEP"
"Find enriched motifs in my ChIP-seq peaks"
"Process my scATAC-seq data with Signac"
"Assemble my bacterial genome with SPAdes"
"Polish my assembly with Pilon"
"Check assembly completeness with BUSCO"
```

The agent will use the skill patterns to generate correct code.

