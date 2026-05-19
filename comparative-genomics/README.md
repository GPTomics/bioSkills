# comparative-genomics

## Overview

Comparative genomics across species for inferring evolutionary history, gene-family dynamics, structural rearrangements, polyploidy, horizontal gene transfer, introgression, positive selection, and ancestral state reconstruction. Covers prokaryote-to-vertebrate scale with explicit treatment of incomplete lineage sorting, ghost lineages, WGD, gBGC, ortholog conjecture, and contamination artifacts.

**Tool type:** mixed | **Primary tools:** OrthoFinder3, PAML, HyPhy, MCScanX, JCVI, GENESPACE, SyRI, Cactus, Minigraph-Cactus, PGGB, PGR-TK, wgd v2, KsRates, ALE, GeneRax, AleRax, Whale.jl, Dsuite, TreeMix, qpAdm, AdmixTools v2, HGTector v2, AvP, TOGA, CESAR 2.0, LiftOff, CAFE5, Count, Panaroo, PPanGGOLiN, PEPPAN, skani, FastANI, GTDB-Tk

## Skills

| Skill | Description |
|-------|-------------|
| ancestral-reconstruction | Sequence, discrete trait, continuous trait, and gene content ASR with PAML codeml, IQ-TREE2 --ancestral, GRASP indel-aware (protein resurrection), corHMM hidden Markov, phytools::make.simmap stochastic mapping, BM/OU/EB continuous-trait models, HiSSE state-dependent diversification |
| hgt-detection | HGT via composition (GC%, codon, tetranucleotide), phylogenetic incongruence (ALE, GeneRax, AleRax, RANGER-DTL), BLAST distribution (HGTector v2, DarkHorse, Alien Index, AvP); contamination filtering (FCS-GX, BlobTools); Lawrence-Ochman amelioration; DGL discrimination |
| ortholog-inference | OrthoFinder3 HOGs (Quest-for-Orthologs leader 2025), SonicParanoid2, Broccoli, ProteinOrtho, OMA / FastOMA; ortholog conjecture caveat; synteny-aware in WGD lineages; TOGA WGA-anchored orthology + gene-loss |
| positive-selection | PAML codeml site/branch/branch-site M0-M8 + HyPhy BUSTED / BUSTED-S / BUSTED-MH / MEME / FEL / FUBAR / aBSREL / RELAX; GARD recombination pre-screen; PREQUAL/HmmCleaner alignment filter; gBGC confounder; asymptotic alpha + polyDFE + GRAPES; RERconverge / CSUBST / PhyloAcc for convergence |
| synteny-analysis | MCScanX, JCVI/MCScan Python, GENESPACE (riparian plots, pan-gene tracks), SyRI + plotsr (structural rearrangements), AnchorWave (WGD-aware), i-ADHoRe (deep), ntSynt (alignment-free macrosynteny); microsynteny vs macrosynteny; repeat-masking critical |
| gene-tree-species-tree-reconciliation | ALE / GeneRax / AleRax probabilistic DTL; Whale.jl Bayesian DTL+WGD; RANGER-DTL / NOTUNG parsimony; ALE-rooting for deep species trees (Williams 2017); ILS vs DTL via ecceTERA / DLCpar |
| whole-genome-alignment | Progressive Cactus (reference-free clade-scale), Minigraph-Cactus (HPRC pangenome), LASTZ chain/net (UCSC), MUMmer4, minimap2 -x asm5/10/20, AnchorWave, Winnowmap2 (centromeres); HAL toolkit for downstream extraction |
| whole-genome-duplication | wgd v2 + KsRates (rate-corrected Ks dating), DupGen_finder, MAPS, POInT, SLEDGe (ML classifier), Whale.jl (Bayesian DTL+WGD); 2R / 3R / Ss4R vertebrate WGDs; angiosperm zeta |
| pangenome-analysis | Bacterial: Panaroo + PPanGGOLiN + PEPPAN (replaces deprecated Roary); Eukaryotic: Minigraph-Cactus + PGGB + vg + PanGenie; PGR-TK / PANGEA (DGI / Diploid Genomics) for repetitive / clinical genes (MHC, DAZ, OPN); pan-GWAS with Scoary / pyseer |
| genome-distance-and-species-delineation | skani (default in GTDB-Tk 2.4+), FastANI, pyani, OrthoANI; GTDB-Tk r220 taxonomy; TYGS / GGDC dDDH for novel-species naming; 95% ANI + AF >= 0.5 species threshold (Jain 2018); AAI for distant comparison |
| introgression-detection | Dsuite (Patterson D + Fbranch + ABBAclustering 2024), AdmixTools v2 (qpAdm / qpGraph), TreeMix + OptM, QuIBL / Twisst (locus-level topology), PhyloNet (explicit networks), sprime (archaic), HyDe; ILS vs introgression vs ghost lineage |
| gene-family-evolution | CAFE5 (gamma rate categories), Count (Csurös ASR), BadiRate; annotation-pipeline normalization mandatory; HGT-affected families -> ALE / GeneRax / AleRax; convergent rate shifts via RERconverge |
| comparative-annotation-projection | TOGA + CESAR 2.0 (Kirilenko 2023; Zoonomia / Bird10000 standard); LiftOff (pairwise); UCSC liftOver (coordinate); Comparative Annotation Toolkit (multi-reference); GeMoMa (evidence integration) |

## Example Prompts

- "Build Progressive Cactus alignment of 30 mammal genomes with TOGA annotation projection across all"
- "Identify whole-genome duplication events in the soybean lineage using wgd v2 + KsRates with multiple outgroups"
- "Run ALE undated to identify horizontal gene transfers across 100 bacterial genomes and report per-branch transfer counts"
- "Classify these MAGs to GTDB taxonomy using GTDB-Tk + skani with 95% ANI + AF >= 0.5 species delineation"
- "Test for introgression between Heliconius species using Dsuite Dtrios + Fbranch + Twisst window-level confirmation"
- "Detect positive selection in immune genes across primates with GARD pre-screen + BUSTED-MH + MEME, cross-validated with PAML M8 vs M8a"
- "Reconstruct ancestral cytochrome c with GRASP for protein resurrection, design 8 alternative constructs at ambiguous sites"
- "Build pangenome of 200 E. coli strains with Panaroo + PPanGGOLiN, compute Heaps law for openness, run Scoary pan-GWAS for antibiotic resistance"
- "Identify gene families with lineage-specific expansion in cetaceans using CAFE5 + RERconverge convergent rate shifts"
- "Apply PGR-TK to MHC class II locus across HPRC haplotypes; check upstream for current PANGEA (DGI) repo pointer"

## Requirements

```bash
# Phylogenetic / sequence ASR
conda install -c bioconda paml iqtree fastml prank macse

# Orthology
conda install -c bioconda orthofinder=3.0 sonicparanoid broccoli proteinortho diamond mmseqs2

# Positive selection
conda install -c bioconda hyphy gard prequal hmmcleaner

# Synteny + WGD
conda install -c bioconda mcscanx anchorwave syri plotsr mummer4 minimap2
pip install jcvi wgd ksrates

# Whole-genome alignment
docker pull quay.io/comparative-genomics-toolkit/cactus:latest
conda install -c bioconda lastz mash pggb vg pangenie

# Pangenome
conda install -c bioconda panaroo ppanggolin peppan anvio bakta
conda install -c bioconda pgr-tk

# Distance / taxonomy
conda install -c bioconda skani fastani gtdbtk checkm2 mash dashing2

# Reconciliation
conda install -c bioconda ale generax

# Introgression
conda install -c bioconda treemix sprime
git clone https://github.com/millanek/Dsuite && cd Dsuite && make
Rscript -e "remotes::install_github('uqrmaie1/admixtools')"

# Gene family
conda install -c bioconda cafe treepl

# Annotation projection
pip install liftoff
conda install -c bioconda ucsc-liftover
# TOGA: conda env from hillerlab/TOGA

# R packages
Rscript -e "install.packages(c('ape', 'phytools', 'geiger', 'corHMM', 'OUwie', 'mclust'))"
Rscript -e "remotes::install_github(c('jtlovell/GENESPACE', 'nclark-lab/RERconverge', 'thej022214/hisse'))"

# Contamination + QC
conda install -c bioconda fcs-gx blobtoolkit kraken2 busco compleasm
```

## Related Skills

- **phylogenetics** - Tree building (modern-tree-inference, bayesian-inference, species-trees, divergence-dating) is required input for nearly all comparative-genomics workflows
- **alignment** - MSA (multiple-alignment, alignment-trimming) precedes most analyses; PRANK / MACSE / PREQUAL / HmmCleaner
- **population-genetics** - association-testing, selection-statistics, linkage-disequilibrium, population-structure
- **causal-genomics** - Heritability partitioning with conservation annotations; trait-MR with phylogenetic context
- **genome-annotation** - Annotation-transfer, prokaryotic-annotation, functional-annotation; comparative-annotation-projection complements
- **genome-assembly** - assembly-qc (BUSCO / Compleasm) gates inclusion
- **variant-calling** - SV-calling complements WGA-based synteny / SyRI
- **metagenomics** - GTDB-Tk for MAG taxonomy; AMR + mobile-element annotation for HGT context
- **clinical-databases** - HLA / KIR analyses link to pangenome-analysis PGR-TK
