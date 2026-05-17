---
name: bio-chipseq-chromatin-state-segmentation
description: Segments the genome into chromatin states from combinatorial histone modification and chromatin factor ChIP-seq data. Uses ChromHMM (multivariate HMM on binarized signal, v1.27), Segway (Dynamic Bayesian Network on continuous signal), EpiSegMix (flexible-distribution HMM with duration modeling, 2024), EpiLogos (multi-biosample visualization), IDEAS (cell-type-aware joint), and full-stack ChromHMM (Vu Ernst 2022) for cross-cell-type segmentations. Handles state-count selection (15 vs 18 vs 25 states), binarization choice, OverlapEnrichment / NeighborhoodEnrichment downstream analysis, and cross-biosample integration. Use when learning chromatin states from a histone mark panel, characterizing learned states by genomic feature enrichment, or comparing chromatin landscapes across cell types.
tool_type: cli
primary_tool: ChromHMM
---

## Version Compatibility

Reference examples tested with: ChromHMM 1.27+, Segway 3.0+, EpiSegMix 1.0+, EpiLogos (Meuleman lab), IDEAS 1.20+, samtools 1.19+, bedtools 2.31+. ChromHMM requires Java 8+; runs as `java -mx<MEMORY> -jar ChromHMM.jar <command>`.

# Chromatin State Segmentation

**"Integrate multiple histone modification ChIP-seq tracks into chromatin states"** -> Learn a small set of recurring combinatorial patterns of histone marks (active promoter, active enhancer, poised enhancer, polycomb-repressed, heterochromatic, transcribed, etc.) and segment the genome by which state each region belongs to. Output: per-state genomic intervals, state-by-mark emission matrix, and state-state transition matrix.

- CLI (canonical): ChromHMM `BinarizeBam` -> `LearnModel` -> `OverlapEnrichment` / `NeighborhoodEnrichment`
- CLI (continuous signal): Segway `train` -> `posterior` -> `identify`
- CLI (flexible distributions): EpiSegMix (2024)
- Visualization across biosamples: EpiLogos (Meuleman lab)
- Cell-type-aware joint: IDEAS

Chromatin state segmentation requires a panel of histone marks; minimum 4-5 marks (e.g., H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3) for meaningful states. With fewer marks, simpler peak-based annotation (chipseq/peak-annotation) is more appropriate.

## Tool Taxonomy

| Tool | Method | Strength | Fails when |
|------|--------|----------|------------|
| **ChromHMM** (Ernst & Kellis 2012; v1.27 current) | Multivariate HMM on binarized 200 bp bins | Canonical; widely used; integrated with Roadmap Epigenomics 15-state model; mature toolchain | Binarization throws away signal quantitation; default 200 bp bins may be too coarse for sharp boundaries |
| **Segway** (Hoffman 2012) | Dynamic Bayesian Network on continuous signal | Higher resolution; uses signal magnitudes not binarized | More complex setup; slower; less standardized output |
| **EpiSegMix** (Sczyrba 2024) | HMM with flexible read-count distributions + duration modeling | Modern; handles both narrow and broad mark distributions in one model | Newer; smaller user base |
| **EpiLogos** (Meuleman lab) | Multi-biosample visualization tool | Built on top of ChromHMM/Segway segmentations; compare ChromHMM states across 100s of biosamples | Visualization tool, not a segmentation method itself |
| **IDEAS** (Zhang 2016) | Cell-type-aware joint inference | Across-cell-type segmentation respecting cell-type identity | Slower; complex parameter tuning |
| **EpiCSeg** (Mammana 2015) | Negative binomial mixture | Read-count-based; doesn't need binarization | Less standardized output |
| **GenoSTAN** | HMM with various emission distributions | Flexible | Less actively developed |
| **Roadmap 25-state model** (Kundaje 2015) | ChromHMM 25-state precomputed model | Reference for cross-cell-type interpretation | Tied to Roadmap mark panel (5 core marks) |
| **Full-stack ChromHMM** (Vu Ernst 2022) | 100-state segmentation across 1032 datasets / 127 reference epigenomes | Comprehensive cross-tissue annotation | Computationally intensive to retrain |

## ChromHMM Workflow

ChromHMM is the de facto standard. The workflow has 4 stages:

### Step 1: Binarize ChIP-seq signal

```bash
# Build cellMarkFileTable: cell_type<TAB>mark<TAB>file<TAB>(optional control)
cat > cellMarkFileTable.txt << EOF
GM12878	H3K4me3	gm12878_h3k4me3.bam	gm12878_input.bam
GM12878	H3K27me3	gm12878_h3k27me3.bam	gm12878_input.bam
GM12878	H3K27ac	gm12878_h3k27ac.bam	gm12878_input.bam
GM12878	H3K4me1	gm12878_h3k4me1.bam	gm12878_input.bam
GM12878	H3K36me3	gm12878_h3k36me3.bam	gm12878_input.bam
EOF

# Binarize BAMs into 200 bp bins; emission = whether mark exceeds Poisson threshold
java -mx16G -jar ChromHMM.jar BinarizeBam \
    -b 200 \
    chromsizes_hg38.txt \
    bam_dir/ \
    cellMarkFileTable.txt \
    binarized_output/
```

Output: per-chromosome `_binary.txt` files, one row per 200 bp bin, columns = marks, values 0/1.

### Step 2: Learn model

```bash
# Train HMM with N states; common choices: 15, 18, 25
# 15 states: Ernst & Kellis 2011 model; canonical
# 18 states: extends with additional regulatory states
# 25 states: Roadmap Epigenomics extended model
java -mx16G -jar ChromHMM.jar LearnModel \
    -p 8 \
    binarized_output/ \
    model_15state/ \
    15 \
    hg38

# Output: model_15state.txt (emission + transition matrices),
# emissions_15.png (visualization), transitions_15.png,
# per-chromosome _segments.bed (state assignments)
# AND automatically runs OverlapEnrichment + NeighborhoodEnrichment
```

### Step 3: Interpret states from emission matrix

| Roadmap 15-state assignments (canonical) |
|-------------------------------------------|
| 1_TssA — Active TSS (high H3K4me3, H3K27ac) |
| 2_TssAFlnk — Flanking TSS (H3K4me3, H3K27ac, lower) |
| 3_TxFlnk — Transcript flanking |
| 4_Tx — Strong transcription (H3K36me3, H3K79me2 if available) |
| 5_TxWk — Weak transcription |
| 6_EnhG — Enhancer in gene body (H3K4me1, H3K27ac) |
| 7_Enh — Generic enhancer (H3K4me1, H3K27ac) |
| 8_ZNF/Rpts — Zinc-finger / repeats |
| 9_Het — Heterochromatin (H3K9me3) |
| 10_TssBiv — Bivalent TSS (H3K4me3 + H3K27me3) |
| 11_BivFlnk — Bivalent flanking |
| 12_EnhBiv — Bivalent enhancer (H3K4me1 + H3K27me3) |
| 13_ReprPC — Polycomb-repressed (H3K27me3) |
| 14_ReprPCWk — Weak Polycomb |
| 15_Quies — Quiescent (no signal) |

### Step 4: Functional enrichment of states

```bash
# OverlapEnrichment: enrichment of each state for external feature sets
java -mx16G -jar ChromHMM.jar OverlapEnrichment \
    model_15state/GM12878_15_segments.bed \
    /path/to/anchor_files/ \
    enrichment_output/GM12878 \
    -labels

# NeighborhoodEnrichment: enrichment relative to anchor positions (e.g., TSS)
java -mx16G -jar ChromHMM.jar NeighborhoodEnrichment \
    model_15state/GM12878_15_segments.bed \
    /path/to/tss_anchors.txt \
    enrichment_output/GM12878_TSS \
    -labels
```

Anchor files: BED files of features (CGIs, repeats, conserved elements, etc.) for OverlapEnrichment; position files for NeighborhoodEnrichment.

## Choosing State Count

| States | Use case | Mark panel size |
|--------|----------|-----------------|
| 8-10 | Initial exploration; small mark panel (3-4 marks) | 3-5 marks |
| 15 | Roadmap Epigenomics canonical | 5 core (H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3) |
| 18 | Roadmap extended (adds bivalent states, fine enhancer subtypes) | 5-7 marks |
| 25 | Roadmap Epigenomics extended; cross-cell-type compatibility | 6+ marks |
| 50+ | Full-stack model (Vu Ernst 2022) | Many marks across many cell types |

**Practical workflow:** Train at N=15, 18, 25; compare emission matrices; choose the smallest N where biology is interpretable. Higher N risks over-segmentation (state splitting random variation).

## Segway Workflow

```bash
# Train Segway model (10 states by default; supports more)
segway train \
    --num-labels 25 \
    --num-instances 3 \
    --resolution 100 \
    chromsizes_hg38.bed \
    h3k4me3.bw h3k27ac.bw h3k4me1.bw h3k36me3.bw h3k27me3.bw \
    --traindir traindir/

# Posterior inference + segmentation
segway posterior --traindir traindir/ --identifydir identifydir/ \
    chromsizes_hg38.bed \
    h3k4me3.bw h3k27ac.bw h3k4me1.bw h3k36me3.bw h3k27me3.bw

# Output: identifydir/segway.bed (state assignments)
```

Segway uses bigWig (continuous signal) vs ChromHMM binarized binary. Trade-off: more information per region (continuous) but more complex training.

## EpiLogos Visualization

EpiLogos doesn't perform segmentation; it visualizes existing ChromHMM/Segway segmentations across many biosamples (epilogos.org).

```bash
# Use precomputed ChromHMM segmentations across multiple cell types
# Web interface: https://epilogos.altius.org/
# Local: github.com/meuleman/epilogos
```

Useful for: cross-cell-type comparison; identifying tissue-specific regulatory states; cohort-level chromatin landscape summaries.

## Full-Stack ChromHMM (Vu Ernst 2022)

The full-stack model trained on 1032 datasets / 127 reference epigenomes:

```bash
# Use precomputed model from Ernst lab
# github.com/ernstlab/full_stack_ChromHMM_annotations
# Annotate new sample by applying model to binarized data
java -mx16G -jar ChromHMM.jar MakeSegmentation \
    full_stack_model_100states.txt \
    binarized_sample/ \
    full_stack_output/
```

Useful for: applying a comprehensive cross-tissue annotation to a new sample; comparing to canonical Roadmap states.

## Per-Tool Failure Modes

### ChromHMM -- Bin size 200 bp too coarse for sharp boundaries

**Trigger:** Studying TF binding boundaries or sharp enhancer transitions at 200 bp resolution.

**Mechanism:** ChromHMM default 200 bp bins; biology may shift within a bin.

**Fix:** Reduce to `-b 100` or `-b 50` (smaller bin); increases memory and compute time but improves boundary resolution. Re-train model at finer resolution.

### ChromHMM -- Binarization throws away signal quantitation

**Trigger:** Distinguishing low- from high-signal regions of the same state.

**Mechanism:** ChromHMM binarizes each 200 bp bin to 0/1 per mark; state assignment uses combinatorial pattern, not magnitude.

**Fix:** Use Segway (continuous signal) or EpiSegMix (flexible distributions) for magnitude-aware segmentation.

### ChromHMM / Segway -- Wrong state count

**Trigger:** Training with N=50 states on a 4-mark panel; or N=10 on a 7-mark panel.

**Mechanism:** Excess states fragment biology; insufficient states force unrelated regions into the same state.

**Symptom:** Emission matrix shows redundant states (multiple states with same emission profile) at high N; or biologically distinct regions lumped together at low N.

**Fix:** Train at N=15, 18, 25; inspect emission matrix similarity; choose the smallest N where states are interpretable as distinct biology.

### Mark panel mismatch with model

**Trigger:** Applying Roadmap 25-state model to a sample with different mark panel.

**Mechanism:** Model was trained on specific marks; emission probabilities are mark-specific. Applying to different mark panel produces nonsensical state assignments.

**Fix:** Either train a new model on the available mark panel; OR ensure the exact same marks (and ordering) as used in the model.

### IgG control instead of input

**Trigger:** Using IgG controls in `BinarizeBam` for histone marks.

**Mechanism:** ChromHMM's binarization compares mark signal to control; histone mark biology assumes input (sonicated chromatin) as background, not IgG.

**Fix:** Use sonicated input as control for histone mark ChIP. IgG is not appropriate for ChromHMM binarization of histone marks.

### Cross-cell-type state mapping

**Trigger:** Training separate models per cell type and trying to compare state assignments.

**Mechanism:** State 5 in cell type A may not correspond to state 5 in cell type B if trained independently.

**Fix:** Train one model on concatenated data from all cell types (joint segmentation); or apply a single precomputed model (Roadmap 15-state, full-stack) to all samples for consistent state labels.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ChromHMM and Segway segments differ | Different bin sizes / binarization vs continuous | Both can be valid; inspect emission matrices; pick the tool matching the resolution needs |
| State assignment varies wildly between replicates | Insufficient marks; over-binned | Increase mark panel; reduce state count |
| Active TSS state overlaps polycomb state at promoters | Bivalent biology (Bernstein 2006) | Expected for ESC-like cells; not an error; consider bivalent-specific state in N=18 model |
| Roadmap 25-state model annotates unknown cell type | Cross-cell-type generalization | Use cautiously; verify against tissue-specific tracks |
| Heterochromatin (H3K9me3) state has too many bins | H3K9me3 covers large fraction of genome | Expected; heterochromatin is genome-wide |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `java.lang.OutOfMemoryError` | Insufficient JVM heap | `java -mx32G -jar ChromHMM.jar ...` |
| `BinarizeBam` very slow | Large BAMs without index | `samtools index` all BAMs first |
| All states have similar emissions | Mark panel too small | Need at least 5 marks for canonical 15-state model |
| Segments file empty for some chromosomes | Chromosome not in chromsizes file | Add or use `-chrom` flag to restrict |
| State labels don't match Roadmap | Trained model independently | Use Roadmap precomputed model OR map states by emission similarity |
| ChromHMM "no signal in marks" | All bins binarized to 0 | Check signal quality; verify control normalization |

## References

- Ernst J & Kellis M 2012 Nat Methods 9:215 (ChromHMM v1)
- Ernst J & Kellis M 2017 Nat Protoc 12:2478 (ChromHMM protocol)
- Hoffman MM et al 2012 Nat Methods 9:473 (Segway)
- Sczyrba A et al 2024 Bioinformatics (EpiSegMix)
- Meuleman W et al 2020 Nature 583:744 (EpiLogos / DHS index)
- Zhang Y & Hardison 2016 Nucleic Acids Res 44:6721 (IDEAS)
- Mammana A & Chung HR 2015 Genome Biol 16:151 (EpiCSeg)
- Roadmap Epigenomics Consortium 2015 Nature 518:317 (Roadmap 25-state model)
- Vu H & Ernst J 2022 Nat Commun 13:2783 (full-stack ChromHMM)
- Kundaje A et al 2015 Nature 518:317 (Roadmap integrative analysis)

## Related Skills

- chip-seq/peak-calling - Peak calling per mark before segmentation
- chip-seq/chipseq-qc - Replicate concordance per mark; QC before integration
- chip-seq/peak-annotation - cCRE classification (PLS/pELS/dELS) complementary to chromatin states
- chip-seq/spike-in-normalization - Spike-in normalize per-mark BAMs before binarization for cross-condition state comparison
- atac-seq/single-cell-atac - scATAC + ChIP integration via multimodal methods
- machine-learning/model-validation - Model selection (state count); cross-validation
- data-visualization/genome-tracks - Visualize state segmentations
- gene-regulatory-networks/coexpression-networks - Cross-reference chromatin states with co-expression modules
