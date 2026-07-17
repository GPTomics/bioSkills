---
name: bio-chipseq-super-enhancers
description: Identifies super-enhancers from H3K27ac, MED1, or BRD4 ChIP-seq using ROSE, ROSE2, LILY, HOMER -style super, and ENCODE dELS cross-referencing. Handles peak stitching parameters, ranking choices, hockey-stick inflection, marker choice (H3K27ac vs MED1/BRD4), and cross-condition comparison with spike-in normalization. Constructs core regulatory circuitry (Saint-Andre 2016) from SE-encoded TFs. Use when identifying cell-identity / cancer-associated regulatory domains, comparing super-enhancers between conditions, identifying master transcription factor networks, or predicting BET-inhibitor responsiveness.
tool_type: mixed
primary_tool: ROSE
---

## Version Compatibility

Reference examples tested with: ROSE (stjude/ROSE, 2018+), ROSE2 (linlabbcm/rose2, 2021+), LILY (BoevaLab/LILY, 2020+), HOMER 4.11+, samtools 1.19+, bedtools 2.31+, GenomicRanges 1.54+.

The original Young-lab ROSE is Python 2; ROSE2 (linlabbcm/rose2) and the stjude/ROSE fork are the Python-3 implementations with the same algorithm. For hg38 data use stjude/ROSE (`python ROSE_main.py`, whose genomeDict includes HG38); rose2's released genomeDict covers only HG18/HG19/MM8/MM9/MM10/RN4/RN6, so `rose2 -g HG38` fails. LILY (Boeva 2017) is a refactored implementation with input-control background subtraction for low-quality H3K27ac data.

# Super-Enhancer Calling

**"Identify super-enhancers driving cell identity / cancer biology"** -> Stitch nearby active enhancer peaks (H3K27ac, MED1, or BRD4) within a stitching window, exclude proximal-promoter signal, rank by total signal, find the hockey-stick inflection point where signal sharply increases, and classify all stitched regions above the inflection as super-enhancers.

- CLI (stjude/ROSE, Python 3): `python ROSE_main.py -g HG38 -i peaks.gff -r h3k27ac.bam -c input.bam -s 12500 -t 2500 -o rose_out/`
- CLI (HOMER): `findPeaks tag_dir/ -style super -i input_tag_dir/`
- CLI (LILY): variant with input-control background subtraction
- R (custom hockey-stick): rank enhancers by signal, find tangent-line inflection

The SE concept (Whyte 2013) is a thresholding heuristic on a continuous signal distribution (Pott & Lieb 2015 *Nat Genet*), not a categorical biological category. Genetic dissection of super-enhancers (Hay 2016; Moorthy 2017) shows constituent elements contribute unequally and many are individually dispensable/redundant; the "SE" label is a useful operational definition for BET-inhibitor responsiveness and cell-identity gene regulation, not an absolute biological property.

## Marker Choice: H3K27ac vs MED1 vs BRD4

| Marker | Captures | When to prefer |
|--------|----------|----------------|
| **H3K27ac** | Active regulatory elements broadly | Most widely available; standard for SE definition since Whyte 2013 |
| **MED1** | Mediator complex accumulation (the defining biology) | Direct readout of SE; less common antibody; lower signal-to-noise |
| **BRD4** | BET cofactor accumulation | Most predictive of BET-inhibitor responsiveness; clinical relevance |
| **H3K27ac + MED1 intersection** | High-confidence SE | Gold standard if both available |
| **dELS from ENCODE cCREs** | Cell-type-agnostic distal enhancer registry | Cross-reference; not SE-specific by itself |

**Operational rule:** H3K27ac for discovery; MED1 or BRD4 ChIP for functional / therapeutic claims. SE called on H3K27ac alone may not respond to BET inhibitors; SE called on BRD4 will.

## Algorithmic Taxonomy

| Tool | Method | Strength | Fails when |
|------|--------|----------|------------|
| **ROSE** (Whyte 2013) | Stitch within 12.5 kb, exclude ±2.5 kb of TSS, rank by signal, hockey-stick inflection | Original; widely cited; canonical reference | Original Young-lab code is Python 2; run the stjude/ROSE Py3 fork instead |
| **ROSE2** (linlabbcm/rose2) | Same algorithm, Python 3 port | Maintained; pip-installable `rose2` console command | Released genomeDict has no HG38 (HG18/HG19/MM8/MM9/MM10/RN4/RN6 only) -> use stjude/ROSE for hg38 |
| **LILY** (Boeva 2017) | ROSE-like with input-control background subtraction | Works on lower-quality H3K27ac data; subtracts input | Adds complexity; less validated; specific to neuroblastoma/glioma in original paper |
| **HOMER `-style super`** | Native ROSE-like in HOMER framework; stitching without TSS exclusion | Integrated with HOMER workflow | Different stitching defaults; not directly comparable to ROSE counts |
| **Custom hockey-stick (R)** | Generic rank-by-signal + tangent inflection | Flexible; works on any signal definition | Reinvents algorithm; verify against ROSE on known dataset |

**Most papers use ROSE/ROSE2 with default stitching (12.5 kb) and TSS exclusion (2.5 kb).** This is the de facto standard for cross-paper comparison. HOMER's `-style super` produces different counts and is not directly comparable.

## Decision Tree: SE Calling Workflow

| Scenario | Recommended pipeline |
|----------|----------------------|
| Standard SE discovery, H3K27ac available | ROSE2 with default `-s 12500 -t 2500`; input control for subtraction |
| Predict BET-inhibitor response | BRD4 ChIP -> ROSE2 (or H3K27ac SE intersected with BRD4 peaks) |
| Compare SE between conditions (drug treatment) | ROSE2 per condition + spike-in normalization (HDACi/BETi/EZH2i need ChIP-Rx) |
| Build core regulatory circuitry | ROSE2 + Saint-Andre 2016 algorithm: identify TFs encoded by SE that bind own SE + cross-bind other SE-encoded TFs |
| Low-quality H3K27ac (low FRiP) | LILY with input subtraction |
| Compare with ENCODE dELS atlas | ROSE2 + intersect with ENCODE cCRE dELS BED |
| Differential SE between conditions | ROSE2 per condition + signal-quantitative differential (DiffBind on SE regions) |

## ROSE / ROSE2 Workflow

**Goal:** Identify super-enhancers by stitching nearby active enhancer peaks within a stitching distance and ranking by total signal.

**Approach:** Convert peaks to GFF, exclude promoter-proximal peaks via `-t` (TSS exclusion window), stitch enhancers within `-s` (default 12.5 kb), rank by total H3K27ac (or MED1/BRD4) signal, find the hockey-stick inflection point, classify regions above as super-enhancers.

```bash
# Install ROSE2 (Python 3 port; unmaintained ROSE Py2 not recommended)
git clone https://github.com/linlabbcm/rose2.git
pip install ./rose2

# Convert peaks BED to GFF (ROSE requires GFF input)
awk 'BEGIN{OFS="\t"} {print $1,"peaks","enhancer",$2,$3,".",$6,".","ID="NR}' \
    peaks.narrowPeak > peaks.gff

# Filter promoter peaks before SE calling (within 2.5 kb of TSS)
# ROSE handles this via -t flag; preferable to pre-filter for clarity
bedtools intersect -a peaks.narrowPeak -b promoters_2kb.bed -v > enhancer_peaks.bed

# Run stjude/ROSE with input control (Python-3 fork; genomeDict includes HG38)
python ROSE_main.py -g HG38 -i peaks.gff \
    -r h3k27ac.bam -c input.bam \
    -o rose_output/ \
    -s 12500 \
    -t 2500
```

ROSE outputs:
- `*_AllEnhancers.table.txt` — all stitched enhancer regions ranked by signal
- `*_SuperEnhancers.table.txt` — SE only (above hockey-stick inflection)
- `*_Enhancers_withSuper.bed` — BED with SE / TE classification
- `*_Plot_points.png` — hockey-stick plot

## Cross-Condition SE Comparison

This is the analysis most often done wrong. SE calling thresholds depend on absolute signal, so any global shift (HDACi, BETi, EZH2i) confounds direct SE-count comparison.

**Wrong approach:** Call ROSE2 on condition A and condition B separately, intersect SE BEDs, report "gained/lost SE."

**Right approach:**
1. Spike-in normalize signal between conditions (see chip-seq/spike-in-normalization)
2. Build a union SE set from both conditions
3. Quantify signal at union SE regions per condition (DiffBind on the union)
4. Apply differential testing with appropriate normalization (background-bin TMM or spike-in)

```r
library(DiffBind)
# Union of SE BED files from condition A and B
union_se <- rtracklayer::import('union_SE.bed')
# Run DiffBind quantification on this region set with spike-in normalization
```

For BET-inhibitor experiments: the biology IS that all SE decrease globally; spike-in is mandatory.

## Core Regulatory Circuitry (Saint-André 2016)

The CRC algorithm identifies master TF networks from SE annotations:

1. List all TFs encoded by SE-associated genes
2. For each such TF, check if its motif appears in its own SE (auto-regulation)
3. Build a graph where TFs encoded by SE-A bind to motifs in SE-B
4. Identify highly-interconnected sub-networks (CRC)

```bash
# Install CRC pipeline (console command is `crc`; -g is the genome BUILD, not a GTF)
pip install git+https://github.com/linlabcode/CRC.git

# Requires: SE enhancer table, subpeak BED, chromosome-FASTA dir
crc -e SE_table.txt -g HG38 -s subpeaks.bed -c chroms/ -o crc_out/ -n SAMPLE
```

CRC outputs the connected components of the regulatory network. Master TFs typically appear in the largest component with high out-degree.

## ENCODE dELS Cross-Reference

ENCODE distal Enhancer-Like Signatures (dELS) are the cell-type-agnostic regulatory atlas (see chip-seq/peak-annotation). Cross-referencing SE against dELS:

- Validates SE constituents are at canonical regulatory elements
- Identifies SE constituents NOT in the dELS registry (potentially cell-type-specific)
- Provides chromatin-state context (DNase + H3K27ac signatures)

```bash
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed
awk -F'\t' '$NF == "dELS"' GRCh38-cCREs.bed > dels.bed

# Fraction of SE constituents overlapping dELS
bedtools intersect -a SuperEnhancers.bed -b dels.bed -u | wc -l
bedtools intersect -a SuperEnhancers.bed -b dels.bed -wa -wb > se_with_dels.tsv
```

## Per-Tool Failure Modes

### ROSE -- Python 2 dependency

**Trigger:** Running the original Young-lab `ROSE_main.py` (Python 2 code) on a modern system.

**Mechanism:** The original ROSE is Python 2 code; `print` statements without parens, `dict.iteritems()`, etc.

**Symptom:** SyntaxError on first import.

**Fix:** Use the stjude/ROSE fork (`python ROSE_main.py`, Python 3, genomeDict includes HG38) or ROSE2 (`rose2`, Python 3, but its released genomeDict has no HG38 -- HG18/HG19/MM8/MM9/MM10/RN4/RN6 only); identical algorithm and output format.

### ROSE / ROSE2 -- Stitching distance default not appropriate for all biology

**Trigger:** Using default `-s 12500` (12.5 kb) on small genomes or compact gene structures.

**Mechanism:** Default was set on human/mouse vertebrate genomes; Drosophila / yeast / plants have different regulatory architecture.

**Fix:** For non-vertebrate genomes, reduce stitching distance proportionally (e.g., -s 2500 for Drosophila, -s 500 for yeast).

### ROSE / ROSE2 -- TSS exclusion can remove promoter-associated enhancers

**Trigger:** Default `-t 2500` (exclude peaks within 2.5 kb of TSS) on promoter-proximal enhancers (e.g., pELS class).

**Mechanism:** TSS-proximal enhancers are filtered out; SE definition becomes distal-only.

**Symptom:** Lower SE counts than expected for cell types with promoter-enhancer architecture (e.g., human ES cells).

**Fix:** Reduce TSS exclusion to `-t 500` or `-t 0` if including promoter-proximal regulatory regions; document the decision.

### H3K27ac SE vs BRD4 SE -- BET-inhibitor mismatch

**Trigger:** Calling SE on H3K27ac and claiming BET-inhibitor responsiveness.

**Mechanism:** H3K27ac marks active enhancers broadly; not all H3K27ac-positive SE have BRD4 accumulation.

**Symptom:** Predicted BET-sensitive genes don't respond to BET inhibitors in cell-based assays.

**Fix:** For BET-inhibitor claims, use BRD4 ChIP for SE calling, or intersect H3K27ac SE with BRD4 peaks.

### Cross-condition SE counting -- Wrong normalization

**Trigger:** Comparing SE counts in HDACi-treated vs DMSO without spike-in normalization.

**Mechanism:** SE calling thresholds depend on absolute signal; HDACi globally increases H3K27ac, raising every region's signal and shifting the hockey-stick inflection.

**Symptom:** Reports "1000 SE in HDACi vs 500 in DMSO" when biology is just global H3K27ac increase.

**Fix:** Spike-in normalize BAMs (ChIP-Rx with Drosophila chromatin), call SE on scaled signal; OR quantify signal at a union peak set rather than calling SE per condition.

### LILY -- Input subtraction artifacts

**Trigger:** Running LILY without high-quality matched input control.

**Mechanism:** LILY subtracts background based on input signal; mismatched input introduces artifactual negative signal.

**Fix:** Use LILY only when input quality is good (same library prep, same depth, same fragmentation); otherwise use ROSE2 with standard input handling.

### Hockey-stick inflection -- Sensitive to peak count

**Trigger:** Calling SE on a small peak set (< 5000 enhancers).

**Mechanism:** Hockey-stick inflection depends on having a long "tail" of typical enhancers; few peaks distort the inflection.

**Symptom:** SE count is unreasonably high (50%+ of all peaks called SE) or unreasonably low (< 50 SE).

**Fix:** Require ≥ 5000 enhancer peaks input to ROSE2; if fewer, use absolute signal cutoff (e.g., top 5% by signal density) rather than hockey-stick.

## Reconciliation: When SE Calls Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ROSE2 vs HOMER -style super differ | Different stitching distance / TSS handling | Use ROSE2 standard for cross-paper comparison; HOMER for HOMER-integrated workflows |
| H3K27ac SE ≠ MED1 SE at same locus | H3K27ac is broad; MED1 marks subset of active SE | MED1 SE is the more functional definition; H3K27ac includes inactive-but-acetylated regions |
| SE called in DMSO but not in BETi (or vice versa) | Global signal shift confounds threshold | Spike-in normalize; compare quantitatively at union SE set |
| SE shifts location between replicates | Marginal calls below inflection; hockey-stick inflection noisy | Use top N SE by rank for robustness; or require SE in ≥ 2/3 replicates |
| LILY and ROSE2 disagree on SE count | LILY's input subtraction differs | Trust ROSE2 unless input quality is poor (low FRiP) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `SyntaxError: invalid syntax` in ROSE | Python 2 codebase | Use ROSE2 (linlabbcm) Python 3 port |
| GFF format error | Wrong column ordering | Use awk template: `chr<TAB>peaks<TAB>enhancer<TAB>start<TAB>end<TAB>.<TAB>strand<TAB>.<TAB>ID=N` |
| ROSE2 reports 0 SE | Hockey-stick inflection failed; too few enhancers | Inspect `_Plot_points.png`; ≥ 5000 peaks input recommended |
| Genome flag error in ROSE2 | Genome not pre-configured | Genome flag must be one of HG18, HG19, HG38, MM8, MM9, MM10 |
| All SE at promoters | TSS exclusion too narrow OR data dominated by promoter signal | Verify `-t 2500`; check input is H3K27ac at enhancers not full chromatin |
| Cross-condition SE gain/loss not reproducible | No spike-in normalization | Spike-in (ChIP-Rx) or quantitative differential on union SE |

## References

- Whyte WA et al 2013 Cell 153:307 (super-enhancers, ROSE)
- Lovén J et al 2013 Cell 153:320 (SE characterization, BET sensitivity)
- Hnisz D et al 2013 Cell 155:934 (SE in cell identity)
- Pott S & Lieb JD 2015 Nat Genet 47:8 (SE as continuum critique)
- Lin CY et al 2016 Nature 530:57-62 (medulloblastoma super-enhancers)
- Saint-André V et al 2016 Genome Res 26:385 (core regulatory circuitry)
- Boeva V et al 2017 Nat Genet 49:1408 (LILY; neuroblastoma SE)
- Hnisz D et al 2017 Cell 169:13 (phase-separation model of transcriptional control)
- Hay D et al 2016 Nat Genet 48:895 (genetic dissection of the alpha-globin super-enhancer)
- Moorthy S et al 2017 Genome Res 27:246 (SE constituents have equivalent/redundant regulatory roles in ESCs)
- Sengupta S & George RE 2017 Trends Cancer 3:269 (SE function review)

## Related Skills

- chip-seq/peak-calling - Generate H3K27ac / MED1 / BRD4 peaks for SE input
- chip-seq/chipseq-qc - Filter hyper-ChIPable peaks before SE calling
- chip-seq/spike-in-normalization - Mandatory for cross-condition SE comparison
- chip-seq/differential-binding - Quantitative differential testing on union SE set
- chip-seq/peak-annotation - Annotate SE-associated genes; cross-reference dELS
- chip-seq/cut-and-run-tag - SE calling on CUT&RUN/CUT&Tag H3K27ac (different spike-in)
- atac-seq/enhancer-gene-linking - ENCODE-rE2G for SE-target gene assignment
- data-visualization/genome-tracks - SE region visualization
