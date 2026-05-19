---
name: bio-clip-seq-stamp-antibody-free
description: Profile RNA-binding protein targets without antibody or UV crosslinking using STAMP (APOBEC1-RBP fusion, C-to-U editing), scSTAMP (single-cell), TRIBE/HyperTRIBE (ADAR-RBP, A-to-I editing), DART-seq (APOBEC1-YTH for m6A), or Bullseye/SAILOR edit-site detection pipelines. Use when antibody is unavailable or specificity is doubtful, when single-cell RBP profiling is needed (scSTAMP), or when in vivo RBP profiling without UV is preferred.
tool_type: mixed
primary_tool: STAMP
---

## Version Compatibility

Reference examples tested with: STAMP / scSTAMP (Brannan 2021/2024 Yeo lab github), Bullseye 1.0+, SAILOR 1.1+, samtools 1.19+, REDItools2 1.3+, JACUSA2 2.0+, scanpy 1.10+, anndata 0.10+, pysam 0.22+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt the example to match the actual CLI rather than retrying.

# STAMP / Antibody-Free RBP Profiling

**"Profile RBP-RNA targets without UV crosslinking or immunoprecipitation"** -> Express a fusion of the RBP-of-interest with a deaminase (APOBEC1 for STAMP, ADAR for TRIBE) in cells; the deaminase edits RNA nucleotides adjacent to where the RBP binds, producing a C-to-U (STAMP) or A-to-I (TRIBE, read as A-to-G) editing signature in standard RNA-seq. The targets are recovered computationally from the editing pattern. Three properties make this approach valuable: (a) no UV crosslinking required (works in tissue/in vivo); (b) no IP step (no antibody needed - the RBP itself targets the deaminase); (c) compatible with single-cell readout because the editing signal exists in standard scRNA-seq (scSTAMP, scTRIBE). Trade-off: editing is offset from the binding site (typically 0-50 nt away); resolution is approximate; off-target editing from deaminase alone must be subtracted.

- CLI (STAMP, bulk): standard RNA-seq pipeline + Bullseye or SAILOR for C-to-U edit detection vs APOBEC1-only control
- CLI (TRIBE, bulk): standard RNA-seq + REDItools2 or JACUSA2 for A-to-I edit detection vs ADAR-only control
- CLI (DART-seq for m6A): same as STAMP, with APOBEC1-YTH fusion (YTH is the m6A reader)
- Python (scSTAMP single-cell): 10x Genomics or Smart-seq2 pipeline + custom editing-rate quantification per cell + per-cell binding-target inference
- CLI (general edit-site detection): `JACUSA2 call-2 -r ref.fa -p 8 -F 1024 -A,B treated.bam,control.bam -t pileup.tsv` then filter for C-to-U or A-to-I

STAMP (Brannan 2021) is the canonical antibody-free RBP profiling method. TRIBE (McMahon 2016) and HyperTRIBE (Xu 2018) are earlier ADAR-based variants. DART-seq (Meyer 2019) is the m6A-specific application using YTH-fused APOBEC1. scSTAMP (Brannan 2024) extends STAMP to single-cell readout.

## Methods Taxonomy

| Method | Deaminase | Edit signature | Cells supported | Single-cell | Strength | Fails when |
|--------|-----------|----------------|-----------------|-------------|----------|------------|
| STAMP (Brannan 2021) | APOBEC1 | C->U in mRNA (reads as C->T) | Any | Yes (scSTAMP) | Antibody-free; no UV; in vivo | APOBEC1 also edits ssDNA off-target; saturated edits at high APOBEC1 expression |
| scSTAMP (Brannan 2024) | APOBEC1 | C->U per cell | Single cell (10x or Smart-seq2) | Yes (native) | Per-cell RBP profiling | Coverage per cell limits sensitivity; ~25% of cytosines accessible per transcript |
| TRIBE (McMahon 2016) | ADAR catalytic domain | A->I (reads as A->G in cDNA) | Drosophila standard; mammalian works | Yes (scTRIBE) | First antibody-free | Edits restricted to certain ADAR consensus; lower edit rate than HyperTRIBE |
| HyperTRIBE (Xu 2018) | ADAR E488Q hyperactive mutant | A->I in much wider context | Drosophila / mammalian | Yes | Higher edit rate than original TRIBE | Hyperactive may edit off-target; needs ADAR-only control |
| TRIBE-DiCo (Erickson 2024) | ADAR catalytic + dimerization | A->I | Mammalian | Yes | Improved specificity | Newer; less validation |
| DART-seq (Meyer 2019) | APOBEC1-YTH | C->U near m6A | Any | Yes (scDART) | m6A reader profiling | Indirect (edits near m6A, not at RBP binding sites) |
| Bullseye | NA (analysis tool) | NA | Any | Yes | STAMP / DART analysis | Just an analysis pipeline |
| SAILOR | NA (analysis tool) | NA | Any | Yes | RNA editing analysis | Just an analysis pipeline |
| REDItools2 | NA (analysis tool) | NA | Any | NA | Generic RNA editing | Generic; not RBP-specific |
| JACUSA2 (Piechotta 2017) | NA (analysis tool) | NA | Any | NA | Multi-sample edit-site detection | Generic; not RBP-specific |
| ADAR-CLIP | NA - this is regular ADAR CLIP | NA | NA | NA | CLIP for ADAR | Not an antibody-free method; just a different CLIP target |

Methodology evolves; the Brannan lab and Yeo lab papers (2021, 2024) are canonical. Verify deaminase fusion expression level (low expression for specificity; saturation degrades specificity).

## STAMP vs m6A-Specific Methods

For m6A profiling specifically, antibody-free choices include:
- **DART-seq (APOBEC1-YTH fusion):** This skill covers the methodology, but for m6A detection see clip-seq/m6a-clip. Only 44% of DART edits fall within DRACH motifs (Liu 2023); strong off-target component.
- **GLORI (Liu 2023):** Antibody-free, chemical, stoichiometric single-base m6A; this is the new (2023) gold standard for m6A. See clip-seq/m6a-clip.
- **m6Anet (Hendra 2022):** Nanopore direct RNA m6A; AUC 0.83 on HEK293T.

If the use case is m6A profiling, the m6a-clip skill is the canonical reference; this skill (stamp-antibody-free) focuses on the broader RBP-editing-fusion paradigm where the target is not m6A but the RBP's RNA targets.

## Critical Choice: STAMP (APOBEC1) vs TRIBE (ADAR)

| Property | STAMP | TRIBE |
|----------|-------|-------|
| Deaminase | APOBEC1 (cytidine -> uridine) | ADAR (adenosine -> inosine) |
| Edit signature | C->U (reads as C->T) | A->I (reads as A->G) |
| ssRNA preference | Yes (APOBEC1 acts on ssRNA + ssDNA) | No (ADAR acts on dsRNA stems by default; ADAR2cd in TRIBE relaxes this) |
| Edit clusters per target | 10-1000 | ~5-50 (lower; HyperTRIBE higher) |
| Off-target | APOBEC1 alone has detectable C->U on ssDNA + RNA | ADAR has weak intrinsic A->I |
| Spatial offset from RBP binding | 0-50 nt | 0-30 nt |
| Cell line tested | HEK293, K562, mouse tissue | Drosophila (original), mouse, human |
| Single-cell | scSTAMP 2024 | scTRIBE 2020 |
| Compatible methods | C->U is rare in mRNA; signal is clean | A->I is common at ALU repeats; baseline ADAR editing competes |
| Cytosine accessibility | ~25-35% of mRNA bases are C; APOBEC1 needs ssRNA | All A residues are potential ADAR targets |

Both work; STAMP has more clusters per target (advantage for low-coverage scenarios) and cleaner background (C->U is rare in mRNA). TRIBE has more flexibility (ADAR variants tunable) and lower off-target. Practical choice often comes down to lab familiarity.

## scSTAMP / scTRIBE Single-Cell Workflow

The defining advantage of antibody-free RBP profiling is compatibility with single-cell readout. scSTAMP processes 10x Genomics or Smart-seq2 libraries.

```bash
# Standard 10x cellranger pipeline produces BAM with per-cell barcodes
cellranger count \
    --id=scstamp_sample \
    --transcriptome=refdata-gex-GRCh38 \
    --fastqs=fastq_dir \
    --localcores=16 --localmem=64

# scSTAMP analysis (Yeo lab github)
# Quantify per-cell C->U editing
python scstamp_analysis.py \
    --bam scstamp_sample/outs/possorted_genome_bam.bam \
    --barcodes scstamp_sample/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
    --control apobec1_only_sample/outs/possorted_genome_bam.bam \
    --output per_cell_edits.h5
```

Per-cell edit-rate matrix can be integrated with standard scRNA-seq clustering. The per-cell binding profile is reconstructed from cells with sufficient coverage (>= 10000 unique reads typically).

## Editing-Site Detection Pipelines

**Goal:** Recover specific (not background) RBP-fusion-induced editing sites from RNA-seq libraries by subtracting the deaminase-only control.

**Approach:** Process fusion-sample BAM and deaminase-only-control BAM in parallel; use Bullseye (STAMP/DART), SAILOR (Yeo), or JACUSA2 (general) to call C-to-U (STAMP/DART) or A-to-I (TRIBE/HyperTRIBE) edit sites at edit rate >= 0.1 and coverage >= 10, requiring fusion-vs-control edit ratio > 3 as the specificity threshold.

**Bullseye** (Meyer 2019) is the canonical DART-seq/STAMP analysis pipeline. It subtracts deaminase-only control and applies edit-rate thresholds.

```bash
# Bullseye for STAMP
bullseye \
    --ip stamp_sample.bam \
    --control apobec1_only.bam \
    --reference genome.fa \
    --edit_type c2t \
    --threshold 0.1 \
    --min_coverage 10 \
    --output stamp_edits.bed
```

**SAILOR** (Yeo lab 2017) is the Yeo lab pipeline for in vivo editing.

```bash
# SAILOR with paired-end RNA-seq
sailor \
    --bam stamp_sample.bam \
    --bg apobec1_only.bam \
    --ref genome.fa \
    --output sailor_out/
```

**JACUSA2** is a general-purpose RNA editing pipeline.

```bash
# JACUSA2 multi-sample edit-site detection
jacusa2 call-2 \
    -r genome.fa \
    -p 8 \
    -F 1024 \
    -A stamp1.bam,stamp2.bam \
    -B control1.bam,control2.bam \
    -o jacusa_edits.tsv

# Post-filter for C->U (STAMP) at edit rate >= 0.1
awk '$5 == "C" && $9 ~ /U/ && $11 >= 0.1' jacusa_edits.tsv > stamp_edits_filtered.tsv
```

## Per-Method Failure Modes

### STAMP -- APOBEC1 over-expression saturation

**Trigger:** Strong APOBEC1-RBP expression (>>10x endogenous).

**Mechanism:** At high APOBEC1 expression, the deaminase saturates editing - every accessible C in mRNA is edited, regardless of RBP binding.

**Symptom:** Edit count per gene >> expected; non-specific editing across mRNAs; APOBEC1-only control has nearly as many edits as the fusion.

**Fix:** Titrate fusion expression with inducible promoter; aim for low-to-moderate expression giving clean fusion-specific edits. Yeo lab convention: doxycycline-inducible with mid-range dox dose. Compare edits in fusion vs APOBEC1-only; require fusion edits / APOBEC1-only edits > 3.

### STAMP -- APOBEC1 off-target on ssDNA

**Trigger:** APOBEC1 expressed in DNA-replicating cells.

**Mechanism:** APOBEC1 has intrinsic ssDNA editing activity; some "C->U" calls are actually genomic ssDNA edits read through transcription.

**Symptom:** Edits cluster at replication-fork regions or LINE-1 elements; non-specific genome-wide.

**Fix:** Bullseye filters genomic SNVs vs RNA edits via strand information; verify mismatch is C->U on the transcribed strand, not the genomic C->T.

### TRIBE -- Editing at ALU repeats

**Trigger:** TRIBE in mammalian cells; many edits at ALU sequences.

**Mechanism:** ADAR has baseline activity at ALU dsRNA structures; this is NOT TRIBE-specific signal. ALU edits dominate the apparent target list.

**Symptom:** Top edited regions are all ALU repeats; target list looks generic.

**Fix:** Subtract ADAR-only control or wild-type ADAR baseline; filter out ALU-overlapping edits unless RBP is known to bind repeats.

### DART-seq -- Spatial offset from m6A

**Trigger:** DART-seq applied with expectation of single-base m6A resolution.

**Mechanism:** APOBEC1-YTH edits Cs 0-50 nt from the YTH-bound m6A site. The exact m6A position is not the edit position.

**Symptom:** DART edits scattered around DRACH motifs; only 44% of edits within DRACH (Liu 2023).

**Fix:** Treat DART edits as "near m6A"; cross-reference with single-base m6A methods (GLORI, miCLIP2). DART is hypothesis-generating, not precise localization.

### scSTAMP -- Coverage limitation per cell

**Trigger:** scSTAMP on 10x library; per-cell coverage limits target detection.

**Mechanism:** Single cells have ~5000-50000 mRNA molecules; editing-rate quantification at any single position needs >= 10 reads. Most positions have 0-3 reads per cell.

**Symptom:** Per-cell binding-target list is sparse; many cells have 0 detected targets.

**Fix:** Aggregate cells into pseudo-bulk by cluster/cell-type; quantify editing at pseudobulk level; or use ultra-deep Smart-seq2 (~1M reads/cell) instead of 10x for higher per-cell coverage.

### No control subtraction

**Trigger:** STAMP/DART/TRIBE run without deaminase-only control.

**Mechanism:** Deaminases have intrinsic baseline editing (APOBEC1 ~3-5% C->U; ADAR ~5-10% A->I at ALUs). Without control, all edits look like signal.

**Symptom:** Edit count enormous; target list non-specific.

**Fix:** Always run deaminase-only (APOBEC1 or ADAR catalytic domain) control in parallel. Bullseye / SAILOR / JACUSA all support control subtraction.

### Strand-specific edit interpretation

**Trigger:** Generic variant caller used instead of edit-aware tool.

**Mechanism:** Variant callers report any C->T mismatch; without strand information it is not possible to distinguish C->U (STAMP signal on transcribed strand) from G->A (the reverse complement of C->T on the genome strand from anti-sense reads).

**Symptom:** Edit count inflated 2x; signal not stranded.

**Fix:** Use editing-specific tool (Bullseye, SAILOR, JACUSA2) that respects strand. Or filter for proper strand: C->U on +sense and G->A on -sense.

## Decision Tree by Use Case

| Scenario | Method | Why |
|----------|--------|-----|
| Antibody for RBP doesn't exist | STAMP or TRIBE | The original use case |
| RBP in tissue / in vivo (no UV possible) | STAMP / TRIBE | No CL needed |
| Single-cell RBP profiling | scSTAMP or scTRIBE | The only practical option |
| m6A reader profiling (YTHDF) | DART-seq (APOBEC1-YTH) | Specific reader fusion |
| Drosophila RBP | TRIBE (original development) | Most validated in Drosophila |
| Mammalian RBP | STAMP (more validated in mammalian) | Yeo lab benchmarks |
| Need precise binding site | Use CLIP / eCLIP not STAMP/TRIBE | Editing is offset from binding |
| Repeat-binding RBP | CLIP + CLAM, not TRIBE | ADAR baseline at ALUs swamps TRIBE |
| Comparison across methods | STAMP + classic eCLIP both | Triangulation increases confidence |
| Low input cell numbers (< 50k) | scSTAMP | Compatible with sparse libraries |
| Time-course binding dynamics | STAMP with inducible expression | Live-cell editing accumulates |
| Cross-link sensitive RBP | STAMP / TRIBE (no UV) | Some RBPs degrade with UV |

## Reconciliation: STAMP vs CLIP

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| STAMP finds targets eCLIP missed | Targets that crosslink poorly; or low-abundance | Validate orthogonally; STAMP often more sensitive for low-abundance targets |
| eCLIP finds targets STAMP missed | RBP-RNA contact too far from accessible C; or APOBEC1 saturated | Check fusion expression level; ssRNA accessibility |
| STAMP edits clustered; eCLIP peaks broader | Spatial offset of editing from binding | Both correct; report at appropriate resolution |
| STAMP top targets generic mRNAs | Saturated APOBEC1; or no control subtraction | Titrate fusion expression; verify APOBEC1-only control |
| TRIBE editing dominated by ALUs | ADAR baseline activity | Subtract ADAR-only; filter ALU repeats |
| Discordant target lists across labs for same RBP | Fusion expression varies; control differs | Standardize protocols; cross-validate |
| scSTAMP pseudobulk = bulk STAMP | Cell aggregation correct | Trust both for low-coverage targets |
| scSTAMP per-cell sparse | Coverage limitation | Pseudobulk by cluster; or use Smart-seq2 |

**Operational rule:** STAMP/TRIBE for hypothesis-generation, antibody-free profiling, or single-cell. CLIP/eCLIP for high-resolution validation. Best paper figure: STAMP + eCLIP concordance for top targets.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Edit count enormous | No control subtraction | Add APOBEC1-only / ADAR-only control |
| Edits in DNA / off-target | APOBEC1 ssDNA activity | Filter genomic SNVs; trust strand-specific RNA edits |
| Saturated editing on every gene | High fusion expression | Titrate down; use inducible promoter |
| TRIBE all edits at ALUs | ADAR baseline | Subtract ADAR-only; filter ALU |
| DART edits not at m6A | Spatial offset (0-50 nt) | Expected; cross-reference single-base m6A |
| scSTAMP per-cell sparse | 10x coverage limit | Pseudobulk by cluster; Smart-seq2 alternative |
| Generic variant caller | No strand awareness | Use Bullseye / SAILOR / JACUSA |
| Edits in non-edited strand | Anti-sense transcription | Filter by strand-specific mate |
| Same target list as RNA-seq abundance | Saturated APOBEC1 | Reduce fusion expression |
| Reproducibility issue across labs | Fusion construct differs | Standardize promoter, tag position, deaminase variant |

## References

- Brannan KW et al 2021 Nat Methods 18:507 (STAMP, APOBEC1-RBP fusion)
- Brannan KW et al 2024 Cell Genomics (scSTAMP single-cell extension)
- McMahon AC et al 2016 Cell 165:742 (TRIBE original, Drosophila)
- Xu W et al 2018 Cell 174:1567 (HyperTRIBE)
- Erickson AW et al 2024 (TRIBE-DiCo)
- Meyer KD 2019 Nat Methods 16:1275 (DART-seq, APOBEC1-YTH)
- Rahman R et al 2018 Bioinformatics 34:1457 (SAILOR)
- Piechotta M et al 2017 BMC Bioinformatics 18:7 (JACUSA)
- Picardi E & Pesole G 2013 Bioinformatics 29:1813 (REDItools)
- Tegowski M et al 2022 Mol Cell 82:868 (scDART single-cell DART-seq)
- Hayashi Y et al 2023 (DART-seq benchmarking)

## Related Skills

- clip-seq/m6a-clip - DART-seq is part of the m6A toolkit
- clip-seq/clip-deep-learning - Computational target prediction
- clip-seq/ago-clip-mirna-targets - AGO-CLIP for comparison
- single-cell/preprocessing - scSTAMP downstream
- single-cell/clustering - scSTAMP per-cluster pseudobulk
- single-cell/markers-annotation - Cell-type-specific targets
- methylation-analysis/methylation-calling - Editing as related to methylation
- epitranscriptomics/m6anet-analysis - Nanopore m6A alternative to DART
