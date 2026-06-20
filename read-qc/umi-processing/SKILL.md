---
name: bio-read-qc-umi-processing
description: Extracts UMIs and collapses reads to original molecules with umi_tools (directional dedup) or builds error-corrected single-strand/duplex consensus reads with fgbio. Use when the library has UMIs and accurate molecule counting or below-sequencer-floor error correction is needed - single-cell, low-input RNA-seq, targeted panels, and ctDNA/liquid-biopsy rare-variant detection. For UMI extraction during QC use fastp-workflow; do not dedup non-UMI bulk RNA-seq.
tool_type: cli
primary_tool: umi_tools
---

## Version Compatibility

Reference examples tested with: umi_tools 1.1+, fgbio 2.1+, samtools 1.19+, STAR 2.7+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# UMI Processing -- count original molecules, or build a consensus below the error floor

Collapse PCR/optical duplicates by (coordinate + UMI) with umi_tools, or call error-corrected consensus reads with fgbio.

**"Deduplicate reads using UMIs"** -> Extract the UMI before alignment, then group reads by UMI + mapping position after alignment to count original molecules.
- CLI: `umi_tools extract` -> align -> `umi_tools dedup` (molecule counting)
- CLI: `fgbio GroupReadsByUmi` -> `fgbio CallMolecularConsensusReads`/`CallDuplexConsensusReads` (error correction)

Scope: this skill OWNS UMI extraction, dedup, and consensus calling. UMI extraction during QC -> read-qc/fastp-workflow. Single-cell matrices -> single-cell/preprocessing. Non-UMI DNA coordinate dedup -> alignment-files/duplicate-handling. OUT OF SCOPE: non-UMI bulk RNA-seq (do NOT dedup it -- read-qc/rnaseq-qc).

## The Single Most Important Modern Insight

1. **UMIs resolve the PCR-vs-biological duplicate confound that coordinates alone cannot, by collapsing on (coordinate + UMI) instead of coordinate -- and this forces a hard pipeline order: extract UMI on the FASTQ, align, THEN dedup.** Two reads at the same coordinate are the same molecule only if they also share a UMI; two independent molecules at one coordinate carry different UMIs. The confound dominates at high coverage, high expression, amplicon (every molecule shares the same primer-defined ends), and low input. Dedup cannot run before alignment because duplicate identity needs mapping COORDINATES; and `extract` must run before alignment so the aligner does not try to map the UMI bases as genomic sequence (extract moves the UMI into the read name / RX tag).

2. **umi_tools' DIRECTIONAL method (default) folds UMI errors back into their parent via a count-gradient rule; naive exact-UMI collapse OVER-counts.** Sequencing/PCR errors inside the UMI mutate a true UMI into a 1-off neighbor that looks like a new molecule. Directional builds a directed graph where an edge a->b exists when they are within edit distance 1 AND n_a >= 2*n_b - 1 (the parent is at least ~twice the error child, because errors are rarer than originals), then collapses each network to one molecule. This is why directional beats `cluster` (single-linkage over-merges, under-counts) and `unique` (no error model, over-counts).

3. **UMI-tools COUNTS molecules; fgbio builds a CONSENSUS read to push the error rate BELOW the sequencer floor -- and only DUPLEX consensus reaches the ctDNA/MRD floor.** Single-strand consensus (CallMolecularConsensusReads) votes within one strand's family and roughly halves errors, but cannot catch a lesion fixed into the molecule before the first copy (oxidative 8-oxo-G, C>T deamination). Duplex consensus (CallDuplexConsensusReads) keeps a base only where BOTH original strands agree -- a real mutation is on both strands, an artifact almost never -- reaching <1e-7 error for sub-0.1% VAF detection, at the cost of ~2x raw reads (families missing one strand are discarded).

Bridges: do NOT dedup non-UMI bulk RNA-seq (high-expression genes make genuine duplicate coordinates; read-qc/rnaseq-qc). CellRanger/STARsolo ALREADY UMI-collapse and emit a final matrix -- do not re-dedup their output. Deep amplicon needs LONGER UMIs because every molecule shares coordinates, so the UMI alone must separate them (4^L space; collisions under-count).

## Tool Taxonomy

| Tool / command | Role | When |
|----------------|------|------|
| umi_tools extract | Move UMI from read into the header (FASTQ stage) | Inline UMIs before alignment |
| umi_tools dedup | Collapse to one read per (coord + UMI) via directional | Molecule counting (bulk, targeted) |
| umi_tools count | Emit a gene x cell molecule matrix | Single-cell from a tagged raw BAM |
| umi_tools group | Tag reads with UG (group id) + BX (representative UMI), no dedup | Inspect grouping / feed consensus |
| fgbio GroupReadsByUmi | Group reads into source-molecule families (MI tag) | First step of consensus calling |
| fgbio CallMolecularConsensusReads | Single-strand consensus | Moderate-VAF error correction |
| fgbio CallDuplexConsensusReads | Duplex consensus (both strands agree) | ctDNA / MRD sub-0.1% VAF |
| fgbio FilterConsensusReads | Filter/mask untrustworthy consensus bases | Mandatory after consensus calling |
| fastp --umi | Extract only (no dedup) | UMI extraction folded into QC (route OUT) |

## Decision Tree by Scenario

| Goal | Use | Why |
|------|-----|-----|
| Count molecules (bulk/targeted RNA or DNA) | umi_tools dedup --method directional | Models UMI errors; the standard |
| Single-cell molecule matrix | umi_tools count (tagged raw BAM) or the aligner's own collapse | per-cell + per-gene |
| Already have a CellRanger/STARsolo matrix | nothing | It is already UMI-deduplicated |
| Moderate-VAF somatic error correction | fgbio single-strand consensus | Halves errors |
| ctDNA / MRD sub-0.1% VAF | fgbio duplex consensus + FilterConsensusReads | Below the single-strand floor |
| Non-UMI bulk RNA-seq | do NOT dedup | Duplicate coordinates are biological |

Default when uncertain: umi_tools directional dedup for counting; fgbio duplex for ctDNA.

## Extraction (FASTQ stage, before alignment)

`--bc-pattern` alphabet (string method): N = UMI base (extracted to the read name), C = cell barcode (extracted), X = a fixed/known base REATTACHED to the read (not discarded). True discard uses the regex method's `(?P<discard_N>...)` group, shown below.

```bash
# Inline 8 nt UMI at the start of R1
umi_tools extract --stdin=R1.fq.gz --read2-in=R2.fq.gz \
    --stdout=R1_umi.fq.gz --read2-out=R2_umi.fq.gz --bc-pattern=NNNNNNNN

# 10x 3' v3: 16 nt cell barcode + 12 nt UMI on R1
umi_tools extract --stdin=R1.fq.gz --read2-in=R2.fq.gz \
    --stdout=R1_umi.fq.gz --read2-out=R2_umi.fq.gz \
    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN

# Variable-position UMI with an anchor (regex method)
umi_tools extract --extract-method=regex --stdin=R1.fq.gz --stdout=R1_umi.fq.gz \
    --bc-pattern='(?P<umi_1>.{8})ATGC(?P<discard_1>.{4})'

# fgbio reads structure (M=UMI, T=template, C=cell, B=sample barcode, S=skip)
fgbio FastqToBam --input R1.fq.gz R2.fq.gz --read-structures 8M+T +T \
    --sample S1 --library L1 --output unmapped.bam      # UMI -> RX tag
```

## umi_tools dedup (molecule counting)

```bash
samtools sort -o sorted.bam aligned.bam && samtools index sorted.bam

# Directional (default), paired, with the diagnostic edit-distance stats
umi_tools dedup -I sorted.bam -S dedup.bam --paired --output-stats=stats

# Single-cell from a RAW aligned BAM whose CB/UB are in tags (NOT a CellRanger BAM)
umi_tools count -I tagged.bam -S counts.tsv \
    --per-gene --gene-tag=XT --per-cell --cell-tag=CB \
    --umi-tag=UB --extract-umi-method=tag
```

| Method | Behavior | Verdict |
|--------|----------|---------|
| directional (default) | Count-gradient graph (n_a >= 2n_b-1); folds UMI errors into parent | Best; the default |
| adjacency | Resolve each component by abundance, one edge out | Reasonable |
| cluster | One molecule per connected component (single-linkage) | Over-merges, under-counts |
| unique | Exact UMI only, no error model | Over-counts; only PCR-free/high-diversity |
| percentile | Drop UMIs below 1% of mean count | Crude denoiser |

`--edit-distance-threshold` default 1; `--output-stats` writes the edit-distance file (observed-vs-null confirms UMI errors were collapsed); `umi_tools group --output-bam` writes UG + BX tags without deduplicating.

## fgbio consensus (error correction)

```bash
# Group reads into source-molecule families (writes MI tag from raw RX)
fgbio GroupReadsByUmi --input mapped.bam --output grouped.bam --strategy adjacency --edits 1

# Single-strand consensus (--min-reads required; raise to >=2-3 when error correction matters)
fgbio CallMolecularConsensusReads --input grouped.bam --output consensus.bam --min-reads 3

# Duplex consensus for ctDNA: group with the paired strategy, then call duplex
fgbio GroupReadsByUmi --input mapped.bam --output grouped.bam --strategy paired --edits 1
fgbio CallDuplexConsensusReads --input grouped.bam --output duplex.bam --min-reads 2 1 1

# Mandatory final step: filter/mask untrustworthy consensus bases
fgbio FilterConsensusReads --input duplex.bam --output filtered.bam --ref ref.fa \
    --min-reads 2 1 1 --max-base-error-rate 0.1 --min-base-quality 40 --max-no-calls 0.2
```

GroupReadsByUmi `--strategy`: identity (exact), edit (cluster by edits), adjacency (umi_tools directional port), paired (DUPLEX -- a read with UMI A-B is the opposite strand of one with B-A, tagged MI .../A and .../B). The consensus pipeline aligns, groups, calls consensus, then RE-aligns the consensus reads (the sequence changed). RX = raw UMI, MI = molecular id (SAM tags).

## Saturation and collision

A fully-random L-mer UMI has 4^L sequences (L=8 -> 65,536; L=12 -> ~16.8M). When the molecules at a locus approach the usable space, independent molecules COLLIDE on the same UMI and are under-counted. For bulk/RNA the key is coordinate+UMI, so the space is 4^L per coordinate and collisions are rare; for AMPLICON every molecule shares coordinates, so the UMI alone separates them and deep panels need longer UMIs (AmpUMI sizes this). UMIs do NOT fix capture/ligation bias upstream of tagging, errors before UMI attachment (only duplex does), or low library complexity.

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Re-running dedup on CellRanger output | CellRanger/STARsolo already UMI-collapse | Use their matrix as-is; do not re-dedup |
| Deduped a non-UMI bulk RNA-seq BAM | Coordinate dups are biological there | Do not dedup; report duplication as a diagnostic |
| Molecule count too high | `--method unique` (no UMI error model) | Use directional (default) |
| Aligner soft-clips/mismaps the UMI | Dedup attempted before extract, or UMI left in read | extract first; UMI must leave the aligned sequence |
| Amplicon molecules under-counted | UMI too short -> collisions at shared coordinates | Use a longer UMI; size with AmpUMI |
| Duplex yields few consensus reads | Many families missing one strand | Expected; duplex needs ~2x raw reads |
| Consensus BAM still noisy | Skipped FilterConsensusReads | Always filter/mask after calling consensus |

## References

Smith T, Heger A, Sudbery I. 2017. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research 27(3):491-499.
Liu D. 2019. Algorithms for efficiently collapsing reads with Unique Molecular Identifiers. PeerJ 7:e8275.
Islam S, Zeisel A, Joost S, et al. 2014. Quantitative single-cell RNA-seq with unique molecular identifiers. Nature Methods 11(2):163-166.
Schmitt MW, Kennedy SR, Salk JJ, et al. 2012. Detection of ultra-rare mutations by next-generation sequencing. PNAS 109(36):14508-14513.
Kennedy SR, Schmitt MW, Fox EJ, et al. 2014. Detecting ultralow-frequency mutations by Duplex Sequencing. Nature Protocols 9(11):2586-2606.
Clement K, Farouni R, Bauer DE, Pinello L. 2018. AmpUMI: design and analysis of unique molecular identifiers for deep amplicon sequencing. Bioinformatics 34(13):i202-i210.

## Related Skills

read-qc/fastp-workflow - UMI extraction folded into preprocessing
read-qc/rnaseq-qc - Why non-UMI bulk RNA-seq must NOT be deduplicated
alignment-files/duplicate-handling - Coordinate dedup for non-UMI DNA
single-cell/preprocessing - scRNA-seq UMI matrices and downstream
liquid-biopsy/ctdna-mutation-detection - Duplex consensus for rare-variant detection
