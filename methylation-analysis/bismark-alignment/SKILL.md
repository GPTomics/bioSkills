---
name: bio-methylation-bismark-alignment
description: Aligns bisulfite-converted (WGBS, RRBS, PBAT) and enzymatic (EM-seq) short reads to an in-silico C->T/G->A-converted reference with Bismark (Bowtie2 or HISAT2), preparing the genome index, choosing the directional vs non-directional vs PBAT strand flag, deduplicating WGBS/EM-seq (never RRBS), and bounding bisulfite conversion efficiency with unmethylated lambda and methylated pUC19 spike-ins. Covers why the library protocol (not the aligner) decides whether calls are meaningful, why incomplete conversion masquerades as methylation, the 3-letter reduced-complexity mapping bias (50-70% efficiency is normal), and M-bias end-clipping. Use when aligning bisulfite or EM-seq reads, preparing a bisulfite genome, choosing the strand flag, or diagnosing low mapping efficiency. For methylation extraction see methylation-calling; for long-read MM/ML modification calling see long-read-sequencing/nanopore-methylation.
tool_type: cli
primary_tool: Bismark
---

## Version Compatibility

Reference examples tested with: Bismark 0.24+, Bowtie2 2.5+, HISAT2 2.2+, Trim Galore 0.6.10+, samtools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags and defaults

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The genome build and the aligner backend ARE the versions that matter. The bisulfite index is built once per genome FASTA with a specific backend (`--bowtie2` vs `--hisat2`); the index must match the backend used at alignment time, and the FASTA build (hg38 vs T2T-CHM13) fixes every downstream coordinate. EM-seq uses the identical aligners and flags as bisulfite - only the upstream chemistry and the coverage/efficiency expectations change.

# Bismark Alignment

**"Align my bisulfite sequencing reads"** -> Confirm the library type to pick the strand flag, trim the chemistry-specific artifacts, then map the C->T-converted reads to a C->T/G->A-converted reference - because the protocol and conversion, not the aligner, decide whether the calls mean anything.
- CLI: `bismark_genome_preparation --bowtie2 genome/` then `bismark --genome genome/ -1 R1.fq.gz -2 R2.fq.gz -o out/`

Scope: short-read bisulfite (WGBS/RRBS/PBAT) and enzymatic (EM-seq) alignment, the genome index, the strand/library flag, deduplication, and conversion QC. Methylation extraction from the BAM (XM tag, MethylDackel, cytosine reports) -> methylation-calling. Per-CpG/DMR statistics -> differential-cpg-testing, dmr-detection. Long-read native MM/ML modification calling -> long-read-sequencing/nanopore-methylation. Adapter trimming mechanics -> read-qc/adapter-trimming.

## The Single Most Important Modern Insight -- Methylation Is Never Sequenced; the Survivors of a Deamination Assay Are

A bisulfite (or EM-seq) run never reads methylation. It reads which cytosines SURVIVED deamination, against a 3-letter genome deliberately depleted of cytosines, as a C-vs-T choice. Every methylation call is two stacked conditional bets, and both fail silently:

1. **Conversion went to completion in BOTH directions.** An unmethylated C that escapes deamination survives as C and is called methylated -> false HYPER-methylation (under-conversion, the dominant fear). A genuinely methylated C deaminated anyway reads T -> false HYPO-methylation (over-conversion). Neither error is visible in the BAM or the mapping rate - only spike-in controls see them, and one control sees only one direction (lambda for under, pUC19 for over).
2. **The read mapped to the right place despite throwing its cytosines away.** The 3-letter alphabet collapses uniqueness, so a wrong library-type flag (PBAT or non-directional run as directional) silently drops half to nearly all reads, and a C/T SNP masquerades as an unmethylated CpG with no alignment penalty.

Organize the work around defending these two bets - chemistry control (both directions) and library/strand correctness - not around listing `bismark` flags. The aligner reports a clean, sorted, indexed BAM whether conversion failed or half the reads went unmapped.

## Why 3-Letter Mapping Is Hard (and Why 50-70% Is Normal)

After conversion, unmethylated Cs become Ts, so the read/genome alphabet collapses toward {A,G,T}. A normal aligner would penalize every C->T as a mismatch, so bisulfite aligners convert all Cs to T in BOTH the reads AND the reference, map in the reduced alphabet, then recover methylation by comparing the original read to the original reference. Bismark builds two converted indices (C->T for OT/CTOT, G->A for OB/CTOB) and aligns each read against both. Reduced complexity means more multi-mapping and a lower mapping efficiency (~50-70% for WGBS vs >95% for ordinary DNA) - this is expected, not a bug. The same collapse means a sample CpG->TpG variant aligns with no extra mismatch and is scored as an unmethylated CpG: methylation at a C/T-polymorphic site is a hypothesis until SNP-aware (Bis-SNP, BISCUIT).

## The Four Strands and the Library-Type Flag

Bisulfite PCR generates four strand species: OT (original top), OB (original bottom), CTOT (complement of OT), CTOB (complement of OB). The library protocol decides which exist, and the flag must match or reads vanish silently:

| Library | Strands sequenced | Bismark flag | Dedup? | Trim Galore special-case |
|---------|-------------------|--------------|--------|--------------------------|
| WGBS (directional) | OT, OB | (default) | YES | M-bias end-clip |
| EM-seq (directional) | OT, OB | (default) | YES | M-bias end-clip (gentler) |
| RRBS | OT, OB | (default) | NO | `--rrbs` (MspI fill-in) |
| PBAT / scBS-seq | CTOT, CTOB | `--pbat` | usually NO | aggressive 5' clip (random priming) |
| non-directional | all four | `--non_directional` | YES | M-bias end-clip (Trim Galore `--non_directional` is RRBS-only, needs `--rrbs`) |

PBAT does bisulfite conversion FIRST then tags by random priming, so its reads originate from CTOT/CTOB - the OPPOSITE of directional. PBAT needs `--pbat` for strand reasons; it is unrelated to RRBS. A non-directional library run as directional silently loses ~half its reads; PBAT run as directional maps near zero.

## Tool Taxonomy

| Tool | Citation | Strategy | When |
|------|----------|----------|------|
| Bismark | Krueger & Andrews 2011 *Bioinformatics* 27:1571 | 3-letter, Bowtie2/HISAT2 backend | de-facto standard; self-contained (index + align + dedup + extractor); teach this |
| bwa-meth | Pedersen 2014 arXiv:1401.1129 | 3-letter, BWA-MEM | lean clinical/cfDNA; handles indels/clipping; pairs with MethylDackel for calling |
| BISCUIT | Zhou 2024 *Nucleic Acids Res* 52:e32 | 3-letter, BWA-derived | when SNPs / allele-specific methylation are needed alongside (joint genetic+epigenetic) |
| gemBS | Merkel 2019 *Bioinformatics* 35:737 | 3-letter, GEM3 | population-scale; the ENCODE WGBS pipeline mapper |
| abismal / methylpy | de Sena Brandine & Smith 2021 *NAR Genom Bioinform* 3:lqab115 | 2-letter (purine/pyrimidine) | memory-constrained, large cohorts |

All produce a BAM whose methylation is recovered by a SEPARATE caller (Bismark extractor, MethylDackel, or the tool's own). Alignment and calling are two steps.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Standard WGBS or EM-seq, mammalian | Bismark default (directional) + dedup | OT/OB only; the common case |
| RRBS | `trim_galore --rrbs` then Bismark default, NO dedup | MspI fixed ends look like (but are not) PCR duplicates |
| PBAT / scBS-seq | `bismark --pbat` | reads come from CTOT/CTOB, not OT/OB |
| Non-directional library | `bismark --non_directional` | all four strands present; default loses half |
| Precious low-input (cfDNA / FFPE / single-cell) | prefer EM-seq or TAPS upstream | bisulfite degrades 84-96% of input; same aligners apply |
| Need SNPs / allele-specific methylation | -> bwa-meth + Bis-SNP, or BISCUIT | C/T SNPs masquerade as methylation in 3-letter space |
| Large mammalian genome, low RAM | `bismark --hisat2` (index must match) | HISAT2 backend is lighter than Bowtie2 |
| Extract per-CpG methylation from the BAM | -> methylation-calling | this skill stops at the deduplicated, M-bias-clipped BAM |
| Long-read ONT/PacBio modBAM (MM/ML tags) | -> long-read-sequencing/nanopore-methylation | native modification calling, not bisulfite |

## Prepare the Genome Index

**Goal:** Build the bisulfite-converted index once per genome, with the backend that alignment will use.

**Approach:** Place the reference FASTA(s) in a folder, run `bismark_genome_preparation` with the chosen backend; it writes `Bisulfite_Genome/` containing the C->T and G->A converted indices.

```bash
bismark_genome_preparation --bowtie2 genome/   # or --hisat2 for large genomes, lower RAM
# genome/ holds the FASTA (e.g. hg38.fa); writes genome/Bisulfite_Genome/{CT_conversion,GA_conversion}
# The backend chosen here MUST match the bismark alignment backend below.
```

## Trim First, with the Chemistry-Specific Flag

**Goal:** Remove adapters and the library-specific end artifacts before alignment so they do not become spurious methylation calls.

**Approach:** Run Trim Galore (Cutadapt wrapper). Add `--rrbs` for RRBS (clips the MspI end-repair fill-in), `--non_directional` for non-directional, or extra 5' clipping for PBAT. Bismark itself does not trim. Mechanics live in read-qc/adapter-trimming.

```bash
trim_galore --paired R1.fq.gz R2.fq.gz                 # WGBS / EM-seq (auto-detect adapter, -q 20)
trim_galore --rrbs --paired R1.fq.gz R2.fq.gz          # RRBS: extra 2 bp off 3' R1 (+ 5' R2) = MspI fill-in
trim_galore --clip_r2 6 --paired R1.fq.gz R2.fq.gz     # PBAT/scBS: random-priming bias at 5' (amount from M-bias)
```

## Align

```bash
bismark --genome genome/ -1 R1_val_1.fq.gz -2 R2_val_2.fq.gz \
    --bowtie2 \         # must match the index backend; --hisat2 if prepared that way
    --parallel 4 \      # instances PER direction; total threads scale up several-fold per instance
    -o out/             # writes *_bismark_bt2_pe.bam + *_PE_report.txt (mapping efficiency, %meth per context)
# Add --pbat for PBAT/scBS, or --non_directional for non-directional libraries (NOT both).
```

## Deduplicate (WGBS/EM-seq Only)

**Goal:** Remove PCR duplicates from random-fragmentation libraries, while leaving RRBS untouched.

**Approach:** `deduplicate_bismark` removes reads sharing mapping coordinate + strand. Run it on the by-name (unsorted) Bismark BAM, before extraction. For RRBS, SKIP it: every fragment starts at an MspI cut site, so identical coordinates are biologically distinct molecules, not PCR copies (apparent duplication ~90-95% is real data).

```bash
deduplicate_bismark --paired --bam out/sample_R1_bismark_bt2_pe.bam   # WGBS/EM-seq ONLY
# RRBS: do NOT run this. UMI-tagged RRBS can dedup by UMI+coordinate; optical dups can still be removed.
samtools sort out/sample_R1_bismark_bt2_pe.deduplicated.bam -o out/sample.sorted.bam   # IGV/downstream
samtools index out/sample.sorted.bam
```

## Conversion QC: Both Directions, and the Spike-In Is an Optimistic Floor

**Goal:** Bound both conversion error directions before believing any methylation level.

**Approach:** Spike unmethylated lambda phage (measures under-conversion -> false hyper) AND CpG-methylated pUC19 (measures over-conversion -> false hypo). Align each spike-in genome separately and read off context methylation. With no spike-in, sample CHH methylation is a weak fallback (somatic tissue only; confounded in ESCs/neurons/plants).

```bash
bismark_genome_preparation --bowtie2 lambda/   # lambda: residual %meth = non-conversion rate (target <=1%)
bismark --genome lambda/ -1 R1.fq.gz -2 R2.fq.gz -o lambda_qc/
# pUC19 (CpG-methylated): fraction of CpGs called UNmethylated = over-conversion (expect ~96-98% methylated)
```

Spike-ins are naked, fully accessible DNA that denature completely, so their conversion is an OPTIMISTIC upper bound. Real genomic conversion is region-dependent: GC-rich CpG islands and structured regions denature less, under-convert more, and inflate apparent methylation exactly where the biology is. Treat the spike-in number as a floor; a rising per-GC-bin CHH rate flags local under-conversion.

## Per-Method Failure Modes

### PBAT or non-directional run with the default flag
**Trigger:** running PBAT/scBS or a non-directional library without `--pbat`/`--non_directional`. **Mechanism:** PBAT reads come from CTOT/CTOB and non-directional from all four strands, but the default tries only OT/OB. **Symptom:** near-zero (PBAT) or ~halved (non-directional) mapping efficiency on a clean-looking run. **Fix:** confirm the kit/protocol directionality, pass the matching flag; never reach for `-N 1` first.

### Incomplete conversion read as methylation
**Trigger:** no conversion control, or only a lambda (under-conversion) control. **Mechanism:** an unmethylated C surviving deamination is indistinguishable from real 5mC. **Symptom:** globally elevated methylation, worst in GC-rich CpG islands. **Fix:** report BOTH a lambda non-conversion rate (<=1%) and a pUC19 over-conversion rate; add per-GC CHH as an internal check.

### RRBS deduplicated by coordinate
**Trigger:** running `deduplicate_bismark` on RRBS. **Mechanism:** MspI cuts give every fragment a fixed start, so distinct molecules share coordinates. **Symptom:** ~90-95% of reads discarded, coverage decimated. **Fix:** skip coordinate dedup for RRBS; use UMIs if dedup is required.

### MspI fill-in not trimmed
**Trigger:** RRBS aligned without `trim_galore --rrbs`. **Mechanism:** end-repair fills MspI overhangs with unmethylated dCTP, creating artificial cytosines at fragment ends. **Symptom:** artificial hypomethylation clustered at MspI sites. **Fix:** `trim_galore --rrbs`; Bismark aligns RRBS fine but does NOT fix this trimming artifact.

### M-bias not clipped before calling
**Trigger:** calling methylation off raw read ends. **Mechanism:** end-repair fills 5' overhangs with unmethylated dCTP, worst at the start of R2. **Symptom:** an M-bias plot (methylation vs read position) shows a dip/spike at the ends instead of a flat line. **Fix:** read the M-bias plot, clip the affected ends; extraction `--ignore`/`--clip` mechanics live in methylation-calling.

### Index/backend mismatch
**Trigger:** index prepared with `--bowtie2`, alignment run with `--hisat2` (or vice versa). **Mechanism:** the two backends use incompatible converted indices. **Symptom:** Bismark errors or fails to find the index. **Fix:** prepare and align with the same backend.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Lambda non-conversion <=1% | manufacturer spec; field standard | residual apparent methylation on unmethylated spike-in = false-positive floor (EM-seq v2 ~<=0.5%) |
| pUC19 ~96-98% methylated | manufacturer spec | bounds over-conversion -> false hypo; lambda alone cannot see this direction |
| WGBS mapping efficiency ~50-70% | Krueger & Andrews 2011; 3-letter complexity | reduced alphabet costs uniqueness; below this, diagnose (library flag > trimming > reference > biology) |
| EM-seq mapping efficiency typically higher | Vaisvila 2021 | no chemical fragmentation -> flatter coverage; WGBS expectations are too pessimistic |
| Bisulfite degrades 84-96% of input | Grunau 2001 *Nucleic Acids Res* 29:e65 | only ~4-16% of molecules survive intact; the reason low-input fails and EM-seq/TAPS exist |
| `-N` = 0 (seed mismatches) | Bismark manual | default; `-N 1` raises sensitivity AND mis-mapping - last resort, not the low-mapping fix |
| `--rrbs` clips 2 bp | Trim Galore guide | the MspI end-repair fill-in length; confirm on the installed version |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Near-zero mapping efficiency | PBAT run as directional | add `--pbat` |
| ~Half the reads unmapped | non-directional run as directional | add `--non_directional` |
| RRBS loses ~90% of reads | deduplicated by coordinate | skip `deduplicate_bismark` for RRBS |
| Globally high methylation | incomplete conversion (no/one-sided control) | lambda + pUC19 spike-ins; check per-GC CHH |
| Artificial hypomethylation at MspI sites | `--rrbs` trimming omitted | `trim_galore --rrbs` |
| FastQC per-base content / GC FAIL | expected for converted libraries (C depleted) | not a defect; do not "fix" a healthy bisulfite library |
| 0% sites at C/T variants | C/T SNP read as unmethylated CpG | SNP-aware calling (Bis-SNP/BISCUIT) or mask known C/T SNPs |
| Bismark cannot find the index | backend mismatch with genome prep | re-prep or align with the matching `--bowtie2`/`--hisat2` |
| Output named "5mC" | standard BS and EM-seq report 5mC+5hmC summed | label the sum; oxBS/TAB pairing is needed to separate (see methylation-calling) |

## References

- Krueger F, Andrews SR. 2011. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. *Bioinformatics* 27:1571-1572.
- Vaisvila R, Ponnaluri VKC, Sun Z, et al. 2021. Enzymatic methyl sequencing detects DNA methylation at single-base resolution from picograms of DNA. *Genome Res* 31:1280-1289.
- Meissner A, Gnirke A, Bell GW, et al. 2005. Reduced representation bisulfite sequencing for comparative high-resolution DNA methylation analysis. *Nucleic Acids Res* 33:5868-5877.
- Miura F, Enomoto Y, Dairiki R, Ito T. 2012. Amplification-free whole-genome bisulfite sequencing by post-bisulfite adaptor tagging. *Nucleic Acids Res* 40:e136.
- Hansen KD, Langmead B, Irizarry RA. 2012. BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions. *Genome Biol* 13:R83.
- Grunau C, Clark SJ, Rosenthal A. 2001. Bisulfite genomic sequencing: systematic investigation of critical experimental parameters. *Nucleic Acids Res* 29:e65.
- Pedersen BS, Eyring K, De S, Yang IV, Schwartz DA. 2014. Fast and accurate alignment of long bisulfite-seq reads. arXiv:1401.1129.
- Zhou W, Johnson BK, Morrison J, et al. 2024. BISCUIT: an efficient, standards-compliant tool suite for simultaneous genetic and epigenetic inference in bulk and single-cell studies. *Nucleic Acids Res* 52:e32.
- Merkel A, Fernandez-Callejo M, Casals E, et al. 2019. gemBS: high throughput processing for DNA methylation data from bisulfite sequencing. *Bioinformatics* 35:737-742.
- de Sena Brandine G, Smith AD. 2021. Fast and memory-efficient mapping of short bisulfite sequencing reads using a two-letter alphabet. *NAR Genom Bioinform* 3:lqab115.

## Related Skills

- methylation-calling - Extract per-CpG methylation from the aligned BAM
- methylkit-analysis - Downstream import, filtering, normalization
- read-qc/adapter-trimming - Trim Galore before Bismark (RRBS/PBAT handling)
- read-qc/quality-reports - FastQC (expect per-base C-depletion FAIL on converted libraries)
- alignment-files/sam-bam-basics - BAM manipulation after alignment
- sequence-io/read-sequences - FASTQ handling before alignment
- long-read-sequencing/nanopore-methylation - Native long-read MM/ML modification calling (out of scope here)
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
