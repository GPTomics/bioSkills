---
name: bio-long-read-sequencing-nanopore-methylation
description: Calls DNA base modifications (5mC, 5hmC, 6mA, 4mC) directly from Oxford Nanopore and PacBio HiFi long reads encoded as MM/ML SAM tags, piles them into per-site bedMethyl with modkit (or pb-CpG-tools for PacBio), and produces phased allele-specific methylation. Covers why methylation is a basecalling decision that cannot be recovered later, the MM/ML tag-drop failure that silently zeroes methylation through alignment, the MM ? vs . no-call semantics, 5mC/5hmC resolution vs bisulfite, modkit's 10th-percentile auto-threshold, and the haplotagged ASM workflow. Use when calling 5mC/5hmC/6mA from a modBAM, generating bedMethyl, preserving methylation tags through alignment, doing allele-specific or differential methylation, or QC-ing a modification BAM.
tool_type: cli
primary_tool: modkit
---

## Version Compatibility

Reference examples tested with: modkit 0.3+, dorado 1.0+, minimap2 2.28+, samtools 1.19+, pb-CpG-tools 2.3+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Inputs that determine what is even possible - record them:
- The basecaller MODIFICATION model (e.g. `5mCG_5hmCG`) fixes which mods can ever be piled up; it must be requested at basecall time and cannot be added later.
- The MM/ML tags must survive every fastq/alignment step or methylation is silently lost.
- modkit auto-estimates the pass threshold from the data (per run); fix it for cross-sample comparisons.

If code throws an error, introspect the installed tool (`modkit pileup --help`, `modkit --help`) and adapt the example to the actual API rather than retrying.

# Nanopore Methylation

**"Call methylation from my long reads"** -> First confirm the MM/ML tags exist and survived alignment, then pile them into per-site bedMethyl - because methylation is a basecalling decision, not something that can be added now.
- CLI: `modkit pileup aligned.bam out.bed --ref ref.fa --cpg --combine-strands`

## The Single Most Important Modern Insight -- Methylation Is a Basecalling Decision, and the Tags Silently Die in Alignment

Two facts gate the entire skill:

1. **If the reads were not basecalled with a modification model, the signal is already gone.** Mods are inferred from raw signal at basecall time (ONT Remora model in Dorado; PacBio kinetics model in jasmine) and written into the unaligned BAM as MM/ML tags. A plain BAM or a FASTQ cannot yield methylation - there is no post-hoc tool. The only fix is re-basecalling from POD5 (`dorado basecaller sup,5mCG_5hmCG pod5/`). The agent's FIRST move is to check the tags exist: `samtools view in.bam | head | grep -o 'MM:Z:[^\t]*'`.
2. **The MM/ML tags silently die in a normal alignment workflow.** `samtools fastq` drops auxiliary tags unless given `-T MM,ML`; minimap2 ignores them unless given `-y`; hard-clipping breaks MM's per-base skip counting unless `-Y` is set. Miss any one and the aligned BAM still sorts, indexes, and looks fine, but `modkit pileup` returns an empty/all-canonical bedMethyl with no error. Use `dorado aligner` (carries tags natively) or `samtools fastq -T MM,ML | minimap2 -y -Y`, and re-grep for MM:Z AFTER alignment.

## End-to-End modBAM Pipeline

```bash
# 1. MODS-BASECALL (from POD5; the only step that can ever produce methylation)
dorado basecaller sup,5mCG_5hmCG pod5/ > calls.bam      # unaligned BAM, has MM/ML

# 2. TAG-PRESERVING ALIGN (route a is simplest)
dorado aligner ref.mmi calls.bam > aligned.bam                                  # a) native
samtools fastq -T MM,ML calls.bam | minimap2 -y -Y -ax lr:hq ref.fa - \
  | samtools sort -o aligned.bam && samtools index aligned.bam                  # b) manual
samtools view aligned.bam | head | grep -q 'MM:Z' && echo 'tags survived'       # verify!

# 3. PILEUP -> bedMethyl (auto-thresholds at the 10th percentile of ML; NOT 0.5)
modkit pileup aligned.bam out.bed --ref ref.fa --cpg --combine-strands
bgzip out.bed && tabix -p bed out.bed.gz

# 4. (optional) DIFFERENTIAL methylation, long-read native
modkit dmr pair -a A.bed.gz -b B.bed.gz --ref ref.fa --regions cpgislands.bed -o dmr.tsv
```

## MM / ML Tag Spec (what the numbers mean)

- `MM:Z` encodes modification positions: `<canonical base><strand><mod code><. or ?>,<skip counts>;`. Mod codes: `m`=5mC, `h`=5hmC, `a`=6mA, `c`=4mC. The `.`/`?` modifier is load-bearing: `.` = skipped bases are implicitly canonical (count toward the unmodified denominator); `?` = skipped bases are no-call/unknown (land in Nnocall, outside the denominator). Misreading `?` as `.` inflates the canonical denominator and deflates methylation.
- `ML:B:C` is a uint8 per call: value N means probability in `[N/256, (N+1)/256)`, so 255 is ~0.998, never exactly 1.0. Do not threshold `== 1.0`.

## bedMethyl Columns (modkit, 18 columns)

Cols 1-9 are tab-delimited BED9; cols 10-18 are space-delimited (a parsing gotcha). The ones that matter:

| Col | Name | Meaning |
|-----|------|---------|
| 10 | Nvalid_cov | Nmod + Ncanonical + Nother_mod (the denominator; this is "coverage" for QC) |
| 11 | percent_modified | (Nmod / Nvalid_cov) * 100 (a percent, 0-100) |
| 12 | Nmod | passing calls of this modification |
| 13 | Ncanonical | passing calls of the canonical base |
| 14 | Nother_mod | passing calls of a different mod on the same base (5hmC in a 5mC row) |
| 16 | Nfail | calls below the pass threshold (excluded from Nvalid_cov) |
| 18 | Nnocall | aligned canonical base with no mod call (e.g. `?`-skipped) |

For count-based DMR (DSS/methylKit) hand over Nmod and Nvalid_cov, never percent_modified.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| ONT 5mC for mammals | `dorado ...sup,5mCG_5hmCG` -> `modkit pileup --cpg --combine-strands` | mammalian 5mC is overwhelmingly CpG |
| Compare ONT to WGBS/array | `modkit pileup --combine-mods` (or `--preset traditional`) | WGBS conflates 5mC+5hmC; combine to match |
| Study 5hmC biology | keep 5mC and 5hmC split; ideally add oxBS/TAB-seq | bisulfite cannot separate them |
| Plants (CHG/CHH) or bacterial 6mA/4mC | all-context model + `--motif` (not `--cpg`) | methylation is not CpG-restricted there |
| Allele-specific methylation / imprinting | phase + haplotag -> `modkit pileup --partition-tag HP` | one read carries SNV phase AND methylation |
| PacBio HiFi 5mC | `ccs --hifi-kinetics` -> `jasmine` -> pb-CpG-tools (or modkit) | primrose is deprecated; Revio does 5mC on-instrument |
| Differential methylation statistics | `modkit dmr` (native) or export to -> methylation-analysis | DSS/methylKit for dispersion modeling |
| RNA modifications (m6A etc.) | -> epitranscriptomics | direct-RNA mods are out of scope here |

## Phased Allele-Specific Methylation

```bash
# Order is strict: align (tags preserved) -> phase+haplotag -> pileup partitioned by HP
# 1-2. Clair3/DeepVariant -> whatshap/longphase phase + haplotag (adds HP:i:1/2) -> haplotype-phasing
modkit pileup aligned.haplotagged.bam asm_out/ --ref ref.fa --cpg --combine-strands --partition-tag HP
# --partition-tag writes one UNCOMPRESSED bedMethyl per HP value into asm_out/, named by the tag
# value (e.g. 1.bed, 2.bed). bgzip + tabix each before dmr, which requires indexed inputs:
bgzip asm_out/1.bed && tabix -p bed asm_out/1.bed.gz
bgzip asm_out/2.bed && tabix -p bed asm_out/2.bed.gz
modkit dmr pair -a asm_out/1.bed.gz -b asm_out/2.bed.gz --ref ref.fa -o asm.tsv
```

Each haplotype gets ~half the coverage, so the per-site 10x floor effectively wants ~20x total. Imprinted loci (one haplotype ~fully methylated) are the canonical positive control.

## Per-Method Failure Modes

### No methylation in a plain BAM
**Trigger:** BAM basecalled without a mods model. **Mechanism:** mods are a basecall-time decision. **Symptom:** no MM:Z tags; modkit returns nothing. **Fix:** re-basecall from POD5 with a mods model; there is no post-hoc tool.

### Tags died in alignment (the #1 silent killer)
**Trigger:** `samtools fastq | minimap2` without `-T MM,ML`/`-y`. **Mechanism:** fastq export and minimap2 drop the tags. **Symptom:** valid aligned BAM, empty bedMethyl, no error. **Fix:** `dorado aligner`, or `samtools fastq -T MM,ML | minimap2 -y -Y`; verify MM:Z after alignment.

### Methylation fraction looks too low
**Trigger:** misreading `?` (no-call) as `.` (canonical). **Mechanism:** unscored bases counted as unmethylated. **Symptom:** deflated percent_modified; large Nnocall. **Fix:** check the MM modifier and Nnocall; `modkit update-tags` to convert styles if needed.

### 5mC understated vs WGBS
**Trigger:** comparing ONT-5mC-only to bisulfite. **Mechanism:** WGBS reads 5mC+5hmC together. **Symptom:** ONT looks lower by the 5hmC fraction. **Fix:** `--combine-mods`/`--preset traditional` to combine before comparing.

### Cross-sample thresholds not comparable
**Trigger:** relying on modkit's auto-threshold per sample. **Mechanism:** the 10th-percentile cut is data-dependent. **Symptom:** sample-specific thresholds confound a DMR. **Fix:** fix a common `--filter-threshold`/`--mod-thresholds` across samples.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Nvalid_cov >= 10 per CpG | field standard | below it, single-site fractions are noisy (~20x total for phased) |
| modkit pass = 10th percentile of ML | modkit docs | discards the lowest-confidence ~10%; improves WGBS concordance |
| R10 ONT vs WGBS r ~ 0.84-0.95 | benchmarks | adequate-depth site-level concordance (R10 > R9) |
| ML 255 ~ 0.998 (not 1.0) | SAM spec | uint8 bin [255/256, 1.0); never test == 1.0 |
| no `--min-coverage` flag on pileup | modkit API | filter bedMethyl on Nvalid_cov post-hoc instead |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Empty bedMethyl / Nvalid_cov 0 | tags dropped or never present | grep MM:Z; re-basecall or re-align preserving tags |
| `modkit pileup --min-coverage` unknown flag | no such flag | filter on Nvalid_cov (col 10) after pileup |
| `modkit extract in.bam out.tsv` errors | needs a subcommand | `modkit extract full` / `modkit extract calls` |
| `modkit dmr` fails on raw .bed | inputs must be indexed | `bgzip` + `tabix -p bed` first |
| Sparse bedMethyl with `--combine-strands` | records lack MN tags (old basecaller or hard-clipped) | basecall with current Dorado and align with `-Y` (no hard-clip) |
| Methylation lower than expected | `?` read as `.`; or 5hmC excluded vs WGBS | check Nnocall; `--combine-mods` to match WGBS |
| primrose not found (PacBio) | deprecated/archived | use `jasmine` (or Revio on-instrument 5mC) |

## References

- Simpson JT, Workman RE, Zuzarte PC, et al. 2017. Detecting DNA cytosine methylation using nanopore sequencing. *Nat Methods* 14:407-410.
- Yuen ZW-S, Srivastava A, Daniel R, et al. 2021. Systematic benchmarking of tools for CpG methylation detection from nanopore sequencing (METEORE). *Nat Commun* 12:3438.
- Tse OYO, Jiang P, Cheng SH, et al. 2021. Genome-wide detection of cytosine methylation by single molecule real-time sequencing. *PNAS* 118(5):e2019768118.
- Cheetham SW, Kindlova M, Ewing AD. 2022. Methylartist: tools for visualizing modified bases from nanopore sequence data. *Bioinformatics* 38(11):3109-3112.
- SAM Optional Fields Specification (SAMtags): MM/ML base-modification tags. samtools/hts-specs.

## Related Skills

- basecalling - The upstream gate: methylation must be requested at basecall time
- long-read-alignment - Carry MM/ML through with `-y -Y` (or use dorado aligner)
- haplotype-phasing - Phase + haplotag the BAM for allele-specific methylation
- clair3-variants - SNVs to phase before allele-specific methylation
- methylation-analysis/dmr-detection - DMR statistics downstream of bedMethyl
- methylation-analysis/methylkit-analysis - methylKit differential methylation
- epitranscriptomics/m6anet-analysis - Direct-RNA m6A (out of scope here)
- workflows/methylation-pipeline - End-to-end methylation pipeline
