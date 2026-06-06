---
name: bio-long-read-sequencing-structural-variants
description: Detects structural variants (deletions, insertions, inversions, duplications, translocations) from Oxford Nanopore and PacBio long-read alignments with Sniffles2, cuteSV, SVIM, and assembly-based callers, joint-genotypes cohorts via the Sniffles2 .snf workflow, and benchmarks with Truvari against GIAB. Covers why an SV call is a representation artifact (the tandem-repeat BED, aligner, and Truvari params set precision/recall as much as the caller), the cuteSV per-platform parameter trap, soft-clipped supplementary alignments as the SV substrate, and the somatic/mosaic boundary to Severus/nanomonsv. Use when calling germline or somatic SVs from ONT/HiFi reads, joint-genotyping a cohort, choosing or tuning an SV caller, or benchmarking SV calls.
tool_type: cli
primary_tool: sniffles
---

## Version Compatibility

Reference examples tested with: Sniffles 2.2+, cuteSV 2.1+, minimap2 2.28+, samtools 1.19+, truvari 4.0+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Results depend on inputs that outlive the binary version - record them:
- The reference-matched tandem-repeat BED supplied to the caller (Sniffles `--tandem-repeats`) drives the FP rate in repeats more than any other setting. Record which TR BED was used.
- Benchmark numbers depend on the region set + TR handling + Truvari params; record all three.
- cuteSV parameters are platform-specific (ONT vs HiFi vs CLR); the defaults are not platform-appropriate.

If code throws an error, introspect the installed tool (`sniffles --help`, `cuteSV --help`) and adapt the example to the actual API rather than retrying.

# Long-Read Structural Variants

**"Find structural variants in my long reads"** -> Map with the SV-ready preset (soft-clipped supplementaries), call with a TR-aware caller, and benchmark stating the region set and Truvari params.
- CLI: `sniffles --input aln.bam --vcf svs.vcf --reference ref.fa --tandem-repeats TR.bed`

Long reads are the killer app for SVs: a single read spans the breakpoint (within-read CIGAR or split alignment) and resolves repeats short reads cannot. By convention SV = >=50 bp; the 30-100 bp range is a VNTR-dominated gray zone where callers disagree most.

## The Single Most Important Modern Insight -- An SV Call Is a Representation Artifact as Much as a Biological Fact

In tandem repeats and segmental duplications, the same biological event has many valid VCF encodings - a deletion can be written as the reciprocal insertion on the other allele, and a VNTR expansion's breakpoints slide freely across repeat units. Consequently:

1. **The tandem-repeat BED, the aligner, and the Truvari parameters decide precision/recall as much as the caller does.** A claim like "caller X has F1 0.95" is meaningless without also stating the region set, the TR BED supplied to the caller, and the Truvari params - change any one and the number moves more than the gap between callers.
2. **Without a TR BED, one event fragments into several false-positive calls** with inconsistent breakpoints. `--tandem-repeats` makes clustering repeat-aware (widening the merge window inside annotated TRs) - the single biggest FP-reduction lever, not a nicety.
3. **`truvari refine` exists precisely to re-harmonize representations** within TR regions; benchmarking TR-dense regions without it systematically understates recall.

## Caller Taxonomy

| Tool | Regime | Best for | Citation |
|------|--------|----------|----------|
| Sniffles2 | germline + population + mosaic | the default germline workhorse; cohort joint genotyping; .snf merge | Smolka 2024 *Nat Biotechnol* 42:1571 |
| cuteSV | germline | high sensitivity, speed; per-platform tuning required | Jiang 2020 *Genome Biol* 21:189 |
| SVIM | germline | scores (not hard-filters) SVs; good INS detection | Heller 2019 *Bioinformatics* 35:2907 |
| pbsv | germline (PacBio) | two-step discover->call; official PacBio tool | PacBio (no journal paper) |
| NanoVar | germline, low-depth | 4-8x ONT clinical | Tham 2020 *Genome Biol* 21:56 |
| dipcall / SVIM-asm / PAV | assembly-based germline | most accurate single sample with phased HiFi; truth-set generation | Li 2018; Heller 2021; Ebert 2021 |
| Severus | somatic (tumor-normal) | cancer T/N, complex/subclonal | Keskus 2026 *Nat Biotechnol* |
| nanomonsv | somatic (tumor-normal) | precise somatic breakpoints, MEI | Shiraishi 2023 *NAR* 51:e74 |
| SVision-pro | de novo + somatic, complex | resolving nested CSVs | Wang 2025 *Nat Biotechnol* 43:181 |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Single ONT/HiFi germline sample | Sniffles2 + `--tandem-repeats` | TR-aware, auto support, fast |
| Cohort germline | Sniffles2 per-sample `.snf` -> merge | re-genotypes from raw signal; true joint genotypes |
| Maximum sensitivity / speed | cuteSV with the platform-matched param set | per-platform tuning is mandatory |
| Phased HiFi, want best per-sample accuracy | assembly-based (dipcall/SVIM-asm) -> hifi-assembly | resolves the alt haplotype directly |
| Tumor-normal somatic SVs | Severus or nanomonsv | paired callers; Sniffles `--mosaic` is single-sample only |
| Low-VAF mosaic in one sample | Sniffles2 `--mosaic` | lowers support, reports VAF (not a T/N caller) |
| Low coverage (4-8x) | NanoVar | designed for low-depth clinical |
| Benchmarking | Truvari (+`refine`) vs GIAB Tier1/CMRG | the field standard; state region + params |

## Alignment for SV Calling

Map with minimap2 (the modern default; NGMLR is a higher-precision/slower legacy niche for Sniffles). Use the platform preset and keep soft-clipped supplementary alignments - split-read callers reconstruct breakpoints from the clipped sequence on those records.

```bash
minimap2 -ax map-ont --MD -Y ref.fa ont.fq.gz | samtools sort -o aln.bam && samtools index aln.bam
#   -Y keeps SEQ on supplementaries (the SV substrate); --MD for cuteSV; map-hifi/map-pb for PacBio
```

## Sniffles2 - germline and the .snf population workflow

```bash
# Single sample (always supply --reference for INS sequence and --tandem-repeats for repeats)
sniffles --input aln.bam --vcf svs.vcf --reference ref.fa --tandem-repeats human_GRCh38_TR.bed

# Cohort: per-sample .snf signature index, then merge + joint-genotype
sniffles --input s1.bam --snf s1.snf --reference ref.fa --tandem-repeats TR.bed
sniffles --input s2.bam --snf s2.snf --reference ref.fa --tandem-repeats TR.bed
sniffles --input s1.snf s2.snf --vcf cohort.vcf --reference ref.fa

# Force-call / regenotype a known SV set in a new sample
sniffles --input new.bam --genotype-vcf known_svs.vcf --vcf genotyped.vcf

# Single-sample low-VAF / mosaic (NOT a tumor-normal caller)
sniffles --input tumor.bam --vcf mosaic.vcf --mosaic
```

The `.snf` is a binary signature index (NOT a VCF - never bcftools it); it retains sub-threshold signatures so the merge re-genotypes an SV even in a sample that did not independently pass support.

## cuteSV - the per-platform parameter trap

cuteSV's defaults are not platform-appropriate; the README gives distinct sets by error rate. `--genotype` is OFF by default. Positional args: `cuteSV <bam> <ref> <out.vcf> <work_dir>`. Force-calling moved to the separate cuteFC tool.

| Platform | --max_cluster_bias_INS | --diff_ratio_merging_INS | --max_cluster_bias_DEL | --diff_ratio_merging_DEL |
|----------|------------------------|--------------------------|------------------------|--------------------------|
| ONT | 100 | 0.3 | 100 | 0.3 |
| PacBio HiFi/CCS | 1000 | 0.9 | 1000 | 0.5 |
| PacBio CLR | 100 | 0.3 | 200 | 0.5 |

```bash
mkdir cutesv_work
cuteSV aln.bam ref.fa cutesv.vcf cutesv_work --genotype \
  --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
  --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3   # ONT set
```

## Benchmarking with Truvari

```bash
truvari bench --base giab_tier1.vcf.gz --comp calls.vcf.gz \
  --includebed tier1_regions.bed --pctseq 0.7 --refdist 500 --passonly -o bench/
truvari refine bench/        # re-harmonize TR-region representations for a fair comparison
```

`--pctseq` (default 0.7) compares the actual inserted/deleted sequence, not just coordinates - set 0 for depth-based callers lacking alt sequence, keep 0.7 for long-read callers. Region set dominates the headline: Tier1 (resolvable INS/DEL >=50 bp) overstates whole-genome performance; CMRG reflects hard clinical loci. Tier1 v0.6 is INS/DEL only - do not report INV recall against it.

## Per-Method Failure Modes

### One VNTR fragments into many false positives
**Trigger:** calling in tandem repeats without a TR BED. **Mechanism:** the breakpoint slides across repeat units, scattering signatures. **Symptom:** several calls with inconsistent breakpoints where one event exists. **Fix:** supply `--tandem-repeats` to the caller; `truvari refine` when benchmarking.

### cuteSV defaults inflate or fragment calls
**Trigger:** running cuteSV with one parameter set across platforms. **Mechanism:** HiFi settings over-merge ONT noise; ONT settings fragment clean HiFi signatures. **Symptom:** FP inflation or split calls. **Fix:** use the platform-matched set; remember `--genotype` is off by default.

### Missing insertion sequence / breakpoints
**Trigger:** Sniffles without `--reference`, or alignment without `-Y`. **Mechanism:** no reference -> no ALT sequence; hard-clipped supplementaries -> lost breakpoint sequence. **Symptom:** INS lack sequence; imprecise breakpoints. **Fix:** add `--reference` and align with `-Y`.

### Treating Sniffles --mosaic as a cancer caller
**Trigger:** somatic SV calling with single-sample `--mosaic`. **Mechanism:** mosaic mode lowers support in one sample; it has no normal to subtract. **Symptom:** germline SVs reported as somatic; FP at low VAF. **Fix:** Severus or nanomonsv (paired tumor-normal).

### Comparing F1 across studies that handled repeats differently
**Trigger:** quoting F1 without region + TR BED + Truvari params. **Mechanism:** representation handling moves the number more than the caller. **Symptom:** apples-to-oranges comparisons. **Fix:** fix the region set, TR BED, and Truvari params; run `truvari refine`.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| SV >= 50 bp | GIAB convention | 30-100 bp is a VNTR gray zone where callers disagree |
| Sniffles `--minsvlen` 35, `--mapq` 25, `--minsupport auto` | Sniffles2 manpage | the actual defaults (support is coverage-derived, not a fixed 3) |
| Coverage ~20-30x germline; >30-60x mosaic/somatic | SV practice | large SVs callable from 5-10x; low-VAF needs depth |
| Truvari `--pctseq 0.7`, `--refdist 500` | English 2022 | sequence-aware INS matching; loosen refdist to 1000 only for fuzzy callers |
| cuteSV params per platform | cuteSV README | error rate sets cluster bias / merge ratio |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Many FP calls in repeats | no TR BED | supply `--tandem-repeats` |
| cuteSV VCF has no GT | `--genotype` off by default | add `--genotype` |
| Cannot bcftools the `.snf` | `.snf` is a binary signature index | use it as Sniffles input, not a VCF |
| INS records lack sequence | `--reference` not supplied | add `--reference ref.fa` |
| Imprecise/missing breakpoints | supplementaries hard-clipped | align with minimap2 `-Y` |
| Looking for cuteSV force-calling flag | moved to cuteFC | use the cuteFC tool |
| Somatic SVs from a single sample | germline/mosaic caller | Severus / nanomonsv (paired) |

## References

- Smolka M, Paulin LF, Grochowski CM, et al. 2024. Detection of mosaic and population-level structural variants with Sniffles2. *Nat Biotechnol* 42:1571-1580.
- Jiang T, Liu Y, Jiang Y, et al. 2020. Long-read-based human genomic structural variation detection with cuteSV. *Genome Biol* 21:189.
- Heller D, Vingron M. 2019. SVIM: structural variant identification using mapped long reads. *Bioinformatics* 35:2907-2915.
- English AC, Menon VK, Gibbs RA, Metcalf GA, Sedlazeck FJ. 2022. Truvari: refined structural variant comparison preserves allelic diversity. *Genome Biol* 23:271.
- Zook JM, Hansen NF, Olson ND, et al. 2020. A robust benchmark for detection of germline large deletions and insertions. *Nat Biotechnol* 38:1347-1355.
- Wagner J, Olson ND, Harris L, et al. 2022. Curated variation benchmarks for challenging medically relevant autosomal genes (CMRG). *Nat Biotechnol* 40:672-680.
- Keskus AG, Bryant A, Ahmad T, et al. 2026. Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing. *Nat Biotechnol* 44:247-257.

## Related Skills

- long-read-alignment - SV-ready mapping (`-Y` soft-clip, platform preset)
- basecalling - Read accuracy/length that gates breakpoint precision
- clair3-variants - Small variants (<50 bp) are Clair3's job, not an SV caller's
- haplotype-phasing - Haplotag the BAM for haplotype-specific / phased SVs
- genome-assembly/hifi-assembly - Phased assembly for assembly-based SV calling
- variant-calling/structural-variant-calling - The variant-calling-side SV view
- variant-calling/vcf-manipulation - Filter/merge the SV VCFs
- genome-intervals/gtf-gff-handling - Annotate SVs against gene models
