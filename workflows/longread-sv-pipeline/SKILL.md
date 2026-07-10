---
name: bio-workflows-longread-sv-pipeline
description: Orchestrates an end-to-end long-read structural-variant pipeline - basecalling to minimap2 alignment (platform-matched preset) to Sniffles2/cuteSV/pbsv calling to optional assembly-based calling (dipcall/PAV) to two-step .snf cohort merging to Truvari benchmarking - chaining ONT and PacBio HiFi runs while handing the SV signal mechanism off to the component skills. Use when running a long-read SV workflow from reads to a benchmarked callset, choosing the minimap2 preset and SV caller by platform and goal, deciding when long reads are worth it for the insertions and repeat-mediated SVs short reads physically miss, building a joint-genotyped cohort with the two-step .snf design, or parameterizing a Truvari benchmark against GIAB HG002 Tier 1 plus CMRG. Not for the SV signal mechanism itself (see variant-calling/structural-variant-calling) or short-read SV.
tool_type: cli
primary_tool: Sniffles
workflow: true
depends_on:
  - long-read-sequencing/basecalling
  - long-read-sequencing/long-read-alignment
  - long-read-sequencing/long-read-qc
  - long-read-sequencing/structural-variants
qc_checkpoints:
  - after_qc: "Read N50 >10kb, mean quality >Q10 (ONT R10) / HiFi rq>0.99"
  - after_alignment: "Mapping rate >90%, coverage >=15x for confident SV calling"
  - after_calling: "SV count in the expected range, INS/DEL ratio sane, genotypes concordant"
  - after_benchmark: "Truvari F1 reported WITH its refdist/pctsize/pctseq, on Tier 1 AND CMRG"
---

## Version Compatibility

Reference examples tested with: minimap2 2.28+, Sniffles 2.2+, cuteSV 2.1+, pbsv 2.9+, dipcall 0.3+, bcftools 1.19+, samtools 1.19+, truvari 4.0+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Use minimap2 >= 2.28: the `lr:hq` accurate-read preset was added in 2.27, and 2.28 fixes the 2.27 `--MD` regression. Supply a reference-matched tandem-repeat BED to the caller - it is the single biggest false-positive lever in repeats. Truvari renamed the alt-sequence-similarity param from `--pctsim` to `--pctseq` at v4; confirm against `truvari bench --help`.

If code throws an error, introspect the installed tool and adapt the example to the actual API rather than retrying.

# Long-Read SV Pipeline

**"Detect structural variants from my long-read sequencing data"** -> Chain basecalling, platform-matched minimap2 alignment, an SV caller selected by platform and goal, an optional assembly-based branch, cohort merging, and a parameterized Truvari benchmark - with the SV mechanism delegated to the component skills.

This is a workflow (orchestration) skill: it makes the stage-to-stage decisions and quality gates that connect the component skills. It does NOT re-teach how a caller sees an SV or the VCF representation minefield - that lives in variant-calling/structural-variant-calling and long-read-sequencing/structural-variants.

## Why long reads for SV: the physics that justifies this pipeline

Run this pipeline instead of a short-read SV workflow for one mechanistic reason: a single long read (PacBio HiFi ~15-25 kb, ONT tens of kb to >Mb ultralong) physically spans the SV and both flanks in one molecule, turning SV detection from an inference-over-fragments problem into near-direct observation. Two consequences decide whether long reads earn their per-sample cost:

- Insertions become tractable. Placing and sizing an insertion needs reads that carry the novel bases; when an INS exceeds Illumina read length (150 bp) no short read spans it, so short-read INS recall is stuck at ~30-50% while long reads reach ~90%+. Ebert 2021 (*Science* 372:eabf7117) found 68% of 107,590 assembly-discovered SVs were missed by short reads. If insertions matter, this pipeline is the answer, not a tuning knob.
- Repeat-mediated junctions resolve. A 20 kb read anchored in unique sequence on both flanks spans a breakpoint buried in a 5 kb repeat that no 150 bp read can straddle, bringing segmental-duplication NAHR, mobile-element insertions, and VNTR/STR expansions into reach.

The governing pipeline principle: an SV call is a *representation artifact* of choices made upstream. The aligner preset, the tandem-repeat BED handed to the caller, and the Truvari matching parameters move precision/recall as much as the caller does. Chaining decisions - not the caller name - are what this skill is about.

## Pipeline map

```
POD5/FAST5 (ONT only)
    |  [Step 0] Dorado basecall  -> long-read-sequencing/basecalling
    v            (model + methylation are IRREVERSIBLE choices; sup model for SV)
FASTQ (ONT / PacBio HiFi)
    |  [QC]     NanoPlot / NanoComp -> long-read-sequencing/long-read-qc
    v            gate: read N50 >10 kb, sane quality, chimera screen
[Step 1] minimap2 alignment       -> long-read-sequencing/long-read-alignment
    |            preset by platform/chemistry; -Y keeps breakpoint seq on split reads
    v            gate: mapping rate >90%, coverage >=15x
[Step 2] SV calling               -> long-read-sequencing/structural-variants
    |            Sniffles2 / cuteSV / pbsv, caller by platform + goal
    |            + tandem-repeat BED (biggest FP lever)      variant-calling/structural-variant-calling
    v
[Step 3, optional] assembly-based SV (dipcall / PAV against a phased diploid assembly)
    |            highest-quality callset; the way truth sets are built
    v
[Step 4] cohort merge             -> two-step Sniffles2 .snf (per-sample -> combine)
    |
    v
[Step 5] benchmark                -> Truvari vs GIAB HG002 Tier 1 + CMRG
                 an F1 is meaningless without refdist/pctsize/pctseq
```

## Platform decision: ONT vs PacBio HiFi

The platform sets the preset, the caller options, and what bonus channels are available. Decide before basecalling.

| Dimension | ONT (R10.4.1) | PacBio HiFi |
|-----------|---------------|-------------|
| Per-base accuracy | ~Q20+ simplex, higher duplex | ~Q30+ (circular consensus) |
| Read length | tens of kb; ultralong >Mb achievable | ~15-25 kb |
| minimap2 preset | `lr:hq` (R10/Q20 accurate) or `map-ont` (older R9) | `map-hifi` |
| Best for | ultralong spans, repeat/centromere traversal, native methylation | highest base accuracy, small variants + SV in one run |
| SV caller | Sniffles2 or cuteSV | Sniffles2, cuteSV, or pbsv (official, TR-aware) |
| Bonus channel | 5mCG/6mA methylation if requested AT basecall time | 5mCG via kinetics; phasing native from HiFi length |

R10.4 chemistry moved ONT simplex to ~Q20, which is why `lr:hq` (not the noisy-read `map-ont`) is the right preset for modern ONT - it rewrites the scoring/chaining model for accurate reads and runs faster at equal accuracy. Older R9 data still needs `map-ont`.

## SV caller selection (by platform and goal)

Do not re-derive the caller mechanism here; pick by goal and hand tuning to the component skill.

| Goal / platform | Caller | Why / cross-reference |
|-----------------|--------|-----------------------|
| ONT/HiFi germline, cohorts, mosaic | Sniffles2 | field standard; two-step `.snf` population merge scales linearly in N; `--mosaic` for low-VAF (Smolka 2024) |
| Highest recall on noisy ONT | cuteSV | signature clustering; MUST pass the per-platform param set and `--genotype` (Jiang 2020) |
| PacBio HiFi, official, TR-aware | pbsv | expects pbmm2 alignments; single- and joint-sample modes |
| tandem-vs-interspersed DUP detail | SVIM | reports origin AND destination of duplications (Heller 2019) |
| Highest-quality callset / truth set | dipcall or PAV | assembly-vs-reference from a phased diploid assembly (Li 2018; Ebert 2021) |
| Somatic (tumor-normal) | Severus / nanomonsv | matched-normal subtraction; do NOT use Sniffles `--mosaic` for somatic (Keskus 2025) |

Methods evolve; verify current best practice against each tool's docs before committing. Deeper caller tuning (cuteSV per-platform params, the tandem-repeat BED, aligner effects) lives in long-read-sequencing/structural-variants.

## Step 0: Basecalling (ONT only)

**Goal:** Convert raw POD5/FAST5 signal into reads suitable for SV calling, capturing methylation if it will ever be needed.

**Approach:** Basecall with Dorado using a chemistry-matched `sup` (super-accuracy) model; request modified bases at basecall time because methylation cannot be recovered later. PacBio HiFi arrives as reads already, so this step is skipped.

```bash
# sup model maximizes accuracy for SV; 5mCG_5hmCG requested now (irreversible if omitted).
# See long-read-sequencing/basecalling for model selection and duplex.
dorado basecaller sup pod5_dir/ --modified-bases 5mCG_5hmCG > reads.bam
samtools fastq -T MM,ML reads.bam | gzip > reads.fastq.gz   # -T carries methylation tags through
```

## Step 1: Alignment

**Goal:** Produce a sorted, indexed BAM whose split (supplementary) alignments retain the breakpoint sequence SV callers reconstruct from.

**Approach:** Align with the platform-matched minimap2 preset; keep `-Y` so supplementary alignments are soft-clipped (not hard-clipped), which preserves the junction bases on split reads. SV calling rides on supplementary, not secondary, alignments.

```bash
# ONT R10/Q20: lr:hq (accurate reads). Older R9: map-ont. HiFi: map-hifi. PacBio CLR: map-pb.
minimap2 -ax lr:hq -t 16 --MD -Y reference.fa reads.fastq.gz | \
    samtools sort -@ 4 -o aligned.bam
samtools index aligned.bam
```

**QC checkpoint** (gate before spending compute on calling):

```bash
samtools flagstat aligned.bam                                  # mapping rate should be >90%
samtools depth -a aligned.bam | awk '{s+=$3} END{print "mean cov:", s/NR}'
# Gate: >=15x for confident SV calling; below ~10x callers drift toward false negatives.
```

## Step 2: SV calling

**Goal:** Call SVs (>=50 bp DEL/INS/DUP/INV/BND) from the aligned reads with a caller matched to the platform.

**Approach:** Run Sniffles2 (the default) with a reference-matched tandem-repeat BED - it clusters the repeat-driven false positives that otherwise dominate the callset. `--minsvlen 50` enforces the GIAB >=50 bp SV convention (Sniffles2 defaults to 35).

```bash
# Sniffles2: the tandem-repeat BED is the single biggest false-positive lever in repeats.
sniffles --input aligned.bam --reference reference.fa \
    --tandem-repeats human_GRCh38_TR.bed \
    --vcf svs.vcf.gz --threads 8 --minsvlen 50 --output-rnames
```

cuteSV as an alternative - its defaults are NOT platform-appropriate, and `--genotype` is off by default:

```bash
# ONT param set shown. HiFi: 1000/0.9/1000/0.5. CLR: 100/0.3/200/0.5. See structural-variants.
mkdir -p work_dir   # cuteSV requires the work dir to pre-exist; it does not create it
cuteSV aligned.bam reference.fa svs.vcf work_dir/ --threads 8 --genotype \
    --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3
```

## Step 3 (optional): Assembly-based SV

**Goal:** Produce the highest-quality SV callset by comparing a phased diploid assembly to the reference, rather than inferring from read alignments.

**Approach:** Assemble the genome (hifiasm/verkko), then call variants from the two haplotype assemblies aligned to the reference. dipcall (the syndip method) and PAV (the HGSVC method) are the assembly-vs-reference callers; this is how the GIAB and HGSVC truth sets themselves are built. Use it when an assembly already exists or when callset quality outranks turnaround.

```bash
# dipcall needs two haplotype assemblies (hap1/hap2) plus minimap2/k8/htsbox on PATH.
run-dipcall prefix reference.fa hap1.fa hap2.fa > prefix.mak
make -j2 -f prefix.mak                                          # emits prefix.dip.vcf.gz + prefix.dip.bed
```

## Step 4: Cohort merging (the two-step .snf design)

**Goal:** Build a joint-genotyped multi-sample SV matrix, not a union of per-sample discovery VCFs.

**Approach:** Sniffles2's population design processes each sample independently into a compact `.snf`, then combines the `.snf` files in a second pass - scaling linearly in N. A union of per-sample discovery VCFs is wrong: a sample recorded 0/0 may simply not have had that event *discovered* in it (a false missing), which corrupts allele frequencies. The two-step `.snf` combine force-genotypes every sample at every merged site.

```bash
# Pass 1: per-sample .snf (each sample processed once, independently).
for s in sample1 sample2 sample3; do
    sniffles --input ${s}.bam --reference reference.fa \
        --tandem-repeats human_GRCh38_TR.bed --snf ${s}.snf
done
# Pass 2: combine into a jointly genotyped cohort VCF (linear in N).
sniffles --input sample1.snf sample2.snf sample3.snf --vcf cohort.vcf.gz
```

For sequence-aware AF work across callsets, prefer Truvari `collapse` over position-only merging (position-only mergers inflate allele frequency by up to 2.2x; English 2022) - see variant-calling/structural-variant-calling for the merger decision table.

## Step 5: Benchmarking - an SV F1 is meaningless without its parameters

**Goal:** Report a defensible, reproducible accuracy figure - not a number that looks good because of loose matching.

**Approach:** Truvari `bench` counts a call as a true positive only if it matches a truth variant under ALL of `--refdist`, `--pctsize`, and `--pctseq` simultaneously. Every one of these moves the score, so a bare F1 is uninterpretable. Report the full parameter set, run `truvari refine` for a harmonized re-comparison, and stratify by region.

```bash
# Report EVERY parameter. --pctseq 0 disables alt-sequence checking and quietly inflates INS scores.
truvari bench -b HG002_SV_Tier1.vcf.gz -c svs.vcf.gz -o bench_tier1/ --passonly \
    -f reference.fa --refdist 500 --pctsize 0.70 --pctseq 0.70 --sizemin 50   # -f persists reference to params.json
truvari refine bench_tier1/                                     # harmonized breakpoint re-comparison (needs the reference bench recorded)

# CMRG is NOT optional: Tier 1 EXCLUDES the medically relevant repetitive genes.
truvari bench -b HG002_CMRG_SV.vcf.gz -c svs.vcf.gz -o bench_cmrg/ --passonly \
    --refdist 500 --pctsize 0.70 --pctseq 0.70 --sizemin 50
```

Three escalating bars are routinely conflated, and a pipeline can pass the first while failing the ones that matter:

- Event detection - something of about the right type/size near the right place (loose refdist, no sequence check). Easy.
- Breakpoint accuracy - POS/END within a few bp (tight `--refdist`, `--pctseq` on). Matters at exon/splice boundaries.
- Genotype accuracy - the sample GT (het/hom) is correct (genotype-aware comparison). A caller can detect an event perfectly and still call het-as-hom, which is fatal for Mendelian analyses.

Region stratification is decisive: Tier 1 (Zook 2020 *Nat Biotechnol* 38:1347) is conservative isolated SVs, while CMRG (Wagner 2022 *Nat Biotechnol* 40:672) covers the repetitive medically relevant genes Tier 1 leaves out - where GRCh38 false duplications cause reference-specific misses that masking raised from 8% to 100% recall. A good Tier 1 F1 certifies nothing about the genes clinicians care about; run both.

## Filtering and annotation

```bash
bcftools view -i 'QUAL>=20 && ABS(SVLEN)>=50' svs.vcf.gz -Oz -o svs.filtered.vcf.gz
bcftools index svs.filtered.vcf.gz          # ABS() is mandatory: DEL SVLEN is negative by convention
bcftools stats svs.filtered.vcf.gz > sv_stats.txt

AnnotSV -SVinputFile svs.filtered.vcf.gz -genomeBuild GRCh38 -outputFile annotated_svs
# gene overlap, DGV/gnomAD-SV population AF, ClinVar pathogenicity
```

## Phased and methylation-aware SV (bonus channels)

Heterozygous variants on the same long read are physically phased, so SVs can be assigned to haplotypes with no statistical phasing. Haplotag the BAM (whatshap/`sniffles --phase`) before or during calling to get haplotype-resolved SVs; see long-read-sequencing/haplotype-phasing. If methylation was requested at basecall time (Step 0), the MM/ML tags ride through alignment (via minimap2 `-y` / `samtools fastq -T`) and give a per-haplotype methylation channel alongside the SV call at no extra sequencing cost - useful for imprinting and allele-specific silencing, but it must be captured at basecall time or it is gone.

## SV types detected

| Type | ALT | Notes for long reads |
|------|-----|----------------------|
| Deletion | DEL | excellent recall; breakpoints base-precise when a read spans the junction |
| Insertion | INS | the reason to use long reads; the read carries the inserted sequence |
| Duplication | DUP | tandem vs interspersed distinguishable (SVIM reports origin + destination) |
| Inversion | INV | resolved when unique anchors flank the repeat-embedded breakpoints |
| Translocation | BND | paired breakend records linked by MATEID; complex events are BND graphs |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Few SVs / missing known INS | coverage <10x or missing tandem-repeat BED | raise depth to >=15x; pass `--tandem-repeats` |
| Many false positives in repeats | no tandem-repeat BED supplied | provide a reference-matched TR BED (biggest FP lever) |
| `map-ont` on R10 data is slow/less accurate | wrong preset for accurate reads | use `lr:hq` for R10/Q20 ONT; `map-ont` only for R9 |
| Split reads lost breakpoint sequence | aligned without `-Y` (hard-clipped supplementaries) | re-align with `-Y` |
| Methylation channel gone | not requested at basecall time | rebasecall with `--modified-bases`; it is irreversible |
| cuteSV recall poor / no genotypes | ran defaults; `--genotype` off | pass the per-platform param set and `--genotype` |
| Cohort "0/0" wrong, AF too low | took a union of per-sample discovery VCFs | use the two-step `.snf` combine (force-genotypes all sites) |
| `ABS(SVLEN)>=50` filter drops all deletions | filtered raw SVLEN (DEL is negative) | always wrap in `ABS()` |
| Truvari F1 not reproducible / suspiciously high | reported without params, or `--pctseq 0` | state refdist/pctsize/pctseq/sizemin; never disable pctseq to look good |
| Passed Tier 1 but clinical genes fail | benchmarked only on Tier 1 | also run CMRG (Tier 1 excludes those genes) |

## Related Skills

- long-read-sequencing/basecalling - Dorado model choice and requesting methylation at basecall time (Step 0)
- long-read-sequencing/long-read-alignment - minimap2 preset selection, `-Y` soft-clipping, MM/ML tag passthrough
- long-read-sequencing/long-read-qc - read-length/quality QC and chimera screening before alignment
- long-read-sequencing/structural-variants - caller tuning (cuteSV per-platform params, tandem-repeat BED, Truvari) - the SV mechanism for long reads
- long-read-sequencing/haplotype-phasing - haplotag the BAM for phased/somatic SVs
- variant-calling/structural-variant-calling - the SV signal model, SVLEN-sign / symbolic-vs-BND / CIPOS representation, force-genotyping, sequence-aware merging (also short-read SV)
- variant-calling/consensus-sequences - why symbolic `<DEL>`/`<INS>` alleles are not directly consensus-able

## References

- Li H. Minimap2: pairwise alignment for nucleotide sequences. 2018 *Bioinformatics* 34:3094-3100.
- Sedlazeck FJ, Rescheneder P, Smolka M, Fang H, Nattestad M, von Haeseler A, Schatz MC. Accurate detection of complex structural variations using single-molecule sequencing. 2018 *Nature Methods* 15:461-468. (Sniffles v1 + NGMLR)
- Smolka M, Paulin LF, Grochowski CM, Horner DW, Mahmoud M, Behera S, et al. Detection of mosaic and population-level structural variants with Sniffles2. 2024 *Nature Biotechnology*. doi:10.1038/s41587-023-02024-y. (two-step .snf population merge; mosaic SVs)
- Jiang T, Liu Y, Jiang Y, Li J, Gao Y, Cui Z, et al. Long-read-based human genomic structural variation detection with cuteSV. 2020 *Genome Biology* 21:189.
- Heller D, Vingron M. SVIM: structural variant identification using mapped long reads. 2019 *Bioinformatics* 35:2907-2915.
- Li H, Bloom JM, Farjoun Y, Fleharty M, Gauthier L, Neale B, MacArthur D. A synthetic-diploid benchmark for accurate variant-calling evaluation. 2018 *Nature Methods* 15:595-597. (dipcall/syndip)
- Ebert P, Audano PA, Zhu Q, Rodriguez-Martin B, Porubsky D, Bonder MJ, et al. Haplotype-resolved diverse human genomes and integrated analysis of structural variation. 2021 *Science* 372:eabf7117. (PAV; 68% of SVs missed by short reads)
- Keskus AG, et al. Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing. 2025 *Nature Biotechnology*. doi:10.1038/s41587-025-02618-8.
- English AC, Menon VK, Gibbs RA, Metcalf GA, Sedlazeck FJ. Truvari: refined structural variant comparison preserves allelic diversity. 2022 *Genome Biology* 23:271. (defaults refdist 500, pctsize 0.70, pctseq 0.70, sizemin 50; up to 2.2x AF inflation from position-only merging)
- Zook JM, Hansen NF, Olson ND, Chapman L, Mullikin JC, Xiao C, et al. A robust benchmark for detection of germline large deletions and insertions. 2020 *Nature Biotechnology* 38:1347-1355. (GIAB HG002 SV Tier 1)
- Wagner J, Olson ND, Harris L, McDaniel J, Cheng H, Fungtammasan A, et al. Curated variation benchmarks for challenging medically relevant autosomal genes. 2022 *Nature Biotechnology* 40:672-680. (GIAB-CMRG; false-duplication masking raises recall 8%->100%)
- pbsv - PacBio structural variant caller (no dedicated publication): github.com/PacificBiosciences/pbsv
