---
name: bio-workflows-cnv-pipeline
description: Orchestrates the copy-number pipeline from BAM to segmented, integer-called, annotated CNVs, forking on germline-vs-somatic - CNVkit (somatic exome/panel: coverage -> assay-matched reference/PoN -> fix -> segment -> purity/ploidy-aware call), GATK gCNV (germline rare-CNV cohort), and allele-specific callers (ASCAT/FACETS/PURPLE) for purity/ploidy. Use when committing the build + target/access BED + PoN once (assay-matched), building the reference from normals BEFORE segmenting, fitting purity/ploidy BEFORE integer calls in tumors, centering on the true (non-diploid) mode before GISTIC2 recurrence, or routing cfDNA to ichorCNA. Hands mechanism to the copy-number component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: CNVkit
goal_approach_exempt: true
workflow: true
depends_on:
  - copy-number/cnvkit-analysis
  - copy-number/gatk-cnv
  - copy-number/copy-ratio-segmentation
  - copy-number/allele-specific-copy-number
  - copy-number/cnv-visualization
  - copy-number/cnv-annotation
  - copy-number/recurrent-cnv
qc_checkpoints:
  - after_coverage: "Uniform coverage across targets; flag systematically low-depth targets (capture dropout)"
  - after_fix: "log2-ratio noise (.cnr spread/MAD) within tolerance; high bin noise -> over-segmentation"
  - after_call: "Integer CN off a fitted purity/ploidy (not defaults); tumor purity above the ~40% death zone"
  - after_recurrent: "GISTIC2 input is diploid-CENTERED (uncentered WGD inverts recurrence)"
---

## Version Compatibility

Reference examples tested with: CNVkit 0.9.10+, GATK 4.5+ (gCNV / ModelSegments), ASCAT/FACETS/PURPLE (allele-specific), GISTIC2 2.0.23 (recurrent), ichorCNA 0.5+ (cfDNA)

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: GATK gCNV runs COHORT mode (model all normals, no prior) vs CASE mode (score a singlet against a prior model) — order is DetermineGermlineContigPloidy -> GermlineCNVCaller -> PostprocessGermlineCNVCalls. Sequenza's `copynumber` dependency was REMOVED from Bioconductor 3.18+ (needs a fork). GATK gCNV/ModelSegments have no single method paper — cite the GATK docs. Confirm in-tool before quoting.

# CNV Pipeline

**"Detect copy number variants from my sequencing data"** -> Fork germline-vs-somatic, commit the build + target/access BED + assay-matched reference, bias-correct against normals, segment, and integer-call off a fitted purity/ploidy.
- CLI: cnvkit target/access/antitarget -> coverage -> reference(normals) -> fix -> segment -> call  (OR GATK gCNV for germline cohorts)

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A CNV callset is decided at three seams, not inside the caller.

1. **The reference build + target/access BED + reference/PoN is one made-once commitment inherited by everything downstream.** The capture-kit target BED, the `access` mappability BED, and the annotation refFlat must all be the SAME build as the BAMs (a GRCh37 BED against GRCh38 BAMs silently produces zero-coverage bins). And the PoN is the identity of the assay: it MUST be built from the same capture kit, chemistry, and (ideally) batch as the cases. A PoN from a different kit imports the wrong bias profile and fabricates CNVs at capture boundaries.
2. **The reference/PoN is built BEFORE anything is segmented, and it absorbs shared signal.** `fix` needs the reference to bias-correct; segmenting raw log2 without normalizing segments the capture bias, not biology. Beware: tangent normalization / a pooled PoN ABSORBS any CNV shared across the normals — a real common CNV becomes invisible; GC correction alone does NOT remove the replication-timing wave.
3. **The diploid baseline is a commitment, not a given — fit purity/ploidy BEFORE integer calls in tumors.** In WGD/hyper-aneuploid tumors the data mode is not diploid; naive centering inverts every call. `cnvkit.py call` with wrong `--purity`/`--ploidy` (or defaults on an impure/WGD tumor) assigns integer copy numbers off the wrong baseline. Fit purity/ploidy (ASCAT/FACETS/PURPLE) first; below ~40% purity calls degrade and below ~20% no bulk caller works.

## Pipeline map

```
BAM (tumor +/- matched normal, OR germline cohort)
  | fork: germline rare-CNV cohort? --> GATK gCNV  (copy-number/gatk-cnv)
  v  else somatic exome/panel:
  | [1] target/access/antitarget BED (build-matched)   (copy-number/cnvkit-analysis)
  v
  | [2] per-sample coverage
  v
  | [3] build reference/PoN from NORMALS first (assay-matched)
  v     ^-- tangent absorbs CNV shared across the PoN
  | [4] fix (bias-correct) -> segment -> call
  v     ^-- purity/ploidy fitted BEFORE integer call (copy-number/allele-specific-copy-number)
  | [5] visualize + gene-level annotate            (copy-number/cnv-visualization, cnv-annotation)
  v
  | [6] (cohort) center on true mode -> GISTIC2 recurrence  (copy-number/recurrent-cnv)
  v
Segmented, integer-called, annotated CNVs
```

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Build + target/access/refFlat BED | Any build mismatch -> zero-coverage bins / shifted annotations |
| Reference / PoN (assay-matched) | A different-kit PoN imports the wrong bias -> false CNVs at capture boundaries; tangent absorbs CNVs shared across the PoN |
| Purity/ploidy (fitted, not default) | Wrong baseline shifts every integer call; WGD inverts calls |
| Diploid centering (cohort) | Uncentered WGD segments into GISTIC2 invert recurrence |

## The canonical order and why

1. **Prepare target/access/antitarget BEDs** on the committed build.
2. **Per-sample coverage** (target + antitarget/off-target bins).
3. **Build the reference/PoN from normals FIRST** — order-trap: `fix` needs the reference; segmenting raw log2 segments capture bias.
4. **fix -> segment -> call**, with purity/ploidy fitted BEFORE the integer call — order-trap: default purity/ploidy on an impure/WGD tumor mis-assigns every integer CN.
5. **Visualize + gene-level annotate** (positive control: known CNVs recovered if present).
6. **(Cohort) center on the true mode, THEN GISTIC2** — order-trap: uncentered WGD inverts recurrence; and do NOT concatenate per-sample `.cns` and call recurrence naively — feed a diploid-centered `.seg` matrix to GISTIC2.

## Choosing the caller (the germline-vs-somatic fork)

Pipeline-level selection only; mechanism lives in the component skills.

| Situation | Lean toward | Hand off to |
|-----------|-------------|-------------|
| Exome/targeted panel, somatic (tumor) CNV | CNVkit (target + antitarget bins) | copy-number/cnvkit-analysis |
| Germline rare-CNV from a cohort of exomes | GATK gCNV (DetermineGermlineContigPloidy -> GermlineCNVCaller -> PostprocessGermlineCNVCalls) | copy-number/gatk-cnv |
| WGS, need allele-specific CN + purity/ploidy | ASCAT / Sequenza / FACETS / PURPLE | copy-number/allele-specific-copy-number |
| Relative copy-ratio segments (research) | GATK ModelSegments/CallCopyRatioSegments | copy-number/copy-ratio-segmentation |
| Cohort recurrent/driver CNV | GISTIC2 (diploid-centered input) | copy-number/recurrent-cnv |
| cfDNA / low-pass tumor fraction | ichorCNA (NOT CNVkit) | workflows/liquid-biopsy-pipeline |

## Primary path: CNVkit (somatic exome/panel)

```bash
# 1. Targets on the committed build (annotate with refFlat, split for WES)
cnvkit.py target capture_targets.bed --annotate refFlat.txt --split -o targets.bed
cnvkit.py access genome.fa -o access.bed
cnvkit.py antitarget targets.bed --access access.bed -o antitargets.bed

# 2-3. Coverage per sample, then build the reference from NORMALS (assay-matched) BEFORE any fix
cnvkit.py coverage $bam targets.bed -o cov/${s}.targetcoverage.cnn
cnvkit.py coverage $bam antitargets.bed -o cov/${s}.antitargetcoverage.cnn
cnvkit.py reference cov/normal*.{,anti}targetcoverage.cnn --fasta genome.fa -o reference.cnn

# 4. fix (bias-correct) -> segment -> call. Fit purity/ploidy first (ASCAT/FACETS) for tumors:
cnvkit.py fix cov/${s}.targetcoverage.cnn cov/${s}.antitargetcoverage.cnn reference.cnn -o ${s}.cnr
cnvkit.py segment ${s}.cnr -o ${s}.cns
cnvkit.py call ${s}.cns --purity 0.6 --ploidy 2 -o ${s}.call.cns   # purity/ploidy from an allele-specific fit
```

A runnable somatic CNVkit script (manual target -> coverage -> reference -> fix -> segment -> call path) is in this skill's examples/; germline cohorts use GATK gCNV (copy-number/gatk-cnv), not CNVkit.

## QC checkpoints between steps

| After | Gate | Interpretation |
|-------|------|----------------|
| Coverage | Uniform depth across targets; flag low-depth targets | Capture dropout -> phantom deletions |
| fix | `.cnr` log2 spread / MAD within tolerance | High bin noise is the #1 CNV false-positive lever (over-segmentation) |
| segment/call | Sane segment count; integer CN consistent with known events; purity plausible | Over-segmentation = noisy reference / low purity; wrong purity shifts every call |
| annotate | Known CNVs recovered (positive control) | Build/BED mismatch surfaces as missing known events |
| recurrent | GISTIC2 input diploid-centered | Uncentered WGD inverts recurrence |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Zero-coverage bins / shifted annotations | Target BED build != BAM build | Pin one build across BED, access, refFlat, BAMs |
| False CNVs at capture boundaries | PoN from a different kit/chemistry | Build the PoN from the same kit/chemistry/batch |
| A real common CNV vanishes | Tangent/pooled PoN absorbed the shared signal | Use a PoN that does not carry the event, or germline-CNV logic |
| Every integer call shifted / inverted | Default purity/ploidy on an impure/WGD tumor | Fit purity/ploidy (ASCAT/FACETS/PURPLE) BEFORE `call` |
| Inverted recurrence in the cohort | Uncentered WGD segments into GISTIC2 | Center on the true (non-diploid) mode first |
| Cohort recurrence looks wrong | Concatenated per-sample `.cns` naively | Feed a diploid-centered `.seg` matrix to GISTIC2 (copy-number/recurrent-cnv) |
| Sequenza install fails | `copynumber` removed from Bioconductor 3.18+ | Use a maintained fork (ShixiangWang/igordot) |

## Related Skills

- copy-number/cnvkit-analysis - CNVkit coverage/fix/segment/call details
- copy-number/gatk-cnv - GATK gCNV (germline cohort) and ModelSegments
- copy-number/copy-ratio-segmentation - segmentation algorithm and depth-bias correction
- copy-number/allele-specific-copy-number - purity/ploidy and integer allele-specific CN (ASCAT/FACETS/PURPLE)
- copy-number/cnv-visualization - scatter/diagram/heatmap plotting
- copy-number/cnv-annotation - gene-level CNV annotation
- copy-number/recurrent-cnv - cohort recurrent/driver CNV with GISTIC2
- copy-number/hrd-scoring - HRD scar score for PARP eligibility
- workflows/liquid-biopsy-pipeline - cfDNA tumor-fraction CNV (ichorCNA)
- workflows/somatic-variant-pipeline - consumes purity/ploidy for VAF-to-CCF

## References

- Steele CD, Abbasi A, Islam SMA, et al (2022) Signatures of copy number alterations in human cancer. *Nature* 606:984-991. DOI 10.1038/s41586-022-04738-6. (copy-number signatures need ABSOLUTE CN.)
- Telli ML, Timms KM, Reid J, et al (2016) Homologous Recombination Deficiency (HRD) score predicts response to platinum-containing neoadjuvant chemotherapy. *Clinical Cancer Research* 22:3764-3773. DOI 10.1158/1078-0432.CCR-15-2477. (GIS >= 42 HRD threshold.)
- GATK gCNV / ModelSegments have no single method paper — cite the GATK/Broad documentation.
