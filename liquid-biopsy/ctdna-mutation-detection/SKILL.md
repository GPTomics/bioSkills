---
name: bio-ctdna-mutation-detection
description: Detects somatic mutations in circulating tumor DNA, treating low-VAF detection as a signal-versus-noise problem set by error suppression and molecules sampled, not by the choice of caller. Distinguishes de novo CALLING (scanning a panel for unknown variants, bounded by per-locus error and multiple testing) from tumor-informed DETECTION (tracking a pre-specified variant set, where panel integration reaches single-ppm). Covers VarDict and Mutect2 for de novo calling, UMI-aware callers, and a pysam-based known-variant VAF tracker, with matched-WBC subtraction as the mandatory defense against clonal hematopoiesis (the dominant false positive). Use when calling or tracking tumor mutations from plasma cfDNA, setting a VAF threshold, or deciding whether a low-VAF call is tumor versus CHIP.
tool_type: mixed
primary_tool: VarDict
---

## Version Compatibility

Reference examples tested with: pysam 0.22+, pandas 2.2+, VarDictJava 1.8+, GATK 4.5+, Ensembl VEP 111+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: VarDict's `-c -S -E -g` are 1-based BED COLUMN INDICES, not genomic coordinates; var2vcf_valid.pl's `-E` suppresses the END tag (opposite meaning to VarDict's `-E`). VEP gnomAD flags are `--af_gnomade` (exomes)/`--af_gnomadg` (genomes); the bare `--af_gnomad` is a legacy alias that returns only exome AF, so prefer the explicit forms.

# ctDNA Mutation Detection

**"Detect mutations in my cfDNA sample"** -> Either scan a panel for unknown low-VAF somatic variants (de novo calling) or quantify a pre-specified mutation set across samples (tumor-informed tracking) — two different statistical problems.
- CLI: `vardict-java | teststrandbias.R | var2vcf_valid.pl` for de novo low-VAF calling on a consensus BAM
- CLI: `gatk Mutect2` with the read-orientation model for de novo calling with artifact filtering
- Python: `pysam` pileup of ref/alt counts at fixed loci for known-variant tracking (MRD)

## The Single Most Important Modern Insight -- detection is a signal-vs-noise problem, and tracking is not the same problem as calling

Low-VAF ctDNA detection is set by two limits that no caller can overcome: the per-base error floor (raw Illumina ~1e-3 caps naive VAF detection near 0.5-1%) and the number of tumor molecules physically present in the tube (1 ng cfDNA ~= 303 haploid genome-equivalents; at 0.01% VAF in 10 ng the expected mutant count is ~0.3 copies — there is nothing to detect at any depth). The achieved limit of detection is the worse of the two. Error suppression (UMI consensus -> ~1e-5, duplex -> <1e-7) and input mass move the floor; swapping VarDict for Mutect2 does not.

Critically, de novo CALLING and known-variant DETECTION are different statistical problems. De novo calling scans every covered position for an unknown alt and pays a multiple-testing tax across 1e5-1e6 loci, so per-locus thresholds must be stringent (practical LoD ~0.1-0.5% on UMI consensus). Tumor-informed detection tests ONE hypothesis — "is tumor present?" — by integrating signal across a pre-specified set of N patient-specific loci; the multiple-testing penalty collapses and per-locus signal that is individually indistinguishable from noise sums into a confident panel-level call. This is why per-locus LoD is poor while panel-integrated LoD reaches single-ppm. Conflating the two is the most common conceptual error in the field.

## Methods Landscape

| Method | Class | Citation | Role | When |
|--------|-------|----------|------|------|
| VarDict / vardict-java | de novo caller | AstraZeneca-NGS | sensitive low-VAF amplicon/capture calling with explicit strand-bias test | de novo panel calling on a UMI-consensus BAM |
| Mutect2 (tumor-only) | de novo caller | GATK | local-assembly somatic caller + learned orientation-bias artifact model | de novo calling needing FFPE/OxoG artifact filtering, PoN, germline resource |
| umi-varcal | UMI-aware caller | Sater 2020 *Bioinformatics* 36(9):2718 | own UMI-aware pileup + per-position Poisson test against local background | UMI-tagged BAM where a consensus-aware caller is wanted (floor ~0.3%) |
| CAPP-Seq / iDES | tumor-informed integration | Newman 2014 *Nat Med* 20:548; 2016 *Nat Biotechnol* 34:547 | hybrid-capture deep panel + molecular barcoding + in-silico background polishing | de novo ctDNA to ~0.02%; iDES stacks ~15x error suppression |
| INVAR | tumor-informed integration | Wan 2020 *Sci Transl Med* 12:eaaz8084 | integrate variant reads across 100s-1000s patient loci, background-weighted | MRD/monitoring with tumor WES; quantifies to ~1e-5, best ~2.5 ppm |
| MRDetect | tumor-informed integration | Zviran 2020 *Nat Med* 26:1114 | shallow WGS vs patient SNV compendium, read-level SVM noise model | MRD trading depth for breadth (~35x WGS, thousands of SNVs); ~1e-5 |
| PhasED-Seq | tumor-informed integration | Kurtz 2021 *Nat Biotechnol* 39:1537 | enrich phased (co-occurring) variants to suppress single-molecule error | sub-ppm MRD where phased variants are available |

Per-locus LoD for any of these is error- and sampling-limited (~0.1-0.5%); the tumor-informed methods reach ppm only by integrating across a known, large, patient-specific variant set. Methodology evolves — verify current best practice against each tool's live docs before committing to one.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Tumor tissue available, MRD/monitoring of a known cancer | Tumor-informed tracking (INVAR/MRDetect/Signatera-class), or the pysam tracker below for a fixed list | Integrating across N pre-specified loci is the only route to ppm; CHIP excluded by construction (CHIP variants are not on the tumor list) |
| No tumor tissue, screening/discovery | de novo panel calling (VarDict or Mutect2) + matched WBC | Must scan for unknown variants; CHIP subtraction is mandatory or most calls are not tumor |
| VAF regime > 1% | Any standard caller on a deduplicated BAM | Above the raw error floor; consensus not strictly required |
| VAF regime 0.1-1% | UMI single-strand consensus + VarDict/umi-varcal/Mutect2 | Below the raw 1e-3 floor; consensus needed to recover real signal from error |
| VAF regime < 0.1% | Duplex consensus + tumor-informed integration | Single-strand consensus cannot remove one-strand deamination/oxidation; only duplex + panel integration reaches this regime |
| No matched WBC available | Do NOT report de novo calls as somatic-tumor | Without WBC subtraction, CHIP (the majority of non-germline cfDNA variants) is indistinguishable from tumor |

## CHIP -- the dominant false positive, not background

Clonal hematopoiesis of indeterminate potential (CHIP) is the single largest source of false-positive somatic calls in plasma, and it is the null hypothesis for any low-VAF cfDNA variant. Razavi 2019 sequenced cfDNA with matched white-blood-cell DNA (508 genes, >60,000x) and found that **53.2% of non-germline cfDNA variants in cancer patients and 81.6% in non-cancer controls** had features consistent with clonal hematopoiesis; only ~24.4% of cfDNA somatic variants in patients were also in the matched tumor (the remainder split between white-cell CHIP and variants of uncertain origin). These are bona fide somatic mutations — in cancer genes — that come from lysed leukocytes, not tumor. No error-suppression tier removes them because they are not errors.

The biology: CHIP arises in hematopoietic stem cells and rises steeply with age (Jaiswal 2014: ~10% prevalence over age 70), enriched for PPM1D/TP53/CHEK2 clones after prior chemo/radiation — exactly the monitored population. The canonical genes are DNMT3A, TET2, ASXL1 (the big three), then PPM1D, TP53, JAK2, SF3B1, SRSF2, GNB1, GNAS, CBL, ATM, CHEK2. TP53 and ATM are both CHIP genes and bona fide tumor suppressors, so a low-VAF TP53 cfDNA call is the ambiguous case par excellence.

The only reliable filter is matched buffy-coat/WBC subtraction: sequence the WBC fraction of the same draw at comparable depth and remove any cfDNA variant also present in WBC. gnomAD filtering removes germline only — CHIP variants are somatic and absent from germline databases, so they sail straight through. A canonical-CHIP-gene list (the example's `CHIP_GENES`) is a heuristic flag for extra scrutiny, NOT a substitute for WBC subtraction. See analytical-validation for the LoB/LoD statistics that quantify how confidently a subtracted call clears background.

## De Novo Calling with VarDict

**Goal:** Scan a target panel for unknown low-VAF somatic variants on a UMI-consensus BAM.

**Approach:** Run vardict-java with a lowered `-f`, pipe through the strand-bias test, then convert to VCF — matching `-f` across both stages so the threshold is not silently re-applied.

```bash
AF_THR=0.005   # 0.5% — practical UMI-consensus de novo floor; below this approaches the per-base error floor
vardict-java -G ref.fa -f $AF_THR -N sample -b consensus.bam \
  -c 1 -S 2 -E 3 -g 4 targets.bed | \
  teststrandbias.R | \
  var2vcf_valid.pl -N sample -E -f $AF_THR > sample.vcf
```

Key flags: `-G` indexed reference; `-f` min VAF (VarDict default 0.01); `-N` sample name; `-b` BAM. `-c 1 -S 2 -E 3 -g 4` are the 1-based BED COLUMN INDICES for chrom/start/end/gene in a standard 4-column BED — they are column positions, not genomic values. On var2vcf_valid.pl, `-E` means "do NOT print the END tag" (unrelated to VarDict's `-E`); its `-f` default is 0.02, so set it to match. For PCR/amplicon data add `-P 0` (positional std is expected to be ~0). For paired tumor/normal use the `testsomatic.R | var2vcf_paired.pl` path instead.

## De Novo Calling with Mutect2 and the Orientation-Bias Model

**Goal:** Call de novo somatic variants while filtering FFPE-deamination (C>T) and OxoG (G>T) strand-biased artifacts that dominate low-VAF false positives.

**Approach:** Collect F1R2/F2R1 counts during calling, learn the orientation-bias prior, then apply it during filtering alongside a panel of normals and germline resource.

```bash
gatk Mutect2 -R ref.fa -I consensus.bam --f1r2-tar-gz f1r2.tar.gz \
  --germline-resource af-only-gnomad.vcf.gz --panel-of-normals pon.vcf.gz \
  -O unfiltered.vcf.gz
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
gatk FilterMutectCalls -R ref.fa -V unfiltered.vcf.gz \
  --ob-priors read-orientation-model.tar.gz -O filtered.vcf.gz
```

Mutect2 is run tumor-only here (no normal sample arg); the orientation model is the load-bearing low-VAF filter. At true ctDNA VAFs Mutect2 is underpowered relative to a dedicated UMI/duplex + background-polishing pipeline, and local assembly can miss extremely low-AF alt support — it is a reasonable de novo caller on consensus reads with the orientation model + PoN (and ideally a matched normal, not shown in this tumor-only command), not a substitute for tumor-informed integration at ppm.

## Track Known Mutations Across Serial Samples

**Goal:** Quantify the VAF of a pre-specified mutation set at fixed loci for MRD monitoring — the detection (not calling) problem.

**Approach:** For each target mutation, pileup reads at the position, count ref/alt/other alleles, and compute VAF with depth; aggregate across loci as the panel-level detection signal. The single-base pileup below tracks SNVs only — indel reporters (e.g. EGFR exon-19 deletions) need `read.indel`/CIGAR-aware counting; a single-base comparison silently scores every indel read as `other` and reports the locus as cleared.

```python
import pysam

def track_known_variants(bam_file, variants):
    '''Pileup ref/alt counts at fixed (chrom, pos, ref, alt) SNV loci; pos is 1-based.
    SNVs only - indel reporters need read.indel/CIGAR handling, not a single-base compare.'''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    rows = []
    for chrom, pos, ref, alt in variants:
        counts = {'ref': 0, 'alt': 0, 'other': 0}
        for col in bam.pileup(chrom, pos - 1, pos, truncate=True):
            for read in col.pileups:
                if read.is_del or read.is_refskip:
                    continue
                base = read.alignment.query_sequence[read.query_position]
                counts['alt' if base == alt else 'ref' if base == ref else 'other'] += 1
        depth = sum(counts.values())
        rows.append({'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                     'depth': depth, 'alt_count': counts['alt'],
                     'vaf': counts['alt'] / depth if depth else 0.0})
    bam.close()
    return rows
```

Annotate calls for interpretation with Ensembl VEP (`--cache --offline --fasta --vcf --everything`); the gnomAD allele-frequency flags are `--af_gnomade` (exomes) and `--af_gnomadg` (genomes) — the bare `--af_gnomad` is a legacy alias returning only exome AF. gnomAD presence separates germline; only WBC presence separates CHIP.

## Per-Method Failure Modes

### CHIP misclassified as tumor
Trigger: de novo calling without matched WBC. Mechanism: leukocyte-derived clonal somatic variants in cancer genes look identical to tumor signal and pass gnomAD filtering. Symptom: low-VAF calls in DNMT3A/TET2/TP53; "tumor" mutations not in the matched tissue. Fix: subtract matched buffy-coat/WBC genotype; never report somatic-tumor without it.

### Strand-biased / deamination artifacts at low VAF
Trigger: FFPE-style C>T or oxidative G>T at VAF near the floor. Mechanism: damage on one template strand is inherited by every PCR copy, so single-strand UMI consensus votes unanimously for the artifact. Symptom: alt support concentrated on one strand. Fix: VarDict strand-bias test or Mutect2 orientation model; for sub-0.1% require duplex consensus.

### Calling below the error floor
Trigger: lowering `-f` to e.g. 0.001 on a non-consensus BAM. Mechanism: the raw ~1e-3 error rate manufactures alt reads at that frequency. Symptom: a flood of low-VAF calls scaling with depth. Fix: do consensus upstream; do not set a VAF threshold below the demonstrated error floor of the input.

### Germline-vs-somatic confusion at low coverage
Trigger: classifying by VAF alone when depth is low. Mechanism: a true 50% het reads 3/12 = 0.25 by chance. Symptom: germline hets mislabeled subclonal somatic. Fix: gnomAD + matched-WBC presence (germline ~0.5 in WBC; CHIP at clone VAF; tumor-only absent from WBC).

### Per-locus LoD quoted as the assay LoD
Trigger: reporting a single-variant sensitivity for a multi-locus tracking assay (or vice versa). Mechanism: panel-integrated LoD is orders of magnitude below per-locus LoD. Symptom: a "0.1%" claim that does not match observed ppm-level tracking. Fix: state per-locus vs panel-integrated explicitly; see analytical-validation.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Raw Illumina error ~1e-3 caps naive VAF near 0.5-1% | Schmitt 2012 *PNAS* 109:14508 | Per-base miscall rate sets the per-locus VAF floor; alt support below it is mostly error |
| UMI single-strand consensus -> ~1e-5; duplex -> <1e-7 | Schmitt 2012; Newman 2016 *Nat Biotechnol* 34:547 | Family consensus erases PCR/sequencing error; duplex strand concordance also catches one-strand damage |
| CAPP-Seq de novo LoD ~0.02% at 96% specificity | Newman 2014 *Nat Med* 20:548 | Deep hybrid-capture + reporter set; demonstrates the de novo panel floor |
| iDES ~15x error suppression (UMI ~3x x polishing ~3x) | Newman 2016 *Nat Biotechnol* 34:547 | Molecular consensus and in-silico background polishing are orthogonal and stack |
| INVAR quantifies to ~1e-5, detects to ~2.5 ppm | Wan 2020 *Sci Transl Med* 12:eaaz8084 | Integrating variant reads across 100s-1000s of patient loci collapses multiple testing |
| MRDetect ~1e-5 tumor fraction at ~35x WGS, 95% spec | Zviran 2020 *Nat Med* 26:1114 | Breadth (thousands of SNVs) + read-level SVM (~14.4x error reduction) substitutes for depth |
| CHIP = 53.2% (cancer pts) / 81.6% (controls) of cfDNA variants | Razavi 2019 *Nat Med* 25:1928 | Most non-germline cfDNA variants are not tumor; matched WBC is mandatory |
| Depth >= 1000-5000x unique consensus for panels | community / Phallen 2017 *Sci Transl Med* 9:eaan2415 (~30,000x) | Detecting <1% VAF needs enough unique molecules sampled at each locus |
| ~303 genome-equivalents per ng cfDNA (3.3 pg/haploid) | standard constant | Input mass sets a hard Poisson ceiling on detectable VAF independent of sequencing |
| LoB / LoD / LoD95 per CLSI EP17 | CLSI EP17-A2 | A bare VAF without input mass + replicate detection rate is not a sensitivity spec |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Flood of low-VAF calls scaling with depth | `-f` set below the input's error floor on non-consensus reads | Do UMI/duplex consensus first; keep `-f` >= demonstrated floor |
| VarDict emits nothing or wrong regions | `-c -S -E -g` read as genomic values | They are 1-based BED column indices; use `-c 1 -S 2 -E 3 -g 4` for a 4-column BED |
| var2vcf re-filters away VarDict calls | var2vcf_valid.pl `-f` default 0.02 mismatched | Set var2vcf `-f` to match VarDict's `-f` |
| Real amplicon calls dropped as positional artifacts | var2vcf `-P` (filter pstd=0) on by default | Add `-P 0` for PCR/amplicon data |
| "Tumor" variants absent from matched tissue | CHIP not subtracted | Sequence and subtract matched WBC; flag CHIP-gene hits |
| `--af_gnomad` returns only exome AF | bare flag is a legacy exome-only alias | Use `--af_gnomade` (exomes) / `--af_gnomadg` (genomes) |
| ppm "LoD" not reproducible | per-locus LoD quoted for a tracking assay | Report panel-integrated LoD with input mass and LoD95 |

## References

- Schmitt MW, Kennedy SR, Salk JJ, Fox EJ, Hiatt JB, Loeb LA. 2012. Detection of ultra-rare mutations by next-generation sequencing. *Proc Natl Acad Sci USA* 109(36):14508-14513. — Duplex Sequencing; DCS error <1e-7.
- Newman AM, Bratman SV, To J, et al. 2014. An ultrasensitive method for quantitating circulating tumor DNA with broad patient coverage. *Nat Med* 20(5):548-554. — CAPP-Seq; ~0.02% LoD.
- Newman AM, Lovejoy AF, Klass DM, et al. 2016. Integrated digital error suppression for improved detection of circulating tumor DNA. *Nat Biotechnol* 34(5):547-555. — iDES; UMI + background polishing.
- Phallen J, Sausen M, Adleff V, et al. 2017. Direct detection of early-stage cancers using circulating tumor DNA. *Sci Transl Med* 9(403):eaan2415. — TEC-Seq deep panel.
- Razavi P, Li BT, Brown DN, et al. 2019. High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. *Nat Med* 25(12):1928-1937. — CHIP is the majority of cfDNA variants; matched WBC.
- Wan JCM, Heider K, Gale D, et al. 2020. ctDNA monitoring using patient-specific sequencing and integration of variant reads. *Sci Transl Med* 12(548):eaaz8084. — INVAR; integration to ~2.5 ppm.
- Zviran A, Schulman RC, Shah M, et al. 2020. Genome-wide cell-free DNA mutational integration enables ultra-sensitive cancer monitoring. *Nat Med* 26(7):1114-1124. — MRDetect; shallow WGS + read SVM.
- Kurtz DM, Soo J, Co Ting Keh L, et al. 2021. Enhanced detection of minimal residual disease by targeted sequencing of phased variants in circulating tumor DNA. *Nat Biotechnol* 39(12):1537-1547. — PhasED-Seq; phased-variant enrichment.
- Sater V, Viailly P-J, Lecroq T, et al. 2020. UMI-VarCal: a new UMI-based variant caller that efficiently improves low-frequency variant detection in paired-end sequencing NGS libraries. *Bioinformatics* 36(9):2718-2724. — UMI-aware pileup + per-position Poisson test.

## Related Skills

- cfdna-preprocessing - UMI/duplex consensus input that sets the error floor
- analytical-validation - LoD/LoB and the panel-integration math behind detection
- longitudinal-monitoring - track detected variants across serial samples
- tumor-fraction-estimation - orthogonal burden estimate to cross-check
- variant-calling/variant-calling - general somatic calling principles
- clinical-databases/variant-prioritization - clinical annotation and interpretation
