---
name: bio-long-read-sequencing-clair3-variants
description: Calls germline small variants (SNPs and indels) from Oxford Nanopore and PacBio HiFi long reads with Clair3, a two-stage (pileup + full-alignment) deep-learning caller, selecting the chemistry- and basecaller-version-matched model, enabling read-based phasing, and benchmarking against GIAB with stratification. Covers why the model string is the experiment (no auto-detection, silent degradation on mismatch), why ONT homopolymer/STR indels are the residual error whole-genome F1 hides, and the somatic/trio/RNA boundary to the ClairS/Clair3-Trio family. Use when calling germline SNVs/indels from ONT or HiFi BAMs, choosing a Clair3 model, phasing variants, or benchmarking long-read calls.
tool_type: cli
primary_tool: Clair3
---

## Version Compatibility

Reference examples tested with: Clair3 2.0+, whatshap 2.0+, bcftools 1.19+, hap.py 0.3.15+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Results depend on inputs that outlive the binary version - record them:
- The Clair3 MODEL must match the platform + chemistry + basecaller tier + basecaller version (e.g. `r1041_e82_400bps_sup_v500`). There is NO auto-detection; `--model_path` is mandatory and a mismatch silently degrades calls.
- Clair3 v2 moved TensorFlow -> PyTorch; models are `pileup.pt`/`full_alignment.pt`. v1 TensorFlow models do NOT load in v2.
- The full ONT model set (every version, hac/fast, `_with_mv` signal-aware) lives in the rerio `clair3_models/` repo; only a subset is bundled.

If code throws an error, introspect the installed tool (`run_clair3.sh --help`) and adapt the example to the actual API rather than retrying.

# Clair3 Variant Calling

**"Call variants from my long reads"** -> Run Clair3 with the model that matches how the reads were basecalled, phase, and benchmark with stratification - because the model string, not the command, determines accuracy.
- CLI: `run_clair3.sh --bam_fn=aln.bam --ref_fn=ref.fa --output=out/ --threads=16 --platform=ont --model_path=/models/r1041_e82_400bps_sup_v500`

Scope: germline diploid SNPs + small indels. NOT structural variants (-> structural-variants), NOT somatic/mosaic (-> ClairS/ClairS-TO), NOT RNA (-> Clair3-RNA).

## The Single Most Important Modern Insight -- The Model String Is the Experiment, and ONT Indels Hide in the Strata

Clair3's accuracy is gated by two facts a naive user misses:

1. **The model is hand-picked and a mismatch fails silently.** There is no auto-detection - the user must point `--model_path` at a specific model folder. Three axes must ALL match: chemistry (`r941` vs `r1041`), basecaller tier (`fast`/`hac`/`sup`), and basecaller version (`g5014`/`v430`/`v500`/`v520`), plus the optional `_with_mv` signal-aware axis if the BAM has Dorado `mv` tags. Wrong model = no crash, no warning, measurably worse calls (indels most). Derive the model from the basecaller string in the run metadata; pick the model version closest to but not above the basecaller version.
2. **ONT indels in homopolymers/STRs are the residual error whole-genome F1 conceals.** Even on R10.4.1 sup, insertions/deletions in homopolymer runs and short tandem repeats are the weak point (G/C homopolymers worst), because the pore cannot reliably count identical consecutive bases. A genome-wide indel F1 of ~99.5% hides much lower performance inside LowComplexity/homopolymer strata - exactly the medically relevant loci. HiFi largely solves this; do not transfer ONT-indel pessimism to HiFi. Always benchmark with GIAB stratification, never a single global number.

## Two-Stage Architecture

Clair3 "symphonizes" two networks: a fast **pileup model** (summarized per-position statistics) that calls the large majority of sites, and a slow **full-alignment model** (haplotype-resolved read tensor) that re-evaluates only the uncertain subset. Internally Clair3 phases the top het-SNP pileup calls with WhatsHap, haplotags the BAM, and feeds the haplotagged reads to the full-alignment model - which is why read-based phasing buys ~6% indel F1, not cosmetics. Output: `merge_output.vcf.gz` (final).

## Model Selection

Model name anatomy (`r1041_e82_400bps_sup_v500`): pore (`r1041`=R10.4.1), flowcell (`e82`), speed (`400bps`), basecaller tier (`sup`/`hac`/`fast`), basecaller version (`v500`=Dorado 5.0.0, `g5014`=Guppy 5.0.14). The `_with_mv` suffix uses Dorado move-table tags for best accuracy when present.

| Data | `--platform` | Model |
|------|--------------|-------|
| ONT R10.4.1 sup, Dorado v5.x, mv tags present | ont | `r1041_e82_400bps_sup_v520_with_mv` |
| ONT R10.4.1 sup, Dorado v5.0.0 | ont | `r1041_e82_400bps_sup_v500` |
| ONT R10.4.1 hac | ont | `r1041_e82_400bps_hac_v500`/`_v520` |
| ONT R9.4.1 (any tier) | ont | `r941_prom_sup_g5014` |
| PacBio HiFi Revio | hifi | `hifi_revio` |
| PacBio HiFi Sequel II | hifi | `hifi_sequel2` |
| Illumina (supported) | ilmn | `ilmn` |
| PacBio CLR | - | not supported -> PEPPER-Margin-DeepVariant |

## Decision Tree by Scenario

| Scenario | Tool | Why |
|----------|------|-----|
| Germline SNV/indel, single sample | Clair3 | this skill |
| Somatic, paired tumor-normal | ClairS | VAF-aware; Clair3 germline priors cannot find low-VAF somatic |
| Somatic, tumor-only | ClairS-TO | tumor-only ensemble |
| De novo / Mendelian trio | Clair3-Nova / Clair3-Trio | family-aware |
| Long-read RNA variants | Clair3-RNA | RNA model |
| ONT R10.4.1, also considering DeepVariant | either | neck-and-neck on R10 sup; native-ONT DeepVariant (Kolesnikov 2024) superseded PEPPER-Margin |
| Non-human / draft / bacterial reference | Clair3 + `--include_all_ctgs` | default calls only chr1-22,X,Y -> empty output otherwise |
| Cohort joint genotyping | Clair3 gVCF -> GLnexus | `bcftools merge` on gVCFs is NOT joint genotyping |

## Core Commands

```bash
# Germline ONT calling (model MUST match the basecaller)
run_clair3.sh \
  --bam_fn=aln.bam --ref_fn=ref.fa --output=clair3_out/ \
  --threads=16 --platform=ont \
  --model_path=/opt/models/r1041_e82_400bps_sup_v500
# final VCF: clair3_out/merge_output.vcf.gz

# Phase the final output VCF (WhatsHap); --longphase_for_phasing swaps only the INTERNAL
# phaser to LongPhase (faster, SV-aware). For a LongPhase-phased final VCF use
# --use_longphase_for_final_output_phasing instead of --enable_phasing.
run_clair3.sh ... --enable_phasing --longphase_for_phasing
# Phased calls go to clair3_out/phased_merge_output.vcf.gz; merge_output.vcf.gz stays UNPHASED.

# Non-human / draft assembly reference - call ALL contigs
run_clair3.sh ... --include_all_ctgs

# Targeted / amplicon panel
run_clair3.sh ... --bed_fn=panel.bed --gvcf

# Benchmark against GIAB with stratification (the step that reveals ONT indel errors)
hap.py giab_truth.vcf.gz clair3_out/merge_output.vcf.gz \
  -f giab_confident.bed -r ref.fa --engine=vcfeval \
  --stratification giab_stratifications.tsv -o bench/hg002
```

## Per-Method Failure Modes

### Silent model mismatch
**Trigger:** `--model_path` pointing at a model that does not match the basecaller chemistry/tier/version. **Mechanism:** no auto-detection; the wrong network runs. **Symptom:** no error, lower F1 (indels most). **Fix:** derive the model from the basecaller string; verify the folder exists (rerio for the full set); for v2 ensure `.pt` models.

### Clair3 found nothing on a non-human reference
**Trigger:** bacterial genome or draft assembly without chr1-22,X,Y names. **Mechanism:** Clair3 calls only standard human contigs by default. **Symptom:** near-empty VCF. **Fix:** `--include_all_ctgs`.

### Global F1 looks great, clinical genes are wrong
**Trigger:** reporting only whole-genome F1. **Mechanism:** ONT indel errors concentrate in homopolymer/STR/low-complexity strata. **Symptom:** ~99.5% global indel F1 but much lower in LowComplexity. **Fix:** stratify with GIAB BEDs (Dwarshuis 2024); use CMRG for medically relevant genes.

### Treating Clair3 as a somatic caller
**Trigger:** lowering `--snp_min_af`/`--indel_min_af` to catch low-VAF variants. **Mechanism:** germline model expects ~0.5/1.0 allele fractions, is not VAF-aware. **Symptom:** germline-model false positives at low AF, missed true somatic. **Fix:** ClairS (paired) / ClairS-TO (tumor-only).

### v1 model with v2 Clair3
**Trigger:** an old TensorFlow model dir with Clair3 v2. **Mechanism:** v2 needs PyTorch `.pt` models. **Symptom:** model load failure. **Fix:** use `pileup.pt`/`full_alignment.pt` models (Converted Rerio).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Recommended depth ~20-60x | Clair3 guidance | sensitivity (hets, indels) falls off below ~20x; `--min_coverage` default 2 is a floor, not a recommendation |
| Phasing buys ~6% indel F1 | Zheng 2022 | haplotagged reads disambiguate indel alleles in repeats |
| ONT R10.4.1 sup: SNP F1 ~99.99%, indel F1 ~99.5% | GIAB benchmarks | indel residual lives in homopolymer/STR strata |
| `--var_pct_full` 0.3 (default) | Clair3 README | fraction of low-quality pileup calls re-run by full-alignment; raise for recall, slower |
| Stratify with GIAB / CMRG | Dwarshuis 2024 | global F1 hides the ONT indel problem |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Empty/near-empty VCF on non-human ref | default calls only chr1-22,X,Y | `--include_all_ctgs` |
| Model fails to load | v1 TF model with v2 Clair3 | use `.pt` (PyTorch) models |
| `--model_path .../models/ont` not found | no generic `ont`/`hifi` model | point at a specific model subfolder |
| Worse-than-expected indels | wrong-version or wrong-tier model | match the basecaller model exactly |
| "joint genotyping" gave odd merges | `bcftools merge` on gVCFs is not joint calling | use GLnexus |
| Looking for somatic/low-VAF variants | germline caller | use ClairS / ClairS-TO |

## References

- Zheng Z, Li S, Su J, Leung AWS, Lam TW, Luo R. 2022. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling (Clair3). *Nat Comput Sci* 2:797-803.
- Zheng Z, He M, Yu X, et al. 2026. Accelerated long-read variant calling with Clair3 for whole-genome sequencing. *Bioinformatics* (advance access) btag181.
- Kolesnikov A, Cook D, Nattestad M, et al. 2024. Local read haplotagging enables accurate long-read small variant calling. *Nat Commun* 15:5907.
- Dwarshuis N, Kalra D, McDaniel J, et al. 2024. The GIAB genomic stratifications resource for human reference genomes. *Nat Commun* 15:9029.
- Lin JH, Chen LC, Yu SC, Huang YT. 2022. LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants. *Bioinformatics* 38(7):1816-1822.
- Chen L, Zheng Z, Su J, et al. 2025. ClairS-TO: a deep-learning method for long-read tumor-only somatic small variant calling. *Nat Commun* 16:9630.

## Related Skills

- basecalling - The basecaller model+version the Clair3 model must match
- long-read-alignment - Produces the BAM (keep `--MD`; use minimap2 >=2.28)
- haplotype-phasing - whatshap/longphase phasing and haplotagging Clair3 uses internally
- medaka-polishing - ONT consensus; medaka diploid variant calling is deprecated in favor of Clair3
- structural-variants - SVs are out of Clair3's scope (Sniffles2/cuteSV)
- variant-calling/deepvariant - DeepVariant native ONT/HiFi models (neck-and-neck on R10)
- variant-calling/vcf-statistics - Summarize/filter the VCF Clair3 emits
- clinical-databases/variant-prioritization - Prioritize the called variants
