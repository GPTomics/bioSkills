---
name: bio-phasing-imputation-haplotype-phasing
description: Estimates haplotype phase from population linkage disequilibrium with SHAPEIT5, SHAPEIT4, Eagle2, or Beagle - turning unphased genotypes (0/1) into phased haplotypes (0|1) for imputation input, compound-heterozygote calls, HLA typing, or population genetics. Covers why statistical phase is an INFERENCE (not a measurement) whose error concentrates at rare variants, why a genome-wide switch-error rate hides catastrophic rare-variant error and must be reported MAC-stratified, the SHAPEIT5 common-scaffold-then-rare design (phase_common, ligate, phase_rare, switch), reference-based vs within-cohort phasing, the build-matched genetic map, chrX male-haploid handling, and the switch-vs-flip-vs-Hamming distinction. Use when phasing genotypes before imputation, for compound-het/ASE/HLA, or benchmarking against trios. Read-backed / molecular phasing (long reads, Hi-C) is long-read-sequencing/haplotype-phasing; panel choice is reference-panels; imputation is genotype-imputation.
tool_type: cli
primary_tool: SHAPEIT5
---

## Version Compatibility

Reference examples tested with: SHAPEIT5 5.1.1, Eagle 2.4.1, Beagle 5.4 (22Jul22), bcftools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

SHAPEIT4 to SHAPEIT5 changed the CLI substantially: SHAPEIT5 is a SUITE of binaries (`phase_common`, `phase_rare`, `ligate`, `switch`), not a single `shapeit` command, and `phase_common` is the engine formerly known as SHAPEIT4. The genetic map and the reference panel must match the data's genome build (GRCh37 vs GRCh38); a build-mismatched map silently degrades phasing. PBWT and Ne defaults have drifted between betas; confirm against the installed `--help`.

# Statistical Haplotype Phasing -- Inferring Phase From Population LD

**"Resolve which alleles sit together on each chromosome"** -> Estimate haplotype phase from population linkage disequilibrium via the Li-Stephens HMM - because phase is INFERRED statistically from how haplotypes are shared across a population, not read off the genotype, so a switch error is a model uncertainty (the rate, not zero, is the deliverable), not a typo.
- CLI: `phase_common --input target.bcf --filter-maf 0.001 --map chr20.b38.gmap.gz --region chr20 --output scaffold.bcf` then `ligate` then `phase_rare` (SHAPEIT5), or Eagle2/Beagle for common-variant phasing

Scope: population/statistical phasing of array or sequence genotypes for imputation input, compound-het/ASE/HLA, and population genetics. Read-backed / molecular single-sample phasing (long reads, Hi-C, 10x linked reads) is a PHYSICALLY DIFFERENT signal -> long-read-sequencing/haplotype-phasing (the two are easily conflated; do not run SHAPEIT on long-read evidence or trust statistical phase for a private clinical variant). Panel choice -> reference-panels. Imputation against a panel -> genotype-imputation. The input VCF and biallelic normalization -> variant-calling/variant-normalization. End-to-end orchestration -> workflows/gwas-pipeline.

## The Single Most Important Modern Insight -- A Phased Haplotype Is a Statistical Estimate, and Its Error Concentrates Exactly Where the Biology of Interest Lives

Statistical phasing reconstructs which alleles are on the same chromosome by borrowing LD across many individuals or a reference panel (Delaneau 2019 *Nat Commun* 10:5436). That works beautifully for common variants in LD with their neighbors and fails, by construction, for rare variants - which are young, carried by few people, and in LD with almost nothing. Three facts drive every decision:

1. **The genome-wide switch-error rate lies, because it is dominated by easy common sites.** A headline "switch error rate 0.3%" is averaged over millions of common heterozygous sites and says nothing about the singleton or doubleton that is most likely to be the compound-het, the de-novo, or the pathogenic allele of interest - those are phased at MAC-dependent accuracy an order of magnitude worse, and a true singleton is essentially a coin flip without special machinery (Hofmeister 2023 *Nat Genet* 55:1243). Report accuracy stratified by minor allele count, never as one number.
2. **The deliverable is a switch-error rate against an independent truth set, not the tool name.** "We used SHAPEIT" is not a switch-error rate. A switch error changes which haplotype an allele sits on without changing any genotype, so it is invisible to every per-site genotype QC; for any phase-dependent claim, measure the rate against a trio (Mendelian truth via `switch --pedigree`) or read-backed truth.
3. **The modern arc is the scaffold design, and rare-variant phasing needs biobank scale to work at all.** SHAPEIT5 phases common variants into a fixed, near-perfect scaffold, then places each rare allele onto it by PBWT/IBD haplotype matching - which depends on finding a long shared haplotype, itself a function of cohort size. This is why rare-variant phasing in a small cohort cannot be trusted for a cis/trans call without orthogonal (trio or read-backed) evidence.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| SHAPEIT5 | Hofmeister 2023 *Nat Genet* 55:1243 | suite (phase_common/phase_rare/ligate/switch); scaffold design for rare/singleton phasing; PBWT | biobank-scale WGS/WES; rare-variant phasing |
| SHAPEIT4 (= phase_common engine) | Delaneau 2019 *Nat Commun* 10:5436 | sub-linear common-variant phasing; integrates panels, scaffolds, read-backed phase | common-variant phasing / pre-phasing; legacy |
| Eagle2 | Loh 2016 *Nat Genet* 48:1443 | HMM + PBWT-derived HapHedge; reference-based (`--vcfRef`) and within-cohort | array data; the classic imputation-server phaser |
| Beagle 5.x | Browning 2021 *Am J Hum Genet* 108:1880 | Java; does BOTH phasing (gt=, no ref=) and imputation; two-stage for sequence | one tool for phase and impute; no compile |
| Trio / pedigree phasing | (Mendelian transmission) | deterministic phase where the trio is informative | gold standard; validating other phasers via `switch` |
| WhatsHap (boundary) | Patterson 2015 *J Comput Biol* 22:498 | read-backed phasing (weighted MEC) from aligned reads | -> long-read-sequencing/haplotype-phasing; can seed SHAPEIT as a scaffold |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Array data, small-to-modest cohort, have a panel | Eagle2 `--vcfRef` or phase_common `--reference` | a panel models LD better than a few thousand samples |
| Array data, large cohort, no panel | Eagle2 or phase_common within-cohort | LD is modeled from the cohort; accuracy rises with N |
| WGS/WES, biobank scale, need rare variants phased | SHAPEIT5: phase_common -> ligate -> phase_rare | the scaffold design is the only route to accurate rare-variant phase |
| Pre-phasing as imputation input | Eagle2 or Beagle 5 | small switch errors largely wash out in imputation -> genotype-imputation |
| One tool for phase and impute, no compile | Beagle 5.x (gt= to phase, add ref= to impute) | pragmatic single tool |
| Trio/pedigree available | trio/pedigree phasing; use `switch` to benchmark | deterministic where informative; the truth ruler |
| Long reads on the same sample | -> long-read-sequencing/haplotype-phasing (then seed SHAPEIT as a scaffold) | read-backed phase is local and deterministic; combine, do not replace |
| Common-variant phasing only, modest data | SHAPEIT4 or Beagle | rare-variant machinery is unnecessary overhead |

## The Common-Scaffold-Then-Rare Design (SHAPEIT5)

Rare variants carry too little LD to phase in a joint model, and a joint HMM over millions of rare sites does not scale, so SHAPEIT5 splits the problem. Use the full pipeline when N > ~2,000; below that, `phase_common` alone suffices (too few rare-allele carriers for the rare step to add value).

1. **phase_common** phases the common variants (e.g. `--filter-maf 0.001`) into accurate haplotypes - the scaffold. Run per chunk for large chromosomes, with OVERLAPPING regions.
2. **ligate** stitches the per-chunk common scaffolds into one chromosome; chunks must overlap so ligate can resolve phase across the seam (a non-overlapping seam is a guaranteed switch).
3. **phase_rare** takes the FULL genotypes plus the fixed scaffold and places each rare allele onto the already-phased common haplotypes by IBD matching. Do not filter rare variants out of the phase_rare input - placing them is the whole point.

## Switch Error vs Flip vs Hamming -- the Metrics

A single rate hides the failure mode. Report more than one, and look at the distribution of switch positions.

| Metric | What it counts | Inflates on |
|--------|----------------|-------------|
| Switch error rate (SER) | fraction of consecutive het-site pairs whose phase relationship is wrong | many small local errors; the standard headline |
| Flip error | an isolated het phased wrong then immediately corrected (two switches one site apart) | noisy single sites; double-counts in raw SER |
| Hamming error | fraction of het sites on the wrong haplotype under the best global alignment | a few LARGE block swaps - high Hamming, low switch count |
| Long switch / block flip | a sustained segment on the wrong haplotype | poor long-range LD; ruinous for cis/trans yet only 2 switches |

SER and Hamming measure different sins: many tiny flips give high SER but modest Hamming; one half-chromosome block swap gives catastrophic Hamming but only two switches. Het density matters too - SER is per-het-pair, so sparse het sites mean the same SER spans more bp. Typical magnitudes (order-of-magnitude, dataset-specific): Eagle2 + HRC reference, European array ~1.36%; Eagle2 within-cohort N~5,000 ~1.5%; within-cohort N~150,000 (UK Biobank) ~0.27-0.35%; SHAPEIT5 for a variant in ~1 of 100,000 < ~5%. The pattern: common-variant phasing in a big cohort is sub-1%; rare-variant phasing is single-digit-percent at best and worsens steeply as MAC approaches 1.

## Reference-Based vs Within-Cohort

Reference-based phasing wins when the cohort is small (a few thousand samples cannot model LD as well as a 32k-100k+ haplotype panel); phase against the biggest ancestry-matched panel available (Eagle2 `--vcfRef`). Within-cohort phasing wins when the cohort is large and ancestry-matched to itself, because accuracy rises monotonically with N; by UK-Biobank scale within-cohort is more accurate than any external panel. The crossover is in the tens of thousands. Ancestry match dominates either way - a mismatched panel phases worse than a smaller matched one or within-cohort -> reference-panels.

## Per-Method Failure Modes

### Genome-wide SER trusted for a rare-variant call
**Trigger:** quoting one switch-error rate and treating all haplotypes as equally trustworthy. **Mechanism:** SER is dominated by easy common sites; rare-variant phase is far worse and MAC-dependent. **Symptom:** a confident compound-het (cis/trans) call from a small-cohort statistical phase that is actually near chance. **Fix:** stratify accuracy by MAC; confirm rare-variant cis/trans with a trio or read-backed phase.

### Wrong-build or flat genetic map
**Trigger:** a GRCh37 map on GRCh38 data, or a uniform map "for simplicity". **Mechanism:** the map sets the HMM's recombination (transition) rates; wrong coordinates or a flat rate mis-place where haplotype breaks are expected. **Symptom:** degraded phasing, more long switches, no error message. **Fix:** use the build-matched per-chromosome map shipped with the tool; the default population map is right.

### Non-overlapping ligate seam
**Trigger:** chunking a chromosome with abutting (non-overlapping) regions. **Mechanism:** ligate needs overlap to resolve the phase relationship across the seam. **Symptom:** a guaranteed switch at every chunk boundary. **Fix:** make `--region` / `--input-region` / `--scaffold-region` overlap between adjacent chunks.

### chrX male coded diploid
**Trigger:** phasing male chrX non-PAR as diploid heterozygous. **Mechanism:** males are haploid outside the PARs; a het call there is biologically impossible. **Symptom:** corrupted male chrX phase. **Fix:** pass the male sample list (SHAPEIT5 `--haploids`; Eagle handles mixed ploidy); keep PAR1/PAR2 as separate diploid regions with build-correct coordinates.

### Multiallelic records fed to a phaser
**Trigger:** phasing raw multiallelic sites. **Mechanism:** phasers expect biallelic records; a multiallelic record is undefined behavior. **Symptom:** tool errors or mis-phased sites. **Fix:** `bcftools norm -m -any` to split and left-align first -> variant-calling/variant-normalization.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `--filter-maf 0.001` defines the common/rare scaffold split | SHAPEIT5 docs | common variants build the accurate scaffold; rarer variants are phased onto it |
| Use phase_common -> ligate -> phase_rare when N > ~2,000 | SHAPEIT5 docs | below that, too few rare-allele carriers for the rare step to help |
| Report SER stratified by MAC, not genome-wide | Hofmeister 2023 *Nat Genet* 55:1243 | phasing quality is a steep function of MAC; a single number hides rare-variant failure |
| Eagle2 `--Kpbwt` default 10000 (raise at large N) | Loh 2016 *Nat Genet* 48:1443 | more conditioning haplotypes raise accuracy at biobank scale |
| phase_rare `--effective-size` ~15000 (verify) | SHAPEIT5 docs | Ne sets expected recombination; often tuned per dataset, confirm with --help |
| Genetic map must match the data build | Delaneau 2019 *Nat Commun* 10:5436 | a build-mismatched map mis-assigns recombination rates silently |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Switch at every chunk boundary | non-overlapping ligate seams | overlap adjacent chunk regions |
| Corrupted male chrX phase | male non-PAR coded diploid | pass `--haploids`; split PAR/nonPAR |
| Phaser errors on some sites | multiallelic records | `bcftools norm -m -any` first |
| Rare-variant cis/trans call does not replicate | small-cohort statistical phase of rare variants | use SHAPEIT5 at scale; confirm with trio/read-backed |
| Phasing mysteriously bad in one region | wrong-build or flat genetic map | build-match the map |
| SHAPEIT4 syntax fails under SHAPEIT5 | SHAPEIT5 split into phase_common/phase_rare/ligate | use the suite binaries, not a single `shapeit` |
| Beagle OutOfMemoryError | JVM heap too small / whole genome in one job | raise `-Xmx`; phase per chromosome |

## References

- Hofmeister RJ, Ribeiro DM, Rubinacci S, Delaneau O. 2023. Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. *Nat Genet* 55:1243-1249.
- Delaneau O, Zagury JF, Robinson MR, Marchini JL, Dermitzakis ET. 2019. Accurate, scalable and integrative haplotype estimation. *Nat Commun* 10:5436.
- Loh PR, Danecek P, Palamara PF, et al. 2016. Reference-based phasing using the Haplotype Reference Consortium panel. *Nat Genet* 48:1443-1448.
- Browning BL, Tian X, Zhou Y, Browning SR. 2021. Fast two-stage phasing of large-scale sequence data. *Am J Hum Genet* 108:1880-1890.
- Durbin R. 2014. Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT). *Bioinformatics* 30:1266-1272.
- Patterson M, Marschall T, Pisanti N, et al. 2015. WhatsHap: weighted haplotype assembly for future-generation sequencing reads. *J Comput Biol* 22:498-509.
- Li N, Stephens M. 2003. Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data. *Genetics* 165:2213-2233.

## Related Skills

- reference-panels - Select the ancestry-matched panel that reference-based phasing copies from
- genotype-imputation - Imputation consumes the phased haplotypes (pre-phasing)
- imputation-qc - Switch-error benchmarking sits alongside imputation quality QC
- long-read-sequencing/haplotype-phasing - Read-backed / molecular single-sample phasing (a different signal)
- variant-calling/variant-normalization - Split multiallelics and left-align before phasing
- causal-genomics/fine-mapping - Phased haplotypes feed haplotype-level fine-mapping
- clinical-databases/hla-typing - HLA typing is a high-stakes consumer of long-range phase
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate
