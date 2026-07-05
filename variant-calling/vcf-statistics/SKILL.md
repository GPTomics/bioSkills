---
name: bio-vcf-statistics
description: Compute and interpret VCF quality-control metrics (Ti/Tv, het/hom, novel/known, missingness, HWE, contamination, relatedness) with bcftools stats, vcftools, plot-vcfstats, and identity tools (somalier, peddy, KING). Use when judging whether a callset is trustworthy, diagnosing a low Ti/Tv or outlier het/hom sample, deciding whether an HWE deviation is error or biology, screening a cohort for sample swaps/contamination/wrong-sex before analysis, or comparing call sets before and after filtering. Not for applying filters (see variant-calling/filtering-best-practices) or normalizing representation (see variant-calling/variant-normalization).
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, vcftools 0.1.16+, somalier 0.2.19+, peddy 0.4.8+, cyvcf2 0.30+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: `bcftools gtcheck` was rewritten around bcftools 1.10; the old `-G/--GTs-only` flag is gone. Modern gtcheck cross-checks all samples in one file when `-g` is omitted, and uses `-E/--error-probability`. Verify flags against the installed build.

# VCF Statistics

**"Is this callset any good, and are the samples who the manifest says they are?"** -> Summarize variant counts and quality distributions, then read each QC metric as a signal whose expected value is context-dependent.

- CLI: `bcftools stats` (+ `plot-vcfstats`), `vcftools`, `somalier`, `peddy`
- Python: `cyvcf2` for custom per-record statistics

## The governing principle

No QC metric has a universal pass value. Each metric's expected range depends on the assay (WGS vs WES), the sample's ancestry, and the cohort it sits in, so QC is not threshold-checking but reading a deviation for its mechanistic meaning and acting on it. Three rules follow. First, always compare a sample against a matched cohort (same assay, same inferred ancestry), never against an absolute number. Second, order of operations matters: apply genotype-level filters (set `./.` where GQ/DP/allele-balance fail) BEFORE computing cohort missingness, call rate, or HWE, because low-quality genotypes left in the matrix drive spurious missingness and HWE deviation. Third, sample-identity QC is a graph problem, not a per-sample check: build the all-pairs relatedness matrix and reconcile it against the declared pedigree/manifest, and do this BEFORE any association or burden analysis, because one undetected swap or contaminated sample can fabricate or erase a genome-wide-significant hit.

## QC metric decision table

Expected values are stated conventions (they shift with capture kit, ancestry, reference build, and caller); treat them as starting points and confirm against a matched cohort.

| Metric | Expected (WGS) | Expected (WES) | A deviation MEANS | Action |
|--------|---------------|---------------|-------------------|--------|
| Ti/Tv (overall) | ~2.0-2.1 | ~3.0-3.3 | Low -> false-positive transversions dilute the signal (random errors have Ti/Tv ~0.5); high -> over-filtering removed transversions | Tighten site filters if low; check filter for transversion bias if high |
| Ti/Tv (novel only) | near overall | near overall | Novel Ti/Tv well below the known-site value -> FP contamination in the novel fraction | Raise stringency; the novel set is where FPs concentrate |
| Het/hom ratio | ancestry-dependent (~1.5-1.6 EUR, ~2.0+ AFR) | same | High vs same-ancestry cohort -> contamination or reference bias; low -> inbreeding/consanguinity/ROH or chromosome loss | Stratify by ancestry first; then flag within-ancestry outliers |
| Novel fraction (vs dbSNP) | low for common; high for rare/singleton | same | Novel COMMON variants -> almost always artifacts; high novel rare -> expected or under-studied population | Stratify novel% by frequency; investigate common novels |
| Call rate (per sample) | >95-98% | >95-98% | Low -> low-coverage/low-quality sample | Drop worst samples, then recompute per-variant missingness (iterate) |
| Call rate (per variant) | >95-99% | >95-99% | Low -> site in a hard-to-genotype region | Drop worst variants after sample QC; rare variants tolerate more missing |
| Het allele-balance (hets) | centered on 0.5 | 0.5 | Shifted away from 0.5 + elevated het count -> contamination signature | Confirm with VerifyBamID2 (on the BAM) or CHARR (on the VCF/gVCF) -- the real test |
| Excess-het / HWE | in equilibrium within ancestry | same | EXCESS het -> collapsed paralog/CNV mapping artifact; het DEFICIT -> often real (Wahlund/inbreeding) | Filter on excess het only, within ancestry, in controls; do NOT blanket-filter HWE |
| Relatedness (kinship) | matches manifest | matches manifest | Unexpected high kinship -> sample swap/duplicate; a "replicate" that is not -> mislabel | Reconcile the all-pairs matrix against the pedigree with KING/somalier/peddy |

## bcftools stats

**Goal:** Generate comprehensive variant statistics (counts, Ti/Tv, per-sample het/hom and singletons, indel and depth distributions).

**Approach:** Run bcftools stats and read the section-tagged output lines; add `-s -` for per-sample metrics.

```bash
bcftools stats input.vcf.gz > stats.txt              # cohort-level
bcftools stats -s - input.vcf.gz > per_sample.txt    # per-sample (PSC/PSI lines)
bcftools stats file1.vcf.gz file2.vcf.gz > cmp.txt   # compare two callsets
plot-vcfstats -p qc_plots/ stats.txt                 # render PDF + PNGs (needs matplotlib)
```

Output sections: `SN` summary numbers, `TSTV` transition/transversion, `SiS` singletons, `AF` allele-frequency spectrum, `QUAL` quality distribution, `IDD` indel-length distribution, `ST` substitution types, `DP` depth distribution, `PSC` per-sample counts (hom-ref, het, hom-alt, transitions, transversions, missing), `PSI` per-sample indels.

```bash
bcftools stats input.vcf.gz | grep "^SN"    | cut -f3-   # counts
bcftools stats input.vcf.gz | grep "^TSTV"  | cut -f5     # Ti/Tv ratio
bcftools stats -s - input.vcf.gz | grep "^PSC"           # per-sample het/hom/missing
```

## Ti/Tv ratio

Transitions (purine<->purine A<->G, pyrimidine<->pyrimidine C<->T) are favored over transversions because CpG deamination (methylated C->T) is the single most common vertebrate point mutation and is a transition, and because transitions are more often synonymous. A random error spectrum gives Ti/Tv ~0.5, so a callset diluted with false positives drifts DOWNWARD. WES runs higher than WGS (~3.0-3.3 vs ~2.0-2.1) because coding regions are CpG- and constraint-enriched and transitions at CpGs plus codon degeneracy push the ratio up.

Ti/Tv is the single fastest gross-error smell test. A WES callset reporting Ti/Tv ~2.1 is telling the analyst the filtering is too loose. Stratify: compute Ti/Tv on the novel fraction separately (below) since that is where false positives concentrate. Ancestry and target design shift the exact number, so compare against a matched cohort, not the absolute.

## Het/hom ratio and ancestry

The per-sample het:non-ref-hom ratio reflects heterozygosity relative to the reference, so it is strongly ancestry-dependent: African-ancestry genomes diverge more from GRCh (a mostly European-ancestry assembly) and carry more het calls (commonly ~2.0+), European ancestry sits ~1.5-1.6, and the value shifts across populations. A single global het/hom cutoff is therefore wrong: it would flag every AFR sample in a EUR-tuned pipeline.

Infer ancestry first (peddy/somalier project onto 1000 Genomes PCs), then flag outliers WITHIN each ancestry group. An elevated het/hom versus same-ancestry peers indicates contamination (foreign reads manufacture spurious hets) or reference bias; a depressed het/hom indicates inbreeding/consanguinity, a long run of homozygosity, or chromosome loss.

## Novel/known ratio via dbSNP

Overlapping the callset with dbSNP gives an orthogonal false-positive signal independent of Ti/Tv. Annotate known/novel and stratify by frequency, because the diagnostic differs by allele frequency: novel COMMON variants are almost always artifacts (real common variants are already catalogued), while novel rare/singleton variants are expected and biologically real. A single novel% without frequency stratification is uninformative.

```bash
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz -Oz -o annotated.vcf.gz
# novel fraction = records with ID "." over total; stratify by INFO/AF
bcftools view -H annotated.vcf.gz | awk '{n++; if($3==".") novel++} END{print "novel:", novel/n}'
```

## Missingness and call rate

Call rate is the fraction of sites with a non-missing genotype for a sample; missingness is its complement. Only `./.` (no-call) counts as missing; `0/0` (confident hom-ref) does NOT, so treating no-call as hom-ref biases allele frequencies. GWAS convention drops samples below ~95-98% call rate and variants below ~95-99% (dataset-dependent conventions).

The iteration trap: sample-level and variant-level missingness are coupled, so applying both thresholds in one pass is wrong. Drop the worst samples, recompute per-variant missingness, drop the worst variants, and repeat. Critically, apply genotype-level filters (`./.` where GQ<20 or DP<8) BEFORE computing cohort missingness or HWE, or low-quality genotypes will drive both.

```bash
vcftools --gzvcf input.vcf.gz --missing-indv --out sample_miss   # .imiss (per-sample)
vcftools --gzvcf input.vcf.gz --missing-site --out site_miss     # .lmiss (per-site)
```

## HWE filtering (excess-het only)

Hardy-Weinberg testing flags genotype frequencies that deviate from p^2:2pq:q^2, but a naive two-sided HWE gate throws away real biology. The discipline:

- **Filter on EXCESS heterozygosity only.** Heterozygote excess is the signature of a mapping artifact: a duplicated/collapsed region piles reads from two paralogous copies onto one locus, manufacturing spurious hets everywhere. Heterozygote DEFICIT, by contrast, is often real (Wahlund effect from substructure, inbreeding, selection, a true null allele). GATK's `ExcessHet` and `InbreedingCoeff` and gnomAD's `ExcessHet` filter the excess side only.
- **Compute within an ancestry-homogeneous subgroup.** Pooling populations with different allele frequencies induces a Wahlund het-deficit that mimics genotyping error; HWE on a mixed cohort filters good variants.
- **Use the EXACT test, not chi-square.** The chi-square approximation is anticonservative for rare variants and small samples; the standard is the Wigginton-Cutler-Abecasis exact test (implemented in vcftools `--hardy` and PLINK).
- **In case/control studies, compute HWE in CONTROLS only.** A true association at a strong-effect locus produces HWE deviation in cases; filtering it removes the hit.

```bash
vcftools --gzvcf controls.vcf.gz --hardy --out hwe               # exact-test P per site
# GATK ExcessHet (phred): larger = more excess het = more suspect
bcftools query -f '%CHROM\t%POS\t%INFO/ExcessHet\n' input.vcf.gz | awk '$3>54.69'
```

The threshold 54.69 is GATK's default ExcessHet cutoff (phred-scaled p ~= 3.4e-6, ~the 1000-sample z=-4.5 boundary); tune it to the cohort size.

## Contamination signatures

Cross-sample contamination is visible in VCF statistics before any dedicated test: the het allele-balance distribution shifts away from 0.5 (foreign reads add minor-allele support at true hom sites and skew true hets), the het count and het/hom ratio rise, and the novel-fraction Ti/Tv drops. These are SIGNALS, not the measurement. The real test runs on the BAM/CRAM: VerifyBamID2 (Zhang et al. 2020) estimates the contamination fraction alpha ancestry-agnostically by modeling observed allele fractions against population frequencies, and CHARR estimates alpha directly from VCF-level reference-read counts at hom-alt sites. An alpha above ~0.02-0.03 is a red flag; somatic pipelines are sensitive to even 1%. GATK pipelines feed `--contamination alpha` from VerifyBamID2. See variant-calling/gatk-variant-calling for wiring contamination estimates into calling.

```bash
# quick het allele-balance sanity check from AD (het genotypes only)
bcftools query -i 'GT="het"' -f '[%AD]\n' input.vcf.gz | \
    awk -F',' '{ab=$2/($1+$2); s+=ab; n++} END{print "mean het AB:", s/n}'   # expect ~0.5
```

## Sample-swap, relatedness, and sex checks (mandatory cohort QC)

Sample swaps are among the most common errors in sequencing studies, so identity QC is not optional. Build the all-pairs relatedness matrix and reconcile it against the manifest; swaps appear as off-diagonal surprises (unexpected relatedness) or on-diagonal failures (a "replicate" that is not).

| Tool | Input | Detects | Notes |
|------|-------|---------|-------|
| `bcftools gtcheck` | VCF | Same-file swaps/duplicates (quick) | Native, no reference panel; discordance score, not a relatedness graph |
| KING (Manichaikul 2010) | PLINK bed | Robust kinship without allele-freq/ancestry assumptions | Reference method; kinship bands below |
| peddy (Pedersen 2017) | VCF + PED | Reported vs inferred sex, relationships, ancestry (PCA on 1000G) | Fast, VCF-only; ideal PED-vs-VCF reconciliation |
| somalier (Pedersen 2020) | BAM/CRAM/VCF | Relatedness/ancestry/sex at scale; cross-checks RNA-seq vs WGS | Tiny per-sample sketches; tens of thousands of samples in seconds |

KING kinship coefficient bands: >0.354 duplicate/MZ twin, [0.177, 0.354] first-degree (parent-child, full sib), [0.0884, 0.177] second-degree, [0.0442, 0.0884] third-degree. A pair not expected to be related at ~0.5 is a swap or duplicate.

```bash
# bcftools: cross-check all samples in one file (no -g), or against a truth VCF (-g)
bcftools gtcheck input.vcf.gz > gtcheck.txt          # DC lines: query, genotyped, discordance, sites
bcftools gtcheck -g reference.vcf.gz query.vcf.gz    # concordance to a genotyping panel

# peddy: PED-vs-VCF sex/relatedness/ancestry, 4 CPUs, HTML + CSVs
python -m peddy -p 4 --plot --prefix cohort_qc input.vcf.gz cohort.ped

# somalier: extract sketches then relate against the pedigree
somalier extract -d extracted/ --sites sites.vcf.gz -f ref.fa input.vcf.gz
somalier relate --ped cohort.ped extracted/*.somalier    # writes an HTML relatedness report
somalier ancestry --labels 1kg-labels.tsv 1kg/*.somalier ++ extracted/*.somalier    # PCA projection (labelled ++ query)

# vcftools: KING-robust kinship directly from a VCF
vcftools --gzvcf input.vcf.gz --relatedness2 --out kin    # .relatedness2 (Manichaikul method)
```

## Stratified evaluation

A caller with 99% overall accuracy may drop to 70% in difficult regions, so single-number accuracy hides where a callset fails. Evaluate stratified by region class using GIAB stratification BED files (github.com/genome-in-a-bottle/genome-stratifications).

```bash
bcftools stats -R easy_regions.bed      input.vcf.gz > easy.txt
bcftools stats -R difficult_regions.bed input.vcf.gz > difficult.txt
```

Key strata and their failure modes: homopolymer runs (systematic indel errors, Illumina/Ion Torrent), tandem repeats / low-complexity (alignment ambiguity inflates FP and FN), segmental duplications (paralogous mapping -> false hets), high-GC >70% / low-GC <25% (coverage-bias missingness), MHC/centromeric (extreme polymorphism or repetitiveness). A drop in Ti/Tv within difficult regions confirms elevated false positives there.

## Quick counts with query

```bash
bcftools view -H input.vcf.gz | wc -l              # total records
bcftools view -v snps   -H input.vcf.gz | wc -l    # SNPs
bcftools view -v indels -H input.vcf.gz | wc -l    # indels
bcftools view -f PASS   -H input.vcf.gz | wc -l    # PASS variants
bcftools query -f '%QUAL\n' input.vcf.gz | awk '{s+=$1;n++} END{print "mean QUAL:", s/n}'
```

See examples/vcf_stats.py for a cyvcf2 script computing counts, Ti/Tv, and mean QUAL in one pass; the usage guide covers per-sample genotype distributions and the allele-frequency spectrum.

## Quick Reference

| Task | Command |
|------|---------|
| Full stats | `bcftools stats input.vcf.gz` |
| Per-sample het/hom/missing | `bcftools stats -s - input.vcf.gz \| grep "^PSC"` |
| Ti/Tv ratio | `bcftools stats input.vcf.gz \| grep "^TSTV" \| cut -f5` |
| Per-sample missingness | `vcftools --gzvcf in.vcf.gz --missing-indv` |
| Exact HWE (controls) | `vcftools --gzvcf controls.vcf.gz --hardy` |
| KING kinship | `vcftools --gzvcf in.vcf.gz --relatedness2` |
| Sex/ancestry/relatedness | `python -m peddy -p 4 --plot --prefix qc in.vcf.gz in.ped` |
| Scalable identity QC | `somalier extract ...` then `somalier relate --ped ...` |
| Plot stats | `plot-vcfstats -p dir stats.txt` |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Every AFR sample flagged as high het/hom | Single global cutoff across ancestries | Infer ancestry, flag outliers within each group |
| Good variants filtered by HWE | Two-sided HWE on a mixed cohort or on cases | Excess-het only, within ancestry, controls-only |
| Spurious HWE/missingness | Cohort metrics computed before genotype filtering | Set `./.` on low GQ/DP first, then recompute |
| `bcftools gtcheck -G 1` errors | `-G` removed after the 1.10 rewrite | Drop `-G`; cross-check runs by default without `-g` |
| Low Ti/Tv only in novel set | False positives concentrate in novel variants | Raise stringency; recheck against dbSNP overlap |
| `plot-vcfstats not found` / no plots | matplotlib missing or not on PATH | `pip install matplotlib`; check `which plot-vcfstats` |

## Related Skills

- variant-calling/filtering-best-practices - Apply the site and genotype filters these metrics motivate
- variant-calling/gatk-variant-calling - VerifyBamID2 contamination estimate feeding `--contamination`
- variant-calling/variant-normalization - Normalize before comparing or annotating call sets
- variant-calling/vcf-basics - View, query, and understand VCF fields
- variant-calling/vcf-manipulation - Compare and merge call sets
- variant-calling/joint-calling - Cohort genotyping where population QC applies
- alignment-files/bam-statistics - Upstream alignment QC that drives variant statistics

## References

- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. 2021 *GigaScience* 10:giab008. (bcftools stats/gtcheck)
- Danecek P, Auton A, Abecasis G, et al. The variant call format and VCFtools. 2011 *Bioinformatics* 27:2156-2158. (vcftools QC)
- Wigginton JE, Cutler DJ, Abecasis GR. A note on exact tests of Hardy-Weinberg equilibrium. 2005 *American Journal of Human Genetics* 76:887-893. (exact HWE test)
- Manichaikul A, Mychaleckyj JC, Rich SS, et al. Robust relationship inference in genome-wide association studies. 2010 *Bioinformatics* 26:2867-2873. (KING kinship)
- Pedersen BS, Quinlan AR. Who's Who? Detecting and Resolving Sample Anomalies in Human DNA Sequencing Studies with Peddy. 2017 *American Journal of Human Genetics* 100:406-413. (peddy)
- Pedersen BS, Bhetariya PJ, Brown J, et al. Somalier: rapid relatedness estimation for cancer and germline studies using efficient genome sketches. 2020 *Genome Medicine* 12:62. (somalier)
- Zhang F, Flickinger M, Gagliano Taliun SA, et al. Ancestry-agnostic estimation of DNA sample contamination from sequence reads. 2020 *Genome Research* 30:185-194. (VerifyBamID2)
