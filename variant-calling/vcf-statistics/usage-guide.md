# VCF Statistics Usage Guide

## Overview

This guide covers computing and interpreting VCF quality-control metrics. The core idea: no QC metric has a universal pass value. Each metric's expected range depends on the assay (WGS vs WES), the sample's ancestry, and the cohort it sits in, so QC means reading a deviation for its mechanistic meaning and acting on it, not checking a fixed threshold. It also covers cohort identity QC (sample swaps, contamination, relatedness, sex), which is mandatory before any association or burden analysis.

## Prerequisites

- bcftools (`conda install -c bioconda bcftools`) with plot-vcfstats (requires python + matplotlib)
- vcftools (`conda install -c bioconda vcftools`) for per-sample missingness, exact HWE, and KING-robust kinship
- somalier (`conda install -c bioconda somalier`) and peddy (`pip install peddy`) for scalable identity/relatedness/sex/ancestry QC
- cyvcf2 (`pip install cyvcf2`) and numpy for custom per-record statistics in Python

## Quick Start

Tell your AI agent what you want to do:
- "Generate comprehensive statistics for my VCF including Ti/Tv ratio"
- "My exome Ti/Tv is only 2.1 -- is the callset over-permissive?"
- "Flag het/hom outliers in my cohort accounting for ancestry"
- "Decide whether this HWE deviation is a genotyping error or real biology"
- "Screen my cohort for sample swaps, contamination, and wrong-sex annotations"
- "Compare variant statistics before and after filtering"

## Reading bcftools stats output

`bcftools stats` is the workhorse. Add `-s -` to enable per-sample metrics.

```bash
bcftools stats input.vcf.gz > stats.txt
bcftools stats -s - input.vcf.gz > per_sample.txt
```

Section codes:

| Code | Description |
|------|-------------|
| `SN` | Summary numbers (samples, records, SNPs, indels, multiallelic) |
| `TSTV` | Transition/transversion counts and ratio (ratio in field 5) |
| `SiS` | Singleton stats |
| `AF` | Allele-frequency distribution |
| `QUAL` | Quality-score distribution |
| `IDD` | Indel-length distribution |
| `ST` | Substitution types |
| `DP` | Depth distribution |
| `PSC` | Per-sample counts: hom-ref, het, hom-alt, transitions, transversions, missing |
| `PSI` | Per-sample indels |

```bash
bcftools stats input.vcf.gz | grep "^SN"   | cut -f3-
bcftools stats input.vcf.gz | grep "^TSTV" | cut -f5
bcftools stats -s - input.vcf.gz | grep "^PSC"
```

## Interpreting the key metrics

### Ti/Tv ratio

Expected: WGS ~2.0-2.1, WES ~3.0-3.3 (conventions; they shift with capture kit, ancestry, reference, and caller). WES runs higher because coding regions are CpG- and constraint-enriched and transitions at CpGs plus codon degeneracy inflate the ratio. A random error spectrum gives Ti/Tv ~0.5, so false positives pull the ratio down. A WES callset reporting ~2.1 signals filtering that is too loose. Compute Ti/Tv on the novel (not-in-dbSNP) fraction separately, since false positives concentrate there.

### Het/hom ratio and ancestry

The per-sample het:non-ref-hom ratio is strongly ancestry-dependent because it reflects divergence from the reference: African ancestry commonly ~2.0+, European ~1.5-1.6. A single global cutoff is wrong. Infer ancestry first (peddy/somalier project onto 1000 Genomes PCs), then flag outliers within each ancestry group. Elevated het/hom versus same-ancestry peers indicates contamination or reference bias; depressed het/hom indicates inbreeding, a run of homozygosity, or chromosome loss.

### Novel/known via dbSNP

Overlap with dbSNP gives a false-positive signal independent of Ti/Tv. Stratify novel% by allele frequency: novel COMMON variants are almost always artifacts, while novel rare/singleton variants are expected. A single unstratified novel% is uninformative.

```bash
bcftools annotate -a dbsnp.vcf.gz -c ID input.vcf.gz -Oz -o annotated.vcf.gz
bcftools view -H annotated.vcf.gz | awk '{n++; if($3==".") novel++} END{print "novel frac:", novel/n}'
```

### Missingness / call rate and the iteration trap

Only `./.` counts as missing; `0/0` is a confident hom-ref call and does not. Sample-level and variant-level missingness are coupled, so drop the worst samples, recompute per-variant missingness, drop the worst variants, and repeat rather than applying both thresholds at once. Apply genotype-level filters (set `./.` where GQ<20 or DP<8) BEFORE computing cohort missingness or HWE, or low-quality genotypes drive both.

```bash
vcftools --gzvcf input.vcf.gz --missing-indv --out sample_miss   # .imiss
vcftools --gzvcf input.vcf.gz --missing-site --out site_miss     # .lmiss
```

### HWE: excess-het only

A naive two-sided HWE gate throws away real biology. Filter on EXCESS heterozygosity only (the signature of collapsed paralogs/CNVs manufacturing spurious hets); heterozygote DEFICIT is often real (Wahlund substructure, inbreeding, selection, null alleles). Compute the exact test within an ancestry-homogeneous subgroup, and in case/control studies use controls only (a true association produces HWE deviation in cases).

```bash
vcftools --gzvcf controls.vcf.gz --hardy --out hwe               # exact test P per site
bcftools query -f '%CHROM\t%POS\t%INFO/ExcessHet\n' input.vcf.gz | awk '$3>54.69'
```

### Contamination signatures

Contamination shows up in VCF stats as a het allele-balance distribution shifted away from 0.5, an elevated het count and het/hom ratio, and a depressed novel-fraction Ti/Tv. These are signals, not the measurement. Run VerifyBamID2 or CHARR to estimate the contamination fraction alpha (alpha > ~0.02-0.03 is a red flag; somatic pipelines are sensitive to 1%).

```bash
bcftools query -i 'GT="het"' -f '[%AD]\n' input.vcf.gz | \
    awk -F',' '{ab=$2/($1+$2); s+=ab; n++} END{print "mean het AB:", s/n}'   # expect ~0.5
```

## Identity QC: sample swaps, relatedness, sex

Sample swaps are among the most common errors in sequencing studies. Build the all-pairs relatedness matrix and reconcile it against the manifest before any analysis. Swaps appear as off-diagonal surprises (unexpected relatedness) or on-diagonal failures (a "replicate" that is not).

### bcftools gtcheck (quick, native)

Modern gtcheck cross-checks all samples in one file when `-g` is omitted; the old `-G/--GTs-only` flag was removed after the 1.10 rewrite. Output DC lines carry the query sample, genotyped sample, an accumulated discordance score (not a rate), and the number of sites compared; smaller discordance means more similar samples.

```bash
bcftools gtcheck input.vcf.gz > gtcheck.txt          # all-pairs cross-check
bcftools gtcheck -g reference.vcf.gz query.vcf.gz    # concordance to a genotyping panel
```

### peddy (PED-vs-VCF)

peddy samples ~25000 sites plus chrX to check reported sex, relationships, and ancestry (PCA projection onto 1000 Genomes) directly from a VCF and its PED file.

```bash
python -m peddy -p 4 --plot --prefix cohort_qc input.vcf.gz cohort.ped
```

### somalier (scalable)

somalier extracts tiny per-sample sketches at informative sites, then relates them; it scales to tens of thousands of samples in seconds and cross-checks RNA-seq against WGS from the same individual.

```bash
somalier extract -d extracted/ --sites sites.vcf.gz -f ref.fa input.vcf.gz
somalier relate --ped cohort.ped extracted/*.somalier
somalier ancestry --labels 1kg-labels.tsv 1kg/*.somalier ++ extracted/*.somalier   # labelled ++ query
```

### KING kinship via vcftools

```bash
vcftools --gzvcf input.vcf.gz --relatedness2 --out kin    # Manichaikul KING-robust method
```

KING kinship bands: >0.354 duplicate/MZ twin, [0.177, 0.354] first-degree, [0.0884, 0.177] second-degree, [0.0442, 0.0884] third-degree. A pair not expected to be related at ~0.5 is a swap or duplicate.

## Stratified and comparative evaluation

### Stratify by region class

A caller at 99% overall accuracy may drop to 70% in hard regions, so evaluate stratified using GIAB stratification BED files (github.com/genome-in-a-bottle/genome-stratifications). A drop in Ti/Tv within difficult regions confirms elevated false positives there.

```bash
bcftools stats -R easy_regions.bed      input.vcf.gz > easy.txt
bcftools stats -R difficult_regions.bed input.vcf.gz > difficult.txt
```

### Compare before and after filtering

```bash
bcftools stats raw.vcf.gz filtered.vcf.gz > comparison.txt
grep "^SN" comparison.txt | head -20
plot-vcfstats -p comparison_plots comparison.txt
```

## Custom statistics in Python (cyvcf2)

The `examples/vcf_stats.py` script computes counts, Ti/Tv, and mean QUAL in one pass. Two common extensions:

### Per-sample genotype distribution

```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
samples = vcf.samples
hom_ref = [0] * len(samples); het = [0] * len(samples)
hom_alt = [0] * len(samples); missing = [0] * len(samples)

for variant in vcf:
    for i, gt in enumerate(variant.gt_types):   # cyvcf2 codes: 0 hom-ref, 1 het, 2 unknown, 3 hom-alt
        if gt == 0:
            hom_ref[i] += 1
        elif gt == 1:
            het[i] += 1
        elif gt == 3:
            hom_alt[i] += 1
        else:
            missing[i] += 1

for i, s in enumerate(samples):
    ratio = het[i] / hom_alt[i] if hom_alt[i] else 0
    print(f'{s}: het/hom={ratio:.2f} HET={het[i]} HOM_ALT={hom_alt[i]} MISS={missing[i]}')
```

### Allele-frequency spectrum

```python
from cyvcf2 import VCF

bins = {'rare(<1%)': 0, 'low(1-5%)': 0, 'common(5-50%)': 0, 'frequent(>50%)': 0}
for variant in VCF('input.vcf.gz'):
    af = variant.INFO.get('AF')
    if af is None:
        continue
    af = af[0] if isinstance(af, tuple) else af
    if af < 0.01:
        bins['rare(<1%)'] += 1
    elif af < 0.05:
        bins['low(1-5%)'] += 1
    elif af < 0.5:
        bins['common(5-50%)'] += 1
    else:
        bins['frequent(>50%)'] += 1
bins
```

## Generating plots

```bash
bcftools stats input.vcf.gz > stats.txt
plot-vcfstats -p output_directory stats.txt
```

Creates a multi-page `summary.pdf` and individual PNGs (substitution types, indel-length distribution, Ts/Tv by quality, per-sample SNPs/indels, quality and depth distributions). Requires matplotlib.

## Example Prompts

### Callset quality
> "Generate comprehensive statistics for my VCF including Ti/Tv ratio and flag whether the ratio is in range for an exome"

> "Compute Ti/Tv on the novel (not-in-dbSNP) fraction separately to check for false positives"

### Cohort QC
> "Flag het/hom ratio outliers in my cohort after inferring ancestry, not against a global cutoff"

> "Decide whether these HWE-failing sites are genotyping errors or real, filtering on excess heterozygosity within ancestry"

> "Screen my cohort for sample swaps, contamination, and wrong-sex annotations with somalier and peddy"

### Comparison
> "Compare variant statistics before and after filtering and quantify how many SNPs and indels were removed"

## What the Agent Will Do

1. Identify the VCF/BCF and check indexing, then run `bcftools stats` (per-sample, region-based, or comparison mode as needed)
2. Extract key metrics -- counts, Ti/Tv (overall and novel), per-sample het/hom, depth, missingness
3. Contextualize each metric against the assay and, where relevant, inferred ancestry rather than an absolute threshold
4. For cohorts, run identity QC (gtcheck/peddy/somalier/KING) and reconcile relatedness against the manifest
5. Flag issues (low Ti/Tv, ancestry-adjusted het/hom outliers, excess-het sites, contamination signatures, unexpected relatedness) with the recommended action

## Tips

- Ti/Tv ~2.0-2.1 (WGS) and ~3.0-3.3 (WES) indicate good calls; below range signals false positives, since random errors have Ti/Tv ~0.5
- Never apply a single global het/hom cutoff across ancestries; stratify first
- Apply genotype-level filters before computing cohort missingness or HWE, or garbage genotypes drive both
- Filter HWE on excess heterozygosity only, within ancestry, in controls; heterozygote deficit is often real
- The het allele-balance distribution shifting off 0.5 is a contamination signal; confirm with VerifyBamID2/CHARR on the BAM
- Run identity QC (swap/relatedness/sex) early -- one undetected swap can create or erase a significant hit
- plot-vcfstats requires matplotlib; install it separately if plots fail

## Related Skills

- variant-calling/filtering-best-practices - Apply the filters these metrics motivate
- variant-calling/gatk-variant-calling - Feed VerifyBamID2 contamination into calling
- variant-calling/variant-normalization - Normalize before comparing or annotating call sets
- variant-calling/vcf-basics - Query and understand VCF fields
- variant-calling/joint-calling - Cohort genotyping where population QC applies
