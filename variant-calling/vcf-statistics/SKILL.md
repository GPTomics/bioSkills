---
name: bio-vcf-statistics
description: Generate variant statistics, sample concordance, and quality metrics using bcftools stats and gtcheck. Use when evaluating variant quality, comparing samples, or summarizing VCF contents.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# VCF Statistics

Generate statistics and quality metrics using bcftools.

## Statistics Tools

| Command | Purpose |
|---------|---------|
| `bcftools stats` | Comprehensive variant statistics |
| `bcftools gtcheck` | Sample concordance and relatedness |
| `bcftools query` | Custom summaries |

## bcftools stats

**Goal:** Generate comprehensive variant statistics including counts, Ti/Tv ratio, and quality distributions.

**Approach:** Run bcftools stats and parse section-tagged output lines (SN, TSTV, AF, QUAL, DP).

**"How many variants are in this VCF?"** → Compute summary counts, substitution types, and quality distributions from variant records.

### Basic Statistics

```bash
bcftools stats input.vcf.gz > stats.txt
```

### View Key Metrics

```bash
bcftools stats input.vcf.gz | grep "^SN"
```

Output sections:
- `SN` - Summary numbers
- `TSTV` - Transitions/transversions
- `SiS` - Singleton stats
- `AF` - Allele frequency distribution
- `QUAL` - Quality distribution
- `IDD` - Indel distribution
- `ST` - Substitution types
- `DP` - Depth distribution

### Summary Numbers (SN)

```bash
bcftools stats input.vcf.gz | grep "^SN" | cut -f3-
```

Reports:
- Number of samples
- Number of records
- Number of SNPs
- Number of indels
- Number of multiallelic sites
- Number of multiallelic SNPs

### Transition/Transversion Ratio

```bash
bcftools stats input.vcf.gz | grep "^TSTV"
```

Transitions (purine-to-purine: A↔G, or pyrimidine-to-pyrimidine: C↔T) are chemically favored over transversions because they preserve the purine/pyrimidine ring structure. CpG deamination (methylated C→T) is the single most common point mutation in vertebrate genomes and is a transition, which further inflates the Ti/Tv ratio. Exomes have higher Ti/Tv than whole genomes because coding regions are enriched for CpG dinucleotides.

Expected Ti/Tv ratio:
- Whole genome: ~2.0-2.1
- Exome: ~2.8-3.3
- Ti/Tv below expected range suggests excess false-positive SNPs (random errors produce Ti/Tv ~0.5)
- Ti/Tv above expected range suggests over-filtering that disproportionately removes transversions

### Per-Sample Statistics

```bash
bcftools stats -s - input.vcf.gz > per_sample.txt
```

### Compare Two VCFs

```bash
bcftools stats input1.vcf.gz input2.vcf.gz > comparison.txt
```

### Region-Specific Stats

```bash
bcftools stats -r chr1:1000000-2000000 input.vcf.gz > region_stats.txt
bcftools stats -R exome.bed input.vcf.gz > exome_stats.txt
```

## Plotting Statistics

**Goal:** Visualize variant statistics as publication-quality plots.

**Approach:** Pipe bcftools stats output to plot-vcfstats to generate PDF and PNG plots.

### Generate Plots

```bash
bcftools stats input.vcf.gz > stats.txt
plot-vcfstats -p output_dir stats.txt
```

Creates:
- `output_dir/summary.pdf`
- Individual PNG files

### Comparison Plots

```bash
bcftools stats file1.vcf.gz file2.vcf.gz > comparison.txt
plot-vcfstats -p comparison_dir comparison.txt
```

## Expected QC Metric Ranges

Interpreting variant statistics requires context-dependent thresholds. A metric that looks normal in WGS may be alarming in WES.

| Metric | WGS Expected | WES Expected | Flag If |
|--------|-------------|-------------|---------|
| Ti/Tv ratio | ~2.0-2.1 | ~2.8-3.3 | WGS: <1.8 or >2.5; WES: <2.5 or >3.5 |
| Het/Hom ratio | 1.5-2.0 | 1.5-2.0 | Outlier relative to cohort |
| Call rate per variant | >95% | >95% | <90% |
| Call rate per sample | >95% | >95% | <90% |
| Singleton rate (WGS) | ~100-200k | Variable | >2x cohort mean |
| Mendelian error rate (trios) | <0.5% | <0.5% | >1% |

Interpretation notes:
- Het/Hom ratio too high suggests contamination (mixed DNA inflates heterozygosity); too low suggests inbreeding, population structure, or excess homozygous reference calls from low coverage
- Excess singletons per sample suggest sample quality issues, contamination, or library artifacts
- Call rate thresholds apply to common variants; rare variants naturally have higher missingness
- These ranges assume human diploid germline calling; somatic, polyploid, or non-model organisms require different expectations

## Stratified Evaluation

A variant caller with 99% overall accuracy may perform at 70% in difficult genomic regions. Always evaluate stratified by region complexity.

```bash
bcftools stats -R easy_regions.bed input.vcf.gz > easy_stats.txt
bcftools stats -R difficult_regions.bed input.vcf.gz > difficult_stats.txt
```

Key stratification categories (using GIAB stratification BED files):
- **Easy/high-confidence regions** - baseline accuracy measurement
- **Homopolymer runs** - systematic indel errors, especially for Illumina and Ion Torrent
- **Tandem repeats / low-complexity** - alignment ambiguity inflates both FP and FN
- **Segmental duplications** - paralogous mapping produces false heterozygous calls
- **High GC (>70%) / Low GC (<25%)** - coverage bias creates systematic missingness
- **MHC / centromeric regions** - extreme polymorphism or repetitiveness defeats standard callers

GIAB stratification BED files are available at github.com/genome-in-a-bottle/genome-stratifications. Comparing Ti/Tv ratio between easy and difficult regions is a quick diagnostic: a drop in Ti/Tv within difficult regions confirms elevated false-positive rates there.

## Population-Scale QC Metrics

For multi-sample VCFs, per-sample and per-site population metrics reveal systematic issues invisible in single-sample analyses.

**Per-sample checks:**
- Variant count should be roughly consistent across samples from the same population (~4.5M SNPs per WGS sample in humans); outliers suggest processing errors or contamination
- Per-sample Het/Hom ratio: outliers suggest contamination (high) or sample swaps between populations (unexpected value)
- Excess singletons in one sample relative to the cohort mean suggests library or sequencing issues

**Per-site checks:**
- Excess heterozygosity (InbreedingCoeff < -0.3 in GATK, or HWE test): suggests genotyping error at sites where reads from paralogous regions pile up
- Hardy-Weinberg equilibrium: sites failing HWE (p < 1e-5) within a homogeneous population may be genotyping artifacts. HWE filtering must be applied within populations, never across mixed-ancestry cohorts where true allele frequency differences violate HWE assumptions

```bash
bcftools query -f '%CHROM\t%POS\t%INFO/InbreedingCoeff\n' input.vcf.gz | \
    awk '$3 < -0.3' > excess_het_sites.txt
```

## bcftools gtcheck

**Goal:** Verify sample identity and detect sample swaps by comparing genotype concordance.

**Approach:** Use bcftools gtcheck to compute pairwise discordance rates; interpret thresholds based on expected relatedness.

```bash
bcftools gtcheck -G 1 input.vcf.gz > relatedness.txt
bcftools gtcheck -g reference.vcf.gz query.vcf.gz
```

Discordance thresholds for interpreting results:
- Same individual (replicates/re-extractions): <0.5% discordance
- First-degree relatives (parent-child, siblings): ~10-15% discordance
- Unrelated individuals: >10% discordance (typically 15-25%)
- For trios: Mendelian error rate should be <0.5%; rates >1% suggest sample swaps or contamination

Output format (DC lines):

```
DC  0  sample1  sample2  0.95  1234  1200
```

Fields: DC tag, index, sample1, sample2, discordance rate, sites compared, discordant sites. Sample swaps are one of the most common errors in genomics studies. Running gtcheck is essential for any multi-sample study and should be performed early in the QC pipeline before downstream analysis.

## Quick Statistics with Query

**Goal:** Compute targeted summary statistics using bcftools query and shell tools.

**Approach:** Extract specific fields with bcftools query/view and aggregate with awk for counts, means, and distributions.

### Count Variants

```bash
bcftools view -H input.vcf.gz | wc -l
```

### Count by Type

```bash
# SNPs
bcftools view -v snps -H input.vcf.gz | wc -l

# Indels
bcftools view -v indels -H input.vcf.gz | wc -l
```

### Count PASS Variants

```bash
bcftools view -f PASS -H input.vcf.gz | wc -l
```

### Quality Distribution

```bash
bcftools query -f '%QUAL\n' input.vcf.gz | \
    awk '{sum+=$1; count++} END {print "Mean QUAL:", sum/count}'
```

### Depth Distribution

```bash
bcftools query -f '%INFO/DP\n' input.vcf.gz | \
    awk '{sum+=$1; count++} END {print "Mean DP:", sum/count}'
```

### Genotype Counts

```bash
# Count heterozygous sites per sample
bcftools query -f '[%GT\t]\n' input.vcf.gz | \
    awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="0/1" || $i=="0|1") het[i]++}
        END {for(i in het) print "Sample", i, "het:", het[i]}'
```

### Allele Frequency Spectrum

```bash
bcftools query -f '%INFO/AF\n' input.vcf.gz | \
    awk '{
        if($1<0.01) rare++
        else if($1<0.05) low++
        else if($1<0.5) common++
        else freq++
    } END {
        print "Rare (<1%):", rare
        print "Low (1-5%):", low
        print "Common (5-50%):", common
        print "Frequent (>50%):", freq
    }'
```

## Sample Statistics

**Goal:** Compute per-sample variant counts, genotype distributions, and missingness rates.

**Approach:** Use bcftools query/view/stats per sample to tabulate sample-level metrics.

### List and Count Samples

```bash
bcftools query -l input.vcf.gz
bcftools query -l input.vcf.gz | wc -l
```

### Per-Sample Variant Counts

```bash
for sample in $(bcftools query -l input.vcf.gz); do
    count=$(bcftools view -s "$sample" -H input.vcf.gz | \
        bcftools view -c 1 -H | wc -l)
    echo "$sample: $count"
done
```

### Missing Genotypes per Sample

```bash
bcftools stats -s - input.vcf.gz | grep "^PSC"
```

## cyvcf2 Statistics

**Goal:** Compute variant statistics programmatically in Python for custom analyses.

**Approach:** Iterate variants with cyvcf2, accumulate counts by type, and compute quality/genotype distributions with numpy.

### Basic Counts

```python
from cyvcf2 import VCF

stats = {'snps': 0, 'indels': 0, 'other': 0}

for variant in VCF('input.vcf.gz'):
    if variant.is_snp:
        stats['snps'] += 1
    elif variant.is_indel:
        stats['indels'] += 1
    else:
        stats['other'] += 1

print(f'SNPs: {stats["snps"]}')
print(f'Indels: {stats["indels"]}')
print(f'Other: {stats["other"]}')
```

### Quality Statistics

```python
from cyvcf2 import VCF
import numpy as np

quals = []
for variant in VCF('input.vcf.gz'):
    if variant.QUAL:
        quals.append(variant.QUAL)

quals = np.array(quals)
print(f'Mean QUAL: {np.mean(quals):.1f}')
print(f'Median QUAL: {np.median(quals):.1f}')
print(f'Min QUAL: {np.min(quals):.1f}')
print(f'Max QUAL: {np.max(quals):.1f}')
```

### Genotype Distribution

```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
samples = vcf.samples

hom_ref = [0] * len(samples)
het = [0] * len(samples)
hom_alt = [0] * len(samples)
missing = [0] * len(samples)

for variant in vcf:
    for i, gt in enumerate(variant.gt_types):
        if gt == 0:
            hom_ref[i] += 1
        elif gt == 1:
            het[i] += 1
        elif gt == 3:
            hom_alt[i] += 1
        else:
            missing[i] += 1

for i, sample in enumerate(samples):
    print(f'{sample}: HOM_REF={hom_ref[i]}, HET={het[i]}, HOM_ALT={hom_alt[i]}, MISS={missing[i]}')
```

### Transition/Transversion Calculation

```python
from cyvcf2 import VCF

transitions = 0
transversions = 0

ti_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

for variant in VCF('input.vcf.gz'):
    if not variant.is_snp:
        continue
    ref = variant.REF
    alt = variant.ALT[0]
    if (ref, alt) in ti_pairs:
        transitions += 1
    else:
        transversions += 1

ratio = transitions / transversions if transversions > 0 else 0
print(f'Transitions: {transitions}')
print(f'Transversions: {transversions}')
print(f'Ti/Tv ratio: {ratio:.2f}')
```

## Common Workflows

### Quality Control Report

```bash
bcftools stats input.vcf.gz > stats.txt
grep "^SN" stats.txt | cut -f3-
grep "^TSTV" stats.txt | cut -f5
plot-vcfstats -p qc_plots stats.txt
```

### Compare Before/After Filtering

```bash
bcftools stats raw.vcf.gz filtered.vcf.gz > comparison.txt
grep "^SN" comparison.txt | head -20
```

## Quick Reference

| Task | Command |
|------|---------|
| Full stats | `bcftools stats input.vcf.gz` |
| Summary only | `bcftools stats input.vcf.gz \| grep "^SN"` |
| Ti/Tv ratio | `bcftools stats input.vcf.gz \| grep "^TSTV"` |
| Per-sample | `bcftools stats -s - input.vcf.gz` |
| Compare VCFs | `bcftools stats file1.vcf.gz file2.vcf.gz` |
| Sample check | `bcftools gtcheck -G 1 input.vcf.gz` |
| Plot stats | `plot-vcfstats -p dir stats.txt` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `No data` | Empty VCF | Check if VCF has variants |
| `plot-vcfstats not found` | Not installed | Install with bcftools |
| `Cannot open` | Invalid VCF | Check file format |

## Related Skills

- variant-calling/vcf-basics - View and query VCF files
- variant-calling/filtering-best-practices - Evaluate filter impact on QC metrics
- variant-calling/vcf-manipulation - Compare and merge call sets
- variant-calling/variant-calling - Assess calling quality
- variant-calling/joint-calling - Cohort-level genotyping where population metrics apply
- alignment-files/bam-statistics - Upstream alignment QC that affects variant statistics
