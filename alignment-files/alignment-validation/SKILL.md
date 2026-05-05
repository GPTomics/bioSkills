---
name: bio-alignment-validation
description: Validate alignment quality with insert size distribution, proper pairing rates, GC bias, strand balance, and other post-alignment metrics. Use when verifying alignment data quality before variant calling or quantification.
tool_type: mixed
primary_tool: samtools
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, picard 3.1+, pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Validation

Post-alignment quality control to verify alignment quality and identify issues.

**"Check alignment quality"** → Compute post-alignment QC metrics (mapping rate, pairing, insert size, strand balance) to identify issues before downstream analysis.
- CLI: `samtools flagstat`, `samtools stats`, Picard `CollectAlignmentSummaryMetrics`
- Python: `pysam.AlignmentFile` iteration with metric calculations

## Two Different Validations

| Concern | Tools | What it catches |
|---------|-------|-----------------|
| **File integrity** | `samtools quickcheck`, `picard ValidateSamFile` | Truncation, missing EOF, malformed records, wrong CIGAR, MAPQ out of range |
| **Sequence dictionary identity** | `samtools dict` + M5 diff | BAM aligned to wrong reference flavor / different decoy / chr vs no-chr |
| **QC metrics** | `samtools stats`, `flagstat`, `mosdepth`, Picard `CollectMultipleMetrics` / `CollectHsMetrics` / `CollectWgsMetrics` | Are the data biologically reasonable for the assay? |
| **Contamination / sample swap** | `verifybamid2`, `somalier`, Picard `CrosscheckFingerprints` | Cross-sample contamination, tumor-normal swap, mislabeled sample |

A file can pass `quickcheck` and still be malformed in ways that crash GATK three hours into HaplotypeCaller. Conversely, a QC-poor BAM can be structurally valid.

### File Integrity
```bash
# Fast: header + EOF block check (misses mid-file truncation, invalid CIGAR)
samtools quickcheck -v in.bam || echo "QUICKCHECK FAILED"
samtools quickcheck -v *.bam > bad_bams.fofn   # one fail-line per bad file

# Slow but thorough: structural validation
picard ValidateSamFile I=in.bam MODE=SUMMARY R=ref.fa

# Production: ignore expected-but-noisy
picard ValidateSamFile I=in.bam MODE=SUMMARY R=ref.fa \
    IGNORE=INVALID_MAPPING_QUALITY \
    IGNORE=MISMATCH_FLAG_MATE_NEG_STRAND
```

CI-safe one-liner:
```bash
test -s in.bam \
  && samtools quickcheck -v in.bam \
  && [ $(samtools view -c -F 2304 in.bam) -gt 1000 ] \
  || { echo "BAM failed integrity"; exit 1; }
```

### Sequence Dictionary Cross-Validation (M5)
```bash
# Compare per-contig MD5 between BAM and reference
diff \
    <(samtools view -H in.bam | grep '^@SQ' | tr '\t' '\n' | grep '^M5:' | sort) \
    <(samtools dict ref.fa | grep '^@SQ' | tr '\t' '\n' | grep '^M5:' | sort)
```

If M5s differ, the BAM was aligned to a different sequence than the current reference (even if contig names match). Concrete failure modes: GRCh38 vs GRCh38.p13 vs GRCh38_no_alt (alt contigs differ); UCSC `chr1` vs Ensembl `1` (names differ, M5s match -- pure renaming); soft-masked vs hard-masked (M5 matches, viewers differ). The M5 tag is the only definitive identity check.

### Contamination and Sample Swap

No alignment QC is complete without these in production:
```bash
# Cross-sample contamination
verifybamid2 --SVDPrefix /resources/1000g.b38.vcf.gz.SVD \
    --Reference ref.fa --BamFile sample.bam --Output sample.contam
# FREEMIX < 0.01 = clean; > 0.03 = problematic; > 0.05 breaks germline calling

# Relatedness, sex check, sample swap detection
somalier extract -d extracted/ -s /resources/sites.GRCh38.vcf.gz \
    -f ref.fa sample.bam
somalier relate --infer extracted/*.somalier

# Tumor/normal pairing verification
picard CrosscheckFingerprints I=tumor.bam I=normal.bam \
    HAPLOTYPE_MAP=Homo_sapiens_assembly38.haplotype_database.txt
# LOD > 5 = same individual; < -5 = different
```

Sample-swap rates of 0.5-1% in production cohorts are typical. Without `somalier` or `CrosscheckFingerprints`, swaps are detected only when a downstream finding contradicts clinical expectation.

## Insert Size Distribution

**Goal:** Verify that the fragment length distribution matches the library preparation protocol.

**Approach:** Extract template_length from properly paired reads and compare the distribution to expected values for the library type.

### samtools stats

```bash
samtools stats input.bam > stats.txt
grep "^IS" stats.txt | cut -f2,3 > insert_sizes.txt
```

### Picard CollectInsertSizeMetrics

```bash
java -jar picard.jar CollectInsertSizeMetrics \
    I=input.bam \
    O=insert_metrics.txt \
    H=insert_histogram.pdf
```

### Expected Insert Sizes by Library

| Library | Mean insert | Distribution shape | Diagnostic |
|---------|-------------|-------------------|-----------|
| TruSeq DNA PCR-free WGS | 400-500 bp | Roughly Gaussian | Sharp peak; bimodality = degraded sample |
| TruSeq DNA Nano (PCR) WGS | 300-400 bp | Gaussian, narrower | |
| Twist / IDT exome capture | 250-350 bp | Gaussian | |
| TruSeq Stranded mRNA | 200-300 bp | Right-skewed (transcript distribution) | Long tail = poor size selection |
| Ribo-Zero rRNA-depleted | 250-400 bp | Right-skewed | |
| Smart-seq2 / Smart-seq3 | 200-700 bp | Broad | |
| 10x Chromium (3') | n/a -- not informative | n/a | |
| TruSeq ChIP | 200-400 bp | Sharp | |
| ATAC-seq (Buenrostro / Omni-ATAC) | Multimodal | Peaks at ~50, ~180, ~340 bp | **Missing multimodal pattern = bad library**; missing ~180 bp = under-digested |
| Hi-C / Micro-C | Multimodal | Peak at ligation-junction size | |
| cfDNA / ctDNA | 160-180 bp | Multimodal; ~167 bp mononucleosomal + ~340 dinuc | Tumor-derived shorter (~145 bp); shape itself is a biomarker |
| FFPE | 100-250 bp | Right-skewed, broad | |
| aDNA | 30-80 bp | Sharp left-skewed | |
| ONT (native) | 1-30 kb | n/a | |
| PacBio HiFi | 10-25 kb | Sharp peak | |

For ATAC, the multimodal pattern *is* the QC. If the mononucleosomal peak (~180 bp) is absent, Tn5 was over-titrated, under-titrated, or DNA was degraded. Use ATACseqQC `fragSizeDist()` for the standard ATAC fragment-size diagnostic.

### Python Insert Size Analysis

```python
import pysam
import numpy as np
import matplotlib.pyplot as plt

def get_insert_sizes(bam_file, max_reads=100000):
    sizes = []
    bam = pysam.AlignmentFile(bam_file, 'rb')
    for i, read in enumerate(bam.fetch()):
        if i >= max_reads:
            break
        if read.is_proper_pair and not read.is_secondary and read.template_length > 0:
            sizes.append(read.template_length)
    bam.close()
    return sizes

sizes = get_insert_sizes('sample.bam')
print(f'Median insert size: {np.median(sizes):.0f}')
print(f'Mean insert size: {np.mean(sizes):.0f}')
print(f'Std dev: {np.std(sizes):.0f}')

plt.hist(sizes, bins=100, range=(0, 1000))
plt.xlabel('Insert Size')
plt.ylabel('Count')
plt.savefig('insert_size_dist.pdf')
```

## Proper Pairing Rate

Percentage of reads correctly paired.

### samtools flagstat

```bash
samtools flagstat input.bam

samtools flagstat input.bam | grep "properly paired"
```

### Calculate Pairing Rate

```bash
proper=$(samtools view -c -f 2 input.bam)
mapped=$(samtools view -c -F 4 input.bam)
rate=$(echo "scale=4; $proper / $mapped * 100" | bc)
echo "Proper pairing rate: ${rate}%"
```

### Expected Rates

| Metric | Good | Marginal | Poor |
|--------|------|----------|------|
| Proper pair | > 90% | 80-90% | < 80% |
| Mapped | > 95% | 90-95% | < 90% |
| Singletons | < 5% | 5-10% | > 10% |

## GC Bias

GC content correlation with coverage.

### Picard CollectGcBiasMetrics

```bash
java -jar picard.jar CollectGcBiasMetrics \
    I=input.bam \
    O=gc_bias_metrics.txt \
    CHART=gc_bias_chart.pdf \
    S=gc_summary.txt \
    R=reference.fa
```

### deepTools computeGCBias

```bash
computeGCBias \
    -b input.bam \
    --effectiveGenomeSize 2913022398 \
    -g hg38.2bit \
    -o gc_bias.txt \
    --biasPlot gc_bias.pdf
```

### Interpret GC Bias

| Issue | Symptom |
|-------|---------|
| Under-representation | Low GC coverage drops |
| Over-representation | High GC coverage elevated |
| PCR bias | Strong correlation |

## Strand Balance

A balanced 0.48-0.52 forward/reverse ratio applies to WGS / WES / generic DNA-seq on autosomes. Expected to deviate for: stranded RNA-seq (deliberately strand-asymmetric -- verify with RSeQC `infer_experiment.py`), bisulfite (CT vs GA), small-RNA / strand-specific RNA-seq, and chrY/chrM regions. Per-chromosome strand imbalance > 5% on autosomes suggests aligner artifacts; on chrX/chrY suggests sex-mismatch.

### Calculate Strand Ratio

```bash
forward=$(samtools view -c -F 16 input.bam)
reverse=$(samtools view -c -f 16 input.bam)
echo "Forward: $forward"
echo "Reverse: $reverse"
ratio=$(echo "scale=4; $forward / $reverse" | bc)
echo "F/R ratio: $ratio"
```

### Check Strand Bias per Chromosome

```bash
for chr in chr1 chr2 chr3; do
    fwd=$(samtools view -c -F 16 input.bam $chr)
    rev=$(samtools view -c -f 16 input.bam $chr)
    echo "$chr: F=$fwd R=$rev ratio=$(echo "scale=2; $fwd/$rev" | bc)"
done
```

## Mapping Quality Distribution

### Extract MAPQ Distribution

```bash
samtools view input.bam | cut -f5 | sort -n | uniq -c | sort -k2 -n
```

### Calculate Mean MAPQ

```bash
samtools view input.bam | awk '{sum+=$5; count++} END {print "Mean MAPQ:", sum/count}'
```

### MAPQ Distribution Is Bimodal and Aligner-Specific

Mean MAPQ is misleading; distributions are bimodal (0 and aligner-max). For aligner-specific scales and "unique mapping" sentinels, see sam-bam-basics. The fraction of primary mapped reads at MAPQ >= 30 is a more informative summary than the mean.

## Chromosome Coverage Balance

### Calculate Per-Chromosome Coverage

```bash
samtools idxstats input.bam | awk '{print $1, $3/$2}' | head -25
```

### Check for Aneuploidy / Sex Chromosome Imbalance

Median-normalized per-autosome coverage (1.0 = expected diploid; 0.5 = monosomy/sex; 1.5 = trisomy):
```bash
samtools idxstats in.bam | awk '$2>0 && $1!~/^chr[XYM]|^GL|^KI|^chrUn|^chrEBV/ {
    cov[$1] = $3 / $2
}
END {
    n = asort(cov, sorted)
    med = sorted[int(n/2)+1]
    for (c in cov) printf "%s\t%.3f\n", c, cov[c]/med
}'
```

For full ancestry / contamination / relatedness checking, use `verifybamid2`, `somalier`, or `peddy` -- they account for population AFs, not just per-contig depth.

## Mismatch Rate

### Picard CollectAlignmentSummaryMetrics

```bash
java -jar picard.jar CollectAlignmentSummaryMetrics \
    I=input.bam \
    R=reference.fa \
    O=alignment_summary.txt
```

### Key Metrics

| Metric | Description | Good Value |
|--------|-------------|------------|
| PCT_PF_READS_ALIGNED | Mapped % | > 95% |
| PF_MISMATCH_RATE | Mismatches | < 1% |
| PF_INDEL_RATE | Indels | < 0.1% |
| STRAND_BALANCE | Strand ratio | ~0.5 |

## Comprehensive Validation Script

**Goal:** Run all key alignment QC checks in a single pass and generate a summary report.

**Approach:** Combine samtools flagstat, stats, idxstats, and strand counts into one script that outputs pass/warn/fail calls.

```bash
#!/bin/bash
BAM=$1
REF=$2
NAME=$(basename $BAM .bam)
OUTDIR=${3:-qc}

mkdir -p $OUTDIR

echo "=== Alignment Validation: $NAME ===" | tee $OUTDIR/report.txt

echo -e "\n--- Flagstat ---" | tee -a $OUTDIR/report.txt
samtools flagstat $BAM | tee -a $OUTDIR/report.txt

echo -e "\n--- Mapping Rate ---" | tee -a $OUTDIR/report.txt
mapped=$(samtools view -c -F 4 $BAM)
total=$(samtools view -c $BAM)
rate=$(echo "scale=2; $mapped / $total * 100" | bc)
echo "Mapping rate: ${rate}%" | tee -a $OUTDIR/report.txt

echo -e "\n--- Proper Pairing ---" | tee -a $OUTDIR/report.txt
proper=$(samtools view -c -f 2 $BAM)
pair_rate=$(echo "scale=2; $proper / $mapped * 100" | bc)
echo "Proper pairing: ${pair_rate}%" | tee -a $OUTDIR/report.txt

echo -e "\n--- Insert Size ---" | tee -a $OUTDIR/report.txt
samtools stats $BAM | grep "insert size average" | tee -a $OUTDIR/report.txt

echo -e "\n--- Strand Balance ---" | tee -a $OUTDIR/report.txt
fwd=$(samtools view -c -F 16 $BAM)
rev=$(samtools view -c -f 16 $BAM)
strand_ratio=$(echo "scale=3; $fwd / $rev" | bc)
echo "Forward: $fwd, Reverse: $rev, Ratio: $strand_ratio" | tee -a $OUTDIR/report.txt

echo -e "\n--- Chromosome Coverage ---" | tee -a $OUTDIR/report.txt
samtools idxstats $BAM | head -25 | tee -a $OUTDIR/report.txt

echo -e "\nReport: $OUTDIR/report.txt"
```

## Python Validation Module

A skeleton; full implementation is in `examples/validate_alignment.py`:
```python
import pysam

class AlignmentValidator:
    def __init__(self, bam_file):
        self.bam = pysam.AlignmentFile(bam_file, 'rb')

    def report(self, sample_size=100000):
        # Sample reads, compute mapping rate, proper-pair rate, MAPQ dist, strand balance
        # See examples/validate_alignment.py for full implementation
        ...
```

The first-N-reads sampling pattern is biased toward chr1 (different GC content and complexity than chrM/chrX/chrY/alt contigs). For unbiased per-chromosome statistics, use `samtools view -s 42.01 input.bam` (seed-prefixed fraction `INT.FRAC` is reproducible; bare `0.01` is not) instead of head-of-file iteration.

## Quality Thresholds Summary

| Metric | Good | Warning | Fail |
|--------|------|---------|------|
| Mapping rate | > 95% | 90-95% | < 90% |
| Proper pairing | > 90% | 80-90% | < 80% |
| Duplicate rate (assay-specific) | see bam-statistics decision table | -- | -- |
| Strand balance | 0.48-0.52 | 0.45-0.55 | Outside |
| Mean MAPQ | > 40 | 30-40 | < 30 |
| GC bias | < 1.2x | 1.2-1.5x | > 1.5x |

## Related Skills

- bam-statistics - Per-assay metric thresholds, depth/coverage tools, mosdepth
- alignment-filtering - Aligner-specific MAPQ thresholds (canonical home)
- duplicate-handling - Library-aware dedup decisions
- sam-bam-basics - MAPQ-by-aligner table
- chip-seq/chipseq-qc - ChIP-specific QC (FRiP, NSC, RSC)
