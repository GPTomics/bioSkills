---
name: bio-splicing-qc
description: Assesses RNA-seq data quality specifically for alternative splicing analysis. QC layers include experimental design audit (library prep, read length, depth, replicates), STAR 2-pass cohort-style alignment, junction saturation curves and discovery plateau detection, novel-vs-known junction ratio diagnostics, junction-overhang distribution, splice-site strength scoring (MaxEntScan intrinsic + SpliceAI context-aware), strandedness verification, GENCODE basic vs comprehensive choice, and rRNA contamination screening. Splicing analysis is more demanding than DGE on read length, depth, library prep, alignment strategy, and annotation choice — failures silently bias PSI estimates and inflate novel-junction false positives. Use when evaluating data suitability for splicing analysis, troubleshooting low event detection, or designing sequencing experiments where AS is a primary endpoint.
tool_type: python
primary_tool: RSeQC
---

## Version Compatibility

Reference examples tested with: RSeQC 5.0+, STAR 2.7.11+, samtools 1.19+, pysam 0.22+, regtools 1.0+, maxentpy 0.0.1+, spliceai 1.3+, matplotlib 3.8+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Splicing-Specific Quality Control

Splicing analysis is more demanding than DGE on read length, depth, library prep, alignment strategy, and annotation choice. Failures in any of these silently bias PSI estimates and inflate novel-junction false positives. The decision sequence is: experimental design → library prep → alignment strategy → annotation → diagnostic metrics. Each layer's failure mode is distinct.

## QC Layer Taxonomy

| Layer | Target | Tool | Fails when |
|-------|--------|------|------------|
| Experimental design | Read length, depth, replicates, library type | Pre-sequencing review | <PE 75nt; n<3 vs n<3; <30M reads/sample |
| Library prep | poly(A) vs rRNA depletion | Pre-sequencing review | poly(A) library used for IR analysis |
| Alignment | STAR 2-pass cohort-style | STAR | 1-pass loses 14% novel junctions; per-sample 2-pass introduces inconsistency |
| Junction discovery | Saturation, novelty | RSeQC `junction_saturation`, `junction_annotation` | Curve still rising = under-sequenced; novel% >40% suggests biology or artifact |
| Strand specificity | Library protocol consistency | RSeQC `infer_experiment` | Wrong `--libType` halves usable junctions |
| Splice site strength | Cryptic vs canonical | MaxEntScan, SpliceAI | Weak splice sites (MaxEnt<5) may indicate cryptic, regulated, or annotation error |
| Junction overhang | Read-junction support quality | pysam CIGAR parsing | Overhang <8nt = high false-positive rate |
| Contamination | rRNA, adapters | fastq_screen | >20% rRNA in "depleted" library = failed depletion |
| Annotation | GENCODE basic vs comprehensive | Annotation choice | Basic for canonical events; comprehensive for DTU |

## Decision Tree by Question

| Question | Recommended QC |
|----------|-----------------|
| Will my planned RNA-seq design support AS analysis? | Pre-sequencing audit: library type, read length, depth, replicates |
| Is my data suitable for cassette exon analysis? | Junction saturation + known/novel ratio + read length |
| Why does my AS analysis call so few events? | Saturation curve, depth, library type, alignment 2-pass |
| Why does my AS analysis call so many novel junctions? | Annotation completeness + novel% + biology check (TDP-43, SF3B1) |
| Are my SpliceAI predictions calibrated for my tissue? | MaxEntScan + SpliceAI concordance for known sites |
| Did STAR 2-pass actually run cohort-style? | Verify SJ.out.tab merging across samples |
| Is intron retention detectable in my data? | Library type (must be rRNA-depleted); strand-specific |
| Are my microexons detectable? | Read length >=100; aligner anchor settings; consider VAST-TOOLS |

## Experimental Design Audit (Before Sequencing)

| Decision | For splicing analysis | Rationale |
|----------|------------------------|-----------|
| **Library prep** | rRNA depletion (Ribo-Zero, RiboCop) | poly(A) selection loses pre-mRNA, nascent transcripts, and detained introns; for IR analysis rRNA depletion is mandatory (Sims 2014 *Genome Res*) |
| **Read length** | PE 100-150 nt (PE 150 strongly preferred) | Junction-spanning reads need >=8 nt overhang on each exon; 50 nt SE biases toward shorter exons (Wang 2008 *Nature*) |
| **Pairing** | Paired-end | Single-end loses fragment-level disambiguation of junctions |
| **Depth** | 50-100M reads/sample | DGE-grade 30M misses low-PSI events; 100M for low-abundance event discovery |
| **Strandedness** | Stranded library (Illumina TruSeq stranded) | Distinguishes overlapping antisense; some tools double-count unstranded junctions |
| **Replicates** | n>=3 per condition | n=2 vs n=2 has poor calibration in most tools (especially SUPPA2) |
| **Annotation** | GENCODE basic for canonical, comprehensive for DTU/discovery | basic = high-confidence; comprehensive includes putative — affects FDR control |
| **Microexons** | PE 100+ with `--alignSJoverhangMin 8`; VAST-TOOLS | Default aligners miss 3-27nt exons |
| **Long-intron genes (TTN, brain)** | Increased `--alignIntronMax` | Default 1Mb may miss >1Mb introns |

## STAR 2-Pass Alignment

**Goal:** Maximize novel-junction sensitivity for downstream AS analysis.

**Approach:** Run STAR once per sample to discover novel junctions (pass 1), merge novel junctions across cohort, then re-align with the augmented junction set (pass 2). Cohort-style 2-pass beats per-sample basic 2-pass for differential splicing because all samples use the same junction reference.

```bash
# Pass 1: per-sample
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix pass1_${sample}_ \
    --outSJtype Standard \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1
```

```bash
# Cohort-style 2-pass: collect all SJ.out.tab from pass 1
cat pass1_*_SJ.out.tab | awk '$5 > 0 && $7 >= 3' | sort -u > cohort_novel_SJ.tab

# Pass 2: re-align with augmented junctions
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbFileChrStartEnd cohort_novel_SJ.tab \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix pass2_${sample}_ \
    --outSJtype Standard \
    --twopassMode None \
    --quantMode GeneCounts \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 3
```

| Approach | Novel-junction recovery | Cohort consistency |
|----------|-------------------------|--------------------|
| 1-pass with annotation | ~80-86% (depends on GENCODE completeness) | High (annotation-based) |
| Per-sample basic 2-pass (`--twopassMode Basic`) | >=94% | Variable (each sample has its own junction set) |
| Cohort-style 2-pass (manual merge) | >=94% | High (shared junction reference) |

Per-sample 2-pass (`--twopassMode Basic`) is simpler but produces inconsistent junction sets across samples; for differential splicing the **cohort-style** version is preferred (Veeneman 2016 *Bioinformatics*).

The pass-1 filter `awk '$5 > 0 && $7 >= 3'` keeps junctions with strand info AND >=3 unique reads — adjust threshold to balance discovery vs noise.

## Junction Saturation

**Goal:** Determine whether sequencing depth is sufficient for comprehensive splicing detection.

**Approach:** Run RSeQC junction saturation; check whether the discovery curve plateaus.

```bash
junction_saturation.py \
    -i sample.bam \
    -r gencode_v45.bed \
    -o sample_junc_sat \
    -m 50 \
    -v 100000
```

```python
import subprocess
import pandas as pd

samples = ['s1.bam', 's2.bam', 's3.bam']
for sample in samples:
    subprocess.run([
        'junction_saturation.py',
        '-i', sample,
        '-r', 'gencode_v45.bed',
        '-o', sample.replace('.bam', '_junc_sat')
    ], check=True)
```

The output `*.junctionSaturation_plot.r` plots known + novel junctions vs subsampled reads.

**Plateau detection rule:** if from 80% to 100% of reads, the junction count rises by <2%, consider it plateaued. Still rising means more sequencing would yield more junctions.

For AS analysis, **plateau on the known junction curve** is the requirement; novel-junction curves often don't plateau even at deep coverage (which is biologically informative — novel junctions are inherently rarer events).

## Novel-vs-Known Junction Ratio

**Goal:** Detect annotation/mapping issues or biologically interesting cryptic splicing.

**Approach:** Classify junctions with RSeQC and compute the novel:known ratio.

```bash
junction_annotation.py -i sample.bam -r gencode_v45.bed -o sample_junc_annot
```

```python
import pandas as pd

junc = pd.read_csv('sample_junc_annot.junction.xls', sep='\t')
total = junc['total_splicing_events'].sum()

by_class = junc.groupby('annotation')['total_splicing_events'].sum()
known_frac = by_class.get('known', 0) / total
novel_frac = (by_class.get('partial_novel', 0) + by_class.get('novel', 0)) / total

print(f'known: {known_frac:.1%}, novel: {novel_frac:.1%}')
```

| Known fraction | Status | Interpretation |
|----------------|--------|----------------|
| >=80% | Healthy | Comprehensive annotation, good alignment |
| 60-80% | Acceptable | Check annotation completeness or organism |
| <60% | Suspect or interesting | Mapping artifacts, contamination, OR biologically informative |

**High novel-junction rate may be biology, not artifact:**
- **TDP-43 loss** (ALS/FTD post-mortem brain): cryptic exon de-repression in UNC13A, STMN2, ATG4B (Brown 2022 *Nature*; Klim 2019 *Nat Neurosci*)
- **SF3B1-mutant** cancer (MDS, CLL, uveal melanoma): cryptic 3'ss ~10-30nt upstream of canonical (Darman 2015 *Cell Rep*)
- **Non-model organism**: GENCODE-grade annotation unavailable; novel junctions reflect annotation gaps not biology
- **Microbial / viral contamination**: reads aligning to host but with unusual junctions

If novel% >40%, drill down: check organism, check spliceosomal mutation status, check known disease signatures.

## Junction Read Overhang and Coverage

**Goal:** Profile per-junction read counts and overhang distribution to identify weakly-supported events.

**Approach:** Parse CIGAR for N (intron) operations; tally per-junction reads and minimum exon overhangs.

```python
import pysam
from collections import defaultdict

def junction_stats(bam_path):
    bam = pysam.AlignmentFile(bam_path, 'rb')
    counts = defaultdict(int)
    min_overhang = defaultdict(lambda: float('inf'))

    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary:
            continue
        ref_pos = read.reference_start
        cumulative_query = 0
        cigar = read.cigartuples
        for i, (op, length) in enumerate(cigar):
            if op == 3:
                left_match = sum(l for o, l in cigar[:i] if o in (0, 7, 8))
                right_match = sum(l for o, l in cigar[i+1:] if o in (0, 7, 8))
                overhang = min(left_match, right_match)
                key = (read.reference_name, ref_pos, ref_pos + length)
                counts[key] += 1
                min_overhang[key] = min(min_overhang[key], overhang)
            if op in (0, 2, 3, 7, 8):
                ref_pos += length

    bam.close()
    return counts, dict(min_overhang)

counts, overhang = junction_stats('sample.bam')
print(f'total junctions: {len(counts)}')
print(f'>= 10 reads: {sum(1 for c in counts.values() if c >= 10)}')
print(f'overhang >= 8 nt: {sum(1 for k, c in counts.items() if overhang[k] >= 8)}')
```

Junction reads with overhang <8 nt are common false positives, especially for novel sites. Most callers default to >=8 nt anchor for this reason. Microexon-aware aligners use overhang as low as 6 nt with explicit configuration.

## Splice Site Strength (MaxEntScan and SpliceAI)

**Goal:** Score donor and acceptor splice sites to flag weak / cryptic sites and to predict variant impact on splicing.

**Approach:** Use MaxEntScan (sequence information content) and SpliceAI (context-aware deep-learning) — they answer different questions.

```python
from maxentpy.maxent import score5, score3

donor = 'CAGGTAAGT'
acceptor = 'TTTTTTTTTTTTTTTTTTTTCAG'
print(f"5'ss MaxEnt: {score5(donor):.2f}")
print(f"3'ss MaxEnt: {score3(acceptor):.2f}")
```

| Score | Interpretation | Source |
|-------|----------------|--------|
| 5'ss MaxEnt > 8 | Strong donor | Yeo & Burge 2004 *J Comput Biol* |
| 5'ss MaxEnt 5-8 | Moderate | |
| 5'ss MaxEnt < 5 | Weak / cryptic | |
| 3'ss MaxEnt > 8 | Strong acceptor | |
| 3'ss MaxEnt < 5 | Weak / cryptic | |
| SpliceAI delta > 0.2 | PP3 supporting (ClinGen SVI 2023) | Walker 2023 *AJHG* |
| SpliceAI delta > 0.5 | PP3 moderate | |
| SpliceAI delta > 0.8 | PP3 strong | |

**MaxEntScan vs SpliceAI:**
- **MaxEntScan** scores sequence information content (intrinsic strength). Captures position-wise dependencies at the consensus.
- **SpliceAI** predicts in-vivo usage probability given full pre-mRNA context (10 kb window).
- A position with **high MaxEnt but low SpliceAI** is intrinsically strong but contextually silenced (chromatin, trans factors).
- A position with **low MaxEnt but high SpliceAI** is intrinsically weak but contextually used (enhancer-driven, e.g. weak donors stabilized by ESEs).
- Report both for variant interpretation; for variant impact see `splice-variant-prediction`.

## Picard CollectRnaSeqMetrics and Gene-Body Coverage

**Goal:** Get integrated RNA-seq QC including intronic / exonic / intergenic mapping rates and gene-body coverage uniformity.

**Approach:** Run picard CollectRnaSeqMetrics for mapping distribution; RSeQC `geneBody_coverage.py` for 5'-3' bias.

```bash
picard CollectRnaSeqMetrics \
    I=sample.bam \
    O=sample.rna_metrics.txt \
    REF_FLAT=refFlat.txt \
    STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
    RIBOSOMAL_INTERVALS=rRNA_intervals.interval_list

# Strandedness conversion (foot-gun):
# Reverse-stranded (Illumina TruSeq Stranded; NEB Ultra II Directional — both dUTP):
#   rMATS  --libType fr-firststrand
#   featureCounts -s 2
#   Picard STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND
# Forward-stranded (Lexogen QuantSeq FWD, certain ligation-based kits):
#   rMATS  --libType fr-secondstrand
#   featureCounts -s 1
#   Picard STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND
# STAR has no library-strand flag; pass --outSAMstrandField intronMotif
# (works for any library) so downstream tools can read XS tags.

geneBody_coverage.py \
    -i sample.bam \
    -r gencode_v45.bed \
    -o sample_geneBody
```

| Metric | Healthy | Concerning |
|--------|---------|------------|
| PCT_CODING_BASES | >=50% | <30% (suggests degradation or mis-priming) |
| PCT_UTR_BASES | 20-40% | >>50% (3' bias) |
| PCT_INTRONIC_BASES | <30% (poly(A)); <60% (rRNA-depleted) | >50% (poly(A)) suggests pre-mRNA contamination |
| PCT_INTERGENIC_BASES | <10% | >20% (genomic DNA contamination) |
| MEDIAN_5PRIME_TO_3PRIME_BIAS | 0.7-1.3 | >2 or <0.5 (severe degradation) |
| Gene body coverage curve | Flat | Strong 3' skew = RIN low or library mis-prep |

3' bias (degraded RNA) directly reduces splicing-event detection because junction reads scatter across the gene body; with 3' bias they concentrate near the 3' end and miss CDS junctions.

## Strandedness Verification

```bash
infer_experiment.py -i sample.bam -r gencode_v45.bed -s 200000
```

Output reports the fraction of reads consistent with each library type:

| Output pattern | Library type | rMATS `--libType` |
|----------------|---------------|---------------------|
| ~50% / ~50% | Unstranded | `fr-unstranded` |
| >=90% "++ , --" | Forward-stranded | `fr-secondstrand` |
| >=90% "+- , -+" | Reverse-stranded (Illumina TruSeq stranded) | `fr-firststrand` |

**Wrong strand setting halves usable junction reads** — always verify before quantification. RSeQC `infer_experiment.py` is fast and authoritative.

## Annotation Choice

| GENCODE level | Contents | Use for |
|---------------|----------|---------|
| Basic | High-confidence canonical isoforms | Standard rMATS, leafcutter, SUPPA2 |
| Comprehensive | All transcripts including putative/predicted | DTU pipelines (DRIMSeq+DEXSeq, satuRn), isoform discovery |
| RefSeq | NCBI curated | Less complete than GENCODE; legacy use |
| Ensembl | Same content as GENCODE in vertebrates | Different attribute conventions |

Comprehensive captures more biology but inflates DTU multiple-testing burden and includes annotation noise. For event-level (rMATS) AS, basic is usually adequate; for transcript-level DTU (DRIMSeq, satuRn), comprehensive may be necessary to capture rare isoforms.

## rRNA Contamination Check

```bash
fastq_screen --conf fastq_screen.conf --threads 8 sample_R1.fq.gz
```

Or post-alignment:

```bash
samtools view -c sample.bam | awk '{print "total:",$0}'
samtools view -c -L rRNA_intervals.bed sample.bam | awk '{print "rRNA:",$0}'
```

| rRNA fraction | Library type | Status |
|----------------|---------------|--------|
| >=20% | "depleted" | Failed depletion; redo |
| 5-20% | "depleted" | Acceptable; some rRNA leakage |
| <5% | poly(A) | Healthy |
| <5% | "depleted" | Excellent depletion |
| 1-3% | poly(A) | Suggests RNA degradation |

>5% rRNA in a poly(A) library suggests degraded RNA; >20% in a "depleted" library indicates failed depletion.

## Per-Tool Failure Modes

### RSeQC `junction_saturation`: Subsampling Behavior

**Trigger:** Running on extremely deep BAM (>200M reads).

**Mechanism:** RSeQC subsamples at 5%, 10%, ..., 100%; with very deep BAMs, the early subsamples are still tens of millions of reads, masking saturation behavior.

**Symptom:** Curve appears flat throughout; uninformative.

**Fix:** Subsample BAM with `samtools view -s 0.1` before running junction_saturation; or use `-s` flag to set custom step intervals.

### STAR 2-Pass: Per-Sample Inconsistency

**Trigger:** Using `--twopassMode Basic` on differential splicing cohorts.

**Mechanism:** Per-sample 2-pass means each sample has its own SJ.out.tab; samples may differ in which novel junctions they re-align against.

**Symptom:** Inconsistent novel junction calls across replicates; rMATS `--novelSS` differential calls don't replicate.

**Fix:** Switch to cohort-style 2-pass (collect all pass-1 SJ.out.tabs, merge, re-align all samples with merged set).

### MaxEntScan: Out-of-Range Sequences

**Trigger:** Sequences with N bases or wrong length.

**Mechanism:** `score5` expects exactly 9 nt (3 exon + 6 intron); `score3` expects 23 nt (20 intron + 3 exon).

**Symptom:** ValueError or silently incorrect score.

**Fix:** Pre-validate sequence length and N-content; use a wrapper that returns NaN for invalid inputs.

### SpliceAI: TensorFlow Memory

**Trigger:** Running spliceai on large VCF without GPU.

**Mechanism:** TensorFlow CPU mode is slow; default batch size may exceed memory.

**Symptom:** OOM kill; very slow runtime (hours per chromosome).

**Fix:** Use `-D 50` for screening (fastest); split VCF by chromosome; use GPU when available.

### `infer_experiment.py`: Sample Size

**Trigger:** Running on very low-coverage region or small subsample (-s).

**Mechanism:** Default sample size is 200,000 reads; with low coverage, this isn't met.

**Symptom:** "0 of 200000 reads" output; cannot infer strand.

**Fix:** Lower `-s` to actual available reads; or use `-q 30` to filter by quality.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `STAR: SJDBoverhang differs from genome` | Index built with different overhang than current run | Rebuild index with `--sjdbOverhang` matching read length - 1 |
| `RSeQC: BED format error` | Annotation BED has wrong column order | Convert with `awk` or `bedtools` |
| `MaxEntScan: invalid sequence character N` | N in input | Filter or replace; document |
| `samtools view: missing index` | BAM not indexed | `samtools index sample.bam` |
| `STAR: too many SJs in cohort merge` | Cohort SJ.out.tab too large after merge | Filter to junctions in >=3 samples or with >=3 unique reads |
| `regtools: invalid CIGAR` | Non-spec read in BAM | Filter with `samtools view -h -F 0x100 -F 0x800` |

## Quality Thresholds

| Metric | Good | Acceptable | Poor | Source |
|--------|------|------------|------|--------|
| Read length (PE) | 150 nt | 100 nt | <75 nt | Wang 2008 *Nature* |
| Sequencing depth | >=100M | 50-100M | <30M | DGE-grade insufficient |
| Junction saturation | Plateau (<2% growth in last 20%) | Near plateau | Still rising | RSeQC convention |
| Known-junction fraction | >=80% | 60-80% | <60% (suspect or interesting) | RSeQC convention |
| Junctions >=10 reads | >=50% | 30-50% | <30% | rMATS reliability cutoff |
| 5'ss / 3'ss MaxEnt | >8 | 5-8 | <5 | Yeo & Burge 2004 |
| Strandedness | >90% one direction | 70-90% | <70% | RSeQC convention |
| rRNA in depleted library | <5% | 5-20% | >20% | Sims 2014 *Genome Res* |
| 2-pass STAR | Cohort-style | Per-sample basic | 1-pass only | Veeneman 2016 *Bioinformatics* |

## Troubleshooting Low Event Detection

| Issue | Possible causes | Solutions |
|-------|-----------------|-----------|
| Few events called | Low depth; short reads; SE; wrong strand | Increase depth; use PE150; verify libType |
| High novel junctions | Annotation gaps; mapping artifacts; biology (TDP-43, SF3B1) | Update annotation; check 2-pass; consider biology |
| Low IR detection | poly(A) library | Use rRNA depletion |
| Microexons missing | Default aligner anchors too long | VAST-TOOLS, MicroExonator, or long-read |
| Many weak splice sites | Cryptic splicing | Validate with MaxEnt + SpliceAI; consider RNA-seq from secondary tissue |
| FDR uncalibrated at low n | n=2 vs n=2 | Use leafcutter or Shiba; avoid SUPPA2 alone |
| PSI variance high across replicates | Library prep / RIN inconsistency | Check RIN; consider RNA degradation |
| Sashimi plot mismatch with PSI | Junction-imbalance bias in rMATS | Run Shiba; or filter by overhang distribution |

## Common Pitfalls

- **Skipping STAR 2-pass** — loses ~14% of novel junctions; matters for any non-canonical organism or condition.
- **Per-sample 2-pass instead of cohort-style** — produces inconsistent junction sets; differential splicing calls don't replicate.
- **poly(A) library for IR analysis** — biases toward mature transcripts; depletes pre-mRNA / nascent / detained intron signal.
- **PE 50nt single-end** — junction-spanning reads need >=8nt overhang on both sides; biases toward shorter exons.
- **Wrong `--libType`** — halves usable junctions; always verify with `infer_experiment.py`.
- **Using basic GENCODE for DTU** — basic excludes putative/rare isoforms; DTU pipelines may underdetect.
- **Using MaxEntScan alone for variant interpretation** — misses context-dependent regulation; pair with SpliceAI.
- **Treating high novel% as artifact reflexively** — could be biology (TDP-43, SF3B1, non-model organism); investigate.

## Related Skills

- splicing-quantification - PSI estimation after QC passes
- read-alignment/star-alignment - STAR 2-pass detail and parameter tuning
- read-qc/quality-reports - General sequencing QC (FastQC, MultiQC)
- read-qc/contamination-screening - rRNA / adapter / cross-species contamination
- splice-variant-prediction - SpliceAI / Pangolin for variant impact
- long-read-splicing - When short-read QC is fundamentally limiting (microexons, complex isoforms)
- differential-splicing - Downstream tool that requires QC pass

## References

- Yeo & Burge 2004 *J Comput Biol* - MaxEntScan
- Jaganathan et al 2019 *Cell* - SpliceAI
- Walker et al 2023 *Am J Hum Genet* - ClinGen SVI splicing thresholds
- Veeneman et al 2016 *Bioinformatics* - STAR 2-pass benchmark
- Wang et al 2008 *Nature* - PE 75+ recommendation for AS
- Sims et al 2014 *Genome Res* - rRNA depletion vs poly(A) for AS
- Brown et al 2022 *Nature* - cryptic exons in TDP-43 loss
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing in ALS
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
- Wang et al 2024 *Nat Protoc* - rMATS-turbo
- Dobin et al 2013 *Bioinformatics* - STAR aligner
