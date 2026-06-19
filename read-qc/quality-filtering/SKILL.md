---
name: bio-read-qc-quality-filtering
description: Filters reads by quality, length, N content, and complexity with Trimmomatic, fastp, and Cutadapt, including sliding-window trimming, per-read unqualified-base filtering, and 2-color poly-G removal. Use when reads have poor-quality tails, when an assembly or k-mer workflow needs clean input, or when a junk read subpopulation must be dropped. For adapter removal use adapter-trimming; for all-in-one preprocessing use fastp-workflow.
tool_type: cli
primary_tool: trimmomatic
---

## Version Compatibility

Reference examples tested with: Trimmomatic 0.39+, fastp 0.23+, Cutadapt 4.4+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Quality Filtering -- trim lightly or not at all, and never without a length filter

Trim low-quality bases and drop low-quality reads with Trimmomatic (sliding window / MAXINFO), fastp (per-read filter + window cut), or Cutadapt (BWA-style quality trim).

**"Filter reads by quality"** -> Remove low-quality bases and/or discard reads below quality/length thresholds.
- CLI: `fastp -i in.fq -o out.fq --cut_right -q 20 -l 36` (window trim + per-read filter + length gate)
- CLI: `trimmomatic SE in.fq out.fq SLIDINGWINDOW:4:20 MINLEN:36`

Scope: this skill OWNS quality/length/N/complexity filtering. Adapter removal -> read-qc/adapter-trimming. Single-pass trim+QC -> read-qc/fastp-workflow. Reading the quality plots -> read-qc/quality-reports. OUT OF SCOPE: contamination removal (read-qc/contamination-screening).

## The Single Most Important Modern Insight

1. **Modern local aligners SOFT-CLIP low-quality tails, so quality trimming is usually unnecessary -- and AGGRESSIVE quality trimming actively harms downstream results.** Williams 2016 showed aggressive trimming changed expression estimates for >10% of genes; Del Fabbro 2013 showed stringent Q>30 DEGRADES de novo assembly; MacManes 2014 found gentle trimming (remove only Phred<2-5) optimal for RNA-seq; GATK discourages quality trimming because BQSR recalibrates qualities itself. Trim ADAPTER always (read-qc/adapter-trimming); quality-trim lightly or not at all before a soft-clipping aligner. The workflows that genuinely need quality trimming are assembly, k-mer/pseudo-alignment, small-RNA, amplicon, and variant calling WITHOUT BQSR.

2. **Quality FILTERING (drop whole reads) and quality TRIMMING (cut bases within a read) are different operations with different tools.** fastp's `-q/-u/-n` filters whole reads by the fraction of unqualified bases; `--cut_right` / Trimmomatic `SLIDINGWINDOW` trims bases from a window scan. Filtering removes a junk subpopulation (a low-Q hump in the per-sequence-quality plot); trimming shortens reads with decayed tails. Choose by whether the problem is some bad reads or bad ends.

3. **A short post-trim read mis-maps, so quality trimming MUST be paired with a minimum-length filter.** Williams 2016 showed that adding a post-trim min-length filter mitigates most of the expression distortion that trimming introduces, because over-trimmed fragments that would map spuriously are dropped instead. `MINLEN` (Trimmomatic, always last), `-l` (fastp), `-m` (cutadapt) are not optional add-ons; they are the safety mechanism that makes trimming safe.

Two-color note: on NextSeq/NovaSeq the quality scores are binned to four values (RTA3: 2, 12, 23, 37), so a sliding-window threshold like 4:15 partitions between the 12 and 23 bins rather than acting on a smooth gradient -- thresholds tuned on HiSeq-era 0-40 qualities behave differently. And poly-G tails are HIGH quality, so a quality filter does not remove them (use poly-G trimming).

## Tool Taxonomy

| Tool | Mechanism | When it wins |
|------|-----------|--------------|
| fastp | Per-read unqualified-base filter (`-q/-u/-n`) + window cut (`--cut_right`) + auto poly-G | DEFAULT; one fast pass, filtering and trimming together |
| Trimmomatic | `SLIDINGWINDOW` / `MAXINFO` window trim; ordered step pipeline; orphan handling | Legacy/reproducibility pipelines; MAXINFO length-vs-quality balance |
| Cutadapt | `-q` BWA running-sum quality trim (combined with adapter removal) | When already running cutadapt for adapters; precise per-end control |

## Decision Tree by Scenario

| Workflow | Quality trimming | Why |
|----------|------------------|-----|
| Alignment-based DNA/RNA (BWA-MEM, STAR, Bowtie2 local, HISAT2) | Light or none | Aligner soft-clips tails; aggressive trim distorts expression |
| GATK variant calling with BQSR | None | BQSR recalibrates; trimming interferes |
| De novo assembly | Moderate (~Q20) + min-length | Low-Q errors corrupt the de Bruijn graph; stringent Q>30 over-trims |
| k-mer / pseudo-alignment (kallisto/salmon) | Light + adapter | Errors create phantom k-mers |
| A junk read subpopulation (bimodal per-seq quality) | FILTER whole reads (`-e`/AVGQUAL) | Trimming cannot fix a globally bad read |
| Variant calling WITHOUT BQSR | Moderate + min-length | No recalibration safety net |

Default when uncertain: trim adapter, apply a light window trim plus a minimum-length filter, then confirm with FastQC.

## Trimmomatic

Steps run in COMMAND-LINE ORDER; put quality steps before MINLEN so the length check reflects all trimming.

```bash
# Single-end: light leading/trailing + window, length-gated
trimmomatic SE -phred33 in.fq.gz out.fq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

# Paired-end (four outputs: paired + orphan)
trimmomatic PE -phred33 -threads 8 \
    R1.fq.gz R2.fq.gz \
    R1_paired.fq.gz R1_unpaired.fq.gz R2_paired.fq.gz R2_unpaired.fq.gz \
    SLIDINGWINDOW:4:20 MINLEN:36

# MAXINFO: adaptive length-vs-quality balance (strictness <0.2 favors length, >0.8 favors correctness)
trimmomatic SE in.fq.gz out.fq.gz MAXINFO:40:0.5 MINLEN:36
```

| Step | Meaning |
|------|---------|
| SLIDINGWINDOW:W:Q | scan 5'->3'; cut from the point where the W-bp window mean drops below Q |
| MAXINFO:L:S | adaptive trim balancing target length L against error rate; strictness S in 0-1 |
| LEADING:Q / TRAILING:Q | cut 5'/3' bases below Q (also removes N) |
| MINLEN:L / AVGQUAL:Q | DROP read if shorter than L / if mean quality below Q |
| CROP:L / HEADCROP:N | cap length / remove first N bases (do NOT HEADCROP random-hexamer bias -- see below) |

Do NOT HEADCROP the first ~12 bp of RNA-seq to "fix" the wavy per-base-content plot: that pattern is random-hexamer priming bias (Hansen 2010), not adapter, and trimming it just discards real data without removing the underlying bias.

## fastp

fastp separates per-read FILTERING from window TRIMMING. Quality filtering is on by default (-q 15).

```bash
# Per-read quality filter: base < Q20 is 'unqualified'; drop read if >40% unqualified or >5 Ns
fastp -i in.fq.gz -o out.fq.gz -q 20 -u 40 -n 5 -l 36

# Window trim from the 3' (Trimmomatic SLIDINGWINDOW analogue) + length gate
fastp -i R1.fq.gz -I R2.fq.gz -o R1.fq.gz -O R2.fq.gz \
      --cut_right --cut_window_size 4 --cut_mean_quality 20 -l 36

# Drop globally low-quality reads by mean quality (filter, not trim)
fastp -i in.fq.gz -o out.fq.gz -e 25

# 2-color poly-G (auto-enabled for NextSeq/NovaSeq from the instrument ID)
fastp -i in.fq.gz -o out.fq.gz --trim_poly_g

# Low-complexity filter (e.g. poly-A / homopolymer-rich reads)
fastp -i in.fq.gz -o out.fq.gz --low_complexity_filter --complexity_threshold 30
```

fastp flags: `-q` qualified quality (default 15), `-u` unqualified percent limit (default 40), `-n` N base limit (default 5), `-e` average-quality filter (default 0 = off), `-l` length required (default 15), `--length_limit` max length (long form only), `--cut_front/--cut_tail/--cut_right` window cut modes (off by default), `--cut_window_size` (4), `--cut_mean_quality` (Q20).

## Cutadapt

`-q` uses the BWA running-partial-sum algorithm, not a fixed cutoff, so a single high-Q base inside a low-Q run does not stop trimming. Quality trimming runs BEFORE adapter removal.

```bash
# 3'-only quality trim with a length gate (5',3' form: -q 15,20)
cutadapt -q 20 -m 36 -o out.fq.gz in.fq.gz

# Combined adapter + light quality trim, paired
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -m 36 \
         -o R1.fq.gz -p R2.fq.gz R1.fq.gz R2.fq.gz
```

## Quantitative Thresholds

| Parameter | Typical | Rationale |
|-----------|---------|-----------|
| Window quality | Q20 (4:20) | 1% error; light. Aggressive (Q25-30) distorts expression/assembly (Williams 2016, Del Fabbro 2013) |
| fastp -q / -u | Q15 / 40% | fastp defaults; a base under Q15 is unqualified, read dropped if >40% unqualified |
| MINLEN / -l / -m | 36 (150 bp reads) | Mandatory after trimming; short reads mis-map. Scale up for longer inserts |
| complexity_threshold | 30 (30%) | fastp default for low-complexity filtering |
| MAXINFO strictness | 0.2-0.8 | <0.2 favors length, >0.8 favors correctness |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Expression estimates shift for many genes | Aggressive quality trimming | Trim lightly; always add a min-length filter (Williams 2016) |
| Variant calling worse after trimming | Trimmed before/around BQSR | Do not quality-trim for GATK BQSR workflows |
| Window threshold behaves oddly on NovaSeq | Binned quality (4 values) makes windows coarse | Expect step-like behavior; do not port HiSeq thresholds blindly |
| Reads mis-map after trimming | No min-length filter, over-trimmed fragments | Add MINLEN / -l / -m |
| Poly-G tails survive quality filtering | Poly-G is high quality on 2-color | Use `--trim_poly_g` / cutadapt `--nextseq-trim` |
| R1/R2 out of sync | Independent SE trimming of mates | Use Trimmomatic paired outputs or fastp/cutadapt paired mode |

## References

Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30(15):2114-2120.
Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17):i884-i890.
MacManes MD. 2014. On the optimal trimming of high-throughput mRNA sequence data. Frontiers in Genetics 5:13.
Del Fabbro C, Scalabrin S, Morgante M, Giorgi FM. 2013. An extensive evaluation of read trimming effects on Illumina NGS data analysis. PLoS ONE 8(12):e85024.
Williams CR, Baccarella A, Parrish JZ, Kim CC. 2016. Trimming of sequence reads alters RNA-Seq gene expression estimates. BMC Bioinformatics 17:103.
Hansen KD, Brenner SE, Dudoit S. 2010. Biases in Illumina transcriptome sequencing caused by random hexamer priming. Nucleic Acids Research 38(12):e131.

## Related Skills

read-qc/adapter-trimming - Remove adapter before quality filtering
read-qc/quality-reports - Read the quality plots that motivate filtering
read-qc/fastp-workflow - All-in-one preprocessing in a single pass
read-alignment/bwa-alignment - Soft-clipping aligner that absorbs low-quality tails
read-alignment/star-alignment - Soft-clipping RNA aligner (light trimming preferred)
