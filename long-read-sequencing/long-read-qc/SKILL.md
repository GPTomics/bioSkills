---
name: bio-long-read-sequencing-long-read-qc
description: Assesses Oxford Nanopore and PacBio long-read quality with NanoPlot, cramino, NanoComp, pycoQC/toulligQC, and seqkit, and filters reads with chopper/Filtlong for the downstream goal. Covers why read-only Qscore is an uncalibrated posterior (real accuracy needs a reference BAM), why the sequencing_summary.txt is required for run-health metrics, intent-conditioned filtering (preserve long reads and small replicons for assembly, filter almost nothing for variant calling), the chimera/internal-adapter trap that fabricates SVs, and PacBio rq-based HiFi QC. Use when judging a long-read run, computing read N50 or percent identity, filtering reads before assembly or variant calling, comparing barcodes/runs, or reading run-health red flags.
tool_type: cli
primary_tool: nanoplot
---

## Version Compatibility

Reference examples tested with: NanoPlot 1.42+ (NanoPack2), cramino 0.14+, chopper 0.7+, Filtlong 0.2+, seqkit 2.5+, pycoQC 2.5+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags (chopper/cramino are fast-moving Rust tools)

Inputs that determine what QC is even possible - record them:
- `sequencing_summary.txt` is produced by the basecaller (Dorado/Guppy), not the FASTQ. pycoQC/toulligQC REQUIRE it for pore activity, yield-over-time, and translocation speed. FASTQ-only hand-off permanently loses the run-health layer.
- Percent identity requires a reference BAM (NanoPlot `--bam` / cramino); it cannot come from FASTQ.

If code throws an error, introspect the installed tool (`NanoPlot --help`, `cramino --help`) and adapt the example to the actual API rather than retrying.

# Long-Read QC

**"Is my long-read run any good?"** -> Read length N50 and yield from FASTQ, real percent identity from a reference BAM, run-health from the sequencing_summary, then filter for the downstream goal.
- CLI: `NanoPlot --fastq reads.fq.gz -o qc/` (overview), `cramino aln.bam` (fast BAM stats + identity), `pycoQC -f sequencing_summary.txt -o run.html` (run health)

## The Single Most Important Modern Insight -- Read-Only Qscore Is a Self-Graded Posterior; Real Accuracy and the Failures That Sink a Run Are Only Visible Against a BAM and the Summary

Three corrections a naive long-read QC misses:

1. **Per-read Qscore is an uncalibrated basecaller posterior, not an empirical error rate.** It is the Phred of the mean per-base error probability (NOT the arithmetic mean of Q values), assigned by the basecaller to its own output. ONT's own data: bases labeled Q20 are empirically ~Q12.5 on older chemistries; R10 sup and HiFi are better calibrated but read-only Q still overstates accuracy. Real accuracy is gap-compressed identity from a reference BAM (cramino, NanoPlot `--bam`). Treat Q thresholds as relative knobs, not accuracy guarantees.
2. **The sequencing_summary.txt is the run-health layer, and it is not in the FASTQ.** Pore/channel activity, yield-over-time, translocation speed, and barcode breakdown come from the basecaller's summary TSV. Hand a collaborator only FASTQ and that layer is gone (re-basecalling from POD5 can regenerate it; FASTQ cannot).
3. **The right filter depends on intent, not a fixed cutoff.** Assembly wants the long reads (which are the lowest-Q) and small replicons preserved - subsample by quality, never hard-length-cut. Variant calling wants depth - filter almost nothing and let the caller model per-base Q. HiFi is already Q20+ - do not Phred-filter it like noisy CLR.

## Tool Roles

| Tool | Input | Reports |
|------|-------|---------|
| NanoPlot | FASTQ / BAM / summary | length dist, length-vs-quality, yield; `--bam` adds percent identity |
| cramino | BAM/CRAM | fast N50, yield, gap-compressed identity, `--phased` block N50, `--karyotype` |
| NanoComp | multiple FASTQ/BAM/summaries | compare runs/barcodes (length, quality, identity) |
| pycoQC / toulligQC | sequencing_summary.txt | run health: pore activity, mux map, yield/speed over time, barcodes |
| seqkit stats -a | FASTA/FASTQ | N50, quartiles, total bases, GC |
| chopper | FASTQ (stdin) | filter/trim by mean Q and length |
| Filtlong | FASTQ | keep best reads by length x identity; subsample to a target depth |

Read N50 = the length where 50% of total bases are in reads at least that long (length-weighted, far above the median); it predicts assembly contiguity. NanoFilt and the rrwick Porechop are deprecated/unmaintained (use chopper and Porechop_ABI).

## Intent-Conditioned Filtering Decision Tree

| Goal | Filter | Why |
|------|--------|-----|
| Bacterial / small-genome assembly | light Q/length, then subsample by quality to ~50-100x (`filtlong --target_bases`) | a hard 10 kb length cut erases small plasmids; quality-subsampling beats length filtering |
| Eukaryotic / large-genome assembly | minimal; keep the long tail | the longest (lowest-Q) reads span repeats; over-filtering loses N50 |
| SV calling | light Q only; trim chimeras | chimeras fabricate SVs; trimming matters more than Q filtering |
| SNV / small-variant calling | almost nothing (`chopper -q 10`) | callers model per-base Q and want depth |
| PacBio HiFi | `rq >= 0.99` only | already Q20+; Phred filtering adds nothing |
| cDNA / direct RNA | orient/trim (pychopper), no hard length cut | transcript length is biology; a length cut biases the expression matrix |

## Core Commands

```bash
# Overview from FASTQ (length + posterior quality only - not real accuracy)
NanoPlot --fastq reads.fq.gz -o qc_fastq/ --N50
seqkit stats -a reads.fq.gz                     # N50 + quartiles, fast

# Real accuracy: fast BAM stats incl. gap-compressed identity (needs a reference BAM)
cramino aln.bam
NanoPlot --bam aln.bam -o qc_bam/               # percent identity scatter

# Run health (requires the basecaller's summary)
pycoQC -f sequencing_summary.txt -o run_qc.html

# Compare barcodes / runs
NanoComp --bam s1.bam s2.bam s3.bam --names s1 s2 s3 -o compare/

# Filter for VARIANT calling: light quality only
chopper -q 10 -i reads.fq.gz | gzip > q10.fq.gz

# Subsample for ASSEMBLY: by quality to ~100x of a 5 Mb genome (never a hard length cut)
filtlong --target_bases 500000000 reads.fq.gz | gzip > subsampled.fq.gz
```

## Per-Method Failure Modes

### Trusting FASTQ Qscore as accuracy
**Trigger:** judging a run from `NanoStat --fastq` mean Q. **Mechanism:** Q is an uncalibrated posterior. **Symptom:** "Q20 reads" that are ~94% accurate. **Fix:** align and read gap-compressed identity (cramino / NanoPlot `--bam`).

### QC without the summary
**Trigger:** only FASTQ/BAM at hand-off. **Mechanism:** run-health metrics live in sequencing_summary.txt. **Symptom:** cannot see pore death, mux map, or yield-over-time. **Fix:** obtain the summary (or re-basecall from POD5 to regenerate it).

### Over-filtering erases assembly value
**Trigger:** a blunt `-q 15` or hard 10 kb length cut before assembly. **Mechanism:** the longest reads are the lowest-Q; small plasmids fall under a length floor. **Symptom:** worse N50; missing plasmids. **Fix:** subsample by quality (Filtlong `--target_bases`), keep the long tail, never length-floor above the smallest replicon.

### Chimeras masquerade as SVs
**Trigger:** undetected internal adapters (two molecules ligated as one read). **Mechanism:** the read's halves map to different loci. **Symptom:** phantom translocations/insertions in the SV VCF. **Fix:** check whether Dorado already trimmed/split; use Porechop_ABI for unknown adapters; suspect a biologically implausible long-read spike.

### Re-filtering HiFi like CLR
**Trigger:** Phred-quality-filtering PacBio HiFi. **Mechanism:** HiFi is Q20+ consensus already. **Symptom:** wasted reads, no accuracy gain. **Fix:** filter on `rq >= 0.99` only.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Q20-labeled bases ~Q12.5 empirically | ONT EPI2ME | read-only Q overstates accuracy; verify by alignment |
| Subsample assembly data to ~50-100x | Wick 2026 | >100x slows assemblers and can propagate systematic errors |
| Pore occupancy <~70% in hour 1 rarely recovers | ONT guidance | run-health red flag for early pore death |
| Translocation ~400 b/s (R10 DNA) | ONT chemistry | drift off target correlates with falling basecall Q |
| HiFi `rq >= 0.99` (Q20); `>= 0.999` for Q30 | PacBio CCS | the canonical HiFi accuracy filter |
| `-q 10` as a light QC floor | convention | a relative knob, not a 90%-accuracy guarantee |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| NanoPlot gives no percent identity | run on FASTQ | use `--bam` (identity needs alignment) |
| pycoQC errors / empty | no sequencing_summary.txt | supply the basecaller summary |
| cramino fails on FASTQ | cramino is BAM/CRAM only | give it the aligned BAM |
| Assembly N50 dropped after filtering | hard length/quality cut removed long reads | subsample by quality instead |
| Missing small plasmids | length floor above the replicon size | lower/remove the length floor |
| Phantom SVs in the VCF | chimeric reads | trim/split internal adapters |

## References

- De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. 2018. NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics* 34(15):2666-2669.
- De Coster W, Rademakers R. 2023. NanoPack2: population-scale evaluation of long-read sequencing data (cramino, chopper). *Bioinformatics* 39(5):btad311.
- Leger A, Leonardi T. 2019. pycoQC, interactive quality control for Oxford Nanopore Sequencing. *J Open Source Softw* 4(34):1236.
- Steinig E, Coin L. 2022. Nanoq: ultra-fast quality control for nanopore reads. *J Open Source Softw* 7(69):2991.
- Bonenfant Q, Noé L, Touzet H. 2023. Porechop_ABI: discovering unknown adapters in Oxford Nanopore sequencing reads. *Bioinform Adv* 3(1):vbac085.
- Shen W, Le S, Li Y, Hu F. 2016. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLoS ONE* 11(10):e0163962.

## Related Skills

- basecalling - Produces the reads and the sequencing_summary.txt this QC needs
- long-read-alignment - Produces the BAM required for real percent identity
- structural-variants - Chimeras flagged here fabricate SVs there
- medaka-polishing - QC/subsample reads before polishing
- genome-assembly/long-read-assembly - Subsample by quality before assembling
- genome-assembly/genome-profiling - K-mer ploidy/size estimate alongside read QC
- read-qc/quality-reports - General (short-read-oriented) read QC
- sequence-io/sequence-statistics - FASTA/FASTQ summary statistics
