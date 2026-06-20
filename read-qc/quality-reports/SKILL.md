---
name: bio-read-qc-quality-reports
description: Generates and interprets per-file and cross-sample QC reports from FASTQ data with FastQC, falco, and MultiQC, covering Phred quality, per-base composition, GC, duplication, overrepresented sequences, and adapter content. Use when performing initial QC on raw sequencing reads, validating preprocessing, or judging a multi-sample cohort for outliers and batch effects. For long reads use NanoPlot; for adapter/quality remediation route to adapter-trimming, quality-filtering, or fastp-workflow.
tool_type: cli
primary_tool: fastqc
---

## Version Compatibility

Reference examples tested with: FastQC 0.12+, MultiQC 1.21+, falco 1.2+, seqkit 2.5+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Quality Reports -- the traffic light is a hypothesis about WGS DNA, not a verdict

Generate per-file QC with FastQC/falco and aggregate the cohort with MultiQC, then READ THE PLOTS against the assay rather than trusting pass/warn/fail.

**"Run quality control on FASTQ files"** -> Compute per-base quality, composition, GC, duplication, and adapter profiles per file, then aggregate across samples to find outliers.
- CLI: `fastqc -t 8 *.fastq.gz` then `multiqc .`
- Long reads: `NanoPlot --fastq reads.fastq.gz` (FastQC assumes fixed-length short reads)

Scope: this skill OWNS raw-FASTQ QC reporting and interpretation, and carries the cross-cutting quality-score / chemistry / duplication concepts the rest of read-qc depends on. Remediation lives elsewhere -> read-qc/adapter-trimming, read-qc/quality-filtering, read-qc/fastp-workflow. Contamination -> read-qc/contamination-screening. Transcriptome QC on the BAM -> read-qc/rnaseq-qc. OUT OF SCOPE: any modification of the reads.

## The Single Most Important Modern Insight

1. **FastQC pass/warn/fail are heuristics calibrated to random whole-genome DNA, so they FALSE-FAIL on every other assay.** RNA-seq fails per-base content (random-hexamer bias) and duplication (high-expression molecules); amplicon fails duplication and GC by design; bisulfite fails base content (C->T conversion); small-RNA fails length distribution; single-cell R1 fails everything (it is barcode+UMI, not biology). A red light is a HYPOTHESIS about a WGS library. For any other protocol, first ask "is this module expected to deviate for this chemistry?" before treating the red as a defect. Read the plot shape; the traffic light is calibration noise.

2. **On 2-color chemistry (NextSeq, NovaSeq, MiniSeq) G is the ABSENCE of signal, so poly-G tails are called at HIGH quality and the quality plot will NOT flag them.** When a cluster runs out of template, dark cycles read as a run of Gs with high confidence. Quality trimming alone does not remove them. They surface as a 3'-end RISE in G content (per-base sequence content) and a spurious high-GC spike, and they mis-map or manufacture false somatic variants if left in. The fix is a chemistry-aware poly-G trim (fastp auto-enables it from the instrument ID; cutadapt `--nextseq-trim`), not a quality cutoff. Read the per-base CONTENT plot on any 2-color run, not just the quality plot.

3. **The duplication percentage is read-level and complexity-blind: it cannot tell a PCR jackpot from genuine high abundance.** Identical reads from a highly expressed transcript, a targeted amplicon, or a ChIP/ATAC peak are counted as duplicates even though they are independent biological molecules. Duplication % is a function of BOTH library complexity AND sequencing depth (a good library sequenced deeply shows high duplication). It is a PROMPT to reason about library complexity (preseq), never an automatic "remove duplicates" -- and removing duplicates in non-UMI RNA-seq is actively wrong (read-qc/umi-processing, read-qc/rnaseq-qc).

Bonus trap: NovaSeq/NextSeq emit BINNED quality scores (RTA3 uses four values: 2, 12, 23, 37), so FastQC box plots look blocky/quantized. This is the instrument's quality table, NOT bad data and NOT something to fix. The bin edges are RTA-version-specific (NovaSeq X / RTA4 differs) -- never hard-code one bin set.

## Tool Taxonomy

| Tool | Role | Mechanism / when |
|------|------|------------------|
| FastQC | Per-file short-read QC (HTML + zip) | Java; the module set below; duplication/overrep from the first 100k distinct reads. The de-facto standard per-file report. |
| falco | Drop-in FastQC re-implementation (C++) | ~3x faster, lower memory, same module names and MultiQC-compatible output. Use when FastQC throughput bottlenecks a large cohort. |
| MultiQC | Cross-sample aggregator (SCRAPER, not a re-analyzer) | Walks directories, regex-matches each tool's log/report, parses the numbers, builds one cohort report. The unit of review for multi-sample studies. |
| seqkit stats | Instant tabular FASTA/FASTQ numbers | `seqkit stats -a`: N50, Q20%, Q30%, GC%, length quartiles. For quick numbers and assembly/long-read contexts where FastQC is the wrong shape. |
| NanoPlot / NanoComp | Long-read (ONT/PacBio) QC | Read-length and quality distributions, yield, N50, length-vs-quality. The correct first pass for long reads; FastQC's fixed-length assumptions break there. |

## Decision Tree by Scenario

| Scenario | Use | Why |
|----------|-----|-----|
| Per-file Illumina short-read QC | FastQC (or falco) | Module-level diagnostics; read the plots by assay |
| Many samples / a study cohort | FastQC/falco then MultiQC | Outlier and batch detection is RELATIVE; only visible overlaid |
| Long reads (ONT/PacBio) | NanoPlot / NanoComp | FastQC is built for fixed-length short reads |
| Instant numbers, assembly input | seqkit stats -a | N50/Q20/Q30/GC in one line; no HTML overhead |
| Large cohort, FastQC too slow | falco then MultiQC | Same output, ~3x faster |

Default when uncertain: FastQC on each file, then MultiQC over the run directory, and judge each sample against the cohort.

## FastQC Modules -- thresholds and the expert read

Thresholds are FastQC's `limits.txt` defaults (calibrated to random WGS DNA). The expert read is what to conclude BEYOND the traffic light.

| Module | Default warn / fail | Expert read |
|--------|--------------------|-------------|
| Per base sequence quality | warn LQ<10 or median<25; fail LQ<5 or median<20 | 3' decay is normal; blocky boxes on NovaSeq are binning; this plot will NOT reveal poly-G on 2-color |
| Per tile sequence quality | spatial deviation (no numeric) | A hot tile band across cycles = a localized flowcell problem (bubble, debris, edge); reason no MultiQC table replaces raw FastQC |
| Per sequence quality scores | distribution of per-read mean Q | A low-Q hump = a junk subpopulation to FILTER (not trim) |
| Per base sequence content | warn dev>10%; fail dev>20% | First ~12 bp skew = random-hexamer priming (Hansen 2010), expected for RNA-seq, do NOT trim it. A 3'-end skew is poly-G / adapter -- act on that |
| Per sequence GC content | warn dev>15%; fail dev>30% | SHAPE matters: bimodal/secondary peak = contamination; sharp spike = adapter dimer / overrepresented; a shifted single peak = wrong-GC reference assumption |
| Per base N content | warn N>5%; fail N>20% | Ns at a fixed position = a failed cycle; rising 3' Ns = dying clusters |
| Sequence length distribution | warn if lengths differ; fail if any length 0 | WARNs trivially after trimming and on long reads -- ignore for those |
| Sequence duplication levels | warn if <70% would remain; fail if <50% | Read-level, complexity-blind (see insight 3); high = think complexity, not dedup |
| Overrepresented sequences | warn >0.1%; fail >1% | Most diagnostic module: it prints the sequence -- BLAST it (adapter dimer, rRNA, primer, poly-G) |
| Adapter content | warn k-mer>5%; fail >10% | A curve climbing toward 3' = read-through from short inserts; this panel IS an insert-size readout (route to adapter-trimming) |
| K-mer content | (deprecated, off by default) | Only appears in old reports; do not build guidance on it |

Algorithm note (why duplication/overrep are estimates): FastQC tracks only the first 100,000 DISTINCT sequences, keys on the first 50 bp for reads >75 bp (so 3' errors do not fragment a duplicate family), counts by exact identity, and extrapolates the "% remaining if deduplicated" headline. It is a sample-based estimate, not a full-library dedup.

## Duplication Taxonomy -- four causes, four actions

| Class | Mechanism | Detected by | Action |
|-------|-----------|-------------|--------|
| Optical | One real cluster mis-segmented (non-patterned flowcell) | Same tile, pixel distance (Picard default 100) | Removable; spatially local artifact |
| ExAmp / patterned | One molecule seeds two nanowells (HiSeq X/4000, NovaSeq) | Spatially clustered, larger radius (Picard 2500 for patterned) | Removable; the reason patterned flowcells need the bigger pixel distance |
| PCR | Same fragment amplified and sequenced twice | Identical 5' coordinates post-alignment (+UMI if present) | Mark/remove for variant calling; NEVER coordinate-dedup amplicon (use UMIs) |
| Natural / biological | Independent identical molecules (high coverage, expressed genes, amplicon start) | Indistinguishable from PCR at read level without UMIs | KEEP -- removing biases quantification (do not dedup non-UMI RNA-seq) |

The read-level duplication % FastQC reports cannot separate these. Use preseq (Daley & Smith 2013) to model the complexity curve and ask "how many NEW molecules would more sequencing buy?" -- that curve, not a single %, judges whether a library is exhausted or just deeply sequenced.

## Quality scores and encoding

Phred Q = -10*log10(P_error): Q20 = 1% error, Q30 = 0.1%, Q40 = 0.01%. Q30 is the routine Illumina target; bulk Q40+ is uncommon on legacy chemistry (phasing, signal decay) and is a tell for re-binned or synthetic data on old runs, though XLEAP-SBS (NovaSeq X, NextSeq 2000) genuinely reaches Q40+. Modern data is universally Phred+33; any Phred+64 file (Illumina 1.3-1.7) feeds 31-too-high scores to a +33-assuming tool and passes garbage silently -- convert it (`seqtk seq -Q64 -V`). A quality byte below ASCII 64 (digits/punctuation) proves +33; detection tools sample reads to break ties.

## MultiQC -- the cohort is the unit of review

MultiQC does NOT re-analyze data; it scrapes tool logs/reports (`search_patterns.yaml`), parses the numbers, and tabulates them per sample. Consequences: it is only as good as the files left on disk and the sample-name parsing (name collisions merge samples -- check `multiqc_data/multiqc_sources.txt`), and it reports whatever the upstream tool wrote (a wrong reference or wrong strandedness shows as a coherent-but-wrong table, not an error). Read the General Statistics table FIRST -- outliers jump out as a column anomaly -- then overlay per-base-quality / GC / duplication and ask whether the low-quality set maps to one lane / prep batch / operator. The batch effect caught here at QC is the one not chased for a month in the DE results.

```bash
# Per-file QC, then aggregate the run
fastqc -t 8 -o qc/raw/ raw_data/*.fastq.gz
multiqc qc/raw/ -o qc/multiqc/ -f

# Compare before vs after trimming in one report
fastqc -t 8 -o qc/trimmed/ trimmed/*.fastq.gz
multiqc qc/ -o qc/compare/ -f          # picks up both raw/ and trimmed/

# Long reads do not go through FastQC
NanoPlot --fastq ont_reads.fastq.gz -o qc/nanoplot/
```

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Every RNA-seq sample fails per-base content | Random-hexamer 5' bias (Hansen 2010) | Expected; do not trim the first bases |
| High-Q reads but a 3' G-content rise on NovaSeq | 2-color poly-G (dark cycles = G) | Chemistry-aware poly-G trim (fastp / cutadapt --nextseq-trim), not -q |
| FastQC quality boxes look quantized/blocky | NovaSeq/NextSeq binned qualities (RTA3) | Expected; not a defect, do not "fix" |
| MultiQC merges two samples into one row | Over-aggressive name cleaning / collision | Check multiqc_sources.txt; use `--fn_as_s_name` or fix names |
| Duplication 60%, urge to dedup RNA-seq | Read-level dup is complexity-blind | Do not dedup non-UMI RNA-seq; assess complexity (preseq) |
| FastQC crashes / huge plot on long reads | Fixed-length short-read assumptions | Use NanoPlot / seqkit stats instead |
| FastQC module missing in MultiQC | The fastqc_data.txt was not on disk / wrong dir | Point MultiQC at the directory holding the zip/data files |

## References

de Sena Brandine G, Smith AD. 2019. Falco: high-speed FastQC emulation for quality control of sequencing data. F1000Research 8:1874.
Ewing B, Hillier L, Wendl MC, Green P. 1998. Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Research 8(3):175-185.
Ewing B, Green P. 1998. Base-calling of automated sequencer traces using phred. II. Error probabilities. Genome Research 8(3):186-194.
Hansen KD, Brenner SE, Dudoit S. 2010. Biases in Illumina transcriptome sequencing caused by random hexamer priming. Nucleic Acids Research 38(12):e131.
Daley T, Smith AD. 2013. Predicting the molecular complexity of sequencing libraries. Nature Methods 10(4):325-327.
Ewels P, Magnusson M, Lundin S, Kaller M. 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32(19):3047-3048.
Shen W, Le S, Li Y, Hu F. 2016. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE 11(10):e0163962.
De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. 2018. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics 34(15):2666-2669.

## Related Skills

read-qc/adapter-trimming - Remove read-through adapter flagged by the adapter-content panel
read-qc/quality-filtering - Drop low-quality reads and trim ends
read-qc/fastp-workflow - All-in-one QC + trim, including 2-color poly-G
read-qc/contamination-screening - Resolve a bimodal-GC or unexpected overrepresented-sequence signal
read-qc/rnaseq-qc - Transcriptome QC (strandedness, gene-body coverage) on the aligned BAM
sequence-io/sequence-statistics - Programmatic per-file sequence summaries
