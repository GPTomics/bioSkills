# BWA-MEM2 Alignment - Usage Guide

## Overview

bwa-mem2 is the maintained, architecture-aware reimplementation of BWA-MEM: near-identical output, ~1.5-3x faster, ~2x the memory, and a different index format. It is the standard DNA short-read aligner for WGS, WES, and germline/somatic variant-calling pipelines. The aligner itself is the easy part of the job; this skill emphasizes the decisions that actually determine whether downstream variant calling works -- read groups, the reference analysis set (decoy/ALT), the -M/-Y output flags, and the strict collate/fixmate/sort/markdup ordering.

## Prerequisites

```bash
conda install -c bioconda bwa-mem2 bwa samtools
```

- A reference FASTA. For human variant calling, choose the analysis set (decoy, optionally ALT; see the SKILL's Build Index) before indexing.
- For ancient-DNA work, the original `bwa` binary (the `aln`/`samse`/`sampe` path), not just bwa-mem2.

## Quick Start

Tell your AI agent what you want to do:
- "Align my paired-end WGS reads to the human reference with read groups for GATK"
- "Build a bwa-mem2 index for my reference genome"
- "Align and mark duplicates in one streaming pipeline"
- "Map reads for structural-variant calling with soft-clipped supplementary alignments"
- "Run bwa-mem2 with reproducible output across thread counts"

## Example Prompts

### Basic alignment
> "Align reads_R1.fq.gz and reads_R2.fq.gz to reference.fa with bwa-mem2, add read groups for sample NA12878, and output a coordinate-sorted, indexed BAM."

### WGS / variant-calling pipeline
> "Run a complete WGS alignment with duplicate marking for sample NA12878, then tell me which QC numbers to check before calling variants."

### Structural variants
> "Map my WGS reads for split-read SV detection -- which flags keep the supplementary alignments that Manta and Delly need?"

### Reproducibility and reference
> "Align with reproducible results across different thread counts, and confirm I am mapping to a decoy-containing GRCh38 analysis set."

### Ancient DNA
> "I have short, damaged ancient-DNA reads with a low mapping rate using bwa mem -- what should I use instead and why?"

## What the Agent Will Do

1. Confirm the reference analysis set (decoy/ALT) is appropriate for the downstream caller and index it with `bwa-mem2 index` if needed.
2. Align with read groups (SM/ID/PL/LB) injected at mapping time, streaming to a coordinate-sorted BAM to avoid a large intermediate SAM.
3. For SV pipelines, add `-Y` (and never `-M`) so supplementary split reads keep their full sequence.
4. Mark duplicates in the strict order collate -> fixmate -m -> sort -> markdup, skipping dedup entirely for amplicon/PCR data.
5. Add `-K 100000000` when bit-stable, thread-count-invariant output is required.
6. Run the QC gate (alignment rate, properly-paired, duplicate/complexity) before handing the BAM to a variant caller, routing the stats and their interpretation to alignment-files/bam-statistics.

## Read Group Fields

| Field | Description | Example |
|-------|-------------|---------|
| ID | unique read-group id; the BQSR error-model unit | flowcell.lane |
| SM | sample name; variant callers group by SM | NA12878 |
| PL | platform; affects error modeling | ILLUMINA |
| LB | library; MarkDuplicates dedups WITHIN a library | lib1 |
| PU | platform unit (optional but recommended) | flowcell.lane.barcode |

## Tips

- Always inject read groups at mapping time; a BAM without SM/LB is rejected or mis-merged by GATK, and fixing it later is a full rewrite.
- Pipe straight to `samtools sort` to avoid writing a huge intermediate SAM; never leave an uncompressed SAM on disk for a cohort.
- Use `-Y` for any pipeline that calls structural variants, and never `-M` there -- `-M` hides the split-read evidence.
- Use `-K 100000000` for reproducible output; without it, per-batch insert-size estimation makes multithreaded runs vary.
- Never run MarkDuplicates on amplicon/multiplex-PCR data -- identical primer-defined ends are by design, so dedup deletes real coverage; use UMIs instead.
- A low mapping rate is a reference/input problem (wrong build, wrong species, un-trimmed adapter, contamination), not something more sequencing fixes; check the reference and trim first.
- For allele-specific work, expect reference bias at heterozygous sites and correct it (WASP) rather than trusting raw allele counts.
- For WGS archival, write CRAM instead of BAM (`samtools view -C -T reference.fa`, ~2x smaller); CRAM is reference-backed and needs the same reference to decode.
- With multiple lanes per library, give each lane its own read group (distinct ID/PU, same SM/LB), merge the per-lane BAMs, then mark duplicates once across the merged library so cross-lane PCR duplicates are caught (-> alignment-files).

## Related Skills

- bowtie2-alignment - ChIP/ATAC DNA mapping with end-to-end vs local modes
- star-alignment - RNA splice-aware alignment (when reads cross junctions)
- read-qc/fastp-workflow - Trim and QC reads before alignment
- alignment-files/duplicate-handling - Mark/remove duplicates; UMI-aware dedup
- alignment-files/sam-bam-basics - SAM flags, CIGAR, the cross-tool MAPQ scale, SA tags
- alignment-files/bam-statistics - flagstat/idxstats/stats QC gate; what a high mapping rate hides
- variant-calling/variant-calling - Call variants from the aligned BAM
- variant-calling/structural-variant-calling - SV calling from split/discordant reads (needs -Y)
