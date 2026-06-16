---
name: bio-read-alignment-bwa-alignment
description: Aligns DNA short reads (paired- or single-end) to a reference genome with bwa-mem2, the maintained successor to BWA-MEM, for WGS/WES and germline/somatic variant-calling pipelines; covers index build, read-group injection, the collate/fixmate/sort/markdup ordering, soft-clipping for SV split reads, ALT/decoy-aware mapping on GRCh38, -K determinism, and streaming straight to a sorted BAM. Use when mapping DNA short reads to a reference for variant calling, coverage, ChIP/ATAC (alongside bowtie2-alignment), or SV detection. RNA-seq spliced alignment is star-alignment/hisat2-alignment; BAM sort/dedup/stats, the QC gate, and the cross-tool MAPQ scale are alignment-files; read trimming is read-qc; counting reads over features is rna-quantification.
tool_type: cli
primary_tool: bwa-mem2
---

## Version Compatibility

Reference examples tested with: bwa-mem2 2.2.1+, bwa 0.7.17+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# BWA-MEM2 Alignment -- Read Groups, the Reference, and the Output Contract Decide Downstream Truth

**"Align my DNA reads"** -> Map reads to a reference with seed-and-extend, inject read groups, and stream to a sorted BAM -- because the aligner is the easy part: the reference analysis set, the @RG metadata, the -M/-Y output flags, and the dedup ordering are the decisions that silently determine whether variant calling works.
- CLI: `bwa-mem2 mem -t 8 -R '@RG\tID:s1\tSM:s1\tPL:ILLUMINA\tLB:lib1' ref.fa R1.fq.gz R2.fq.gz | samtools sort -@4 -o aligned.bam -`

Scope: DNA short-read mapping with bwa-mem2 (and original BWA-MEM / bwa-aln), the GRCh38 analysis-set choice (below), and the output contract it must satisfy for GATK/DeepVariant/SV callers. BAM sort/dedup/index/stats, the QC gate (flagstat/idxstats interpretation), and contig-naming reconciliation -> alignment-files. The cross-tool MAPQ scale -> alignment-files/sam-bam-basics. Read trimming/QC -> read-qc. Variant calling -> variant-calling. Structural-variant calling consumes the split reads produced here -> variant-calling/structural-variant-calling. OUT OF SCOPE: RNA (use star-alignment/hisat2-alignment), long reads (long-read-sequencing/long-read-alignment), bisulfite (methylation-analysis/bismark-alignment).

## The Single Most Important Modern Insight

1. **Read groups are a hard contract, not metadata decoration.** GATK (HaplotypeCaller/Mutect2/BQSR) and Picard refuse to run or behave wrongly without `@RG`. SM names the sample callers group by; ID is the BQSR error-model unit; PL sets the error model; LB is the unit MarkDuplicates dedups within (two reads at the same coordinate from different libraries are NOT duplicates). Inject them at mapping time with `-R`; adding them later (Picard AddOrReplaceReadGroups) is a full BAM rewrite. The catastrophic error is a clean-looking BAM that GATK rejects or mis-merges because SM/LB are missing.
2. **On GRCh38, ALT/decoy handling is a correctness decision and a blind +ALT mapping is actively worse than no ALT.** ALT contigs are alternate haplotypes of hyperpolymorphic loci (MHC/HLA); a read matches the primary copy AND the ALT copy, becomes a multimapper, and MAPQ collapses to 0 -- so variant callers drop it and HLA variants vanish. bwa-mem rescues this only if the `<idxbase>.alt` file is present (it then scores non-ALT hits against non-ALT hits only); adding ALT contigs to the FASTA WITHOUT the `.alt` imports the ambiguity and discards the fix. Use a decoy (hs38d1) always -- it absorbs reads from sequence missing from the primary assembly that would otherwise mismap and make recurrent false SNPs. The analysis-set layering is below.
3. **The -M/-Y flags and the dedup ordering are downstream contracts that fail silently.** `-M` marks split (chimeric) pieces as secondary (0x100) for legacy Picard -- but modern tools read supplementary (0x800) natively, and `-M` hides the split-read evidence SV callers need, so do NOT use it for SV work; use `-Y` (soft-clip supplementary) so every split piece keeps its full sequence. Duplicate marking has a strict order: `collate (name) -> fixmate -m -> sort (coordinate) -> markdup`; the `-m` adds the mate-score/MC tags markdup requires, fixmate needs mates adjacent, markdup needs coordinate order. Any other order silently produces wrong duplicate flags.

## What MAPQ Means Here (and what it does not)

bwa-mem MAPQ runs 0-60 and is an ordinal confidence rank derived from the gap between the best and second-best alignment score and the seed coverage -- it is NOT a calibrated `-10*log10 P(wrong)`. It clusters bimodally at 0 (a competing locus scores as well -> multimapper) and 60 (no competitor, good seed support); the middle is sparse. It models only repeat ambiguity, not contamination or heuristic error, so MAPQ 60 means "no competing locus found," not "P(wrong)=1e-6". A `-q 20` pre-filter drops most multimappers; rely on the caller's own MAPQ handling (GATK ignores MAPQ-0 and duplicates) rather than double-filtering aggressively. This scale is bwa-specific -- a `-q 60` "unique" filter empties a Bowtie2 (cap 42/44) or STAR (255=unique) BAM; see alignment-files/sam-bam-basics for the cross-tool scale.

## Tool Taxonomy

| Tool / subcommand | Citation | Mechanism / role | When |
|-------------------|----------|------------------|------|
| `bwa-mem2 mem` | Vasimuddin 2019 IEEE IPDPS 314-324 | architecture-aware reimplementation of BWA-MEM; near-identical output, ~1.5-3x faster, ~2x RAM, different index | the default DNA aligner today |
| `bwa mem` | Li 2013 arXiv:1303.3997 | SMEM seed + chain + banded affine SW; the reference implementation defining "correct" | when bwa-mem2 is unavailable or to match a legacy pipeline |
| `bwa aln` + `samse`/`sampe` | Li & Durbin 2009 Bioinformatics 25:1754 | bounded backtracking, no exact-seed requirement | ancient DNA / very short (<70 bp) damaged reads where the 19-mer seed fails |
| `bwa bwasw` | Li & Durbin 2010 Bioinformatics 26:589 | SW over BWT for long/divergent queries | legacy, superseded by mem; rarely used |
| DRAGEN / Parabricks fq2bam | -- | hardware-accelerated BWA-MEM (FPGA/GPU); concordant, not bit-identical | high-throughput production where a validated pipeline accepts the small delta |
| STAR / HISAT2 | -- | splice-aware (route OUT) | any RNA library -> star-alignment, hisat2-alignment |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Human WGS/WES germline or somatic | `bwa-mem2 mem` + GRCh38 + decoy (+ALT/postalt if HLA matters); add `-Y` | seed-and-extend standard; GATK/DeepVariant expect its soft-clipped, ALT-aware output |
| SV / split-read detection | `bwa-mem2 mem -Y` (no `-M`) | supplementary alignments with full soft-clipped sequence are the split-read signal |
| Reproducible / functional-equivalence pipeline | add `-K 100000000` (and `-Y`) | pins per-batch insert-size estimation so output is thread-count-invariant |
| ChIP-seq / CUT&RUN | bwa-mem2 or bowtie2 | both work; bowtie2 is the ENCODE peak-assay default -> bowtie2-alignment |
| ATAC-seq | bowtie2 (local/dovetail) | adapter read-through + short fragments favor local mode -> bowtie2-alignment |
| Ancient / very short damaged DNA | `bwa aln -l 1024 -n ~0.01-0.03 -o 2` + samse/sampe | mem's mandatory 19-mer seed fails on short damaged reads, biasing toward reference |
| Allele-specific (ASE / allelic binding) | bwa-mem2 then WASP-filter | reads carrying the alt allele align worse -> reference bias at het sites |
| RNA-seq | route OUT to star-alignment / hisat2-alignment | spliced reads need an N-CIGAR-capable aligner |

Default when uncertain: `bwa-mem2 mem` with read groups, streamed to a coordinate-sorted BAM, on a decoy-containing GRCh38 analysis set; add `-Y` for any pipeline that calls SVs.

## Build Index

```bash
bwa-mem2 index reference.fa
# emits reference.fa.0123 .amb .ann .bwt.2bit.64 .pac  (NOT interchangeable with `bwa index`'s .bwt/.sa)
```

For human DNA, obtain the decoy/ALT analysis set rather than a bare GRCh38: bwakit's `run-gen-ref hs38DH` downloads `hs38DH.fa` (primary + ALT + hs38d1 decoy + HLA) AND the `hs38DH.fa.alt` lift file, which must sit next to the index basename for ALT-aware mapping; the GATK/Broad resource bundle ships an equivalent decoy reference. Pick the layer by what the downstream caller can handle: no-alt + decoy (hs38) when the caller is not ALT-aware (the safe minimum), or full + decoy (hs38DH) with ALT-aware mapping plus `bwa-postalt.js` when MHC/HLA accuracy matters -- mismatching mapper-awareness to the analysis set silently degrades MHC/segdup accuracy. Use GRCh38 over GRCh37, and never mix builds across a cohort.

## Align with Read Groups, Stream to a Sorted BAM

```bash
bwa-mem2 mem -t 8 \
    -R '@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1' \
    reference.fa reads_1.fq.gz reads_2.fq.gz | \
    samtools sort -@ 4 -o aligned.sorted.bam -
samtools index aligned.sorted.bam
# single-end: drop reads_2.fq.gz. The literal \t in -R must survive the shell (single quotes do this).
```

## Mark Duplicates (the strict ordering)

```bash
# collate -> fixmate -m -> sort -> markdup. -m adds the ms/MC tags markdup requires.
bwa-mem2 mem -t 8 -R '@RG\tID:s1\tSM:s1\tPL:ILLUMINA\tLB:lib1' reference.fa R1.fq.gz R2.fq.gz | \
    samtools collate -@ 4 -O -u - | \
    samtools fixmate -m -@ 4 -u - - | \
    samtools sort -@ 4 -u - | \
    samtools markdup -@ 4 - aligned.markdup.bam
samtools index aligned.markdup.bam
# NEVER run this on amplicon/multiplex-PCR data: identical primer-defined ends are by design, so
# markdup deletes almost all real coverage. Use UMIs (fgbio/UMI-tools) for PCR-dup removal there.
```

## SV / Split-Read Mapping

```bash
# -Y soft-clips supplementary alignments so every split piece keeps its full sequence (SA-tag chain).
# Do NOT add -M (it demotes split pieces to secondary and hides the SV evidence).
bwa-mem2 mem -t 8 -Y -R '@RG\tID:s1\tSM:s1\tPL:ILLUMINA' reference.fa R1.fq.gz R2.fq.gz | \
    samtools sort -@ 4 -o aligned.sv.bam -
```

## Reproducible Output Across Thread Counts

```bash
# Per-batch insert-size estimation makes default multithreaded output vary; -K pins the batch size.
bwa-mem2 mem -t 8 -K 100000000 -Y \
    -R '@RG\tID:s1\tSM:s1\tPL:ILLUMINA\tLB:lib1' \
    reference.fa R1.fq.gz R2.fq.gz | samtools sort -@ 4 -o aligned.bam -
```

## Ancient / Very Short Damaged DNA

```bash
# mem's 19-mer seed fails on short deaminated reads; bwa aln backtracks with no exact-seed requirement.
bwa index reference.fa
bwa aln -l 1024 -n 0.02 -o 2 -t 8 reference.fa reads.fq.gz > reads.sai
bwa samse reference.fa reads.sai reads.fq.gz | samtools sort -o ancient.bam -
# -l 1024 disables the seed (seed longer than the read); -n relaxes edit distance; -o 2 allows gaps.
# The -n value is benchmark-dependent (commonly ~0.01-0.03); treat as a tunable starting point.
```

## Key Parameters (verified against bwa/bwa-mem2 source defaults)

| Flag | Default | Effect |
|------|---------|--------|
| -t | 1 | threads |
| -k | 19 | min seed (SMEM) length; lower = more sensitive on short/divergent reads, slower; the reason aDNA defeats mem |
| -r | 1.5 | re-seed a MEM longer than k*1.5; lower = more sensitive, slower |
| -c | 500 | discard a seed occurring > N times (NOT 10000; the old sourceforge page is stale) |
| -A / -B / -O / -E | 1 / 4 / 6 / 1 | match / mismatch / gap-open / gap-extend; lower -B for divergent data, lower -O/-E to permit longer indels |
| -L | 5 | clip penalty; higher discourages soft-clipping (toward end-to-end) |
| -T | 30 | min alignment score to OUTPUT (below it the read is reported unmapped) |
| -M | off | mark split hits secondary (legacy Picard); harmful for SV -- prefer -Y |
| -Y | off | soft-clip supplementary alignments (keep full sequence for SV callers) |
| -K | auto | input bases per batch; fix it (e.g. 100000000) for thread-count-invariant output |
| -R | -- | the @RG header line (literal \t) |

## Per-Method Failure Modes

### Missing read groups
**Trigger:** mapping without `-R '@RG...'`. **Mechanism:** GATK groups reads by SM and models error by ID/LB. **Symptom:** GATK errors ("no read group"), or MarkDuplicates/BQSR misbehave; samples cannot be told apart. **Fix:** inject SM/ID/PL/LB at mapping time; Picard AddOrReplaceReadGroups is a rewrite if missed.

### GRCh38 + ALT without the .alt file
**Trigger:** ALT contigs in the FASTA but no `<idxbase>.alt`. **Mechanism:** ALT copies become ordinary contigs, so every MHC/HLA read multimaps and MAPQ collapses to 0. **Symptom:** a MAPQ-0 spike and missing variants in MHC/HLA. **Fix:** supply the `.alt` (alt-aware mapping) + `bwa-postalt.js`, or use a no-alt + decoy analysis set if not ALT-aware (see Build Index above).

### -M used in an SV pipeline
**Trigger:** `bwa-mem2 mem -M ...` then split-read SV calling. **Mechanism:** `-M` marks split pieces secondary (0x100), which SV callers and MarkDuplicates skip. **Symptom:** SV caller finds little split-read support. **Fix:** drop `-M`, add `-Y`, then dedup, then call -> variant-calling/structural-variant-calling.

### fixmate without -m, or wrong dedup order
**Trigger:** `samtools markdup` on a name-sorted BAM, or `fixmate` without `-m`. **Mechanism:** markdup needs coordinate order plus the ms/MC tags `-m` writes. **Symptom:** silently wrong duplicate flags; over- or under-marking. **Fix:** `collate -> fixmate -m -> sort -> markdup` -> alignment-files/duplicate-handling.

### Index built with the wrong tool
**Trigger:** pointing bwa-mem2 at a `bwa index` directory or vice versa. **Mechanism:** bwa-mem2 uses `.bwt.2bit.64`; bwa uses `.bwt`/`.sa`. **Symptom:** "fail to locate the index" / format error. **Fix:** rebuild with the matching tool (`bwa-mem2 index`).

### Reference bias at heterozygous sites
**Trigger:** allele-specific analysis (ASE, allelic ChIP/ATAC, low-VAF somatic) off a raw BAM. **Mechanism:** the alt-allele read carries an extra mismatch, sometimes failing the 19-mer seed -> reference allele over-counted. **Symptom:** false allelic imbalance toward reference; depressed alt-allele VAF. **Fix:** WASP-filter (re-map allele-swapped reads, keep only if placement is unchanged), N-mask, a personalized reference, or a SNP-graph/pangenome.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| min seed -k 19 | bwamem.c::mem_opt_init | sets the sensitivity floor; the reason short damaged aDNA reads fail mem |
| seed cap -c 500 | bwamem.c::mem_opt_init | above this a seed is uninformative and explodes extension; verified default (not 10000) |
| scoring -A1 -B4 -O6 -E1 | bwamem.c::mem_opt_init | the substitution-rate the scheme tolerates is ~0.75*exp(-log(4)*B/A) |
| min output score -T 30 | bwamem.c::mem_opt_init | reads scoring below 30 are reported unmapped |
| -K 100000000 for reproducibility | CCDG / functional-equivalence pipelines | pins per-batch insert-size estimation independent of thread count |
| bwa-mem2 human index/runtime RAM ~10 GB | bwa-mem2 README (approximate) | ~2x original bwa; plan node memory accordingly |
| GRCh38 + hs38d1 decoy for human WGS | GATK/CCDG guidance | decoy absorbs missing-assembly reads that otherwise make recurrent false variants |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "fail to locate the index" | bwa vs bwa-mem2 index mismatch | rebuild with the matching tool (`bwa-mem2 index`) |
| GATK "no read group" | missing `-R '@RG...'` | inject SM/ID/PL/LB at mapping time |
| Low mapping rate | wrong reference/species, un-trimmed adapter, contamination | confirm species with read-qc/contamination-screening; check the reference build; trim first -> read-qc; interpret stats -> alignment-files/bam-statistics |
| MAPQ-0 spike, missing HLA variants | GRCh38+ALT without `.alt`, or no decoy | use decoy + ALT-aware mapping/postalt, or no-alt+decoy (see Build Index) |
| Wrong/unstable BAM across runs | default per-batch insert estimation with multithreading | add `-K 100000000` and fix input order |
| Almost all reads marked duplicate | MarkDuplicates run on amplicon/PCR data | do not dedup amplicon; use UMIs -> read-qc/umi-processing, alignment-files/duplicate-handling |
| markdup mis-marks duplicates | name-sorted input or fixmate without `-m` | `collate -> fixmate -m -> sort -> markdup` |

## References

- Li H. 2013. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997.
- Li H, Durbin R. 2009. Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics* 25:1754-1760.
- Li H, Durbin R. 2010. Fast and accurate long-read alignment with Burrows-Wheeler transform. *Bioinformatics* 26:589-595.
- Vasimuddin M, Misra S, Li H, Aluru S. 2019. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *IEEE IPDPS* 2019:314-324.
- Li H, Handsaker B, Wysoker A, et al. 2009. The Sequence Alignment/Map format and SAMtools. *Bioinformatics* 25:2078-2079.
- van de Geijn B, McVicker G, Gilad Y, Pritchard JK. 2015. WASP: allele-specific software for robust molecular quantitative trait locus discovery. *Nat Methods* 12:1061-1063.

## Related Skills

- bowtie2-alignment - ChIP/ATAC DNA mapping with end-to-end vs local modes
- star-alignment - RNA splice-aware alignment (when reads cross junctions)
- read-qc/fastp-workflow - Trim and QC reads before alignment
- read-qc/umi-processing - UMI extraction/dedup for amplicon and low-input libraries (do not coordinate-dedup amplicon)
- alignment-files/duplicate-handling - Mark/remove duplicates; UMI-aware dedup
- alignment-files/sam-bam-basics - SAM flags, CIGAR, the cross-tool MAPQ scale, SA tags
- alignment-files/bam-statistics - flagstat/idxstats/stats QC gate; what a high mapping rate hides
- variant-calling/variant-calling - Call variants from the aligned BAM
- variant-calling/structural-variant-calling - SV calling from split/discordant reads (needs -Y)
