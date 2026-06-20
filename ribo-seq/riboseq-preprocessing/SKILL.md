---
name: bio-ribo-seq-riboseq-preprocessing
description: Preprocess ribosome profiling reads with UMI handling, adapter trimming, contaminant/rRNA depletion, and footprint-aware alignment. Use when preparing Ribo-seq FASTQ for periodicity QC, ORF detection, translation efficiency, or stalling analysis, or when deciding how to deduplicate, which aligner to use, or how to size-select ribosome-protected fragments.
tool_type: cli
primary_tool: STAR
---

## Version Compatibility

Reference examples tested with: cutadapt 4.4+, umi_tools 1.1+, STAR 2.7.11+, bowtie2 2.5.3+, SortMeRNA 4.3+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Ribo-seq Preprocessing

**"Preprocess my ribosome profiling data"** -> Extract UMIs, trim the 3' linker, deplete rRNA/tRNA contaminants, align footprints with end-to-end (non-soft-clipped) settings, deduplicate only when UMIs allow it, and QC the read-length distribution.
- CLI: `umi_tools extract` -> `cutadapt` -> `bowtie2`/`SortMeRNA` (contaminant removal) -> `STAR` (genome + transcriptome projection) -> `umi_tools dedup` -> `samtools`

The canonical modern order (nf-core/riboseq, McGlincy & Ingolia 2017) is UMI-extract FIRST (the UMI lives in the read and must move to the read name before the linker is cut), then trim, then contaminant removal (before the expensive aligner), then align, then dedup on the BAM.

## Upstream context that changes the analysis (ask before trusting the data)

- **How were cells harvested, and with which drug?** Cycloheximide (CHX) pre-treatment of live cells lets initiation continue while elongation arrests, fabricating start-codon and 5'-ramp density and distorting downstream dwell-time work (Hussmann 2015). Flash-freeze with no drug (or CHX only in the lysis buffer) is the gold standard. Harvest method is recorded at preprocessing because it gates which downstream conclusions are valid (see ribosome-stalling).
- **Which nuclease?** RNase I (eukaryotes) trims close to the ribosome with little sequence bias, giving sharp ~28-30 nt footprints and crisp periodicity. RNase I is inhibited by the E. coli ribosome and FAILS in bacteria, so bacterial protocols use micrococcal nuclease (MNase), which has sequence bias, broader footprints, and forces 3'-end P-site anchoring (Mohammad 2019). A eukaryote-tuned pipeline silently misanalyzes MNase/bacterial data.
- **Are there UMIs?** The dedup decision depends entirely on this (table below).

## The decisions that shape preprocessing

### Deduplication: with-UMI vs without-UMI (the load-bearing choice)

| Situation | What to do | Why |
|-----------|-----------|-----|
| Library has UMIs (McGlincy & Ingolia design or kit) | `umi_tools extract` before trim, `umi_tools dedup` on the BAM (`--method directional`) | UMI separates a true PCR duplicate (same position + length + UMI) from two independent ribosomes on the same codon (same position + length, different UMI) |
| No UMIs | Do NOT position-deduplicate; keep all reads | Many distinct ribosomes give identical 5' position AND identical footprint length; `markdup`/Picard would delete real footprints and flatten high-occupancy codons |
| Low input (single cells, scarce tissue, selective/IP profiling) | UMIs are essential | Few input molecules force heavy PCR; without UMIs amplified-once and amplified-1000x are indistinguishable |

### Alignment: genome (STAR, spliced) vs transcriptome (bowtie2, unspliced)

| Axis | Genome (STAR) | Transcriptome (bowtie2) |
|------|---------------|-------------------------|
| Splicing / novel junctions | Handles introns; required for junction-spanning footprints | Cannot span genomic introns; only annotated transcript cDNA |
| Multimapping | Lower (isoforms collapse to one locus) | High (every shared isoform + paralog multiplies hits) |
| Novel/uORF discovery | Strong (ribotricer/Ribo-TISH work off genome BAM + GTF) | Limited to annotated transcripts |
| P-site / periodicity coords | Project with `--quantMode TranscriptomeSAM` | Native transcript coords (convenient for riboWaltz) |
| Recommended | DEFAULT for mammals: STAR genome + transcriptome projection in one pass | Compact genomes (yeast) or when transcript-coordinate counts are the explicit goal |

### Contaminant removal approach

| Approach | Tool | Tradeoff |
|----------|------|----------|
| Combined-index depletion | bowtie2/STAR vs an rRNA+tRNA+snoRNA+snRNA FASTA, keep unmapped | Fast, full control of the contaminant set; the de-facto standard |
| Dedicated rRNA filter | SortMeRNA v4 (rRNA HMM/k-mer DBs) | rRNA-specialized but covers only rRNA; often paired with a separate ncRNA index |
| Layered (nf-core/riboseq) | BBSplit (broad) then SortMeRNA (rRNA) | Production-grade; most thorough |

rRNA is the dominant contaminant: commonly 50-90% (often >80%) of a Ribo-seq library, because nuclease digestion of the ribosome itself produces abundant rRNA fragments in the footprint size range. Wet-lab depletion (RiboZero/RiboCop/biotinylated subtraction oligos) reduces but never eliminates it, so in-silico removal is mandatory. Effective mRNA depth is a small fraction of raw reads.

## Extract UMIs

**Goal:** Move the UMI from the read sequence into the read name so it survives every later step and can deduplicate the final BAM.

**Approach:** Run `umi_tools extract` FIRST, before adapter trimming, with the barcode pattern matching the library's read structure (N = random UMI base extracted to the name, X = fixed base kept).

```bash
# Only when the library has UMIs. Pattern is library-specific.
# McGlincy & Ingolia 2017 split the 7-nt UMI (5 nt in the linker + 2 nt from circularization)
umi_tools extract \
    --bc-pattern=NNNNN \
    --stdin reads.fastq.gz \
    --stdout reads.umi.fastq.gz \
    --log umi_extract.log
```

When the UMI is split across the read (an inline 5' portion plus a portion inside the 3' linker, as in McGlincy & Ingolia 2017), the linker-embedded part is otherwise lost at trimming: extract it from the 3' end too (a second `umi_tools extract` with a `--3prime` pattern, or cutadapt's `{N}` linker capture) rather than discarding it. A pattern matching only the 5' inline bases recovers half the UMI and under-collapses duplicates.

## Trim the 3' linker

**Goal:** Remove the 3' adapter that is always read through because footprints (~28-30 nt) are far shorter than the read.

**Approach:** Run cutadapt with the known adapter and a PERMISSIVE length floor, and discard reads where no adapter was found.

```bash
# --discard-untrimmed: a footprint without read-through adapter is almost never a real footprint
# -m 15: permissive floor (do NOT narrow to 28-32 yet; inspect the length distribution first)
cutadapt \
    -a CTGTAGGCACCATCAAT \
    --discard-untrimmed \
    -m 15 -M 40 \
    -j 0 \
    -o reads.trimmed.fastq.gz \
    reads.umi.fastq.gz
```

The classic Ingolia linker `CTGTAGGCACCATCAAT` is an example only; the real sequence is protocol/kit-specific and McGlincy-Ingolia linkers embed the UMI and sample barcode, so the trimmed "adapter" region may include them.

## Remove rRNA and other contaminants

**Goal:** Discard rRNA/tRNA/snoRNA reads before the expensive spliced aligner runs.

**Approach:** Align to a combined contaminant index and keep only the unmapped reads, OR use a dedicated rRNA filter.

```bash
# Option A: combined contaminant index (rRNA + tRNA + snoRNA + snRNA), keep unmapped
bowtie2 -x contaminant_index \
    -U reads.trimmed.fastq.gz \
    --un-gz reads.noncontam.fastq.gz \
    -S /dev/null -p 8

# Option B: SortMeRNA v4 (use a per-sample --workdir; a shared kvdb collides across runs)
sortmerna \
    --ref rRNA_db/silva-euk-18s-id95.fasta \
    --ref rRNA_db/silva-euk-28s-id98.fasta \
    --reads reads.trimmed.fastq.gz \
    --aligned rRNA_hits --other reads.noncontam \
    --fastx --workdir sortmerna_sampleA --threads 8
```

## Align footprints (STAR, Ribo-seq-tuned)

**Goal:** Map cleaned footprints with settings appropriate for 28-30 nt reads, preserving the exact ends needed for P-site assignment.

**Approach:** Use STAR end-to-end (no soft-clipping), short-read seeding, a low mismatch cap, and transcriptome projection in one pass.

```bash
# --alignEndsType EndToEnd: the single most important Ribo-seq STAR flag.
#   STAR defaults to Local, which soft-clips footprint ends and corrupts P-site offsets.
# --seedSearchStartLmax 15: STAR's default 50 is wrong for ~30 nt reads.
# Do NOT set --alignIntronMax 1 on a genome (that forbids splicing and defeats STAR).
STAR --runMode alignReads \
    --genomeDir STAR_index \
    --readFilesIn reads.noncontam.fastq.gz \
    --readFilesCommand zcat \
    --alignEndsType EndToEnd \
    --seedSearchStartLmax 15 \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 10 --outSAMmultNmax 1 --outMultimapperOrder Random \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix sampleA_ --runThreadN 8

samtools index sampleA_Aligned.sortedByCoord.out.bam
```

Multimapping is higher in Ribo-seq than RNA-seq (paralogs, ncRNA, repeats). `--outFilterMultimapNmax 1` (unique-only) is simplest but silently drops translated paralogs/repeats; keeping a few multimappers with one random primary, or resolving by EM (RSEM) downstream, retains that signal. STAR's default `--outFilterScoreMinOverLread`/`--outFilterMatchNminOverLread` (0.66) are tuned for ~100 nt reads; very short footprints occasionally need these relaxed if good alignments are rejected.

## Deduplicate (only with UMIs)

**Goal:** Collapse PCR duplicates without destroying genuine co-occupancy.

**Approach:** Run `umi_tools dedup` on the aligned, sorted, indexed BAM; the `directional` method tolerates UMI sequencing errors.

```bash
# Run ONLY if the library has UMIs. Without UMIs, skip this entirely.
umi_tools dedup \
    --stdin sampleA_Aligned.sortedByCoord.out.bam \
    --stdout sampleA.dedup.bam \
    --method directional --log umi_dedup.log
samtools index sampleA.dedup.bam
```

Deduplicate WHICHEVER BAM the downstream step counts on. RiboCode and riboWaltz consume the transcriptome-projected BAM (`Aligned.toTranscriptome.out.bam`), so with UMIs that BAM must be deduplicated too (`umi_tools dedup --per-contig`, because reads sit on transcript "chromosomes"); deduplicating only the genome BAM leaves the ORF/periodicity inputs PCR-inflated. Note also that these blocks assume single-end reads (the Ribo-seq norm); paired-end kits that place the UMI on R2 need a different extract pattern.

## QC the preprocessing

**Goal:** Confirm the library captured real footprints before trusting any downstream analysis.

**Approach:** Plot the read-length distribution, report the contaminant fraction and mapping rate, and (with UMIs) the post-dedup complexity.

```bash
# Read-length distribution is THE key plot: expect a sharp mammalian peak ~28-30 nt,
# sometimes a ~20-22 nt shoulder (the open-A-site footprint population, Lareau 2014).
samtools view sampleA.dedup.bam | awk '{print length($10)}' | sort -n | uniq -c
samtools flagstat sampleA.dedup.bam
```

A permissive trim floor matters here: a tight 28-32 nt gate applied before this plot discards the ~20-22 nt population and hides QC problems. Size-select narrowly only after inspecting the distribution, and prefer per-read-length analysis downstream (riboWaltz, RiboFlow assign per-length P-site offsets).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Flat 33/33/33 frame downstream; weak periodicity | Footprint ends soft-clipped by STAR Local mode | Add `--alignEndsType EndToEnd`; never rely on STAR defaults for footprints |
| Junction-spanning footprints all lost | `--alignIntronMax 1` set on a genome alignment | Remove it (or only use it when the "genome" is a transcriptome FASTA) |
| High-occupancy codons look flattened after "dedup" | Position-based dedup on a library WITHOUT UMIs | Do not deduplicate without UMIs; same position + length is mostly real biology |
| SortMeRNA errors on the second sample | Shared default kvdb workdir collides across runs | Give each sample a fresh `--workdir` |
| Very few reads survive trimming | `--discard-untrimmed` plus a wrong adapter sequence | Confirm the actual linker (kit/protocol-specific); inspect a few raw reads |
| Mammalian peak missing, broad smear instead | Over-digestion, wrong size gate too early, or MNase data analyzed as RNase I | Plot length distribution first; for bacteria expect MNase breadth and 3'-anchoring |

## Related Skills

- ribosome-periodicity - Validate 3-nt periodicity and calibrate P-site offsets on the aligned BAM
- orf-detection - Detect translated ORFs once footprints are aligned and offsets known
- translation-efficiency - Needs matched RNA-seq processed consistently with the footprints
- read-qc/quality-reports - General read quality control before footprint-specific steps
- read-alignment/star-alignment - General STAR alignment background

## References

- Ingolia NT, Ghaemmaghami S, Newman JRS, Weissman JS. 2009. Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling. Science 324(5924):218-223. doi:10.1126/science.1168978
- McGlincy NJ, Ingolia NT. 2017. Transcriptome-wide measurement of translation by ribosome profiling. Methods 126:112-129. doi:10.1016/j.ymeth.2017.05.028
- Mohammad F, Green R, Buskirk AR. 2019. A systematically-revised ribosome profiling method for bacteria reveals pauses at single-codon resolution. eLife 8:e42591. doi:10.7554/eLife.42591
- Lareau LF, Hite DH, Hogan GJ, Brown PO. 2014. Distinct stages of the translation elongation cycle revealed by sequencing ribosome-protected mRNA fragments. eLife 3:e01257. doi:10.7554/eLife.01257
- Smith T, Heger A, Sudbery I. 2017. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Res 27(3):491-499. doi:10.1101/gr.209601.116
