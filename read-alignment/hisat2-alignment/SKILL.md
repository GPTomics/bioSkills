---
name: bio-read-alignment-hisat2-alignment
description: Aligns RNA-seq reads to a genome with HISAT2, the splice-aware aligner whose hierarchical graph FM-index runs at roughly a quarter of STAR's memory (~7 GB for human), whose SNP/haplotype graph index reduces reference bias in the index itself, and whose MAPQ is GATK-friendly (60 for unique, no 255 problem). Use when RNA alignment must fit a memory-constrained machine, when feeding StringTie/Cufflinks transcript assembly via --dta, or when a SNP-aware graph index is wanted for allele-robust mapping. Feature-rich/high-RAM RNA alignment and fusion detection are star-alignment; DE on known transcripts only should skip alignment for rna-quantification/alignment-free-quant; the QC gate and contig-naming reconciliation are alignment-files; counting is rna-quantification.
tool_type: cli
primary_tool: HISAT2
---

## Version Compatibility

Reference examples tested with: hisat2 2.2+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# HISAT2 Alignment -- Graph-Indexed Spliced Mapping at a Quarter of STAR's Memory

**"Align my RNA-seq reads with low memory"** -> Map reads across exon-exon junctions with a hierarchical graph FM-index that fits a small machine -- because HISAT2 buys splice-aware alignment at ~7 GB instead of STAR's ~30 GB, its MAPQ is GATK-friendly, and its SNP-graph index can remove reference bias before a single read is mapped.
- CLI: `hisat2 -p 8 -x index -1 R1.fq.gz -2 R2.fq.gz | samtools sort -@4 -o aligned.bam -`

Scope: low-memory RNA splice-aware mapping with HISAT2 -- index building (plain / annotation-aware / SNP-graph), strandedness, the --dta transcript-assembly mode, and manual two-pass. Contig naming and the QC gate -> alignment-files. Feature-rich/high-RAM RNA alignment, native gene counts, and fusion detection -> star-alignment. Counting reads over genes -> rna-quantification. DE without a BAM -> rna-quantification/alignment-free-quant. OUT OF SCOPE: DNA (bwa-alignment/bowtie2-alignment), long reads (long-read-sequencing/long-read-alignment), HLA typing (HISAT-genotype, a separate tool).

## The Single Most Important Modern Insight

1. **The hierarchical graph FM-index is why HISAT2 exists: near-STAR spliced alignment at ~1/4 the RAM.** HISAT2 uses one global FM-index to anchor a read plus ~55,000 small local graph FM-indexes (each ~56 kb), and extends a spliced read within the relevant local index rather than stitching genome-wide as STAR does. Most introns fit inside one local window, so spliced extension is a cheap local operation -- the resident human index is ~4-7 GB vs STAR's ~30 GB. That memory win is the reason to choose HISAT2; the cost is slightly lower novel-junction sensitivity than STAR two-pass and no native gene counts or fusion output.
2. **The SNP/haplotype graph index removes reference bias in the index, and the MAPQ is GATK-friendly.** A `hisat2-build --snp --haplotype` (or the prebuilt `grch38_snp` index) encodes millions of known variants as alternate graph nodes, so a read carrying a known alt allele traverses the alt node with no mismatch penalty -- the bias that over-counts the reference allele is removed structurally, for all those sites at once, without a per-sample personalized reference. (Private/novel variants still cause bias, so rigorous ASE still needs WASP or a personalized reference.) HISAT2 also assigns unique reads MAPQ 60 (not STAR's 255), so its output goes into GATK without the reassignment STAR needs.
3. **--dta is for transcript assembly only, and using it for plain counting throws away reads.** `--dta` raises the minimum anchor length required to report a de-novo spliced alignment, deliberately suppressing short-anchor junction reads -- because StringTie/Cufflinks cannot reliably assemble a transcript from a 3-5 bp anchor and such reads produce spurious isoforms. That trades junction sensitivity for assembly cleanliness, so `--dta` belongs only in a transcript-assembly pipeline; for plain gene counting it just discards usable junction reads. Strandedness (`--rna-strandness RF` for the common dUTP/TruSeq case) must also be set, or sense reads land in "no feature" and counts roughly halve.

## How HISAT2 Splices (the mechanism in brief)

A read is seeded by the global FM-index, then the relevant ~56 kb local FM-index is selected and the read is extended across the junction within it: the unaligned remainder is anchored in the local index and extended by repeated FM-index extension. Because the spliced extension is a narrow, local operation rather than a genome-wide seed-cluster-stitch, HISAT2 needs far less RAM than STAR -- and evaluates a narrower set of candidate splice configurations, which is the source of both its speed/memory advantage and its slightly lower novel-junction sensitivity.

## Tool Taxonomy

| Mode / index | Citation | Mechanism / role | When |
|--------------|----------|------------------|------|
| `hisat2-build` (plain) | Kim 2019 Nat Biotechnol 37:907 | genome-only HGFM | quick index; junctions supplied at align time |
| `hisat2-build --ss --exon` | Kim 2019 | annotation-aware HGFM (better short-anchor placement) | when build RAM allows; or use prebuilt `*_tran` indexes |
| `hisat2-build --snp --haplotype` | Kim 2019 | SNP/haplotype graph (reference-bias reduction) | allele-robust mapping; the `grch38_snp` index |
| `hisat2` alignReads | Kim 2019 | spliced alignment via local FM-index extension | the default RNA-to-genome mapping |
| `--dta` / `--dta-cufflinks` | HISAT2 manual | longer-anchor reporting for assemblers | StringTie / Cufflinks transcript assembly ONLY |
| manual two-pass (`--novel-splicesite-*`) | HISAT2 manual | discover then reuse novel junctions | novel-junction sensitivity (cohort: merge across samples) |
| STAR | Dobin 2013 Bioinformatics 29:15 | higher RAM, native counts, fusions, 2-pass | feature-rich RNA (route OUT) -> star-alignment |
| Salmon / kallisto | Patro 2017 Nat Methods 14:417 | alignment-free quantification | DE on known transcripts only (route OUT) -> rna-quantification/alignment-free-quant |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| RNA-seq on a memory-constrained machine (<32 GB) | HISAT2 | ~7 GB graph index vs STAR's ~30 GB |
| StringTie/Cufflinks transcript assembly | HISAT2 `--dta` | longer-anchor reporting the assemblers need |
| Allele-robust mapping / known-variant-aware | HISAT2 SNP-graph index (`grch38_snp`) | alt-allele reads traverse graph nodes without penalty |
| RNA variant calling | HISAT2 (MAPQ 60) then GATK SplitNCigarReads | GATK-friendly MAPQ, no 255 reassignment |
| Need native gene counts, fusions, or top novel-junction sensitivity | route OUT to star-alignment | HISAT2 has no GeneCounts/chimeric output |
| DE on known transcripts only | route OUT to rna-quantification/alignment-free-quant | Salmon/kallisto are faster and model multimapping |
| Plain gene-level counting | HISAT2 without `--dta` | `--dta` discards short-anchor junction reads |

Default when uncertain: HISAT2 with `--rna-strandness RF` (verify the strand), streamed to a coordinate-sorted BAM; add `--dta` only for transcript assembly.

## Build Index

```bash
# Plain genome-only index (cheap; supply junctions at align time with --known-splicesite-infile).
hisat2-build -p 8 reference.fa hisat2_index

# Annotation-aware (better short-anchor placement). NOTE: a full human --ss --exon build needs a LOT of RAM;
# prefer the prebuilt grch38_tran / grch38_snp_tran indexes, or pass junctions at align time instead.
hisat2_extract_splice_sites.py annotation.gtf > splice_sites.txt
hisat2_extract_exons.py        annotation.gtf > exons.txt
hisat2-build -p 8 --ss splice_sites.txt --exon exons.txt reference.fa hisat2_index
```

## Basic Alignment with Strandedness

```bash
# RF = reverse-stranded (dUTP / Illumina TruSeq Stranded mRNA -- the common case). Verify, do not assume.
hisat2 -p 8 -x hisat2_index --rna-strandness RF \
    --rg-id sample1 --rg SM:sample1 --rg PL:ILLUMINA \
    -1 reads_1.fq.gz -2 reads_2.fq.gz \
    --new-summary --summary-file sample.summary.txt | \
    samtools sort -@ 4 -o aligned.sorted.bam -
samtools index aligned.sorted.bam
# Single-end stranded: --rna-strandness R (reverse) or F (forward). Unstranded: omit the flag.
```

## For StringTie / Cufflinks (transcript assembly)

```bash
# --dta reports longer anchors the assemblers need; use ONLY for assembly, not for plain counting.
hisat2 -p 8 -x hisat2_index --rna-strandness RF --dta \
    -1 r1.fq.gz -2 r2.fq.gz | samtools sort -@ 4 -o aligned.bam -
```

## Manual Two-Pass (cohort novel-junction discovery)

```bash
# Pass 1: discover novel junctions per sample.
for r1 in *_R1.fq.gz; do
    base=$(basename "$r1" _R1.fq.gz); r2=${r1/_R1/_R2}
    hisat2 -p 8 -x hisat2_index --novel-splicesite-outfile "${base}.novel.txt" \
        -1 "$r1" -2 "$r2" -S /dev/null
done
# Merge across the cohort so every sample sees the same junction set (avoids a per-sample junction batch effect).
cat *.novel.txt | sort -u > cohort.novel.txt
# Pass 2: re-align every sample with the shared novel-junction set.
for r1 in *_R1.fq.gz; do
    base=$(basename "$r1" _R1.fq.gz); r2=${r1/_R1/_R2}
    hisat2 -p 8 -x hisat2_index --rna-strandness RF --novel-splicesite-infile cohort.novel.txt \
        -1 "$r1" -2 "$r2" | samtools sort -@ 4 -o "${base}.bam" -
done
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| -x | -- | index BASENAME |
| -1 / -2 / -U | -- | paired / single-end reads |
| --rna-strandness | unstranded | FR / RF / F / R (dUTP/TruSeq = RF / R) |
| --dta / --dta-cufflinks | off | longer anchors for StringTie / Cufflinks (assembly only) |
| --known-splicesite-infile | -- | supply junctions at align time (cheap-index alternative to --ss build) |
| --novel-splicesite-outfile / -infile | -- | manual two-pass |
| --max-intronlen | 500000 | shorter than STAR's effective ~1 Mb; raise for long-intron genes |
| -k | 5 (HFM) / 10 (HGFM) | max alignments reported per read |
| --no-softclip / --no-spliced-alignment | off | force end-to-end / disable splicing (DNA mode) |

## Per-Method Failure Modes

### --dta used for plain counting
**Trigger:** `--dta` on a run whose downstream is featureCounts/htseq, not StringTie. **Mechanism:** `--dta` suppresses short-anchor junction reads. **Symptom:** lower junction-read recovery and counts than a non-dta run. **Fix:** drop `--dta` for counting; keep it only for transcript assembly.

### Wrong strandedness
**Trigger:** omitting or mis-setting `--rna-strandness`. **Mechanism:** the XS strand tag is mislabeled and sense reads are assigned to "no feature." **Symptom:** counts ~halved; StringTie builds transcripts on the wrong strand. **Fix:** infer strand (RSeQC infer_experiment.py, or STAR GeneCounts) and set RF for dUTP/TruSeq.

### --ss --exon human build runs out of RAM
**Trigger:** a full human annotation-aware build on a small machine. **Mechanism:** building the annotation-aware HGFM needs far more RAM than a plain build. **Symptom:** the build is killed (OOM). **Fix:** use a prebuilt `grch38_tran`/`grch38_snp_tran` index, or build plain and pass junctions at align time via `--known-splicesite-infile`.

### max-intronlen too small for long-intron genes
**Trigger:** the default `--max-intronlen 500000` on genes with introns near or above ~1 Mb. **Mechanism:** junctions longer than the cap are not formed. **Symptom:** long-gene junction reads soft-clipped or mismapped. **Fix:** raise `--max-intronlen` for organisms/genes with very long introns.

### Genome/GTF contig-naming mismatch
**Trigger:** the BAM uses `chr1`/`chrM` but the counting GTF uses `1`/`MT`. **Mechanism:** no overlapping features. **Symptom:** zero counts despite a high alignment rate. **Fix:** reconcile naming (same source/release) -> alignment-files.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| HISAT2 human graph index RAM ~4.3 GB plain / ~6.7 GB SNP | Kim 2019 (approximate) | the ~1/4-of-STAR footprint that motivates choosing HISAT2 |
| --max-intronlen 500000 default | HISAT2 manual | shorter than STAR's ~1 Mb; raise for long-intron genes |
| --rna-strandness RF for dUTP/TruSeq | library-prep chemistry | the overwhelmingly common stranded protocol |
| unique-read MAPQ 60 (since v2.0.4) | HISAT2 manual / changelog | GATK-friendly; no 255 reassignment needed |
| -k 5 (HFM) / 10 (HGFM) | HISAT2 manual | max reported alignments differs by index type |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Counts ~halved, wrong-strand transcripts | missing/incorrect `--rna-strandness` | infer strand; set RF for dUTP/TruSeq |
| Lower counts than expected | `--dta` used for plain counting | drop `--dta` unless assembling transcripts |
| `--ss --exon` build killed (OOM) | full human annotation-aware build | use a prebuilt index or `--known-splicesite-infile` at align time |
| Long-gene junction reads clipped | `--max-intronlen` too small | raise it for long-intron genes |
| 0 counts despite high alignment rate | genome/GTF contig-naming mismatch | reconcile `chr1` vs `1` (same source/release) -> alignment-files |
| "Could not locate a HISAT2 index" | `-x` given a `.ht2` file | pass the index basename |
| htseq-count miscounts HISAT2 output | htseq wants name-sorted input | pipe to `samtools sort -n` for htseq; featureCounts accepts coordinate order |

## References

- Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. 2019. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol* 37:907-915.
- Kim D, Langmead B, Salzberg SL. 2015. HISAT: a fast spliced aligner with low memory requirements. *Nat Methods* 12:357-360.
- Dobin A, Davis CA, Schlesinger F, et al. 2013. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29:15-21.
- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. 2017. Salmon provides fast and bias-aware quantification of transcript expression. *Nat Methods* 14:417-419.

## Related Skills

- star-alignment - Feature-rich, higher-RAM splice-aware alternative (native counts, fusions)
- bwa-alignment - DNA short-read mapping (when reads do not cross junctions)
- read-qc/rnaseq-qc - RNA destination metrics: rRNA, gene-body coverage, strandedness
- read-qc/fastp-workflow - Trim adapters/poly-A before alignment
- alignment-files/bam-statistics - flagstat/idxstats QC gate; what a high mapping rate hides; contig naming
- rna-quantification/featurecounts-counting - Count aligned reads over genes
- rna-quantification/alignment-free-quant - Salmon/kallisto when only known-transcript DE is needed
- differential-expression/deseq2-basics - Downstream DE from the count matrix
