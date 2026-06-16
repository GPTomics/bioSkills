---
name: bio-read-alignment-star-alignment
description: Aligns RNA-seq reads to a genome with STAR, the fast splice-aware aligner whose splice-junction database (built from a GTF at sjdbOverhang = readlength-1) and two-pass mode set junction sensitivity, whose 255-for-unique MAPQ breaks GATK, and whose GeneCounts output reveals library strandedness. Use when RNA reads must be placed on the genome for novel-isoform discovery, fusion detection, RNA variant calling, coverage tracks, splicing QC, or single-cell (STARsolo). Memory-constrained RNA alignment is hisat2-alignment; DE on known transcripts only should skip alignment for rna-quantification/alignment-free-quant; the QC gate and contig-naming reconciliation are alignment-files; counting is rna-quantification; DNA is bwa-alignment/bowtie2-alignment.
tool_type: cli
primary_tool: STAR
---

## Version Compatibility

Reference examples tested with: STAR 2.7.11+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# STAR RNA-seq Alignment -- The Junction Database, the 255 MAPQ, and the Strand Column Decide the Result

**"Align my RNA-seq reads"** -> Map reads across exon-exon junctions to the genome, building or reusing a splice-aware index from a GTF -- because an RNA read spans introns the genome does not contain, so the junction database (and its sjdbOverhang), the two-pass choice, the 255-unique MAPQ, and the strand column STAR reports are what actually determine the downstream counts and variant calls.
- CLI: `STAR --runMode alignReads --genomeDir idx/ --readFilesIn R1.fq.gz R2.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate`

Scope: RNA splice-aware mapping with STAR -- index generation, two-pass, GeneCounts, chimeric/fusion output, and the MAPQ/strandedness traps. Contig-naming reconciliation and the QC gate -> alignment-files. Low-memory RNA alignment -> hisat2-alignment. Counting reads over genes/transcripts -> rna-quantification. Single-cell droplet/plate counting (STARsolo) -> single-cell. DE on known transcripts without a BAM -> rna-quantification/alignment-free-quant. OUT OF SCOPE: DNA (bwa-alignment/bowtie2-alignment), long reads (long-read-sequencing/long-read-alignment).

## The Single Most Important Modern Insight

1. **The splice-junction database is part of the index/run config, and a wrong sjdbOverhang or per-sample two-pass silently changes the answer.** STAR builds short artificial junction-flank sequences from the GTF so reads with a short overhang on one exon can still be placed; `--sjdbOverhang` sets the flank length and should equal `max(readlength) - 1` (default 100 is fine near 100 bp reads but degrades junction sensitivity for very short reads). Two-pass (`--twopassMode Basic`) discovers novel junctions and re-aligns, raising novel-junction sensitivity -- but it is PER-SAMPLE: each sample is re-aligned against its own augmented index, so a junction found only in the deep/disease sample is rescued asymmetrically, a batch effect that confounds junction/splicing/sQTL comparisons. The cohort-correct recipe is to pool every sample's pass-1 `SJ.out.tab`, filter, and feed one common `--sjdbFileChrStartEnd` to a uniform second pass. Per-sample two-pass is fine for plain gene-level DE; it bites splicing analyses.
2. **STAR's MAPQ is 255-for-unique and a multiplicity code, so it breaks GATK and a copied MAPQ filter deletes every multimapper.** Unique reads get MAPQ 255 -- the SAM "mapping quality unavailable" value -- which GATK treats as missing and drops, the classic silently-empty RNA VCF; fix it at align time with `--outSAMmapqUnique 60`. Multimappers get 3 / 1 / 0 for 2 / 3-4 / >=5 loci (pure locus count, no score information), so a generic `samtools view -q 10` or `featureCounts -Q 10` copied from a DNA pipeline DELETES every multimapper while keeping every unique -- a directional "discard paralog/gene-family/rRNA/pseudogene reads" filter that biases against recently-duplicated gene families. A hard MAPQ filter is almost never what an RNA analysis wants.
3. **STAR reports library strandedness for free in GeneCounts, and getting strand wrong roughly halves counts.** `--quantMode GeneCounts` emits a 4-column `ReadsPerGene.out.tab`: gene, unstranded, forward, reverse. Summing columns 3 and 4 reveals the protocol -- roughly equal is unstranded (use col 2), col 3 dominant is forward, col 4 dominant is reverse (the common dUTP/TruSeq case, use col 4). Never assume the strand: feeding the wrong column (or the wrong `-s` to a counter) sends sense reads to "no feature" and roughly halves the counts, distorting DE.

## How STAR Places a Splice (the mechanism in brief)

STAR seeds with Maximal Mappable Prefix search on an uncompressed suffix array: it finds the longest exact prefix, and when that prefix ends (at a junction, mismatch, or read end) it restarts from the next base -- so a junction-spanning read naturally decomposes into seeds on two exons that STAR then stitches across the intron with an `N` (skipped-region) CIGAR. The stitch is scored with splice priors that penalize non-canonical motifs (`--scoreGapNoncan -8`), long introns (`--scoreGenomicLengthLog2scale -0.25`), and reward annotated junctions (`--sjdbScore 2`), plus a hard anchor floor (`--alignSJoverhangMin`): a novel junction on a tiny anchor carries almost no information, so it must clear both the soft score prior and the hard overhang minimum. This is why annotation (the sjdb) and two-pass matter for short-overhang and novel junctions.

## Tool Taxonomy

| Mode / output | Citation | Mechanism / role | When |
|---------------|----------|------------------|------|
| `STAR --runMode alignReads` | Dobin 2013 Bioinformatics 29:15 | MMP seed + stitch; spliced genomic BAM | the default RNA-to-genome alignment |
| `--twopassMode Basic` | Dobin 2013 (STAR manual) | discover novel junctions, re-align | novel-isoform / variant / splicing work (per-sample); pool for cohorts |
| `--quantMode GeneCounts` | STAR manual | per-gene counts + the strand-detection 3 columns | quick counts and strandedness inference |
| `--quantMode TranscriptomeSAM` | STAR manual | transcriptome-coord BAM for RSEM / Salmon-aln | isoform-level EM quantification -> rna-quantification |
| `--chimSegmentMin` -> Chimeric.out.junction | STAR manual | split-across-loci reads for fusion calling | STAR-Fusion / Arriba fusion detection |
| STARsolo (`--soloType`) | STAR manual | cell-barcode + UMI single-cell counting | scRNA-seq (route OUT) -> single-cell |
| HISAT2 | Kim 2019 Nat Biotechnol 37:907 | graph FM-index, ~1/4 the RAM | memory-constrained RNA (route OUT) -> hisat2-alignment |
| Salmon / kallisto | Patro 2017 Nat Methods 14:417 | alignment-free transcript quantification | DE on known transcripts only (route OUT) -> rna-quantification/alignment-free-quant |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| RNA-seq, ample RAM (>=32 GB), need a genomic BAM | STAR | fastest splice-aware aligner; native counts, fusions, 2-pass |
| RNA variant calling | STAR 2-pass + `--outSAMmapqUnique 60` then GATK SplitNCigarReads | 2-pass splices novel junctions; 60 avoids the 255 drop |
| Novel-isoform / splicing / sQTL across a cohort | cohort 2-pass (pool SJ.out.tab -> common `--sjdbFileChrStartEnd`) | per-sample 2-pass is a junction batch effect |
| Fusion detection | `--chimSegmentMin 12 --chimOutType ...` -> STAR-Fusion / Arriba | chimeric junctions are the fusion signal |
| Single-cell RNA | STARsolo (route OUT) | barcode+UMI counting -> single-cell |
| Memory-constrained (<32 GB) | route OUT to hisat2-alignment | STAR needs ~30 GB for human |
| DE on known transcripts only | route OUT to rna-quantification/alignment-free-quant | Salmon/kallisto are faster and model multimapping better |
| Small genome (bacterial/viral/plasmid) | STAR with reduced `--genomeSAindexNbases` | the default 14 silently builds a bad index / segfaults |

Default when uncertain: STAR with `--outSAMtype BAM SortedByCoordinate`, `--quantMode GeneCounts` (to also read strandedness), and `--twopassMode Basic` for per-sample novel-junction work; set `--outSAMmapqUnique 60` for any STAR -> GATK path.

## Generate the Genome Index

```bash
# sjdbOverhang = max(readlength) - 1 (149 for 2x150 reads). Default 100 degrades junctions for short reads.
STAR --runMode genomeGenerate --runThreadN 8 \
    --genomeDir star_index/ \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile annotation.gtf \
    --sjdbOverhang 149
# Small genome (e.g. 5 Mb): add --genomeSAindexNbases <= min(14, log2(GenomeLength)/2 - 1), or STAR segfaults.
```

## Basic Alignment

```bash
STAR --runThreadN 8 \
    --genomeDir star_index/ \
    --readFilesIn reads_1.fq.gz reads_2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix sample_ \
    --outSAMtype BAM SortedByCoordinate
samtools index sample_Aligned.sortedByCoord.out.bam
# STAR already coordinate-sorts -- a subsequent `samtools sort` is redundant. Single-end: one file in --readFilesIn.
```

## Two-Pass + Gene Counts + Strandedness

```bash
STAR --runThreadN 8 --genomeDir star_index/ \
    --readFilesIn r1.fq.gz r2.fq.gz --readFilesCommand zcat \
    --outFileNamePrefix sample_ \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --outSAMattrRGline ID:sample1 SM:sample1 PL:ILLUMINA LB:lib1 \
    --outSAMmapqUnique 60        # so a downstream GATK RNA-variant step does not drop the 255 uniques
# STAR read groups use --outSAMattrRGline (SPACE-separated tags), NOT bwa's tab-delimited -R '@RG\t...'.
# GATK requires read groups; comma-with-spaces separates groups for multiple --readFilesIn files.

# Detect strandedness from ReadsPerGene.out.tab (skip the 4 N_* summary rows, sum cols 3 vs 4):
awk 'NR>4 {f+=$3; r+=$4} END {printf "fwd(col3)=%d  rev(col4)=%d -> use col %s\n", f, r, (f>2*r?"3 fwd": r>2*f?"4 rev":"2 unstranded")}' sample_ReadsPerGene.out.tab
```

## ENCODE Long-RNA-seq Parameter Set

```bash
STAR --runThreadN 8 --genomeDir star_index/ --readFilesIn r1.fq.gz r2.fq.gz --readFilesCommand zcat \
    --outFileNamePrefix sample_ --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
    --sjdbScore 1 --outSAMattributes NH HI AS NM MD
# BySJout keeps only reads whose junctions passed the dataset-wide collapse; multimapNmax 20 retains real
# multi-locus genes; the 0.04 mismatch ratio scales with read length; intron caps at ~1 Mb cover human genes.
```

## Fusion Detection

```bash
# Chimeric junctions for STAR-Fusion (params per the STAR-Fusion wiki).
STAR --runThreadN 8 --genomeDir star_index/ --readFilesIn r1.fq.gz r2.fq.gz --readFilesCommand zcat \
    --outFileNamePrefix sample_ --outSAMtype BAM SortedByCoordinate \
    --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
    --chimOutType Junctions
# Arriba instead reads chimeric alignments from the BAM: --chimSegmentMin 10 --chimOutType WithinBAM SoftClip.
```

## Key Parameters (STAR defaults unless noted)

| Parameter | Default | Description |
|-----------|---------|-------------|
| --sjdbOverhang | 100 | junction-flank length at index build; set to readlength-1 |
| --twopassMode | None | `Basic` for per-sample novel-junction discovery |
| --outSAMmapqUnique | 255 | MAPQ for unique reads; set 60 for GATK |
| --outFilterMultimapNmax | 10 | max loci to report (ENCODE uses 20 for RNA) |
| --alignIntronMin / Max | 21 / 0 (auto) | gap < min is a deletion; cap Max ~1 Mb for human |
| --alignSJoverhangMin / SJDBoverhangMin | 5 / 3 | novel / annotated junction anchor floor (ENCODE 8 / 1) |
| --quantMode | -- | `GeneCounts` and/or `TranscriptomeSAM` |
| --chimSegmentMin | 0 (off) | turn on chimeric/fusion detection |
| --genomeSAindexNbases | 14 | reduce to min(14, log2(L)/2-1) for small genomes |
| --genomeLoad / --limitBAMsortRAM | NoSharedMemory / -- | shared-memory reuse; explicit sort RAM |

## Per-Method Failure Modes

### STAR 255 MAPQ into GATK
**Trigger:** STAR BAM (uniques at 255) fed to GATK RNA variant calling. **Mechanism:** GATK reads 255 as "MAPQ unavailable" and drops the read. **Symptom:** a silently empty or near-empty RNA VCF. **Fix:** `--outSAMmapqUnique 60` at align time (then SplitNCigarReads in the GATK RNA workflow).

### MAPQ filter deletes every multimapper
**Trigger:** a `-q 10` / `-Q 10` MAPQ filter (copied from DNA) on STAR output. **Mechanism:** STAR multimappers are MAPQ <= 3, uniques 255, so the filter keeps only uniques. **Symptom:** systematic under-counting of paralog/gene-family/rRNA/pseudogene loci; DE driven by multimapper fraction. **Fix:** do not MAPQ-filter RNA for counting; handle multimappers in the counter (NH-aware) or via EM -> rna-quantification/featurecounts-counting.

### Per-sample two-pass as a batch effect
**Trigger:** `--twopassMode Basic` per sample for a splicing/junction comparison. **Mechanism:** each sample is aligned against its own novel-junction-augmented index. **Symptom:** junction recovery correlates with depth/condition, confounding sQTL/differential splicing. **Fix:** pool pass-1 SJ.out.tab across the cohort, filter, feed one `--sjdbFileChrStartEnd` to a uniform second pass.

### sjdbOverhang mismatched to read length
**Trigger:** index built with default 100 for 36-50 bp reads, or a mismatch between index build and align values. **Mechanism:** junction flanks far longer than the reads degrade short-overhang sensitivity; a mismatch errors at align time. **Symptom:** reduced novel-junction-spanning sensitivity, or "present sjdbOverhang not equal to genome generation step." **Fix:** rebuild with `sjdbOverhang = max(readlength)-1`.

### Small genome with default SAindexNbases
**Trigger:** indexing a bacterial/viral/plasmid genome with `--genomeSAindexNbases 14`. **Mechanism:** the SA pre-index string is too long for the genome. **Symptom:** STAR silently builds a bad index or segfaults at align time. **Fix:** set `--genomeSAindexNbases min(14, log2(GenomeLength)/2 - 1)` (1 Mb -> 9, 100 kb -> 7).

### Index built with a different STAR version
**Trigger:** an index built months ago loaded by a bumped STAR module. **Mechanism:** STAR refuses an index whose versionGenome differs. **Symptom:** "Genome version is INCOMPATIBLE with running STAR version," or a subtly different version that loads and behaves differently across a cohort. **Fix:** rebuild the index with the exact aligning STAR version; pin the version for the whole cohort.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| sjdbOverhang = max(readlength) - 1 | STAR manual | sets the junction-flank length for short-overhang reads |
| --outSAMmapqUnique 60 for GATK | STAR manual / GATK RNA best practice | 255 is "unavailable" and is dropped by GATK |
| --outFilterMultimapNmax 20, --outFilterMismatchNoverReadLmax 0.04 | ENCODE long-RNA-seq pipeline | retains real multi-locus genes; mismatch budget scales with read length |
| --genomeSAindexNbases <= min(14, log2(L)/2 - 1) | STAR manual | the default 14 corrupts/segfaults small-genome indexes |
| STAR human-genome index RAM ~30 GB | STAR docs (approximate) | the reason to route memory-constrained jobs to HISAT2 |
| --chimSegmentMin 12 (STAR-Fusion) | STAR-Fusion wiki | minimum chimeric-segment length for fusion calling |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Empty / tiny RNA VCF after GATK | STAR 255 MAPQ dropped as "unavailable" | `--outSAMmapqUnique 60` (and SplitNCigarReads) |
| GATK "no read group" on a STAR BAM | STAR omits @RG unless asked | add `--outSAMattrRGline ID:.. SM:.. PL:ILLUMINA LB:..` (space-separated, not bwa's `-R`) |
| htseq-count silently miscounts STAR output | htseq wants name-sorted (or `-r pos`) input; STAR coordinate-sorts | emit `--outSAMtype BAM Unsorted` for htseq, or `samtools sort -n`; featureCounts accepts either |
| Counts ~halved, antisense artifacts | wrong strandedness column / `-s` | infer from GeneCounts cols 3 vs 4 (or RSeQC); TruSeq/dUTP = reverse |
| Paralog/gene-family genes under-counted | a `-q 10` MAPQ filter deleted multimappers | drop the MAPQ filter; handle NH>1 in the counter -> rna-quantification |
| "not enough memory for BAM sorting" | `--limitBAMsortRAM` too small | set it explicitly (e.g. 10000000000) |
| Segfault / bad index on a small genome | `--genomeSAindexNbases 14` | reduce per min(14, log2(L)/2 - 1) |
| "Genome version INCOMPATIBLE" | index built with another STAR version | rebuild with the running version; pin it |
| 0 counts despite high mapping rate | genome/GTF contig-naming mismatch | reconcile `chr1` vs `1`, `chrM` vs `MT` (same source/release) -> alignment-files |
| STAR reads `.gz` FASTQ as garbage / fails | STAR does not auto-detect gzip | add `--readFilesCommand zcat` (bwa/bowtie2/HISAT2 auto-detect gzip) |

## References

- Dobin A, Davis CA, Schlesinger F, et al. 2013. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29:15-21.
- Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. 2019. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol* 37:907-915.
- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. 2017. Salmon provides fast and bias-aware quantification of transcript expression. *Nat Methods* 14:417-419.
- Burset M, Seledtsov IA, Solovyev VV. 2000. Analysis of canonical and non-canonical splice sites in mammalian genomes. *Nucleic Acids Res* 28:4364-4375.

## Related Skills

- hisat2-alignment - Low-memory splice-aware alternative to STAR
- bwa-alignment - DNA short-read mapping (when reads do not cross junctions)
- read-qc/rnaseq-qc - RNA destination metrics: rRNA, gene-body coverage, strandedness
- read-qc/fastp-workflow - Trim adapters/poly-A before alignment
- alignment-files/bam-statistics - flagstat/idxstats QC gate; what a high mapping rate hides; contig naming
- rna-quantification/featurecounts-counting - Count aligned reads over genes (NH-aware multimapper handling)
- rna-quantification/alignment-free-quant - Salmon/kallisto when only known-transcript DE is needed
- differential-expression/deseq2-basics - Downstream DE from the count matrix
- single-cell/data-io - STARsolo single-cell counts into a single-cell workflow
