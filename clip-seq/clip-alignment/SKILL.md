---
name: bio-clip-seq-clip-alignment
description: Align preprocessed CLIP-seq reads (eCLIP, iCLIP, iCLIP2, PAR-CLIP) to genome with STAR or bowtie2 using crosslink-preserving parameters, choosing between unique-mapper-only and multi-mapper-aware alignment for repeat-binding RBPs, deciding STAR vs HISAT2 memory trade-offs, and applying ENCODE-compatible filters. Use when turning preprocessed CLIP FASTQ into a deduplicated, MAPQ-filtered BAM ready for peak calling or crosslink-site detection.
tool_type: cli
primary_tool: STAR
---

## Version Compatibility

Reference examples tested with: STAR 2.7.11b+, bowtie2 2.5.3+, HISAT2 2.2.1+, samtools 1.19+, CLAM 1.2+, umi_tools 1.1.5+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed binary (`<tool> -h`) and adapt the example to match the actual CLI rather than retrying.

# CLIP-seq Alignment

**"Align preprocessed CLIP reads to genome with crosslink-preserving parameters"** -> Map UMI-extracted, adapter-trimmed reads to the genome (NOT transcriptome) with end-to-end alignment, strict mismatch ceiling, and unique-mapper-only filtering by default. The 5' end of the read (R2 5' in paired-end eCLIP; R1 5' in iCLIP) carries the reverse-transcriptase truncation = crosslink site -1; any soft-clipping or 5' trimming during alignment destroys nucleotide resolution.

- CLI (eCLIP / iCLIP / iCLIP2, ENCODE pattern): `STAR --runMode alignReads --genomeDir STAR_index --readFilesIn R1.trim.fq.gz R2.trim.fq.gz --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 1 --alignEndsType EndToEnd --outFilterMismatchNoverReadLmax 0.04 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFileNamePrefix sample_`
- CLI (PAR-CLIP, mismatch ceiling raised for T->C signal): same as above but `--outFilterMismatchNoverReadLmax 0.07`
- CLI (low memory, no splicing): `bowtie2 -x genome_index -U R1.trim.fq.gz --very-sensitive -p 8 | samtools view -bS - | samtools sort -o aligned.bam`
- CLI (repeat-binding RBP, multi-mapper rescue): `STAR ... --outFilterMultimapNmax 100 --outSAMmultNmax -1` then process with `CLAM` (see repeat-element section)

End-to-end alignment is non-negotiable. `--alignEndsType Local` (STAR default for some pipelines) soft-clips low-quality 5' bases and discards exactly the truncation signal. Mismatch ceiling 0.04 (4% of read length) excludes most sequencing-error reads; 0.07 is the PAR-CLIP override to retain T->C reads.

## Algorithmic Taxonomy

| Aligner | Splice-aware | Memory (human) | CLIP-suitable | Strength | Fails when |
|---------|--------------|----------------|---------------|----------|------------|
| STAR | Yes | ~30 GB peak (sjdbOverhang 100) | Yes (ENCODE eCLIP standard) | Splice-aware, fast on deep libraries, excellent multi-mapper logs | RAM hungry; small clusters or laptops cannot run human genome index |
| HISAT2 | Yes | ~8 GB | Yes (good alternative when STAR memory infeasible) | Low memory; comparable splice accuracy | Slightly lower multi-mapper precision; smaller community adoption for CLIP |
| bowtie2 | No (no splice) | ~3 GB | Limited (loses intronic and read-through reads at splice junctions) | Fast, low memory, mature | Misses spliced reads; not suitable for mRNA-binding RBPs (PTBP1, U2AF2) |
| bwa-mem2 | No | ~30 GB index | Not recommended | Fast for DNA; no splice support | Same splicing issue as bowtie2; no CLIP-specific advantage |
| chromap | No | ~5 GB | Not recommended for CLIP | Very fast for ATAC/ChIP | No splice support; pre-applies fragment shift inappropriate for CLIP |
| novoalign | Yes (splice with -X) | ~10 GB | Acceptable | Old standard; some labs still use | Commercial; community moved to STAR/HISAT2 |
| Salmon (transcriptome) | N/A | low | Not suitable | Pseudo-alignment | Discards intronic and unannotated reads, both common in CLIP |

Methodology evolves; verify the current ENCODE eCLIP pipeline (encodeproject.org/eclip) before locking parameters. As of 2025, ENCODE 4 still uses STAR 2.4.0j (legacy version pinning); modern reanalyses use STAR 2.7.x with the same parameter set.

## Critical Choice: Unique-Mapper-Only vs Multi-Mapper Rescue

Two valid alignment strategies depend on the RBP biology:

**Strategy A -- Unique-mapper-only (default, ENCODE standard):** `--outFilterMultimapNmax 1` discards every read that maps to more than one genomic location. Loses 5-15% of reads but produces unambiguous positions for downstream peak calling and crosslink-site detection.

**Strategy B -- Multi-mapper rescue (for repeat-binding RBPs):** `--outFilterMultimapNmax 100 --outSAMmultNmax -1` retains up to 100 alignments per read; downstream EM-based assignment (CLAM, Xinglab) probabilistically allocates them. Required for RBPs that bind Alu, LINE-1, LTR, or other repetitive elements (e.g., MATR3, ILF3, FUS at LINE-1, HNRNPK at SINEs, PUM2 in some repeat contexts).

| RBP class | Strategy | Why |
|-----------|----------|-----|
| Splicing factors (PTBP1, U2AF2, RBFOX, SRSF1) | A (unique) | Bind defined intronic/exonic motifs, not repeats |
| 3' UTR stability (HuR, PUM2, AUF1) | A (unique) | Binding sites are usually in unique 3' UTR sequence |
| Ribosomal proteins (RPS19, RPL35A) | A (unique) | Bind mRNA bodies and snoRNAs in unique sequence |
| Translation initiation (EIF3J, EIF2S2) | A (unique) | Bind 5' UTRs and snoRNAs |
| Repeat / TE binders (MATR3, ZFP36 isoforms, HNRNPK, LINE-1 ORFs) | B (multi-mapper) | Genuine biology is in repeat regions |
| Y-RNA / 7SK / vault RNAs (TROVE2, LARP7) | A (unique) but with `--outFilterMultimapNmax 50` for the ncRNA itself | These small RNAs have a few unique copies; pure unique discards them |
| Mitochondrial RBPs (FASTKD2, LRPPRC, TFAM) | A (unique) with chrM retained | chrM has unique sequence but must not be excluded by blacklists |
| Histone mRNA (SLBP) | A (unique) but DO NOT pre-filter rRNA index that includes histone | Replication-dependent histone genes are repetitive in some indices |

## Per-Aligner Failure Modes

### STAR -- Local alignment soft-clips the truncation site

**Trigger:** Pipeline copied from RNA-seq tutorial uses `--alignEndsType Local` (or omits the flag, letting STAR default).

**Mechanism:** Local alignment soft-clips up to 12% of the read end if scoring improves. In CLIP, the 5' read end is the truncation = CL site -1 base; if it carries even one mismatched base from RT errors, STAR soft-clips it.

**Symptom:** Crosslink-site density looks "smoothed"; PureCLIP / CTK CITS detect 5-30% fewer single-nt sites than expected; peak boundaries look fuzzy.

**Fix:** Always pass `--alignEndsType EndToEnd`. To confirm: `samtools view dedup.bam | awk '{ if ($6 ~ /S/) print }' | wc -l` should report < 1% of reads with soft-clip CIGAR operations.

### STAR -- Mismatch ceiling discards PAR-CLIP signal

**Trigger:** PAR-CLIP library aligned with `--outFilterMismatchNoverReadLmax 0.04` (the iCLIP/eCLIP default).

**Mechanism:** PAR-CLIP T->C conversion rate is 20-50% of T positions in crosslinked reads. For a 30 nt read with 8 Ts and 50% conversion, ~4 T->C mismatches occur (13% of read length). The 4% mismatch ceiling drops these reads.

**Symptom:** 40-70% read loss at alignment for PAR-CLIP; downstream PARalyzer / wavClusteR find few clusters.

**Fix:** For PAR-CLIP only, set `--outFilterMismatchNoverReadLmax 0.07` (7%). Verify post-alignment: `samtools view dedup.bam | awk '{ for(i=12;i<=NF;i++) if($i ~ /^MD:/) print $i }' | head` should show many T->C-indicating MD tags.

### STAR -- Multi-mappers silently discarded when needed

**Trigger:** Repeat-binding RBP (MATR3, LINE-1 binders) aligned with default `--outFilterMultimapNmax 1`.

**Mechanism:** Reads mapping to more than one Alu/LINE/LTR instance are removed; 10-30% of true binding sites are lost.

**Symptom:** Compared to published RBP literature for repeat binders, the analysis peak count is 3-10x lower; repeat-overlap fraction is < 5% when literature suggests 15-30%.

**Fix:** Raise `--outFilterMultimapNmax` to 50-100; emit all alignments with `--outSAMmultNmax -1`; downstream use CLAM (Xinglab) to probabilistically assign multi-mappers to a single best location via expectation-maximization. Or restrict to a non-repeat analysis and note the limitation.

### bowtie2 -- Splice junctions missed silently

**Trigger:** Used `bowtie2` instead of STAR for an mRNA-binding RBP (typical of PTBP1, U2AF2 intronic binding studies); BAM looks normal but introns/exons map poorly.

**Mechanism:** bowtie2 has no splice model; reads spanning exon-intron junctions soft-clip 5-20 nt or fail to align entirely.

**Symptom:** Lower read recovery (~70-80% vs STAR ~90-95%); peaks at splice sites under-called; intronic peaks over-represented (reads that would have aligned across an exon now align fully in the intron upstream).

**Fix:** Use STAR or HISAT2 for any RBP that binds mRNA, pre-mRNA, or splice signals. Reserve bowtie2 for non-coding RNA targets (snoRNA, 7SK, Y RNA) where short reads do not span junctions.

### STAR -- Memory exhausted on small machine

**Trigger:** Human genome index loaded into < 32 GB RAM machine.

**Mechanism:** STAR's suffix-array index requires ~30 GB RAM for human/mouse. Smaller machines crash or swap.

**Symptom:** "STAR EXITED" with bus error; or alignment takes > 24h on a single sample.

**Fix:** Switch to HISAT2 (~8 GB peak); or build a STAR index with `--genomeSAsparseD 2` (halves RAM at the cost of ~30% alignment slowdown); or use a public cluster with > 64 GB. Do NOT use bowtie2 as a memory workaround for splice-aware needs.

### Read-2 5' end trimmed inadvertently

**Trigger:** A `--clip5pNbases` flag carried over from RNA-seq quality trimming; or a `fastp --trim_front2 5` step in preprocessing.

**Mechanism:** Removes the eCLIP truncation base on R2 5'.

**Symptom:** Crosslink sites called from R2 5' positions cluster artificially at uniform offsets across genome; motif enrichment around crosslink sites collapses.

**Fix:** Verify the BAM with `samtools view dedup.bam | awk '$2 ~ /83|163/ { print $4 }'` (R2 5' positions on minus and plus strand) and confirm they map to the EXPECTED truncation sites, not shifted by 5 nt.

## ENCODE 4 eCLIP STAR Parameters (Reference)

```bash
STAR --runMode alignReads \
    --runThreadN 16 \
    --genomeDir /path/to/STAR_hg38_index \
    --genomeLoad NoSharedMemory \
    --readFilesIn R1.trim.fq.gz R2.trim.fq.gz \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --outFilterMultimapNmax 1 \
    --alignEndsType EndToEnd \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes All \
    --outFileNamePrefix sample_ \
    --outFilterScoreMinOverLread 0.66 \
    --outFilterMatchNminOverLread 0.66
```

The `--outFilterScoreMinOverLread 0.66` / `--outFilterMatchNminOverLread 0.66` (matched score / matched bases >= 66% of read length) are ENCODE-specific stringency added on top of the mismatch ceiling. These two flags are why ENCODE eCLIP BAMs are smaller than naive STAR output.

## Post-Alignment Filtering

```bash
# Sort and index
samtools index sample_Aligned.sortedByCoord.out.bam

# MAPQ filter (255 in STAR = unique; bowtie2 uses different scheme)
samtools view -b -q 10 sample_Aligned.sortedByCoord.out.bam > sample_q10.bam
samtools index sample_q10.bam

# UMI deduplication (see clip-preprocessing for `--method=unique` rationale)
umi_tools dedup \
    --stdin=sample_q10.bam \
    --stdout=sample_dedup.bam \
    --method=unique \
    --paired \
    --log=sample_dedup.log
samtools index sample_dedup.bam
```

**MAPQ thresholds differ by aligner:**
- STAR: 255 = uniquely mapped (== `--outFilterMultimapNmax 1` already filtered); lower MAPQ values mean multi-mappers
- bowtie2: 42 = unique; 0-1 = multi-mapper; intermediate = ambiguous
- HISAT2: 60 = unique; similar to bowtie2 scheme but tool-specific

For STAR `-q 10` is conventional but redundant if `--outFilterMultimapNmax 1` was already set. For bowtie2/HISAT2, `-q 30` is a stricter unique-mapper proxy.

## Multi-Mapper Rescue with CLAM (Repeat-Binding RBPs)

```bash
# 1. Re-align permitting multi-mappers
STAR --runMode alignReads \
    --genomeDir STAR_index --readFilesIn R1.fq.gz R2.fq.gz \
    --readFilesCommand zcat --outFilterMultimapNmax 100 \
    --outSAMmultNmax -1 --alignEndsType EndToEnd \
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix mm_

samtools index mm_Aligned.sortedByCoord.out.bam

# 2. CLAM preprocessing splits unique vs multi-mapper reads
CLAM preprocessor -i mm_Aligned.sortedByCoord.out.bam -o clam_out/ --read-tagger-method median

# 3. EM-based multi-mapper realignment
CLAM realigner -i clam_out/unique.sorted.bam -o clam_out/ --winsize 50 --max-tags 0

# 4. Downstream peakcaller (CLAM peakcaller is a wrapper for Piranha+EM)
CLAM peakcaller -i clam_out/unique.sorted.bam clam_out/realigned.sorted.bam \
    -o clam_peaks.bed -p 8 --gtf gencode.v38.annotation.gtf
```

CLAM rescues 10-30% additional peaks in repeat regions (Zhang & Xing 2017). It is the only EM-based multi-mapper solution actively maintained for CLIP as of 2025.

## HISAT2 Low-Memory Alternative

```bash
hisat2 --rna-strandness FR --no-softclip \
    -p 8 -x hisat2_grch38_index \
    -1 R1.trim.fq.gz -2 R2.trim.fq.gz \
    --no-unal \
    --score-min L,0,-0.2 \
    -S sample.sam 2> hisat2.log

samtools view -bS sample.sam | samtools sort -o sample.bam -
samtools index sample.bam
```

`--no-softclip` is the HISAT2 equivalent of STAR's `--alignEndsType EndToEnd`. `--score-min L,0,-0.2` is roughly equivalent to STAR's `--outFilterMismatchNoverReadLmax 0.04` (penalty scales linearly with read length, with -0.2 per mismatch slope).

## Decision Tree by Use Case

| Scenario | Recommended aligner + parameters | Why |
|----------|----------------------------------|-----|
| eCLIP, ENCODE-comparable, 32+ GB RAM | STAR ENCODE pattern (above) | Reproducible against ENCODE peak calls |
| iCLIP / iCLIP2, single-end | STAR same params, `--readFilesIn R1.fq.gz` | Single-end variant of ENCODE block |
| PAR-CLIP | STAR with `--outFilterMismatchNoverReadLmax 0.07` | T->C is signal, not error |
| Repeat-binding RBP (MATR3, LINE-1 binders) | STAR `--outFilterMultimapNmax 100 --outSAMmultNmax -1` + CLAM | EM rescue of multi-mappers |
| Low memory (< 16 GB) | HISAT2 with `--no-softclip` | STAR index unloadable |
| Bacterial / yeast (small genome, no splicing) | bowtie2 `--very-sensitive` | Splice-awareness wasted |
| snoRNA-only / 7SK-only analysis | bowtie2 with custom index | Avoid genome-scale splicing tangle |
| Allele-specific CLIP | STAR ENCODE + WASP filter (`--waspOutputMode SAMtag`) | Reference-allele mapping bias must be removed |
| Long-read CLIP (dirCLIP, nanopore) | minimap2 `-ax splice -uf -k14` | STAR cannot align long reads |

## Allele-Specific Alignment with WASP

CLIP-seq inherits reference-allele mapping bias from RNA-seq. For allele-specific binding analyses (BEAPR, ASPRIN), STAR's WASP integration removes reads where the alternative-allele version of the read would not have aligned at the same position.

```bash
# STAR with WASP filter; --varVCFfile is the heterozygous SNP VCF
STAR --runMode alignReads \
    --genomeDir STAR_index --readFilesIn R1.fq.gz R2.fq.gz --readFilesCommand zcat \
    --alignEndsType EndToEnd --outFilterMultimapNmax 1 --outFilterMismatchNoverReadLmax 0.04 \
    --outSAMtype BAM SortedByCoordinate \
    --varVCFfile sample.het.vcf.gz \
    --waspOutputMode SAMtag \
    --outSAMattributes vA vG vW \
    --outFileNamePrefix wasp_

samtools view -b -e '[vW]==1' wasp_Aligned.sortedByCoord.out.bam > wasp_pass.bam
```

WASP filtering is MANDATORY for allele-specific binding analyses; reference-allele bias inflates REF allele frequency 1-5% in unfiltered CLIP data.

## Reconciliation: When Alignment Outputs Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| STAR retains more reads than bowtie2 | STAR handles splice junctions; bowtie2 fails on spliced reads | Trust STAR for mRNA-binding RBPs |
| Read count after STAR < expected | `--outFilterMultimapNmax 1` discarding multi-mappers | If RBP binds repeats, switch to multi-mapper mode + CLAM |
| Crosslink-site density looks smoothed | Soft-clip ON (alignEndsType Local) | Re-align with `--alignEndsType EndToEnd` |
| 40-70% read loss in PAR-CLIP | T->C mismatches exceed 4% ceiling | Raise `--outFilterMismatchNoverReadLmax` to 0.07 |
| HISAT2 calls peaks STAR misses | HISAT2 lower-stringency soft-clip behaviour | Re-run HISAT2 with `--no-softclip`; differences should narrow to < 5% |
| Two STAR versions (2.4 vs 2.7) give different BAMs | Index sjdb format changed; parameter defaults drift | Pin STAR version for cross-study comparison; document |

**Operational rule:** For ENCODE-comparable analysis, use STAR 2.4.0j or 2.7.x with the ENCODE eCLIP parameter block; use `--alignEndsType EndToEnd`; use `--outFilterMultimapNmax 1` unless the biology requires multi-mappers (then add CLAM); UMI-dedupe with `umi_tools dedup --method=unique`. Document any deviation in methods.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `STAR ERROR ... Could not allocate memory` | Index too large for RAM | Switch to HISAT2; or rebuild STAR index with `--genomeSAsparseD 2` |
| `samtools sort: not enough memory` after STAR | sort -m default 768M too small | `samtools sort -m 4G -@ 8` |
| Crosslink density at uniform 5 nt offsets across genome | R2 5' trimmed by mistake | Recheck preprocessing; cutadapt `-g` and fastp `--trim_front2` are banned for CLIP |
| 90% of reads MAPQ < 10 | Genome index built for different species | Verify with `samtools view -h sample.bam \| head` against expected chromosome names |
| HISAT2 splice junctions miss-detected | Used `--no-splice` (HISAT2 splice OFF) | Remove `--no-splice` for mRNA studies |
| umi_tools dedup OOM | Too many UMIs at one position; deep library | Switch to `--method=unique` (less RAM than directional) |
| Multi-mapper count after STAR `--outFilterMultimapNmax 100` = 0 | `--outSAMmultNmax` not set | Add `--outSAMmultNmax -1` to emit all alignments |
| Reads MAPQ 0 dominate after bowtie2 | `--very-sensitive` ON but multi-mappers retained | Add `samtools view -q 30` post-filter |
| CLAM EM never converges | Background regions empty; sparse coverage | Lower `--max-tags` to 5; provide gene model GTF; verify multi-mapper BAM exists |

## References

- Van Nostrand EL et al 2016 Nat Methods 13:508 (eCLIP / ENCODE alignment parameter block)
- Konig J et al 2010 Nat Struct Mol Biol 17:909 (iCLIP truncation-as-CL principle)
- Dobin A et al 2013 Bioinformatics 29:15 (STAR aligner)
- Kim D et al 2019 Nat Biotechnol 37:907 (HISAT2)
- Langmead B & Salzberg SL 2012 Nat Methods 9:357 (bowtie2)
- Zhang Z & Xing Y 2017 Nucleic Acids Res 45:9260 (CLAM multi-mapper assignment)
- van de Geijn B et al 2015 Nat Methods 12:1061 (WASP allele-specific alignment)
- Hafner M et al 2010 Cell 141:129 (PAR-CLIP T->C mismatch tolerance need)
- Chakrabarti AM et al 2023 Genome Biol 24:235 (nf-core/clipseq alignment defaults)

## Related Skills

- clip-seq/clip-preprocessing - UMI extraction and adapter trimming before alignment
- clip-seq/clip-qc - Post-alignment QC (library complexity, read distribution, FRiP)
- clip-seq/crosslink-site-detection - Why the 5' end must be preserved
- clip-seq/clip-peak-calling - Downstream peak calling on the dedup BAM
- read-alignment/star-alignment - General STAR usage and indexing
- read-alignment/bowtie2-alignment - General bowtie2 usage
- alignment-files/duplicate-handling - Picard MarkDuplicates as UMI-less fallback
- single-cell/scatac-analysis - Cross-reference for single-cell variants
