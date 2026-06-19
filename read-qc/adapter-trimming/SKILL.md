---
name: bio-read-qc-adapter-trimming
description: Removes sequencing adapters from FASTQ reads with Cutadapt and Trimmomatic, including paired-end read-through, small-RNA 3' adapters, amplicon primers, and anchored/linked adapters. Use when FastQC shows adapter content climbing toward the 3' end, when inserts are shorter than the read length (small-RNA, cfDNA, FFPE), or before assembly/k-mer analysis. For all-in-one trimming use fastp-workflow; for quality/length filtering use quality-filtering.
tool_type: cli
primary_tool: cutadapt
---

## Version Compatibility

Reference examples tested with: Cutadapt 4.4+, Trimmomatic 0.39+, fastp 0.23+, FastQC 0.12+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Adapter Trimming -- adapter content IS the insert-size distribution

Remove adapter sequence that the polymerase read INTO once it ran off the end of a short insert, using Cutadapt (precise, the correctness reference) or Trimmomatic (palindrome mode for paired read-through).

**"Trim adapters from my reads"** -> Detect and remove 3' adapter introduced by read-through, then length-filter the survivors.
- CLI: `cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 20 -o R1.fq -p R2.fq in_R1.fq in_R2.fq`
- All-in-one alternative: `fastp` (PE overlap analysis needs no adapter sequence) -> read-qc/fastp-workflow

Scope: this skill OWNS adapter/primer removal. Quality and length filtering -> read-qc/quality-filtering. Single-pass trim+QC -> read-qc/fastp-workflow. Contaminant/PhiX k-mer removal -> read-qc/contamination-screening. OUT OF SCOPE: quality-score trimming as a standalone goal (usually unnecessary before soft-clipping aligners; see insight 2).

## The Single Most Important Modern Insight

1. **Adapter appears only when the insert is shorter than the read, so adapter content is a direct readout of the insert-size distribution -- and adapter trimming is 3'-only for standard Illumina.** The library is `[P5]-[insert]-[P7]`; a read primes at the insert boundary and reads 5'->3' into the insert, running into the 3'/P7-side adapter only if it runs out of insert. Short-insert libraries (small-RNA ~22 nt, cfDNA ~167 bp, FFPE, degraded RNA, ancient DNA) are read-through-dominated; long-insert WGS may show almost none. The FastQC adapter-content curve climbing toward the 3' end IS that insert-size signal.

2. **Adapter trimming is the one near-universal preprocessing step; quality trimming usually is not.** Local aligners (BWA-MEM, STAR, Bowtie2 local, HISAT2) SOFT-CLIP low-quality tails, so quality trimming is redundant or harmful for alignment-based DNA/RNA (MacManes 2014, Williams 2016; GATK discourages it before BQSR). But aligners do NOT reliably remove ADAPTER -- adapter is foreign sequence with genuine base quality, so the aligner may try to align it and anchor a wrong placement. Trim adapter; leave quality trimming to the cases that need it (assembly, k-mer/pseudo-alignment, small-RNA, amplicon, no-BQSR variant calling).

3. **Small-RNA inverts the logic: the adapter is on EVERY read, so DISCARD reads with no adapter.** A ~22 nt miRNA insert is far shorter than a 50-75 nt read, so read-through is universal; a read with no detectable adapter is an adapter dimer, a too-long contaminant, or junk. Use `--discard-untrimmed` plus a tight length gate (`-m 18 -M 30`). This is the OPPOSITE of genomic DNA, where the no-adapter reads are the good full-length inserts.

Two-color trap: on NextSeq/NovaSeq, a high-quality poly-G tail is NOT adapter and is not removed by adapter trimming -- it needs a chemistry-aware poly-G trim (`cutadapt --nextseq-trim=20`, or fastp's auto poly-G). See read-qc/quality-reports and read-qc/fastp-workflow.

## Verified Adapter Sequences

| Kit | Read | Sequence |
|-----|------|----------|
| Illumina TruSeq | R1 3' | AGATCGGAAGAGCACACGTCTGAACTCCAGTCA |
| Illumina TruSeq | R2 3' | AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT |
| TruSeq (shared stem -- catches both) | -- | AGATCGGAAGAGC |
| Nextera / Tn5 | transposase | CTGTCTCTTATACACATCT |
| TruSeq small-RNA | 3' | TGGAATTCTCGGGTGCCAAGG |

The R1 3' adapter is the reverse complement of the R2-side region; trimming the shared 13 bp stem `AGATCGGAAGAGC` on both mates catches TruSeq read-through.

## Tool Taxonomy

| Tool | Mechanism | When it wins |
|------|-----------|--------------|
| Cutadapt | Error-tolerant semiglobal alignment of a supplied adapter | PRECISION: small-RNA 3' adapter, amplicon/16S primers, anchored/linked adapters, demultiplexing. The correctness reference. |
| Trimmomatic | ILLUMINACLIP simple + palindrome modes; ordered step pipeline | Legacy/reproducibility pipelines; palindrome PE read-through detection |
| fastp | PE overlap analysis (no adapter sequence needed) + auto poly-G | DEFAULT general-purpose trim; one fast pass (route OUT -> fastp-workflow) |
| Trim Galore | Cutadapt + FastQC wrapper with adapter auto-detect | Bisulfite/RRBS (`--rrbs`), Bismark pipelines |
| BBDuk | k-mer match against an adapter/contaminant reference | Contaminant/PhiX removal in the same pass (route OUT -> contamination-screening) |

## Decision Tree by Scenario

| Scenario | Use | Why |
|----------|-----|-----|
| General Illumina PE WGS/WES/RNA | fastp, or cutadapt with the TruSeq stem | Overlap analysis needs no sequence; cutadapt for explicit control |
| Small-RNA / miRNA | cutadapt `-a TGGAATTCTCGGGTGCCAAGG -m 18 -M 30 --discard-untrimmed` | Adapter on every read; gate length and drop no-adapter reads |
| Amplicon / 16S primers | cutadapt linked/anchored adapters | Primers are at fixed positions; needs precise placement |
| PE read-through, no adapter sequence known | fastp overlap, or Trimmomatic palindrome | Both detect read-through from the R1/R2 overlap |
| Bisulfite / RRBS | Trim Galore `--rrbs` | Handles MspI fill-in and Bismark conventions |
| NextSeq/NovaSeq with poly-G tails | fastp (auto) or cutadapt `--nextseq-trim` | Poly-G is high-Q; quality trim alone misses it |

Default when uncertain: fastp for bulk PE, cutadapt with the TruSeq stem for explicit single-tool control.

## Cutadapt

The algorithm is semiglobal (overlap) alignment, so a partial 3' adapter at the read end is detected. Two defaults drive behavior: `-e` (error rate, default 0.1) is computed against the LENGTH OF THE MATCHED REGION, not the whole adapter (an 8 bp match with 1 error is rate 0.125 and is rejected at the default); `-O` (minimum overlap, default 3) costs only ~0.07 bases lost per read by chance.

```bash
# Single-end 3' adapter
cutadapt -a AGATCGGAAGAGC -m 20 -o trimmed.fq.gz in.fq.gz

# Paired-end TruSeq (shared stem on both mates); both reads of a pair are discarded together
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 20:20 \
         -o R1.fq.gz -p R2.fq.gz in_R1.fq.gz in_R2.fq.gz

# Small-RNA: adapter on every read -> discard untrimmed, gate length
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 18 -M 30 --discard-untrimmed -j 8 \
         -o mirna.fq.gz raw.fq.gz

# Amplicon: linked 5'...3' primers (anchor with ^ to require the 5' primer)
cutadapt -g ^FWDPRIMER...REVPRIMER -o trimmed.fq.gz in.fq.gz

# 2-color poly-G aware (treats G as low quality so high-Q poly-G is trimmed)
cutadapt --nextseq-trim=20 -a AGATCGGAAGAGC -m 20 -o out.fq.gz in.fq.gz

# Higher error tolerance / longer required overlap when matches are missed / spurious
cutadapt -a ADAPTER -e 0.15 -O 5 -m 20 -o out.fq.gz in.fq.gz
```

Key flags: `-a/-g/-b` (3'/5'/anywhere, R1), `-A/-G/-B` (R2), `-q` (quality trim, BWA running-sum, runs BEFORE adapter removal), `--pair-filter {any,both,first}` (default any), `--max-n`, `--action {trim,mask,lowercase,none}`. When a filtering option discards reads in PE mode, both files MUST be processed together or they fall out of sync.

## Trimmomatic

`ILLUMINACLIP:<adapters.fa>:<seedMismatches>:<palindromeClip>:<simpleClip>:<minAdapterLen>:<keepBothReads>`

- seedMismatches: mismatches tolerated in the initial seed (commonly 2).
- palindromeClip (~30): log-odds threshold for the PE palindrome alignment; ~30 needs ~50 matched bases.
- simpleClip (~10): log-odds threshold for an adapter-vs-read match; ~10 needs ~16 bases.
- keepBothReads: DEFAULT False -- after palindrome detects read-through, R2 is redundant (reverse complement of R1) and is DROPPED; set True if a downstream tool needs both mates.

SIMPLE mode tests each adapter against each read. PALINDROME mode (PE-only) aligns R1+adapter against the reverse complement of R2+adapter, so it detects read-through even when only a few adapter bases remain or the adapter is entirely past the read end. Steps run in COMMAND-LINE ORDER; put ILLUMINACLIP first and MINLEN last so the length check reflects all prior trimming.

```bash
# Paired-end, palindrome-capable adapter file, MINLEN last
trimmomatic PE -phred33 -threads 8 \
    in_R1.fq.gz in_R2.fq.gz \
    R1_paired.fq.gz R1_unpaired.fq.gz R2_paired.fq.gz R2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads MINLEN:36

# Built-in adapter files ship with the install
ls $CONDA_PREFIX/share/trimmomatic-*/adapters/
```

PE mode emits FOUR files: paired (both mates survived) and unpaired/orphan (mate dropped). Feed the paired files to the aligner; the orphans stay synchronized out of the way.

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| FastQC still shows adapter after trimming | Wrong adapter, too-low `-e`, or only partial stem used | Use the shared stem AGATCGGAAGAGC; raise `-e` to 0.15; BLAST the overrepresented sequence |
| Reads truncated / many lose a few bp | `-O` too low -> random 3-mer matches | Raise `-O` (e.g. 5); the default loses ~0.07 bp/read by chance |
| Aligner reports R1/R2 out of sync | Mates trimmed/filtered independently | Process pairs together (cutadapt `-p`; Trimmomatic paired outputs) |
| Small-RNA yields huge "reads" | Forgot `--discard-untrimmed` / length gate | Add `--discard-untrimmed -m 18 -M 30` |
| 3' G-content rise persists after trimming | 2-color poly-G is high-quality, not adapter | `cutadapt --nextseq-trim` or fastp auto poly-G |
| Half of R2 disappears in Trimmomatic | keepBothReads default False drops redundant R2 | Add `keepBothReads` (True) if both mates are needed downstream |

## References

Martin M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17(1):10-12.
Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30(15):2114-2120.
MacManes MD. 2014. On the optimal trimming of high-throughput mRNA sequence data. Frontiers in Genetics 5:13.
Williams CR, Baccarella A, Parrish JZ, Kim CC. 2016. Trimming of sequence reads alters RNA-Seq gene expression estimates. BMC Bioinformatics 17:103.
Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17):i884-i890.

## Related Skills

read-qc/quality-reports - Read the adapter-content panel that triggers trimming
read-qc/quality-filtering - Quality and length filtering after adapter removal
read-qc/fastp-workflow - All-in-one adapter + quality trim with auto poly-G
read-qc/contamination-screening - k-mer removal of PhiX/vector/contaminant sequence
small-rna-seq/smrna-preprocessing - Full small-RNA adapter + length workflow
read-alignment/bwa-alignment - Soft-clipping aligner that handles low-quality tails without trimming
