---
name: bio-alignment-amplicon-clipping
description: Trim PCR primers from aligned reads in amplicon-panel BAMs using samtools ampliconclip. Use when processing SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA hot-spot panels, or any amplicon assay where primer-derived bases would falsely confirm reference at primer footprints.
tool_type: cli
primary_tool: samtools
---

## Version Compatibility

Reference examples tested with: samtools 1.19+, pysam 0.22+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Amplicon Clipping

**"Trim primer-derived bases from amplicon BAM"** -> Soft- or hard-clip the 5' primer footprint after alignment using a primer BED, then repair fixmate/MD/NM tags.
- CLI: `samtools ampliconclip -b primers.bed input.bam -o clipped.bam` (since samtools 1.11)
- Alternative: `iVar trim`, `BAMClipper`, `fgbio ClipBam`

## Why Primer Trimming After Alignment

Amplicon panels (SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA hot-spot panels, fusion panels, 16S rRNA) use designed PCR primers for enrichment. Primer-derived bases at read 5' ends do NOT reflect biological sequence -- they reflect the primer template. Without trimming:
- False reference confirmation at primer footprint positions.
- Variant allele frequency suppressed at variants under primers (the primer sequence cannot capture the variant base).
- Strand bias artifacts (primers are typically one-strand).

Standard amplicon BAMs should NEVER be processed by `samtools markdup` -- by design every read at a primer location is a "duplicate" by coordinate. See duplicate-handling for the assay-aware decision.

## Tool Selection

| Tool | When | Notes |
|------|------|-------|
| `samtools ampliconclip` | Default for amplicon panels (since 1.11) | Soft- or hard-clip from BED; modifies CIGAR; invalidates MD/NM |
| `iVar trim` | SARS-CoV-2 ARTIC pipelines | Coordinates by primer name/position; applies hard or soft clip |
| `BAMClipper` | Capture / hybrid panels with primer overlap | 5'-end clipping with overlap handling |
| `fgbio ClipBam` | When read-pair coordination matters | Soft/hard-clip with mate-aware end adjustment |
| `cutadapt` (pre-alignment) | Legacy / when alignment is downstream | Trims at FASTQ stage; less precise for amplicon |

## Soft-Clip vs Hard-Clip

**Goal:** Decide whether trimmed bases are kept in the BAM (reversible) or discarded (irreversible).

**Approach:** Soft-clip is the safe default; hard-clip only when archiving and disk is constrained.

| Mode | Flag | What it does | Reversible? |
|------|------|--------------|-------------|
| Soft-clip | (default) / `--soft-clip` | Bases kept in SEQ; CIGAR uses `S`; bases not aligned | Yes (CIGAR can be re-extended) |
| Hard-clip | `--hard-clip` | Bases removed from SEQ; CIGAR uses `H` | **No** (bases lost) |

Soft-clip is the recommended default. Hard-clip is irreversible -- once applied, the trimmed bases cannot be recovered for re-analysis with different primer coordinates.

## Basic ampliconclip Workflow

**Goal:** Trim primers from a coordinate-sorted, indexed amplicon BAM and produce a downstream-ready BAM.

**Approach:** Run `ampliconclip` with primer BED, optionally `--both-ends` and `--strand`, then re-fixmate (CIGAR changed) and re-calmd (MD/NM tags invalidated by clip).

```bash
# 1. Soft-clip primers (default; reversible)
samtools ampliconclip --both-ends --strand --soft-clip \
    -b primers.bed input.bam -o clipped.bam

# 2. Re-pair tags (CIGARs changed -- mate info needs refresh)
samtools sort -n clipped.bam | \
    samtools fixmate -m - - | \
    samtools sort -o sorted.bam -

# 3. Repair MD/NM tags (invalidated by clip; required by mpileup BAQ and IGV)
samtools calmd -b sorted.bam reference.fa > clipped_final.bam
samtools index clipped_final.bam
```

### Strand-Aware Clipping

`--strand` clips primer bases only on the strand the primer is designed for. Without `--strand`, both strands are clipped at the primer site, removing valid biological sequence on the off-strand.

### Both-End Clipping

`--both-ends` allows clipping at both 5' and 3' positions of the read (some primers can appear at either end after alignment). Necessary for amplicon designs where reads can read through the entire amplicon.

## Primer BED Format

```
# tab-separated, 0-based half-open like all BED
chr1   100   125   primer_1_F    +
chr1   500   525   primer_1_R    -
chr1   600   625   primer_2_F    +
chr1   1000  1025  primer_2_R    -
```

Tools that consume the BED: column 1-3 (region), column 6 (strand) is required for `--strand`. ARTIC primer schemes ship pre-built BEDs (e.g., `primer.bed` from artic-network/primer-schemes).

## SARS-CoV-2 ARTIC Comparison

| Tool | Approach | When |
|------|----------|------|
| `samtools ampliconclip` | Soft-clip from BED, post-alignment | nf-core/viralrecon, modern ARTIC workflows |
| `iVar trim` | Soft- or hard-clip with primer-position parsing | Original ARTIC field bioinformatics |

Modern viral consensus pipelines tend to use ampliconclip then `samtools consensus --config illumina --ambig` for IUPAC heterozygote handling. See reference-operations for consensus generation.

## After Clipping: Required Re-Processing

Clipping invalidates several tags and CIGAR-derived fields:

| Field | Impact | Repair |
|-------|--------|--------|
| CIGAR | New `S` or `H` operations added | Automatic from ampliconclip |
| MD:Z | Mismatch positions now wrong | `samtools calmd -b in.bam ref.fa` |
| NM:i | Edit distance recomputed | `samtools calmd` |
| TLEN | Template length changes when both mates clipped | `samtools fixmate -m` |
| MS, MC | Mate score / CIGAR | `samtools fixmate -m` |

A clipped BAM that bypasses fixmate + calmd causes silent failures in `bcftools mpileup` BAQ (which depends on MD), IGV mismatch coloring, and any tool using NM for filtering.

## Why Not Markdup

Amplicon reads at primer locations are by design coordinate-degenerate -- every read mapped to the same amplicon shares the same start/end coordinates because they all come from the same primer pair. `samtools markdup` would mark essentially every read as a duplicate and erase the dataset. For amplicon panels:

- WITHOUT UMIs: skip dedup entirely; rely on coverage uniformity from amplicon design.
- WITH UMIs (deep panels): use `fgbio GroupReadsByUmi` -> `CallMolecularConsensusReads` instead of markdup. See duplicate-handling.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `MD tag mismatch` after clipping | calmd not run | Run `samtools calmd -b clipped.bam ref.fa` |
| Variant calls with strand bias at every amplicon end | Forgot `--strand` | Re-run with strand-aware clipping |
| Markdup output shows ~100% duplicates | Amplicon BAM was processed with markdup | Restart from raw alignment; use ampliconclip; skip markdup |
| Unexpected reference confirmation at primer-overlapping variants | ampliconclip not run | Run before variant calling |

## Quick Reference

| Task | Command |
|------|---------|
| Soft-clip primers | `samtools ampliconclip --both-ends --strand -b primers.bed in.bam -o clipped.bam` |
| Hard-clip (irreversible) | `samtools ampliconclip --both-ends --strand --hard-clip -b primers.bed in.bam -o clipped.bam` |
| Repair MD/NM after clip | `samtools calmd -b clipped.bam ref.fa > final.bam` |
| Repair mate info | `samtools sort -n - \| samtools fixmate -m - - \| samtools sort -o out.bam -` |

## Related Skills

- duplicate-handling - Why amplicon BAMs should not be markdup'd; UMI-aware alternatives
- alignment-filtering - Post-clip filtering for amplicon variant calling
- alignment-sorting - Re-sort after fixmate
- pileup-generation - mpileup flags for amplicon (`-aa -A -d 600000 -B`)
- reference-operations - Consensus generation from amplicon BAMs (samtools consensus)
- read-qc/quality-reports - Pre-alignment adapter/quality trimming
