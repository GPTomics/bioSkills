---
name: bio-long-read-sequencing-long-read-alignment
description: Aligns Oxford Nanopore and PacBio long reads (and assemblies) to a reference with minimap2 using the error-rate-matched preset (map-ont, lr:hq, map-hifi, map-pb, splice/splice:hq, asm5/10/20, ava), producing a sorted/indexed BAM for variant, SV, methylation, or isoform analysis. Covers why the preset rewrites the scoring/chaining model, why SV calling rides on supplementary not secondary alignments, carrying MM/ML methylation tags through with -y, the multi-part-index MAPQ trap, and when to swap in Winnowmap/VACmap/lra/pbmm2. Use when mapping ONT or PacBio reads, choosing a minimap2 preset by platform/chemistry, preparing input for Clair3/medaka/Sniffles/modkit, aligning into repeats/centromeres, or spliced-aligning cDNA/Iso-Seq.
tool_type: cli
primary_tool: minimap2
---

## Version Compatibility

Reference examples tested with: minimap2 2.28+, samtools 1.19+, winnowmap 2.03+, pbmm2 1.13+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Version-driven behavior to record:
- `lr:hq` and `map-iclr` were added in minimap2 2.27; `lr:hqae` in 2.28. Use >=2.28.
- `--MD` was broken by the 2.27 `--ds` addition and fixed in 2.28; use >=2.28 for any MD-dependent caller.
- A prebuilt `.mmi` index bakes in k/w/H/I - it must be built with the same preset used for alignment.

If code throws an error, introspect the installed tool (`minimap2 --help`, man page) and adapt the example to the actual API rather than retrying.

# Long-Read Alignment with minimap2

**"Align my long reads to the reference"** -> Map with the preset that matches the reads' ERROR RATE (not just platform), keeping the supplementary alignments and tags that downstream callers need.
- CLI: `minimap2 -ax lr:hq --MD -Y ref.fa reads.fq | samtools sort -o aln.bam` (accurate ONT/R10), `minimap2 -ax map-ont` (noisy R9 ONT), `minimap2 -ax map-hifi` (PacBio HiFi)

## The Single Most Important Modern Insight -- The Preset Rewrites the Scoring Model, So the Wrong One Fabricates or Erases Variants

`-x <preset>` is not a label. The man page defines each preset as a literal bundle that rewrites k-mer/window AND the entire scoring model (match `-A`, mismatch `-B`, gap-open `-O`, gap-extend `-E`), Z-drop `-z`, and chaining bandwidth `-r`. So the wrong preset does not merely "align worse" - it changes which gaps the chainer will span, and thereby fabricates or erases the exact insertions, deletions, introns, and SV breakpoints the downstream caller is built to find. Three corollaries an expert holds:

1. **Preset = read error rate, not platform.** "ONT" is no longer one regime: noisy R9/fast/hac = `map-ont`; accurate Q20+/duplex/R10-sup = `lr:hq` (2.27+, ~4x fewer CPU-hours, equal/better accuracy). `map-hifi` is literally `lr:hq` + HiFi scoring.
2. **Supplementary alignments ARE the SV signal.** SV callers read the split-read pattern (primary + supplementary chimeric pieces), not the tidy primary. Feeding them secondaries, or hard-clipping supplementaries, silently degrades SV sensitivity.
3. **A tag absent at alignment time is unrecoverable.** MM/ML, MD, cs - if minimap2 did not write them, no downstream tool can reconstruct them; the pipeline succeeds and produces empty/wrong results.

## Preset Taxonomy

| Preset | Read type / when correct | Notes |
|--------|--------------------------|-------|
| `map-ont` | ONT noisy genomic (R9, fast/hac) | the historic default; ~10% error scoring |
| `lr:hq` | accurate long reads <1% err (ONT Q20+/duplex/R10 sup) | 2.27+; the modern accurate-ONT default |
| `map-hifi` | PacBio HiFi/CCS genomic | = `lr:hq` + HiFi scoring (2.27+) |
| `map-pb` | PacBio CLR (legacy, ~15% err) | homopolymer-compressed minimizers; NEVER for HiFi |
| `splice` | noisy long RNA (ONT cDNA/direct RNA) | add `-uf` for stranded direct RNA |
| `splice:hq` | accurate long RNA (PacBio Iso-Seq, R10 cDNA) | |
| `asm5` / `asm10` / `asm20` | assembly-to-ref at ~0.1% / ~1% / ~5% divergence | PAF output; `--cs` for paftools call |
| `ava-ont` / `ava-pb` | all-vs-all read overlap (miniasm) | overlaps only, no base alignment |
| `lr:hqae` | accurate reads back to THEIR OWN assembly | 2.28+; fixes centromere self-mapping mismaps |

## Aligner Decision Tree

| Situation | Aligner | Why |
|-----------|---------|-----|
| Standard ONT/HiFi to a normal reference (SNV/SV/general) | minimap2 | the de-facto standard; default for Sniffles2, cuteSV, Clair3 |
| Accurate ONT (Q20+/duplex/R10 sup) | minimap2 `-x lr:hq` | ~4x faster than map-ont, equal/better |
| Centromeres / satellite arrays / segmental dups / T2T reference | Winnowmap2 | minimap2 minimizer-masking mismaps long tandem repeats; Winnowmap down-weights via meryl repetitive k-mers |
| Complex/nested SVs, inversions, tandem dups | VACmap (or lra) | variant-aware nonlinear chaining resolves CSVs minimap2 splits |
| Accurate reads -> a diploid assembly built from them | minimap2 `-x lr:hqae` (2.28+) | avoids self-assembly centromere mismaps |
| PacBio-native (.bam/.xml, want sorted+indexed in one call) | pbmm2 | minimap2 + PacBio plumbing; presets SUBREAD/CCS/HIFI/ISOSEQ |
| Legacy Sniffles1 reproduction | NGMLR | the 2018 standard, now superseded by minimap2+Sniffles2 |

## Tag Requirements by Downstream Tool

A missing tag is a silent failure. Add the tag at alignment time.

| Tag / flag | What it does | Needed for |
|------------|--------------|-----------|
| `--MD` | mismatch positions vs ref | many small-variant callers, IGV mismatch coloring (use minimap2 >=2.28) |
| `-Y` | soft-clip supplementary (default hard-clips) | SV callers: keeps breakpoint/insertion SEQ on the split read |
| `-y` | copy MM/ML (and other) tags from the input | methylation: carries Dorado MM/ML through alignment |
| `--cs` | minimap2 difference string | `paftools.js call` (assembly/long-read variant calling) requires it |
| `--eqx` | `=`/`X` CIGAR instead of `M` | tools that read match/mismatch from CIGAR |
| `-L` | move >65535-op CIGAR to CG:B tag | ultra-long ONT reads (else unrepresentable in BAM) |

Supplementary (flag 0x800) = split piece of one read across loci = the SV substrate, controlled by chaining + `-Y`. Secondary (flag 0x100) = multi-mapping alternative, controlled by `--secondary`/`-N`/`-p`. SV work keeps primary+supplementary and is fine with `--secondary=no`.

## Core Commands

```bash
# Accurate ONT (Q20+/R10 sup) -> genome, SV+variant ready, sorted+indexed
minimap2 -ax lr:hq -t 16 --MD -Y -R '@RG\tID:s1\tSM:s1' ref.fa reads.fq.gz \
  | samtools sort -@4 -o aln.bam && samtools index aln.bam

# Noisy ONT (R9 / fast / hac)
minimap2 -ax map-ont -t 16 --MD -Y ref.fa r9.fq.gz | samtools sort -o ont.bam

# PacBio HiFi (minimap2, or pbmm2 in one sorted+indexed call)
minimap2 -ax map-hifi -t 16 --MD -Y ref.fa hifi.fq.gz | samtools sort -o hifi.bam
pbmm2 align --preset HIFI --sort -j 16 ref.fa hifi.bam hifi.aligned.bam

# Methylation passthrough: carry Dorado MM/ML through alignment (the -y trap). -Y soft-clips
# supplementary records so hard-clipping does not break the MM per-base skip counting.
samtools fastq -T MM,ML dorado.mod.bam \
  | minimap2 -ax lr:hq -y -Y --MD ref.fa - \
  | samtools sort -o meth.bam        # then modkit pileup meth.bam ...

# Direct RNA (ONT): stranded forward-only, small k for terminal-exon sensitivity
minimap2 -ax splice -uf -k14 -G500k ref.fa dRNA.fq.gz | samtools sort -o drna.bam
#   -G500k raises max-intron above the 200k default only for genes with long introns

# Assembly-to-reference: PAF is correct here; --cs enables paftools variant calling
minimap2 -cx asm5 --cs ref.fa asm.fa > asm.paf
paftools.js call asm.paf > asm.var.vcf

# Repeats / centromeres / T2T: Winnowmap (precompute repetitive k-mers)
meryl count k=15 output merylDB ref.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
winnowmap -W repetitive_k15.txt -ax map-ont ref.fa reads.fq.gz | samtools sort -o wm.bam

# Prebuild index - bake the SAME preset's k/w in (else the preset's k/w is ignored)
minimap2 -x lr:hq -d ref.lrhq.mmi ref.fa
```

## Per-Method Failure Modes

### Hard-clipped supplementaries break SV insertion calls
**Trigger:** mapping for SV calling without `-Y`. **Mechanism:** minimap2 hard-clips supplementary records, discarding the breakpoint-spanning bases. **Symptom:** imprecise/missing insertions and translocations. **Fix:** add `-Y` (soft-clip) so split reads keep full SEQ.

### Methylation tags silently dropped
**Trigger:** aligning a Dorado mod BAM without preserving tags. **Mechanism:** `samtools fastq` strips MM/ML unless `-T MM,ML`; minimap2 ignores them unless `-y`. **Symptom:** aligned BAM has no MM/ML; modkit produces empty bedMethyl, no error. **Fix:** `samtools fastq -T MM,ML | minimap2 -y -Y`, or use `dorado aligner`.

### Multi-part index destroys MAPQ
**Trigger:** reference larger than `-I` (default 8G) - large plant/polyploid or concatenated refs. **Mechanism:** minimap2 builds a multi-part index and scores batches independently, so cross-batch best hits are invisible and MAPQ is wrong. **Symptom:** "no @SQ lines ... use --split-prefix"; spurious MAPQ. **Fix:** `-I <bigger-than-ref>` or `--split-prefix`.

### Wrong preset on accurate reads
**Trigger:** `map-ont` on Q20/R10/duplex, or `map-pb` on HiFi. **Mechanism:** noisy-read scoring on accurate reads (or CLR scoring on HiFi). **Symptom:** ~4x slower for no gain (map-ont case), or spurious clips/indels (map-pb-on-HiFi). **Fix:** `lr:hq` for accurate ONT, `map-hifi` for HiFi.

### Direct-RNA junctions on the wrong strand
**Trigger:** `-ax splice` on direct RNA without `-uf`. **Mechanism:** splice defaults to `-ub` (GT-AG on both strands), but dRNA is stranded. **Symptom:** invented/misplaced introns. **Fix:** add `-uf` (and usually `-k14`).

### Centromere/SD mismapping looks fine in flagstat
**Trigger:** plain minimap2 into long tandem repeats. **Mechanism:** minimizer masking collapses minimizer density, so reads map to the wrong paralog/copy. **Symptom:** reads still "map" (flagstat clean) but produce false SVs/heterozygosity in repeats. **Fix:** Winnowmap2 with a meryl repetitive-k-mer set.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `lr:hq` for reads <1% error | minimap2 2.27 NEWS / Li | accurate-read preset; ~4x fewer CPU-hours than map-ont |
| `-I 8G` default index batch | minimap2 man page | refs above it split into a MAPQ-breaking multi-part index |
| `distinct=0.9998` meryl k-mer cutoff | Winnowmap2 (Jain 2022) | flags the most-frequent k-mers to down-weight in repeats |
| `-G 200k` default max intron (splice) | minimap2 man page | raise only to the real longest intron; excess slows and invents alignments |
| minimap2 >= 2.28 | minimap2 NEWS | `lr:hq`/`lr:hqae` present and the 2.27 `--MD` regression fixed |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| SV insertions imprecise/missing | supplementaries hard-clipped | add `-Y` |
| modkit bedMethyl empty after alignment | MM/ML dropped | `samtools fastq -T MM,ML | minimap2 -y -Y` |
| "no @SQ lines ... use --split-prefix" | ref exceeds `-I`, multi-part index | `-I <bigger>` or `--split-prefix` |
| Preset k/w seems ignored | `.mmi` built with a different preset | rebuild index with the same `-x` preset |
| `paftools.js call` fails on PAF | missing base CIGAR / cs | `minimap2 -cx asm5 --cs` |
| Reads mismap in centromeres/SDs | minimizer masking | Winnowmap2 with meryl repetitive k-mers |
| Spurious wrong-strand introns (direct RNA) | `splice` default `-ub` | add `-uf` |

## References

- Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094-3100.
- Li H. 2021. New strategies to improve minimap2 alignment accuracy. *Bioinformatics* 37(23):4572-4574.
- Jain C, Rhie A, Hansen NF, et al. 2022. Long-read mapping to repetitive reference sequences using Winnowmap2. *Nat Methods* 19:705-710.
- Ren J, Chaisson MJP. 2021. lra: a long read aligner for sequences and contigs. *PLoS Comput Biol* 17(6):e1009078.
- Ding H, et al. 2026. VACmap: an accurate long-read aligner for unraveling complex genomic rearrangements. *Nat Commun* 16:11198.
- Sedlazeck FJ, et al. 2018. Accurate detection of complex structural variations using single-molecule sequencing (NGMLR/Sniffles). *Nat Methods* 15:461-468.

## Related Skills

- basecalling - The basecaller chemistry/error rate that picks the preset; carries MM/ML to pass with `-y`
- long-read-qc - Read length/quality before mapping; % identity from the aligned BAM
- structural-variants - Consumes the supplementary (split-read) signal this preserves with `-Y`
- clair3-variants - Small-variant calling on this BAM (needs the matched basecaller model)
- nanopore-methylation - Pileup of the MM/ML tags carried through with `-y`
- isoseq-analysis - Spliced alignment of full-length cDNA/Iso-Seq
- alignment-files/sam-bam-basics - Sort/index/inspect the BAM this produces
- alignment-files/alignment-filtering - Filter by MAPQ and secondary/supplementary flags
- genome-assembly/long-read-assembly - Assemble the reads instead of reference-mapping
