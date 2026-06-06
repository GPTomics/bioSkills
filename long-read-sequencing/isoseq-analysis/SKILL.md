---
name: bio-long-read-sequencing-isoseq-analysis
description: Discovers, classifies, filters, and quantifies full-length transcript isoforms from PacBio Iso-Seq/Kinnex (HiFi) and Oxford Nanopore (cDNA/direct-RNA) long reads, using the isoseq+pigeon pipeline, SQANTI3, and ONT tools (IsoQuant, FLAIR, Bambu, StringTie2). Covers why a novel isoform is an artifact until proven otherwise (RT template-switching, intra-priming, and 5' degradation manufacture junctions and truncations), the SQANTI3 structural categories and their trust order, the Kinnex skera-split step, orthogonal CAGE/poly-A/short-read-junction validation, and why long-read isoform quantification needs EM. Use when building a full-length isoform catalog, classifying/filtering long-read transcripts, running Iso-Seq or ONT cDNA/dRNA analysis, or judging novel-isoform reliability.
tool_type: mixed
primary_tool: SQANTI3
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: isoseq 4.3+, pigeon 1.2+, SQANTI3 5.2+, pbmm2 1.13+, minimap2 2.28+, IsoQuant 3.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python/R: `pip show <pkg>` / `packageVersion('<pkg>')` for SQANTI3/IsoQuant/Bambu

Results depend on inputs that outlive the binary version - record them:
- The reference annotation + genome version drive SQANTI3/pigeon classification; record them.
- Orthogonal support files (CAGE refTSS BED, poly-A motif/atlas, short-read STAR SJ.out) determine which novels survive; record their provenance.
- The Iso-Seq binary was renamed `isoseq3` -> `isoseq` in v4; the classifier `pigeon` is a separate binary.

If code throws an error, introspect the installed tool (`isoseq --help`, `pigeon --help`, `sqanti3_qc.py --help`) and adapt the example to the actual API rather than retrying.

# Full-Length Isoform Analysis

**"Find the isoforms in my long-read RNA data"** -> Build a full-length isoform catalog, then classify and filter it against the reference with orthogonal end/junction support - because discovery without curation is a catalog of artifacts.
- CLI: `isoseq refine ... && isoseq collapse ... && pigeon classify ... && pigeon filter ...` (PacBio), `IsoQuant`/`FLAIR`/`Bambu` (ONT)

## The Single Most Important Modern Insight -- A Novel Isoform Is an Artifact Until Proven Otherwise

RT template-switching, intra-priming on genomic poly-A, and 5' RNA degradation actively MANUFACTURE novel junctions and truncated isoforms. So the classification + filter + orthogonal validation IS the analysis, not a QC postscript. Invert the posture from "I discovered N novel isoforms" to "I curated N novel isoforms that survived artifact filtering." Three consequences:

1. **A high novel-isoform fraction is a RED FLAG, not a success** - it usually means an under-powered filter or degraded RNA, not unusually rich biology.
2. **ISM (incomplete-splice-match) is the RNA-degradation thermometer, not a discovery.** ISMs are 5'-truncated FSMs; a high ISM fraction signals bad RNA integrity. Do not report ISMs as novel isoforms without CAGE 5' support.
3. **The orthogonal validation triad is mandatory:** CAGE peaks for the 5' TSS (catches 5' degradation), poly-A atlas/motif for the 3' TES (catches intra-priming), and short-read STAR junctions for splice sites (catches RT-switch/NNC junk).

## SQANTI3 Structural Categories (trust order)

Reference comparison is junction-chain based. NIC > NNC in trust, always; ISM is a diagnostic, not a discovery.

| Category (field value) | Meaning | Trust |
|------------------------|---------|-------|
| FSM (`full-splice_match`) | every internal junction matches a reference transcript; ends may differ | highest (known); ends still need CAGE/polyA |
| ISM (`incomplete-splice_match`) | junction subset of a reference (fewer 5' exons) | low - the 5'-degradation/RT-dropoff signature; trust only with CAGE |
| NIC (`novel_in_catalog`) | novel combination of KNOWN splice sites | high among novels - RT-switching cannot fake a NIC |
| NNC (`novel_not_in_catalog`) | >=1 genuinely novel splice site | lower - where junction artifacts concentrate; needs canonical/short-read support |
| genic / genic_intron | overlaps introns/exons; within an intron | low - pre-mRNA / gDNA carryover |
| fusion | spans >=2 genes | RT-chimera until proven by short-read split reads |
| intergenic / antisense | no gene overlap / antisense | novel-gene candidate or artifact; needs ORF/CAGE/conservation |

Mono-exon transcripts have no junctions to validate and are the false-discovery sink (intra-priming + gDNA run unchecked) - require ORF + CAGE + polyA + conservation before belief.

## Platform / Tool Decision Tree

| Data / goal | Tool | Why |
|-------------|------|-----|
| PacBio Iso-Seq/Kinnex, turnkey | isoseq + pigeon | native PacBio collapse + SQANTI-style classify/filter, SMRT Link integrated |
| Any long-read transcriptome, full curation | SQANTI3 | structural classification + ~50 QC descriptors + rules/ML filter + rescue; PacBio and ONT |
| ONT bulk discovery + quantification | IsoQuant | intron-graph; lowest novel FP rate among ONT tools |
| ONT, want built-in differential splicing | FLAIR | align -> correct junctions -> collapse -> diffSplice |
| Quantification with a precision knob | Bambu | NDR (novel discovery rate) calibrates precision; R/Bioconductor |
| Genome-guided assembly / hybrid short+long | StringTie2 `-L` (`--mix`) | fast long-read transcript assembly |
| ONT single-cell long-read isoforms | FLAMES | single-cell/spatial full-length isoforms |
| Differential isoform usage (DTU/DTE) | -> alternative-splicing | this skill yields the filtered set + counts and hands off |

## cDNA vs Direct-RNA and Spliced Alignment

PacBio Iso-Seq and ONT cDNA sequence reverse-transcribed cDNA (modifications erased; strand from primers); ONT direct-RNA sequences native RNA (true strand, poly-A length, modifications preserved, lower accuracy). Match the minimap2 preset to the chemistry:

```bash
minimap2 -ax splice ref.fa ont_cdna.fq        # ONT cDNA (orient first with pychopper)
minimap2 -ax splice -uf -k14 ref.fa drna.fq   # ONT direct RNA (stranded -> -uf, small k)
minimap2 -ax splice:hq -uf ref.fa hifi.fa     # PacBio HiFi (or pbmm2 --preset ISOSEQ)
```

`-uf` forces the forward transcript strand - correct for stranded dRNA/Iso-Seq, wrong for unoriented ONT PCR-cDNA (orient with pychopper first).

## PacBio Iso-Seq / Kinnex Pipeline

```bash
# 0. Kinnex (MAS-seq) ONLY: deconcatenate the array into segmented reads FIRST
skera split movie.hifi_reads.bam mas_adapters.fasta movie.segmented.bam   # skip for classic Iso-Seq

# 1. Remove cDNA primers; 2. produce FLNC (full-length non-chimeric)
lima movie.segmented.bam primers.fasta movie.fl.bam --isoseq --peek-guess
isoseq refine movie.fl.5p--3p.bam primers.fasta movie.flnc.bam --require-polya

# 3. cluster (reference-free) or skip and align FLNC directly; 4. map; 5. collapse to isoforms
isoseq cluster2 movie.flnc.bam clustered.bam                              # cluster2 scales to large sets
pbmm2 align --preset ISOSEQ --sort ref.fa clustered.bam mapped.bam
isoseq collapse --do-not-collapse-extra-5exons mapped.bam movie.flnc.bam collapsed.gff
#   collapsed.flnc_count.txt = FLNC molecules per isoform = the real DEPTH metric

# 6. classify + filter with pigeon (needs the collapsed.sorted.gff after prepare, NOT a BAM)
pigeon prepare collapsed.gff            # sorts the transcript GFF
pigeon prepare annotation.gtf ref.fa    # sorts the annotation -> annotation.sorted.gtf, indexes genome
pigeon classify collapsed.sorted.gff annotation.sorted.gtf ref.fa \
    --fl collapsed.flnc_count.txt --cage-peak cage.refTSS.bed --poly-a polyA.motif.list
pigeon filter collapsed_classification.txt --isoforms collapsed.sorted.gff
pigeon report --exclude-singletons collapsed_classification.filtered_lite_classification.txt saturation.txt
```

pigeon is PacBio's productized SQANTI3 (classify/filter, NOT a quantifier). Substitute SQANTI3 itself for the full descriptor set, ML filter, rescue module, and ONT support:

```bash
sqanti3_qc.py collapsed.gff annotation.gtf ref.fa --CAGE_peak cage.bed --polyA_motif_list polyA.txt \
    --short_reads short_reads_fofn.txt    # isoforms positional defaults to GTF/GFF; add --fasta for FASTA input
sqanti3_filter.py rules collapsed_classification.txt   # or: sqanti3_filter.py ml ...
```

## Per-Method Failure Modes

### Counting ISMs as novel isoforms
**Trigger:** reporting incomplete-splice-match transcripts as discoveries. **Mechanism:** 5' RNA degradation truncates FSMs into ISMs. **Symptom:** inflated novel/ISM fraction tracking RNA quality, not biology. **Fix:** treat ISM fraction as an integrity QC; keep ISMs only with CAGE 5' support.

### Intra-priming false 3' ends
**Trigger:** trusting 3' ends without poly-A validation. **Mechanism:** oligo-dT primes on a genomic internal A-stretch. **Symptom:** spurious short/mono-exon transcripts; `perc_A_downstream_TTS` >59%. **Fix:** SQANTI3/pigeon filter on downstream genomic A-content and poly-A motif; `--require-polya` alone does NOT catch this.

### Believing NNC novels without scrutiny
**Trigger:** treating NNC like NIC. **Mechanism:** novel splice sites are where RT template-switching and mapping artifacts land. **Symptom:** novel junctions absent from short-read data. **Fix:** require canonical junctions or short-read SJ coverage; prefer NIC.

### Feeding pigeon a BAM
**Trigger:** `pigeon classify mapped.bam ...`. **Mechanism:** pigeon classifies the collapsed.sorted.gff after `pigeon prepare`, not an alignment. **Symptom:** wrong-input error. **Fix:** `isoseq collapse` -> `pigeon prepare` -> `pigeon classify`.

### Comparing isoform counts across libraries of different depth
**Trigger:** raw isoform counts as abundance. **Mechanism:** discovery is depth-unsaturated; truncated reads are multi-isoform-compatible. **Symptom:** deeper libraries "have more isoforms"; double-counted abundance. **Fix:** rarefaction curve (`--exclude-singletons`); EM quantification (Bambu/IsoQuant/NanoCount), not raw FLNC counts.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `perc_A_downstream_TTS` > 59-60% = intra-priming | SQANTI (Tardaguila 2018) | genomic A-rich window means the poly-A was internal, not the real tail |
| novel junction trusted if canonical OR short-read cov >= 3 | SQANTI3 rules filter | a single criterion for RT-switch/NNC artifacts |
| ML filter needs >= 250 Reference-Match FSM | SQANTI3 | enough true-positive labels to train; else falls back to rules |
| exclude singletons (1-FLNC) for saturation | pigeon report | singletons are the dominant unreliable novel bucket |
| FLNC count = depth metric | isoseq collapse | independently sequenced full-length molecules, before clustering/dedup |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `isoseq3: command not found` | renamed in v4 | use `isoseq` (subcommands unchanged) |
| pigeon classify wrong input | fed a BAM | give the collapsed.sorted.gff after `pigeon prepare` |
| Huge novel-isoform count | filter skipped/underpowered | run pigeon/SQANTI3 filter with CAGE/polyA/short-read support |
| Many mono-exon novels | intra-priming / gDNA carryover | filter on poly-A; require ORF/CAGE for mono-exon |
| Wrong-strand spliced alignment | `-uf` on unoriented cDNA | orient with pychopper, or drop `-uf` for cDNA |
| Isoform counts not comparable across samples | depth-unsaturated discovery | EM quantification + rarefaction curve |

## References

- Tardaguila M, de la Fuente L, Marti C, et al. 2018. SQANTI: extensive characterization of long-read transcript sequences for quality control in full-length transcriptome identification and quantification. *Genome Res* 28(3):396-411.
- Pardo-Palacios FJ, Arzalluz-Luque A, Kondratova L, et al. 2024. SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms. *Nat Methods* 21(5):793-797.
- Prjibelski AD, Mikheenko A, Joglekar A, et al. 2023. Accurate isoform discovery with IsoQuant using long reads. *Nat Biotechnol* 41(7):915-918.
- Tang AD, Soulette CM, van Baren MJ, et al. 2020. Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia (FLAIR). *Nat Commun* 11:1438.
- Chen Y, Sim A, Wan YK, et al. 2023. Context-aware transcript quantification from long-read RNA-seq data with Bambu. *Nat Methods* 20(8):1187-1195.
- Al'Khafaji AM, Smith JT, Garimella KV, et al. 2024. High-throughput RNA isoform sequencing using programmed cDNA concatenation (MAS-ISO-seq/Kinnex). *Nat Biotechnol* 42(4):582-586.

## Related Skills

- long-read-alignment - Spliced alignment of cDNA/direct-RNA (splice/splice:hq, `-uf`)
- basecalling - Direct-RNA (RNA004) basecalling; cDNA vs direct-RNA chemistry
- nanopore-methylation - Direct-RNA modifications are separate from isoform structure
- alternative-splicing/long-read-splicing - Long-read splicing analysis (define the boundary)
- alternative-splicing/isoform-switching - Differential isoform usage (DTU) downstream
- alternative-splicing/differential-splicing - Differential splicing downstream
- rna-quantification/tximport-workflow - Transcript-level quantification downstream
- genome-annotation/eukaryotic-gene-prediction - Long-read isoforms as annotation evidence
