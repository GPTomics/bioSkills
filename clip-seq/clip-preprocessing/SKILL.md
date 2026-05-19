---
name: bio-clip-seq-clip-preprocessing
description: Preprocess CLIP-seq reads (eCLIP, iCLIP, iCLIP2, iCLIP3, irCLIP, PAR-CLIP, FLASH) with protocol-specific UMI extraction, adapter trimming, length filtering, and post-alignment PCR-duplicate collapse. Use when raw CLIP FASTQ must be turned into deduplicated, crosslink-preserving BAM input for peak calling; choosing between two-pass and single-pass adapter trimming; deciding minimum read length; or mapping UMI patterns to specific eCLIP/iCLIP/iCLIP2/iCLIP3 library preps.
tool_type: cli
primary_tool: umi_tools
---

## Version Compatibility

Reference examples tested with: umi_tools 1.1.5+, cutadapt 4.6+, fastp 0.23.4+, samtools 1.19+, pysam 0.22+, picard 3.1+, preseq 3.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# CLIP-seq Preprocessing

**"Preprocess raw CLIP reads into UMI-deduplicated, alignable FASTQ"** -> Extract random barcodes, trim adapters without disturbing the 5' truncation site, length-filter to remove unmappable shorts, and (post-alignment) collapse PCR duplicates by UMI + position. The 5' end of the read carries the iCLIP/eCLIP truncation signature one base downstream of the protein-RNA crosslink; preserving this base is the single most important constraint of CLIP preprocessing.

- CLI (eCLIP, paired-end): `umi_tools extract --bc-pattern=NNNNNNNNNN --stdin R1.fq.gz --read2-in R2.fq.gz --stdout R1_umi.fq.gz --read2-out R2_umi.fq.gz`
- CLI (iCLIP/iCLIP2, single-end): `umi_tools extract --bc-pattern=NNNXXXXNN --extract-method=string --stdin R1.fq.gz --stdout R1_umi.fq.gz` (3+2 random Ns flanking a 4 nt library barcode; demultiplex by the X positions first if multiplexed)
- CLI (PAR-CLIP): `umi_tools extract --bc-pattern=NNNN ...` (most protocols use 4 nt random barcodes; verify the lab's exact prep)
- CLI (3' trim only, eCLIP convention): `cutadapt -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --quality-base 33 --quality-cutoff 6 -m 18 -o R1.trim.fq.gz -p R2.trim.fq.gz R1_umi.fq.gz R2_umi.fq.gz`

The eCLIP convention is: trim quality and adapter from the 3' end ONLY. The 5' end of read 2 (in paired-end eCLIP) is the truncation site of the RT enzyme at the protein-RNA adduct, located one nucleotide downstream of the crosslink. Trimming the 5' end discards that exact base. `cutadapt -g`, `fastp --trim_front1`, and aggressive quality trimming of the 5' end are all banned for CLIP unless a documented protocol-specific reason exists.

## Read Structure by Protocol

| Protocol | UMI length and location | Truncation-site read | Adapter set | Notes |
|----------|-------------------------|----------------------|-------------|-------|
| eCLIP (Van Nostrand 2016 / ENCODE) | 10 nt at R1 5' end (random nucleotides preceding insert) | R2 5' end (after R1 UMI is stripped) = crosslink -1 | Illumina TruSeq R1 + R2 + inline X1A/X1B inverted adapter for two-pass | ENCODE pipeline normative |
| seCLIP (single-end eCLIP) | 10 nt at R1 5' end | R1 5' end | TruSeq R1 only | ENCODE accepts both eCLIP and seCLIP |
| iCLIP (Konig 2010) | 5 nt random (NNNXXXXNN: 3 N + 4 X library barcode + 2 N), single-end | R1 5' end after barcode strip | L3 adapter at 3' | Multiplexed - demultiplex by the 4 X bases |
| iCLIP2 (Buchbender 2020) | 5 or 9 nt random (NNNXXXXNN or longer), single-end | R1 5' end | L3 adapter at 3' | Increased complexity vs iCLIP; same UMI pattern |
| iCLIP3 (Buchbender et al, bioRxiv 2026.03.01.708747) | 10 nt random + dual sample index, single-end | R1 5' end | TruSeq | Silica-column RNA isolation; non-radioactive; streamlined low-input protocol. Preprint - verify final published version before pinning a pipeline. |
| irCLIP (Zarnegar 2016) | 5 nt random + barcode (similar to iCLIP) | R1 5' end | IR700 adapter + standard TruSeq sequencing adapter | Infrared replaces 32P; otherwise iCLIP-like |
| PAR-CLIP (Hafner 2010) | 0-4 nt depending on prep | T->C transitions within reads (NOT a truncation method) | TruSeq R1 | UMI optional; rely on T->C signature for CL |
| FLASH (Aktas 2020) | Sample barcode + UMI in custom adapter | R1 5' | Custom L3 design | 1.5 day protocol; adapter design proprietary to MPI |
| miCLIP / miCLIP2 (Linder 2015 / Kortel 2021) | iCLIP-style barcodes | R1 5' = m6A -1 (truncation OR C-to-T) | iCLIP L3 | m6A-specific |
| STAMP (Brannan 2021) | NA (no UV) | NA (C-to-U editing) | 10x or bulk RNA-seq adapters | Antibody-free editing-based; preprocess as RNA-seq |

If the read structure differs from the lab's documentation, run `seqkit head -n 1000 R1.fq.gz | seqkit stats -a` and inspect the first 12 bases of 100 random reads. Random-barcode positions show ~25% base composition per position; library/sample barcodes will be fixed across reads.

## Critical Choice: One Adapter Pass vs Two

**One pass (single-end iCLIP, PAR-CLIP):** Single 3' adapter; cutadapt with `-a <L3>` and `-m 18` is sufficient.

**Two passes (eCLIP, paired-end):** The eCLIP library design ligates an inverted X1A/X1B inline adapter that can appear at either end of a short fragment due to read-through. The ENCODE pipeline does pass 1 with the standard 3' adapter, then pass 2 trims a residual 5' adapter that read-through events leave on read 2. Trimming a 5' adapter from R2 is a special case: cutadapt's `-G <ADAPTER>` (uppercase) anchored at R2's 5' end is the right invocation. Do NOT use `-g` (lowercase) for R1 5' trimming - that destroys the truncation site.

**Cutadapt full eCLIP-style invocation:**

```bash
# Pass 1 - 3' adapter (both reads)
cutadapt \
    -a AGATCGGAAGAGCACACGTCT \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --quality-base 33 -q 6 \
    -m 18 \
    -o R1.p1.fq.gz -p R2.p1.fq.gz \
    R1.umi.fq.gz R2.umi.fq.gz

# Pass 2 - read-through 5' adapter on R2 only (ENCODE eCLIP)
cutadapt \
    -G GATCGTCGGACTGTAGAACTCTGAAC \
    --quality-base 33 -q 6 \
    -m 18 \
    -o R1.p2.fq.gz -p R2.p2.fq.gz \
    R1.p1.fq.gz R2.p1.fq.gz
```

The `-q 6` is intentionally permissive. Aggressive `-q 20` quality trimming chews back the 5' end of R2 and breaks truncation-based crosslink-site detection. The `-m 18` minimum-length cutoff is non-negotiable: reads shorter than 18 nt are functionally unmappable (multi-map prevalence > 50%, see CIMS analysis discussions in the CLIP review literature) and must be discarded before alignment.

## Per-Protocol Failure Modes

### eCLIP -- 5' end trimmed by mistake

**Trigger:** Pipeline written for generic RNA-seq applied to eCLIP; `fastp` defaults trim both ends; `cutadapt -g` invoked on R2.

**Mechanism:** The 5' end of R2 in eCLIP is the truncation site = crosslink -1. Quality-trim from 5' or untargeted 5' adapter trim discards that exact base.

**Symptom:** Downstream PureCLIP / iCount / CTK CITS analysis fails to call crosslink sites; peak calling still works but the peaks lose nucleotide precision.

**Fix:** Use `-q 6` (3' only by default in cutadapt) and never `--trim_front2` in fastp. Validate by inspecting BAM read 2 5' positions: 60-90% of unique R2 5' positions should map within 100 nt windows around known RBP binding motifs, not be uniformly distributed.

### iCLIP / iCLIP2 -- Demultiplex confusion with UMI

**Trigger:** Multiplexed iCLIP library; user runs `umi_tools extract` with `NNNNNNNNN` (9 N) when the actual prep is `NNNXXXXNN` (5 random + 4 sample barcode).

**Mechanism:** umi_tools treats the entire 9-base prefix as UMI, losing the sample identity in the middle 4 bases. Reads from different samples are merged.

**Symptom:** Library complexity inflated artificially; per-sample read counts seem high but binding profiles look averaged across samples; sample-specific motifs absent.

**Fix:** Demultiplex BEFORE UMI extraction. Use `umi_tools extract --bc-pattern=NNNXXXXNN --extract-method=string --filter-cell-barcode` with a whitelist of the 4 X-base barcodes; or split the FASTQ with `je demultiplex` against the library barcode table first, then `umi_tools extract --bc-pattern=NNNNN` (the 5 surviving random Ns).

### PAR-CLIP -- T->C mistaken for sequencing error

**Trigger:** Standard variant-calling pipeline applied to PAR-CLIP without recognizing the T->C signature.

**Mechanism:** PAR-CLIP's diagnostic mutation is T->C (4SU adduct pairs with G during RT). Upon crosslinking the per-position T->C rate jumps from ~0.5% baseline to 20-50% (see PAR-CLIP literature; Hafner 2010, Spitzer 2014). Pipelines that filter "high-error" reads or apply STAR `--outFilterMismatchNoverLmax 0.04` (4% mismatch ceiling) discard the very reads carrying the signal.

**Symptom:** Loss of 40-70% of PAR-CLIP reads; downstream T->C site calling (PARalyzer, wavClusteR) finds almost nothing.

**Fix:** For PAR-CLIP only, raise the STAR mismatch ceiling to `--outFilterMismatchNoverLmax 0.07` (downstream in clip-alignment); also raise quality-trim tolerance in cutadapt to `-q 6` (already the CLIP default but worth re-confirming for PAR-CLIP). Track T->C rate per nucleotide position with `samtools mpileup` to confirm the signature is preserved post-alignment. **Critical exception summary:** iCLIP/eCLIP/iCLIP2 keep 0.04; PAR-CLIP only raises to 0.07. See clip-seq/clip-alignment for the alignment-stage override.

### Quality trimming too aggressive

**Trigger:** Inherited fastp/Trimmomatic pipeline with `-q 20` or `MINLEN 36`.

**Mechanism:** CLIP fragments are short (20-75 nt insert) and the 3' end carries adapter. Aggressive quality trimming and length filtering destroy 30-60% of usable reads.

**Symptom:** Per-replicate unique-mapped read count < 500k from a ~30M raw library; library complexity calculation crashes from too few reads.

**Fix:** Use `-q 6 -m 18` (cutadapt) or `--qualified_quality_phred 6 --length_required 18` (fastp). Compare retained fraction: a properly preprocessed CLIP library retains ~70-85% of raw reads through trimming.

### Library complexity below threshold

**Trigger:** PCR duplication rate >> 50%; deduplicated unique fragment count < 1M.

**Mechanism:** Low input cells (< 5M), over-amplification (> 25 PCR cycles), or failed IP all produce libraries dominated by a small number of PCR-amplified molecules.

**Symptom:** `preseq lc_extrap` predicts plateau well below 5M unique reads at infinite sequencing depth; `picard CollectLibraryComplexity` reports `ESTIMATED_LIBRARY_SIZE` < 1M; UMI families have median size > 8.

**Fix:** No analytic rescue. Re-prep the library with more input cells and fewer PCR cycles (target 14-18 cycles for eCLIP, 16-20 for iCLIP2). If the dataset must be salvaged, downsample to the unique-fragment fraction and acknowledge the loss of statistical power in differential analyses.

## UMI Deduplication Decision

After alignment, collapse PCR duplicates by `(UMI, position, strand)`. The choice between `unique` and `directional` methods is a precision/recall tradeoff:

| Method | Match rule | Behaviour | When to use |
|--------|-----------|-----------|-------------|
| `--method=unique` | Exact UMI match | Strictest; two UMIs differing by 1 base treated as independent molecules | ENCODE convention; reproducibility against published peaks |
| `--method=directional` | Network of UMIs differing by hamming-1; pick most-abundant | Collapses UMI sequencing errors; slightly fewer unique fragments reported | Highest precision; preferred when UMI sequencing error rate > 1% |
| `--method=cluster` | Connected components within edit distance | Most aggressive collapse | Default umi_tools behaviour pre-2017; not recommended for CLIP |
| `--method=adjacency` | Adjacency clustering | Between unique and directional | Rarely used for CLIP |

For eCLIP and iCLIP, the **ENCODE convention is `--method=unique`** (Van Nostrand 2016); the Yeo lab pipeline uses this exclusively. Directional adds 5-15 minutes runtime for human eCLIP at typical depth and is more conservative on rare UMI sequencing errors but produces slightly different absolute counts.

```bash
samtools index aligned.bam
umi_tools dedup \
    --stdin=aligned.bam \
    --stdout=dedup.bam \
    --method=unique \
    --paired \
    --log=dedup.log
samtools index dedup.bam
```

Without UMIs (e.g., older eCLIP with only sample barcodes, some custom protocols), fall back to `picard MarkDuplicates --REMOVE_DUPLICATES true`. The trade-off: position-only dedup over-collapses (genuinely independent fragments that share start positions are lost) at ~5-15% in deep CLIP libraries; UMI dedup recovers them.

## Pre-Mapping rRNA Filter (Optional but Standard)

eCLIP libraries are 5-30% rRNA reads even after polyA depletion. ENCODE's pipeline pre-maps to a rRNA + RepeatMasked repeat index with bowtie2, then aligns unmapped reads to the genome. This is purely a performance optimization (avoiding STAR's multi-mapper tangle on rRNA) and does not change which reads survive deduplication; it only reorders the alignment stages.

```bash
# Pre-map to rRNA + repeats
bowtie2 -x repbase_repeats -U R1.trim.fq.gz \
    --un-gz R1.norep.fq.gz \
    -p 8 -S /dev/null
# Then align unmapped reads to genome with STAR
```

For non-eCLIP protocols (iCLIP, PAR-CLIP) the rRNA pre-map is optional; the alternative is to filter rRNA-overlapping reads after STAR alignment with `bedtools intersect -v -a aligned.bam -b rRNA.bed`.

## Library Complexity Assessment

```bash
# preseq lc_extrap predicts unique fragments at deeper sequencing
preseq lc_extrap -B -P aligned.bam -o complexity.txt
# Output: TOTAL_READS, EXPECTED_DISTINCT, LOWER_0.95CI, UPPER_0.95CI
# At 100M reads, EXPECTED_DISTINCT >= 10M = good complexity; < 3M = library failed

# picard direct estimate
picard EstimateLibraryComplexity \
    I=aligned.bam \
    O=picard_complexity.txt
# ESTIMATED_LIBRARY_SIZE > 5M = healthy
```

ENCODE eCLIP requires >= 1M unique fragments per replicate (after UMI dedup). Anything below 500k is functionally unusable for genome-wide peak calling. Between 500k and 1M, restrict analysis to high-expression transcripts.

## QC Checkpoints After Preprocessing

| Metric | Target | If below target |
|--------|--------|-----------------|
| Adapter trim retention | >= 70% reads | Adapter pattern wrong; recheck library prep docs |
| Mean read length post-trim | >= 25 nt | RNA fragmentation too aggressive in IP step |
| % reads >= 18 nt | >= 80% | Drop the prep; libraries with > 30% < 18 nt indicate degraded RNA |
| UMI-extract success rate | >= 95% | UMI pattern wrong; reads do not start with random Ns |
| PCR duplication rate | 30-70% (CLIP normal) | < 30% means under-sequenced; > 90% means over-amplified |
| Library complexity (preseq) | >= 1M unique at sequenced depth | Library failed; cannot rescue analytically |

CLIP libraries have HIGH duplication rates (40-70%) by design - the IP enriches a small pool of molecules. This is NOT a problem if UMIs collapse duplicates correctly. A "30% duplication rate" CLIP library typically means either (a) the IP failed (no enrichment, so library looks like RNA-seq) or (b) the library was undersequenced and unique molecules dominate.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| umi_tools extract: "Read does not match pattern" for > 10% of reads | UMI pattern wrong (5 vs 9 vs 10 nt) | Inspect first 12 bases of 100 reads; rerun extract with correct `--bc-pattern` |
| cutadapt: 95% reads "Too short, filtered" | Adapter sequence wrong; reads are adapter-only after trim | Verify adapter sequence in library prep documentation; check both R1 and R2 adapters |
| All reads aligning to chrM after preprocessing | Pre-map to rRNA skipped; rRNA reads dominate | Add bowtie2 rRNA pre-map OR `samtools view -F 4 -L exclude_chrM_rRNA.bed` post-align |
| umi_tools dedup very slow (> 4h on 30M BAM) | Default `--method=directional` on dense libraries | Switch to `--method=unique` (ENCODE convention) |
| eCLIP truncation positions look uniform across genome | 5' adapter trim destroyed R2 5' end | Re-preprocess; `-g` should never touch R2 5' in eCLIP |
| 30% T->C mismatches in PAR-CLIP fail mapping | STAR mismatch ceiling 0.04 too strict | Raise to `--outFilterMismatchNoverLmax 0.07` for PAR-CLIP only |
| iCLIP2 reads after dedup << expected unique | Demultiplex done with whole UMI; samples merged | Demultiplex by 4 X bases of NNNXXXXNN first, then extract UMI |

## Reconciliation: Preprocessing Tools

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| fastp gives more reads than cutadapt | fastp default polyG trimming kept short reads; cutadapt min-length filtered them | Both correct; pick one and document |
| umi_tools and je-suite give different dedup counts | Different UMI distance metric (unique vs hamming-1 vs network) | Use ENCODE convention `umi_tools --method=unique` for cross-study comparability |
| Library complexity (preseq) vs (picard) disagree | preseq extrapolates non-linearly; picard estimates at sequenced depth only | Trust preseq at sequenced depth, picard for absolute library size |
| Read count after Yeo lab pipeline >> nf-core/clipseq | Yeo includes rRNA reads; nf-core pre-filters | Both correct; downstream peak callers handle this differently |

**Operational rule:** For ENCODE comparability, follow the Yeo lab eCLIP pipeline exactly: `umi_tools extract` (10 N R1) -> `cutadapt -q 6 -m 18` (two-pass) -> `bowtie2 pre-map rRNA` (unmapped to genome) -> `STAR --alignEndsType EndToEnd --outFilterMultimapNmax 1` -> `umi_tools dedup --method=unique`. For iCLIP/iCLIP2, follow the iCount preprocessing tutorial. Document any deviation in methods.

## References

- Van Nostrand EL et al 2016 Nat Methods 13:508 (eCLIP protocol, ENCODE standard)
- Konig J et al 2010 Nat Struct Mol Biol 17:909 (original iCLIP, truncation principle)
- Buchbender A et al 2020 Methods 178:33 (iCLIP2 protocol, library complexity gain)
- Lee FCY et al 2021 bioRxiv 2021.08.27.457890 (iiCLIP / improved iCLIP, motif specificity)
- Hafner M et al 2010 Cell 141:129 (PAR-CLIP, 4SU labeling, T->C signature)
- Zarnegar BJ et al 2016 Nat Methods 13:489 (irCLIP, non-radioactive)
- Aktas T et al 2020 Nucleic Acids Res 48:e15 (FLASH, fast protocol)
- Smith T et al 2017 Genome Res 27:491 (UMI-tools, network-based dedup)
- Daley T & Smith AD 2013 Nat Methods 10:325 (preseq library complexity)
- Chakrabarti AM et al 2023 Genome Biol 24:235 (nf-core/clipseq pipeline)

## Related Skills

- clip-seq/clip-alignment - Downstream STAR/bowtie2 alignment with ENCODE parameters
- clip-seq/clip-qc - Library complexity, FRiP, IDR, read-distribution QC after preprocessing
- clip-seq/crosslink-site-detection - Why preserving the 5' R2 base is critical
- clip-seq/clip-peak-calling - Downstream peak callers consume the dedup BAM
- clip-seq/stamp-antibody-free - STAMP/DART-seq use RNA-seq preprocessing instead of UMI-based CLIP preprocessing
- read-qc/umi-processing - General UMI handling concepts
- read-qc/adapter-trimming - General adapter trimming
- alignment-files/duplicate-handling - Picard MarkDuplicates fallback when UMIs unavailable
