---
name: bio-clip-seq-clip-peak-calling
description: Call protein-RNA binding sites from CLIP-seq BAM with CLIPper, PureCLIP, Skipper, Piranha, omniCLIP, CTK, CLAM, or Paraclu. Use when choosing between coverage-based, HMM-based, beta-binomial window-based, and crosslink-site-based peak callers; applying ENCODE eCLIP thresholds (log2 IP/SMInput >= 3, -log10 p >= 3); deciding when SMInput is mandatory; or reconciling peak-set discordance between callers for the same RBP.
tool_type: cli
primary_tool: CLIPper
---

## Version Compatibility

Reference examples tested with: CLIPper 2.0+, PureCLIP 1.3.1+, Piranha 1.2.1+, omniCLIP 0.2.0+, CTK 1.1.4+, CLAM 1.2+, Paraclu 9+, Skipper (commit 2023.05+), MACS3 3.0+, bedtools 2.31+, samtools 1.19+, idr 2.0.4+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed binary and adapt the example to match the actual CLI rather than retrying. PureCLIP 2.x changed several flag names; Skipper is distributed as a Snakemake workflow with frequently-evolving paths.

# CLIP-seq Peak Calling

**"Call protein-RNA binding sites from my deduplicated CLIP BAM"** -> Identify regions where read pile-up (and, for iCLIP/eCLIP, single-nucleotide truncations) exceed background from a size-matched input (SMInput) control. The choice of peak caller depends on (a) the CLIP variant (HITS-CLIP, iCLIP, eCLIP, PAR-CLIP), (b) whether SMInput is available, (c) the RBP binding mode (narrow motif vs broad zones vs repeat-binding), and (d) the goal (publication-comparable ENCODE peaks vs single-nucleotide crosslink sites vs high-recall site lists).

- CLI (ENCODE eCLIP canonical): `clipper -b dedup.bam -s hg38 -o peaks.bed --save-pickle --FDR-alpha 0.05 --superlocal` then `eclip-norm peaks.bed -i sminput.bam` for log2 fold-change against SMInput (`--FDR-alpha` is the CLIPper flag for the FDR cutoff; older docs sometimes shorten to `--FDR`)
- CLI (single-nucleotide crosslink sites, iCLIP/eCLIP): `pureclip -i dedup.bam -bai dedup.bam.bai -g genome.fa -ibam sminput.bam -ibai sminput.bam.bai -o sites.bed -or regions.bed -nt 8 -dm 8`
- CLI (high-recall, beta-binomial windowed; needs SMInput): `Skipper` Snakemake workflow with config matching cell type and SMInput BAM
- CLI (no SMInput, no truncation): `Piranha -b 50 -p 0.01 -d ZeroTruncatedNegativeBinomial -s -o peaks.bed dedup.bam` (Piranha takes the BAM as a positional argument, not via `-s`)
- CLI (CIMS/CITS single-nt from CTK): `tag2cluster.pl dedup.bed cluster.bed --multi-tag-method coverage`; then `bedExtractCIMS.pl cluster.bed cims.bed`
- CLI (multi-mapper rescue for repeat-binding RBPs): `CLAM peakcaller -i unique.bam multimap.bam -o clam_out_dir/ --gtf gencode.gtf` (CLAM peakcaller writes peaks into an OUTPUT DIRECTORY, not a single file)

The ENCODE eCLIP gold standard for "stringent" peaks is: log2(IP / SMInput) >= 3 AND -log10(p-value) >= 3. "Lenient" peaks use log2 >= 1 AND -log10 >= 2. Both filters operate on CLIPper peak output normalized against SMInput. Without SMInput, neither stringency level can be reproduced - the algorithm has no background term.

## Algorithmic Taxonomy

| Caller | Model | Resolution | SMInput required | Strength | Fails when |
|--------|-------|------------|------------------|----------|------------|
| CLIPper (Yeo) | Poisson with 500 bp local lambda, cubic-spline interpolation | Peak (10-500 nt) | Not for calling; required for ENCODE normalization | ENCODE eCLIP standard; reproducible against published ENCODE peak sets | Assumes most reads not from binding; fails when IP signal dominates a gene; sensitive to highly expressed transcripts |
| PureCLIP (Krakau 2017) | Non-homogeneous HMM jointly modeling enrichment + truncation + sequence biases | Single-nucleotide CL + broad region | Optional via -ibam; recommended | Single-nt resolution; only caller that explicitly models the iCLIP truncation pattern | Misses broad binding zones; very focal (Skipper 2023 benchmark: only ~4 sites per CLIP on test set; F1 ~0.2) |
| Skipper (Boyle 2023) | GC-stratified beta-binomial, 100 bp feature-respecting windows | Window (~100 nt) | Mandatory | 210-320% more sites than CLIPper; 8x faster; properly normalizes against input | Loses single-nt resolution; relatively new (2023); fewer published peak sets to compare |
| Piranha (Smith lab) | Zero-truncated negative binomial regression with optional covariates | Bin (50-200 nt) | Optional; pass as covariate | Mature, widely cited, handles count overdispersion | Biased toward high-expression transcripts; convergence fails with large covariate values (use `-l` log-space) |
| omniCLIP (Drewe-Boss 2018) | Non-homogeneous HMM with Dirichlet-multinomial of variants | Peak | Required | Models replicate variance; integrates background; can call peaks on any CLIP variant | Slow on deep libraries; F1 ~0.4 on standard benchmarks; fails for mitochondrial transcripts (Skipper benchmark: 0 chrM windows for FASTKD2) |
| CTK CIMS/CITS (Shah 2017) | Empirical FDR on crosslink-induced mutations or truncations | Single-nucleotide | No (uses background mutation rate from non-bound transcripts) | Single-nt resolution; works on HITS-CLIP deletions, PAR-CLIP T->C, iCLIP truncations | Empirical FDR less principled than HMM/beta-binomial; perl-based pipeline harder to integrate |
| CLAM peakcaller (Zhang & Xing 2017) | EM-assigned multi-mapper count + Piranha-like negative binomial | Peak | Optional | Only solution for repeat-binding RBPs; 10-30% additional sites in repeats | Inherits Piranha limitations on coverage-based stats; slower than Piranha |
| Paraclu (Kawaji) | Parametric clustering with min-density and max-density thresholds | Variable | No | Simple, parameter-tunable, works on bedgraph | Heuristic; no statistical significance; cluster boundaries sensitive to thresholds |
| PIPE-CLIP | Online CLIP pipeline with custom peak caller | Peak | Optional | Web interface; integrated end-to-end | Slow on cloud; less customizable; community has moved on |
| CLIPick (Park 2018) | Expression-deconvolved peak caller | Peak | No (RNA-seq used instead) | Models RNA-seq abundance as background | RNA-seq cannot capture nonspecific binding; less popular post-Skipper |
| Flipper (Tu 2024) | Negative-binomial differential test downstream of Skipper | Window | Yes | Companion differential tool to Skipper | Only meaningful in differential context; see clip-seq/differential-clip |
| MACS3 callpeak | Local Poisson | Peak | Optional (treats SMInput as ChIP-seq input) | Familiar, fast | Not designed for CLIP; misses truncation signal; produces wider-than-typical peaks |

Methodology evolves; verify the current ENCODE eCLIP standard operating procedure (encodeproject.org/eclip) and the nf-core/clipseq pipeline configuration before locking on a single caller. The 2023 Skipper benchmark (Boyle Cell Genomics) is the most recent comprehensive comparison and is the rationale for many recent eCLIP reanalyses.

## Critical Choice: Coverage-Based vs Crosslink-Site-Based vs Window-Based

Three fundamentally different statistical frameworks:

**Coverage-based (CLIPper, Piranha, MACS3, Paraclu):** Tally reads in a window; significance from local Poisson or negative binomial. Produces multi-nt peaks. Best when binding zones are broader than the read footprint (PUM2 3' UTR, HuR ARE elements).

**Crosslink-site-based (PureCLIP, CTK CIMS/CITS, PARalyzer):** Identify single nucleotides enriched in truncations (iCLIP/eCLIP) or specific mutations (PAR-CLIP T->C, HITS-CLIP deletions). Produces single-nt sites that can be aggregated into footprints. Best when single-nt resolution is essential (motif registration with mCross; allele-specific binding).

**Window-based (Skipper):** Tile transcriptome a priori into fixed-size feature-respecting windows; test each window's IP/IN ratio with beta-binomial. Produces window-level calls. Best when SMInput is available and the goal is high-recall, comparable-across-RBPs site sets.

| Goal | Recommended caller | Why |
|------|-------------------|-----|
| Reproduce ENCODE published peaks | CLIPper + IP/SMInput normalization | The Yeo lab pipeline is the canonical reference |
| Maximum sensitivity, modern statistics, SMInput available | Skipper | 210-320% more sites; explicit input normalization |
| Single-nucleotide crosslink-site map (motif registration, ASB) | PureCLIP (HMM) or CTK CITS (empirical) | True single-nt resolution |
| Repeat-binding RBP (Alu, LINE, LTR) | STAR multi-mapper + CLAM peakcaller | Only solution that uses multi-mappers |
| No SMInput control | Piranha or PureCLIP with no -ibam | Both work without input but lose specificity |
| PAR-CLIP T->C signature | PARalyzer or CTK CIMS (-substitution) | Designed for the T->C signal |
| HITS-CLIP deletion signature | CTK CIMS (-deletion) | Designed for the deletion signal |
| Custom statistical model needed | omniCLIP | Most flexible HMM among the options |
| Online interactive analysis | Galaxy CLIP-Explorer or PIPE-CLIP | Web interface |

## ENCODE eCLIP Stringency Thresholds

The Van Nostrand et al. ENCODE eCLIP standards define two stringency levels for CLIPper peaks normalized against SMInput:

| Stringency | log2(IP / SMInput) | -log10(p-value) | Use case |
|------------|--------------------|-----------------|----------|
| Stringent (publication) | >= 3 | >= 3 | High-confidence binding sites; differential analysis; motif discovery |
| Lenient (discovery) | >= 1 | >= 2 | Catalogue of all enriched windows; sensitive analyses |

Both thresholds are applied to the IDR-passing reproducible peak set. ENCODE requires:
- >= 2 biological replicates
- >= 1M unique fragments per replicate (or saturated peak detection)
- IDR rescue ratio and self-consistency ratio both < 2
- Narrow-binding RBPs: FRiP >= 0.005 (atypical-binding RBPs exempt)

```bash
# CLIPper -> SMInput normalization -> ENCODE thresholds
clipper -b dedup.bam -s hg38 -o peaks.bed --save-pickle --FDR-alpha 0.05 --superlocal

# Normalize against SMInput; eclip-norm is the Yeo lab tool
python overlap_peakfi_with_bam_PE.py peaks.bed dedup.bam sminput.bam dedup.bam.readnum.txt sminput.bam.readnum.txt peaks.normed.bed
python compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.py peaks.normed.bed peaks.compressed.bed

# Stringent filter: log2 FC >= 3, -log10 p >= 3
awk 'BEGIN{FS=OFS="\t"} $4 >= 3 && $5 >= 3' peaks.compressed.bed > peaks.stringent.bed

# Lenient
awk 'BEGIN{FS=OFS="\t"} $4 >= 1 && $5 >= 2' peaks.compressed.bed > peaks.lenient.bed
```

The Yeo lab scripts (`overlap_peakfi_with_bam_PE.py`, `compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.py`) live in the eclip-pipeline repository; they implement the exact log2 fold-change computation used for ENCODE data downloads.

## Per-Caller Failure Modes

### CLIPper -- High-expression transcript bias

**Trigger:** RBP that binds rare transcripts (snoRNAs, tRNAs, mitochondrial mRNAs) such as TROVE2 (Y RNA), NSUN2 (tRNA), or FASTKD2 (chrM).

**Mechanism:** CLIPper's Poisson local-lambda model assumes only a minority of reads come from binding. For rare-transcript binders, the IP reads ARE the majority signal on that transcript; the model treats them as background.

**Symptom:** Per-transcript peak count plummets compared to qualitatively-visible IP signal; chrM peaks called by Skipper but missed by CLIPper for FASTKD2.

**Fix:** Use Skipper for rare-transcript binders (beta-binomial against SMInput respects the actual IP/IN ratio independently of local lambda). Or pre-restrict CLIPper to the target transcript with `--gene custom.bed`.

### CLIPper -- Requires SMInput downstream

**Trigger:** Submitting CLIPper output as "peaks" without IP/SMInput log2 normalization.

**Mechanism:** CLIPper's own p-value is uncorrected for SMInput background. ENCODE peak BED files distributed on the portal are CLIPper + SMInput log2 normalization combined.

**Symptom:** Peak count from the analysis is 5-10x higher than ENCODE values for the same RBP; the peak set includes obvious housekeeping-gene noise.

**Fix:** Always normalize CLIPper output against SMInput. The Yeo lab scripts (`overlap_peakfi_with_bam_PE.py` + `compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.py`) produce the canonical ENCODE-style peak.

### PureCLIP -- Too focal; misses broad binding zones

**Trigger:** RBP with broad binding mode (PUM2 in 3' UTR clusters; SR proteins across exonic enhancer regions); using PureCLIP and disappointed.

**Mechanism:** PureCLIP's HMM emits single-nt crosslink-site state at high stringency. Skipper 2023 benchmark: PureCLIP F1 ~0.2 on RBFOX2 and similar bulk RBPs (only ~4 sites called per CLIP on test data).

**Symptom:** Site count 100x lower than expected; sites cluster within known binding regions but vast majority of region is "background" in PureCLIP output.

**Fix:** Use PureCLIP for single-nt CL maps (motif registration, ASB) but pair with CLIPper or Skipper for the broader binding-site list. PureCLIP's `-or` regions file gives broader output but is still focal compared to CLIPper.

### Piranha -- Top-expression-decile bias

**Trigger:** Highly expressed transcript (GAPDH, ACTB, rRNA-flanking) appears in most-significant Piranha peaks even when RBP biology is elsewhere.

**Mechanism:** ZTNB regression weights peak significance by absolute count, not relative enrichment. Genes with 10x baseline expression dominate the top of the peak list.

**Symptom:** Top-100 Piranha peaks include housekeeping genes; motif enrichment in top peaks is weak; functional GO term analysis of peak-bearing genes returns generic "cellular metabolism."

**Fix:** Use Piranha with SMInput as a covariate (`-c sminput_counts.bed -l`) so the regression accounts for input pile-up. Or switch to Skipper for principled input normalization.

### Piranha -- Convergence failure with covariates

**Trigger:** `-c rna_seq.bed` for expression normalization; large per-bin RNA-seq counts (> 10000).

**Mechanism:** ZTNB optimization fails on extreme covariate values; the iterative algorithm diverges.

**Symptom:** "Maximum iterations reached" or "Singularity in regression" error; or output BED file empty.

**Fix:** Pass `-l` to convert covariates to log-space; or pre-filter covariate BED to remove top-1% outlier bins; or switch to Skipper which handles input normalization internally.

### omniCLIP -- Mitochondrial transcript blind spot

**Trigger:** Mitochondrial RBP (FASTKD2, LRPPRC) submitted to omniCLIP.

**Mechanism:** omniCLIP's HMM is trained on autosomal/X transcripts; chrM has different background statistics that the model treats as outlier (or excludes if blacklisted).

**Symptom:** chrM peaks return 0; FASTKD2 peak file has no peaks at all.

**Fix:** Use Skipper for chrM-binding RBPs (GC-stratified beta-binomial is per-bin and indifferent to chromosome). Or call peaks on chrM separately with PureCLIP.

### CTK CIMS -- HITS-CLIP only; iCLIP/eCLIP misuse

**Trigger:** CTK CIMS (`-mutation` flag) applied to iCLIP/eCLIP data.

**Mechanism:** iCLIP/eCLIP use truncation (RT stops at adduct), not mutations. CIMS expects deletions/substitutions inside the read; CITS is the truncation-based companion.

**Symptom:** CIMS returns very few sites (mutation rate of CL is ~3-7% for HITS-CLIP, but only ~1-2% for iCLIP since most reads truncate before reaching the adduct).

**Fix:** For iCLIP/eCLIP, use CTK CITS (`tag2cluster.pl ... -cs5 5` truncation-mode). For HITS-CLIP and PAR-CLIP, CIMS is correct (HITS deletions, PAR T->C).

### CLAM peakcaller -- Multi-mapper BAM absent

**Trigger:** Ran STAR with `--outFilterMultimapNmax 1` then tried CLAM; CLAM fails or returns identical output to Piranha.

**Mechanism:** CLAM's EM rescue requires the multi-mapper BAM from STAR (`--outFilterMultimapNmax 100 --outSAMmultNmax -1`). Without it, CLAM falls back to unique-only Piranha-style peak calling.

**Symptom:** CLAM output looks identical to Piranha; repeat peaks not rescued.

**Fix:** Re-align with multi-mapper-aware STAR parameters (see clip-seq/clip-alignment). Verify the multi-mapper BAM has alignments at non-unique positions before running CLAM.

### MACS3 callpeak -- Wide peaks miss footprint

**Trigger:** MACS3 used because the user is familiar from ATAC/ChIP work.

**Mechanism:** MACS3 model is for DNA-protein binding at ~100-500 nt scale; CLIP-seq RBP footprints are 8-30 nt. MACS3 wide peaks (300+ nt) merge multiple distinct binding sites.

**Symptom:** MACS3 peak count is much lower than CLIPper/Skipper; mean peak width >> 100 nt; multiple known motif instances inside one MACS peak.

**Fix:** Use CLIP-specific callers. If MACS3 is forced (institutional pipeline), pass `--nomodel --shift -1 --extsize 50` for narrower peaks and accept the limitation.

## IDR for CLIP-seq

IDR (Li 2011) measures reproducibility of peak rankings across biological replicates. ENCODE eCLIP convention is per-rep + pooled + pseudoreplicates, identical to the ChIP-seq IDR pattern. The CLIP-specific consideration: rank by `signalValue` (column 7) NOT by `score` (column 5), because CLIPper's score is sparse and tied across many peaks.

```bash
# Per-replicate CLIPper + SMInput normalization yields .compressed.bed files
sort -k7,7nr rep1_peaks.compressed.bed > rep1.sorted
sort -k7,7nr rep2_peaks.compressed.bed > rep2.sorted

idr --samples rep1.sorted rep2.sorted \
    --input-file-type bed --rank 7 \
    --output-file idr.out --idr-threshold 0.05 \
    --plot --log-output-file idr.log
```

ENCODE IDR rules for eCLIP:
- Nt (true replicates) and Nself (pseudoreplicates): max(Nt, Nself) / min(Nt, Nself) <= 2 -> pass
- Both ratios > 2 -> library rejected
- IDR threshold 0.05 for true reps; 0.10 for pseudoreplicates

## Decision Tree by Scenario

| Scenario | Caller + parameters | Why |
|----------|---------------------|-----|
| Standard eCLIP with SMInput, ENCODE-comparable | CLIPper + SMInput log2 norm + IDR | Reproducible against ENCODE peak files |
| Standard eCLIP with SMInput, maximum sensitivity | Skipper Snakemake workflow | 210-320% more sites, beta-binomial input-norm |
| iCLIP/iCLIP2, want single-nt CL map | PureCLIP + SMInput (or empty bam if no SMI) | HMM jointly models enrichment + truncation |
| HITS-CLIP, want CIMS deletions | CTK `parseAlignment.pl` + `tag2cluster.pl` + `bedExtractCIMS.pl` | Empirical FDR on deletion-induced mutations |
| PAR-CLIP, T->C signal | PARalyzer or CTK CIMS `-substitution T C` | Designed for the T->C signature |
| Repeat-binding RBP (MATR3, LINE-1) | CLAM peakcaller on multi-mapper BAM | Only solution that recovers repeat peaks |
| No SMInput, mostly hopeless | Piranha `-d ZeroTruncatedNegativeBinomial` or PureCLIP w/o -ibam | Lose ENCODE comparability but get a peak set |
| Mitochondrial RBP (FASTKD2) | Skipper (window-based, chrM-friendly) | omniCLIP/CLIPper miss chrM |
| Highly multiplexed (SPIDR) | Custom Snakemake; no single-tool solution | Each barcode's reads are deep enough for CLIPper |
| Long-read CLIP (dirCLIP) | Custom analysis from BAM; long-read tools not yet mature | Isoform-resolved binding-site calls |
| Need maximum precision for motif analysis | PureCLIP (single-nt) + mCross downstream | mCross requires registered CL positions |
| Replicates with unequal library complexity (>2x diff in unique fragments) | Skipper (more robust to depth imbalance) OR `samtools view -s` down-sample to common depth then CLIPper + IDR | Skipper's GC-stratified beta-binomial handles per-sample depth; CLIPper without down-sampling lets the deeper replicate dominate joint analysis |

## Reconciliation: When Callers Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Skipper >> CLIPper site count | Skipper better recovers low-abundance and rare-transcript binders | Trust Skipper for rare transcripts; trust CLIPper for ENCODE-comparable counts |
| PureCLIP << CLIPper | PureCLIP is single-nt focal; CLIPper is broad peak | Both correct; report at the resolution the downstream analysis needs |
| Piranha top peaks dominated by housekeeping | No input normalization | Add SMInput as `-c` covariate or switch to Skipper |
| CLAM > Piranha by 10-30% | CLAM rescued multi-mappers in repeats | Trust CLAM if RBP biologically binds repeats |
| omniCLIP wider peaks than CLIPper | HMM smoothes binding state across nearby positions | Both correct; choose by downstream analysis needs |
| Per-rep CLIPper calls peak; pooled does not | Replicates inconsistent; one rep dominates | Inspect with IGV; trust IDR-passing peaks only |
| MACS3 < CLIPper peak count | MACS3 wide peaks merge multiple CLIP sites | Use CLIP-specific caller |
| FASTKD2 chrM peaks: Skipper finds; CLIPper / omniCLIP find none | Coverage-based callers blind to chrM lambda | Trust Skipper for organelle-binding RBPs |
| chimeric eCLIP peaks for snoRNA: standard tools find few | snoRNA-targeting eCLIP needs custom processing | Use VanNostrandLab/snoRNA-chimeric pipeline |

**Operational rule for high-confidence reporting:** Require a peak to pass (a) IDR <= 0.05 on true replicates, (b) log2(IP/SMInput) >= 3 AND -log10 p >= 3 (ENCODE stringent), (c) overlap with at least one independent caller within 100 nt. For motif analysis, additionally require PureCLIP single-nt sites within the peak. For ASB, require WASP-filtered alignment + PureCLIP CL sites overlapping heterozygous SNPs.

## Strand-Specific Peak Calling

CLIP libraries are strand-specific by design. Most callers handle strands correctly when the BAM is stranded (read flags 99/147/83/163 for paired). If a caller bug or unstranded library forces manual split:

```bash
# Plus strand reads (F flag NOT 0x10)
samtools view -h -b -F 16 dedup.bam > plus.bam
samtools index plus.bam
# Minus strand reads (F flag 0x10)
samtools view -h -b -f 16 dedup.bam > minus.bam
samtools index minus.bam

# Call peaks on each strand independently, then merge
clipper -b plus.bam -s hg38 -o peaks_plus.bed
clipper -b minus.bam -s hg38 -o peaks_minus.bed
cat peaks_plus.bed peaks_minus.bed | sort -k1,1 -k2,2n > peaks_stranded.bed
```

For paired-end CLIP, mate orientation determines strand: RF (R1 reverse) means R2 5' is the truncation = CL site -1 on plus-strand transcripts. ENCODE eCLIP convention: read 2 of pair is the truncation-bearing read.

## Multi-Mapper Rescue with CLAM (Repeat Binders)

```bash
# Assume STAR was run with --outFilterMultimapNmax 100 --outSAMmultNmax -1
samtools sort -o multimap.bam aligned_multimap.bam
samtools index multimap.bam

# Split unique vs multi
CLAM preprocessor -i multimap.bam -o clam_out/ --read-tagger-method median

# EM-based reassignment
CLAM realigner -i clam_out/unique.sorted.bam -o clam_out/ --winsize 50 --max-tags 0

# Peak calling with reassigned multi-mappers
CLAM peakcaller -i clam_out/unique.sorted.bam clam_out/realigned.sorted.bam \
    -o clam_peaks.bed -p 8 --gtf gencode.v38.annotation.gtf
```

CLAM rescues 10-30% additional peaks in repeat-rich regions. The trade-off: peak coordinates in repeat instances are probabilistic (single best from EM), so motif analysis on CLAM repeat peaks requires careful interpretation.

## Quality Metrics for Peak Sets

```python
import pandas as pd

def peak_qc(peaks_bed):
    peaks = pd.read_csv(peaks_bed, sep='\t', header=None,
                        names=['chrom','start','end','name','score','strand'])
    peaks['width'] = peaks['end'] - peaks['start']
    return {
        'n_peaks': len(peaks),
        'mean_width': peaks['width'].mean(),
        'median_width': peaks['width'].median(),
        'width_5p_95p': (peaks['width'].quantile(0.05), peaks['width'].quantile(0.95)),
        'chroms_with_peaks': peaks['chrom'].nunique(),
        'peaks_chrM': (peaks['chrom'] == 'chrM').sum()
    }
```

| Metric | Expected | Interpretation |
|--------|----------|----------------|
| n_peaks (CLIPper stringent) | 5k - 100k | <5k = under-IP or stringency too high; >>100k = noisy library |
| Mean peak width (CLIPper) | 40-200 nt | Wider = caller over-merging; narrower = under-coverage |
| Peak width 5-95% range | 20-500 nt | Tight range = consistent peaks; wide = mixed sharp + broad |
| FRiP (fraction reads in peaks) | >= 0.005 | ENCODE narrow-binding RBP minimum; many RBPs higher |
| chrM peaks | 0 unless mt-RBP | FASTKD2, LRPPRC expected to have chrM enrichment |
| TPM-rank distribution | Spread across deciles | Concentrated in top decile = caller has expression bias (Piranha symptom) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| CLIPper "no peaks found" | BAM not sorted/indexed; gene annotation missing | `samtools sort && index`; verify `-s hg38` (or `--gene custom.bed`) |
| CLIPper out-of-memory on 50M BAM | `--save-pickle` keeps intermediates in RAM | Remove `--save-pickle`; or use machine with > 64 GB |
| PureCLIP slow (>24h) | `-iv` interval file missing | Restrict to one chromosome with `-iv chr1` for testing; provide GTF |
| Piranha returns no peaks | ZTNB convergence failed silently | Add `-l`; or reduce `-b` bin size; or switch to Skipper |
| omniCLIP segmentation fault | Replicate BAMs differ in chromosome order | Re-sort all BAMs with identical reference dictionary |
| CTK CITS finds nothing on PAR-CLIP | Used CITS (truncation) for substitution data | Use CIMS `-substitution T C` for PAR-CLIP |
| Multi-mapper rescue returns 0 extra peaks | STAR was run unique-only; no multi-mapper BAM | Re-align with `--outFilterMultimapNmax 100 --outSAMmultNmax -1` |
| ENCODE stringent peak count = 0 | Forgot SMInput normalization step | Run Yeo lab `overlap_peakfi_with_bam_PE.py` first; thresholds apply to normalized peaks |
| Strand information lost | `samtools view -F 16` filter applied to merged BAM | Strand is encoded in BAM flags; do not pre-split |
| IDR returns no reproducible peaks | Wrong ranking column | Sort by signalValue / log2FC, not score; for narrowPeak format use `--rank p.value` |

## References

- Van Nostrand EL et al 2016 Nat Methods 13:508 (eCLIP, SMInput, ENCODE pipeline)
- Lovci MT et al 2013 Nat Struct Mol Biol 20:1434 (CLIPper algorithm)
- Krakau S et al 2017 Genome Biol 18:240 (PureCLIP HMM)
- Boyle EA et al 2023 Cell Genomics 3:100317 (Skipper, GC-stratified beta-binomial)
- Uren PJ et al 2012 Bioinformatics 28:3013 (Piranha, ZTNB)
- Drewe-Boss P et al 2018 Genome Biol 19:183 (omniCLIP)
- Shah A et al 2017 Bioinformatics 33:566 (CTK / CIMS / CITS)
- Zhang Z & Xing Y 2017 Nucleic Acids Res 45:9260 (CLAM multi-mapper)
- Frith MC et al 2008 Genome Res 18:1-12 (paraclu)
- Li Q et al 2011 Ann Appl Stat 5:1752 (IDR framework)
- Van Nostrand EL et al 2020 Nature 583:711 (Principles of RNA processing from 150 eCLIP maps)
- ENCODE eCLIP Standards (encodeproject.org/eclip) - canonical thresholds

## Related Skills

- clip-seq/clip-preprocessing - Upstream preprocessing producing dedup BAM
- clip-seq/clip-alignment - Upstream alignment with crosslink-preserving parameters
- clip-seq/clip-qc - FRiP, library complexity, IDR before/after peak calling
- clip-seq/crosslink-site-detection - PureCLIP / CTK CITS / PARalyzer for single-nt CL sites
- clip-seq/differential-clip - DEWSeq / Flipper for cross-condition differential binding
- clip-seq/binding-site-annotation - Annotate called peaks to 5' UTR / CDS / 3' UTR / intron
- clip-seq/clip-motif-analysis - Motif discovery on called peaks
- clip-seq/ago-clip-mirna-targets - Chimeric eCLIP peak calling
- clip-seq/m6a-clip - miCLIP2-specific peak calling
- chip-seq/peak-calling - DNA-protein peak calling for comparison
