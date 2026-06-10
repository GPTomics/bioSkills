---
name: bio-methylation-calling
description: Extracts per-cytosine methylation calls from aligned bisulfite/EM-seq reads with bismark_methylation_extractor (Bismark BAM) or the aligner-agnostic MethylDackel/BISCUIT (bwa-meth BAM), producing the beta value M/(M+U) as a coverage file, bedGraph, or genome-wide cytosine report across CpG/CHG/CHH context. Covers conversion-rate QC as the first gate, the 5mC vs 5hmC summed caveat, variant-aware calling so a C/T SNP does not masquerade as unmethylation, paired-end --no_overlap double-counting, symmetric CpG dyad collapse, and the 0-based vs 1-based coordinate trap. Use when extracting methylation levels from a bisulfite/EM-seq alignment, choosing an extractor for a non-Bismark BAM, QC-ing conversion efficiency, or producing coverage/cytosine-report input for testing. For long-read MM/ML modification calling see long-read-sequencing/nanopore-methylation; for the upstream BAM see bismark-alignment; for per-CpG statistics see differential-cpg-testing.
tool_type: cli
primary_tool: Bismark
---

## Version Compatibility

Reference examples tested with: Bismark 0.24+, MethylDackel 0.6+, samtools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The extractor must match the aligner: bismark_methylation_extractor reads Bismark's own XM call string and CANNOT process a bwa-meth BAM (no XM tag); MethylDackel and BISCUIT recompute the call from BAM + reference and are aligner-agnostic. The genome FASTA build (hg38 vs T2T-CHM13) defines every coordinate downstream, and BISCUIT/Bis-SNP need a SNP-aware reference workflow. Always run `<tool> --help` on the installed build before quoting a default; MethylDackel `--maxVariantFrac` and EM-seq end-trim recommendations have no universal value.

# Methylation Calling

**"Get methylation levels from my bisulfite BAM"** -> Count converted-vs-unconverted bases at each reference cytosine and divide - because a methylation level is a per-cytosine COUNT RATIO, not a measurement, and every number is hostage to conversion completeness, SNPs, overlap, and the 5mC/5hmC conflation.
- CLI (Bismark BAM): `bismark_methylation_extractor -p --comprehensive --bedGraph --cytosine_report --genome_folder genome/ sample_pe.bam`
- CLI (bwa-meth BAM): `MethylDackel extract --mergeContext ref.fa sample.bam`

Scope: per-cytosine M/U counting from a short-read bisulfite/EM-seq/TAPS alignment, the conversion/M-bias/variant QC gates, and the coverage/cytosine-report handoff. Native long-read MM/ML modification calling -> long-read-sequencing/nanopore-methylation. The aligned BAM -> bismark-alignment. Per-CpG statistics on the counts -> differential-cpg-testing. Region calling -> dmr-detection. Array (EPIC/450K) beta-from-intensity is a different readout entirely and is not this skill.

## The Single Most Important Modern Insight -- A Methylation Level Is a Count Ratio, Not a Measurement

The extractor does not measure methylation. Bisulfite/EM-seq turns an epigenetic state into a sequence state (an unmethylated C converts to T; a methylated C stays C), and the extractor tallies, at each reference cytosine, M reads that show C and U reads that show T, then reports `beta = M / (M + U)`. Every beta is therefore hostage to four upstream choices the extractor cannot fix:

1. **Whether conversion was complete.** An unconverted unmethylated C is byte-for-byte identical to a methylated C, so incomplete conversion inflates EVERY beta and is invisible per-site. It must be measured globally before any beta is trusted (CHH rate or a lambda spike-in).
2. **Whether a C/T SNP is hiding as unmethylation.** The reference says C, the sample carries a T allele, the read shows T, and the extractor scores "converted" = unmethylated. Polymorphic CpGs are silently miscalled and can fabricate SNP-driven DMRs.
3. **Whether the paired-end overlap was double-counted.** A cytosine in the R1/R2 overlap is one molecule observed twice; counting both votes inflates coverage and biases beta when the mates disagree.
4. **Whether "5mC" is actually 5mC+5hmC.** 5hmC is also protected from conversion, so a standard WGBS/EM-seq beta is the SUM 5mC+5hmC. In brain/ESC/liver this is a real misattribution that no extractor flag can undo.

Organize the analysis around defending these four (conversion QC -> overlap -> M-bias -> variant-awareness -> what context/chemistry am I even calling), not around listing flags. And one deeper caveat: beta is itself a lossy average - collapsing reads to a `.cov` matrix discards the read-level epiallele, the strand (hemimethylation), and the allele (ASM); see the last sections.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| bismark_methylation_extractor | Krueger & Andrews 2011 *Bioinformatics* 27:1571 | reads Bismark's own XM call string; ecosystem standard | input is a Bismark BAM; want the genome-wide cytosine report / methylKit `.cov`; plant `--CX` |
| MethylDackel | github.com/dpryan79/MethylDackel (no journal) | recomputes calls from BAM + reference; fast; aligner-agnostic | bwa-meth BAM; want built-in M-bias `--OT/--OB` bounds and `--maxVariantFrac` SNP guard |
| BISCUIT | github.com/huishenlab/biscuit (no journal) | aligner + JOINT methylation/SNP/ASM caller; VCF + epiBED | human population / allele-specific work needing SNP-aware betas |
| Bis-SNP | Liu 2012 *Genome Biol* 13:R61 | Bayesian joint genotype + methylation | the original SNP-aware caller; mask C/T SNPs from betas |
| methylpy | github.com/yupenghe/methylpy (no journal) | allc format + binomial methylated-flag | the allc / ALLCools ecosystem; built-in per-site significance |
| asTair | bitbucket.org/bsblabludwig/astair (no journal) | polarity-aware caller (mCtoT / CtoT) | TAPS data (inverted polarity); explicit per-context stats |

EM-seq extraction arithmetic is identical to bisulfite (unmodified C reads as T either way) - the same Bismark/MethylDackel command works. TAPS INVERTS the polarity (modified C->T), so a bisulfite extractor reports 1-beta; use ASTAIR `--mod_mapping mCtoT`.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bismark-aligned BAM, mammalian CpG | `bismark_methylation_extractor -p --comprehensive --cytosine_report` | reads the XM tag; cytosine report is the bsseq/methylKit input |
| bwa-meth-aligned WGBS/EM-seq | `MethylDackel extract --mergeContext` | bismark extractor cannot read a bwa-meth BAM (no XM tag) |
| Human population / polymorphic CpGs / ASM | BISCUIT (or Bis-SNP) | joint SNP+methylation; C/T SNP recognized as a variant, not unmethylation |
| Plant methylome (CHG/CHH real) | `--CX` (Bismark) / `--CHG --CHH` (MethylDackel) + lambda spike-in | non-CpG is biology; CHH-as-conversion-proxy fails |
| TAPS data | ASTAIR `--mod_mapping mCtoT` | polarity inverted; a bisulfite extractor reports 1-beta |
| Need 5hmC resolved from 5mC | second chemistry (oxBS/TAB/ACE) -> route, not a flag | WGBS/EM-seq sums them; the extractor cannot split |
| Read-level heterogeneity / clonality / cfDNA | keep reads (epiread/epiBED), do NOT collapse to `.cov` | the epiallele dies in the per-CpG average |
| Long-read modBAM with MM/ML tags | -> long-read-sequencing/nanopore-methylation | native modification tags, not conversion counting |

## Conversion-Rate QC -- the First Gate

**Goal:** Reject a methylome whose betas are inflated by incomplete conversion before computing anything.

**Approach:** Measure non-conversion globally, two ways. In mammals, genome-wide CHH methylation is near-zero biology, so the observed CHH rate IS the non-conversion rate (~2% CHH implies ~98% conversion). Where non-CpG methylation is real (plants, ESCs, neurons - Schultz 2015), spike in unmethylated lambda phage DNA and align to it: all its cytosines are unmethylated, so its observed "methylation" is the non-conversion rate. Require conversion >=99% (community convention; Shirane 2013 is a representative source). Extract CHH at least once even for a CpG-only mammalian study, purely for this number.

```bash
# Bismark: extract all contexts, then read the CHH methylation rate from the splitting report
bismark_methylation_extractor -p --comprehensive --CX --genome_folder genome/ sample_pe.bam
# The *_splitting_report.txt prints C methylated in CHH context %; (100 - that) ~ conversion efficiency in mammals.

# Lambda spike-in (the gold standard, mandatory in plants/ESC/neurons):
# align to the lambda genome, extract, and treat its global methylation as the non-conversion rate.
bismark --genome lambda_genome/ -1 R1.fq.gz -2 R2.fq.gz -o lambda_qc/
```

MethylDackel `--minConversionEfficiency` can additionally drop individual reads whose own conversion (from non-CpG Cs) is too low - a per-read complement to the global gate, not a replacement for it.

## Extract from a Bismark BAM

```bash
bismark_methylation_extractor \
    -p \                                # paired-end (--no_overlap is ON BY DEFAULT for -p)
    --comprehensive \                   # merge OT/OB/CTOT/CTOB into one file per context
    --bedGraph --cytosine_report \      # coverage/bedGraph + genome-wide every-C report
    --genome_folder genome/ \           # MANDATORY for --cytosine_report (scans the FASTA for all Cs)
    --ignore_r2 2 \                     # trim R2 5' end-repair artifact (set from M-bias; near-universal for EM-seq/PBAT)
    --parallel 4 --gzip \
    -o methylation/ sample_pe.bam
# Outputs: sample_pe.bismark.cov.gz (1-BASED), sample_pe.bedGraph.gz (0-BASED), sample_pe.CpG_report.txt.gz (1-BASED).
```

Collapse the symmetric CpG dyad with `coverage2cytosine --merge_CpG` (NOT `--merge_non_CpG`, which merges the CHG+CHH files). `--merge_CpG` adds the + strand C at position p and the - strand C at p+1 into one dyad entry, doubling effective coverage; it is CpG-only and incompatible with `--CX`.

## Extract from a bwa-meth BAM (MethylDackel)

```bash
MethylDackel mbias ref.fa sample.bam mbias_prefix   # inspect the SVGs; it SUGGESTS --OT/--OB bounds (do NOT accept blindly)
MethylDackel extract \
    --mergeContext \                   # collapse the symmetric CpG dyad (the MethylDackel equivalent of --merge_CpG)
    --maxVariantFrac 0.25 \            # exclude a C if the opposite-strand non-G fraction exceeds this (cheap SNP guard)
    --OT 3,0,0,98 --OB 3,0,0,98 \     # inclusion bounds from mbias: first/last bp to keep on read1,read2
    ref.fa sample.bam
# Output sample_CpG.bedGraph is 0-BASED, half-open: chr start end round(%meth) count_M count_U. CpG ONLY by default;
# add --CHG --CHH for plants. Default --minDepth 1, -q (MAPQ) 10, -p (Phred) 5.
```

## Coverage, Precision, and the Handoff -- Pass Counts, Not Betas

A beta is a binomial proportion: at coverage n its granularity is 1/n and its SE is ~sqrt(beta(1-beta)/n). At n=1 beta is only 0 or 1; at n=4 it lands on {0,.25,.5,.75,1}. A 1/2 site and a 50/100 site both read beta=0.5 but carry wildly different evidence. The extractors impose no biological floor (MethylDackel `--minDepth` and Bismark `.cov` both default to >=1). Do NOT pre-threshold to a single beta and t-test it - that discards the coverage information. Hand the COUNT data (M and total) forward; the downstream DMR callers (DSS/methylKit/bsseq) model the beta-binomial and USE n. A common per-CpG floor is >=10x AFTER symmetric merge, but it is a tradeoff (site count vs precision), not a magic constant.

## 5mC vs 5hmC and the Oxidation Cascade -- a Chemistry Decision, Not a Flag

Standard bisulfite AND standard EM-seq protect BOTH 5mC and 5hmC from conversion, so the "methylated" bin is 5mC+5hmC, never 5mC alone. Worse, bisulfite DEAMINATES the downstream TET-oxidation products 5fC and 5caC, so they read as T and land in the "unmethylated" bin - both bins are biochemically impure at TET-active loci (ESC, early embryo, neurons, some tumors). The full cascade is 5mC -> 5hmC -> 5fC -> 5caC. Resolving any derivative requires a SECOND wet-lab chemistry, not an extractor flag: 5hmC via oxBS-seq (Booth 2012; 5hmC = BS - oxBS by subtraction), TAB-seq (Yu 2012; direct 5hmC), or ACE-seq (Schutsky 2018; enzymatic, low-input); the bisulfite-free TAPS (Liu 2019) sidesteps the harsh chemistry but inverts polarity. Never let a plain WGBS/EM-seq beta be labeled "5mC" or its complement "unmodified C."

## Read-Level Heterogeneity -- the Epiallele Dies in the Average

beta=0.5 is consistent with three opposite biologies that beta cannot distinguish: a 50/50 mixture of fully-methylated and fully-unmethylated cells, every cell ~50% methylated with CpGs scattered differently per molecule (stochastic disorder), or a true uniform intermediate. The discriminator lives in the JOINT CpG pattern on a single read - intramolecular co-methylation - which the per-CpG average integrates out. Read-level metrics (PDR, Landau 2014; epipolymorphism, Landan 2012; methylation entropy; MHL for cfDNA tissue-of-origin, Guo 2017; FDRP/qFDRP) measure clonal/epigenetic instability and power liquid-biopsy deconvolution, and they require a READ-PRESERVING format (BISCUIT epiread/epiBED or the raw BAM) plus tools like Metheor or methclone - none of which are extractors. The moment the pipeline collapses to a `.cov` beta matrix, the epiallele is gone. If the question is heterogeneity, clonality, or cfDNA deconvolution, do NOT collapse to beta; this deeper analysis may warrant its own workflow, but at minimum the calling stage must flag that beta discards it.

## Hemimethylation and Allele-Specific Methylation -- two more things the average hides

**Hemimethylation (the strand axis).** `--merge_CpG` / `--mergeContext` is not a neutral coverage optimization: it bakes in the assumption that the two strands of a dyad agree and silently zeroes the hemimethylation channel (one strand methylated, the other not). Hemimethylation is real biology - the obligate post-replication maintenance intermediate, and a stable heritable mark at CTCF/cohesin sites required for chromatin looping (Xu & Corces 2018). Default destranding ON for bulk symmetric-CpG DMR work; turn it OFF (strand-specific extraction) and budget high per-strand coverage the moment strand asymmetry is the question. Per-molecule dyad state needs hairpin-bisulfite (Laird 2004).

**Allele-specific methylation (the allele axis).** The same C/T SNP that corrupts a beta (insight #2) becomes the measurement axis once reads are phased: assign each read's methylation to the SNP allele it carries and compare beta per allele (BISCUIT `epiread -B snps.bed` then `biscuit asm`; or Bis-SNP). The headline trap: a stable ~50% beta at a KNOWN imprinted control region (H19/IGF2, KCNQ1OT1, SNRPN, GNAS, MEG3) is NOT intermediate methylation and NOT a conversion artifact - it is two superimposed monoallelic states (one allele ~100%, one ~0%) averaged into a deceptive midpoint, a signature to PHASE, not a value to model. Sequence-dependent ASM is the mQTL bridge to causal-genomics/mendelian-randomization; produce the allele-phased betas here, do the colocalization there.

## Per-Method Failure Modes

### No conversion-rate gate
**Trigger:** reporting betas without checking CHH rate or a lambda spike-in. **Mechanism:** an unconverted unmethylated C is identical to a methylated C; non-conversion inflates every beta. **Symptom:** uniformly elevated methylation, fabricated low-methylation regions. **Fix:** require >=99% conversion (CHH-proxy in mammals; lambda spike-in in plants/ESC/neurons).

### bismark extractor on a bwa-meth BAM
**Trigger:** `bismark_methylation_extractor` on a non-Bismark BAM. **Mechanism:** it reads Bismark's XM tag, absent from bwa-meth output. **Symptom:** error or empty/garbage calls. **Fix:** use MethylDackel or BISCUIT, which recompute from BAM + reference.

### C/T SNP read as unmethylation
**Trigger:** extracting human population data with no variant-awareness. **Mechanism:** a T allele at a reference C is scored as converted. **Symptom:** spurious hypomethylation at polymorphic CpGs; SNP-driven false DMRs. **Fix:** BISCUIT/Bis-SNP joint calling, or MethylDackel `--maxVariantFrac`, or mask dbSNP/sample-VCF CpGs.

### Paired-end overlap double-counted
**Trigger:** single-end mode on paired data, or a non-default tool. **Mechanism:** the R1/R2 overlap counts one molecule twice. **Symptom:** inflated coverage, biased beta on mate disagreement. **Fix:** `-p` (Bismark `--no_overlap` is on by default for paired-end); MethylDackel handles overlap automatically.

### M-bias not trimmed
**Trigger:** extracting without inspecting the M-bias plot. **Mechanism:** end-repair fill-in introduces unmethylated Cs at fragment ends; the per-position methylation deviates near read ends. **Symptom:** a non-flat M-bias curve; biased calls in the trimmable region. **Fix:** trim with `--ignore`/`--ignore_r2` (Bismark) or `--OT/--OB` bounds (MethylDackel). The diagnostic is FLATNESS/positional stability, not any particular global level.

### Coordinate base mismatch on a join
**Trigger:** joining a MethylDackel bedGraph (0-based) to a Bismark `.cov` (1-based). **Mechanism:** the two formats index differently. **Symptom:** every site shifts by one; strands silently mismatch. **Fix:** pick one format end-to-end or convert explicitly; the cytosine report and allc are 1-based, both bedGraphs are 0-based.

### --merge_non_CpG mistaken for dyad collapse
**Trigger:** using `--merge_non_CpG` to merge the symmetric CpG strands. **Mechanism:** `--merge_non_CpG` merges the CHG+CHH output files, NOT the CpG dyad. **Symptom:** strand-specific CpG report unchanged; CpG coverage not doubled. **Fix:** symmetric CpG collapse is `coverage2cytosine --merge_CpG` / MethylDackel `--mergeContext`.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Conversion efficiency >=99% | Shirane 2013 *PLoS Genet* 9:e1003439; community | below it, residual unconverted Cs inflate every beta |
| CHH rate ~= (1 - conversion) in mammals | Schultz 2015 *Nature* 523:212 | true mammalian CHH is near-zero, so it reads out non-conversion |
| Per-CpG coverage >=10x (after merge) | community | granularity 1/n; SE of an intermediate beta ~0.16 at 10x |
| `--ignore_r2` ~2 for EM-seq/PBAT | M-bias plot | R2 5' end-repair fill-in adds unmethylated Cs; set from the plot |
| MethylDackel `-q` 10 / `-p` 5 | MethylDackel defaults | MAPQ/base-quality minima; defaults, confirm with `--help` |
| `--maxVariantFrac` no universal value | MethylDackel docs | set per-experiment against known-SNP density |
| `--merge_CpG` doubles dyad coverage | Bismark docs | + and - strand of a symmetric CpG are co-methylated; merge then threshold |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Empty/garbage calls from bismark extractor | bwa-meth BAM (no XM tag) | use MethylDackel/BISCUIT |
| Uniformly high methylation | incomplete conversion | check CHH rate / lambda; require >=99% |
| Spurious hypomethylation at SNPs | C/T SNP read as conversion | `--maxVariantFrac`; BISCUIT/Bis-SNP; mask dbSNP |
| Inflated coverage on overlaps | single-end mode on paired data | `-p`; MethylDackel handles it automatically |
| Sites shifted by one after a join | 0-based vs 1-based format mix | reconcile coordinate base; one format end-to-end |
| Plant CHH/CHG methylome missing | CpG-only default | `--CX` (Bismark) / `--CHG --CHH` (MethylDackel) |
| `coverage2cytosine: option --merge_CpG` ignored | used `--merge_non_CpG` instead, or paired with `--CX` | `--merge_CpG` is CpG-only, incompatible with `--CX` |
| Inverted (1-beta) methylome | bisulfite extractor on TAPS data | ASTAIR `--mod_mapping mCtoT` |

## References

- Krueger F, Andrews SR. 2011. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. *Bioinformatics* 27:1571-1572.
- Liu Y, Siegmund KD, Laird PW, Berman BP. 2012. Bis-SNP: combined DNA methylation and SNP calling for Bisulfite-seq data. *Genome Biol* 13:R61.
- Vaisvila R, Ponnaluri VKC, Sun Z, et al. 2021. Enzymatic methyl sequencing detects DNA methylation at single-base resolution from picograms of DNA. *Genome Res* 31:1280-1289.
- Booth MJ, Branco MR, Ficz G, et al. 2012. Quantitative sequencing of 5-methylcytosine and 5-hydroxymethylcytosine at single-base resolution. *Science* 336:934-937.
- Yu M, Hon GC, Szulwach KE, et al. 2012. Base-resolution analysis of 5-hydroxymethylcytosine in the mammalian genome. *Cell* 149:1368-1380.
- Schutsky EK, DeNizio JE, Hu P, et al. 2018. Nondestructive, base-resolution sequencing of 5-hydroxymethylcytosine using a DNA deaminase. *Nat Biotechnol* 36:1083-1090.
- Liu Y, Siejka-Zielinska P, Velikova G, et al. 2019. Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution. *Nat Biotechnol* 37:424-429.
- Schultz MD, He Y, Whitaker JW, et al. 2015. Human body epigenome maps reveal noncanonical DNA methylation variation. *Nature* 523:212-216.
- Shirane K, Toh H, Kobayashi H, et al. 2013. Mouse oocyte methylomes at base resolution reveal genome-wide accumulation of non-CpG methylation and role of DNA methyltransferases. *PLoS Genet* 9:e1003439.
- Landau DA, Clement K, Ziller MJ, et al. 2014. Locally disordered methylation forms the basis of intratumor methylome variation in chronic lymphocytic leukemia. *Cancer Cell* 26:813-825.
- Landan G, Cohen NM, Mukamel Z, et al. 2012. Epigenetic polymorphism and the stochastic formation of differentially methylated regions in normal and cancerous tissues. *Nat Genet* 44:1207-1214.
- Guo S, Diep D, Plongthongkum N, et al. 2017. Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples and tumor tissue-of-origin mapping from plasma DNA. *Nat Genet* 49:635-642.
- Xu C, Corces VG. 2018. Nascent DNA methylome mapping reveals inheritance of hemimethylation at CTCF/cohesin sites. *Science* 359:1166-1170.
- Laird CD, Pleasant ND, Clark AD, et al. 2004. Hairpin-bisulfite PCR: assessing epigenetic methylation patterns on complementary strands of individual DNA molecules. *PNAS* 101:204-209.
- MethylDackel. github.com/dpryan79/MethylDackel (no associated journal publication).
- BISCUIT. github.com/huishenlab/biscuit (no associated journal publication).

## Related Skills

- bismark-alignment - Produces the aligned BAM consumed here
- differential-cpg-testing - Per-CpG statistical testing on the counts
- dmr-detection - Region-level methods downstream
- methylkit-analysis - methylKit import of the coverage/cytosine report
- long-read-sequencing/nanopore-methylation - Native long-read MM/ML modification calling (the wall)
- causal-genomics/mendelian-randomization - mQTL / causal follow-up of allele-specific methylation
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
