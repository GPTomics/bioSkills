---
name: bio-clip-seq-binding-site-annotation
description: Annotate CLIP-seq peaks or crosslink sites to RNA features (5'UTR, CDS, 3'UTR, intron, splice junction, snoRNA, tRNA, ncRNA, repeat elements) with ChIPseeker, RCAS, RBP-Maps (Yeo splicing regulatory maps), and bedtools, applying feature-priority hierarchies, transcript-context resolution, and metagene aggregation. Use when characterizing where in transcripts an RBP binds, comparing peak distribution across regions, generating splicing-regulatory maps relative to alternative-splicing events, or distinguishing exonic vs intronic vs UTR binding.
tool_type: mixed
primary_tool: ChIPseeker
---

## Version Compatibility

Reference examples tested with: ChIPseeker 1.40+, RCAS 1.30+, GenomicFeatures 1.56+, GenomicRanges 1.56+, rbp-maps (Yeo github), bedtools 2.31+, pybedtools 0.10+, pyranges 0.0.129+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt the example to match the actual API rather than retrying. ChIPseeker 1.40+ changed the priority defaults; verify the priority vector before pipelining.

# Binding Site Annotation

**"Annotate where in transcripts my RBP binds"** -> Map CLIP peaks or single-nucleotide crosslink sites to RNA features and report the per-feature distribution. The interpretation is RBP-class-specific: splicing factors (PTBP1, U2AF2, RBFOX) bind intron-exon junctions; mRNA-stability regulators (HuR, PUM2) bind 3' UTRs; translation factors (EIF3J, RPS19) bind 5' UTRs and CDS; and small-ncRNA-binding RBPs (NSUN2 tRNAs, LARP7 7SK, TROVE2 Y-RNAs) bind specific non-coding transcripts. A correct annotation pipeline (a) resolves overlapping features by priority, (b) preserves transcript-isoform context, (c) generates metagene distributions, and (d) flags repeat-element overlap separately.

- R (peak-level, fast): `ChIPseeker::annotatePeak(peaks, TxDb=txdb, level='gene', tssRegion=c(-100,100))` then `plotAnnoPie(anno)`
- R (transcript-level with RNA-specific regions): RCAS `runReport(queryRegions=peaks, gffData=gencode_gtf, genomeVersion='hg38')`
- CLI (regions only): `bedtools intersect -s -wa -wb -a peaks.bed -b features.bed` (manual hierarchy)
- R (splicing-regulatory map for splice factors): RBP-Maps (`yeolab/rbp-maps`) generates the 1400 nt vectorized cassette-exon map used in Yeo lab ENCODE papers
- CLI (metagene aggregation): `deepTools computeMatrix` or RSeQC `geneBody_coverage.py` for read profiles

The Yeo lab convention for ENCODE eCLIP: ChIPseeker for global region distribution + RBP-Maps for position-resolved splicing context + custom analysis for repeat-element overlap. ChIPseeker alone over-counts intronic binding because TSS-region default extends too far for RNA features.

## Algorithmic Taxonomy

| Tool | Input | Annotation level | Feature priority | Output | Strength | Fails when |
|------|-------|------------------|------------------|--------|----------|------------|
| ChIPseeker (R) | Peak BED + TxDb | Gene-level or transcript-level | Promoter > 5' UTR > 3' UTR > Exon > Intron > Intergenic | Per-peak annotation + pie chart + distance-to-TSS | Mature, widely used, fast | TSS-region default `c(-3000, 3000)` over-extends; promoter category meaningless for RNA |
| RCAS (R) | Peak BED + GFF | Transcript-level | Customizable | HTML report + per-region table | RNA-specific design; ncRNA-aware | Slow; HTML-heavy; less actively maintained |
| RBP-Maps (Yeo) | eCLIP BAM + alternative splicing event tables | Splice junction-level | NA (positional metagene) | Splicing regulatory map (1400nt vector around cassette exon) | The standard for splicing-factor CLIP analysis | Splicing-only; not for 3' UTR or ncRNA binders |
| ASTERIX (Yeo) | Peak BED | Region | NA | Region distribution + region-specific motifs | Companion to RBP-Maps for non-splicing regions | Less granular than ChIPseeker |
| HOMER annotatePeaks.pl | Peak BED | Gene-level | Promoter > UTR > Exon > Intron > Intergenic | Annotation TSV | Familiar from ChIP-seq | Same TSS issue as ChIPseeker; less RNA-aware |
| bedtools intersect | Peak BED + feature BED | User-defined | User-defined | Per-feature overlap | Most flexible | Manual hierarchy logic must be written |
| pyranges (Python) | Peak BED + GTF | Customizable | User-defined | DataFrame | Fast, pythonic, scriptable | Hand-written priority logic |
| pybedtools | Peak BED + GTF | Customizable | User-defined | iterators | Flexible, scripted | Same as pyranges |
| Peakhood (Reuter 2022) | Peak BED + transcripts | Transcript context | NA | Per-peak transcript context | Resolves transcript-isoform ambiguity | Specialized; not a global region tool |

Methodology evolves; verify the GTF source (GENCODE preferred over Ensembl for human; both have feature differences) and the priority hierarchy against published RBP literature.

## Critical Choice: Annotation Hierarchy

RNA features are nested and overlapping. A 3' UTR peak is ALSO in an exon; an intronic peak near a splice site is ALSO in a transcribed-but-spliced region. Without a priority hierarchy the same peak gets counted multiple times.

**Standard CLIP hierarchy (most-specific first):**

1. **3' UTR** (specific to mature mRNA biology)
2. **5' UTR**
3. **CDS** (coding sequence)
4. **Exon** (catch-all for non-UTR exonic)
5. **Splice site** (within 50 nt of donor/acceptor; key for splicing factors)
6. **Intron**
7. **Promoter** (within 1 kb of TSS; usually not relevant for RNA-binding)
8. **5' flanking** / **3' flanking** (within 1 kb of gene)
9. **ncRNA** (lncRNA, snoRNA, snRNA, miRNA host)
10. **Repeat element** (Alu, LINE, LTR, SINE) - separate axis from above
11. **Intergenic**

ChIPseeker's default `level='gene'` reports the FIRST hit by priority. For CLIP, set `level='transcript'` and explicitly customize `tssRegion=c(-100, 100)` (default `c(-3000, 3000)` is for ChIP).

**Why CLIP needs the tight TSS window:** ChIPseeker was designed for ChIP-seq where the "Promoter" category captures transcription-factor binding within 3 kb upstream of a gene. For CLIP-seq, the binding substrate is RNA - the RBP cannot bind upstream of the TSS because there is no RNA upstream of the TSS. A 6 kb TSS window catches 5'-UTR peaks and unrelated upstream sequence, mislabeling 30-50% of CLIP peaks as "Promoter (2-3kb)" - implausible for an RNA-binding study. Tightening to `c(-100, 100)` constrains the "Promoter" category to the mRNA 5' end immediate vicinity and pushes deeper-5'-UTR peaks into the "5UTR" category where they belong.

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile('peaks.stringent.bed')

# CLIP-appropriate annotation
anno <- annotatePeak(
    peaks,
    TxDb = txdb,
    level = 'transcript',     # RNA-isoform-aware
    tssRegion = c(-100, 100), # tight; not a meaningful promoter for RNA binding
    genomicAnnotationPriority = c(
        'Promoter', '5UTR', '3UTR', 'Exon',
        'Intron', 'Downstream', 'Intergenic'
    )
)

print(anno)
plotAnnoPie(anno)
plotAnnoBar(anno)
plotDistToTSS(anno)  # CLIP-specific: dist to TSS is dist to mRNA 5' end
```

## RBP Class -> Expected Annotation

| RBP class | Expected dominant region | RBP examples |
|-----------|--------------------------|--------------|
| Splicing factors | Intron (within 500 nt of splice junction) | PTBP1, U2AF2, RBFOX1/2, SRSF1-9, HNRNPC, MBNL1 |
| 3' UTR mRNA stability/decay | 3' UTR | HuR (ELAVL1), PUM2, AUF1, ZFP36, TIA1, NUDT21 |
| 5' UTR translation regulation | 5' UTR, CDS start | EIF4A3, EIF3J, eIF2 subunits, LARP1, PCBP1/2 |
| Ribosome-associated | CDS | RPL/RPS proteins, RACK1 |
| Nuclear export | CDS, 3' UTR | NXF1, ALYREF, SRSF3 |
| m6A reader | 3' UTR, stop codon region | YTHDF1/2/3, IGF2BP1/2/3 |
| miRNA effector | 3' UTR (seed-matched sites) | AGO1/2/3/4, DICER1, TNRC6A/B/C |
| snoRNA-targeting RBP | snoRNA, rRNA | DKC1, FBL, NOP56, NHP2 |
| ncRNA-binding | Specific ncRNA target | TROVE2 (Y RNA), LARP7 (7SK), SSB (vault) |
| Mitochondrial | chrM (mt-mRNA + 12S/16S rRNA) | FASTKD2, LRPPRC, TFAM, MTPAP |
| Histone mRNA | Histone 3' end stem-loop | SLBP |
| Repeat-binding (TE) | Alu, LINE-1, LTR, SINE | MATR3, ZFP36, HNRNPK, HNRNPA2B1, FUS (LINE-1) |
| LCD/IDR / phase separation | Diverse (Alu, intron, 3' UTR) | FUS, TDP-43, EWSR1, TAF15 |

Mismatch between observed and expected suggests: (a) failed IP (cross-reactive antibody binding wrong RBP), (b) artifact (high-expression contamination), (c) genuinely new biology - investigate.

## RBP-Maps: Splicing Regulatory Maps

For splicing factors (PTBP1, U2AF2, RBFOX, SRSF1, etc.), the Yeo lab RBP-Maps tool aggregates eCLIP signal in a 350 nt window flanking the upstream, cassette, and downstream exons, producing a 1400 nt vectorized regulatory map. This map shows the position-dependent enrichment of RBP binding relative to alternative splicing events and is THE standard for splicing-CLIP analysis.

```bash
# Yeo lab rbp-maps Snakemake workflow
git clone https://github.com/YeoLab/rbp-maps
cd rbp-maps

# Inputs:
#   1. eCLIP BAM (rep1 + rep2)
#   2. Cassette exon BED (from RNA-seq KD experiment; see clip-seq/differential-clip)
#   3. Background "native" exons (constitutively spliced)

python rbpmaps/RBPMaps.py \
    --ip rep1.bam rep2.bam \
    --input sminput.bam \
    --positive cassette_inc.bed \
    --negative cassette_exc.bed \
    --background native_constitutive.bed \
    --output rbpmaps_out/ \
    --normalization input_normalization
```

Output: 1400 nt vector with eCLIP signal at every base, plotted as a metagene for cassette inclusion vs exclusion vs constitutive. A peak in the upstream-intron region (right of the cassette exon's 5' splice site) indicates intronic regulation; a peak in the cassette exon body indicates exonic regulation.

## Per-Tool Failure Modes

### ChIPseeker -- TSS-region default over-extends

**Trigger:** Used ChIPseeker default `tssRegion=c(-3000, 3000)` for CLIP peaks.

**Mechanism:** The 6 kb TSS window catches CLIP peaks deep in the 5' UTR or upstream introns and labels them "Promoter (2-3kb)". For RNA biology, the "promoter" category is meaningless - the RBP cannot bind pre-mRNA upstream of the TSS.

**Symptom:** 30-50% of peaks labeled "Promoter (2-3kb)" - implausible for an RNA-binding study.

**Fix:** Set `tssRegion=c(-100, 100)` or `c(-50, 50)`. Remove the "Promoter" category from the priority vector entirely if it gives 0 hits.

### ChIPseeker -- Gene-level loses isoform context

**Trigger:** Used `level='gene'` (default) on peaks that should be transcript-specific.

**Mechanism:** Gene-level annotation picks the canonical/longest transcript and assigns the peak to that transcript's features. A peak in an alternatively-spliced exon may be labeled "intron" because the gene-level reference picks an isoform where the exon is skipped.

**Symptom:** Splicing factor peaks heavily "intronic" when they should be exonic at alternative-cassette exons.

**Fix:** Use `level='transcript'`. Match isoforms to RBP-Maps cassette-exon coordinates for splicing factors.

### RCAS -- Slow on large peak sets

**Trigger:** RCAS report on > 100k peaks.

**Mechanism:** RCAS generates an HTML report with many plots; rendering is slow.

**Symptom:** Run time > 1 h; out-of-memory; HTML failing to render.

**Fix:** Sub-sample peaks to top 10k; or use ChIPseeker for the global view and reserve RCAS for per-region deep dives.

### RBP-Maps -- Requires SE event table

**Trigger:** Splicing-factor CLIP with no companion RNA-seq KD or no SE event BED.

**Mechanism:** RBP-Maps's regulatory metagene only makes sense relative to cassette exons differentially regulated by the RBP. Without the SE event table, the tool falls back to plotting at all native exons (no differential signal).

**Symptom:** Flat metagene; no position-dependent enrichment.

**Fix:** Generate the cassette-exon table from companion RNA-seq KD experiments (rMATS / MAJIQ / LeafCutter; see alternative-splicing skills) OR use a published table from ENCODE shRNA RNA-seq.

### bedtools -- Manual hierarchy errors

**Trigger:** Hand-written `bedtools intersect` chain for region annotation.

**Mechanism:** Forgetting to enforce priority leads to double-counting: a 3' UTR peak intersects 3UTR.bed AND exon.bed.

**Symptom:** Per-region counts sum to > total peak count by 20-50%.

**Fix:** Process in priority order; remove already-annotated peaks from the remaining set with `-v` before next intersect. Or use ChIPseeker / RCAS which encode the hierarchy internally.

### Repeat-element overlap missed

**Trigger:** Repeat-binding RBP (MATR3, HNRNPK) annotated without RepeatMasker BED.

**Mechanism:** Genomic feature BEDs (UTR/CDS/intron) do not flag repeat instances within them. An MATR3 peak in an intronic LINE-1 is labeled "intron" when "LINE-1 within intron" is the biology.

**Symptom:** Repeat-binding RBP looks like an intronic binder; functional analysis (GO terms) returns generic "RNA processing".

**Fix:** Pass RepeatMasker BED (UCSC) as a separate axis; annotate each peak with BOTH region AND repeat-class. Cross-check repeat-class fraction against literature: > 15% Alu / LINE / LTR overlap indicates a repeat binder.

### Mitochondrial peaks excluded by default

**Trigger:** Standard TxDb annotation drops chrM transcripts; FASTKD2 / LRPPRC peaks all "intergenic".

**Mechanism:** Some TxDb objects (older releases) exclude mitochondrial chromosome; ChIPseeker reports those peaks as having no annotation.

**Symptom:** chrM peaks labeled "Distal Intergenic" or simply missing from the annotation table.

**Fix:** Verify TxDb includes chrM (`genes(txdb, filter=list(tx_chrom='chrM'))`); add custom mt-mRNA + mt-rRNA + mt-tRNA BED if missing.

### Strand information ignored

**Trigger:** `bedtools intersect -wa -wb` without `-s` flag.

**Mechanism:** CLIP is strand-specific. Without `-s`, peaks on one strand can be annotated to features on the other strand - especially for overlapping transcripts (sense/antisense pairs).

**Symptom:** Antisense transcript peaks contaminate the sense annotation; per-region distribution shifted.

**Fix:** Always pass `-s` to bedtools intersect for CLIP. ChIPseeker handles strand automatically when peaks have BED column 6.

## Metagene Analysis

Metagene = aggregated profile of CLIP signal at each base of a normalized feature (e.g., 5' UTR, CDS, 3' UTR scaled to common length). Reveals position-dependent binding.

**Goal:** Generate a position-resolved CLIP-signal profile aggregated across all 3' UTRs (or other features) to identify position-dependent binding modes (e.g., stop-codon-proximal m6A readers, polyA-proximal CPEB).

**Approach:** Convert dedup BAM to strand-stranded bigWig, then use deepTools `computeMatrix scale-regions` to align signal across feature instances of variable length, and `plotProfile` / `plotHeatmap` for visualization.

```bash
# deepTools computeMatrix for metagene
# Step 1: peaks BED -> bedgraph
bedtools genomecov -bg -strand + -ibam dedup.bam > clip_plus.bg
bedtools genomecov -bg -strand - -ibam dedup.bam > clip_minus.bg
# combine to bigwig
bedGraphToBigWig clip_plus.bg chrom.sizes clip_plus.bw

# Step 2: metagene over 3' UTRs
computeMatrix scale-regions \
    -S clip_plus.bw clip_minus.bw \
    -R three_prime_UTRs.bed \
    --regionBodyLength 500 \
    --beforeRegionStartLength 100 \
    --afterRegionStartLength 100 \
    -o matrix.gz

plotProfile -m matrix.gz -o metagene.png --perGroup
plotHeatmap -m matrix.gz -o metagene_heatmap.png
```

| Metagene shape | Biology |
|----------------|---------|
| Peak at stop codon region | m6A reader (YTHDF1/2/3); IGF2BP1; some 3' UTR factors |
| Peak in proximal 3' UTR (within 200 nt of stop) | PUM2, IGF2BP, miRNA effectors |
| Peak in distal 3' UTR (within 200 nt of poly-A) | NUDT21, CPEB1, polyadenylation regulators |
| Peak at 5' splice site flanking intron | RBFOX, MBNL, HNRNPC |
| Peak at 3' splice site / branch point | U2AF2, PTBP1 |
| Peak at 5' UTR (within 50 nt of TSS) | EIF3, EIF4A3, LARP1 |
| Flat across CDS | Cytoplasmic ribosome-associated; less position-specific |

## Decision Tree by Scenario

| Scenario | Tool + parameters | Why |
|----------|-------------------|-----|
| Global region distribution (3' UTR / CDS / intron) | ChIPseeker `tssRegion=c(-100,100) level='transcript'` | Mature, fast, ENCODE-comparable |
| Splicing factor regulatory map | RBP-Maps (Yeo) with RNA-seq KD cassette table | The standard for splicing-CLIP |
| RNA-feature aware (snoRNA, lncRNA-specific) | RCAS or custom GTF + bedtools | RCAS handles ncRNA features |
| Transcript-isoform context resolved | Peakhood + ChIPseeker level='transcript' | Resolves which isoform the peak supports |
| Metagene 3' UTR / 5' UTR | deepTools computeMatrix + plotProfile | Standard metagene visualization |
| Repeat-element overlap (MATR3, HNRNPK) | bedtools intersect with RepeatMasker | Repeat axis is independent of region axis |
| Mitochondrial RBP (FASTKD2) | Custom chrM-aware annotation | Standard TxDb may drop chrM |
| Bulk RBP overview | ChIPseeker plotAnnoPie + plotDistToTSS | Quick orienting plots |
| Compare two RBPs' region distributions | ChIPseeker annotatePeakList + side-by-side bar | Direct comparability |
| miRNA target site (AGO-CLIP) | bedtools intersect with 3' UTR + TargetScan seed scan | miRNA biology specific |
| m6A reader (YTHDF, IGF2BP) | ChIPseeker + custom stop-codon-distance plot | m6A biology stops-codon-proximal |

## Reconciliation: When Annotations Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ChIPseeker "Promoter" >> 5' UTR for CLIP | TSS region too wide | Tighten `tssRegion=c(-100,100)` |
| ChIPseeker gene-level intronic; RBP-Maps shows exonic at cassettes | Gene-level picks canonical isoform; RBP-Maps uses cassettes | Trust RBP-Maps for splicing factors |
| Many peaks "intergenic" near a known gene | TxDb missing the transcript | Verify GENCODE vs Ensembl GTF; add lincRNA / unannotated transcripts |
| chrM peaks missing | TxDb excluded chrM | Add custom mt-mRNA BED |
| Sum of per-region counts > total peaks | Hierarchy not enforced | Use ChIPseeker; or `bedtools intersect -v` chain |
| Repeat-binding RBP looks "intronic" | RepeatMasker axis not added | Annotate peaks with RepMask intersect separately |
| Metagene flat for known position-specific RBP | Wrong feature BED (e.g., gene-body instead of transcript-isoform) | Use isoform-aware feature BED |
| RBP-Maps signal flat at cassette exons | No RNA-seq KD reference table | Generate cassettes from companion KD; use ENCODE shRNA tables |

**Operational rule:** Report (a) ChIPseeker global region pie, (b) RBP-Maps metagene for splicing factors, (c) repeat-element fraction separately, (d) transcript-isoform context with Peakhood if ambiguous. Three views together prevent misinterpretation.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "Promoter" >> 30% of peaks in CLIP | TSS region default too wide | `tssRegion=c(-100,100)` |
| Strand-mixing artifacts | Bedtools intersect without `-s` | Always pass `-s` for CLIP |
| 3' UTR fraction < 10% for HuR | Wrong txdb version (UTR poorly annotated) | Use GENCODE v38+ TxDb |
| chrM peaks lost | TxDb excludes mitochondrial | Add mt-transcript BED |
| RBP-Maps flat metagene | No cassette exon BED | Generate from RNA-seq KD; or use ENCODE table |
| Repeat fraction over-counted | Same repeat instance counted multiple times | `bedtools merge` repeat BED first |
| ChIPseeker `level='gene'` collapses isoforms | Default level=gene | Set `level='transcript'` |
| Annotation runs slow for large peak set | RCAS HTML generation | Use ChIPseeker for global; RCAS for top-K |
| GO enrichment yields generic terms | Repeat-binding RBP annotated as intronic | Re-annotate with RepeatMasker axis |
| Peak in unannotated lincRNA labeled intergenic | Older GENCODE missing the lincRNA | Update GENCODE; add manual transcript BED |

## References

- Yu G et al 2015 Bioinformatics 31:2382 (ChIPseeker)
- Yusuf D et al 2017 BMC Bioinformatics 18:227 (RCAS)
- Yee BA et al 2019 RNA 25:193 (RBP-Maps splicing regulatory maps)
- Van Nostrand EL et al 2020 Nature 583:711 (ENCODE 150 RBP eCLIP, region distribution analysis)
- Quinlan AR & Hall IM 2010 Bioinformatics 26:841 (bedtools)
- Hentze MW et al 2018 Nat Rev Mol Cell Biol 19:327 (RBP function review)
- Reuter JS et al 2022 Bioinformatics 38:1139 (Peakhood transcript-context)
- ENCODE eCLIP standards (encodeproject.org/eclip) - region distribution conventions

## Related Skills

- clip-seq/clip-peak-calling - Upstream peak calls
- clip-seq/crosslink-site-detection - Single-nt CL sites for fine-grained metagene
- clip-seq/clip-motif-analysis - Motif discovery within annotated regions
- clip-seq/differential-clip - Cross-condition annotated peaks for regulatory maps
- clip-seq/ago-clip-mirna-targets - 3' UTR seed-matched site annotation
- clip-seq/m6a-clip - Stop-codon-proximal m6A annotation
- genome-intervals/gtf-gff-handling - GTF preparation for annotation
- genome-intervals/interval-arithmetic - bedtools intersect patterns
- alternative-splicing/differential-splicing - Cassette exon tables for RBP-Maps
- chip-seq/peak-annotation - DNA-protein annotation analogue
