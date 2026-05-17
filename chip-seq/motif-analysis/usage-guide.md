# Motif Analysis - Usage Guide

## Overview

Discovers de novo motifs and tests known motif enrichment in ChIP-seq, ATAC-seq, or other peak sequences using HOMER, MEME-ChIP (STREME, CentriMo, TOMTOM, FIMO), monaLisa (regression-based binned analysis), and AME (differential enrichment). Embeds background-selection theory (GC-matched, dinucleotide-shuffled, Markov order-2, peak-flanks), summit-centering requirements, and recognizes when to use deep-learning attribution methods (TF-MoDISco on BPNet/chromBPNet) for soft motif syntax.

## Prerequisites

```bash
# HOMER
conda install -c bioconda homer
perl $CONDA_PREFIX/share/homer/configureHomer.pl -install hg38 mm10

# MEME suite (5.4+ has STREME, not DREME)
conda install -c bioconda meme

# Motif databases
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
wget https://hocomoco12.autosome.org/final_bundle/hocomoco12/H12CORE/H12CORE_meme_format.meme
```

```r
# Optional: monaLisa for regression-based binned analysis
BiocManager::install(c('monaLisa', 'JASPAR2024', 'BSgenome.Hsapiens.UCSC.hg38'))
```

## Quick Start

Tell the agent what to do:
- "Run HOMER findMotifsGenome.pl on summit-centered narrowPeak with repeat masking and 8 cores"
- "Run MEME-ChIP with JASPAR 2024 motifs on peaks resized to ±100 bp around summit"
- "Find motifs central to my CTCF peaks (test direct vs tethered binding via CentriMo)"
- "Compare motif enrichment between gained-in-treatment and gained-in-control peaks (AME differential)"
- "Use monaLisa to find motifs discriminating top-quintile log2FC peaks from bottom-quintile"
- "Scan my motif of interest genome-wide with FIMO at p < 1e-5"
- "My top motif is the AluY consensus; diagnose and fix the background"

## Example Prompts

### Standard discovery (TF or sharp histone)
> "Run HOMER de novo + known motif discovery on FOXA1 ChIP-seq peaks. Recenter to ±100 bp around summit, mask repeats, use 8 cores. Report top 10 known motifs and top 5 de novo."

### MEME-ChIP comprehensive
> "Run MEME-ChIP on H3K27ac peaks with JASPAR 2024 vertebrates database. Include STREME for short motifs, MEME for longer motifs, CentriMo for central enrichment, TOMTOM for known-motif matching."

### Differential motifs
> "I have peaks gained in EZH2-inhibitor treatment vs DMSO control. Run AME to find motifs enriched in gained-only peaks."

### Binned monaLisa
> "Split my peaks into quintiles by log2FC from differential analysis. Run monaLisa to find which TF motifs discriminate top quintile from bottom, controlling for GC content."

### Central enrichment test
> "Run CentriMo on my CTCF peaks to verify direct sequence-specific binding (motif should be centrally enriched within ±50 bp of summit) vs tethered/indirect binding (no central enrichment)."

### FIMO scanning
> "Scan the human genome (hg38) for FOXA1 motif occurrences at p < 1e-5 with HOCOMOCO v12 PWM."

### Diagnosing artifacts
> "My top motif is GC-rich, doesn't match any TF in TOMTOM, and looks like an Alu fragment. Check for repeat contamination and re-run with `-mask` plus GC-matched background."

### Deep-learning motifs
> "Extract motifs from a trained chromBPNet model using TF-MoDISco from the attribution scores."

## What the Agent Will Do

1. **Validate peak input**: confirm peaks pass FRiP and aren't dominated by hyper-ChIPable artifacts (see chipseq-qc)
2. **Recenter peaks on summit**: narrowPeak column 10 gives summit offset; resize to ±100-250 bp (sharp marks) or use full width (variable marks)
3. **Extract peak sequences**: `bedtools getfasta` with the matching genome FASTA
4. **Choose background appropriately**: GC-matched genomic (HOMER default with `-mask`), Markov order-2 (STREME / MEME-ChIP default), peak-flanks (when peak GC differs from genome), or differential set (AME)
5. **Run discovery + enrichment:**
   - HOMER: de novo + known in one command
   - MEME-ChIP: STREME + MEME + CentriMo + TOMTOM + FIMO
   - monaLisa: binned regression with GC control
   - AME: differential between two sequence sets
6. **TOMTOM motif comparison**: verify de novo motifs against JASPAR / HOCOMOCO / known TF databases
7. **Verify against expected biology**: known TF should rank highly for validated ChIP; failure to recover known motif is a QC failure indicator
8. **Output**: HTML report (HOMER / MEME-ChIP) + ranked motif table + sequence logos
9. **Document**: exact peak input (count, width, summit-centered or not), background choice, motif DB version, significance threshold

## Decision Tree: Which Tool

| Goal | Tool |
|------|------|
| Quick discovery + known + report | HOMER findMotifsGenome.pl |
| Comprehensive (de novo + central + DB + scan) | MEME-ChIP |
| Long or gapped motifs | MEME (classical, 5.4+ still works) |
| Short motifs (3-15 bp) at scale | STREME |
| Central enrichment of known motifs | CentriMo |
| Differential between two peak sets | AME |
| Binned regression (correlate with continuous variable) | monaLisa |
| Genome-wide motif scanning | FIMO (tighten p ≤ 1e-5) |
| Soft motif syntax / cooperativity | TF-MoDISco on BPNet/chromBPNet attribution (see chip-deep-learning) |
| Position-specific motifs (TF dimer arrangements) | STREME + CentriMo + ME-Aria |

## Tips

- **Always summit-center sequences.** Motif enrichment improves dramatically when sequences are 200 bp around summit vs full peak width. CentriMo requires summit centering.
- **Mask repeats unless explicitly studying TE-derived enhancers.** Repeat-derived motifs (AluY, LINE, LTR) are common false positives.
- **Use Markov order-2 background for de novo** (STREME/MEME default). Order-0 (mononucleotide) misses dinucleotide biases like CpG.
- **HOMER's auto background can include peaks.** For large peak sets (>5% of genome), supply explicit `-bg`.
- **FIMO genome-wide at default p ≤ 1e-4 = millions of false positives.** Tighten to 1e-5 for whole-genome scans; restrict to peaks for fine-grained p.
- **Differential motif analysis is more interpretable than absolute enrichment.** Compare gained-in-A vs gained-in-B (AME) rather than each set vs genome.
- **monaLisa for continuous variables.** With a per-peak measurement like log2FC, accessibility, or methylation, binned regression is more powerful than absolute enrichment.
- **CentriMo central enrichment is the gold standard for direct binding.** Tethered (indirect) binding shows motif enrichment but not central enrichment.
- **Failure to recover the known TF motif for a validated ChIP is a QC failure.** Re-check antibody validation, peak quality, and background choice before claiming "novel motif."
- **For broad histone marks (H3K27me3, H3K9me3), motif analysis is generally not informative**: Polycomb and heterochromatin do not depend on sequence-specific binding directly.

## Troubleshooting

### No significant motifs

Causes in order of frequency:
1. **Peak quality / hyper-ChIPable contamination** → see chipseq-qc, filter blacklist + custom hyper-ChIPable
2. **Too few peaks** (< 500) → motif discovery needs ≥ 500-1000 high-confidence peaks
3. **Wrong sequence window** → recenter to summit ±100-250 bp
4. **Wrong background** → for unusual GC composition, supply matched background explicitly
5. **Broad histone mark** → motif analysis may not apply

### Top motif is GC-rich or matches no TF database

1. Repeat contamination → use `-mask` (HOMER) or pre-mask sequences
2. CpG island bias → matched background; consider order-2 Markov
3. Hyper-ChIPable artifacts → filter peaks; see chipseq-qc
4. Compositional control failed → switch tool (STREME order-2 typically better)

### Known TF motif not recovered

1. ChIP quality issue → check FRiP, NSC, antibody validation
2. Tethered (indirect) binding → CentriMo may show no central enrichment despite overall enrichment
3. Cofactor binding (motif belongs to partner) → check known interaction partners

### Too many motif hits (FIMO)

Default `--thresh 1e-4` generates millions of false positives genome-wide. Use `--thresh 1e-5` for whole-genome OR restrict to peaks `fimo motif.meme peaks.fa`.

### MEME-ChIP "sequences too short"

MEME requires sequences ≥ motif min width. Resize peaks to ≥ 200 bp.

### MEME-ChIP / STREME memory failure

Long sequences (> 1 kb each) at large counts exceed RAM. Resize to ±100-250 bp around summit; downsample peak count if needed.

## Related Skills

- chip-seq/peak-calling - Generate input peaks
- chip-seq/chipseq-qc - Filter hyper-ChIPable peaks before motif analysis
- chip-seq/chip-deep-learning - TF-MoDISco on BPNet/chromBPNet for sequence syntax
- chip-seq/peak-annotation - Separate promoter vs enhancer peaks for stratified motif analysis
- atac-seq/motif-deviation - chromVAR per-cell motif activity (ATAC-specific)
- atac-seq/footprinting - TOBIAS footprint analysis
- sequence-manipulation/motif-search - General motif scanning
- genome-intervals/proximity-operations - Extract peak FASTA via bedtools getfasta
