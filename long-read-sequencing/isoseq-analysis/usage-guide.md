# Full-Length Isoform Analysis - Usage Guide

## Overview
This skill discovers, classifies, filters, and quantifies full-length transcript isoforms from PacBio Iso-Seq/Kinnex (HiFi) and Oxford Nanopore (cDNA and direct-RNA) long reads. The central principle is that a long-read novel isoform is an artifact until proven otherwise: RT template-switching, intra-priming on genomic poly-A, and 5' RNA degradation manufacture novel junctions and truncated isoforms, so the SQANTI3/pigeon classification plus artifact filter plus orthogonal end/junction support (CAGE for 5' TSS, poly-A atlas for 3' TES, short-read STAR junctions) is the analysis, not a postscript. It covers the PacBio isoseq+pigeon pipeline (including the Kinnex skera-split step), SQANTI3 for any platform, and the ONT tools (IsoQuant, FLAIR, Bambu, StringTie2, FLAMES). Differential isoform usage is handed off to the alternative-splicing category.

## Prerequisites
```bash
# PacBio path
conda install -c bioconda isoseq pbmm2 pigeon lima
# Cross-platform curation
conda install -c bioconda sqanti3
# ONT path
conda install -c bioconda isoquant flair bambu stringtie pychopper minimap2
# Reference annotation + genome, plus orthogonal support: CAGE refTSS BED, poly-A motif list,
# short-read STAR SJ.out.tab for junction validation
```

## Quick Start
Tell your AI agent what you want to do:
- "Run the PacBio Iso-Seq pipeline and classify the isoforms with pigeon"
- "Deconcatenate my Kinnex data and build a filtered isoform catalog"
- "Discover and quantify isoforms from my ONT cDNA with IsoQuant"
- "Tell me which of my novel isoforms are trustworthy"

## Example Prompts

### PacBio Iso-Seq / Kinnex
> "I have Kinnex HiFi reads. Run skera, lima, refine, collapse, then classify and filter with pigeon using my CAGE peaks and poly-A motifs, and report the FSM/ISM ratio and saturation."

### ONT discovery and quantification
> "These are ONT direct-RNA reads. Spliced-align them correctly, discover and quantify isoforms with IsoQuant, then curate with SQANTI3 and short-read junctions."

### Artifact triage
> "My catalog has 40% novel isoforms. Tell me whether that is real biology or an artifact problem, and which classes (ISM, NNC, mono-exon) to distrust."

### Validation
> "Validate my novel isoforms against CAGE 5' peaks, poly-A sites, and short-read junctions, and drop the intra-priming and RT-switching artifacts."

## What the Agent Will Do
1. Identify the platform and chemistry (PacBio Iso-Seq/Kinnex vs ONT cDNA vs ONT direct-RNA) and choose the pipeline.
2. For Kinnex, deconcatenate with skera first; for ONT cDNA, orient with pychopper.
3. Build the isoform set (isoseq collapse, or IsoQuant/FLAIR/Bambu for ONT).
4. Classify against the reference (pigeon or SQANTI3) and apply the artifact filter with CAGE/poly-A/short-read support.
5. Report the FSM/ISM ratio, novel-class breakdown, and a saturation curve excluding singletons.
6. Quantify with an EM tool and hand differential isoform usage to alternative-splicing.

## Tips
- A high novel-isoform fraction is a red flag, not a success - it usually means an under-powered filter or degraded RNA.
- ISM fraction is an RNA-integrity thermometer; do not report ISMs as novel without CAGE 5' support.
- NIC (novel combination of known junctions) is far more trustworthy than NNC (novel splice site); mono-exon novels are the false-discovery sink.
- The Iso-Seq binary is `isoseq` (renamed from `isoseq3` in v4); `pigeon` is the classifier/filter, not a quantifier, and it consumes the collapsed.sorted.gff after `pigeon prepare`, not a BAM.
- Kinnex/MAS-seq data must be deconcatenated with `skera split` before lima.
- Match the minimap2 preset to chemistry; `-uf` is for stranded direct-RNA/Iso-Seq, not unoriented ONT cDNA.
- Long-read isoform quantification needs EM (Bambu/IsoQuant/NanoCount); raw FLNC counts are not short-read-equivalent abundances, and discovery is depth-unsaturated.

## Related Skills

- long-read-alignment - Spliced alignment presets for cDNA/direct-RNA
- basecalling - Direct-RNA chemistry and cDNA vs direct-RNA tradeoffs
- alternative-splicing/long-read-splicing - Long-read splicing (define the boundary)
- alternative-splicing/isoform-switching - Differential isoform usage downstream
- rna-quantification/tximport-workflow - Transcript-level quantification downstream
- genome-annotation/eukaryotic-gene-prediction - Isoforms as annotation evidence

## Resources
- [PacBio Iso-Seq docs](https://isoseq.how/)
- [SQANTI3](https://github.com/ConesaLab/SQANTI3)
- [IsoQuant](https://github.com/ablab/IsoQuant)
- [pychopper](https://github.com/epi2me-labs/pychopper)
