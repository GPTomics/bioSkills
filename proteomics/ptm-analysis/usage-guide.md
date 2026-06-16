# PTM Analysis - Usage Guide

## Overview
Localize and quantify post-translational modifications (phosphorylation, acetylation, ubiquitination, glycosylation) as three stacked inferences: which sites the enrichment chemistry actually captured, where on the peptide the modification sits (with its own false localization rate), and whether an abundance change survives subtracting the protein-level change. The central decision is never trusting a phospho-only fold-change as "regulation" without protein-level adjustment.

## Prerequisites
```bash
pip install pandas numpy scipy
# R packages: BiocManager::install(c("MSstatsPTM"))
# Kinase activity: KSEAapp (CRAN) or PTM-SEA / ssGSEA2.0 (GitHub broadinstitute/ssGSEA2.0)
# CLI search engines (upstream): MaxQuant, FragPipe/MSFragger, DIA-NN, Spectronaut
```

## Quick Start
Tell your AI agent what you want to do:
- "Load my MaxQuant Phospho (STY)Sites.txt, expand multiplicity, and keep class I sites"
- "Adjust phosphosite changes for protein abundance using MSstatsPTM and a paired global proteome"
- "Build a kinase-motif logo using an experiment-matched background, not the whole proteome"
- "Infer which kinases are active with KSEA from my site fold-changes"
- "Check whether my diGly sites are confounded by NEDD8/ISG15 or an iodoacetamide artifact"

## Example Prompts

### Site Identification and Localization
> "Filter Phospho (STY)Sites.txt to class I (localization probability >= 0.75) and report the residue distribution"

> "Expand the MaxQuant site-table multiplicity into Intensity___1/___2/___3 before any quantification"

> "Estimate an empirical global false localization rate for my phosphosites"

### Protein-Adjusted Quantification
> "Use MSstatsPTM groupComparisonPTM and call only ADJUSTED.Model hits as regulated"

> "Show me which apparent site changes are actually driven by protein abundance"

> "Compare phosphosite changes after drug treatment, adjusting for protein-level stabilization"

### Other PTMs
> "Treat my K-GG data as ubiquitin plus NEDD8 plus ISG15 and flag the chemistry confounds"

> "Map acetylation sites allowing four or more missed cleavages because acetyl-K blocks trypsin"

> "Disambiguate glycosite N->D from spontaneous deamidation"

### Motif and Kinase Activity
> "Run motif analysis with an experiment-matched S/T/Y background and render a sequence logo"

> "Run KSEA or PTM-SEA to infer active kinases and report z-scores with substrate counts"

## What the Agent Will Do
1. Frame the question as enrichment -> localization -> quantification and identify which layer the user is asking about
2. Load the search-engine site output, drop Reverse/contaminant, and filter localization probability to class I
3. Expand MaxQuant multiplicity to a long, multiplicity-resolved site matrix before any statistics
4. Run MSstatsPTM with a paired global proteome and report PTM, PROTEIN, and ADJUSTED models, calling only ADJUSTED hits regulated
5. Perform motif analysis against an experiment-matched background and kinase-activity inference with a curated prior
6. Report three numbers (peptide FDR, localization probability, global FLR) and triage hits by functional evidence

## Tips
- Always acquire a paired global (unenriched) proteome; without it, "regulated site" claims are unfalsifiable.
- A between-method or between-lab difference is a chemistry hypothesis first (TiO2 vs Fe-IMAC mono/multi bias), biology second.
- Quantify on multiplicity columns (three underscores), never the collapsed base Intensity.
- For ubiquitinomes, confirm chloroacetamide alkylation and remember K-GG is not ubiquitin-specific.
- Build motif backgrounds from the matched dataset; whole-proteome backgrounds just rediscover disordered-region composition bias.
- Kinase-activity inference is prior-limited; simple z-score matches sophisticated methods, so invest in the substrate prior.
- Most sites have no annotated function; triage by conservation, stoichiometry, and Ochoa functional score before claiming biology.

## Related Skills
- peptide-identification - Identify modified peptides and run open/variable-mod search
- quantification - Underlying protein-level quant feeding the MSstatsPTM PROTEIN dataset
- differential-abundance - Moderated testing on the protein-level intensity matrix
- pathway-analysis/gsea - Enrichment scoring of regulated-site protein lists and PTM-SEA-style signatures
- data-visualization/sequence-logos - Render motif logos from the foreground/background windows
