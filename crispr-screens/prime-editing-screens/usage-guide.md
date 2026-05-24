# Prime-Editing Screens - Usage Guide

## Overview

Decision-grade design and analysis of pooled prime-editor screens. Covers pegRNA design with PRIDICT (Mathis 2023 Nat Biotechnol 41:1151) and PRIDICT2 (Mathis 2025 Nat Biotechnol 43:712, online June 2024); PE chemistry selection (PE2, PE3, PEmax, PEAR, PE5max); pegRNA architecture (spacer + scaffold + PBS + RTT); chromatin context as the primary determinant; PRIME pooled-screen methodology (Erwood 2023); MOSAIC in situ saturation mutagenesis; CRISPResso2 quantification; and cross-modal validation with base-editor screens.

## Prerequisites

```bash
# PRIDICT2 (web + CLI)
git clone https://github.com/uzh-dqbm-cmi/PRIDICT2
# OR use web interface https://pridict.it/

# CRISPResso2 for amplicon analysis
pip install CRISPResso2

# Helpers
pip install pandas numpy biopython
```

Required inputs:
- Intended variant list with genomic coordinates
- Target genome reference (e.g. GRCh38)
- Cell line for PE chemistry validation
- Chromatin accessibility profile (ATAC-seq) for chromatin-aware filtering

## Quick Start

Tell the AI agent what to do:
- "Design a pegRNA library for 320 specific ClinVar MLH1, MSH2, MSH6, PMS2 variants. Use PRIDICT2 to score; select top 3 pegRNAs per variant"
- "Compare BE vs PE for installing a C->T variant at TP53 R273H. Recommend chemistry based on editing window + bystander"
- "Diagnose: my PE screen has 80% scaffold incorporation. Re-design pegRNAs with longer RTT and verify with PRIDICT2"
- "Run CRISPResso2 in prime-editor mode on amplicon FASTQ. Report intended-edit %, scaffold incorporation, indel byproduct per pegRNA"
- "Cross-validate PE screen hits with parallel BE screen at the same variants; identify PE-only hits in transversion variants"

## Example Prompts

### pegRNA Library Design

> "Build pegRNA library for 500 ClinVar MMR gene variants (MLH1, MSH2, MSH6, PMS2). Score each candidate pegRNA with PRIDICT2. Select 3 best pegRNAs per variant with predicted efficiency >50%. Output peg_library.csv with spacer, PBS, RTT, predicted_efficiency."

> "For installing MLH1 c.677A>G as a single intended variant: PE (no bystander) vs BE (no PAM-restriction). Run PRIDICT2 on PE candidates; check editing window for BE; recommend best chemistry."

### Screen Analysis

> "Run CRISPResso2 in prime-editor mode on my 50-pegRNA pilot. Report intended_edit_pct, scaffold_incorp_pct, indel_pct per pegRNA. Flag pegRNAs with intended_edit <5% for filtering or re-design."

> "Apply MAGeCK MLE on my PE screen counts after filtering to efficiency-passing pegRNAs. Aggregate per-pegRNA LFCs to per-variant scores."

### Chromatin-Aware Design

> "Cross-reference my pegRNA library with K562 ATAC-seq accessibility. Flag pegRNAs targeting closed chromatin (signal <0.5) for empirical pilot validation before library order. Mathis 2024 documented chromatin as the dominant determinant beyond PRIDICT sequence prediction."

### MOSAIC / Saturation

> "Design MOSAIC saturation library tiling TP53 DNA-binding domain (amino acids 94-294). For each residue, design pegRNAs covering all 19 possible amino acid changes (or at least all functionally relevant). Include variants from ClinVar."

### Cross-Validation

> "I have parallel BE and PE screens at the same 200 variants. Intersect hits at FDR <0.05; identify high-confidence variants (both methods, same direction). Flag PE-only hits in non-BE-coverable variants (transversions) as genuine PE-unique."

### Diagnostics

> "PE screen has 30% scaffold incorporation across library. Diagnose: RTT too short relative to PBS, or RT processivity issue with this cell line? Recommend pegRNA re-design or chemistry switch."

> "PRIDICT2 predicts 65% efficiency for my top pegRNAs, but pilot CRISPResso2 shows 8%. Identify chromatin context as the cause; cross-reference ATAC-seq."

## What the Agent Will Do

1. Confirm intended variants are PE-installable (NGG PAM within 30 nt; achievable transversion/multi-base if BE-uncoverable)
2. Decide chemistry: PE2 standard; PEmax for high-efficiency; PE3 with second nick if expression supports
3. Design pegRNA candidates: spacer + PBS (11-13 nt, 40-55% GC) + RTT (10-20 nt) + scaffold
4. Run PRIDICT2 to predict per-pegRNA efficiency, indel rate, scaffold incorporation
5. Filter library to top 3 pegRNAs per variant with predicted efficiency >50%
6. Cross-reference with ATAC-seq accessibility; flag closed-chromatin loci for pilot
7. Pilot the library at 20-50 representative loci; verify empirical efficiency matches PRIDICT2
8. Scale to full screen at MOI 0.3 in PE-expressing cell line
9. Endpoint amplicon sequencing of each target locus
10. CRISPResso2 quantification: intended-edit %, scaffold incorp %, indel %
11. Filter library to >5% intended-edit pegRNAs
12. Per-variant hit calling with MAGeCK MLE / drugZ
13. Cross-validate with BE screen if applicable
14. Output: per-variant fitness, BE-PE concordance, validation strategy

## Tips

- PRIDICT2 sequence prediction is a starting point, not a guarantee. Chromatin context is the dominant determinant beyond sequence; always pilot at representative sites.
- For high-stakes variant-function calls, validate with both BE (where applicable) and PE. Concordant hits are high confidence; PE-only at BE-coverable variants are suspect.
- PE is the gold standard for transversions (C->G, C->A, G->C, G->A) where BE doesn't work; for these variants, PE is the only choice.
- Scaffold incorporation is a design issue, not a chemistry issue. Always check post-PE; if >5%, re-design pegRNAs with longer RTT.
- PE3 adds a nicked sgRNA on the opposite strand; this increases editing efficiency but doubles the indel risk. Use PE2 unless PEmax + PRIDICT2 doesn't reach target efficiency.
- For genome-scale saturation editing (BRCA1 / MMR / TP53), MOSAIC is the canonical methodology. 1000s of variants per cell line; smaller per-variant cell counts; suitable for drug-resistance variant scanning.
- Pilot PRIDICT2 predictions in your specific cell line. The deep-learning model was trained on HEK293T / HCT116; other lines may have different chromatin and RT-activity profiles.
- For BE-PE comparison, design parallel libraries targeting the same variants. Concordant hits represent strong biology.

## Decision Cheat Sheet

| Question | Answer |
|----------|--------|
| C->T or A->G at editing window pos 4-8 | BE (preferred; higher efficiency) |
| Multi-base edit | PE |
| Transversion (C->G, C->A) | PE |
| Out-of-window single base | PE |
| Cancer-line variant scanning | PE or BE; cross-validate |
| iPSC / primary cell variant | PE (no bystander confounding) |
| Variant for drug-resistance | PE or BE; depends on chemistry |
| LoF / KO only | Cas9 (no template needed) |

## Thresholds

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| PRIDICT2 prediction threshold | >50% efficiency | Mathis 2024 |
| Intended edit % for screen | >5% | Anzalone 2019 |
| Scaffold incorporation | <2% (clean); <5% (acceptable) | Empirical |
| Indel byproduct (PE2) | <3% | Anzalone 2019 |
| PBS GC content | 40-55% | PRIDICT2 |
| PBS length | 11-13 nt | PRIDICT2 |
| RTT length | 10-20 nt | PRIDICT2 |
| Edit position from cut | 1-30 nt | Anzalone 2019 |

## Validation Checklist

- [ ] Pre-synthesis PRIDICT2 filter (>50%)
- [ ] Chromatin accessibility cross-referenced
- [ ] Pilot empirical efficiency match >70% of PRIDICT2 prediction
- [ ] CRISPResso2 PE mode for quantification
- [ ] Scaffold incorporation <5% library-wide
- [ ] Per-variant 3+ pegRNAs in library
- [ ] BE cross-validation (if BE-coverable variants)

## Related Skills

- crispr-screens/library-design - pegRNA library design
- crispr-screens/base-editing-analysis - Orthogonal BE for variant attribution
- crispr-screens/crispresso-editing - CRISPResso2 PE mode
- crispr-screens/hit-calling - Per-variant hit calling
- crispr-screens/screen-qc - Editing-efficiency QC
- variant-calling/variant-annotation - Annotate edited variants
- clinical-databases/clinvar-lookup - Pathogenicity annotation
