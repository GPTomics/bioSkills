# Splice Variant Prediction - Usage Guide

## Overview
Predict whether a DNA variant alters mRNA splicing using sequence-based deep-learning tools - SpliceAI (10kb context CNN), Pangolin (multi-tissue), MMSplice (modular per-region), SpliceTransformer/TrASPr (tissue-aware transformers), SpliceVault (empirical 300K-RNA lookup of likely mis-splicing outcomes), and CADD-Splice. Applies the ClinGen SVI 2023 framework for ACMG/AMP variant interpretation (PVS1, PP3, BP4 evidence codes), HGVS splicing nomenclature, extended-window scoring for deep-intronic pseudoexons, and tissue-specific predictions.

## Prerequisites
```bash
# Python
pip install spliceai tensorflow pyfaidx pyensembl pysam

# Pangolin (splicing) from GitHub - NOT the PyPI 'pangolin' (SARS-CoV-2 lineage tool)
pip install git+https://github.com/tkzeng/Pangolin.git
# Plus a gffutils annotation database
python -c "import gffutils; gffutils.create_db('gencode.v45.gff3', 'gencode.db')"

# MMSplice (specific dependencies)
pip install mmsplice

# SpliceVault: web at kidsneuro.shinyapps.io/splicevault or R/Python at github.com/kidsneuro-lab/SpliceVault

# Reference files
# - Genome FASTA (GRCh38 or GRCh37)
# - GENCODE GTF / GFF3
```

## Quick Start
Tell your AI agent what you want to do:
- "Predict SpliceAI delta scores for variants in my VCF"
- "Score a panel of clinical variants for splice impact with concordance across SpliceAI, Pangolin, and MMSplice"
- "Identify deep-intronic pseudoexon-creating variants using extended-window SpliceAI"
- "Classify a candidate variant under ClinGen SVI 2023 splicing rules (PVS1/PP3/BP4)"
- "Predict tissue-specific splicing impact with Pangolin for a brain-disease variant"

## Example Prompts

### Single Variant
> "For variant chr17:g.7676154A>G (TP53), compute SpliceAI delta score with -D 50 and Pangolin tissue-specific predictions for brain."

### VCF Annotation
> "Run SpliceAI on my VCF of clinical variants with grch38 reference, then filter to variants with delta_max >= 0.2 and apply ClinGen PP3 thresholds."

### Deep-Intronic
> "I have an unsolved Mendelian case; re-run SpliceAI with -D 2000 to check for pseudoexon-creating deep-intronic variants."

### Concordance
> "Run SpliceAI + Pangolin + MMSplice on a candidate VUS panel; flag discordant predictions for RNA validation."

### ASO Design
> "Predict the splicing impact of occluding a specific ESE region with a 22-mer ASO using SpliceAI on the masked sequence."

### Empirical Outcome
> "For a canonical 5'ss-disrupting variant, query SpliceVault to predict whether it causes exon skipping vs cryptic site activation."

## What the Agent Will Do
1. Validate variant nomenclature (HGVS) with VariantValidator or Mutalyzer
2. Run SpliceAI with appropriate distance window (50 default; 500-2000 for deep-intronic)
3. Run Pangolin if tissue-specific predictions are needed
4. Run MMSplice for calibrated ΔPSI in cassette-exon contexts
5. Query SpliceVault for empirical mis-splicing outcomes
6. Apply ClinGen SVI 2023 thresholds for PP3/BP4/PVS1 evidence
7. Recommend RNA validation for clinical-grade reporting

## Tips
- ClinGen SVI 2023: SpliceAI >=0.2 = PP3, <=0.1 = BP4, both applied at SUPPORTING weight only (0.5/0.8 raise SpliceAI precision but are not ClinGen evidence-strength upgrades; moderate/strong needs RNA/functional PS3 evidence)
- SpliceAI default window is 50nt; extend to 500-2000 for unsolved Mendelian cases (5-15% are deep-intronic pseudoexons)
- Tissue-agnostic prediction (SpliceAI) is fine for screening; use Pangolin/SpliceTransformer when disease tissue is known
- RNA validation (PS3) supersedes computational evidence (PP3); run RT-PCR or RNA-seq when possible
- Branchpoint variants are weakly predicted by all current tools; use BPHunter (Zhang 2022) as supplement
- Concordance across SpliceAI + Pangolin + MMSplice strengthens evidence; discordance flags need experimental work
- Log SpliceAI version and distance window with every clinical call
- SpliceAI alone is NOT sufficient for PVS1; canonical site disruption requires gene-level LoF context

## Related Skills

- splicing-qc - MaxEntScan + library QC for confirming predicted impact
- splicing-quantification - Empirical PSI from RNA-seq to validate predictions
- outlier-splicing-detection - FRASER2/DROP for RNA-seq confirmation in clinical samples
- variant-calling/clinical-interpretation - Broader ACMG/AMP variant interpretation
- variant-calling/variant-annotation - VEP plugin integration for SpliceAI
