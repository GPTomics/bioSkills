# m6Anet Analysis - Usage Guide

## Overview

End-to-end m6A detection from Oxford Nanopore direct-RNA-sequencing (ONT DRS) signal data using m6Anet. Covers the required upstream pipeline (Dorado / Guppy basecalling -> minimap2 transcriptome alignment -> nanopolish eventalign -> m6anet dataprep -> m6anet inference), interpretation of the three per-site output columns (`probability_modified` posterior, `mod_ratio` stoichiometry, `n_reads` coverage) and the per-read modification probability output, the DRACH-only modeling constraint, minimum-coverage thresholds (20-50 reads per site), multi-condition comparison via xPore / Nanocompore / ELIGOS / Dorado native modification calling, reference-transcriptome version pinning, the cDNA-vs-DRS chemistry distinction, and orthogonal validation against MeRIP / GLORI.

## Prerequisites

```bash
pip install m6anet ont-fast5-api xpore nanocompore
conda install -c bioconda minimap2 samtools nanopolish
# Dorado: download binary from https://github.com/nanoporetech/dorado
```

Reference inputs:

- POD5 / FAST5 raw signal files from a direct-RNA-sequencing run (kit SQK-RNA002 or SQK-RNA004; NOT cDNA-Nanopore)
- Transcriptome FASTA (GENCODE / Ensembl transcripts.fa; version pinned per project)
- Dorado modification model name + version (e.g., `m6A_DRACH@v1`); pin per project

## Quick Start

- "Run the full m6Anet pipeline from POD5 to per-site m6A probabilities"
- "Filter m6Anet results for high-confidence sites (n_reads >= 20, prob >= 0.9)"
- "Compare m6A between two conditions with xPore Bayesian diffmod"
- "Run Nanocompore for modification-agnostic comparison between WT and KO"
- "First-pass m6A discovery with Dorado native modification calling on RNA004"
- "Cross-validate my m6Anet calls against published GLORI sites"

## Example Prompts

### Pipeline

> "Basecall RNA004 POD5 files with Dorado rna004_130bps_sup@v5.0.0; align to GENCODE v44 transcriptome with minimap2 -ax map-ont -uf -k14 --secondary=no; run nanopolish eventalign with the m6Anet-required flag set; run m6anet dataprep and inference."

> "Run m6anet dataprep and inference on a pre-computed eventalign.txt file from a HEK293T DRS run."

### Filtering

> "Filter m6anet results/data.site_proba.csv for high-confidence m6A sites: n_reads >= 20 AND probability_modified >= 0.9; export top 20 most-modified transcripts."

> "Compute per-read modification rate per site by aggregating per-read CSV with threshold 0.5; report alongside per-site probability."

### Differential Modification

> "Run xPore diffmod between WT and METTL3-KO with 3 replicates per condition; report diff_mod_rate >= 0.1 AND pval < 0.05."

> "Run Nanocompore sampcomp between two conditions on the same transcriptome; primary statistic GMM_logit_pvalue < 0.05."

### Dorado Native

> "Use Dorado modification calling on RNA004 data with m6A_DRACH@v1 model; report per-read modification probabilities for first-pass discovery."

> "Cross-check Dorado high-recall site set against m6Anet for high-precision filtering."

### Cross-Validation

> "Compare my m6Anet site calls against published GLORI HEK293T m6A sites at gene level."

> "Cross-reference m6Anet calls with MeRIP exomePeak2 peaks at common transcripts; report concordance."

### Sanity Check

> "Confirm my sequencing run is DRS (SQK-RNA002 / SQK-RNA004) and NOT cDNA-Nanopore; modification detection only works on DRS."

## What the Agent Will Do

1. Verify DRS protocol (kit SQK-RNA002 / SQK-RNA004) before any modification-detection pipeline
2. Basecall POD5 / FAST5 with Dorado (modern) or Guppy (legacy); pin basecaller + model version
3. Align to TRANSCRIPTOME (not genome) with minimap2 `-ax map-ont -uf -k14 --secondary=no`
4. Sort and index BAM with samtools
5. Run nanopolish eventalign with the m6Anet-required flag pair `--scale-events --signal-index` (plus `--summary` and `--threads` for housekeeping; `--samples` and `--print-read-names` are needed only by other downstream tools)
6. Run m6anet dataprep + inference
7. Apply minimum-coverage filter (`n_reads >= 20` conservative; `>= 50` stringent)
8. Apply probability filter (`probability_modified >= 0.9` high-confidence; `>= 0.7` discovery)
9. For stoichiometry estimates, compute per-read modification rate
10. For differential between conditions, run xPore (Bayesian, no matched IVT needed) or Nanocompore
11. For RNA004 first-pass screening, run Dorado native modification calling alongside m6Anet
12. Cross-validate top sites against published GLORI / m6A-Atlas / MeRIP peaks
13. Report per-site probability, per-read modification rate, coverage, and orthogonal validation status

## Tips

- m6Anet is DRACH-only by design. Non-DRACH sites are invisible. For non-DRACH discovery, use Dorado native or xPore / Nanocompore.
- cDNA-Nanopore data CANNOT be used for modification detection. PCR erases the signal. Only DRS (SQK-RNA002 / SQK-RNA004) works.
- RNA002 and RNA004 use different chemistry; models do NOT transfer. Pin chemistry version per project.
- Minimum coverage is non-negotiable. Sites with <20 reads at high probability are likely false positives.
- Per-site `probability_modified` is the model's posterior; per-read modification rate is the stoichiometry estimate.
- minimap2 for ONT DRS: use `-ax map-ont -uf -k14 --secondary=no` against TRANSCRIPTOME. Genome alignment breaks the m6Anet pipeline.
- nanopolish eventalign needs `--scale-events --signal-index` for m6Anet; `--samples --print-read-names` are needed only by other downstream tools (yanocomp, f5c-pipeline interop).
- m6Anet v1 used hyphenated CLI (`m6anet-dataprep`); v2 uses subcommands (`m6anet dataprep`). Pin version.
- Dorado modification models are versioned independently of the basecaller. Pin `m6A_DRACH@vX` explicitly.
- For absolute stoichiometry claims at named loci, cross-validate with GLORI (Liu 2023). Per-read modification rate has 5-15% error.
- Zhong 2023 *Nat Commun* 14:1906 benchmark documents tool tradeoffs; CHEUI is the only tool covering m6A AND m5C simultaneously.

## Related Skills

- merip-preprocessing - MeRIP fragment-level cross-validation
- m6a-peak-calling - MeRIP peak set for orthogonal comparison
- m6a-differential - MeRIP differential at the loci where m6Anet has direct stoichiometry
- modification-visualization - Metagene and browser-track rendering of m6Anet site calls
- long-read-sequencing/basecalling - Dorado / Guppy basecalling fundamentals
- long-read-sequencing/long-read-alignment - General minimap2 alignment patterns
- long-read-sequencing/long-read-qc - DRS yield, length distribution, basecall accuracy
- long-read-sequencing/nanopore-methylation - DNA methylation calling sibling
- variant-calling/vcf-basics - Per-site variant call framework analogue
- rna-quantification/featurecounts-counting - For site-vs-MeRIP-peak comparison
