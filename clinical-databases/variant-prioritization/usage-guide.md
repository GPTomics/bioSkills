# Variant Prioritization - Usage Guide

## Overview

Run rare-disease variant prioritization pipelines on trio/quad WES/WGS: cascade filters (QC, population frequency, functional consequence, inheritance), call de novo variants (DeNovoGear, Triodenovo, DeNovoCNN) with IGV review, phase compound heterozygotes (WhatsHap), rank by phenotype concordance (Exomiser hiPHIVE, Phen2Gene, AMELIE), gate by ClinGen gene-disease validity, and flag ACMG SF v3.2 (Miller 2023; 81 genes). Variant *classification* (PVS1 decision tree, Pejaver 2022 PP3/BP4 thresholds, Tavtigian point system) is deferred to `clinical-databases/acmg-classification`.

## Prerequisites

```bash
pip install cyvcf2 pyhgvs pandas myvariant

# De novo callers
conda install -c bioconda denovogear
# Triodenovo: https://genome.sph.umich.edu/wiki/Triodenovo
# DeNovoCNN: https://github.com/Genomics-CISDC/DeNovoCNN

# Phenotype-driven
# Exomiser: https://www.sanger.ac.uk/tool/exomiser/
# Phen2Gene: pip install phen2gene
# AMELIE: https://amelie.stanford.edu/

# Compound het phasing
conda install -c bioconda whatshap

# HPO + ClinGen
# Download HPO: wget http://purl.obolibrary.org/obo/hp.obo
# ClinGen validity: https://search.clinicalgenome.org/kb/gene-validity
```

## Quick Start

Tell the agent what to do:
- "Run rare-disease prioritization on this trio WES VCF; output ranked candidates by inheritance pattern"
- "Identify de novo candidates with DeNovoGear; flag for IGV review"
- "For singleton WES, find compound heterozygous candidates with WhatsHap read-based phasing"
- "Rank candidates by phenotype concordance via Exomiser hiPHIVE; HPO terms: HP:0001250, HP:0001263, HP:0000175"
- "Cross-check candidates against ACMG SF v3.2 (81 genes, Miller 2023); flag CALM1/2/3 + other newly-added"

## Example Prompts

### Trio Pipeline

> "Run full trio pipeline on this Mendelian-suspected family. Filter to grpmax_faf95 < 0.0001 (dominant) or 0.005 (recessive). Apply DeNovoGear; IGV-review reportable DNVs. Output ranked candidates with inheritance pattern."

> "For trio with consanguineous parents, prioritize autosomal-recessive homozygous candidates; relax frequency to grpmax_faf95 < 0.01 for AR; cross-check ClinGen gene-disease validity (Moderate+)."

### Singleton Pipeline

> "Singleton WES, suspect AR; run WhatsHap read-based phasing for compound het candidates within 500 bp; flag for parental confirmation."

> "Singleton from consanguineous proband; prioritize hom-alt + compound het in genes with ClinGen Moderate+ validity."

### De Novo Analysis

> "DeNovoGear + IGV review on this trio. Posterior probability >= 0.9, parental coverage >= 10x, manual IGV review for all reportable DNVs."

> "Flag candidate DNVs in tandem repeat regions, segmental duplications, or with parental coverage < 10x as low-confidence."

### Phenotype-Driven

> "Submit Exomiser hiPHIVE with 8 specific HPO terms for this patient. Top-5 ranks should capture diagnosis 94% (Cipriani 2020)."

> "Run Phen2Gene as secondary ranking; compare top-5 to Exomiser. Manual review if discordant."

> "For phenotype with sparse HPO, run GADO (HPO-free option) and AMELIE (literature-mining); cross-check."

### ACMG SF Reporting

> "For each candidate variant, flag if gene is in ACMG SF v3.2 (Miller 2023; 81 genes). Note v3.2 additions: CALM1, CALM2, CALM3 (calmodulinopathies)."

> "Generate separate SF-only report for opt-in disclosure; reportable only if P/LP per `acmg-classification`."

### Mosaic Suspected

> "Patient has suspected mosaic PIK3CA / GNAS / NRAS disorder. Lower VAF threshold to 2-5%; require depth >= 200x; report affected tissue if available."

### VUS Reclassification

> "Re-evaluate VUS labeled prior to 2022 in this active diagnostic variant list; cross-check Genome Alert! (Yauy 2022) for ClinVar changes; flag for re-curation."

## What the Agent Will Do

1. Filter cascading: QC -> population frequency (grpmax_faf95) -> functional consequence -> inheritance pattern.
2. For trios, run DeNovoGear / Triodenovo / DeNovoCNN; require parental coverage >= 10x; IGV-review reportable DNVs.
3. For singletons, run WhatsHap read-based phasing for compound het candidates.
4. Rank candidates by phenotype concordance via Exomiser hiPHIVE (requires 5-10 specific HPO terms).
5. Apply ClinGen gene-disease validity gate (Moderate+ for diagnostic reporting; Limited/Disputed -> low confidence).
6. Cross-check ACMG SF v3.2 (81 genes); flag P/LP variants in SF genes for separate disclosure report.
7. For variant CLASSIFICATION (PVS1 / PP3 / BS1 / etc.), defer to `clinical-databases/acmg-classification`.
8. For VUS, recommend annual re-review with Genome Alert! monitoring monthly ClinVar releases.

## Tips

- Frequency filter: grpmax_faf95 < 0.0001 (dominant) or < 0.005 (recessive); use Whiffin gene-specific for known-disease genes.
- DeNovoGear / DeNovoCNN required for production DNV calling; raw Mendelian-violation analysis has 10-30% false-positive rate.
- All reportable DNVs need IGV review; tandem repeats, low parental coverage, parental mosaicism, mapping errors are common artifacts.
- WhatsHap read-based phasing works within ~500 bp; longer distances need parents or long-read.
- Exomiser hiPHIVE requires 5-10 specific HPO terms; sparse generic terms degrade accuracy significantly.
- ClinGen gene-disease validity (Limited / Moderate / Strong / Definitive) is mandatory upstream gate; many commercial panels include Limited-validity genes.
- ACMG SF v3.2 = 81 genes (Miller 2023); v3.2 added CALM1/2/3 (calmodulinopathies). Update from v3.1 (78 genes) needed.
- VUS reclassification cycle median ~5 years for active genes; up to 10 for orphan. Use Genome Alert! (Yauy 2022) to detect monthly changes.
- Mosaic variants: VAF 2-30%; require deep coverage >=200x and affected-tissue sampling when possible.
- Compound het cis vs trans: trio phasing is gold standard; without trio, WhatsHap or long-read.
- This skill produces CANDIDATE LISTS; ACMG CLASSIFICATION (PVS1 / PP3 / BS1 / etc.) is in `clinical-databases/acmg-classification`.
- For pharmacogenomic variants in candidate list, defer to `clinical-databases/pharmacogenomics`.

## Related Skills

- clinical-databases/acmg-classification - ACMG/AMP framework, PVS1 decision tree, Pejaver calibration
- clinical-databases/clinvar-lookup - Variant pathogenicity database
- clinical-databases/gnomad-frequencies - Population frequency for BS1/BA1
- clinical-databases/myvariant-queries - Aggregated annotation
- clinical-databases/pharmacogenomics - PGx-specific variants
- variant-calling/clinical-interpretation - Clinical reporting workflow
- variant-calling/filtering-best-practices - Upstream QC
