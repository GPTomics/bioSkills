---
name: bio-clinical-databases-variant-prioritization
description: Prioritizes rare-disease variants from trio/quad WES/WGS with de novo (DeNovoGear, Triodenovo), compound-heterozygous phasing (WhatsHap), mosaic VAF tiering, phenotype-driven ranking (Exomiser, Phen2Gene, AMELIE), ClinGen gene-disease validity gating, and ACMG SF v3.2 secondary findings reporting. Use when running diagnostic exome / genome pipelines, identifying candidate Mendelian disease genes, screening for incidental findings, or auditing VUS reclassification cycles. The ACMG/AMP classification framework (PVS1 decision tree, Pejaver PP3/BP4 calibration, Tavtigian point system) is in clinical-databases/acmg-classification.
tool_type: python
primary_tool: pandas
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, cyvcf2 0.30+, pyhgvs 0.12+, Exomiser 14.0+ (Smedley 2015), Phen2Gene 1.2+ (Zhao 2020), DeNovoGear 1.1.1+ (Ramu 2013), WhatsHap 2.0+ (Patterson 2015), HPO 2024+ (Human Phenotype Ontology). ACMG Secondary Findings list is v3.2 (Miller 2023): 81 genes.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. Phenotype-driven prioritization REQUIRES high-quality HPO terms; without rich phenotypic input Exomiser/AMELIE degrade significantly.

# Rare-Disease Variant Prioritization Pipeline

**'Prioritize candidate disease-causing variants from this trio exome'** -> Filter to rare + functional + inheritance-consistent variants; rank by phenotype concordance; flag ACMG SF v3.2 incidental findings; report tiers with classification logic deferred to `acmg-classification`.

- Python (filtering pipeline): pandas + cyvcf2 + myvariant.info aggregation
- CLI (phenotype-driven ranking): `exomiser --analysis hiPHIVE-prioritised.yml`
- Python (de novo calling): DeNovoGear / Triodenovo / PossibleDeNovo
- CLI (compound het phasing): `whatshap phase --indels` for singletons; trio-based for families
- Python (HPO concordance): Phen2Gene / AMELIE / Phenolyzer
- VCEP curations: `https://cspec.genome.network/cspec/ui/svi/all`

## Pipeline Architecture: The Standard Rare-Disease Funnel

Typical trio exome enters as 40,000-100,000 variants per individual; reaches diagnostic candidate list of 1-10 variants through cascading filters:

| Stage | Filter | Variant count (typical trio) |
|-------|--------|------------------------------|
| Raw joint-called | -- | 100k-150k |
| QC filter (PASS, depth, GQ, missingness) | GATK best practices + Hail QC | 80k-120k |
| Population frequency | gnomAD grpmax_faf95 < 0.0001 (or disease-specific Whiffin max-credible-AF) | 5k-15k |
| Functional consequence | Coding / splice / regulatory | 1k-3k |
| Inheritance pattern | de novo / AR-hom / AR-compoundhet / X-linked / mosaic | 50-500 |
| Phenotype concordance | Exomiser hiPHIVE / Phen2Gene / AMELIE score | 5-50 |
| ACMG classification | Defer to `acmg-classification` | 1-10 |
| ACMG SF v3.2 cross-check | Miller 2023 (81 genes) | Separate output |

## Inheritance-Based Filtering

| Pattern | Filter |
|---------|--------|
| **De novo (DNV)** | Variant in proband, absent in both parents; needs trio | Apply DeNovoGear / Triodenovo / GATK PossibleDeNovo; visual IGV inspection (~10-30% false-positive rate without) |
| **Autosomal recessive; homozygous** | Hom-alt in proband; het in both parents | gnomAD grpmax_faf95 < 0.005 recessive threshold (Whiffin formula) |
| **Autosomal recessive; compound het** | Two het variants in same gene on opposite alleles | Trio-phased OR read-based phasing via WhatsHap (works within ~500 bp; longer needs parents or long-read) |
| **X-linked recessive** | Male proband hemizygous; carrier mother het | chrX coords; check Klinefelter / mosaic XXY |
| **X-linked dominant** | Het in affected; consider XCI skewing in females | Report XCI status if relevant |
| **Mitochondrial heteroplasmy** | mtDNA variant present at varying heteroplasmy across tissues | Use MITOMAP + HmtVar; ACMG criteria do not apply directly |
| **Mosaic** | Sub-clonal VAF in proband; absent in inherited transmissions | VAF 5-30% suggestive; tissue-dependent (blood vs buccal vs affected tissue) |

## De Novo Calling: Trio Analysis

**Goal:** Identify variants present in proband but absent in both parents with high specificity.

**Approach:** Use specialized DNV callers; supplement with manual IGV inspection.

| Tool | Approach | Use case |
|------|----------|----------|
| **DeNovoGear** (Ramu 2013 *Nat Methods*) | Bayesian, considers parent-of-origin | Standard for trio WES |
| **Triodenovo** (Wei 2015) | Bayesian + family-aware | Alternative |
| **GATK PossibleDeNovo annotation** | Hard filter | Quick prefilter; not standalone |
| **DeNovoCNN** (2024) | Deep learning trio caller | Most accurate as of 2024-2026 |

**False-DNV rate:** ~10-30% without manual IGV inspection; concentrated in:
- Tandem repeat regions (DNM rate inflated)
- Heterozygous parent with low coverage
- Mosaic parents (parental mosaicism transmitted to >1 offspring)
- Mapping errors in segmental duplications

## Phenotype-Driven Prioritization

| Tool | Approach | Performance (typical benchmark) | Fails when |
|------|----------|--------------------------------|-----------|
| **Exomiser** (Smedley 2015 *Nat Protoc*) | hiPHIVE: phenotype + interactome + sequence damage | 74% top-1; 94% top-5 (Cipriani 2020) | Sparse HPO (< 5 specific terms); novel-disease gene |
| **Phen2Gene** (Zhao 2020 *NARGAB*) | HPO-to-gene mapping; faster than Exomiser | Similar top-5 | Phenotype-only filtering insufficient |
| **AMELIE** (Birgmeier 2020 *Sci Transl Med*) | Literature-mining + phenotype | Best when literature is rich | New / rare disease without literature; specific patient HPO unmatched |
| **Phenolyzer** (Yang 2015 *Nat Methods*) | Phenotype-based gene scoring | Legacy | Modern multi-feature tools (Exomiser, AMELIE) preferred |
| **GADO** (Deelen 2019 *Nat Commun*) | Gene Network-based; HPO-free option | When HPO is sparse | Phenotype-rich cases where Exomiser hiPHIVE wins |
| **CADA** (Peng 2021) | Cross-species gene prioritization | Animal model integration | Genes without orthologs; rare-disease without animal model |

**Critical requirement:** all phenotype-driven tools degrade significantly with sparse HPO terms. Capture 5-10 specific HPO terms; avoid generic "intellectual disability" alone.

## ClinGen Gene-Disease Validity: Mandatory Gating

Strong et al. 2017 *AJHG* + ClinGen ongoing curation: **Limited / Moderate / Strong / Definitive** evidence per gene-disease pair.

| Category | When to apply |
|----------|---------------|
| **Definitive** | Strong literature evidence + functional / population genetic evidence | Apply full ACMG framework |
| **Strong** | -- | Apply full framework |
| **Moderate** | -- | Apply framework but flag |
| **Limited** | Single case report or weak segregation | Treat candidate cautiously; PP2 / BP1 should not apply |
| **Disputed** | Contradicting evidence | Do not call pathogenic without VCEP curation |
| **No Known Disease Relationship** | Gene not associated with the queried disease | Do not call |

**Many commercial panels include genes with only Limited validity.** ClinGen-curated `https://search.clinicalgenome.org/kb/gene-validity` is the authoritative directory.

## ACMG Secondary Findings v3.2 (Miller 2023 *Genet Med* 25:100866)

81 genes for opt-in/opt-out reporting on clinical exome/genome. Growth: 56 -> 59 -> 73 -> 78 -> 81. **v3.2 additions: CALM1, CALM2, CALM3** (calmodulinopathy; long QT / CPVT; high actionability via beta-blockade + ICD).

Inclusion criteria: ClinGen Strong or Definitive gene-disease validity + ClinGen ADWG actionability scoring.

```python
ACMG_SF_V3_2_GENES = [
    # Cardiomyopathies
    'ACTA2', 'ACTC1', 'COL3A1', 'DES', 'FBN1', 'FLNC', 'GLA', 'LMNA', 'MYBPC3',
    'MYH11', 'MYH7', 'MYL2', 'MYL3', 'PRKAG2', 'PKP2', 'RBM20', 'SCN5A', 'SMAD3',
    'TGFBR1', 'TGFBR2', 'TMEM43', 'TNNI3', 'TNNT2', 'TPM1', 'TTN',
    # CALM v3.2 additions (calmodulinopathies)
    'CALM1', 'CALM2', 'CALM3',
    # Arrhythmias and channelopathies
    'CACNA1S', 'KCNH2', 'KCNQ1', 'RYR1', 'RYR2',
    # Vascular
    'ACVRL1', 'ENG',
    # Cancer predisposition
    'APC', 'ATM', 'BAP1', 'BMPR1A', 'BRCA1', 'BRCA2', 'BRIP1', 'CDH1', 'CDKN2A',
    'CHEK2', 'GREM1', 'HOXB13', 'MAX', 'MEN1', 'MLH1', 'MSH2', 'MSH6', 'MUTYH',
    'NF2', 'PALB2', 'PMS2', 'PTEN', 'RAD51C', 'RAD51D', 'RB1', 'RET', 'SDHAF2',
    'SDHB', 'SDHC', 'SDHD', 'SMAD4', 'STK11', 'TMEM127', 'TP53', 'TSC1', 'TSC2',
    'VHL', 'WT1',
    # Other
    'FH', 'GAA', 'HFE', 'HNF1A', 'LDLR', 'NTRK1', 'OTC', 'PCSK9', 'TTR'
]
# Note: above list is illustrative; pin to Miller 2023 supplement for exact set.
```

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Trio WES, suspected Mendelian | Full pipeline with DeNovoGear + Exomiser + HPO | Standard rare-disease workflow |
| Singleton WES | WhatsHap read-based phasing + AR-hom + AR-compoundhet candidates | Compound het hard without trio |
| Suspected mosaic | Lower VAF threshold (2-30%); deep coverage (>200x) | Standard tools miss mosaic |
| Long-read genome | Add SV calling + STR repeat expansion | SVs miss in short-read |
| Newborn screening (BabyScreen+) | 605-gene Mendelian panel with current ACMG SF v3.2 | Lynch / Brett 2025 *Nat Med* |
| Cancer predisposition | ClinGen Hereditary Cancer VCEPs + ACMG SF cancer subset | Use VCEP CSpec |
| Cardiomyopathy / arrhythmia | ClinGen HCM / DCM / LQT VCEPs | Strict gene-disease validity |
| Population screening | ACMG SF v3.2 (81 genes) opt-in/opt-out | Miller 2023 |

## Standard Pipeline Workflow

**Goal:** From a trio joint-called VCF, output ranked candidate variants with inheritance pattern, phenotype concordance, and ACMG SF flags.

**Approach:** Cascading filters with QC, population frequency, functional consequence, inheritance, phenotype.

```python
from cyvcf2 import VCF
import pandas as pd
from pathlib import Path

# Quality + population frequency filter (apply first)
def filter_qc_and_frequency(vcf_path, max_grpmax_faf95=0.0001, min_dp=10, min_gq=20):
    '''Stage 1: QC + frequency filter. Reduces 100k -> ~5-15k variants.'''
    vcf = VCF(vcf_path)
    samples = vcf.samples  # e.g., [proband, mother, father]
    rows = []
    for v in vcf:
        if v.FILTER is not None:
            continue
        if min(v.gt_depths) < min_dp:
            continue
        if v.QUAL is not None and v.QUAL < min_gq:
            continue
        gnomad = (v.INFO.get('grpmax_faf95') or v.INFO.get('AF_grpmax') or
                  v.INFO.get('AF_popmax') or 0)
        if gnomad > max_grpmax_faf95:
            continue
        rows.append({
            'chrom': v.CHROM, 'pos': v.POS, 'ref': v.REF, 'alt': v.ALT[0],
            'genotypes': dict(zip(samples, v.gt_types.tolist())),
            'depth': dict(zip(samples, v.gt_depths.tolist())),
            'gnomad_faf95': gnomad,
            'consequence': v.INFO.get('CSQ', '').split('|')[1] if v.INFO.get('CSQ') else None
        })
    return pd.DataFrame(rows)


def call_de_novo(df, proband, mother, father):
    '''Stage 2: identify DNV candidates: hom-ref both parents, het/hom-alt proband.

    Implements Mendelian-violation logic; supplement with DeNovoGear or DeNovoCNN
    for production (this implementation has 10-30% false-positive rate without IGV).
    '''
    is_dnv = []
    for _, row in df.iterrows():
        gts = row['genotypes']
        if gts[mother] == 0 and gts[father] == 0 and gts[proband] in (1, 3):
            # Mother hom-ref AND father hom-ref AND proband het OR hom-alt
            # Confidence boost: depth at parent sites should be >= 10 to trust hom-ref
            if row['depth'][mother] >= 10 and row['depth'][father] >= 10:
                is_dnv.append(True)
                continue
        is_dnv.append(False)
    df['is_de_novo_candidate'] = is_dnv
    return df


def call_compound_het(df, proband, mother, father, gene_col='gene'):
    '''Stage 3: identify compound het: two het variants in same gene, one from each parent.

    Trio phasing is gold standard; singletons require WhatsHap read-based phasing.
    '''
    het_in_proband = df[df['genotypes'].apply(lambda gts: gts[proband] == 1)]
    candidate_genes = []
    for gene in het_in_proband[gene_col].unique():
        if pd.isna(gene):
            continue
        gene_variants = het_in_proband[het_in_proband[gene_col] == gene]
        # Need >= 2 variants; one inherited from each parent
        maternal_het = gene_variants[gene_variants['genotypes'].apply(
            lambda gts: gts[mother] == 1 and gts[father] == 0)]
        paternal_het = gene_variants[gene_variants['genotypes'].apply(
            lambda gts: gts[father] == 1 and gts[mother] == 0)]
        if len(maternal_het) >= 1 and len(paternal_het) >= 1:
            candidate_genes.append(gene)
    df['is_compound_het_candidate'] = df[gene_col].isin(candidate_genes)
    return df


def flag_acmg_sf(df, acmg_sf_genes, gene_col='gene', clnsig_col='clinvar_sig'):
    '''Stage: flag ACMG Secondary Findings (Miller 2023 v3.2; 81 genes).

    Only P/LP variants in SF genes are reportable as secondary findings.
    '''
    df['is_acmg_sf_candidate'] = (
        df[gene_col].isin(acmg_sf_genes) &
        df[clnsig_col].astype(str).str.contains('athogenic', na=False)
    )
    return df


def filter_by_clingen_validity(df, validity_table, gene_col='gene',
                                min_validity='Moderate'):
    '''Gate on ClinGen gene-disease validity. Limited or Disputed -> low confidence.

    validity_table: DataFrame from `https://search.clinicalgenome.org/kb/gene-validity`
    '''
    rank = {'No Known Disease Relationship': 0, 'Disputed': 0, 'Limited': 1,
            'Moderate': 2, 'Strong': 3, 'Definitive': 4}
    min_rank = rank[min_validity]
    df_merged = df.merge(validity_table, on=gene_col, how='left')
    df_merged['validity_rank'] = df_merged['gene_validity'].map(rank).fillna(0)
    df_merged['pass_validity'] = df_merged['validity_rank'] >= min_rank
    return df_merged


def phenotype_score_with_exomiser_yml(yml_path, vcf_path, hpo_terms, output_dir):
    '''Emit Exomiser command for phenotype-driven ranking.

    HPO terms (e.g., HP:0001250 for seizures) must be SPECIFIC.
    Sparse generic HPO degrades Exomiser hiPHIVE accuracy significantly.
    '''
    return (f'java -jar exomiser-cli-14.0.0.jar --analysis {yml_path} '
            f'--vcf {vcf_path} --hpo {",".join(hpo_terms)} '
            f'--output-dir {output_dir}')
```

## Per-Operation Failure Modes

**1. De novo with false-positive rate 10-30%**
- Trigger: Report DNV candidates from Mendelian-violation analysis without IGV inspection.
- Mechanism: Tandem-repeat regions, low-coverage parents, parental mosaicism, mapping errors in segmental duplications all produce false DNVs.
- Symptom: 10-30% of reported DNVs are artifacts.
- Fix: Use DeNovoGear / DeNovoCNN (Bayesian frameworks); manually inspect candidates in IGV; check parental coverage at site.

**2. Compound het without phasing**
- Trigger: Report two hets in same gene as compound het without confirming phase.
- Mechanism: Trans (compound het) vs cis (same chromosome) is critical for AR mechanism.
- Symptom: False-positive compound het when both variants are in cis.
- Fix: Trio phasing if available; WhatsHap read-based phasing for variants within ~500 bp; consider long-read for broader phasing.

**3. Limited-validity gene reported as diagnostic**
- Trigger: Gene appears on commercial panel; variant labeled disease-causing.
- Mechanism: Commercial panels often include Limited or Disputed validity genes.
- Symptom: False-positive diagnostic report.
- Fix: Cross-check ClinGen gene-disease validity; reject Limited / Disputed without VCEP curation.

**4. Sparse HPO terms degrading Exomiser**
- Trigger: Submit Exomiser with single generic HPO (e.g., HP:0001250 "Seizure" only).
- Mechanism: Phenotype-driven prioritization relies on HPO-to-gene network; sparse terms reduce discriminative power.
- Symptom: Top-5 rank includes implausible genes; correct diagnosis sub-rank.
- Fix: Capture 5-10 specific HPO terms (e.g., "infantile spasms with hypsarrhythmia", "facial dysmorphism with hypertelorism").

**5. ACMG SF v3.1 used instead of v3.2**
- Trigger: Pipeline reports SF based on 78-gene v3.1 list; misses CALM1/2/3 calmodulinopathies.
- Mechanism: v3.2 (Miller 2023) added CALM1, CALM2, CALM3.
- Symptom: Misses calmodulinopathy SF; high-actionability long-QT/CPVT not flagged.
- Fix: Use Miller 2023 v3.2 list (81 genes); re-run prior cohorts.

**6. Mosaic variants below standard VAF threshold**
- Trigger: Filter at VAF >= 30% on standard pipeline.
- Mechanism: Mosaic variants frequently 2-30% VAF; below threshold filters them out.
- Symptom: Mosaic disease missed (e.g., Proteus syndrome PIK3CA, McCune-Albright GNAS).
- Fix: For suspected mosaic disorders, deep coverage (>= 200x); VAF threshold 2-5%; sample affected tissue when possible.

**7. ClinVar P variant in Limited-validity gene**
- Trigger: Variant labeled P in ClinVar; gene-disease validity is Limited.
- Mechanism: ClinVar P is variant-level assertion; gene-disease validity is the upstream question.
- Symptom: Reported P variant in non-disease-associated gene.
- Fix: Apply ClinGen gene-disease validity gate BEFORE variant-level interpretation.

**8. VUS reclassification gaps**
- Trigger: VUS labeled 2017 still in active diagnostic report 2025.
- Mechanism: Median VUS reclassification cycle ~5 years for actively-curated genes (Harrison 2017 follow-up).
- Symptom: Stale classifications drive incorrect clinical decisions.
- Fix: Annual VUS re-review for active diagnostic variants; tools like Genome Alert! (Yauy 2022) automate detection of monthly ClinVar changes.

**9. Inheritance pattern assumed wrong**
- Trigger: Assume AD inheritance for a gene with variable expressivity / incomplete penetrance.
- Mechanism: AD genes can have AR variants in functionally significant compound het pattern.
- Symptom: Miss AR mechanism in mostly-AD gene.
- Fix: Allow multi-inheritance candidate generation; cross-check ClinGen gene-disease inheritance.

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| Exomiser ranks low; ClinVar says P | Sparse or wrong HPO terms; rare disease in atypical gene | Re-run with full HPO; manual review |
| ClinVar P + ClinGen Limited validity | Variant-level vs gene-disease tension | Treat as candidate; require VCEP curation or functional evidence |
| DeNovoGear high posterior; trio coverage uneven | Parental mosaicism or mapping error | IGV review; consider parent-of-origin testing |
| Compound het in phasing-ambiguous gene | Distance > 500 bp; can't phase from reads | Trio phasing; long-read confirmation |
| SF gene with V3.1 list; missing CALM | Miller 2023 v3.2 update | Re-run with v3.2 (81 genes) |
| Phenotype tool disagrees with clinical | Tool-specific phenotype model; literature gap | Cross-check with AMELIE for literature-mining alternative |
| Mosaic suspected but standard pipeline negative | VAF below 30% threshold | Deep targeted sequencing or affected tissue |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Rare-disease frequency filter | grpmax_faf95 < 0.0001 | ClinGen SVI |
| Recessive disease filter | grpmax_faf95 < 0.005 | ClinGen SVI |
| Whiffin gene-specific max-credible-AF | Computed per gene + disease | Whiffin 2017 |
| DNV minimum parental coverage | >= 10x both parents | Standard |
| DNV manual IGV review | Required for all reportable DNVs | Standard |
| Compound het phasing | <= 500 bp read-based; trio gold standard | WhatsHap |
| Exomiser top-1 diagnostic rank | 74%; top-5 94% (with rich HPO) | Cipriani 2020 |
| ACMG SF v3.2 genes | 81 (Miller 2023) | Miller 2023 *Genet Med* |
| VUS reclassification cycle | Median 5 years for active genes; up to 10 for orphan | Harrison 2017 follow-up |
| Mosaic VAF threshold | 2-30% | Convention |
| ClinGen gene-disease validity gate | Moderate or Strong minimum for diagnostic reporting | ClinGen SVI |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Too many candidate variants (>50) | Frequency filter too loose | Tighten to grpmax_faf95 < 0.0001 (dominant) or 0.005 (recessive) |
| No DNV candidates in obvious DNV phenotype | False-negative DNV calling | DeNovoGear / DeNovoCNN; check parental sample swap |
| Compound het in gene known AD only | Phasing not validated | Confirm phase via trio or long-read |
| Exomiser top hit unrelated to phenotype | HPO too generic or wrong | Add specific HPO; check ontology version |
| Mosaic disease missed | VAF threshold too high | Deep coverage; affected tissue sampling; VAF 2-5% |
| SF gene match flagged but variant benign | Wrong variant classification | Apply ACMG framework via `acmg-classification` skill |
| Genotype-phenotype discordance | Locus heterogeneity OR multi-gene contribution | Run digenic / oligogenic analysis tools |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why grpmax_faf95 instead of AF?" | grpmax_faf95 is the Whiffin 2017 ClinGen-recommended frequency; excludes bottleneck groups; per ACMG SVI specifications. |
| "Compound het without phase confirmation" | Trio phased; if singleton, WhatsHap read-based for variants within 500 bp; long-read otherwise. |
| "DNV call without IGV review?" | All reportable DNVs underwent IGV inspection; we report posterior probability + parental coverage. |
| "ClinGen Limited validity gene" | Excluded per gate; we require Moderate or higher for reportable diagnostic candidates. |
| "Why ACMG SF v3.2 not v3.1?" | v3.2 (Miller 2023) added CALM1/2/3 calmodulinopathies (high actionability). We use current. |
| "Phenotype-driven prioritization with single HPO term?" | We submit 5-10 specific HPO terms; sparse input degrades Exomiser. |
| "ACMG classification logic?" | Variant prioritization (this skill) outputs candidates; ACMG classification (PVS1 / PP3 / BS1 / etc.) is in `acmg-classification` skill. |
| "Why not VarSome / Franklin automated ACMG?" | We report aggregated annotations via myvariant.info; ACMG classification per `acmg-classification` skill using Tavtigian point system + Pejaver 2022 calibration. |

## References

- Richards S et al. 2015. Standards and guidelines for the interpretation of sequence variants. *Genet Med* 17:405. (ACMG/AMP)
- Miller DT et al. 2023. ACMG SF v3.2 list for reporting of secondary findings in clinical exome and genome sequencing. *Genet Med* 25:100866.
- Smedley D et al. 2015. Next-generation diagnostics and disease-gene discovery with the Exomiser. *Nat Protoc* 10:2004.
- Zhao M et al. 2020. Phen2Gene: rapid phenotype-driven gene prioritization for rare diseases. *NARGAB* 2:lqaa032.
- Birgmeier J et al. 2020. AMELIE accelerates Mendelian patient diagnosis directly from the primary literature. *Sci Transl Med* 12:eaau9113.
- Cipriani V et al. 2020. An improved phenotype-driven tool for rare Mendelian variant prioritization. *Genes* 11:460.
- Ramu A et al. 2013. DeNovoGear: de novo indel and point mutation discovery and phasing. *Nat Methods* 10:985.
- Patterson M et al. 2015. WhatsHap: weighted haplotype assembly for future-generation sequencing reads. *J Comput Biol* 22:498.
- Strong A et al. 2017. Gene-disease validity framework. *AJHG* 100:895.
- Whiffin N et al. 2017. Using high-resolution variant frequencies to empower clinical genome interpretation. *Genet Med* 19:1151.
- Lynch F, Brett T et al. 2025. BabyScreen+ implementation results from genomic newborn screening. *Nat Med*.
- Yauy K et al. 2022. Genome Alert! *Genet Med* 24:S1098. (VUS reclassification monitoring)
- ClinGen gene-disease validity: `https://search.clinicalgenome.org/kb/gene-validity`
- HPO: `https://hpo.jax.org/`
- ACMG SF v3.2 supplement: `https://www.gimjournal.org/article/S1098-3600(23)00879-1/fulltext`

## Related Skills

- clinical-databases/acmg-classification - PVS1 / PP3 / BS1 / PM2 calibration and Tavtigian point system
- clinical-databases/clinvar-lookup - Variant pathogenicity database query
- clinical-databases/gnomad-frequencies - Population frequency filtering
- clinical-databases/myvariant-queries - Aggregated annotation
- clinical-databases/pharmacogenomics - PGx variant handling
- variant-calling/clinical-interpretation - Clinical reporting workflow
- variant-calling/filtering-best-practices - Upstream QC
