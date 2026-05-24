---
name: bio-clinical-databases-pharmacogenomics
description: Queries PharmGKB / CPIC / DPWG for drug-gene interactions; calls CYP2D6/CYP2C9/CYP2C19/DPYD/TPMT/NUDT15/UGT1A1/SLCO1B1 star alleles and phenotype with PharmCAT, Cyrius (CYP2D6 structural variants), Aldy, Stargazer; applies Caudle 2020 activity-score translation. Use when implementing pharmacogenomic-guided prescribing, applying CPIC vs DPWG guidance, screening HLA risk alleles for ICI / antiepileptics / abacavir, or interpreting compound TPMT+NUDT15 thiopurine risk.
tool_type: mixed
primary_tool: PharmCAT
---

## Version Compatibility

Reference examples tested with: PharmCAT 2.13+, Cyrius 1.1+ (Chen 2021), Aldy 4.0+, Stargazer 2.0+, StarPhase 1.0+ (PacBio HiFi), HIBAG 1.40+, requests 2.31+, pandas 2.2+. CPIC guideline versions are gene-specific; PharmVar releases are quarterly. The 2024 DPYD update (Lam et al. *Clin Pharmacol Ther*) replaced single-variant logic with the activity-score system; the 2025 TPMT/NUDT15 update (Maillard 2026) refines compound-IM dosing.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. PharmVar is the authoritative star-allele source (`https://www.pharmvar.org`); the older Human CYP Allele Nomenclature Database was deprecated in 2017.

# Pharmacogenomics; Star Alleles, Activity Scores, and CPIC/DPWG Guidance

**'What is my patient's CYP2D6 metabolizer status and should I adjust their tamoxifen dose?'** -> Call star alleles (haplotype-level), translate diplotype -> activity score -> phenotype, apply CPIC + DPWG dosing.

- CLI (recommended): `pharmcat -vcf input.vcf.gz -o pharmcat_out`; CPIC-recommended, single-tool reporting
- CLI (CYP2D6 SV-aware): `cyrius -m sample.bam -o cyrius_out`; mandatory addition for CYP2D6
- CLI (multi-gene CN-aware): `aldy genotype -p illumina sample.bam`; alternative
- CLI (long-read 8-field): PacBio HiFi `starphase`; transplant-grade including HLA
- R (SNP-array): HIBAG for HLA-B*57:01/B*15:02/B*58:01/A*31:01 imputation
- API: `requests.get('https://api.pharmgkb.org/v1/data/clinicalAnnotation', ...)`

## Governance: CPIC vs DPWG vs PharmGKB vs FDA

These four authorities are routinely conflated. They differ in scope, scale, and recommendations:

| Authority | Scope | Output | Anchors |
|-----------|-------|--------|---------|
| **CPIC** (US Clinical Pharmacogenetics Implementation Consortium) | Once you have a result, what to prescribe | Level A/B/C/D gene-drug pair + strength of recommendation per phenotype + evidence quality | ~26 guidelines, ~25 genes, 100+ drugs as of 2026 |
| **DPWG** (Dutch Pharmacogenetics Working Group) | Whether to test AND what to prescribe | 5-pt (0-4) evidence + 7-pt (AA-F) clinical-relevance scale | G-Standaard (Dutch EHR-integrated); RCT-validated via PREPARE |
| **PharmGKB clinical annotation levels** | Evidence cataloguing | 1A/1B/2A/2B/3/4 | 1A = guideline OR medical-society OR PGRN/eMERGE implementation; NOT pure evidence |
| **FDA Table of Pharmacogenomic Biomarkers** | Drug label info | ~300 drugs (informational) | NOT an actionability list; many entries are dosing-suggestion-only |
| **FDA Table of Pharmacogenetic Associations** | Actionable subset | Closer to CPIC | Compare head-to-head with CPIC |

**Bank et al 2019** *Clin Pharmacol Ther* 105:951 (DOI 10.1002/cpt.762; first online 2018, print 2019) is the canonical CPIC-vs-DPWG comparison. Notable disagreements:
- CYP2D6 IM + multiple antidepressants: DPWG actionable; CPIC says insufficient evidence.
- HLA-B*15:11 carbamazepine: DPWG actionable; CPIC silent.
- CYP2C19 IM + voriconazole: dosing magnitudes differ 25-50%.

**Donnelly 2024** critiques: (1) EUR over-representation in discovery cohorts; (2) most PGx RCTs are open-label / prescriber-unblinded; (3) publication bias in antiseizure PGx may overstate effects ~2x; (4) subjective composite endpoints.

## PharmGKB Clinical Annotation Levels: What 1A Actually Means

| Level | Requirement |
|-------|-------------|
| **1A** | Variant-drug pair appears in CPIC guideline OR medical-society guideline OR is implemented at a PGRN/eMERGE site |
| **1B** | Replication in multiple cohorts; preponderance of evidence; no formal guideline yet |
| **2A** | Replicated association in a VIP (Very Important Pharmacogene) |
| **2B** | Replicated association in non-VIP gene |
| **3** | Single significant association OR mixed-evidence variant-drug pair |
| **4** | In vitro / case report / molecular evidence only |

1A does NOT require RCT evidence; mechanism + guideline status suffices.

## Star Allele Nomenclature (PharmVar)

PharmVar (`https://www.pharmvar.org`) is authoritative for: CYP1A1, CYP1A2, CYP1B1, CYP2A6, CYP2A13, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP2F1, CYP2J2, CYP2R1, CYP2S1, CYP2W1, CYP3A4, CYP3A5, CYP3A7, CYP3A43, CYP4A11, CYP4F2, CYP19A1, CYP26A1, DPYD, NUDT15, SLCO1B1, TPMT.

**A star allele is a haplotype, not a single variant.** *Suballeles* (*1.001, *1.002, etc.) encode the exact SNV+indel pattern within a defined functional haplotype.

**The *1 reference is the PharmVar consensus reference, NOT biological wild type.** Defined as the absence of all known functional variants at the locus.

### CYP2D6 Activity Scores (Caudle 2020 *Clin Transl Sci*; DOI 10.1111/cts.12692)

| Phenotype | Activity score (AS) range |
|-----------|--------------------------|
| **PM (Poor Metabolizer)** | 0 |
| **IM (Intermediate Metabolizer)** | 0 < AS < 1.25 |
| **NM (Normal Metabolizer)** | 1.25 <= AS <= 2.25 |
| **UM (Ultra-rapid)** | AS > 2.25 |

Key per-allele activity values (selected):

| Allele | Activity | Notes |
|--------|----------|-------|
| \*1, \*2, \*35 | 1.0 | Normal |
| \*3, \*4, \*5 (gene deletion), \*6, \*7, \*8, \*11, \*12, \*15, \*19, \*20, \*36, \*40, \*42 | 0 | No function |
| \*9, \*41, \*17, \*29 | 0.5 | Decreased function (substrate-specific caveats for \*17) |
| **\*10** | **0.25** | **Caudle 2020 RESET from 0.5 to 0.25**; reclassified large fractions of East-Asian populations to IM |
| \*68 | 0 | Hybrid; non-functional |

**\*4xN is clinically silent:** a no-function allele multiplied by N is still no-function. Reporting *4xN as UM is the most-common reportable error in clinical PGx.

### CYP2D6 Structural Complexity

CYP2D6 on 22q13.2 sits adjacent to the highly-similar CYP2D7 pseudogene. Four classes of structural variant that no SNV-only caller can resolve:

1. **Gene deletion (\*5):** ~13 kb deletion; activity 0; diagnostic *REP6/REP7* breakpoint.
2. **Gene duplication/multiplication (\*1xN, \*2xN, \*4xN, \*10xN, \*17xN, \*35xN, \*36xN):** Tandem copies; clinical impact depends on which allele is amplified; **\*4xN is clinically silent**.
3. **CYP2D7 -> CYP2D6 hybrids (\*13):** Pseudogene fused 5'; non-functional.
4. **CYP2D6 -> CYP2D7 hybrids (\*36, \*61, \*63, \*68, \*83):** 5' CYP2D6 with 3' pseudogene exon 9 conversion; typically embedded in duplications upstream of *10 (East Asian) or upstream of *4 (European).

**GATK / DeepVariant alone cannot call any of these.** They operate on multi-mapper-filtered BAMs; 97%+ identity between CYP2D6 and CYP2D7 produces silent miscalls of every \*5, \*13, \*36, \*68, \*4xN sample.

## Algorithmic Taxonomy: Star Allele Callers

| Tool | CYP2D6 SV | CYP2D6 CN | Other PGx genes | Phased | Validation | Fails when |
|------|-----------|-----------|-----------------|--------|------------|-----------|
| **PharmCAT** (Sangkuhl 2020 *Clin Pharmacol Ther*) | No (consumes outside SV calls) | No | 21 CPIC genes; full clinical reporting | Phased or unphased VCF | High; CPIC reference | CYP2D6 SV-rich samples need Cyrius/StellarPGx upstream |
| **Cyrius** (Chen 2021 *Genome Med*) | **Yes (99.3% concordance)** | **Yes** | CYP2D6 only | Phased haplotypes | GeT-RM 99.3% | Other genes (single-purpose tool) |
| **BCyrius** (PubMed 39901590, 2025) | Yes (extended) | Yes | CYP2D6 only | Phased | Extended SV diversity | Other genes |
| **Aldy v4** (Numanagic 2018 *Nat Commun*) | Yes | Yes | CYP2D6, CYP2A6, CYP2B6, etc. | Phased | GeT-RM 82-87% (CYP2D6) | Less accurate than Cyrius for CYP2D6 |
| **Stargazer** (Lee 2019 *Genet Med*) | Limited | Yes | ~50 PGx genes | Statistical phasing | ~84% (CYP2D6) | Fails on rare alleles; statistical phasing is unstable |
| **StellarPGx** | Yes (~99%) | Yes | CYP2D6 + others | Phased | GeT-RM ~99% | Less widely deployed than Cyrius |
| **Astrolabe** (proprietary, formerly Constellation) | Yes | Yes | Multi-gene | Proprietary | Industry-validated | License required |
| **StarPhase** (PacBio HiFi 2024+) | Yes | Yes | All CPIC Level A genes + HLA | Native phasing | Long-read gold standard | Requires PacBio HiFi |

**Canonical clinical workflow 2024-2026:** PharmCAT for the panel + Cyrius (or StellarPGx) for CYP2D6 SVs + dedicated HLA typer (T1K, OptiType, HLA-LA) for HLA.

Twesigomwe 2020 *npj Genom Med*: inter-tool discordance 10-18% on CYP2D6; nearly all in samples carrying SVs.

## HLA-Drug Associations: Mechanistically Distinct from CYP

HLA associations are **idiosyncratic immune reactions**, not dose-response phenomena. Effect sizes (OR 50-1000+) far exceed any CYP polymorphism. Testing rationale is **screen-and-avoid**, not dose-adjust.

| Allele | Drug | Reaction | Population | Landmark |
|--------|------|----------|------------|----------|
| **HLA-B\*57:01** | Abacavir | HSS | All ancestries (5-8% NFE) | Mallal 2008 *NEJM* (PREDICT-1) |
| **HLA-B\*15:02** | Carbamazepine, oxcarbazepine, phenytoin, lamotrigine (weaker) | SJS/TEN | Han Chinese, Thai, Malay, Indian (>=5%) | Chung 2004 *Nature*; FDA black-box 2007 |
| **HLA-A\*31:01** | Carbamazepine | DRESS, MPE, SJS/TEN | Europeans (2-5%), Japanese | McCormack 2011 *NEJM* |
| **HLA-B\*58:01** | Allopurinol | SJS/TEN, DRESS | Han Chinese (10-15%), Thai, Korean | Hung 2005 *PNAS* (OR ~580) |
| **HLA-B\*13:01** | Dapsone | DDS | Han Chinese, SE Asian | Zhang 2013 *NEJM* |
| **HLA-B\*35:02** (NOT \*35:01) | Minocycline | DILI | All | Urban 2017 *J Hepatol* |
| **HLA-B\*35:01** | TMP-SMX | DILI, DRESS-like | Mixed | Li 2021 *Hepatology* |
| **HLA-B\*14:01** | TMP-SMX | DILI | African | Li 2021 |
| **HLA-A\*33:01/03** | Terbinafine | DILI | Multi-ancestry | Nicoletti 2017 |
| **HLA-DRB1\*15:01-DQB1\*06:02 haplotype** | Amoxicillin-clavulanate | DILI | Europeans | Stephens 2013 |
| **HLA-B\*15:13** | Phenytoin | SJS | Malaysian | Chang 2017 |

**Critical:** HLA screening requires 4-field resolution. \*57:01 (abacavir risk) vs \*57:03 (no risk); \*35:02 (minocycline DILI) vs \*35:01 (TMP-SMX DILI). See `clinical-databases/hla-typing` for typing.

## Non-CYP Pharmacogenes: Variant-Level Detail

### DPYD (5-FU / Capecitabine / Tegafur); 2024 Activity Score Framework

CPIC 2024 update (Lam et al. *Clin Pharmacol Ther*) moves DPYD to a full **gene activity score** system. Activity values: normal-function = 1.0, decreased = 0.5, no function = 0.

| Variant | rsID | Allele | Activity |
|---------|------|--------|----------|
| c.1905+1G>A | rs3918290 | DPYD*2A | 0 (splice disruption) |
| c.1679T>G | rs55886062 | DPYD*13 (p.I560S) | 0 |
| c.2846A>T | rs67376798 | (p.D949V) | 0.5 |
| c.1129-5923C>G / c.1236G>A (HapB3) | rs56038477 / rs75017182 | HapB3 | 0.5 |

Gene AS = sum of two lowest activities. Recommended dose: AS 2 = full dose; AS 1.5 = 50% start + TDM; AS 1.0 = 50% start + TDM; AS 0 = avoid.

c.85T>C (DPYD\*9A) is NOT in the CPIC 2024 actionable set despite frequent commercial reporting; evidence does not support clinical decrement.

EU universal pre-treatment testing standard since Henricks 2018 *Lancet Oncol* (70% reduction in grade >=3 fluoropyrimidine toxicity) and EMA 2020 endorsement. US lags; ASCO/NCCN moved 2022-2024.

### TPMT + NUDT15 (Thiopurines); 2025 Update

Maillard 2026 *Clin Pharmacol Ther* update emphasizes greater dose reduction for **compound TPMT/NUDT15 IM**.

| Gene | Variant | Activity | Population |
|------|---------|----------|-----------|
| TPMT *2 | c.238G>C | 0 | -- |
| TPMT *3A | c.460G>A + c.719A>G | 0 | EUR-common |
| TPMT *3B | c.460G>A | 0 | -- |
| TPMT *3C | c.719A>G | 0 | AFR / EAS dominant |
| NUDT15 *3 | c.415C>T (rs116855232) | 0 | 9.8% Han Chinese; <1% EUR |

NUDT15 *3 is the **dominant thiopurine determinant in East Asians**; TPMT-alone testing misses these patients (Yang 2014 *Nat Genet*).

### UGT1A1 (Irinotecan, Atazanavir)

- \*28 (TA7 promoter repeat vs \*1 = TA6, \*37 = TA8); EUR-common
- \*6 (c.211G>A, p.G71R); East Asian dominant
- Severe neutropenia in \*28/\*28 at irinotecan >=180 mg/m^2

### CYP2C19 + Clopidogrel; The Most-Litigated Pair

- **Pare 2010** *NEJM*: no benefit of clopidogrel in \*2 carriers in CURE/ACTIVE-A.
- **TAILOR-PCI** (Pereira 2020 *JAMA*): 5,302 patients post-PCI; primary endpoint MACE @12mo HR 0.66, **p=0.06 (negative by pre-specified alpha)** but positive in sensitivity analyses.
- **Pereira NL et al 2021 meta-analysis** (7 RCTs, 6,409 patients): ~30% MACE reduction. Likely venue *JAMA Cardiology* 6:e215028 — verify exact citation in current literature.
- **Consensus 2024 (ACC/AHA/ESC):** genotype-guided therapy reasonable; strongest in post-PCI ACS.

### Warfarin (CYP2C9 + VKORC1 + CYP4F2)

- **EU-PACT 2013** *NEJM*: PGx dosing positive (European).
- **COAG 2013** *NEJM*: PGx dosing negative; worse in African Americans because algorithm omitted CYP2C9 \*5, \*6, \*8, \*11 alleles common in African ancestry. **Paradigmatic ancestry-algorithm failure** (Daneshjou 2014 *Blood*).
- IWPC algorithm explains 47-55% of dose variance.

### SLCO1B1 + Simvastatin

- rs4149056 (c.521T>C, p.V174A); OR 4.5 per C allele for myopathy on 80 mg simvastatin (SEARCH 2008 *NEJM*).
- 2022 CPIC update broadened to all statins with SLCO1B1 substrate behavior.

### Other Actionable

- **CYP2B6 *6** (c.516G>T + c.785A>G): efavirenz dose 600 -> 400 mg in *6/*6 (ENCORE1).
- **CYP3A5 *3** (rs776746): non-expressers (\*3, \*6, \*7) are the *common* state in non-AFR; expressers need 1.5-2x higher tacrolimus dose.
- **G6PD** (CPIC 2022 Gammal 2023): X-linked; female heterozygotes have mosaic activity that single-timepoint assay misclassifies.

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Multi-gene PGx panel from VCF | PharmCAT | CPIC-recommended; 21 genes + full clinical reporting |
| CYP2D6 with structural variants | Cyrius (or StellarPGx) | Only tools with reliable SV calling from short-read |
| All CPIC Level A + HLA from one sample | PacBio HiFi + StarPhase | Long-read single-pass typing |
| Pre-emptive panel for cohort | PREPARE-style 12-gene panel | Swen 2023 RCT-validated |
| HLA-B\*57:01 abacavir screen | T1K or OptiType (4-field); HIBAG if SNP-array | Need 4-field specificity |
| African-ancestry warfarin | IWPC algorithm + CYP2C9 *5/*6/*8/*11 explicit | COAG failure paradigm |
| East Asian thiopurine | NUDT15 + TPMT | NUDT15 *3 is dominant in EAS |
| Compound IM (TPMT + NUDT15) | Apply 2025 update | More aggressive dose reduction than single-gene IM |
| Activity score interpretation | Caudle 2020 thresholds for CYP2D6; gene-specific for others | Per CPIC |

## PharmCAT Workflow (Recommended Multi-Gene Pipeline)

**Goal:** Generate CPIC-compliant pharmacogenomic report from a phased or unphased VCF covering 21 PGx genes.

**Approach:** Run PharmCAT on the VCF; supplement CYP2D6 with Cyrius output if SVs suspected; cross-reference HLA from separate typing.

```bash
# PharmCAT (CPIC-recommended; covers 21 genes including CYP2C19, CYP2C9, CYP2D6,
# DPYD, TPMT, NUDT15, UGT1A1, SLCO1B1, CYP3A5, CYP4F2, VKORC1, IFNL3/IFNL4, etc.)

# 1. Preprocess VCF (ensures correct ref allele alignment + chr formatting)
pharmcat_vcf_preprocessor.py \
    -vcf input.vcf.gz \
    -refFna GRCh38.fa \
    -o pharmcat_input/

# 2. Run PharmCAT
java -jar pharmcat.jar \
    -vcf pharmcat_input/input.preprocessed.vcf.bgz \
    -o pharmcat_output/

# Output: <sample>.report.html with phenotype, activity score, dosing recommendations
```

For CYP2D6 SV-rich samples, run Cyrius separately and pass outside calls to PharmCAT:

```bash
# Cyrius for CYP2D6 (99.3% concordance vs Aldy 82-87%, Stargazer 84%)
cyrius -m sample.bam -o cyrius_out --threads 8
# Output: cyrius_out/sample.tsv with diplotype + activity score

# Pass outside calls to PharmCAT
java -jar pharmcat.jar \
    -vcf pharmcat_input/input.preprocessed.vcf.bgz \
    -po cyrius_out/cyrius_for_pharmcat.tsv \
    -o pharmcat_output_with_cyrius/
```

## CYP2D6 Activity Score Calculation

**Goal:** Convert CYP2D6 diplotype to activity score and phenotype with Caudle 2020 conventions.

**Approach:** Look up per-allele activity values; handle copy-number duplications; apply Caudle 2020 phenotype bins.

```python
# Caudle 2020 activity values; *10 reset from 0.5 to 0.25 in 2020
CYP2D6_ACTIVITY = {
    '*1': 1.0, '*2': 1.0, '*35': 1.0,
    '*3': 0.0, '*4': 0.0, '*5': 0.0, '*6': 0.0, '*7': 0.0, '*8': 0.0,
    '*11': 0.0, '*12': 0.0, '*15': 0.0, '*19': 0.0, '*20': 0.0,
    '*36': 0.0, '*40': 0.0, '*42': 0.0, '*68': 0.0,
    '*9': 0.5, '*41': 0.5, '*17': 0.5, '*29': 0.5,
    '*10': 0.25,
    '*13': 0.0,
}


def cyp2d6_activity(diplotype):
    '''Convert CYP2D6 diplotype to activity score.

    Accepts e.g. '*1/*4' or '*2xN/*10' or '*4xN/*10'. Copy-number-aware:
    - *4xN is clinically silent (no-function * N = 0)
    - *1xN, *2xN multiply functional activity
    '''
    left, right = diplotype.split('/')
    return _allele_activity(left) + _allele_activity(right)


def _allele_activity(allele_str):
    '''Handle copy-number suffix xN. *4xN remains 0 (the most common mis-classification).'''
    if 'x' in allele_str:
        base, n = allele_str.split('x')
        copies = int(n) if n != 'N' else 2  # 'N' usually >=2; clinical assumes 2 unless quantified
        return CYP2D6_ACTIVITY.get(base, 1.0) * copies
    return CYP2D6_ACTIVITY.get(allele_str, 1.0)


def cyp2d6_phenotype(activity_score):
    '''Caudle 2020 phenotype bins.'''
    if activity_score == 0:
        return 'Poor Metabolizer'
    if activity_score < 1.25:
        return 'Intermediate Metabolizer'
    if activity_score <= 2.25:
        return 'Normal Metabolizer'
    return 'Ultrarapid Metabolizer'


# Example: *4xN/*10; the classic clinical-silence footgun
diplotype = '*4xN/*10'
score = cyp2d6_activity(diplotype)  # 0 (from *4xN) + 0.25 (from *10) = 0.25
print(f'{diplotype}: AS={score}, phenotype={cyp2d6_phenotype(score)}')  # IM, NOT UM
```

## DPYD Activity Score (CPIC 2024)

```python
DPYD_2024_ACTIVITY = {
    'c.1905+1G>A': 0.0,    # *2A; splice donor
    'c.1679T>G': 0.0,      # *13; p.I560S
    'c.2846A>T': 0.5,      # p.D949V
    'HapB3': 0.5,          # c.1129-5923C>G linked with c.1236G>A
}


def dpyd_activity(variants):
    '''Compute DPYD gene activity score from observed variants.

    Sum the two lowest activities across the two alleles. CPIC 2024 dosing:
    - AS 2.0: full dose
    - AS 1.5: 50% start + TDM
    - AS 1.0: 50% start + TDM
    - AS 0.0: avoid
    '''
    activities = sorted([DPYD_2024_ACTIVITY.get(v, 1.0) for v in variants])
    return sum(activities[:2])


def dpyd_dosing(activity_score):
    if activity_score >= 1.99:
        return 'Full dose'
    if activity_score >= 1.0:
        return '50% starting dose + therapeutic drug monitoring'
    return 'Avoid fluoropyrimidines'
```

## PharmGKB API for Drug-Gene Pair Lookup

```python
import requests

PHARMGKB = 'https://api.pharmgkb.org/v1'


def clinical_annotation(gene_symbol):
    '''Query PharmGKB clinical annotations by gene.'''
    r = requests.get(f'{PHARMGKB}/data/clinicalAnnotation',
                     params={'view': 'base', 'location.genes.symbol': gene_symbol},
                     timeout=30)
    return r.json().get('data', [])


def cpic_guideline(gene_symbol):
    '''Query CPIC guidelines via PharmGKB.'''
    r = requests.get(f'{PHARMGKB}/data/guideline',
                     params={'view': 'base', 'relatedGenes.symbol': gene_symbol, 'source': 'CPIC'},
                     timeout=30)
    return r.json().get('data', [])
```

## Per-Operation Failure Modes

**1. *4xN -> "Ultrarapid Metabolizer"**
- Trigger: Pipeline reports CYP2D6 \*4xN as UM.
- Mechanism: \*4 has activity 0; \*4 x N = still 0. Only functional alleles (\*1, \*2, \*35) become UM when amplified.
- Symptom: Patient labeled as needing dose reduction when they should be PM/IM.
- Fix: Look up per-allele activity BEFORE multiplying by N; \*4xN = 0; AS depends entirely on the other allele.

**2. Calling CYP2D6 from short-read without SV-aware tool**
- Trigger: Use GATK + PharmCAT only on CYP2D6.
- Mechanism: 97%+ CYP2D6/CYP2D7 identity; SVs (deletion, duplications, hybrids) silently miscalled.
- Symptom: ~10-18% of samples miscalled (Twesigomwe 2020); concentrated in samples with SVs.
- Fix: Add Cyrius (or StellarPGx) for CYP2D6; pass outside calls to PharmCAT.

**3. Pre-2020 \*10 activity value**
- Trigger: Use activity = 0.5 for CYP2D6 \*10.
- Mechanism: Caudle 2020 reset \*10 from 0.5 to 0.25 based on metabolic-ratio evidence.
- Symptom: East-Asian samples mis-classified as NM (when should be IM).
- Fix: Use Caudle 2020 activity table; \*10 = 0.25.

**4. EUR-only DPYD panel**
- Trigger: Pre-treat fluoropyrimidine using CPIC-core 4-variant panel only.
- Mechanism: 4-variant panel captures EUR DPD-deficient carriers; misses ~30% of AFR carriers (Offer 2014).
- Symptom: African-ancestry patients suffer severe toxicity despite "negative" PGx.
- Fix: Use extended panel for AFR cohorts; supplement with phenotype testing (uracil/dihydrouracil plasma ratio).

**5. TPMT testing without NUDT15**
- Trigger: Pre-treat thiopurines using TPMT-only PGx in East Asian patient.
- Mechanism: NUDT15 *3 (9.8% EAS, <1% EUR) is the dominant determinant in EAS.
- Symptom: EAS patients TPMT-wildtype suffer severe myelosuppression.
- Fix: Always test NUDT15 alongside TPMT; apply Maillard 2026 compound-IM rules.

**6. HLA-B\*57 -> "abacavir risk" (4-field underspecified)**
- Trigger: Screen reports "B*57 present" as contraindication.
- Mechanism: B*57:01 (HSS risk), B*57:02, B*57:03 (no HSS risk).
- Symptom: False contraindication; patient denied effective therapy.
- Fix: Report 4-field; B*57:01 specifically.

**7. CYP3A5 *3 / non-expresser confusion**
- Trigger: Apply "CYP3A5 normal metabolizer" to *3/*3 in tacrolimus dosing.
- Mechanism: *3/*3 are NON-EXPRESSERS (most common state in non-AFR); expressers (any *1) need 1.5-2x higher dose.
- Symptom: Tacrolimus over-dosing in expressers; under-dosing in non-expressers.
- Fix: Apply CPIC 2015 (Birdwell) tacrolimus dosing; flag expresser status.

**8. Activity-based vs allele-based confusion**
- Trigger: Sum activities across substrate-non-specific assumption for *17.
- Mechanism: CYP2D6 *17 metabolizes codeine poorly but moderately for other substrates (Twesigomwe 2020).
- Symptom: Substrate-specific dose recommendations applied generically.
- Fix: Use substrate-specific guidance where available; flag *17 in AFR cohorts.

## Reconciliation: When Tools Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| Cyrius vs Aldy CYP2D6 disagree | SV-rich sample; Aldy less accurate | Trust Cyrius |
| PharmCAT vs CPIC website disagree on phenotype | PharmCAT version lag or *10 activity value drift | Update PharmCAT to current release |
| CPIC vs DPWG dosing differ | Independent guideline bodies | Cite both; use jurisdiction-appropriate one |
| Patient phenotype doesn't match genotype | Drug-drug interaction; clearance physiology; non-pharmacogenetic factor | Consider phenoconversion; clinical reassessment |
| TPMT-only test vs IM phenotype | Missed NUDT15 in EAS | Re-test with NUDT15 |
| *4xN reported as UM | Tool bug | Use SV-aware tool and Caudle 2020 activity table |
| HLA-B*57 reported without 4-field | Insufficient resolution | Re-type at 4-field minimum |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| Cyrius CYP2D6 accuracy | 99.3% on GeT-RM reference samples | Chen 2021 *Genome Med* |
| Aldy CYP2D6 accuracy | 82-87% on GeT-RM | Twesigomwe 2020 |
| Stargazer CYP2D6 accuracy | ~84% on GeT-RM | Twesigomwe 2020 |
| Inter-tool CYP2D6 discordance | 10-18% (concentrated in SV samples) | Twesigomwe 2020 |
| PREPARE ADR reduction | OR 0.70 (95% CI 0.54-0.91) for actionable interactions | Swen 2023 *Lancet* |
| PREPARE actionable variant rate | 93.5% of patients had >=1 actionable variant | Swen 2023 |
| TAILOR-PCI primary endpoint | HR 0.66 (95% CI 0.43-1.02), p=0.06 (negative) | Pereira 2020 *JAMA* |
| Pereira 2021 meta-analysis | ~30% MACE reduction with genotype-guided clopidogrel | Pereira NL et al 2021 (cardiology meta-analysis; verify exact venue against current literature — likely *JAMA Cardiology*) |
| Henricks 2018 DPYD outcome | 70% reduction in grade >=3 fluoropyrimidine toxicity | Henricks 2018 *Lancet Oncol* |
| NUDT15 *3 frequency | 9.8% Han Chinese vs <1% EUR | Yang 2014 *Nat Genet* |
| HLA-B*57:01 OR for abacavir HSS | ~100 (NNT 13) | Mallal 2008 *NEJM* (PREDICT-1) |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| CYP2D6 reported as UM in samples with *4xN | Tool not SV-aware OR Caudle 2020 not applied | Use Cyrius; check *4xN handling |
| East-Asian patient labeled CYP2D6 NM | *10 still at activity 0.5 | Update activity table to Caudle 2020 (*10 = 0.25) |
| African patient suffers warfarin bleeding despite "wildtype" CYP2C9 | Panel omits *5/*6/*8/*11 (AFR-common) | Use ancestry-aware panel; supplement with INR-guided dosing |
| Severe thiopurine toxicity in TPMT-wildtype EAS patient | NUDT15 not tested | Always pair TPMT + NUDT15 |
| Patient with CYP2C19 *2/*2 and clopidogrel failure | Expected; no genotype-guided alternative chosen | Switch to prasugrel/ticagrelor per CPIC |
| DPYD AS = 0 but no dose adjustment | Pre-2024 single-variant rule used instead of activity score | Update to CPIC 2024 activity-score framework |
| HLA-B*57:01 false positive | 2-field B*57 result misinterpreted | Re-type at 4-field |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "TAILOR-PCI missed primary endpoint; why genotype clopidogrel?" | Sensitivity analyses positive; Pereira 2021 meta-analysis (7 RCTs) +30% MACE reduction; ESC 2023 endorses; ACC 2022 weaker. |
| "DPYD universal screening is expensive" | Henricks 2018 70% toxicity reduction + Knikman 2021 cost-effective; EU standard since 2020; US ASCO/NCCN updated 2022-2024. |
| "CYP2D6 SV calling is unreliable" | Cyrius 99.3% on GeT-RM (Chen 2021); not unreliable; the prior tools were. |
| "*10 = 0.25 disagrees with old paper" | Caudle 2020 *Clin Transl Sci* consensus reset based on substrate-metabolic-ratio evidence. |
| "GeneSight is approved by my hospital" | GUIDED trial (Greden 2019) missed primary endpoint; physician-unblinded; literature shows modest effects inseparable from expectancy bias. |
| "Why pair TPMT + NUDT15?" | NUDT15 *3 is the dominant thiopurine determinant in East Asians (9.8% vs TPMT *3C ~2%); compound IM (TPMT + NUDT15) requires more aggressive dose reduction per Maillard 2026. |
| "HLA imputation from SNP array reliable?" | EUR-trained panel on EUR samples ~95%; cross-ancestry drops to 70-80%; for HSCT use sequencing-based typing. |

## References

- Sangkuhl K et al. 2020. Pharmacogenomics Clinical Annotation Tool (PharmCAT). *Clin Pharmacol Ther* 107:203.
- Chen X et al. 2021. Cyrius: accurate CYP2D6 genotyping from WGS. *Genome Med* 13:35.
- Numanagic I et al. 2018. Allelic decomposition and exact genotyping of highly polymorphic and structurally variant genes. *Nat Commun* 9:828. (Aldy)
- Lee SB et al. 2019. Stargazer: a tool for calling star alleles. *Genet Med* 21:361.
- Twesigomwe D et al. 2020. Systematic comparison of pharmacogene haplotype callers. *npj Genom Med* 5:36.
- Caudle KE et al. 2020. Standardizing CYP2D6 genotype to phenotype translation. *Clin Transl Sci* 13:116. (Activity-score reset for *10)
- Bank PCD et al. 2018. Comparison of dosing recommendations of pharmacogenomic-based guidelines (CPIC vs DPWG). *Clin Pharmacol Ther* 105:951.
- Swen JJ et al. 2023. PREPARE: A pre-emptive pharmacogenetic testing strategy. *Lancet* 401:347.
- Henricks LM et al. 2018. DPYD-guided dose individualization to fluoropyrimidines. *Lancet Oncol* 19:1459.
- Pereira NL et al. 2020. Effect of genotype-guided oral P2Y12 inhibitor selection vs conventional clopidogrel therapy on ischemic outcomes after PCI. *JAMA* 324:761. (TAILOR-PCI)
- Pereira NL et al. 2021. CYP2C19 genetic testing for oral P2Y12 inhibitor therapy: a meta-analysis of 7 RCTs. *Lancet* (or *J Am Coll Cardiol*; cited).
- Mallal S et al. 2008. HLA-B*5701 screening for hypersensitivity to abacavir. *NEJM* 358:568. (PREDICT-1)
- Chung WH et al. 2004. Medical genetics: a marker for Stevens-Johnson syndrome. *Nature* 428:486.
- McCormack M et al. 2011. HLA-A*3101 and carbamazepine-induced hypersensitivity reactions in Europeans. *NEJM* 364:1134.
- Hung SI et al. 2005. HLA-B*5801 allele as a genetic marker for severe cutaneous adverse reactions caused by allopurinol. *PNAS* 102:4134.
- Yang JJ et al. 2014. Inherited NUDT15 variant is a genetic determinant of mercaptopurine intolerance. *Nat Genet* 46:1017.
- PharmCAT documentation: `https://pharmcat.org`
- PharmVar: `https://www.pharmvar.org`
- CPIC: `https://cpicpgx.org`
- DPWG: `https://www.knmp.nl/dpwg`

## Related Skills

- clinical-databases/hla-typing - HLA-B*57:01, B*15:02, B*58:01, A*31:01 typing
- clinical-databases/clinvar-lookup - Variant pathogenicity for non-PGx context
- clinical-databases/variant-prioritization - Rare-disease pipeline
- clinical-databases/myvariant-queries - Aggregated PGx variant annotation
- chemoinformatics/admet-prediction - Drug metabolism prediction
