---
name: bio-clinical-biostatistics-cdisc-data
description: Reads, validates, and prepares CDISC SDTM and ADaM clinical trial data for analysis. Covers SDTM domain joins (DM, AE, EX, VS, LB, DS), ADaM architecture (ADSL, BDS, OCCDS, ADTTE) with traceability, treatment-emergent AE conventions, baseline derivation, SUPPQUAL/NSV handling, Define-XML 2.1, and Pinnacle 21 / CORE validation. Use when working with clinical trial datasets in CDISC SDTM/ADaM format, preparing analysis-ready data, or validating for regulatory submission.
tool_type: python
primary_tool: pyreadstat
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: pyreadstat 1.2+, pandas 2.1+, numpy 1.26+. CDISC standards referenced: SDTM 2.0 / SDTMIG 3.4 (SDTM 3.0 / SDTMIG 4.0 in public review through April 2026); ADaMIG v1.3 (2021); OCCDS v1.1 (Nov 2021); BDS-for-TTE v1.0; Define-XML 2.1 (FDA-recommended for studies starting on/after March 15, 2023); Dataset-JSON v1.1 (Dec 2025; FDA Federal Register notice April 2025); Pinnacle 21 Community 4.0+; CORE (CDISC Open Rules Engine, 2024).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R packages cited (essential for ADaM derivation): admiral (Roche/openpharma), metacore, metatools, xportr

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# CDISC SDTM and ADaM Data Handling

**"Load clinical trial data"** -> Parse CDISC SDTM domain files; build or consume ADaM analysis-ready datasets; preserve subject-level and event-level structure; respect traceability and validation expectations for regulatory submission.
- Python: `pyreadstat.read_xport()`, `pd.read_sas()`, `pd.merge()`
- R: `haven::read_xpt()`, `admiral` for ADaM derivation, `Pinnacle21` or `CORE` for validation

## Aggregation Strategy Taxonomy -- Choose the Right Question

| Strategy | Scientific question answered | Example endpoint | Fails when |
|----------|------------------------------|------------------|------------|
| Any event (binary) | Does treatment change probability of experiencing the event at all? | Had any serious AE: Yes/No | Treatment changes event burden but not anyone-event probability |
| Event count | Does treatment change burden of events per patient? | Total AE count per subject | Subjects with 1 vs 10 events treated equivalently |
| Maximum severity | Does treatment shift patients toward more severe manifestations? | Worst AESEV per subject | Confounded with event count (more events -> higher chance of severe) |
| First event + time | Does treatment delay onset of the event? | Time to first serious AE (TTE) | Multiple events per subject ignored |
| Rate (events per person-time) | What is the per-time-unit rate? | AEs per subject-year | Requires exposure-time tracking; differential dropout biases rates |
| Composite (per ICH E9 R1) | Event becomes part of endpoint definition | Death = treatment failure | Direction of components conflict; needs hierarchy |

These are NOT interchangeable. A drug might not change the proportion with AEs (binary: no effect) but increase events per patient (count: harmful). **The choice must be pre-specified in the SAP based on the scientific question, not analytic convenience.**

## Decision Tree by Scenario

| Scenario | Recommended aggregation | Why |
|----------|------------------------|-----|
| Primary safety endpoint, single SAE event | Any event (binary); analyse with logistic | Standard regulatory; cite FDA Safety Reporting Guidance |
| Total adverse-event burden across study | Event count per subject; analyse with Poisson or negative binomial | Captures all events; sandwich SE recommended |
| Toxicity grade comparison across arms | Max severity per subject; ordinal logistic with PO check | Preserves grade ordering; cite Brant test for PO |
| Time-to-first AE (Kaplan-Meier visualisation) | First event + time; censor non-events | See clinical-biostatistics/survival-analysis |
| Rate of exacerbations per patient-year | Rate via negative binomial with offset for exposure time | Standard in COPD/asthma trials |
| Composite endpoint (e.g., MACE) | Component-level definition with hierarchy | Pre-specify per ICH E9(R1) composite strategy |
| Stratification factor extraction | Use STRATA1, STRATA2 from RANDB or DM SUPP | Must appear in analysis (Kahan-Morris 2012) |
| Baseline value derivation | VSBLFL='Y' / LBBLFL='Y'; derive latest pre-dose only if flag missing | Trust SDTM flag when present |

## SDTM vs ADaM -- The Regulatory Layer Cake

| Layer | Standard | Purpose | Granularity | Examples |
|-------|----------|---------|-------------|----------|
| Source CRF | EDC system | Raw data capture | Form/page | Rave, Medidata, Veeva |
| SDTM | SDTM 2.0 / SDTMIG 3.4 | Tabulation; "what happened" | One row per observation | DM, AE, EX, VS, LB, DS |
| ADaM | ADaMIG v1.3 (2021); v3.0 in development | Analysis-ready; "one-PROC-away from CSR table" | Subject (ADSL), parameter-timepoint (BDS), occurrence (OCCDS) | ADSL, ADAE, ADLB, ADTTE, ADRS |
| TLF | Sponsor SAS / R / Python | Tables, listings, figures for CSR | Output | Statistical methods section, demographic table, primary efficacy |

**The ADaM Fundamental Principles (the "ROT" document):**

1. Analysis-ready (one procedure call -> the analysis result)
2. Traceability (every value links back to SDTM via metadata)
3. Clear/unambiguous communication via Define-XML
4. Naming conventions (PARAM, PARAMCD, AVAL, AVALC, BASE, CHG, PCHG, ABLFL, ANL01FL, ...)
5. Structural rules (ADSL one row per subject; BDS one row per subject/parameter/timepoint/analysis flag)

**Postdoc reading:** the ADaM IG v1.3 PDF (cdisc.org), ADaM ROT, FDA Study Data Technical Conformance Guide (current 2024 version), Pinnacle 21 validation rule catalog, PHUSE Connect 2023-2025 conference proceedings.

## SDTM Domain Overview

| Domain | Level | Description | Key Variables |
|--------|-------|-------------|---------------|
| DM | Subject | Demographics (one row per subject) | USUBJID, ARM, ARMCD, ACTARM, ACTARMCD, AGE, SEX, RACE, RFSTDTC, RFXSTDTC, RFENDTC |
| AE | Event | Adverse events (multiple per subject) | USUBJID, AETERM, AEDECOD, AEBODSYS, AESEV, AESER, AESTDTC, AEENDTC |
| EX | Event | Drug exposure/dosing | USUBJID, EXTRT, EXDOSE, EXSTDTC, EXENDTC |
| VS | Event | Vital signs | USUBJID, VSTESTCD, VSSTRESN, VSBLFL, VISIT |
| LB | Event | Lab results | USUBJID, LBTESTCD, LBSTRESN, LBSTRESC, LBORRES, LBBLFL, LBSPEC |
| DS | Event | Disposition | USUBJID, DSDECOD, DSSTDTC |
| SE | Event | Subject elements (treatment epochs) | USUBJID, ETCD, SESTDTC, SEENDTC |
| MH | Event | Medical history | USUBJID, MHDECOD, MHCAT |
| CM | Event | Concomitant medications | USUBJID, CMDECOD, CMSTDTC, CMENDTC |

USUBJID = STUDYID-SITEID-SUBJID is the universal merge key. Subject-level domains (DM) have one row per USUBJID; event-level domains have multiple.

**ARM vs ACTARM:** ARM is planned treatment from randomisation; ACTARM is actual treatment received. **In crossover designs, ARM differs from ACTARM by definition;** in parallel-arm trials, they diverge when subjects are randomised to one arm but receive another (per-protocol violations). Primary analyses use ARM (ITT); safety uses ACTARM.

**RFSTDTC vs RFXSTDTC:** RFSTDTC is "first study activity date" (typically screening start); RFXSTDTC is "first treatment date." For treatment-emergent adverse event (TEAE) calculations, ALWAYS use RFXSTDTC (per ICH E2A) — RFSTDTC includes screening AEs which are not treatment-emergent.

## Reading .xpt Files

```python
import pyreadstat
import pandas as pd

# pyreadstat (recommended -- handles SAS metadata)
dm, meta = pyreadstat.read_xport('dm.xpt')
# meta.column_names, meta.column_labels, meta.variable_value_labels

# pandas built-in (SAS XPORT v5)
dm = pd.read_sas('dm.xpt', format='xport', encoding='utf-8')

# CSV fallback (common in academic datasets)
dm = pd.read_csv('DM.csv')
```

When pyreadstat is available, the metadata object provides column labels, value labels, and format information lost with other readers. **Critical for analysis-dataset derivation:** the metadata carries the controlled-terminology codelist, essential for handling values like AESEV ('MILD'/'MODERATE'/'SEVERE') with semantic ordering.

## SAS XPT v5 vs Dataset-JSON -- The 2025-2026 Transition

**SAS XPT v5** is the current FDA submission format but dates to 1995, with constraints:

- 8-character variable names (so `LBSTRESN` is a max-length name)
- 200-character text values
- No UTF-8 (ASCII only) -> problematic for multilingual trials
- Single dataset per file

**Dataset-JSON v1.1 (CDISC, December 2025; FDA Federal Register notice April 2025)** is the modern replacement. PHUSE-CDISC-FDA pilot has demonstrated drop-in feasibility. FDA adoption timeline pending as of mid-2026; EMA and PMDA exploring in parallel.

**Pragmatic position:** for the next ~2 years, SAS XPT v5 will remain the de facto submission format; sponsors should architect for Dataset-JSON migration but maintain XPT compliance.

## Joining Domains -- The Right Way

```python
import pandas as pd

dm = pd.read_csv('DM.csv')
ae = pd.read_csv('AE.csv')

# WRONG: merging event-level directly onto subject-level inflates rows
# RIGHT: aggregate first, then merge
any_serious = ae.groupby('USUBJID')['AESER'].apply(lambda x: (x == 'Y').any()).reset_index()
any_serious.columns = ['USUBJID', 'HAD_SERIOUS_AE']

analysis = dm.merge(any_serious, on='USUBJID', how='left')
analysis['HAD_SERIOUS_AE'] = analysis['HAD_SERIOUS_AE'].fillna(False)
```

**Always use `how='left'` when merging onto DM** to preserve all randomised subjects, even those with no events. Fill missing event indicators with 0 or False.

### Aggregation strategy must follow the scientific question

| Strategy | Scientific question | Example |
|----------|---------------------|---------|
| Any event (binary) | Does treatment increase probability of experiencing the event at all? | Had any serious AE: yes/no |
| Event count | Does treatment increase event burden per patient? | Total AE count per subject |
| Maximum severity | Does treatment shift toward more severe manifestations? | Worst AESEV per subject |
| First event + time | Does treatment delay onset? | Time to first serious AE |
| Rate (events per person-time) | What is the per-time-unit rate? | AEs per subject-year |

These are NOT interchangeable. A drug might not change the proportion with AEs (binary: no effect) but increase events per patient (count: harmful). The choice must follow the SAP, not analytic convenience.

```python
# Count events per subject
ae_counts = ae.groupby('USUBJID').size().reset_index(name='AE_COUNT')

# Maximum severity per subject (map to numeric first -- string max is unreliable)
severity_map = {'MILD': 1, 'MODERATE': 2, 'SEVERE': 3}
ae['AESEV_NUM'] = ae['AESEV'].map(severity_map)
max_severity = ae.groupby('USUBJID')['AESEV_NUM'].max().reset_index()

# Specific event: COVID-19 adverse event
covid_ae = ae[ae['AEDECOD'] == 'COVID-19']
covid_ae['AESEV_NUM'] = covid_ae['AESEV'].map(severity_map)
had_covid = covid_ae.groupby('USUBJID')['AESEV_NUM'].max().reset_index()
had_covid.columns = ['USUBJID', 'COVID_SEVERITY']

analysis = dm.merge(had_covid, on='USUBJID', how='left')
analysis['HAD_COVID'] = analysis['COVID_SEVERITY'].notna().astype(int)
```

## ADaM Architecture -- The Postdoc Deep Dive

### ADSL (Subject-Level) -- The Spine

**Exactly one row per subject.** Every other ADaM dataset must merge to ADSL on USUBJID. Standard variables:

- **USUBJID** -- universal subject ID
- **TRT01A / TRT01P / ACTARMCD / ARMCD** -- planned and actual treatment, period 1
- **TRTSDT / TRTEDT** -- treatment start/end dates (derived from EX, not SDTM)
- **AGE, SEX, RACE, ETHNIC** -- demographics from DM
- **RANDDT** -- randomisation date
- **DCSREAS / DCSREASP / DCSREASCD** -- discontinuation reason (coded + verbatim)
- **Population flags:** ITTFL, FASFL, SAFFL, PPROTFL, EFFFL (Y/N flags for analysis populations)
- **Stratification factors:** STRATA1, STRATA2 (from randomisation)
- **Baseline covariates** that will be used as model covariates downstream

### BDS (Basic Data Structure) -- Long Format Analysis Data

**One row per subject per parameter per analysis timepoint per analysis flag.** Used for ADVS, ADLB, ADEFF, ADQS, ADTTE.

Required variables:

- **USUBJID, STUDYID** -- merge keys
- **PARAM, PARAMCD, PARAMN** -- parameter name, code, number
- **AVISIT, AVISITN** -- analysis visit name, number
- **ADT, ADY** -- analysis date, analysis day (relative to TRTSDT)
- **AVAL, AVALC** -- analysis value (numeric, character)
- **BASE** -- baseline value (replicated per subject/parameter)
- **CHG, PCHG** -- change from baseline, percent change
- **ABLFL** -- 'Y' for the baseline record
- **ANL01FL, ANL02FL** -- analysis flags for primary/secondary analyses
- **DTYPE** -- derivation type ('LOCF', 'WOCF', 'AVERAGE', 'BOCF', or null for original)
- **BASETYPE** -- when multiple baselines per subject/parameter (crossover)
- **EPOCH** -- study period (SCREENING, TREATMENT, FOLLOW-UP)

```python
# Example: derive ADLB BDS structure from LB SDTM
import pandas as pd

lb = pd.read_csv('LB.csv')
adsl = pd.read_csv('ADSL.csv')

# Filter to active tests
adlb = lb[lb['LBTESTCD'].isin(['ALT', 'AST', 'CREAT', 'HGB'])].copy()
adlb['AVAL'] = adlb['LBSTRESN']
adlb['PARAM'] = adlb['LBTEST']
adlb['PARAMCD'] = adlb['LBTESTCD']

# Merge subject-level treatment from ADSL
adlb = adlb.merge(adsl[['USUBJID', 'TRT01A', 'TRTSDT']], on='USUBJID')

# Compute analysis day
adlb['ADT'] = pd.to_datetime(adlb['LBDTC'], errors='coerce')
adlb['TRTSDT'] = pd.to_datetime(adlb['TRTSDT'], errors='coerce')
adlb['ADY'] = (adlb['ADT'] - adlb['TRTSDT']).dt.days + 1  # Day 1 = first dose

# Set ABLFL from LBBLFL
adlb['ABLFL'] = (adlb['LBBLFL'] == 'Y').map({True: 'Y', False: None})

# Derive BASE per subject/parameter
baselines = adlb[adlb['ABLFL'] == 'Y'][['USUBJID', 'PARAMCD', 'AVAL']]
baselines.columns = ['USUBJID', 'PARAMCD', 'BASE']
adlb = adlb.merge(baselines, on=['USUBJID', 'PARAMCD'], how='left')

# Compute CHG and PCHG
adlb['CHG'] = adlb['AVAL'] - adlb['BASE']
adlb['PCHG'] = 100 * adlb['CHG'] / adlb['BASE']
```

### OCCDS (Occurrence Data Structure) -- One Row per Event

OCCDS v1.1 (Nov 2021) handles AE, CM, MH. Variables: AEDECOD, AEBODSYS, AESEV, AESER, ASTDT (analysis start date), AENDT, TRTEMFL.

**OCCDS v1.1 added TRTEM01FL through TRTEM##FL** for multi-period treatment-emergent flags — essential for crossover and multi-phase studies where a single TRTEMFL is ambiguous.

### ADTTE (Time-to-Event) -- The CNSR Convention Trap

The **ADaM BDS for TTE v1.0** uses BDS structure with extra variables for survival analysis:

- **STARTDT** -- time origin (typically TRTSDT for OS; RANDDT for PFS; response date for DOR)
- **ADT** -- analysis date (event date if event, censoring date if censored)
- **AVAL** = ADT - STARTDT (+1 if "first day = day 1" convention)
- **AVALU** = 'DAYS' (or 'MONTHS' for some endpoints)
- **CNSR** -- censoring indicator. **CONVENTION: CNSR = 0 for events; positive integers for censoring**, integer encodes censoring reason
- **EVNTDESC** -- text description ('Death due to disease', 'Last alive contact')
- **CNSDTDSC, SRCDOM, SRCVAR, SRCSEQ** -- traceability back to SDTM source

**The CNSR convention is OPPOSITE to most statistical packages**, which use 1 = event. R `survival::Surv(time, event)` expects event=1; SAS PROC LIFETEST takes CENSORED= statement that's opposite to CNSR convention. **This is a perpetual bug source.** When passing ADTTE to analysis:

```python
# Convert CDISC ADTTE CNSR to R/Python statistical convention
adtte['event'] = (adtte['CNSR'] == 0).astype(int)  # 1 = event for survival packages
```

### Define-XML 2.1

Every ADaM dataset requires variable-level metadata in Define-XML 2.1 (FDA-recommended for studies starting on/after March 15, 2023). Fields per variable:

- **Origin** -- CRF, derived, predecessor SDTM variable
- **Derivation rule** -- free text or controlled algorithm
- **Codelist** -- linked controlled terminology
- **Length, datatype, label**

The FDA reviewer's Analysis Data Reviewer's Guide (ADRG) is now expected in every NDA/BLA — walks reviewer through how each analysis dataset was built.

**Two-level traceability expectation:** SDTM raw -> ADaM analysis-ready, with no orphan derivations. FDA reviewers explicitly trace AE counts in CSR table -> ADAE rows -> AE SDTM rows. Any break is a flag.

## Treatment-Emergent AE -- The Convention Variation

**ICH E2A (1995)** defines an AE generically. **TEAE is sponsor-defined:**

```
TRTEMFL = 'Y' if AE.ASTDT >= TRTSDT AND AE.ASTDT <= TRTEDT + X days
```

X = post-treatment follow-up window. Common values:

- Small molecules: 28 or 30 days
- Biologics with extended half-life: longer (e.g., 60-90 days for mAbs)
- Cell/gene therapy: indefinite (lifelong monitoring expected)

**Sponsor variation:**

- Day-of-first-dose AE included as TEAE (FDA preference) vs excluded (some EMA reviewers)
- Partial-date imputation: impute day 15 if only month/year known, vs censor as missing
- Worsening of pre-existing AE: flagged via SEV change vs requires new PT (preferred term)

**MedDRA SOC/PT hierarchy:** AEs coded to MedDRA Preferred Terms (PT), grouped by System Organ Class (SOC). Clinically related PTs (e.g., 'Diarrhea' / 'Frequent bowel movements' / 'Loose stools') often combined via Standardized MedDRA Queries (SMQs) or sponsor-defined groupings. ADAE typically carries both AEDECOD (PT) and SMQ/group flags.

## Baseline Derivation

**ABLFL = 'Y'** marks the record whose AVAL becomes BASE for all other records of the same subject/parameter.

**Standard rule:** last non-missing assessment on or before first dose (TRTSDT). If protocol mandates a specific baseline visit ('Day 1 pre-dose'), that visit's record is flagged.

```python
# Derive ABLFL when SDTM baseline flag is missing/inconsistent
import pandas as pd

vs['VSDTC_dt'] = pd.to_datetime(vs['VSDTC'], errors='coerce')
vs = vs.merge(adsl[['USUBJID', 'TRTSDT']], on='USUBJID')
vs['is_pre_treatment'] = vs['VSDTC_dt'] <= pd.to_datetime(vs['TRTSDT'])

# Latest pre-treatment value per subject/parameter
baseline_records = (vs[vs['is_pre_treatment'] & vs['VSSTRESN'].notna()]
                    .sort_values('VSDTC_dt')
                    .groupby(['USUBJID', 'VSTESTCD'])
                    .tail(1))
baseline_records['derived_ABLFL'] = 'Y'
```

**Critical detail:** filter on VSBLFL='Y' (or LBBLFL='Y') as the primary source. Only fall back to derivation when the flag is missing. Trust the SDTM flag when present; CRF-level baseline designation embeds clinical judgement the analyst cannot reconstruct.

**BASETYPE** required when more than one baseline exists per subject/parameter (crossover studies, multi-period trials). Distinguishes "Period 1 Baseline" vs "Period 2 Baseline."

**DTYPE** values per CDISC controlled terminology: 'LOCF' (last observation carried forward), 'WOCF' (worst), 'AVERAGE', 'BOCF' (baseline observation carried forward), null for original. **Pinnacle 21 flags any DTYPE value not in CT.**

## SUPPQUAL and the NSV Transition

**SUPPQUAL (supplemental qualifiers)** is the legacy mechanism for sponsor-defined variables that don't fit standard SDTM domain columns. Long-format QNAM/QVAL pairs:

```python
supp = pd.read_sas('suppae.xpt', format='xport', encoding='utf-8')
supp_pivot = supp.pivot_table(
    index='USUBJID', columns='QNAM', values='QVAL', aggfunc='first'
).reset_index()
ae_enriched = ae.merge(supp_pivot, on='USUBJID', how='left')
```

For record-level SUPPQUAL (where IDVAR and IDVARVAL identify specific rows):

```python
supp_record = supp[supp['IDVAR'] == 'AESEQ'].copy()
supp_record['AESEQ'] = supp_record['IDVARVAL'].astype(float)
supp_pivot_record = supp_record.pivot_table(
    index=['USUBJID', 'AESEQ'], columns='QNAM', values='QVAL', aggfunc='first'
).reset_index()
ae_enriched = ae.merge(supp_pivot_record, on=['USUBJID', 'AESEQ'], how='left')
```

**The 2024-2026 SUPP transition:** Therapeutic Area User Guides (TAUGs) increasingly use NS-- domain extensions or **Non-Standard Variables (NSV) Registry**-listed variables directly in the parent domain, instead of QNAM/QVAL pairs in SUPP--. Not a hard deprecation but the direction is clear. The Non-Standard Variables Registry at cdisc.org is the new canonical place to look up sponsor-extension variables.

## Date Handling -- The Partial-Date Reality

```python
dm['RFSTDT'] = pd.to_datetime(dm['RFSTDTC'], errors='coerce')
ae['AESTDT'] = pd.to_datetime(ae['AESTDTC'], errors='coerce')
ae['AEENDT'] = pd.to_datetime(ae['AEENDTC'], errors='coerce')

# Days from randomization to AE onset
ae_with_ref = ae.merge(dm[['USUBJID', 'RFSTDT']], on='USUBJID')
ae_with_ref['AE_ONSET_DAY'] = (ae_with_ref['AESTDT'] - ae_with_ref['RFSTDT']).dt.days
```

**SDTM dates are ISO 8601 strings.** Partial dates (e.g., '2023-03' without day) are common. `errors='coerce'` converts these to NaT rather than raising errors. For analysis requiring complete dates, CDISC conventions impute missing day as the 1st for start dates and the last day of the month for end dates, but imputation rules should match the SAP.

**SDTM records include EPOCH** (SCREENING, TREATMENT, FOLLOW-UP). For TEAEs, filter AEs to onset during or after the treatment epoch. Including pre-treatment AEs confounds the treatment effect estimate.

## Validation -- Pinnacle 21 and CORE

**Pinnacle 21 (Certara, formerly OpenCDISC)** is the de facto FDA submission validation standard. Validates against FDA Validation Rules + CDISC IG conformance + Define-XML schema. Severity tiers:

- **Reject** -- submission will not be accepted
- **Error** -- must justify
- **Warning** -- should investigate

**FDA Validation Rules** are published quarterly by FDA Office of Translational Sciences; Pinnacle 21 wraps these into its rule engine.

**CORE (CDISC Open Rules Engine, 2024)** is a newer open-source alternative using YAML-defined rules from the CDISC Rules Catalog. Gaining traction but not yet at Pinnacle-21 parity for confirmatory submissions.

```bash
# Pinnacle 21 Community (free; appropriate for non-pivotal trials)
p21-community validate --rules sdtmig-3.4 --output-dir validation_output study_data/
```

## Population Flags

| Flag | Source | Purpose |
|------|--------|---------|
| ITTFL | DM all randomised | Primary efficacy population (ICH E9 default) |
| FASFL | ITT minus eligibility failures + no post-baseline | Practical primary (FAS = Full Analysis Set) |
| SAFFL | EX (received at least one dose) | Safety analysis (AE reporting) |
| PPROTFL | SE + DS + protocol-violation list | Per-protocol; sensitivity only |
| EFFFL | Sponsor-defined | Modified ITT variants |

**FAS vs ITT subtlety:** FAS may exclude post-randomisation subjects (ineligibility, no post-baseline efficacy); ITT cannot. Many SAPs equate them; FDA may insist on stricter ITT at submission. Pre-specify both with explicit FAS exclusion criteria in the protocol.

## Missing Data Considerations -- The Clinical Reasoning Layer

Before any imputation/complete-case decision, **examine the DS (Disposition) domain to tabulate reasons for discontinuation by treatment arm.** If discontinuation rates or reasons differ between arms, missing data is likely informative (MNAR) and standard MMRM-MAR is questionable.

```python
ds = pd.read_csv('DS.csv')
dropouts = ds[ds['DSDECOD'] != 'COMPLETED']
dropouts_by_arm = dropouts.merge(adsl[['USUBJID', 'ARM']], on='USUBJID')
discontinuation_reasons = dropouts_by_arm.groupby(['ARM', 'DSDECOD']).size().unstack(fill_value=0)
```

This is the data-quality precursor to choosing the estimand strategy in trial-reporting (see ICH E9(R1)) — missing patterns informed by DS drive the choice between treatment-policy, hypothetical, or composite ICE strategies.

## Common Pitfalls

| Pitfall | Symptom | Solution |
|---------|---------|----------|
| Event-level merged onto subject-level without aggregation | Row count inflates after merge | Aggregate first, then merge |
| First chronological record used as baseline | Misclassified baseline | Filter on VSBLFL='Y' / LBBLFL='Y'; derive only if missing |
| Character (xxORRES) used for analysis | Inconsistent numeric coercion | Use xxSTRESN (numeric standardised); missing xxSTRESN with present xxORRES means 'NOT DONE' or '<LLOQ' |
| ARM used in safety analysis | Crossover or actual-treatment differs | Use ACTARM for safety; ARM for ITT efficacy |
| RFSTDTC used as TEAE reference | Includes screening AEs | Use RFXSTDTC (first treatment); cite ICH E2A |
| ADTTE CNSR confused with stat-pkg convention | Wrong event/censoring assignment | CDISC: CNSR=0 means event; convert: `event = (CNSR == 0).astype(int)` |
| Partial date parsing error | NaT in date column | `pd.to_datetime(..., errors='coerce')` |
| SUPPQUAL granularity confusion | Wrong rows merged | Check IDVAR before choosing subject vs record-level merge |
| Non-standard column names treated as standard SDTM | Missing variables | Inspect actual columns; map to semantic roles |
| Pinnacle 21 not run before submission | FDA reject | Always validate against current SDTMIG and FDA Validation Rules before submission |

### Common non-standard column mappings

| Standard SDTM | Common alternatives | Role |
|---------------|--------------------|------|
| ARM / ARMCD | TRTGRP, TRT01P, treatment, group | Treatment assignment |
| AEDECOD | AEPT, ae_term, preferred_term | AE preferred term |
| AESEV (text) | AESEV (numeric 1-4), severity, AETOXGR | Severity / toxicity grade |
| USUBJID | SUBJID, subject_id, patient_id | Subject identifier |
| RFXSTDTC | trt_start, first_dose_date | First treatment date |
| LBSTRESN | lab_value_num, result_numeric | Lab numeric result |

## Quantitative Thresholds and Conventions

| Threshold/Convention | Source | Rationale |
|----------------------|--------|-----------|
| TEAE window: 28-30 days post-treatment for small molecules | ICH E2A; sponsor convention | Mode of action, half-life inform window |
| RFXSTDTC for TEAE reference (not RFSTDTC) | ICH E2A | RFSTDTC includes screening; TEAE is post-treatment |
| ABLFL='Y' for last non-missing pre-dose | CDISC ADaM IG v1.3 | Standard baseline definition; trust SDTM flag |
| ARM (planned) for ITT efficacy; ACTARM for safety | ICH E9 | Crossover/PP-violation handling |
| Pinnacle 21 validation before submission | FDA Study Data Technical Conformance Guide | Standard quality gate; reject errors block acceptance |
| Define-XML 2.1 for studies starting >=March 15, 2023 | FDA Study Data Standards Catalog | Older 2.0 still accepted for prior studies |
| Dataset-JSON v1.1 (Dec 2025; FDA notice April 2025) | CDISC + FDA Federal Register | Modern replacement for XPT v5; timeline pending |
| CNSR=0 for events, positive integers for censoring | ADaM BDS-for-TTE v1.0 | OPPOSITE of R `survival` and most stat packages |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "RFXSTDTC or RFSTDTC for TEAE?" | RFXSTDTC per ICH E2A; RFSTDTC would include screening AEs |
| "Baseline from VSBLFL or derived?" | VSBLFL when present; documented derivation rule when missing |
| "Pinnacle 21 errors?" | All errors resolved or justified; warnings reviewed and documented |
| "Define-XML 2.1 ADRG provided?" | Yes — analysis dataset traceability documented to variable level |
| "ITT vs FAS reconciliation?" | Pre-specified in protocol with explicit FAS exclusion criteria |
| "OCCDS v1.1 multi-period flags?" | TRTEM01FL...TRTEM##FL pre-specified for crossover periods |
| "ADTTE CNSR convention conversion documented?" | Explicit: CDISC CNSR=0 means event; convert to event=1 for downstream R/Python |

## References

- CDISC. 2021. Analysis Data Model Implementation Guide (ADaMIG) v1.3.
- CDISC. 2021. Occurrence Data Structure (OCCDS) v1.1.
- CDISC. 2012/2024. ADaM Basic Data Structure for Time-to-Event Analyses v1.0.
- CDISC. 2025. Dataset-JSON v1.1.
- FDA. 2024. Study Data Technical Conformance Guide.
- FDA Federal Register Notice. April 2025. Dataset-JSON Pilot Comment Request.
- ICH. 1995. E2A: Clinical Safety Data Management -- Definitions and Standards for Expedited Reporting.
- ICH. 1998. E9: Statistical Principles for Clinical Trials.
- ICH. 2019. E9(R1) Addendum on Estimands and Sensitivity Analysis.
- Pinnacle 21 (Certara). 2024. Community Edition Validation Rules.
- PHUSE/CDISC. 2024. Dataset-JSON Pilot Reports.

## Related Skills

- clinical-biostatistics/logistic-regression - Model binary outcomes from prepared ADaM/SDTM data
- clinical-biostatistics/trial-reporting - Use prepared analysis datasets for ICH E9(R1) estimands and CONSORT 2025 reporting
- clinical-biostatistics/missing-data-sensitivity - DS-domain reasoning informs estimand choice
- clinical-biostatistics/survival-analysis - ADTTE CNSR convention; time-to-event preparation
- expression-matrix/metadata-joins - General metadata joining patterns
