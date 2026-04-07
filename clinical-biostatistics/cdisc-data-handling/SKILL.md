---
name: bio-clinical-biostatistics-cdisc-data
description: Reads and prepares CDISC SDTM clinical trial data for analysis. Handles domain tables (DM, AE, EX, VS, LB), USUBJID-based joins, event-to-subject aggregation, and SUPPQUAL pivoting. Use when working with clinical trial datasets in CDISC/SDTM format or .xpt files.
tool_type: python
primary_tool: pyreadstat
---

## Version Compatibility

Reference examples tested with: pyreadstat 1.2+, pandas 2.1+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# CDISC SDTM Data Handling

**"Load clinical trial data"** -> Parse CDISC SDTM domain files and prepare subject-level analysis datasets by joining and aggregating across domains.
- Python: `pyreadstat.read_xport()`, `pd.read_sas()`, `pd.merge()`

## Domain Overview

| Domain | Level | Description | Key Variables |
|--------|-------|-------------|---------------|
| DM | Subject | Demographics (one row per subject) | USUBJID, ARM, ARMCD, AGE, SEX, RACE, RFSTDTC |
| AE | Event | Adverse events (multiple per subject) | USUBJID, AETERM, AEDECOD, AEBODSYS, AESEV, AESER |
| EX | Event | Drug exposure/dosing | USUBJID, EXTRT, EXDOSE, EXSTDTC, EXENDTC |
| VS | Event | Vital signs | USUBJID, VSTESTCD, VSSTRESN, VSBLFL |
| LB | Event | Lab results | USUBJID, LBTESTCD, LBSTRESN, LBBLFL |
| DS | Event | Disposition | USUBJID, DSDECOD, DSSTDTC |

SDTM is the raw tabulation standard. Many organizations also distribute ADaM (Analysis Data Model) datasets, which are analysis-ready. The key ADaM dataset is ADSL (subject-level), which already contains many aggregations this skill constructs manually. If ADSL is available, check whether the needed derived variables already exist before rebuilding them from SDTM domains.

USUBJID (typically STUDYID-SITEID-SUBJID) is the universal merge key across all domains. Subject-level domains (DM) have one row per USUBJID; event-level domains (AE, EX, VS, LB, DS) have multiple rows per subject.

## Reading .xpt Files

**Goal:** Load CDISC SDTM domain data from SAS transport (.xpt) files into pandas DataFrames.

**Approach:** Use pyreadstat for full metadata support, with pandas or CSV fallbacks for simpler cases.

```python
import pyreadstat
import pandas as pd

# pyreadstat (recommended - maintained by Roche, handles metadata)
dm, meta = pyreadstat.read_xport('dm.xpt')
# meta.column_names, meta.column_labels, meta.variable_value_labels

# pandas built-in (SAS XPORT v5)
dm = pd.read_sas('dm.xpt', format='xport', encoding='utf-8')

# CSV fallback (common in academic datasets)
dm = pd.read_csv('DM.csv')
```

When pyreadstat is available, the metadata object provides column labels, value labels, and format information that are lost with other readers.

## Joining Domains

**Goal:** Combine treatment, demographic, and outcome information into a single subject-level analysis dataset.

**Approach:** Merge event-level data back to DM using USUBJID, aggregating events to one-row-per-subject before merging.

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

### Additional Aggregation Patterns

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

### Choosing the Aggregation Strategy

Different aggregation strategies answer fundamentally different clinical questions:

| Strategy | Scientific question | Example |
|----------|-------------------|---------|
| Any event (binary) | Does treatment increase the probability of experiencing the event at all? | Had any serious AE: yes/no |
| Event count | Does treatment increase the burden of events per patient? | Total AE count per subject |
| Maximum severity | Does treatment shift patients toward more severe manifestations? | Worst AESEV per subject |
| First event + time | Does treatment delay onset of the event? | Time to first serious AE |

These are not interchangeable. A drug might not change the proportion of patients with AEs (binary: no effect) but increase the number of events per patient (count: harmful). The choice must be driven by the scientific question in the statistical analysis plan, not by analytic convenience.

### Multi-Domain Merge

**Goal:** Build an analysis dataset that combines demographics, adverse events, and baseline lab values.

**Approach:** Aggregate each event-level domain independently, then merge all summaries onto the DM backbone.

```python
ae_summary = ae.groupby('USUBJID').agg(
    ae_count=('AETERM', 'count'),
    had_serious=('AESER', lambda x: (x == 'Y').any()),
    max_severity=('AESEV_NUM', 'max')
).reset_index()

lb_baseline = lb[lb['LBBLFL'] == 'Y'].pivot_table(
    index='USUBJID', columns='LBTESTCD', values='LBSTRESN', aggfunc='first'
).reset_index()

analysis = dm.merge(ae_summary, on='USUBJID', how='left')
analysis = analysis.merge(lb_baseline, on='USUBJID', how='left')
analysis['ae_count'] = analysis['ae_count'].fillna(0).astype(int)
analysis['had_serious'] = analysis['had_serious'].fillna(False)
```

## SUPPQUAL Pivoting

**Goal:** Enrich domain data with supplemental qualifier variables stored in SUPP-- datasets.

**Approach:** Pivot SUPPQUAL long format (QNAM/QVAL pairs) to wide format and merge back to the parent domain.

```python
supp = pd.read_sas('suppae.xpt', format='xport', encoding='utf-8')
supp_pivot = supp.pivot_table(
    index='USUBJID', columns='QNAM', values='QVAL', aggfunc='first'
).reset_index()
ae_enriched = ae.merge(supp_pivot, on='USUBJID', how='left')
```

For record-level SUPPQUAL data (where IDVAR and IDVARVAL identify specific rows), merge on both USUBJID and the record identifier:

```python
supp_record = supp[supp['IDVAR'] == 'AESEQ'].copy()
supp_record['AESEQ'] = supp_record['IDVARVAL'].astype(float)
supp_pivot_record = supp_record.pivot_table(
    index=['USUBJID', 'AESEQ'], columns='QNAM', values='QVAL', aggfunc='first'
).reset_index()
ae_enriched = ae.merge(supp_pivot_record, on=['USUBJID', 'AESEQ'], how='left')
```

## Baseline Values

**Goal:** Extract baseline measurements from event-level domains for use as covariates.

**Approach:** Filter on the baseline flag (xxBLFL='Y') and pivot test codes to columns.

```python
vs_baseline = vs[vs['VSBLFL'] == 'Y'].pivot_table(
    index='USUBJID', columns='VSTESTCD', values='VSSTRESN', aggfunc='first'
).reset_index()
# Columns: USUBJID, SYSBP, DIABP, PULSE, TEMP, etc.

lb_baseline = lb[lb['LBBLFL'] == 'Y'].pivot_table(
    index='USUBJID', columns='LBTESTCD', values='LBSTRESN', aggfunc='first'
).reset_index()
# Columns: USUBJID, ALT, AST, CREAT, HGB, etc.
```

## Date Handling

**Goal:** Convert SDTM ISO 8601 date strings to datetime objects for time-based calculations.

**Approach:** Parse with `pd.to_datetime` using `errors='coerce'` to handle partial dates gracefully.

```python
dm['RFSTDT'] = pd.to_datetime(dm['RFSTDTC'], errors='coerce')
ae['AESTDT'] = pd.to_datetime(ae['AESTDTC'], errors='coerce')
ae['AEENDT'] = pd.to_datetime(ae['AEENDTC'], errors='coerce')

# Days from randomization to AE onset
ae_with_ref = ae.merge(dm[['USUBJID', 'RFSTDT']], on='USUBJID')
ae_with_ref['AE_ONSET_DAY'] = (ae_with_ref['AESTDT'] - ae_with_ref['RFSTDT']).dt.days
```

Partial dates (e.g., '2023-03' without day) are common in SDTM. `errors='coerce'` converts these to NaT rather than raising errors. For analysis requiring complete dates, CDISC conventions impute missing day as the 1st for start dates and the last day of the month for end dates, but imputation rules should match the SAP.

SDTM records include EPOCH (SCREENING, TREATMENT, FOLLOW-UP). For treatment-emergent adverse events (TEAEs), filter AEs to those with onset during or after the treatment epoch. Including pre-treatment AEs in a treatment effect analysis confounds the estimate. Use RFXSTDTC (first treatment date) rather than RFSTDTC (first study activity, which may be screening) as the reference for TEAE determination.

## Missing Data Considerations in CDISC

Before choosing an imputation or complete-case strategy, consider why data is missing. In clinical trials, the reason for missingness often has clinical meaning: a missing lab value because the patient discontinued due to adverse events is fundamentally different from a missing value because a blood draw was missed at a visit. The former is informative (MNAR) and related to outcome; the latter may be ignorable (MAR). Examining the DS (Disposition) domain and linking discontinuation reasons to missing observations in AE/LB/VS is essential before any statistical handling decision.

## Common Pitfalls

- **Subject-level vs event-level confusion:** DM has one row per USUBJID. AE/LB/VS have multiple. Always aggregate events before subject-level merges to avoid inflating denominators.
- **Baseline flag (xxBLFL):** Filter on VSBLFL='Y' or LBBLFL='Y' for baseline values. Do not assume the first chronological record is baseline.
- **Character vs numeric results:** SDTM stores both xxORRES (character) and xxSTRESN (float). Always use xxSTRESN for analysis. Missing xxSTRESN with present xxORRES typically means 'NOT DONE' or '<LLOQ'.
- **Treatment coding:** ARMCD is the short code, ARM is the label. Use ARMCD for programmatic comparisons. In crossover designs, actual treatment at a visit may differ from ARM.
- **Date handling:** SDTM dates are ISO 8601 strings. Partial dates like '2023-03' are common. Parse with `pd.to_datetime(col, errors='coerce')`.
- **Non-standard column names:** Academic and CSV datasets often deviate from SDTM naming. Always inspect actual columns and map to their semantic roles before analysis. Common deviations:

| Standard SDTM | Common alternatives | Role |
|---------------|-------------------|------|
| ARM / ARMCD | TRTGRP, TRT01P, treatment, group | Treatment assignment |
| AEDECOD | AEPT, ae_term, preferred_term | AE preferred term |
| AESEV (text) | AESEV (numeric 1-4), severity, AETOXGR | Severity / toxicity grade |
| USUBJID | SUBJID, subject_id, patient_id | Subject identifier |
| RFXSTDTC | trt_start, first_dose_date | First treatment date |
- **SUPPQUAL granularity:** Some SUPPQUAL records are subject-level (IDVAR is empty), others are record-level (IDVAR='AESEQ'). Check IDVAR before choosing the merge strategy.

## Related Skills

- clinical-biostatistics/logistic-regression - Model binary outcomes from prepared analysis datasets
- expression-matrix/metadata-joins - General metadata joining patterns
