# CDISC Data Handling - Usage Guide

## Overview

Reads and prepares CDISC SDTM clinical trial data for downstream statistical analysis. Covers loading .xpt and CSV domain files, joining subject-level demographics with event-level outcomes, aggregating adverse events and lab values per subject, and pivoting supplemental qualifier datasets.

## Prerequisites

```bash
pip install pyreadstat pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Load my clinical trial .xpt files and merge demographics with adverse events"
- "Create a subject-level analysis dataset from CDISC SDTM domains"
- "Aggregate adverse event severity per subject and join with treatment arms"
- "Pivot my SUPPQUAL dataset and merge it with the AE domain"

## Example Prompts

### Data Loading

> "Read the dm.xpt and ae.xpt files from my SDTM submission package and show me the column names and first few rows"

> "Load all domain CSVs in my data folder and tell me how many subjects are in each"

### Non-Standard Data

> "My clinical trial CSVs use TRTGRP instead of ARM and AEPT instead of AEDECOD. Help me map these to the standard roles and build the analysis dataset."

### Domain Joining

> "Merge demographics with adverse events so I have one row per subject with their treatment arm and whether they had a serious AE"

> "Build an analysis dataset combining DM demographics, baseline vital signs from VS, and baseline lab values from LB"

> "Join the disposition domain DS with DM to create a dataset showing which subjects completed the study"

### Specific Analyses

> "For each subject, count the number of adverse events, find the maximum severity, and flag whether they had any serious AEs"

> "Extract baseline systolic blood pressure and hemoglobin from VS and LB domains for my covariate-adjusted analysis"

> "Calculate the number of days from randomization to first adverse event onset for each subject"

### SUPPQUAL

> "Pivot my SUPPAE supplemental qualifiers and merge them with the AE domain"

## What the Agent Will Do

1. Load domain files (.xpt or CSV) using pyreadstat or pandas
2. Inspect column names and verify expected SDTM variables are present
3. Identify subject-level vs event-level domains
4. Aggregate event-level data to one row per subject using appropriate summary statistics
5. Merge aggregated summaries onto the DM backbone using USUBJID
6. Handle missing values from left joins (subjects with no events)
7. Report the final analysis dataset dimensions and a summary of key variables

## Tips

- Always check that USUBJID values are unique in your demographics (DM) domain before using it as the merge backbone. Duplicate USUBJIDs in DM indicate a data quality issue.
- Verify column names before applying SDTM conventions. Academic datasets distributed as CSV may use different naming (e.g., 'subject_id' instead of 'USUBJID').
- After every merge, check the resulting row count. If it increased unexpectedly, event-level data was not properly aggregated first.
- Use `how='left'` when merging onto DM to preserve all randomized subjects, even those with no events. Fill missing event indicators with 0 or False.
- For baseline values, always filter on the baseline flag (VSBLFL='Y', LBBLFL='Y') rather than selecting the first chronological record.
- When working with dates, expect partial dates (year-month only) and use `errors='coerce'` to avoid parsing failures.
- Keep the ARMCD (short code) for programmatic filtering and ARM (full label) for display purposes.

## Related Skills

- clinical-biostatistics/logistic-regression - Model binary outcomes from prepared analysis datasets
- expression-matrix/metadata-joins - General metadata joining patterns
