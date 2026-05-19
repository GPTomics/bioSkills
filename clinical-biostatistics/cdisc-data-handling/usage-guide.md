# CDISC Data Handling - Usage Guide

## Overview

Reads, validates, and prepares CDISC SDTM and ADaM clinical trial data for downstream statistical analysis. Covers SDTM domain joins (DM, AE, EX, VS, LB, DS), ADaM architecture (ADSL, BDS, OCCDS, ADTTE) with traceability and Define-XML 2.1, treatment-emergent AE conventions per ICH E2A, baseline derivation, SUPPQUAL pivoting and NSV Registry transition, Pinnacle 21 / CORE validation, and the SAS XPT v5 → Dataset-JSON v1.1 transition (FDA Federal Register notice April 2025).

## Prerequisites

```bash
pip install pyreadstat pandas numpy
```

R is recommended for production ADaM derivation (Roche/openpharma):

```r
install.packages(c('haven', 'admiral', 'metacore', 'metatools', 'xportr'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Load my CDISC .xpt files via pyreadstat with full metadata"
- "Build a subject-level analysis dataset by joining DM, AE, EX, VS, LB"
- "Derive ADaM BDS for labs (ADLB) with baseline flag and change-from-baseline"
- "Convert ADTTE CNSR=0/1 convention to standard event indicator for R survival"
- "Validate my submission package against Pinnacle 21 Community + FDA Validation Rules"

## Example Prompts

### Reading and inspecting

> "Read all .xpt files in my SDTM submission package. Show me column labels, value labels, and basic dimensions per domain."

> "Inspect non-standard column names in my academic CSV dataset (e.g., 'TRTGRP' instead of 'ARM') and map to standard SDTM roles."

### SDTM merging

> "Merge DM with AE aggregated by USUBJID to create subject-level dataset with 'had any serious AE' and 'max AE severity' columns."

> "Multi-domain merge: DM + AE summary + baseline VS + baseline LB into single analysis-ready DataFrame."

### ADaM derivation

> "Derive ADaM BDS for labs from LB: PARAM/PARAMCD/AVAL/AVALC/BASE/CHG/PCHG/ABLFL/EPOCH. Pre-compute baseline from LBBLFL='Y' records."

> "Derive ADAE (OCCDS) from AE with TRTEMFL='Y' if AESTDTC >= TRTSDT AND AESTDTC <= TRTEDT + 28 days."

> "Derive ADTTE for OS: PARAMCD='OS', STARTDT=TRTSDT, ADT = death date or last alive contact, CNSR=0 for event/1+ for censoring, AVAL = ADT - STARTDT in days."

### Treatment-emergent AE

> "Implement TEAE flag with X=30 days post-treatment for small molecule. Use RFXSTDTC (first treatment) NOT RFSTDTC (screening) per ICH E2A."

> "For my biologic with extended half-life, set X=90 days post-treatment for TEAE window."

### Baseline derivation

> "Derive baseline from LB when LBBLFL='Y' is present; fall back to latest non-missing pre-dose value when flag missing."

> "Multi-period crossover trial: use BASETYPE to distinguish 'Period 1 Baseline' vs 'Period 2 Baseline'."

### SUPPQUAL

> "Pivot SUPPAE (subject-level QNAM/QVAL pairs) and merge with AE."

> "Pivot SUPPAE record-level (IDVAR='AESEQ'); merge on both USUBJID AND AESEQ."

### Validation

> "Run Pinnacle 21 Community against my SDTM submission. Resolve all Reject errors; justify Warnings in submission appendix."

### Dataset-JSON

> "Convert my submission from SAS XPT v5 to Dataset-JSON v1.1 (CDISC Dec 2025) using the PHUSE-CDISC pilot tooling. FDA Federal Register notice April 2025."

## What the Agent Will Do

1. Load SDTM and/or ADaM domain files via pyreadstat (preferred) or pandas
2. Inspect actual column names; map non-standard to SDTM roles if needed
3. Aggregate event-level data to subject level BEFORE merging onto DM (avoid row inflation)
4. Apply CDISC conventions: ARM (planned) vs ACTARM (actual); RFXSTDTC (first treatment) vs RFSTDTC (first activity)
5. For ADaM, derive PARAM/PARAMCD/AVAL/BASE/CHG/PCHG structure with proper flags
6. Convert ADTTE CNSR convention (CNSR=0 is event) for downstream R/Python (event=1 is event)
7. Tabulate DS (Disposition) by arm to inform missing-data strategy in trial-reporting
8. Validate via Pinnacle 21 against current SDTMIG and FDA Validation Rules

## Tips

- **ADaM dataset types:** ADSL (subject; one row per USUBJID), BDS (parameter-timepoint; long), OCCDS (occurrence; one row per AE/CM/MH), ADTTE (TTE BDS variant).
- **The ADaM Fundamental Principles** require analysis-ready, traceable, Define-XML-documented datasets. "One PROC away from the CSR table."
- **ARM vs ACTARM:** ITT efficacy uses ARM (planned at randomisation); safety uses ACTARM (actual treatment received).
- **RFXSTDTC vs RFSTDTC:** for TEAE calculations, ALWAYS use RFXSTDTC (first treatment) per ICH E2A. RFSTDTC includes screening AEs.
- **ADTTE CNSR convention trap:** CDISC ADaM uses CNSR=0 for events (positive integers for censoring reasons). Most statistical packages use 1=event. ALWAYS convert: `event = (CNSR == 0).astype(int)`.
- **Always check expected counts in merges.** Event-level domains have multiple rows per subject. Aggregate FIRST, then merge to DM.
- **Use `how='left'` for merging onto DM** to preserve all randomised subjects, even those with no events.
- **Filter on VSBLFL='Y' / LBBLFL='Y' for baseline.** Trust the SDTM flag when present. Only derive when missing.
- **Use xxSTRESN (numeric standardised) for analysis**, NOT xxORRES (character). Missing xxSTRESN with present xxORRES means 'NOT DONE' or '<LLOQ'.
- **Partial dates ('2023-03' without day) are common in SDTM.** Parse with `pd.to_datetime(col, errors='coerce')`.
- **DS (Disposition) domain is the gateway** to missing-data strategy. Tabulate dropout reasons by arm; if differential, MAR is suspect.
- **OCCDS v1.1 added TRTEM##FL** multi-period flags for crossover/multi-phase studies. Use these instead of single TRTEMFL.
- **Define-XML 2.1 is FDA-recommended** for studies starting on/after March 15, 2023. Older 2.0 still accepted for prior studies.
- **Pinnacle 21 Community is free** and standard for non-pivotal trials; Enterprise edition for sponsor-pivotal. All errors must be resolved or justified.
- **CORE (CDISC Open Rules Engine, 2024)** is the modern open-source alternative to Pinnacle 21; uses YAML rules from CDISC Rules Catalog. Not yet at Pinnacle parity for confirmatory.
- **Dataset-JSON v1.1 (Dec 2025)** is the modern replacement for SAS XPT v5 (which dates to 1995 with 8-char varname limits). FDA Federal Register notice April 2025; adoption timeline pending.
- **ADaMIG v3.0 is in development** (originally 2025, slipping to 2026+); will consolidate ADaMIG + OCCDS + BDS-TTE + ADAE supplement into single unified IG.
- **SUPPQUAL is being de-emphasised in favor of NSV Registry**-listed variables in parent domain. Therapeutic Area User Guides (TAUGs) increasingly use NS-- domains or NSV directly.

## Related Skills

- clinical-biostatistics/logistic-regression - Model from prepared ADaM/SDTM data
- clinical-biostatistics/trial-reporting - Use prepared data for ICH E9(R1) estimands + CONSORT 2025
- clinical-biostatistics/missing-data-sensitivity - DS reasoning informs estimand choice
- clinical-biostatistics/survival-analysis - ADTTE CNSR convention; time-to-event preparation
- expression-matrix/metadata-joins - General metadata joining patterns
