---
name: bio-metabolomics-msdial-preprocessing
description: Runs the MS-DIAL preprocessing workflow (peak picking, MS2Dec spectral deconvolution, alignment, gap-filling) and imports the alignment-result table into R or Python with honest filtering. Use when preprocessing LC-MS DDA/DIA (SWATH) raw data with MS-DIAL, deciding MS-DIAL vs XCMS, configuring the MsdialConsoleApp console run, or parsing an MS-DIAL export into a clean feature matrix. For programmatic R peak detection and the feature-table-as-artifact framing see metabolomics/xcms-preprocessing; for lipid annotation mode see metabolomics/lipidomics; for MSI-level confidence honesty see metabolomics/metabolite-annotation; for drift correction and QC see metabolomics/normalization-qc.
tool_type: mixed
primary_tool: msdial
---

## Version Compatibility

Reference examples tested with: MS-DIAL 5.x (LC-MS) / MS-DIAL 4.x (GC-MS), pandas 2.2+, R 4.3+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: run `MsdialConsoleApp` with no arguments to print the current subcommand/flag list
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show pandas` then `help(module.function)` to check signatures

The MS-DIAL GUI runs only on Windows; the console (MsdialConsoleApp) is the cross-platform headless entry. Which build supports a task is itself a constraint: MS-DIAL 5-alpha covers DI-MS, IM-MS, LC-MS, LC-IM-MS but NOT GC-MS - GC-EI stays in the MS-DIAL 4 lineage. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt rather than retrying.

# MS-DIAL Preprocessing

**"Process my LC-MS run with MS-DIAL and give me a feature table"** -> Pick peaks per file, deconvolve chimeric MS/MS into clean component spectra (MS2Dec), align across samples, gap-fill, then import the alignment result and filter it honestly.
- CLI: `MsdialConsoleApp lcmsdda|lcmsdia|gcms -i <in> -o <out> -m <param.txt>`
- R: `read.csv(..., skip = 4, check.names = FALSE)` to parse the alignment export
- Python: `pandas.read_csv(..., skiprows=4)` for the same export

## The Single Most Important Insight -- Preprocessing Software Is Not Neutral

The same raw files through MS-DIAL versus XCMS yield different feature tables and different marker lists. Li 2018 benchmarked five tools on a 1,100-compound standard and found that while feature *detection* was broadly similar, *quantification* and the set of selected discriminating markers differed by tool. A metabolomics "hit" is conditional on (raw data + software + version + every parameter + fill/filter order), not on the raw files alone. MS-DIAL's specific differentiator is **MS2Dec deconvolution**: it reconstructs clean, library-matchable MS/MS spectra from chimeric DDA/DIA fragment data, which is what makes wide-window DIA (SWATH) tractable at all. Report the full processing specification as part of the result, and treat a finding that survives only one pipeline as a candidate, not a result.

## MS-DIAL vs XCMS

| Axis | MS-DIAL | XCMS |
|---|---|---|
| Interface | Windows GUI + cross-platform console | R package (scriptable everywhere) |
| Core differentiator | MS2Dec MS/MS deconvolution (DDA + DIA) | centWave peak picking, full programmatic control |
| Annotation | Built-in (library + MS-FINDER + LipidBlast) | Separate (CAMERA, downstream tools) |
| Lipidomics | Strong (predicted-CCS / EAD structural elucidation in v5) | Manual |
| Reproducibility unit | Param file + GUI choices | Versioned R script |
| Best when | DIA data, lipidomics, GUI workflow, built-in IDs | Scripted pipelines, custom parameters, cohort scale |

Use MS-DIAL when DIA deconvolution or built-in lipid annotation is the point; use metabolomics/xcms-preprocessing for fully scripted, version-pinned cohort processing. The strongest untargeted claims replicate across both.

## Decision Tree by Scenario

| Situation | Do | Why |
|---|---|---|
| LC-MS, top-N MS/MS (DDA) | `lcmsdda` console / GUI LC-MS DDA | Cleaner per-precursor MS2, but intensity-biased, stochastic coverage |
| LC-MS, wide-window MS/MS (DIA / SWATH) | `lcmsdia` (ABF input only) | Complete MS2 coverage; chimeric spectra REQUIRE MS2Dec to be usable |
| GC-EI run | `gcms` (MS-DIAL 4 build), or AMDIS/eRah | EI fragments every co-eluting compound; deconvolution IS detection (see below) |
| Headless / Linux cluster | MsdialConsoleApp with a `-m` param file | GUI is Windows-only; console is the reproducible batch path |
| Lipid-focused study | MS-DIAL + LipidBlast | -> metabolomics/lipidomics for lipid annotation mode |
| Already have an alignment CSV | skip processing, parse + filter | See import + honest-filter sections below |

## Why GC-EI Is Different (and stays in MS-DIAL 4)

In GC-EI, 70 eV ionization fragments every compound reproducibly, so the trace at any retention time is a superposition of fragments from several co-eluting molecules. Naive peak picking conflates them; **deconvolution into component spectra IS the feature-detection step**, then each component is matched against EI+RI libraries (NIST, FiehnLib). Cross-run/cross-lab alignment uses **retention index** (Kovats n-alkanes, or Fiehn FAME markers giving diagnostic m/z 74/87) rather than raw RT, because RT drifts with column aging. MS-DIAL 5-alpha explicitly excludes GC-MS; use the `gcms` token in a MS-DIAL 4 build, or AMDIS/eRah, for GC-EI work.

## Run MS-DIAL Headless (console)

**Goal:** Process a folder of converted spectra into an alignment table without the GUI.

**Approach:** Pick the analysis-type token, point `-i`/`-o`/`-m` at input dir, output dir, and a method (parameter) file; keep `-p` only if the project should reopen in the GUI.

```bash
# DDA LC-MS: accepts netCDF/mzML/ABF. Output is *.msdial in the output dir.
MsdialConsoleApp lcmsdda -i ./LCMS_DDA/ -o ./LCMS_DDA_out/ -m ./Msdial-lcms-dda-Param.txt

# DIA/SWATH LC-MS: accepts ABF ONLY (convert vendor raw -> ABF first). MS2Dec is the point.
MsdialConsoleApp lcmsdia -i ./LCMS_DIA/ -o ./LCMS_DIA_out/ -m ./Msdial-lcms-dia-Param.txt

# GC-EI (MS-DIAL 4 build): retention-index alignment, quant-mass quantification.
MsdialConsoleApp gcms -i ./GCMS/ -o ./GCMS_out/ -m ./Msdial-GCMS-Param.txt -p
```

The parameter file is plain text (one `Key=Value` per line). The `Minimum peak height` key is the direct analog of an intensity floor and is instrument-dependent: the GUI default is tuned for a TOF and is often far too high (or its baseline assumption wrong) for an Orbitrap. Set the alignment reference to a pooled QC, never to file #1 by default.

## Import the Alignment Result into R

**Goal:** Split the MS-DIAL alignment export into a feature-metadata frame and an intensity matrix.

**Approach:** The export carries four header rows above the real column header (sample class / file type / injection order / batch), so skip them; metadata columns precede the per-sample Area columns.

```r
# MS-DIAL alignment export: real column header is on row 5, so skip the first 4 rows.
msdial <- read.csv('AlignResult.txt', sep = '\t', skip = 4, check.names = FALSE)

# Metadata columns appear before the per-sample intensity columns. Common ones:
# 'Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name', 'Adduct type',
# 'Fill %', 'MS/MS assigned', 'Reference RT', 'Formula', 'Ontology', 'INCHIKEY',
# 'SMILES', 'Annotation tag (VS1.0)'. Sample columns are everything after these.
meta_cols <- c('Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name',
               'Adduct type', 'Fill %', 'MS/MS assigned', 'Annotation tag (VS1.0)')
meta_cols <- intersect(meta_cols, colnames(msdial))
sample_cols <- setdiff(colnames(msdial), colnames(msdial)[seq_len(max(match(meta_cols, colnames(msdial))))])

feature_info <- msdial[, meta_cols]
intensity <- as.matrix(msdial[, sample_cols])
rownames(intensity) <- msdial[['Alignment ID']]
```

## Import the Alignment Result into Python

**Goal:** Same split, in pandas.

**Approach:** `skiprows=4` to land on the real header; slice metadata vs sample columns by position after the last known metadata column.

```python
import pandas as pd

msdial = pd.read_csv('AlignResult.txt', sep='\t', skiprows=4)
meta_cols = ['Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name', 'Adduct type', 'Fill %', 'MS/MS assigned', 'Annotation tag (VS1.0)']
meta_cols = [c for c in meta_cols if c in msdial.columns]
last_meta = max(msdial.columns.get_loc(c) for c in meta_cols)
sample_cols = msdial.columns[last_meta + 1:]

feature_info = msdial[meta_cols].copy()
intensity = msdial[sample_cols].set_axis(msdial['Alignment ID']) if False else msdial[sample_cols].copy()
intensity.index = msdial['Alignment ID']
```

## Filter the Table Honestly

**Goal:** Keep features supported by real signal and known confidence, without overtrusting annotation tags.

**Approach:** Filter on Fill% (cross-sample presence), require MS/MS support for any feature called identified, and tie the annotation tag to a real MSI confidence level rather than treating a name as proof.

```r
# Fill% is the fraction of samples with a DETECTED (not gap-filled) peak. Low Fill% means
# the feature exists mostly as gap-filled noise-floor integrals, which fabricate intensity
# (an honest 'below detection' becomes a positive number). 70% is a common floor.
keep_fill <- feature_info[['Fill %']] >= 70

# An annotated name without MS/MS is at best an MSI Level 2/3 putative ID (accurate mass
# only). Require 'MS/MS assigned == TRUE' before trusting any identity downstream.
has_msms <- feature_info[['MS/MS assigned']] == 'TRUE'

# Annotation tag confidence (do NOT treat a name as an identification). The exact tag
# vocabulary is MS-DIAL-version-dependent, so inspect unique(feature_info[['Annotation tag (VS1.0)']])
# and map the strings the build actually emits rather than hard-coding them:
#   Metabolite / Lipid  with MS/MS  -> MSI Level 2 (spectral library match)
#   Suggested*          mass-only   -> MSI Level 3 (putative, no MS/MS)
#   Unknown                         -> unannotated feature
feature_info$msi_level <- ifelse(feature_info[['Annotation tag (VS1.0)']] %in% c('Metabolite', 'Lipid') & has_msms, 2,
                          ifelse(grepl('^Suggested', feature_info[['Annotation tag (VS1.0)']]), 3, NA))

filtered <- intensity[keep_fill, ]
```

Confidence-level honesty and orthogonal-evidence identification belong to metabolomics/metabolite-annotation; this skill only routes the tag to the right level. Fill% / blank / drift filtering interacts with normalization-qc - process blanks and pooled QCs through the SAME run, then filter the aligned table.

## Per-Method Failure Modes

### DIA processed as DDA (wrong console token)
- **Trigger:** Running SWATH/DIA data through `lcmsdda`.
- **Mechanism:** `lcmsdda` does not deconvolve wide-isolation chimeric MS/MS, so fragments from co-isolated precursors stay mixed.
- **Symptom:** Library matches to the wrong compound; "clean-looking" spectra that fail orthogonal confirmation.
- **Fix:** Use `lcmsdia` (ABF input only); MS2Dec deconvolution is the entire reason to run DIA in MS-DIAL.

### Over-trusting the annotation tag
- **Trigger:** Filtering on `Annotation tag != Unknown` and calling the survivors "identified."
- **Mechanism:** A `Suggested*` tag is an accurate-mass guess with no MS/MS; a named hit without MS/MS is MSI Level 3.
- **Symptom:** A marker list full of confident-sounding names that do not validate against standards.
- **Fix:** Require `MS/MS assigned == TRUE` for any identity claim; map tags to MSI levels (see filtering section) and defer to metabolomics/metabolite-annotation.

### Gap-fill masquerading as measurement
- **Trigger:** Treating low-Fill% features as quantitative.
- **Mechanism:** Gap-filling integrates whatever signal sits in the m/z-RT box even when no peak exists, turning a true below-detection (MNAR, left-censored) value into a positive number.
- **Symptom:** "Significant" features that are mostly gap-filled in one group; shrunken fold-changes for on/off markers.
- **Fix:** Report per-feature filled fraction; gate on Fill%; for inferential stats prefer MNAR-aware imputation over naive fill (see metabolomics/normalization-qc).

### GC-EI run through an LC pipeline / MS-DIAL 5
- **Trigger:** Sending GC-EI data to MS-DIAL 5-alpha or treating it like LC peak-pick-then-group.
- **Mechanism:** MS-DIAL 5-alpha excludes GC-MS; EI needs component deconvolution, not adduct-style peak picking, and RI (not RT) alignment.
- **Symptom:** No GC mode available; or conflated co-eluting compounds and cross-lab RT misalignment.
- **Fix:** Use the `gcms` token in a MS-DIAL 4 build (or AMDIS/eRah); align on Kovats/FAME retention index.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Fill% >= 70% | Common untargeted practice | Below this, the feature is mostly gap-filled noise-floor integrals, not measurements |
| QC CV (RSD) < 20-30% | Broadhurst 2018 | Technical reproducibility floor; drop features noisier than this in pooled QCs |
| D-ratio (sd_QC/sd_sample) < 0.5 | Broadhurst 2018 | Keeps features whose technical variance is well below biological variance |
| Blank filter: sample mean > 3-5x blank mean | Broadhurst 2018 | Removes background/contaminant features present in process blanks |
| ~10x more features than compounds | Mahieu 2017 | One metabolite makes adducts/isotopes/fragments; counting features over-counts hypotheses |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| All columns land in one field on import | Header offset wrong; tab-separated export read as CSV | `skip=4` (R) / `skiprows=4` (Python), set `sep='\t'` |
| `lcmsdia` rejects mzML input | DIA mode accepts ABF only | Convert vendor raw to ABF (Reifycs ABF converter) before `lcmsdia` |
| `Annotation tag` column not found | Header changes across versions (e.g. `Annotation tag (VS1.0)`) | Match by prefix / inspect `colnames()`; do not hard-code the suffix |
| No GC-MS option in MS-DIAL 5 | 5-alpha excludes GC-MS | Use a MS-DIAL 4 build's `gcms` token, or AMDIS/eRah |
| Console command not found on Linux | Expecting the GUI executable | The GUI is Windows-only; run `MsdialConsoleApp` (cross-platform) |
| Few features detected | `Minimum peak height` default too high for the instrument | Lower it toward the real baseline; defaults are TOF-tuned |

## References

- Tsugawa H, Cajka T, Kind T, Ma Y, Higgins B, Ikeda K, Kanazawa M, VanderGheynst J, Fiehn O, Arita M. MS-DIAL: data-independent MS/MS deconvolution for comprehensive metabolome analysis. *Nat Methods.* 2015; 12(6):523-526.
- Tsugawa H, Ikeda K, Takahashi M, et al. A lipidome atlas in MS-DIAL 4. *Nat Biotechnol.* 2020; 38(10):1159-1163.
- Takeda H, Takahashi M, Ikeda K, et al. MS-DIAL 5 multimodal mass spectrometry data mining unveils lipidome complexities. *Nat Commun.* 2024; 15:9903.
- Li Z, Lu Y, Guo Y, Cao H, Wang Q, Shui W. Comprehensive evaluation of untargeted metabolomics data processing software in feature detection, quantification and discriminating marker selection. *Anal Chim Acta.* 2018; 1029:50-57.
- Mahieu NG, Patti GJ. Systems-level annotation of a metabolomics data set reduces 25,000 features to fewer than 1,000 unique metabolites. *Anal Chem.* 2017; 89(19):10397-10406.
- Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR, Dunn WB. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics.* 2018; 14(6):72.
- Stein SE. An integrated method for spectrum extraction and compound identification from gas chromatography/mass spectrometry data (AMDIS). *J Am Soc Mass Spectrom.* 1999; 10(8):770-781.

## Related Skills

- metabolomics/xcms-preprocessing - Programmatic R preprocessing and the feature-table-as-artifact framing
- metabolomics/lipidomics - Lipid annotation mode and LipidBlast workflows
- metabolomics/metabolite-annotation - MSI confidence levels and orthogonal-evidence identification
- metabolomics/normalization-qc - Drift correction, QC/CV/D-ratio filtering, MNAR-aware imputation
