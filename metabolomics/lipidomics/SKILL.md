---
name: bio-metabolomics-lipidomics
description: Assigns honest lipid annotation levels, designs class-based internal-standard quantification, and runs lipid-aware differential and enrichment analysis with lipidr, guarding against in-source-fragment phantoms, sn-position over-claims, and invalid cross-class quantification. Use when naming or canonicalizing lipid species (shorthand separators, Goslin), deciding shotgun vs RP vs HILIC LC-MS, picking internal standards (SPLASH/EquiSPLASH), interpreting MS-DIAL/LipidSearch output, or comparing lipid classes. For general feature detection see metabolomics/xcms-preprocessing and metabolomics/msdial-preprocessing; for non-lipid annotation confidence see metabolomics/metabolite-annotation; for normalization/QC see metabolomics/normalization-qc; for multivariate stats see metabolomics/statistical-analysis.
tool_type: r
primary_tool: lipidr
---

## Version Compatibility

Reference examples tested with: lipidr 2.16+, pygoslin 2.0+, MS-DIAL 5+

The achievable annotation level is fixed by the acquired evidence, not the software: sn-position and double-bond localization require EAD/OzID/PB/UVPD data that routine CID never produces, and class-resolved quantification requires one isotope-labeled internal standard per class. Verify both before trusting a name or a number.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('lipidr')` then `?function_name` to verify parameters
- Python: `pip show pygoslin` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Lipidomics Analysis

**"Analyze my lipidomics data"** -> Canonicalize names to the resolution level the evidence supports, quantify each class against its own standard, then run class/chain-aware differential and enrichment analysis.
- R: `lipidr::read_skyline()` / `as_lipidomics_experiment()`, `de_analysis()`, `lsea()`
- Nomenclature: `pygoslin` (Python) or `rgoslin` (R) for parsing/canonicalization
- Identification: MS-DIAL 5 (open) or LipidSearch (commercial) upstream

## The Single Most Important Insight -- A Lipid Name Is a Structural-Resolution Claim the Software Usually Overstates

The Liebisch/LIPID MAPS shorthand encodes, in its punctuation, exactly how much structure was measured: `PC 34:1` (space, sum composition) < `PC 16:0_18:1` (underscore, chains known) < `PC 16:0/18:1` (slash, sn-resolved) < `PC 16:0/18:1(9Z)` (double-bond position+geometry). The resolution level is a property of the evidence, not of the string. Tools manufacture overstatement three ways: a formatter that only knows `/`, an in-silico library entry authored at sn-level that a species-level match inherits, and "annotate to the nearest database structure" silently promoting a sum composition to a full structure. sn-position is almost never genuinely measured under CID, so treat every `/` as an unproven `_` until EAD/UVPD/derivatization evidence is in hand. The default rule is: when in doubt, drop a level.

## Structural-Resolution Hierarchy (Separator Semantics)

| Notation | Separator | What was measured | What may NOT be claimed |
|----------|-----------|-------------------|-------------------------|
| `PC 34:1` | space | class + total carbons:double-bonds (accurate mass + isotope + class diagnostic) | the two chains; sn; C=C position |
| `PC 16:0_18:1` | underscore `_` | the two acyl chains (MS/MS acyl losses, RT/ECN-consistent, not an in-source fragment) | which chain is sn-1 vs sn-2 |
| `PC 16:0/18:1` | slash `/` | sn-1/sn-2 assignment (EAD/UVPD/enzymatic - not a CID acyl-loss intensity guess) | C=C position/geometry |
| `PC 16:0/18:1(9Z)` | parentheses | exact double-bond position + cis/trans (OzID/PB/EAD/UVPD) | (full structure) |
| `PC O-34:1` / `PC P-34:1` | `O-` ether / `P-` plasmalogen | ether vs vinyl-ether linkage (diagnostic ion or acid-lability) | a sum composition alone cannot distinguish `P-34:1` from `O-34:2` (vinyl ether = ether + one C=C) |
| `Cer 18:1;O2/16:0` | `;O2` | sphingoid hydroxyl count (old `d18:1`) - measured, not assumed | backbone unsaturation if `d18:1` was a default rather than fragment-confirmed |

Canonicalize every name through Goslin before merging tables or querying LIPID MAPS; never string-match lipid names by hand. Goslin preserves a false `/` faithfully - it is necessary but not sufficient.

## Decision Tree by Question

| Question / situation | Approach | Why |
|----------------------|----------|-----|
| Accurate class-level quantification, high throughput | Shotgun (direct infusion) or HILIC-LC-MS | constant concentration / class bands -> clean ratio to a co-eluting class IS |
| Resolve isobars/isomers, deep low-abundance coverage | RP-LC-MS (± ion mobility) | RT axis adds an identity coordinate; co-elution flags in-source fragments |
| Double-bond position, sn-position, ether/plasmalogen | LC-MS + EAD/OzID/PB/UVPD (± IM) | only these break C=C / glycerol backbone; CID is blind to them |
| Spatial localization | MS-imaging (MS-DIAL 5 spatial mode) | tissue context with predicted-CCS database |
| Need PC acyl chains | negative-mode formate/acetate adduct -> `[M-CH3]-` | `[M+H]+` gives only the m/z 184 head-group ion (class, no chains) |
| Neutral lipids (TG/DG) chains | `[M+NH4]+` adduct | drives neutral-loss-of-fatty-acid fragmentation |
| Suspicious elevated LPC / DG / FA pool | RT co-elution test vs the parent class | an LPC eluting at a PC's RT is an in-source fragment, not biology |
| An apparent odd-chain species (`PC 33:1`) | require MS/MS chain confirmation | usually an in-source fragment or 13C-isotope artifact of an even neighbor |
| Merge names across tools / before a DB lookup | Goslin canonicalization first | abbreviations and separators are tool-specific; hand string-matching corrupts merges |
| Untargeted oxidized-lipid claim | escalate to a targeted, standard-anchored oxylipin panel | untargeted oxidized-lipid IDs are hypotheses; auto-oxidation in the tube fabricates them |

## Load, Normalize, and Run Differential Analysis (lipidr)

**Goal:** Import a quantified lipid table, normalize within class, and find lipids that differ between groups with class/chain-aware output.

**Approach:** Read a Skyline/matrix export into a `LipidomicsExperiment`, attach sample groups, normalize (PQN or class internal standard), then `de_analysis` with an explicit contrast; visualize as a class-faceted volcano.

```r
library(lipidr)

# data_normalized ships with lipidr (PQN-normalized, log2); substitute a real import:
#   d <- read_skyline(list.files(datadir, 'data.csv', full.names = TRUE))
#   d <- add_sample_annotation(d, 'clinical.csv')
#   d <- normalize_pqn(d, measure = 'Area', exclude = 'blank', log = TRUE)
data(data_normalized)

# Contrast references sample-group labels directly; group_col defaults to the first annotation
de_results <- de_analysis(data_normalized, HighFat_water - NormalDiet_water, measure = 'Area')

# logFC.cutoff is on the log2 scale used by limma's topTable inside de_analysis
sig <- significant_molecules(de_results, p.cutoff = 0.05, logFC.cutoff = 1)

plot_results_volcano(de_results, show.labels = FALSE)
```

## Class-Based Internal-Standard Quantification (the non-negotiable)

**Goal:** Convert per-class signal to comparable abundances without baking in class-dependent ionization error.

**Approach:** Ratio each species to a stable-isotope-labeled standard of its OWN class, spiked before extraction so it shares the class's recovery loss; never quantify one class with another class's standard.

```r
# normalize_istd divides each lipid by the internal standard of its matched class.
# Requires one labeled IS per class present in the data (e.g. SPLASH/EquiSPLASH covers ~13 classes).
d_istd <- normalize_istd(data_normalized, measure = 'Area', exclude = 'blank', log = TRUE)

# Class-level summary is only valid WITHIN a class unless per-class response factors were calibrated:
# cross-class molar ratios (e.g. 'PE is 3x PC') carry head-group response bias and are not licensed here.
plot_lipidclass(d_istd, 'sd')
```

## Honest Annotation-Level Assignment (Goslin)

**Goal:** Downgrade any name to the level the evidence supports and verify the claimed level is internally consistent.

**Approach:** Parse with Goslin, read the perceived level, and re-emit at SPECIES (or MOLECULAR_SPECIES) unless sn/C=C evidence exists.

```python
from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel

parser = LipidParser()
lipid = parser.parse('PC 16:0/18:1')      # a slash-claimed name from a tool export

claimed_level = lipid.lipid.info.level    # LipidLevel enum the string asserts
# Without EAD/UVPD evidence, re-emit at the honest molecular-species level (drops the unproven sn):
honest_name = lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES)   # 'PC 16:0_18:1'
sum_name = lipid.get_lipid_string(LipidLevel.SPECIES)                # 'PC 34:1'
```

## Per-Method Failure Modes

### In-source-fragment phantom lyso-/DG-lipidome
- **Trigger:** A labile lipid (PC, TG, plasmalogen) clips an acyl chain in the ESI source before MS1.
- **Mechanism:** The fragment is recorded as an intact precursor; PC->LPC, PE->LPE, TG->DG->MG. The fragment can also be isobaric with a free fatty acid or another class, fabricating phantom signal in several bins; extent is instrument- and tune-dependent.
- **Symptom:** Inflated LPC:PC, DG:TG, or FA pools; an "LPC" eluting at a PC's retention time.
- **Fix:** RT co-elution test (a real LPC elutes at its own ECN position); soften the source (lower in-source CID/transfer energy); treat any large lyso/DG/FA pool as suspect until RT-cleared. Shotgun has no RT axis to run this test - never report elevated lyso-lipids from direct infusion without the in-source-fragment caveat.

### sn-position over-claim
- **Trigger:** A tool exports `/` from CID-only data, or a library back-fills its authored sn arrangement onto a species-level match.
- **Mechanism:** CID acyl-loss intensity bias toward sn-2 is real but small, condition-dependent, and biological samples contain both regioisomers, so the ratio is a blend, not a structure readout.
- **Symptom:** `/`-formatted names with no EAD/UVPD/derivatization evidence file attached.
- **Fix:** Canonicalize through Goslin and re-emit at `MOLECULAR_SPECIES` (`_`); at most state "dominant sn-2 likely X" while reporting `_`.

### Invalid cross-class quantification
- **Trigger:** One global internal standard, or comparing molar abundances across classes after only within-class normalization.
- **Mechanism:** ESI response is head-group-dominated; a PC and a PE at equal moles give signal differing by factors that can exceed an order of magnitude.
- **Symptom:** "Class A is N-fold class B" statements; a single IS used for the whole lipidome.
- **Fix:** One isotope-labeled IS per class; report semi-quantitative within-class unless per-class (and per-adduct) response factors were independently calibrated.

### Ether vs plasmalogen (O-/P-) mis-call
- **Trigger:** Reporting `P-` (plasmalogen) from a sum composition.
- **Mechanism:** `P-34:1` and `O-34:2` share elemental composition (vinyl ether = ether + one C=C); mass cannot distinguish them.
- **Symptom:** Plasmalogen calls with no vinyl-ether diagnostic ion or acid-lability evidence.
- **Fix:** Require a diagnostic fragment or acid-lability test; otherwise report at the level that cannot distinguish them.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| One isotope-labeled IS per lipid class | Köfeler 2021 (good practice); SPLASH/EquiSPLASH | ESI response is head-group-dominated; one global IS miscalibrates every other class |
| EquiSPLASH = 13 deuterated IS at equal 100 µg/mL | Avanti product spec | equimolar comparative use; SPLASH LIPIDOMIX uses unequal physiological concentrations |
| Spike IS before extraction | Köfeler 2021 | only a co-extracted IS corrects class-biased recovery (Folch/Bligh-Dyer/MTBE differ for polar minor classes) |
| MS-DIAL 5 EAD ~14 eV; 96.4% standards delineated, 78.0% sn/OH/C=C correct >1 µM | Takeda 2024 | structural lipidomics yield even with the modern method is incomplete and concentration-dependent |
| ~half of single-software species-level IDs need orthogonal evidence | Köfeler 2021 (Nat Commun) | 510/1108 features, 130/301 PCs & 55/171 TGs violated the ECN/RT model in an audited published set |
| LipidSearch grades: keep A/B/C, drop D | LipidSearch grade definitions | D = mass-only; A = class + all chains = molecular-species level, NOT sn/C=C resolved |
| Shotgun infusion below the aggregation regime | Han/Gross protocol literature | above it lipids aggregate, ESI response goes nonlinear, the IS-ratio assumption collapses |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `could not find function "read_lipidomes"` | non-existent function name | use `read_skyline()` or `as_lipidomics_experiment()` |
| `plot_enrichment` rejects an `enrich.results` argument | wrong signature | `plot_enrichment(de.results, significant.sets, annotation = 'class', measure = 'logFC')`; get sets from `significant_lipidsets()` |
| `lsea(type = 'chain')` errors | no `type` argument | `lsea` tests class/length/unsat sets automatically; rank with `rank.by = c('logFC','P.Value','adj.P.Val')` |
| `de_results$FDR` is NULL | wrong column name | `de_analysis` returns limma columns: `adj.P.Val`, `P.Value`, `logFC` |
| pygoslin `LipidLevel.MOLECULAR_SUBSPECIES` AttributeError | pre-2.0 enum name | current enum is `SPECIES` / `MOLECULAR_SPECIES` / `SN_POSITION` / `STRUCTURE_DEFINED` / `FULL_STRUCTURE` / `COMPLETE_STRUCTURE` |
| Elevated LPC reported from shotgun data | in-source fragmentation with no RT to flag it | add the in-source-fragment caveat; confirm with LC-MS RT co-elution before claiming lyso biology |

## References

- Liebisch G, Vizcaíno JA, Köfeler H, et al. 2013. Shorthand notation for lipid structures derived from mass spectrometry. *J Lipid Res* 54:1523-1530.
- Liebisch G, Fahy E, Aoki J, et al. 2020. Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived lipid structures. *J Lipid Res* 61:1539-1555.
- Fahy E, Subramaniam S, Brown HA, et al. 2005. A comprehensive classification system for lipids. *J Lipid Res* 46:839-861.
- Kopczynski D, Hoffmann N, Peng B, Ahrends R. 2020. Goslin: A Grammar of Succinct Lipid Nomenclature. *Anal Chem* 92:10957-10960.
- Kind T, Liu KH, Lee DY, et al. 2013. LipidBlast in silico tandem mass spectrometry database for lipid identification. *Nat Methods* 10:755-758.
- Takeda H, Takahashi M, Ikeda K, et al. 2024. MS-DIAL 5 multimodal mass spectrometry data mining unveils lipidome complexities. *Nat Commun* 15:9903.
- Mohamed A, Molendijk J, Hill MM. 2020. lipidr: A Software Tool for Data Mining and Analysis of Lipidomics Datasets. *J Proteome Res* 19:2890-2897.
- Köfeler HC, Eichmann TO, Ahrends R, et al. 2021. Quality control requirements for the correct annotation of lipidomics data. *Nat Commun* 12:4771.
- Köfeler HC, Ahrends R, Baker ES, et al. 2021. Recommendations for good practice in MS-based lipidomics. *J Lipid Res* 62:100138.
- McDonald JG, Ejsing CS, Kopczynski D, et al. 2022. Introducing the Lipidomics Minimal Reporting Checklist. *Nat Metab* 4:1086-1088.
- Matyash V, Liebisch G, Kurzchalia TV, et al. 2008. Lipid extraction by methyl-tert-butyl ether for high-throughput lipidomics. *J Lipid Res* 49:1137-1146.
- Bowden JA, Heckert A, Ulmer CZ, et al. 2017. Harmonizing lipidomics: NIST interlaboratory comparison exercise for lipidomics using SRM 1950-Metabolites in Frozen Human Plasma. *J Lipid Res* 58:2275-2288.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream peak detection and feature extraction
- metabolomics/msdial-preprocessing - MS-DIAL alignment and deconvolution upstream of lipid annotation
- metabolomics/metabolite-annotation - General (non-lipid) annotation and confidence levels
- metabolomics/normalization-qc - Sample normalization and QC framing
- metabolomics/statistical-analysis - Multivariate stats on the lipid abundance matrix
