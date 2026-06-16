---
name: bio-proteomics-peptide-identification
description: Peptide-spectrum matching from MS/MS with target-decoy FDR control, framing identification confidence as a property of a ranked list (q-value/PEP) rather than a raw engine score (XCorr, hyperscore, Andromeda, SpecEValue). Covers sequence-database search engines (Comet, MS-GF+, MSFragger, Sage, MaxQuant, MetaMorpheus), concatenated vs separate target-decoy competition, PEP vs q-value, the multi-level FDR cascade, open/mass-tolerant search, rescoring (Percolator, mokapot, MS2Rescore), and pyOpenMS SimpleSearchEngineAlgorithm + FalseDiscoveryRate. Use when identifying peptides from tandem mass spectra and deciding what FDR threshold to act on. Protein grouping and protein-level FDR are protein-inference; PTM site localization is ptm-analysis; DIA peptide-centric scoring is dia-analysis; intensity quant is quantification.
tool_type: mixed
primary_tool: pyOpenMS
---

## Version Compatibility

Reference examples tested with: pyOpenMS 3.1+, pandas 2.2+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Peptide Identification -- Confidence Is a Property of a Ranked List, Not a Single PSM

**"Identify peptides from my MS/MS spectra"** -> Match tandem mass spectra against a protein database, then control false discovery rate by target-decoy competition and act on a q-value -- because a raw match score is meaningless in isolation; only the list-level error rate is interpretable.
- Python: `pyopenms.SimpleSearchEngineAlgorithm().search(...)` for in-process database search, `FalseDiscoveryRate` for q-values
- CLI: `comet`, `msfragger`, `sage`, `MSGFPlus` for high-throughput database searching, `percolator`/`mokapot` for rescoring
- R: `mzID::mzID()` + `flatten()` or `mzR::openIDfile()` + `psms()` to read mzIdentML search results

Scope: this skill owns spectrum-to-peptide matching and PSM/peptide-level FDR. Protein grouping and protein-level (picked) FDR -> protein-inference. PTM site localization and open-search mod discovery follow-up -> ptm-analysis. DIA peptide-centric extraction and scoring -> dia-analysis. FDR-filtered IDs feeding intensities -> quantification. mzML/raw loading -> data-import. OUT OF SCOPE: protein inference, PTM localization scoring, DIA peptide-centric pipelines, label-free/TMT quantification.

## The Single Most Important Modern Insight -- A q-value Is a Verdict on the List, a Raw Score Is Not Even Comparable

1. **Identification confidence is a property of a ranked LIST controlled by target-decoy competition, never a property of one PSM.** The number to act on is a q-value (list-level) or PEP (per-PSM), NOT the engine's raw score. XCorr (Comet), hyperscore (MSFragger/X!Tandem), Andromeda score (MaxQuant), and SpecEValue (MS-GF+) live on different scales, are charge- and length-dependent, and are frequently not even monotone in true probability within a single engine -- which is exactly why rescoring (Percolator/mokapot) exists. "1% FDR" answers "what fraction of the list I keep is wrong," NOT "I am 99% sure of this one ID." The catastrophic error is thresholding on a raw score, or comparing scores across engines.

2. **A q-value is valid only if (a) the decoy DB is a faithful null, (b) targets and decoys competed in ONE concatenated search, and (c) there are enough PSMs for the decoy count to be stable.** Generate decoys at the PROTEIN level then digest (so decoy peptides obey the same enzyme rules), matching the target in size and composition. Concatenated competition gives FDR = (#decoys above threshold) / (#targets above threshold) -- one decoy above threshold estimates one false target. Separate target/decoy searches instead need either the simple Elias-Gygi 2x-decoy estimator FDR = 2 * #decoy / (#target + #decoy) or the more refined mix-max estimator (Keich, Kertesz-Farkas & Noble 2015) -- two distinct options for the separate-search setting, NOT the same formula. Mixing the concatenated and separate forms up is the most common silent FDR error.

3. **PEP and q-value answer different questions; filtering at "PEP <= 0.01" is far stricter than "q <= 0.01."** PEP (posterior error probability, local FDR) is the probability that THIS PSM is wrong; q-value is the FDR of the list cut at this PSM. FDR is the average of PEP over the accepted set (Kall 2008). The worst PSM in a 1%-FDR list typically has a PEP of 10-50%. Use q-value for list cutoffs; use PEP only for per-ID decisions (e.g. picking one PTM site). And PSM-FDR at 1% does NOT give 1% peptide-FDR or 1% protein-FDR -- each level needs its own estimation; hand protein-level control to protein-inference.

## The FDR Vocabulary, Precisely

- **FDR**: the expected proportion of false positives among ALL accepted items at a threshold -- a property of the whole list.
- **q-value**: the minimum FDR at which a given PSM is still accepted; monotone after taking the running minimum from the bottom of the ranked list. Filter on q <= 0.01.
- **PEP (local FDR)**: the probability that THIS PSM is wrong given its score. Local, per-PSM; FDR is the integral of PEP over the accepted set (Kall 2008, "two sides of the same coin").
- **The estimator must match the search mode.** Concatenated target-decoy competition (TDC): FDR = #decoy / #target (no factor 2 -- one best hit per spectrum already resolves the competition). Separate target and decoy searches: either the simple Elias-Gygi 2x-decoy estimator FDR = 2 * #decoy / (#target + #decoy), or the more refined mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). Mix-max is a distinct, calibrated-score procedure for the separate-search setting -- it is NOT a rename of the 2x formula.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---|---|---|---|
| Comet | Eng 2013 | XCorr + E-value; SEQUEST lineage, open-source | Robust default, TPP pipelines; pairs with Percolator |
| X!Tandem | -- | hyperscore + refinement passes | Legacy/free; semi-tryptic refinement niche |
| MS-GF+ | Kim & Pevzner 2014 | SpecEValue via generating-function DP | Calibrated cross-instrument E-value; ETD/CID, low-res, non-standard enzymes |
| MaxQuant / Andromeda | Cox 2011 | binomial probability score; integrated MBR/LFQ/TMT | All-in-one quant pipeline (LFQ, TMT, SILAC); GUI |
| MSFragger | Kong 2017 | hyperscore via fragment-ion indexing (~100x faster) | Open/mass-tolerant search, PTM discovery, huge datasets; core of FragPipe |
| Sage | Lazear 2023 | hyperscore-style, Rust, rescoring-native | Modern scalable open-source pipelines; emits Percolator-ready features |
| MetaMorpheus | Solntsev 2018 | calibration + G-PTM-D multinotch | PTM discovery with built-in calibration; proteoform-aware |
| pFind 3 | Chi 2018 | open-search engine | Maximal unrestricted-PTM/mutation discovery |
| Percolator | Kall 2007 | semi-supervised SVM re-rank on decoy negatives | Boost IDs at fixed FDR; non-tryptic/PTM/large search spaces |
| mokapot | Fondrie & Noble 2021 | Percolator in Python; swappable XGBoost classifier | Python pipelines, Sage output, custom features |
| MS2Rescore + DeepLC + MS2PIP | Declercq 2022; Bouwmeester 2021; Gabriels 2019 | predicted-RT + predicted-intensity rescoring features | Sharpen target/decoy separation; immunopeptidomics |
| Spectral-library search | -- | match empirical reference spectra (intensity + RT) | Faster/more specific for known peptides -> spectral-libraries |
| Protein grouping / protein FDR | Savitski 2015; The 2016 | picked / picked-group FDR | route OUT -> protein-inference |
| PTM site localization | -- | per-site PEP, localization scoring | route OUT -> ptm-analysis |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|---|---|---|
| Standard DDA, clean FDR, scriptable | Comet or Sage + Percolator/mokapot at q <= 0.01 | well-validated; rescoring boosts IDs at fixed FDR |
| Cross-instrument / varied fragmentation / odd enzyme | MS-GF+ | SpecEValue is calibrated so a threshold means the same everywhere |
| Discover unknown PTMs / mass shifts | MSFragger open search (-150..+500 Da) | fragment indexing makes wide-window search feasible; then closed search on discovered mods -> ptm-analysis |
| Huge dataset, reproducible, cloud-scale | Sage (rescoring-native) | Rust speed; emits Percolator features directly |
| All-in-one with quant in the same tool | MaxQuant/Andromeda | integrated LFQ/TMT/SILAC and MBR |
| Non-tryptic (immunopeptidomics, degradomics) | any engine + Percolator/MS2Rescore | rescoring gains are largest where search space explodes |
| Few PSMs (single-protein pulldown) | do NOT trust decoy FDR; inspect spectra manually | decoy counts too noisy below ~hundreds of PSMs |
| Need per-site / per-ID confidence | act on PEP, not q-value | q-value is list-level; PEP is local |

Default when uncertain: concatenated target-decoy search with Comet or Sage, rescore with Percolator/mokapot, filter at q <= 0.01, and hand protein-level FDR to protein-inference.

### Database Search with pyOpenMS

**Goal:** Match tandem mass spectra in an mzML file against a protein FASTA and produce scored PSMs as idXML.

**Approach:** `SimpleSearchEngineAlgorithm` actually scores spectra (the hand-rolled `ProteaseDigestion` loop only digests, it never matches a spectrum). The FASTA must already contain target + decoy sequences concatenated for downstream FDR; decoys carry a recognizable prefix.

```python
from pyopenms import SimpleSearchEngineAlgorithm, IdXMLFile

protein_ids = []
peptide_ids = []
search = SimpleSearchEngineAlgorithm()
# spectra are scored against in-silico fragment ions of every candidate peptide
search.search('sample.mzML', 'human_target_decoy.fasta', protein_ids, peptide_ids)

# protein_ids FIRST in load/store -- the OpenMS argument order is fixed
IdXMLFile().store('search_results.idXML', protein_ids, peptide_ids)
```

### Annotate Target/Decoy and Estimate FDR with pyOpenMS

**Goal:** Convert raw PSM scores into q-values and keep only PSMs at 1% FDR.

**Approach:** `PeptideIndexing` maps each PSM back to proteins and flags target vs decoy from the decoy prefix; `FalseDiscoveryRate.apply` runs the concatenated competition; `IDFilter` keeps q <= 0.01. This is the real pyOpenMS path -- not a hand-rolled decoy/target ratio of unknown provenance.

```python
from pyopenms import PeptideIndexing, FalseDiscoveryRate, IDFilter, FASTAFile

fasta = []
FASTAFile().load('human_target_decoy.fasta', fasta)
indexer = PeptideIndexing()
params = indexer.getParameters()
params.setValue('decoy_string', 'DECOY_')      # must match the decoy prefix in the FASTA
params.setValue('decoy_string_position', 'prefix')
indexer.setParameters(params)
indexer.run(fasta, protein_ids, peptide_ids)   # sets target/decoy flags on every hit

FalseDiscoveryRate().apply(peptide_ids)         # concatenated competition -> per-PSM q-value as the new score
IDFilter().filterHitsByScore(peptide_ids, 0.01) # 0.01 = 1% FDR, the community list-level standard
IDFilter().removeDecoyHits(peptide_ids)
```

### FDR from a Results Table (concatenated competition, made explicit)

**Goal:** Compute q-values from any engine's PSM table when the search was a single concatenated target-decoy search.

**Approach:** Rank by score, walk down accumulating target and decoy counts, FDR = decoys/targets, then take the running minimum from the bottom to get monotone q-values. The decoy/target form is correct ONLY for concatenated competition; separate searches need either the Elias-Gygi 2x-decoy form or the mix-max estimator (Keich, Kertesz-Farkas & Noble 2015).

```python
import pandas as pd

psms = pd.read_csv('search_results.tsv', sep='\t')
psms['is_decoy'] = psms['protein'].str.startswith(('DECOY_', 'REV_', 'XXX_'))
psms = psms.sort_values('score', ascending=False).reset_index(drop=True)

# concatenated target-decoy competition: each decoy above threshold estimates one false target
targets = (~psms['is_decoy']).cumsum()
decoys = psms['is_decoy'].cumsum()
psms['fdr'] = decoys / targets
psms['qvalue'] = psms['fdr'][::-1].cummin()[::-1]   # running min from the bottom -> monotone q-values

kept = psms[(psms['qvalue'] <= 0.01) & (~psms['is_decoy'])]   # 1% list-level FDR
```

## Per-Method Failure Modes

### Concatenated vs separate FDR formula mismatch
**Trigger:** applying #decoy/#target to separately-searched targets and decoys, or 2*decoy/(target+decoy) to concatenated competition.
**Mechanism:** the factor of 2 accounts for false hits that could land in either independent database; concatenated competition already resolves that by a single best hit per spectrum.
**Symptom:** systematically under- or over-estimated FDR; irreproducible ID counts.
**Fix:** confirm the search mode; concatenated -> #decoy/#target; separate -> Elias-Gygi 2x-decoy or the mix-max estimator (Keich, Kertesz-Farkas & Noble 2015). In Percolator, mix-max is the default for separate-search input and `-Y`/`--post-processing-tdc` selects target-decoy competition instead; concatenated input forces TDC automatically.

### Thresholding on raw engine score
**Trigger:** filtering on XCorr/hyperscore/Andromeda score, or comparing scores from two engines.
**Mechanism:** scores are uncalibrated, charge/length-dependent, and not monotone in true probability.
**Symptom:** different cutoffs admit different real FDRs; cross-engine merges nonsensical.
**Fix:** always convert to q-value (or SpecEValue/PEP) first; rescore with Percolator/mokapot.

### Decoy FDR on too few PSMs
**Trigger:** reporting "0% FDR" from a single-protein pulldown or tiny PSM list.
**Mechanism:** the decoy count is a noisy Poisson-like estimate; zero observed decoys does not mean zero false targets.
**Symptom:** spuriously confident IDs from small experiments.
**Fix:** below ~hundreds of PSMs, inspect spectra manually; do not act on the decoy q-value.

### Open-search results used for clean FDR or quant
**Trigger:** taking IDs from a wide-window (-150..+500 Da) search as final, FDR-controlled results.
**Mechanism:** wide windows admit "free" mass shifts that inflate random matches; the target-decoy null differs per mass-shift bin.
**Symptom:** inflated, unreliable FDR on open-search output.
**Fix:** treat open search as discovery; follow with a closed search restricted to the discovered mods -> ptm-analysis.

### Rescoring overfitting
**Trigger:** custom features that leak label information, or training without proper cross-validation.
**Mechanism:** the model learns the decoys, making rescored FDR optimistic.
**Symptom:** ID counts jump but downstream validation fails.
**Fix:** use Percolator/mokapot default cross-validation; predicted-feature rescoring (DeepLC/MS2PIP) is safer; validate with entrapment for high-stakes claims (Wen 2025).

### DIA tool FDR taken at face value
**Trigger:** trusting a DIA tool's reported 1% peptide/protein FDR.
**Mechanism:** entrapment shows several DIA tools do not reliably control FDR (Wen 2025).
**Symptom:** real error rate exceeds the reported FDR.
**Fix:** validate with entrapment for high-stakes DIA claims -> dia-analysis.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Precursor tolerance 10-20 ppm (high-res Orbitrap) | -- | matches FT mass accuracy; tighter = fewer random candidates at fixed FDR |
| Precursor tolerance -150..+500 Da (open search) | Kong 2017 | captures arbitrary PTM/mutation shifts; feasible only with fragment indexing |
| Fragment tolerance 0.02 Da (HCD Orbitrap) / 0.6 Da (ion-trap CID) | -- | instrument-dependent; 0.6 Da on Orbitrap discards resolving power |
| Missed cleavages 2 | -- | covers incomplete trypsin digestion without exploding search space |
| PSM/peptide FDR 1% (q <= 0.01) | Elias & Gygi 2007 | community standard; list-level error, not per-PSM |
| Decoy:target ratio 1:1 | Elias & Gygi 2007 | standard; unequal ratios need formula correction |
| Min PSMs for trustworthy decoy FDR: hundreds+ | -- | below this the decoy count is too noisy |
| Variable mods per peptide <= 2-3 | -- | each variable mod multiplies search space and random-match rate |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| pyOpenMS "search" returns peptides but never scores spectra | used `ProteaseDigestion`, which only digests a FASTA | use `SimpleSearchEngineAlgorithm().search(mzML, fasta, protein_ids, peptide_ids)` |
| `IdXMLFile().load/store` argument error | wrong order | protein_ids FIRST: `IdXMLFile().load(path, protein_ids, peptide_ids)` |
| FDR ignores decoys / all q-values 0 | decoys not annotated before `FalseDiscoveryRate` | run `PeptideIndexing` with matching `decoy_string` first |
| R: `MSnbase::readMzIdData` not found | that function name does not exist | use `mzID::mzID(file)` + `flatten()`, or `mzR::openIDfile()` + `psms()` (PSMatch/Spectra is the modern path) |
| Percolator q-method mismatched to search mode | mix-max is the default for separate-search input | for separate searches, mix-max (default) or `-Y`/`--post-processing-tdc` for target-decoy competition; concatenated input forces TDC automatically; use `--picked-protein` for protein FDR |
| 1% PSM FDR assumed to give 1% protein FDR | each level needs its own estimation | estimate protein-level (picked) FDR -> protein-inference |
| "PEP <= 0.01" returns far fewer IDs than expected | PEP is per-PSM and far stricter than q-value | filter list cutoffs on q-value; reserve PEP for per-ID decisions |

## References

- Elias, J.E. & Gygi, S.P. 2007. Target-decoy search strategy for increased confidence in large-scale protein identifications by mass spectrometry. *Nature Methods* 4(3):207-214.
- Keich, U., Kertesz-Farkas, A. & Noble, W.S. 2015. Improved false discovery rate estimation procedure for shotgun proteomics. *Journal of Proteome Research* 14(8):3148-3161.
- Kall, L., Canterbury, J.D., Weston, J., Noble, W.S. & MacCoss, M.J. 2007. Semi-supervised learning for peptide identification from shotgun proteomics datasets. *Nature Methods* 4(11):923-925.
- Kall, L., Storey, J.D., MacCoss, M.J. & Noble, W.S. 2008. Posterior error probabilities and false discovery rates: two sides of the same coin. *Journal of Proteome Research* 7(1):40-44.
- Eng, J.K., Jahan, T.A. & Hoopmann, M.R. 2013. Comet: an open-source MS/MS sequence database search tool. *Proteomics* 13(1):22-24.
- Kim, S. & Pevzner, P.A. 2014. MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications* 5:5277.
- Cox, J., Neuhauser, N., Michalski, A., Scheltema, R.A., Olsen, J.V. & Mann, M. 2011. Andromeda: a peptide search engine integrated into the MaxQuant environment. *Journal of Proteome Research* 10(4):1794-1805.
- Kong, A.T., Leprevost, F.V., Avtonomov, D.M., Mellacheruvu, D. & Nesvizhskii, A.I. 2017. MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry-based proteomics. *Nature Methods* 14(5):513-520.
- Lazear, M.R. 2023. Sage: an open-source tool for fast proteomics searching and quantification at scale. *Journal of Proteome Research* 22(11):3652-3659.
- Solntsev, S.K., Shortreed, M.R., Frey, B.L. & Smith, L.M. 2018. Enhanced global post-translational modification discovery with MetaMorpheus. *Journal of Proteome Research* 17(5):1844-1851.
- Chi, H., Liu, C., Yang, H. et al. 2018. Comprehensive identification of peptides in tandem mass spectra using an efficient open search engine. *Nature Biotechnology* 36:1059-1061.
- Fondrie, W.E. & Noble, W.S. 2021. mokapot: fast and flexible semisupervised learning for peptide detection. *Journal of Proteome Research* 20(4):1966-1971.
- Bouwmeester, R., Gabriels, R., Hulstaert, N., Martens, L. & Degroeve, S. 2021. DeepLC can predict retention times for peptides that carry as-yet unseen modifications. *Nature Methods* 18:1363-1369.
- Gabriels, R., Martens, L. & Degroeve, S. 2019. Updated MS2PIP web server delivers fast and accurate MS2 peak intensity prediction for multiple fragmentation methods, instruments and labeling techniques. *Nucleic Acids Research* 47(W1):W295-W299.
- Declercq, A., Bouwmeester, R., Hirschler, A., Carapito, C., Degroeve, S., Martens, L. & Gabriels, R. 2022. MS2Rescore: data-driven rescoring dramatically boosts immunopeptide identification rates. *Molecular & Cellular Proteomics* 21(8):100266.
- Wen, B., Freestone, J., Riffle, M., MacCoss, M.J., Noble, W.S. & Keich, U. 2025. Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment. *Nature Methods* 22:1454-1463.

## Related Skills

- protein-inference - Group peptides to protein groups and control protein-level (picked) FDR
- ptm-analysis - Open/variable-mod search follow-up and per-site PTM localization
- dia-analysis - DIA peptide-centric extraction and scoring; entrapment FDR validation
- quantification - FDR-filtered IDs feed label-free/TMT intensity quantification
- spectral-libraries - Empirical and predicted spectral-library search as an ID alternative
- data-import - Load mzML/raw MS data before identification
- database-access/uniprot-access - Build the target FASTA (canonical vs isoform, contaminants)
