---
name: bio-metabolomics-metabolite-annotation
description: Turns untargeted LC-MS/MS features (m/z, RT, MS/MS) into confidence-stratified metabolite annotations using spectral-library matching (matchms), in-silico tools (SIRIUS/CSI:FingerID, MetFrag) and molecular networking, and assigns a defensible MSI/Schymanski confidence level to each. Use when naming detected features, scoring MS/MS against a reference library, running SIRIUS, or deciding what confidence level an evidence set actually supports. For upstream feature extraction see metabolomics/xcms-preprocessing and metabolomics/msdial-preprocessing; for downstream enrichment that must respect these levels see metabolomics/pathway-mapping; for lipid-specific structural annotation see metabolomics/lipidomics.
tool_type: mixed
primary_tool: matchms
---

## Version Compatibility

Reference examples tested with: matchms 0.33+, SIRIUS 6.x, MetFrag 2.5+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Spectral matching needs precursor m/z on every MS/MS spectrum (`add_precursor_mz` filter) or ModifiedCosine silently returns zeros. Level 1 needs an authentic standard run in the same lab under the same method; no software output can substitute for it.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolite Annotation

**"Annotate my metabolomics features with compound identities"** -> Map each feature's m/z and MS/MS to candidate structures, then attach an explicit confidence level to every name.
- Python: `matchms.calculate_scores()` for library matching (matchms)
- CLI: `sirius ... formulas fingerprints structures canopus` for in-silico formula/structure/class (SIRIUS)

## The Single Most Important Insight -- An Annotation Is a Hypothesis Carrying a Confidence Level, Not an Identification

A metabolite name without a stated MSI/Schymanski level is scientifically incomplete. The inference chain m/z -> formula -> structure -> isomer-resolved identity is three separate lossy steps, each needing its own orthogonal evidence axis. A database hit supplies a name, not evidence: with no MS/MS or RT to back it, it is Schymanski Level 4 (formula) at best, often Level 5 (a feature of interest). A high cosine score ranks candidates; it never proves one. Only an in-house authentic standard, same method, with MS, MS/MS, and RT all matching reaches Level 1 ("identification") -- everything else is an honest hypothesis. The field's recurring sin is laundering Level 2/3 hypotheses into Level-1 prose; the canonical worked example is phenylacetylglutamine being reported as phenylacetylglycine in nearly half of NMR studies (Theodoridis 2023). Assign the lowest level the evidence honestly supports and report which database/version was searched.

## Confidence-Level Taxonomy (MSI and Schymanski)

| Schymanski | MSI | Name | Evidence required |
|---|---|---|---|
| Level 1 | 1 | Confirmed structure | In-house authentic standard, same method: MS + MS/MS + RT all match. The only "identification". |
| Level 2a | 2 | Probable structure (library) | MS/MS matches a reference library spectrum; no in-house standard. |
| Level 2b | 2 | Probable structure (diagnostic) | Diagnostic fragments / RT / ionization consistent with exactly one structure; no reference spectrum. |
| Level 3 | 3 | Tentative candidate(s) | Evidence narrows to a structure class or candidate set but isomers remain unresolved. |
| Level 4 | -- | Unequivocal formula | MS1 accurate mass + isotope pattern + adduct logic assign one formula; no structure. |
| Level 5 | 4 | Exact mass | A feature of interest; nothing assigned. |

Promote one level per orthogonal evidence axis that survives scrutiny; cap at Level 2 unless an in-house standard exists. CSI:FingerID and library matching recover constitution only -- no stereochemistry, so enantiomer/regiochemistry claims cannot come from MS/MS.

## Tool Roles

| Tool | Core idea | Output | Best for |
|---|---|---|---|
| matchms (CosineGreedy / ModifiedCosine / spectral entropy) | Score query MS/MS against library spectra | Ranked library hits + matched-peak count | Level 2a when a library spectrum exists |
| SIRIUS + ZODIAC | Fragmentation trees + isotope pattern, dataset-wide formula re-ranking | Ranked molecular formula | Formula (Level 4); the reliable part of SIRIUS |
| CSI:FingerID + COSMIC | Predict fingerprint, search structure DB, calibrated confidence | Ranked structures + FDR-controllable score | Level 2b/3 structure when COSMIC FDR is set |
| CANOPUS | Predict compound class directly from MS2 | ClassyFire + NPClassifier class | Level 3 class for unknowns; often the most honest output |
| MetFrag | Bond-disconnection scoring of candidate list | Explainable fragment-supported ranks | Transparent, scriptable, custom DBs, RT term |
| FBMN (GNPS2) + MS2Query | Modified-cosine network / ML analogue search | Edges = "related to" | Analogue propagation (Level 3 scaffold hypothesis) |

## Decision Tree: Evidence Available -> Tool -> Achievable Level

| Situation | Do | Achievable level |
|---|---|---|
| In-house authentic standard, same method, MS+MS/MS+RT match | Confirm against standard | Level 1 |
| MS/MS available, library spectrum likely exists | matchms library match (entropy or modified cosine) | Level 2a |
| MS/MS available, no library spectrum | SIRIUS formulas + CSI:FingerID + CANOPUS, or MetFrag | Level 2b/3 (formula Level 4) |
| Need class only / compound absent from all DBs | CANOPUS (class); MSNovelist (de novo SMILES) | Level 3 |
| Find analogues / propagate across a network | FBMN on GNPS2 + MS2Query | Level 3 (scaffold hypothesis) |
| Only MS1 m/z + isotopes + clean adduct | Formula assignment (SIRIUS / seven golden rules) | Level 4 |
| Bare m/z, no orthogonal evidence | Report as a feature | Level 5 |
| Biology hinges on a specific isomer / stereocenter | Demand a standard or orthogonal method (NMR, chiral assay) | MS alone insufficient |

## Match MS/MS Against a Spectral Library

**Goal:** Rank library candidates for each query spectrum and attach the matched-peak count, not just the score.

**Approach:** Harmonize metadata, normalize intensities, add precursor m/z, score with ModifiedCosine (analogue-aware) or spectral entropy (identity), then keep only hits above both a score and a matched-peak floor.

```python
from matchms import calculate_scores
from matchms.filtering import default_filters, normalize_intensities, add_precursor_mz
try:
    from matchms.similarity import ModifiedCosineGreedy as ModifiedCosine  # matchms 0.33+
except ImportError:
    from matchms.similarity import ModifiedCosine          # matchms <= 0.32

def prepare(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = add_precursor_mz(spectrum)  # required for ModifiedCosine or scores are zero
    return normalize_intensities(spectrum)

queries = [prepare(s) for s in queries_raw]
references = [prepare(s) for s in references_raw]

scores = calculate_scores(references, queries, ModifiedCosine(tolerance=0.005))

# CosineGreedy/ModifiedCosine return a structured array; the field names are
# class-prefixed and version-dependent (e.g. 'ModifiedCosineGreedy_score' in 0.33),
# so derive them from the dtype rather than hard-coding.
for query in queries:
    pairs = scores.scores_by_query(query)
    score_field, match_field = pairs[0][1].dtype.names
    ref, hit = max(pairs, key=lambda pair: pair[1][score_field])
    if hit[score_field] >= 0.7 and hit[match_field] >= 6:  # score floor + peak-count floor (GNPS defaults)
        print(ref.get('compound_name'), hit[score_field], hit[match_field])  # Level 2a candidate
```

## Run SIRIUS for Formula, Structure, and Class

**Goal:** Annotate features that have no library spectrum, reporting formula and class with more trust than top-1 structure.

**Approach:** Run the SIRIUS subcommand chain on one project space; trust ZODIAC-refined formula over CSI:FingerID structure, and only report a structure as confident when a COSMIC FDR threshold is set.

```bash
# SIRIUS 6 is a multi-command pipeline on one line. A free academic account/license
# is required (since v5); log in once, then the project space persists across runs.
# Credential flags vary by version; run `sirius login --help` to confirm (commonly `-u <email>`).
sirius login -u "$SIRIUS_USER"

sirius --input features.mgf --project ./sirius_project \
    formulas --profile orbitrap \
    fingerprints \
    structures --database bio \
    canopus \
    write-summaries --output ./sirius_summary
# Verify exact subcommand spelling with `sirius <command> --help`: formulas/fingerprints/
# structures/canopus changed plural/singular and options between v5 and v6.
# --database (on structures) is a scientific choice: 'bio' raises plausibility but cannot
# return a novel metabolite; 'pubchem' maximizes recall but floods implausible isomers.
```

## Assemble an Evidence-to-Level Call

**Goal:** Collapse a feature's evidence set into a single defensible confidence level.

**Approach:** Start at Level 5 and promote per surviving orthogonal axis; an authentic standard is the only path to Level 1.

```python
def assign_level(evidence):
    if evidence.get('authentic_standard_same_method'):
        return 1
    if evidence.get('library_match') and evidence['library_match']['score'] >= 0.7 and evidence['library_match']['matches'] >= 6:
        return '2a'  # reference library spectrum, no in-house standard
    if evidence.get('diagnostic_fragments') and evidence.get('single_structure_consistent'):
        return '2b'
    if evidence.get('candidate_set') or evidence.get('canopus_class') or evidence.get('network_propagated'):
        return 3  # isomers unresolved, class only, or "related to" an annotated node
    if evidence.get('unambiguous_formula'):
        return 4  # MS1 + isotopes + adduct logic, no structure
    return 5
```

## Per-Method Failure Modes

### Cosine score is not identity
- **Trigger:** Reporting a name because a single high cosine/modified-cosine score came back.
- **Mechanism:** Cosine rewards shared fragment peaks, and fragments are substructures many distinct molecules share; a high score on few peaks aligns with thousands of unrelated compounds.
- **Symptom:** Confident name that an isomer or scaffold-sharing compound would have produced identically.
- **Fix:** Require a matched-peak floor (>=6) alongside the score (>=0.7); prefer spectral entropy for identity; report Level 2a, not Level 1.

### The isomer wall
- **Trigger:** Claiming a specific positional/stereo/regio isomer from MS/MS.
- **Mechanism:** Constitutional isomers frequently fragment identically; enantiomers have near-identical CID spectra; CSI:FingerID is constitution-only.
- **Symptom:** A specific structure reported where multiple isomers fit the data equally.
- **Fix:** Report Level 3 unless RT or CCS breaks the tie (CCS needs ~0.5-0.6% separation); for biology hinging on the isomer, use NMR or a co-eluting standard.

### In-source fragments and adduct cascades corrupt the input
- **Trigger:** Annotating and counting features before collapsing ion families.
- **Mechanism:** In-source fragmentation creates phantom MS1 features; assuming the wrong adduct shifts the neutral mass and corrupts every downstream candidate, producing a confident, internally consistent, wrong answer.
- **Symptom:** Over-counted "compounds", the same molecule named several ways, invented biology.
- **Fix:** Group ion families (CAMERA / Ion Identity Molecular Networking / khipu) before annotation; never quote feature counts as compound counts.

### Database-mapping inflation poisons pathway analysis
- **Trigger:** Feeding all candidate IDs of an ambiguous feature into enrichment.
- **Mechanism:** One ambiguous m/z maps to many compound IDs across different pathways, so a single uncertain feature lights up several pathways (phantom enrichment).
- **Symptom:** Inflated pathway significance traceable to Level-3 features voting as if they were several confirmed compounds.
- **Fix:** Carry annotation uncertainty (candidate sets, levels) into enrichment; prefer mass-level or probabilistic methods that do not multiply ambiguous IDs (see metabolomics/pathway-mapping; mummichog deliberately avoids prior ID).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| Cosine/modified-cosine >= 0.7 AND >= 6 matched peaks | GNPS defaults (Wang 2016) | Suppresses promiscuous low-complexity spectra that hairball the network. |
| Spectral entropy >= 0.75 -> FDR < 10% | Li 2021 (natural-products benchmark) | Dataset-dependent, NOT a universal constant; entropy beats dot product for identity. |
| MS1 mass error <= 5 ppm (HRMS) | HRMS convention | Tighter than the 10 ppm older default; pairs with isotope-pattern filter. |
| Isotope-pattern ~2% abundance accuracy | Kind & Fiehn 2006 | Removes >95% of false formula candidates even at 3 ppm -- orthogonal info, not better mass accuracy, fixes formula. |
| COSMIC 0.94 / 0.64 / 0.34 ~ 5 / 10 / 20% FDR | Hoffmann 2022 | Calibrated confidence on CSI:FingerID structures; raw top-1 with no COSMIC is Level 3. |
| Predicted CCS within ~3-5% of measured | AllCCS / IMS benchmarks (Zhou 2020) | Use CCS as a falsifier (rejects candidates), not as positive proof of identity. |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| ModifiedCosine scores all zero | Missing precursor m/z on spectra | Apply `add_precursor_mz` filter to both references and queries first. |
| `AttributeError: 'Scores' has no attribute 'scores'` | Indexing `scores.scores[...]` (old tutorials) | Use `scores.scores_by_query(query)` or `scores.to_array(name=...)`. |
| `ValueError: no field of name <X>_score` | Field names are class-prefixed and version-dependent | Read `pair[1].dtype.names` for the score/matches field names rather than hard-coding. |
| `ImportError: cannot import name 'ModifiedCosine'` | Renamed to `ModifiedCosineGreedy` in matchms 0.33 | Try the new name with an ImportError fallback to the old. |
| `sirius formula` not found | v5 used singular subcommands; v6 uses `formulas` | Run `sirius --help`; verify plural/singular per installed version. |
| SIRIUS exits at login | Account/license required since v5 | `sirius login` once with a free academic account before the chain. |
| Pathway enrichment lights up everywhere | Ambiguous features mapped to many DB IDs | Collapse ion families and carry levels into enrichment (metabolomics/pathway-mapping). |

## References

- Sumner LW, et al. 2007. Proposed minimum reporting standards for chemical analysis (CAWG MSI). *Metabolomics* 3:211-221.
- Schymanski EL, Jeon J, Gulde R, Fenner K, Ruff M, Singer HP, Hollender J. 2014. Identifying small molecules via high resolution mass spectrometry: communicating confidence. *Environ Sci Technol* 48:2097-2098.
- Li Y, Kind T, Folz J, Vaniya A, Mehta SS, Fiehn O. 2021. Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. *Nat Methods* 18:1524-1531.
- Dührkop K, Fleischauer M, Ludwig M, Aksenov AA, Melnik AV, Meusel M, Dorrestein PC, Rousu J, Böcker S. 2019. SIRIUS 4: a rapid tool for turning tandem mass spectra into metabolite structure information. *Nat Methods* 16:299-302.
- Dührkop K, Shen H, Meusel M, Rousu J, Böcker S. 2015. Searching molecular structure databases with tandem mass spectra using CSI:FingerID. *PNAS* 112:12580-12585.
- Dührkop K, et al. 2021. Systematic classification of unknown metabolites using high-resolution fragmentation mass spectra (CANOPUS). *Nat Biotechnol* 39:462-471.
- Hoffmann MA, et al. 2022. High-confidence structural annotation of metabolites absent from spectral libraries (COSMIC). *Nat Biotechnol* 40:411-421.
- Ruttkies C, Schymanski EL, Wolf S, Hollender J, Neumann S. 2016. MetFrag relaunched: incorporating strategies beyond in silico fragmentation. *J Cheminform* 8:3.
- Wang M, Carver JJ, Phelan VV, et al. 2016. Sharing and community curation of mass spectrometry data with GNPS. *Nat Biotechnol* 34:828-837.
- Nothias LF, Petras D, Schmid R, et al. 2020. Feature-based molecular networking in the GNPS analysis environment. *Nat Methods* 17:905-908.
- Kind T, Fiehn O. 2006. Metabolomic database annotations via query of elemental compositions: mass accuracy is insufficient even at less than 1 ppm. *BMC Bioinformatics* 7:234.
- Zhou Z, et al. 2020. Ion mobility collision cross-section atlas for known and unknown metabolite annotation in untargeted metabolomics (AllCCS). *Nat Commun* 11:4334.
- Theodoridis G, Gika H, Raftery D, Goodacre R, Plumb RS, Wilson ID. 2023. Ensuring fact-based metabolite identification in LC-MS-based metabolomics. *Anal Chem* 95:3909-3916.
- Huber F, Verhoeven S, Meijer C, et al. 2020. matchms - processing and similarity evaluation of mass spectrometry data. *J Open Source Softw* 5:2411.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream feature extraction (m/z, RT, intensity table)
- metabolomics/msdial-preprocessing - Alternative feature extraction and deconvolution
- metabolomics/pathway-mapping - Downstream enrichment that must respect these confidence levels
- metabolomics/lipidomics - Lipid-specific annotation and structural resolution
- proteomics/spectral-libraries - Related spectral-matching concepts (closed-world peptide search)
