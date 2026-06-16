---
name: bio-proteomics-protein-inference
description: Groups proteins from peptide identifications and controls protein-level FDR, framing inference as a chosen explanation (parsimony or a probability model) of underdetermined peptide evidence rather than a measurement. Reports protein GROUPS (proteins indistinguishable by observed peptides) with a leading protein, not flat lists. Covers shared-vs-unique peptides, indistinguishable/subsumable proteins, parsimony vs probabilistic (ProteinProphet, EPIFANY) vs razor inference, picked-protein and picked-group FDR, and why the two-peptide rule is wrong. Use when resolving which proteins are present from a peptide list, building protein groups, or estimating protein-level FDR. PSM/peptide FDR and search engines are peptide-identification; razor-vs-unique quant consequences are quantification; isoform/proteoform resolution is top-down and out of scope.
tool_type: mixed
primary_tool: pyOpenMS
---

## Version Compatibility

Reference examples tested with: pyOpenMS 3.1+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The pyOpenMS protein-inference class names have varied across releases. Confirm the exact spelling at the installed version with `help(pyopenms.EpifanyAlgorithm)` and `help(pyopenms.BasicProteinInferenceAlgorithm)` before relying on the reference code.

# Protein Inference -- A Chosen Explanation of Peptide Evidence, Reported as Groups

**"Tell me which proteins are present from my identified peptides"** -> Assign the observed peptides to a minimal or probability-weighted set of proteins, reported as groups of indistinguishable proteins with a leading accession -- because bottom-up MS measures peptides, and the protein set behind them is inferred, not observed.
- Python: `pyopenms.BasicProteinInferenceAlgorithm().run(peptide_ids, protein_ids)` for parsimony grouping
- Python: `pyopenms.EpifanyAlgorithm` (TOPP tool `Epifany`) for Bayesian belief-propagation inference
- CLI: `ProteinProphet` (TPP) for EM-based probabilistic inference; `Philosopher filter` for FragPipe FDR

Scope: this skill OWNS peptide-to-protein grouping, the indistinguishable/subsumable distinction, the leading-protein convention, inference-method choice, and protein/protein-group FDR. PSM-level and peptide-level FDR plus the search engines that produce the peptide list -> peptide-identification. The quantitative fallout of razor vs unique peptides on protein abundance -> quantification. OUT OF SCOPE: resolving splice isoforms, single-AA variants, or PTM-defined proteoforms (bottom-up groups cannot separate them; that is top-down / proteoform work).

## The Single Most Important Modern Insight -- Protein Inference Is Underdetermined, So the Honest Unit Is a Group, Not a List

1. **The protein set is not uniquely recoverable from peptides, so a protein group -- not a flat protein list -- is the only honest reporting unit.** Many peptides are shared across paralogs, gene families, and isoforms, so distinct protein sets can explain the same peptide evidence equally well. The inference picks ONE explanation under an assumption (parsimony, or a probability model); proteins that the observed peptides cannot tell apart (indistinguishable) MUST be reported as one group with a designated leading protein. A flat list double-counts indistinguishable proteins and breaks target/decoy symmetry at the protein level, silently corrupting FDR.

2. **Protein FDR is its own estimation problem that INFLATES on large data; the fix is PICKED FDR, not the PSM formula reused.** Controlling PSM-FDR at 1% does not give 1% protein-FDR. A deep run has many false PSMs in absolute terms, and each can nucleate a one-hit-wonder false protein; because true proteins accumulate many peptides while false proteins are hit once, the naive protein-FDR balloons to 10-30% on deep datasets. Savitski 2015 picked-protein FDR pairs each target protein with its decoy and keeps only the higher-scoring of the pair before counting, removing the target/decoy asymmetry; The & Kall 2016 extends this to the group level (picked-group FDR), which is required because parsimony grouping is anticonservative otherwise.

3. **The two-peptide rule is wrong -- it increases protein FDR and discards real proteins.** Requiring >=2 peptides per protein (Gupta & Pevzner 2009, "A strike against the two-peptide rule") removes MORE target proteins than decoy proteins, so it raises protein-level FDR rather than lowering it, while throwing away legitimate low-abundance single-peptide IDs. Replace the blanket rule with: control protein-level (picked) FDR, then judge single-peptide IDs by their score, not their peptide count.

## Vocabulary the Rest of This Depends On

- Shared (degenerate) peptide: maps to >1 protein in the searched database. Cannot, alone, distinguish which protein is present.
- Unique peptide: maps to exactly one protein -- the only direct evidence for a specific protein. "Unique" is DATABASE-RELATIVE: a peptide unique against SwissProt may be shared against TrEMBL+isoforms+contaminants. Always document the exact database (isoforms, contaminants, decoys included).
- Indistinguishable proteins: explained by the SAME set of observed peptides -> one group, never two confident IDs.
- Subset / subsumable protein: its observed peptides are a subset of another protein's -> parsimony drops it (the larger protein explains everything it would).
- Leading / representative protein: the group's reported accession. Convention: most peptides, then highest score, then SwissProt canonical over TrEMBL. Downstream tables key on this accession but must retain group membership -- "protein P12345" usually means "the group led by P12345".
- Protein group vs proteoform: a group is an inference artifact (proteins lumped because peptides cannot separate them); a proteoform is a real molecular species (one gene product with a specific sequence + PTM + cleavage state). Bottom-up groups DO NOT resolve proteoforms -- claiming "isoform X present" from a shared-peptide group is overreach.

## Tool Taxonomy

| Tool / method | Citation | Mechanism / role | When |
|---------------|----------|------------------|------|
| Parsimony (Occam) | -- | Greedy minimal protein set explaining all peptides | Fast default; ties broken arbitrarily; anticonservative group-FDR on large data unless picked |
| ProteinProphet | Nesvizhskii 2003 | EM APPORTIONS shared peptides across candidate proteins, weighted by other evidence | TPP / FragPipe pipelines; the classic probabilistic standard |
| EPIFANY | Pfeuffer 2020 | Bayesian network over the peptide-protein graph, loopy belief propagation + convolution trees | OpenMS-recommended modern inference; strong at controlled protein-group FDR |
| Fido | -- | Bayesian generative model (Percolator `--protein`) | Percolator pipelines; superseded by picked-protein for FDR |
| Razor peptide | -- | Shared peptide assigned winner-take-all to the group with most evidence (MaxQuant) | MaxQuant default; ID-fine but distorts QUANT (route to quantification) |
| Picked-protein FDR | Savitski 2015 | Pair target with its decoy, keep the higher-scoring of the pair, then count decoys | Protein-level FDR on any non-trivial dataset |
| Picked-group FDR | The & Kall 2016 | Picking applied at the protein-GROUP level | When the inference unit is the group (the correct unit on deep data) |
| All-proteins / inclusive | -- | Report every protein any peptide could come from | Almost never; massive false-positive protein inflation |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Standard DDA run, OpenMS-based pipeline | `EpifanyAlgorithm` (or `BasicProteinInferenceAlgorithm` for parsimony) + picked-group FDR | Modern, group-FDR aware; well-calibrated on benchmarks |
| MaxQuant output (`proteinGroups.txt`) | Parse groups as-is; quantify on UNIQUE peptides | Groups already inferred; razor quant is the trap, not the inference |
| FragPipe / TPP pipeline | ProteinProphet inference + Philosopher/Philosopher-style FDR filtering | Native EM apportionment + 2-level FDR |
| Deep dataset (many thousands of proteins) | Picked-GROUP FDR, NOT naive decoy/target | Naive protein-FDR inflates to 10-30% from one-hit-wonders |
| Sensitive differential abundance downstream | Quantify on unique peptides only -> quantification | Razor assignment can flip between conditions and fake DE |
| Want isoform-level answers | Stop -- route to top-down / proteoform methods | Bottom-up groups cannot resolve proteoforms |
| Few PSMs (single-protein pulldown) | Report evidence, do not trust a "0% protein FDR" | Target-decoy FDR is meaningless at tiny counts |

Default when uncertain: run parsimony grouping (`BasicProteinInferenceAlgorithm` with `annotate_indistinguishable_groups`), report protein GROUPS with a leading accession, and control protein-GROUP FDR with picked-group FDR at 1%. Do NOT impose a two-peptide rule.

### Group Proteins by Parsimony with pyOpenMS

**Goal:** Turn an FDR-filtered peptide identification list into protein groups with a leading protein, resolving shared-peptide ambiguity.

**Approach:** Load the idXML from peptide identification, run the parsimony algorithm with indistinguishable-group annotation on, then read the inferred groups off the protein identification run.

```python
from pyopenms import IdXMLFile, BasicProteinInferenceAlgorithm

protein_ids = []
peptide_ids = []
# protein_ids is FIRST in both load() and store() for IdXMLFile
IdXMLFile().load('peptides_1pct_fdr.idXML', protein_ids, peptide_ids)

inference = BasicProteinInferenceAlgorithm()
params = inference.getParameters()
# annotate_indistinguishable_groups reports indistinguishable proteins as ONE group
params.setValue('annotate_indistinguishable_groups', 'true')
inference.setParameters(params)
inference.run(peptide_ids, protein_ids)

# indistinguishable groups live on the protein identification run
for prot_id in protein_ids:
    for group in prot_id.getIndistinguishableProteins():
        leading = group.accessions[0]  # convention: highest-evidence accession first
        print(leading, group.probability, list(group.accessions))
```

### Bayesian Inference + Group FDR with EPIFANY

**Goal:** Assign calibrated protein/group posteriors and control protein-group FDR with a probability model rather than greedy parsimony.

**Approach:** EPIFANY consumes idXML whose PSMs already carry posterior error probabilities (from Percolator or IDPosteriorErrorProbability), then propagates belief over the peptide-protein graph. The TOPP tool is reliably named `Epifany`; the pyOpenMS class spelling has varied across releases, so introspect first.

```python
import pyopenms
from pyopenms import IdXMLFile

# CONFIRM the class name at the installed version before use:
#   help(pyopenms.EpifanyAlgorithm)
algo_cls = getattr(pyopenms, 'EpifanyAlgorithm')

protein_ids = []
peptide_ids = []
IdXMLFile().load('peptides_with_pep.idXML', protein_ids, peptide_ids)

algo = algo_cls()
# EPIFANY expects PSM posteriors as input; greedy_group_resolution controls
# whether shared peptides are razor-resolved after inference
algo.inferPosteriorProbabilities(protein_ids, peptide_ids, False)

for prot_id in protein_ids:
    for group in prot_id.getIndistinguishableProteins():
        print(group.accessions[0], group.probability)
```

### Picked Protein-Group FDR

**Goal:** Estimate protein-group FDR without the inflation that the reused PSM formula causes on large data.

**Approach:** For each target group, find its decoy counterpart (same accessions with the decoy prefix); keep only the higher-scoring member of each target/decoy PAIR; rank the picked set and count decoys as the FDR estimate. This is the operation the reference example demonstrates end to end.

```python
def picked_group_fdr(groups, decoy_prefix='DECOY_'):
    # groups: list of dicts with 'accessions', 'score', 'is_decoy'
    by_base = {}
    for g in groups:
        base = frozenset(a.replace(decoy_prefix, '') for a in g['accessions'])
        # keep only the higher-scoring of the target/decoy pair (the 'pick')
        if base not in by_base or g['score'] > by_base[base]['score']:
            by_base[base] = g
    picked = sorted(by_base.values(), key=lambda g: g['score'], reverse=True)

    targets = decoys = 0
    for g in picked:
        if g['is_decoy']:
            decoys += 1
        else:
            targets += 1
        g['fdr'] = decoys / targets if targets else 1.0
    running_min = 1.0
    for g in reversed(picked):  # monotone q-values from the bottom up
        running_min = min(running_min, g['fdr'])
        g['qvalue'] = running_min
    return [g for g in picked if not g['is_decoy'] and g['qvalue'] <= 0.01]
```

## Per-Method Failure Modes

### Naive (non-picked) protein/group FDR
**Trigger:** Reusing the PSM-level `decoys/targets` formula at the protein level on a deep dataset.
**Mechanism:** False target proteins (one-hit-wonders) and decoy proteins are not symmetric once peptides are mapped to proteins; true proteins absorb many peptides, false ones do not.
**Symptom:** Reported 1% protein FDR, actual 10-30%; reviewer or entrapment check exposes it.
**Fix:** Picked-protein FDR (Savitski 2015) or picked-group FDR (The & Kall 2016); validate with a two-species or entrapment search.

### Two-peptide rule
**Trigger:** Filtering to proteins with >=2 (unique) peptides "for confidence".
**Mechanism:** The rule removes more target proteins than decoy proteins, inverting the FDR effect, and deletes real low-abundance single-peptide proteins.
**Symptom:** Fewer proteins AND higher true FDR than picked FDR at the same nominal cutoff.
**Fix:** Drop the rule; control picked protein-level FDR and score single-peptide IDs individually.

### Razor-peptide quantification
**Trigger:** Quantifying on MaxQuant's default unique+razor peptides for a sensitive comparison.
**Mechanism:** A shared peptide's full intensity is credited to one group; that razor assignment can flip between conditions when peptide counts shift, so a protein's quantity changes for inference reasons, not biology.
**Symptom:** Spurious differential abundance concentrated on proteins sharing peptides with paralogs.
**Fix:** Quantify on unique peptides only for sensitive comparisons -> quantification.

### Parsimony tie-breaking
**Trigger:** Multiple minimal protein sets explain the peptides equally well.
**Mechanism:** Greedy parsimony breaks ties arbitrarily; minimality is a heuristic, not truth, and a real protein with only shared peptides is silently dropped.
**Symptom:** Reported lead protein differs run-to-run or pipeline-to-pipeline on the same data.
**Fix:** Prefer a probabilistic method (EPIFANY/ProteinProphet) that apportions shared evidence; retain group membership.

### Proteoform overreach
**Trigger:** Reporting "isoform X is present" from a group whose evidence is shared peptides.
**Mechanism:** Splice isoforms, variants, and PTM forms collapse into groups in bottom-up data; the group cannot separate them.
**Symptom:** Isoform-specific claim with no isoform-unique peptide behind it.
**Fix:** Require an isoform-unique peptide for any isoform claim, or use top-down / proteoform methods.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Protein / protein-group FDR 1% (sometimes 5% for discovery) | community standard | SEPARATE estimation from PSM FDR; never assume 1% PSM implies 1% protein |
| Picked FDR (target/decoy pairing) | Savitski 2015; The & Kall 2016 | Removes target/decoy asymmetry; dataset-size-independent, unlike naive decoy/target |
| Decoy:target ratio 1:1 | community standard | Standard null; unequal ratios require formula correction |
| Min PSMs for trustworthy protein FDR | hundreds+ | Below ~100s of items decoy counts are too noisy; "0% FDR" from zero decoys is luck, not control |
| Two-peptide rule | DO NOT USE (Gupta & Pevzner 2009) | Increases protein FDR and drops real proteins; replaced by picked FDR + per-ID score |
| Single-peptide IDs | judge by score, not count | A high-confidence unique peptide can be a legitimate ID |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Protein FDR much higher than nominal on deep data | Naive decoy/target reused from PSM level | Picked-protein or picked-group FDR |
| Real low-abundance proteins missing | Two-peptide rule applied | Remove the rule; control picked FDR |
| AttributeError on `EpifanyAlgorithm` / `infer_proteins` | Class name varies by version; the R `ProteinInference::infer_proteins` could not be confirmed to exist | `help(pyopenms.EpifanyAlgorithm)` to find the real name; use pyOpenMS, not an unverified R package |
| Indistinguishable proteins reported as separate IDs | Flat protein list instead of groups | Enable `annotate_indistinguishable_groups`; report groups with a leading protein |
| Spurious DE on paralog-sharing proteins | Razor-peptide quant flipped between conditions | Quantify on unique peptides -> quantification |
| "Unique" peptide count changed when DB changed | Uniqueness is database-relative | Fix and document the database (isoforms, contaminants, decoys) |

## References

- Nesvizhskii, A.I., Keller, A., Kolker, E. & Aebersold, R. (2003). A statistical model for identifying proteins by tandem mass spectrometry. *Analytical Chemistry* 75(17):4646-4658.
- Gupta, N. & Pevzner, P.A. (2009). False discovery rates of protein identifications: a strike against the two-peptide rule. *Journal of Proteome Research* 8(9):4173-4181.
- Savitski, M.M., Wilhelm, M., Hahne, H., Kuster, B. & Bantscheff, M. (2015). A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Molecular & Cellular Proteomics* 14(9):2394-2404.
- The, M., Tasnim, A. & Kall, L. (2016). How to talk about protein-level false discovery rates in shotgun proteomics. *Proteomics* 16(18):2461-2469.
- Pfeuffer, J., Sachsenberg, T., Dijkstra, T.M.H., Serang, O., Reinert, K. & Kohlbacher, O. (2020). EPIFANY: a method for efficient high-confidence protein inference. *Journal of Proteome Research* 19(3):1060-1072.

## Related Skills

- peptide-identification - Produces the FDR-filtered peptide list that feeds inference and shares the target-decoy machinery
- quantification - Consumes inferred groups; razor-vs-unique peptide choice lives here
- data-import - Loads idXML/mzML identification files
- database-access/uniprot-access - Canonical-vs-isoform databases drive uniqueness and the leading-protein convention
