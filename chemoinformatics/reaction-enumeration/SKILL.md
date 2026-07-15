---
name: bio-reaction-enumeration
description: Enumerates virtual chemical libraries via reaction SMARTS transformations using RDKit and reaction templates, with explicit handling of atom mapping, RDChiral template extraction, product validation, RECAP/BRICS fragmentation, R-group decomposition, matched molecular pair analysis (MMPA), and Free-Wilson analysis. Use when generating combinatorial libraries from building blocks, enumerating analog series, deriving structure-activity rules, or extracting transformations from reaction data.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, RDChiral 1.1+, mmpdb 3.1+, scikit-learn 1.4+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Reaction Enumeration

Generate virtual libraries by applying reaction SMARTS to building blocks, enumerate analog series via matched molecular pairs, decompose into R-groups for SAR modeling, or extract transformations from reaction data. Reaction enumeration sits at the intersection of medicinal chemistry, lead optimization, and de novo design; DOGS is one published example of reaction-driven de novo design (Hartenfeller et al. 2012). The two key operations are **transform** (apply known reactions to make new compounds) and **mine** (extract rules from observed analog series). RDKit's reaction SMARTS handles the former; mmpdb / Free-Wilson handle analog-series analysis, while mapped-reaction template extraction requires separate tooling such as RDChiral.

For retrosynthetic planning (target-to-starting-material decomposition), see `chemoinformatics/retrosynthesis`. For ML-driven design, see `chemoinformatics/generative-design`. For scaffold-based design, see `chemoinformatics/scaffold-analysis`.

## Operation Taxonomy

| Operation | Goal | Tool | Fails when |
|-----------|------|------|------------|
| Forward enumeration | Apply reaction to building blocks -> products | RDKit `ReactionFromSmarts` + `RunReactants` | Wrong atom mapping; missing connectivity |
| Reverse enumeration (retrosynthesis) | Product -> starting materials | AiZynthFinder, Chemformer | See retrosynthesis skill |
| Template mining | Reaction database -> reaction SMARTS templates | RXNMapper + RDChiral | Atom mapping ambiguous; mechanism unclear |
| RECAP fragmentation | Molecule -> retro-synthetic fragments | RDKit `Chem.Recap` | Inflexible bond rules |
| BRICS fragmentation | Molecule -> retro-synthetic fragments | RDKit `BRICS` module | Many false fragments |
| R-group decomposition | Set of mols + scaffold -> R-group table | RDKit `Chem.rdRGroupDecomposition` | Multiple scaffolds; ambiguous attachment |
| Matched Molecular Pairs (MMPA) | Set of mols -> transformation rules | mmpdb | Sparse or context-confounded matched pairs |
| Free-Wilson | Compounds + activities -> additive R-group contributions | scikit-learn linear regression | Strict additivity assumption |

## Reaction SMARTS Basics

A reaction SMARTS is `reactants >> products` with atom maps `[atom:idx]` tracking atoms through the transformation:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

amide = AllChem.ReactionFromSmarts(
    '[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]'
)

errors = amide.Validate()
print(errors)
```

**Atom mapping rules:**
- Atoms with the same map index `[C:1]` in both reactant and product are tracked
- Maps must be unique within each reactant/product
- Atoms present in a reactant template but absent from the product template are removed; atoms present only in the product template are created
- Atoms outside the matched reaction-center template are generally carried through with their reactant molecule
- Bond orders may change; map index preserves identity

**Common error:** Missing or inconsistent maps for atoms intended to survive within the reaction center can delete atoms, create duplicates, or obscure which reactant atom a product atom represents. Mapping alone does not define the transformation; the reactant and product templates do.

## Common Reaction Templates

```python
REACTIONS = {
    'amide_coupling': '[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]',
    'reductive_amination': '[C:1](=O).[NH2:2]>>[CH:1][NH:2]',
    'suzuki': '[c:1][Br].[c:2][B](O)O>>[c:1][c:2]',
    'buchwald_hartwig': '[c:1][Br].[NH:2]>>[c:1][N:2]',
    'sn2_substitution': '[CH:1][Br].[N:2]>>[CH:1][N:2]',
    'sonogashira': '[c:1][Br].[CH:2]#[C:3]>>[c:1][C:2]#[C:3]',
    'click_chemistry': '[N-:1]=[N+:2]=[N:3][CH2:4].[CH:5]#[C:6]>>[N:3]1[N:2]=[N:1][C:6]=[C:5]1[CH2:4]',
    'esterification': '[C:1](=[O:2])O.[OH:3][C:4]>>[C:1](=[O:2])[O:3][C:4]',
    'urea_formation': '[N:1]=C=O.[NH:2]>>[N:1]C(=O)[N:2]',
    'sulfonamide': '[S:1](=O)(=O)Cl.[NH:2]>>[S:1](=O)(=O)[N:2]',
}
```

These are illustrative templates; real reactions need stereo, protecting-group, and chemoselectivity considerations. For production library enumeration, use a separately curated and validated template collection, such as a vendor catalog. RXNMapper maps atoms in reaction records; it does not by itself supply or validate reaction templates.

## Combinatorial Library Enumeration

**Goal:** Generate every (R1, R2, ..., Rn) product combination from sets of building blocks.

**Approach:** Cartesian product of reactant lists; apply reaction SMARTS; sanitize + deduplicate.

```python
from itertools import product
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def enumerate_library(rxn_smarts, reactant_lists, mw_max=600):
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    num_warnings, num_errors = rxn.Validate()
    if num_errors:
        raise ValueError(f'Invalid reaction: {rxn_smarts}')

    seen = set()
    products = []
    for combo in product(*reactant_lists):
        mols = [Chem.MolFromSmiles(s) for s in combo]
        if None in mols:
            continue

        for prod_tuple in rxn.RunReactants(tuple(mols)):
            for prod in prod_tuple:
                try:
                    Chem.SanitizeMol(prod)
                    smi = Chem.MolToSmiles(prod)
                    if smi in seen:
                        continue
                    if Descriptors.MolWt(prod) > mw_max:
                        continue
                    seen.add(smi)
                    products.append(smi)
                except Exception:
                    continue
    return products
```

**Scaling:** For a 100x100x100 enumeration (1M products), parallelize with multiprocessing. For 1k x 1k x 1k (1B products), use a streaming approach + filter before materializing.

## RECAP Fragmentation

RECAP (Lewell 1998) breaks molecules at retrosynthetically reasonable bonds into reusable fragments.

```python
from rdkit.Chem import Recap

mol = Chem.MolFromSmiles('c1ccc(C(=O)Nc2ccc(F)cc2)cc1')
hier = Recap.RecapDecompose(mol)
fragments = list(hier.GetLeaves().keys())
```

RDKit's current `Recap.reactionDefs` contains 12 cleavage definitions. Inspect the installed definitions when exact coverage matters instead of relying on a shortened functional-group list. Use cases include building-block library generation and scaffold-decoration enumeration.

## BRICS Fragmentation

BRICS (Degen 2008) is an extension of RECAP with more bond types. Better fragment coverage; more fragments per molecule.

```python
from itertools import islice
from rdkit.Chem import BRICS

mol = Chem.MolFromSmiles('CCN(CC)c1ccc(C(=O)NC2CCCC2)cc1')
fragments = BRICS.BRICSDecompose(mol)

builder = BRICS.BRICSBuild([Chem.MolFromSmiles(f) for f in fragments])
new_mols = list(islice(builder, 10))
```

`BRICSDecompose` produces SMILES with isotope/environment-labeled dummy atoms such as `[1*]`, `[5*]`, and `[16*]`; `BRICSBuild` uses those labels when recombining compatible fragments.

## R-Group Decomposition

**Goal:** Given a set of compounds sharing a scaffold, extract the R-group at each attachment point into a tabular SAR matrix.

**Approach:** Define scaffold with `[*:1]`, `[*:2]` placeholders; RDKit matches each compound and extracts R-groups.

```python
from rdkit.Chem import rdRGroupDecomposition as rgd
from rdkit import Chem

scaffold = Chem.MolFromSmiles('c1ccc(-[*:1])cc1-[*:2]')

mols = [Chem.MolFromSmiles(smi) for smi in [
    'c1ccc(C)cc1F',
    'c1ccc(CC)cc1Cl',
    'c1ccc(CCC)cc1Br',
]]

decomp, _ = rgd.RGroupDecompose([scaffold], mols, asSmiles=True)
```

`decomp` is a list of dicts such as `{'Core': scaffold_smi, 'R1': r1_smi, 'R2': r2_smi}`. It does not contain assay values. Preserve compound identifiers and explicitly join the decomposition to the activity table before Free-Wilson analysis:

```python
import pandas as pd

compound_ids = ['cmpd-1', 'cmpd-2', 'cmpd-3']
activities = pd.DataFrame({
    'compound_id': compound_ids,
    'pIC50': [6.2, 6.8, 7.1],  # example measurements
})
decomp, unmatched = rgd.RGroupDecompose([scaffold], mols, asSmiles=True)
unmatched = set(unmatched)
matched_ids = [cid for i, cid in enumerate(compound_ids) if i not in unmatched]
decomp_df = pd.DataFrame(decomp)
decomp_df.insert(0, 'compound_id', matched_ids)
sar_table = decomp_df.merge(
    activities, on='compound_id', how='inner', validate='one_to_one'
)
```

## Matched Molecular Pairs Analysis (MMPA)

MMPA (Hussain & Rea 2010) extracts SAR rules from compound pairs differing by a single transformation.

```bash
mmpdb fragment data.smi -o data.fragments
mmpdb index data.fragments -o data.mmpdb
mmpdb transform --smiles 'COc1ccccc1' data.mmpdb
```

`mmpdb` produces a database of transformations + statistics on activity changes. The values below are synthetic examples showing the output schema; they are not observations from a cited dataset.

| Transformation | Avg delta(pIC50) | N pairs | Confidence |
|----------------|-------------------|---------|------------|
| Me -> F | +0.5 | 152 | high |
| OMe -> OH | -0.3 | 89 | moderate |
| Ph -> 4-pyridine | +1.2 | 23 | moderate |

**Use case:** Lead optimization. Given a hit, ask "what transformations have improved similar series?" Apply top-ranked transformations to generate analog suggestions.

**Context-based MMPA** conditions transformation statistics on local chemical context (for example, "Me -> F adjacent to an amide"). Raut and Dixit (2025) applied this approach to identify transformations associated with reduced CYP1A2 inhibition; that endpoint-specific result should not be generalized as universal superiority over classical MMPA.

## Free-Wilson Analysis

**Goal:** Decompose activity into additive R-group contributions.

**Approach:** Linear regression with R-group identity as binary features.

```python
import pandas as pd
from sklearn.linear_model import Ridge

def free_wilson(decomp_results, activity_col='pIC50'):
    df = pd.DataFrame(decomp_results)
    r_groups = pd.get_dummies(df[['R1', 'R2']], prefix=['R1', 'R2'])
    X = r_groups.values
    y = df[activity_col].values
    model = Ridge(alpha=0.1).fit(X, y)
    contributions = dict(zip(r_groups.columns, model.coef_))
    return contributions, model.intercept_
```

**Trade-off:** Free-Wilson assumes additivity (R1 contribution independent of R2). Real SAR has interactions; Free-Wilson predictions for un-synthesized combinations are biased when synergy exists. Use as a *first-pass model* for analog prioritization; validate with QSAR.

## Template Extraction from Reaction Data

**Goal:** Given an atom-mapped reaction SMILES, extract a generalizable SMARTS template.

**Approach:** Use `rxnmapper` (Schwaller et al. 2021) for atom mapping, then a template extractor such as RDChiral (Coley et al. 2019).

```python
from rxnmapper import RXNMapper

mapper = RXNMapper()
rxns = ['CCO.OC(=O)c1ccccc1>>CCOC(=O)c1ccccc1']
results = mapper.get_attention_guided_atom_maps(rxns)
mapped_smiles = results[0]['mapped_rxn']
```

RDKit does not provide a `ChemicalReaction.GetReactionTemplateFromMappedReaction` method. Pass the mapped reaction to RDChiral's published template-extraction workflow, checking the installed package's interface and expected reaction-record schema, or use a separately implemented and validated extractor.

## Per-Tool Failure Modes

### Reaction SMARTS -- atom mapping mismatch

**Trigger:** An atom intended to survive the reaction center is absent from, or inconsistently mapped in, the product template.

**Mechanism:** RDKit constructs products from the reaction templates. Reactant-template atoms omitted from the product template are deleted, product-only atoms are created, and inconsistent maps can prevent intended atom identity from being carried across the transformation.

**Symptom:** Products missing expected atoms; valences wrong; sanitize fails.

**Fix:** Validate with `rxn.Validate()`, inspect warnings separately from errors, and manually verify that every reaction-center atom intended to survive has one consistent map number on both sides.

### RECAP/BRICS -- over-fragmentation

**Trigger:** Highly substituted molecule with many breakable bonds.

**Mechanism:** Default bond list breaks at every retrosynthetic position; one molecule yields tens of fragments.

**Symptom:** Building-block enumeration explodes; many small irrelevant fragments.

**Fix:** Filter fragments by MW (>=80 Da), heavy atom count (>=4); use only meaningful fragments downstream.

### MMPA -- insufficient pair count

**Trigger:** The dataset yields few matched pairs for a transformation in the relevant chemical context.

**Mechanism:** Sparse or heterogeneous pairs give imprecise, context-confounded estimates. There is no universal minimum dataset size or pair count that guarantees a meaningful effect.

**Symptom:** Transformations report with N=1-3 pairs; effect sizes erratic.

**Fix:** Report pair counts and uncertainty, examine local contexts, use a project-justified precision threshold, and supplement with experimental or literature SAR knowledge.

### Free-Wilson -- non-additive interactions

**Trigger:** R1 and R2 interact through hydrogen bonding, steric clash, or electronic effects.

**Mechanism:** Free-Wilson is purely additive; cannot capture R1+R2 synergy.

**Symptom:** Predicted activities for un-synthesized combinations are biased low for synergistic pairs.

**Fix:** Use Free-Wilson as first-pass screen; validate predictions with QSAR (random forest, chemprop) which captures interactions.

### R-group decomposition -- multiple scaffolds

**Trigger:** Compound matches multiple scaffold templates.

**Mechanism:** RDKit accepts cores ordered from most to least specific and exposes parameters for multi-core matching and alignment. Ambiguous or inconsistently specified cores can change the resulting labels and SAR table.

**Symptom:** Same compound's R-groups differ between runs.

**Fix:** Order cores from most to least specific, label attachment points explicitly, inspect unmatched compounds, and keep a compound only when its selected core assignment matches the intended SAR series.

### Reaction enumeration -- combinatorial explosion

**Trigger:** Large building-block sets (1k x 1k = 1M products).

**Mechanism:** Cartesian product * RunReactants is O(N^d) where d is reactant count.

**Symptom:** Memory blowup, multi-hour runtime.

**Fix:** Pre-filter building blocks; stream products to file rather than list; use mmpdb-style sparse enumeration only for valid pairings.

## Reconciliation: Free-Wilson vs MMPA

Both methods derive R-group rules but from different perspectives:
- Free-Wilson: linear regression on assembled SAR table; gives R-group contributions
- MMPA: transformation-based; gives delta(activity) for each substitution

If they agree on direction (Me->F improves activity), high confidence. If they disagree, investigate non-additive interactions or look for context dependence in MMPA.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `rxn.Validate()` reports a nonzero error count | Bad atom mapping or invalid SMARTS | Unpack `(num_warnings, num_errors)` and reject on `num_errors`; inspect warnings separately |
| Products contain unexpected fragments | Reactants matched in unintended way | Use more specific SMARTS; constrain with explicit ring members |
| Sanitize fails on products | Reaction breaks valence | Filter via `Chem.SanitizeMol(prod, catchErrors=True)` |
| Duplicate products | Same product from different reactant orientations | Deduplicate by canonical SMILES |
| RECAP produces single fragment | Molecule has no retrosynthetic bonds | Try BRICS for more aggressive fragmentation |
| mmpdb empty output | No pairs satisfy the fragmentation, property, and context criteria | Inspect input parsing and fragmentation output; relax justified filters or obtain relevant analogues |
| R-group decomposition wrong R | Scaffold dummy not aligned | Re-check `[*:1]` / `[*:2]` placement |

## References

- Hartenfeller M et al. "DOGS: Reaction-Driven de novo Design of Bioactive Compounds." *PLoS Comput. Biol.* 8:e1002380 (2012). DOI: 10.1371/journal.pcbi.1002380.
- Lewell XQ, Judd DB, Watson SP, Hann MM. "RECAP—Retrosynthetic Combinatorial Analysis Procedure." *J. Chem. Inf. Comput. Sci.* 38:511–522 (1998). DOI: 10.1021/ci970429i.
- Degen J, Wegscheid-Gerlach C, Zaliani A, Rarey M. "On the Art of Compiling and Using 'Drug-Like' Chemical Fragment Spaces." *ChemMedChem* 3:1503–1507 (2008). DOI: 10.1002/cmdc.200800178.
- Hussain J, Rea C. "Computationally Efficient Algorithm to Identify Matched Molecular Pairs (MMPs) in Large Data Sets." *J. Chem. Inf. Model.* 50:339–348 (2010). DOI: 10.1021/ci900450m.
- Dossetter AG, Griffen EJ, Leach AG. "Matched molecular pair analysis in drug discovery." *Drug Discov. Today* 18:724–731 (2013). DOI: 10.1016/j.drudis.2013.03.003.
- Free SM, Wilson JW. "A Mathematical Contribution to Structure-Activity Studies." *J. Med. Chem.* 7:395–399 (1964). DOI: 10.1021/jm00334a001.
- Schwaller P et al. "Unsupervised attention-guided atom-mapping." *Sci. Adv.* 7:eabe4166 (2021). DOI: 10.1126/sciadv.abe4166.
- Raut JA, Dixit VA. "A context-based matched molecular pair analysis identifies structural transformations that reduce CYP1A2 inhibition." *RSC Med. Chem.* 16:3281–3290 (2025). DOI: 10.1039/D4MD01012D.
- Coley CW, Green WH, Jensen KF. "RDChiral: An RDKit Wrapper for Handling Stereochemistry in Retrosynthetic Template Extraction and Application." *J. Chem. Inf. Model.* 59:2529–2537 (2019). DOI: 10.1021/acs.jcim.9b00286.
- Dalke A, Hert J, Kramer C. "mmpdb: An Open-Source Matched Molecular Pair Platform for Large Multiproperty Data Sets." *J. Chem. Inf. Model.* 58:902–910 (2018). DOI: 10.1021/acs.jcim.8b00173.
- RDKit Book, reaction SMARTS and reaction handling: https://www.rdkit.org/docs/RDKit_Book.html.
- RDKit R-group decomposition API: https://www.rdkit.org/docs/source/rdkit.Chem.rdRGroupDecomposition.html.

## Related Skills

- chemoinformatics/molecular-io - Read/write reaction SMILES
- chemoinformatics/substructure-search - SMARTS pattern matching
- chemoinformatics/scaffold-analysis - Bemis-Murcko scaffolds for R-decomp
- chemoinformatics/molecular-descriptors - Featurize products
- chemoinformatics/admet-prediction - Filter enumerated products
- chemoinformatics/retrosynthesis - Reverse direction (target -> starting materials)
- chemoinformatics/generative-design - Generative alternatives to template enumeration
- chemoinformatics/qsar-modeling - Validate Free-Wilson predictions
