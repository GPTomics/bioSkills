---
name: bio-reaction-enumeration
description: Enumerates virtual chemical libraries via reaction SMARTS transformations using RDKit and Reaction templates, with explicit handling of atom mapping, template extraction (RDKit reaction mining), product validation, RECAP/BRICS fragmentation, R-group decomposition, matched molecular pair analysis (MMPA), and Free-Wilson analysis. Use when generating combinatorial libraries from building blocks, enumerating analog series, deriving structure-activity rules, or extracting transformations from reaction data.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, mmpdb 3.1+, scikit-learn 1.4+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Reaction Enumeration

Generate virtual libraries by applying reaction SMARTS to building blocks, enumerate analog series via matched molecular pairs, decompose into R-groups for SAR modeling, or extract transformations from reaction data. Reaction enumeration sits at the intersection of medicinal chemistry, lead optimization, and de novo design. The two key operations: **transform** (apply known rxn to make new compounds) and **mine** (extract rules from observed analog series). RDKit's reaction SMARTS handles the former; mmpdb / Free-Wilson handle the latter.

For retrosynthetic planning (target-to-starting-material decomposition), see `chemoinformatics/retrosynthesis`. For ML-driven design, see `chemoinformatics/generative-design`. For scaffold-based design, see `chemoinformatics/scaffold-analysis`.

## Operation Taxonomy

| Operation | Goal | Tool | Fails when |
|-----------|------|------|------------|
| Forward enumeration | Apply reaction to building blocks -> products | RDKit `ReactionFromSmarts` + `RunReactants` | Wrong atom mapping; missing connectivity |
| Reverse enumeration (retrosynthesis) | Product -> starting materials | AiZynthFinder, Chemformer | See retrosynthesis skill |
| Template mining | Reaction database -> reaction SMARTS templates | RDKit reaction mining; rxnmapper | Atom mapping ambiguous; mechanism unclear |
| RECAP fragmentation | Molecule -> retro-synthetic fragments | RDKit `Chem.Recap` | Inflexible bond rules |
| BRICS fragmentation | Molecule -> retro-synthetic fragments | RDKit `BRICS` module | Many false fragments |
| R-group decomposition | Set of mols + scaffold -> R-group table | RDKit `Chem.rdRGroupDecomposition` | Multiple scaffolds; ambiguous attachment |
| Matched Molecular Pairs (MMPA) | Set of mols -> transformation rules | mmpdb | Need ≥1k compound dataset |
| Free-Wilson | Compounds + activities -> additive R-group contributions | scikit-learn linear regression | Strict additivity assumption |

## Reaction SMARTS Basics

A reaction SMARTS is `reactants >> products` with atom maps `[atom:idx]` tracking atoms through the transformation:

```python
from rdkit.Chem import AllChem, Chem

amide = AllChem.ReactionFromSmarts(
    '[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]'
)

errors = amide.Validate()
print(errors)
```

**Atom mapping rules:**
- Atoms with the same map index `[C:1]` in both reactant and product are tracked
- Maps must be unique within each reactant/product
- Unmapped atoms are added to or removed from the product
- Bond orders may change; map index preserves identity

**Common error:** Leaving an atom unmapped causes RDKit to either lose or duplicate it.

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

These are templates; real reactions need stereo, protecting-group, and chemoselectivity considerations. For production library enumeration, use validated templates from `rxnmapper` or vendor catalogs.

## Combinatorial Library Enumeration

**Goal:** Generate every (R1, R2, ..., Rn) product combination from sets of building blocks.

**Approach:** Cartesian product of reactant lists; apply reaction SMARTS; sanitize + deduplicate.

```python
from itertools import product
from rdkit import Chem
from rdkit.Chem import AllChem

def enumerate_library(rxn_smarts, reactant_lists, mw_max=600):
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    if rxn.Validate()[0] != 0:
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
                    if Chem.Descriptors.MolWt(prod) > mw_max:
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

RECAP bond types: amide, ester, ether, amine, urea, olefin, quaternary nitrogen, sulfonamide. Use cases: building-block library generation, scaffold-decoration enumeration.

## BRICS Fragmentation

BRICS (Degen 2008) is an extension of RECAP with more bond types. Better fragment coverage; more fragments per molecule.

```python
from rdkit.Chem import BRICS

mol = Chem.MolFromSmiles('CCN(CC)c1ccc(C(=O)NC2CCCC2)cc1')
fragments = BRICS.BRICSDecompose(mol)

builder = BRICS.BRICSBuild([Chem.MolFromSmiles(f) for f in fragments])
new_mols = [next(builder) for _ in range(10)]
```

`BRICSDecompose` produces SMILES with `[<dummy>]` attachment points; `BRICSBuild` recombines fragments at these dummies.

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

`decomp` is a list of dicts `{'Core': scaffold_smi, 'R1': r1_smi, 'R2': r2_smi}`. Combined with activity column, enables Free-Wilson.

## Matched Molecular Pairs Analysis (MMPA)

MMPA (Hussain & Rea 2010) extracts SAR rules from compound pairs differing by a single transformation.

```bash
mmpdb fragment data.smi -o data.fragments
mmpdb index data.fragments -o data.mmpdb
mmpdb transform --smiles 'COc1ccccc1' data.mmpdb
```

`mmpdb` produces a database of transformations + statistics on activity changes.

| Transformation | Avg delta(pIC50) | N pairs | Confidence |
|----------------|-------------------|---------|------------|
| Me -> F | +0.5 | 152 | high |
| OMe -> OH | -0.3 | 89 | moderate |
| Ph -> 4-pyridine | +1.2 | 23 | moderate |

**Use case:** Lead optimization. Given a hit, ask "what transformations have improved similar series?" Apply top-ranked transformations to generate analog suggestions.

**Context-based MMPA** (Awale 2024): condition rules on local chemical context (e.g., "Me->F adjacent to amide"). Outperforms classical MMPA on CYP1A2 inhibition reduction.

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

**Approach:** Use `rxnmapper` (Schwaller 2021) for atom mapping, then RDKit reaction template extraction.

```python
from rxnmapper import RXNMapper

mapper = RXNMapper()
rxns = ['CCO.OC(=O)c1ccccc1>>CCOC(=O)c1ccccc1']
results = mapper.get_attention_guided_atom_maps(rxns)
mapped_smiles = results[0]['mapped_rxn']
```

After atom-mapping, RDKit can extract a template via `ChemicalReaction.GetReactionTemplateFromMappedReaction` (custom implementation; see Coley 2019).

## Per-Tool Failure Modes

### Reaction SMARTS -- atom mapping mismatch

**Trigger:** Map index appears on reactant but not product.

**Mechanism:** RDKit treats unmapped atoms as deleted from product; an atom that was meant to be preserved disappears if its map index is missing.

**Symptom:** Products missing expected atoms; valences wrong; sanitize fails.

**Fix:** Validate with `rxn.Validate()`; manually inspect mapping; use Reaction Atom Mapping Number column 2.

### RECAP/BRICS -- over-fragmentation

**Trigger:** Highly substituted molecule with many breakable bonds.

**Mechanism:** Default bond list breaks at every retrosynthetic position; one molecule yields tens of fragments.

**Symptom:** Building-block enumeration explodes; many small irrelevant fragments.

**Fix:** Filter fragments by MW (>=80 Da), heavy atom count (>=4); use only meaningful fragments downstream.

### MMPA -- insufficient pair count

**Trigger:** mmpdb on small dataset (<500 compounds, <50 actives).

**Mechanism:** MMPA needs ≥10 pairs per transformation to yield a statistically meaningful delta(activity).

**Symptom:** Transformations report with N=1-3 pairs; effect sizes erratic.

**Fix:** Filter to transformations with N>=10; supplement with literature SAR knowledge.

### Free-Wilson -- non-additive interactions

**Trigger:** R1 and R2 interact through hydrogen bonding, steric clash, or electronic effects.

**Mechanism:** Free-Wilson is purely additive; cannot capture R1+R2 synergy.

**Symptom:** Predicted activities for un-synthesized combinations are biased low for synergistic pairs.

**Fix:** Use Free-Wilson as first-pass screen; validate predictions with QSAR (random forest, chemprop) which captures interactions.

### R-group decomposition -- multiple scaffolds

**Trigger:** Compound matches multiple scaffold templates.

**Mechanism:** `RGroupDecompose` returns the first matching scaffold; ambiguous SAR series.

**Symptom:** Same compound's R-groups differ between runs.

**Fix:** Specify scaffold unambiguously; use only compounds matching one scaffold.

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
| `rxn.Validate()` returns errors | Bad atom mapping or invalid SMARTS | Re-check map indices; valences |
| Products contain unexpected fragments | Reactants matched in unintended way | Use more specific SMARTS; constrain with explicit ring members |
| Sanitize fails on products | Reaction breaks valence | Filter via `Chem.SanitizeMol(prod, catchErrors=True)` |
| Duplicate products | Same product from different reactant orientations | Deduplicate by canonical SMILES |
| RECAP produces single fragment | Molecule has no retrosynthetic bonds | Try BRICS for more aggressive fragmentation |
| mmpdb empty output | Insufficient dataset size or no matched pairs | Need >=1000 compounds |
| R-group decomposition wrong R | Scaffold dummy not aligned | Re-check `[*:1]` / `[*:2]` placement |

## References

- Hartenfeller et al., *J. Cheminformatics* 4:38 (2012) -- DOGS rule-based library design.
- Lewell et al., *J. Chem. Inf. Comput. Sci.* 38:511 (1998) -- RECAP.
- Degen et al., *ChemMedChem* 3:1503 (2008) -- BRICS fragmentation.
- Hussain & Rea, *J. Chem. Inf. Model.* 50:339 (2010) -- MMPA.
- Dossetter et al., *Drug Discov. Today* 18:724 (2013) -- Practical MMPA in lead optimization.
- Free & Wilson, *J. Med. Chem.* 7:395 (1964) -- Original Free-Wilson.
- Schwaller et al., *Sci. Adv.* 7:eabe4166 (2021) -- rxnmapper.
- Coley et al., *Chem. Sci.* 10:370 (2019) -- Reaction template extraction.

## Related Skills

- chemoinformatics/molecular-io - Read/write reaction SMILES
- chemoinformatics/substructure-search - SMARTS pattern matching
- chemoinformatics/scaffold-analysis - Bemis-Murcko scaffolds for R-decomp
- chemoinformatics/molecular-descriptors - Featurize products
- chemoinformatics/admet-prediction - Filter enumerated products
- chemoinformatics/retrosynthesis - Reverse direction (target -> starting materials)
- chemoinformatics/generative-design - Generative alternatives to template enumeration
- chemoinformatics/qsar-modeling - Validate Free-Wilson predictions
