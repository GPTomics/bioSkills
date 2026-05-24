---
name: bio-retrosynthesis
description: Performs retrosynthetic planning using AiZynthFinder (MCTS, template-based), Chemformer (template-free transformer), ASKCOS, and emerging RetroSynFormer with explicit handling of route scoring, building-block availability (eMolecules, Enamine, Mcule), forward prediction validation (Molecular Transformer), and disconnection-aware multi-objective search (MO-MCTS). Use when assessing synthetic feasibility of generated or selected molecules, planning multi-step syntheses, building synthesis-aware design pipelines, or screening libraries for retro-route feasibility.
tool_type: python
primary_tool: AiZynthFinder
---

## Version Compatibility

Reference examples tested with: AiZynthFinder 4.4+, Chemformer 1.3+, RDKit 2024.09+, RDChiral 1.1+, Aizynthtrain 1.0+, ASKCOS Lite 0.5+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `aizynthcli --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Retrosynthesis

Plan synthetic routes from a target molecule back to commercially-available building blocks. AiZynthFinder 4.0 (Saigiridharan, Hassen, Lai, Torren-Peraire, Engkvist, Genheden 2024 *J Cheminform* 16:57, AstraZeneca) is the open-source production-grade tool: Monte Carlo Tree Search (MCTS) + template-based expansion + multi-objective scoring (MO-MCTS). Chemformer (Irwin 2022) is template-free transformer alternative. ASKCOS (MIT) is the academic reference. Modern best practice combines retrosynthesis with **forward validation** (the predicted route should also predict the target from starting materials via Molecular Transformer) and building-block availability (eMolecules, Enamine, Mcule, ZINC catalog).

For generative design pipelines that need synthetic feasibility, see `chemoinformatics/generative-design`. For reaction enumeration (forward direction), see `chemoinformatics/reaction-enumeration`.

## Retrosynthesis Method Taxonomy

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| AiZynthFinder 4.0 | Template-based MCTS | Open, scalable, well-validated | Beyond template coverage |
| Chemformer | Template-free transformer | Novel disconnections | Less interpretable; harder to debug |
| ASKCOS | Template-based + neural | MIT-quality academic standard | Setup complexity |
| Molecular Transformer | Forward + retro transformer | Single SMILES-to-SMILES | Less robust to non-training distribution |
| RetroSynFormer | Decision transformer | Modern method | Limited adoption |
| IBM RXN | Cloud service | High quality, easy interface | API access required |
| BKMS_MTHRO / RetroPath | Pathway-based | Metabolic / biochemical | Not for general medchem |
| SyntheMol (StanfordGarbage) | Specialized for medchem | Public domain alternative | Limited tooling |

**Decision:** For most users, **AiZynthFinder 4.4 with USPTO + USPTO-50k templates** is the open-source standard. For high-stakes routes, validate with Molecular Transformer forward prediction.

## Decision Tree by Scenario

| Scenario | Tool | Notes |
|----------|------|-------|
| Standard medchem target | AiZynthFinder default templates | USPTO/Reaxys templates |
| Novel chemotype | AiZynthFinder + Chemformer template-free fallback | Combine both |
| Generated molecules (REINVENT output) | AiZynthFinder batch | Filter to feasible routes |
| Multi-step synthesis planning | AiZynthFinder + manual review | Top-K routes |
| Validate generated route | Molecular Transformer forward | Check round-trip |
| Cost-aware synthesis | AiZynthFinder + custom building-block pricing | Score weight |
| Disconnection-aware design (DAD) | AiZynthFinder MO-MCTS | Multi-objective |
| Patent-aware routes | Custom template exclusion | Specialized |

## AiZynthFinder Setup

**Goal:** Configure AiZynthFinder with USPTO templates + a building-block stock and run MCTS retrosynthesis planning on a target SMILES.

**Approach:** Build a configuration dict pointing to policy templates (ONNX + CSV) and a stock H5, instantiate `AiZynthFinder`, set the target SMILES, then call `tree_search()` followed by `build_routes()`.

```python
from aizynthfinder.aizynthfinder import AiZynthFinder

config_dict = {
    'policy': {
        'files': {
            'uspto': ['policy/uspto_model.onnx', 'templates/uspto_templates.csv'],
        }
    },
    'stock': {
        'files': {
            'zinc': 'stock/zinc.h5',
        }
    },
    'finder': {
        'algorithm': 'mcts',
        'iteration_limit': 100,
        'time_limit': 120,
    }
}

finder = AiZynthFinder(configdict=config_dict)
finder.target_smiles = 'CC(=O)Nc1ccc(C(=O)Nc2cccc(C(F)(F)F)c2)cc1'
finder.tree_search()
finder.build_routes()
```

Output: list of routes, each with depth, building blocks, score, leaf nodes.

## Route Output Analysis

```python
for route in finder.routes:
    print(f'Depth: {route.depth}, Score: {route.score:.2f}')
    print(f'In-stock: {sum(node.in_stock for node in route.leafs())}')
    print(f'Building blocks: {[node.smiles for node in route.leafs()]}')
```

Critical metrics:
- **Depth**: number of synthetic steps. 1-3 typical for medchem.
- **Score**: AiZynthFinder route score (0-1, higher = better)
- **In-stock**: how many leaf nodes are commercially available
- **Stock origin**: ZINC, Enamine, Mcule, eMolecules

## Route Scoring (MO-MCTS)

AiZynthFinder 4.0 supports multi-objective scoring (Saigiridharan et al 2024 *J Cheminform* 16:57):

```python
config_dict['finder']['algorithm'] = 'mo_mcts'
config_dict['finder']['mo_mcts'] = {
    'objectives': [
        {'name': 'state_score', 'weight': 0.5},  # default state score
        {'name': 'broken_bonds_score', 'weight': 0.3},  # complexity reduction
        {'name': 'route_length', 'weight': 0.2, 'maximize': False},  # shorter
    ]
}
```

State score: probability the current state can be solved. Broken bonds: each step should reduce molecular complexity. Route length: shorter is better.

## Building Block Stocks

| Stock | Compounds | Source | Cost-tier |
|-------|-----------|--------|-----------|
| ZINC clean leads | 250k | ZINC22 catalog | Various commercial |
| Enamine Building Blocks | 200k+ | Enamine | $$ |
| Enamine REAL | 29B (make-on-demand) | Enamine | $$$ |
| Mcule | 25M | Mcule | $$ |
| eMolecules | 16M | eMolecules | $$ |
| ChemBridge | 1M | ChemBridge | $$ |

AiZynthFinder accepts stocks as HDF5 (built via `aizynthtrain`):
```bash
aizynthtrain build-stock --input zinc_building_blocks.smi --output zinc.h5
```

## Forward Validation with Molecular Transformer

AiZynthFinder predicts retrosynthesis (target → precursors); Molecular Transformer predicts forward (precursors → product). Validating the round-trip:

```python
from molecular_transformer import predict_forward

precursors = route.leafs()  # building blocks from retro
predicted_product = predict_forward(precursors)
match = (Chem.CanonSmiles(predicted_product) == 
         Chem.CanonSmiles(finder.target_smiles))
```

Routes where the forward prediction reproduces the target are highest confidence. ~30-50% of AiZynthFinder routes pass forward validation (Saigiridharan, Genheden et al 2024 *J Cheminform* 16:57).

## Template-Free with Chemformer

Chemformer uses a Transformer (BART) trained on USPTO reactions for SMILES-to-SMILES:

```python
from chemformer import Chemformer

cf = Chemformer.load_pretrained('USPTO_RETROSYNTHESIS_TEMPLATE_FREE')
predictions = cf.predict('CC(=O)Nc1ccc(C(=O)Nc2cccc(C(F)(F)F)c2)cc1',
                         beam_search=10)
```

Output: 10 predicted precursor SMILES. No templates required; can predict novel disconnections.

**Trade-off:** Template-free is more flexible but harder to debug. Combining with AiZynthFinder template MCTS gives best of both.

## Disconnection-Aware Design (DAD)

Modify generative design to also score retrosynthetic feasibility. AiZynthFinder batch mode for 1000+ molecules.

**Goal:** Add retrosynthetic feasibility scoring to generative design pipelines for hundreds-to-thousands of candidate molecules.

**Approach:** Batch-process generated SMILES through `aizynthcli`, classify each compound by route depth and in-stock leaf count, and feed feasibility back into the generative scoring function.

```bash
aizynthcli --smiles compounds.smi --output routes.json \
           --config config.yaml --policy uspto --stock zinc
```

For each compound, returns top-K routes. Score-feasibility for generative design:
- "Synthesizable" = in-stock leaves >= 2 in best route
- "Routable" = at least one route depth <= 5
- "Easy" = at least one route depth <= 3 with all leaves in-stock

## Cost-Aware Synthesis

Add building-block pricing as objective:

```python
def route_cost(route, price_db):
    total = 0
    for leaf in route.leafs():
        smi = Chem.CanonSmiles(leaf.smiles)
        if smi in price_db:
            total += price_db[smi]
    return total
```

Combine with step-cost estimate (typical: ~$500-2000 per synthesis step).

## Per-Tool Failure Modes

### AiZynthFinder -- template coverage gap

**Trigger:** Target molecule uses bond formation not in training reactions.

**Mechanism:** USPTO templates are biased toward common transformations; novel chemistry (organometallics, exotic heterocycles) missing.

**Symptom:** No solved route or route uses unsuitable simplifications.

**Fix:** Augment templates from Reaxys; combine with Chemformer; manual review.

### Chemformer -- non-canonical SMILES output

**Trigger:** Default Chemformer output.

**Mechanism:** Transformer can produce non-canonical SMILES variants.

**Symptom:** SMILES round-trip fails; validation tools confused.

**Fix:** Canonicalize Chemformer output via RDKit before comparing.

### Route uses non-stock building block

**Trigger:** Leaf node not in stock database.

**Mechanism:** AiZynthFinder tree may end on non-purchasable molecules.

**Symptom:** Route "complete" but route has non-stock leaves.

**Fix:** Filter routes by `in_stock_only=True`; or expand stock to include theoretical building blocks (Enamine REAL).

### MCTS iteration limit too low

**Trigger:** Complex target requiring deep tree search.

**Mechanism:** MCTS may not find route in default 100 iterations.

**Symptom:** No routes returned despite plausible target.

**Fix:** Increase `iteration_limit=500` and `time_limit=600`; consider divide-and-conquer for complex targets.

### Forward validation fails

**Trigger:** Retro route uses chemistry that doesn't actually work in forward.

**Mechanism:** Template-based retro lacks reaction conditions / catalysts; forward prediction more conservative.

**Symptom:** Forward predicts different product than target.

**Fix:** Use as confidence signal, not rejection; many routes don't round-trip but are still valid synthesis-wise.

### Building block stock obsolete

**Trigger:** Old ZINC catalog used; building blocks no longer purchasable.

**Mechanism:** Commercial catalogs change quarterly.

**Symptom:** Routes recommend unavailable building blocks.

**Fix:** Use Enamine REAL or recent ZINC22 for current stock; verify with vendor before synthesis.

## Reconciliation: AiZynthFinder vs Chemformer

| Aspect | AiZynthFinder | Chemformer |
|--------|---------------|------------|
| Approach | Templates + MCTS | Transformer encoder-decoder |
| Speed | Fast for shallow trees | Single-pass per target |
| Interpretability | High (template + atom mapping) | Low (black box) |
| Novel disconnections | Limited | Better |
| Production maturity | High | Medium |
| Cost | CPU | GPU recommended |

For comprehensive coverage, run both and merge unique routes.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `tree_search()` returns no routes | Target outside template coverage | Increase iterations; try Chemformer |
| All routes depth > 8 | Complex target | Likely correct; review manually |
| Route says "solved" but leaves not in stock | Stock incomplete | Update stock; or set `in_stock_only=True` |
| Building block price not found | Compound not in pricing DB | Use Enamine quote or vendor inquiry |
| Chemformer truncates SMILES | Token limit | Increase max_length |
| Forward prediction wrong | Out-of-distribution reaction | Use as confidence signal only |
| MCTS slow on simple target | Default config | Reduce time_limit; use smaller template set |

## References

- Saigiridharan L, Hassen AK, Lai J, Torren-Peraire P, Engkvist O, Genheden S 2024 *J Cheminform* 16:57 -- AiZynthFinder 4.0 (multi-objective MCTS retrosynthesis).
- Irwin et al., *Mach. Learn.: Sci. Technol.* 3:015022 (2022) -- Chemformer.
- Schwaller et al., *ACS Cent. Sci.* 5:1572 (2019) -- Molecular Transformer.
- Coley et al., *ACS Cent. Sci.* 3:434 (2017) -- ASKCOS template extraction.
- *Digital Discovery* (2026) -- RetroSynFormer decision transformer.

## Related Skills

- chemoinformatics/molecular-io - Parse target and route SMILES
- chemoinformatics/molecular-standardization - Standardize before retrosynthesis
- chemoinformatics/generative-design - Add synthetic feasibility to scoring
- chemoinformatics/reaction-enumeration - Forward direction (template enumeration)
- chemoinformatics/admet-prediction - Filter targets before retrosynthesis
