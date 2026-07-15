---
name: bio-retrosynthesis
description: Performs retrosynthetic planning using AiZynthFinder (template-based MCTS), maintained or version-pinned template-free models, ASKCOS, and emerging RetroSynFormer with explicit handling of route scoring, configurable MCTS rewards, building-block availability, and forward-prediction checks. Use when assessing synthetic feasibility of generated or selected molecules, planning multi-step syntheses, building synthesis-aware design pipelines, or screening libraries for retro-route feasibility.
tool_type: python
primary_tool: AiZynthFinder
---

## Version Compatibility

Reference examples tested with: AiZynthFinder 4.4+, RDKit 2024.09+, RDChiral 1.1+, and ASKCOS Lite 0.5+. Chemformer is archived legacy research code; if reproducing it, pin an exact repository commit, checkpoint, and configuration rather than assuming a current pip package or stable Python API.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `aizynthcli --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Retrosynthesis

Plan synthetic routes from a target molecule back to commercially available building blocks. AiZynthFinder combines Monte Carlo Tree Search (MCTS), template-based expansion, configurable search rewards, and route scorers (Saigiridharan et al. 2024). Chemformer is a published template-free transformer baseline whose official repository is now archived; use its exact historical environment for reproduction or select a maintained model with a documented interface. ASKCOS is another open-source synthesis-planning platform. A useful workflow combines retrosynthesis with building-block availability and an independently configured forward-prediction check, while recognizing that a round-trip model match is not experimental validation.

For generative design pipelines that need synthetic feasibility, see `chemoinformatics/generative-design`. For reaction enumeration (forward direction), see `chemoinformatics/reaction-enumeration`.

## Retrosynthesis Method Taxonomy

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| AiZynthFinder 4.4+ | Template-based MCTS | Maintained, configurable open-source planner | Beyond selected policy coverage |
| Chemformer (archived) | Template-free transformer | Reproducing the published baseline | Archived code and checkpoint/config coupling |
| ASKCOS | Template-based + neural | Open-source synthesis-planning platform | Setup complexity |
| Molecular Transformer | Forward + retro transformer | Single SMILES-to-SMILES | Less robust to non-training distribution |
| RetroSynFormer | Decision transformer | Modern method | Limited adoption |

**Decision:** For most users, **AiZynthFinder with its documented public USPTO expansion policy and a current stock** is a practical open-source starting point. For high-stakes routes, apply expert review and an independently configured forward-prediction check; neither model agreement nor a solved search route is experimental validation.

## Decision Tree by Scenario

| Scenario | Tool | Notes |
|----------|------|-------|
| Standard medchem target | AiZynthFinder configured expansion policy | Record the public USPTO policy or licensed template source actually loaded |
| Novel chemotype | AiZynthFinder + maintained or exactly version-pinned template-free comparison | Validate each route independently |
| Generated molecules (REINVENT output) | AiZynthFinder batch | Filter to feasible routes |
| Multi-step synthesis planning | AiZynthFinder + manual review | Top-K routes |
| Validate generated route | Molecular Transformer forward | Check round-trip |
| Cost-aware synthesis | AiZynthFinder + custom building-block pricing | Score weight |
| Disconnection-aware design (DAD) | AiZynthFinder MCTS + verified rewards or post-search reranking | Compare route objectives explicitly |
| Patent-aware routes | Custom template exclusion | Specialized |

## AiZynthFinder Setup

**Goal:** Configure AiZynthFinder with USPTO templates + a building-block stock and run MCTS retrosynthesis planning on a target SMILES.

**Approach:** Create a version-appropriate YAML configuration with the expansion policy and stock, instantiate `AiZynthFinder` from that file, set the target SMILES, then call `tree_search()` followed by `build_routes()`. Use the schema documented for the installed release rather than copying legacy `policy` / `finder` dictionaries.

```python
from aizynthfinder.aizynthfinder import AiZynthFinder

finder = AiZynthFinder(configfile='config.yml')
finder.expansion_policy.select('uspto')  # key defined in config.yml
finder.stock.select('zinc')              # key defined in config.yml
# finder.filter_policy.select('uspto')   # optional configured filter
finder.target_smiles = 'CC(=O)Nc1ccc(C(=O)Nc2cccc(C(F)(F)F)c2)cc1'
finder.tree_search()
finder.build_routes()
```

Output: `finder.routes`, a `RouteCollection` containing ranked `reaction_trees`, initial `scores`, serialized route dictionaries, and route metadata.

## Route Output Analysis

```python
for tree, score in zip(finder.routes.reaction_trees, finder.routes.scores):
    leaves = list(tree.leafs())
    n_steps = len(list(tree.reactions()))
    print(f'Steps: {n_steps}, Score: {score}, Solved: {tree.is_solved}')
    print(f'In-stock: {sum(tree.in_stock(node) for node in leaves)} / {len(leaves)}')
    print(f'Building blocks: {[node.smiles for node in leaves]}')
```

Critical metrics:
- **Number of reactions**: count `tree.reactions()`; do not assume graph depth and reaction count are interchangeable for branched routes
- **Score**: interpret according to the configured scorer; scale and direction are scorer-specific
- **In-stock**: how many leaf nodes are commercially available
- **Solved state**: `tree.is_solved` is true only when all leaf nodes satisfy the configured stock criterion
- **Stock origin**: record the selected stock name, source snapshot, and access date

## Search Rewards and Route Scoring

AiZynthFinder retains the `mcts` search algorithm. It can combine configured search rewards through `search.algorithm_config.search_rewards` and corresponding weights, and it can rank completed routes with loaded scorers. There is no built-in `mo_mcts` algorithm or `finder.mo_mcts` configuration block. Available scorer names depend on the installed version, configuration, and plugins, so inspect the loaded scorers and use only documented names. If the desired objective is not available during search, export routes and rerank them explicitly after search.

## Building Block Stocks

| Stock source | Use | Required provenance |
|--------------|-----|---------------------|
| ZINC-derived snapshot | Publicly reproducible stock baseline | Download/source URL, filters, and snapshot date |
| Vendor building blocks | Purchase-oriented route termination | Vendor catalog version, region, and availability date |
| Make-on-demand catalog | Broader route termination | Catalog release, synthesis/lead-time assumptions, and access date |
| Custom internal stock | Organization-specific availability | Inclusion rules, identifiers, prices, and refresh date |

AiZynthFinder accepts HDF5 stocks built from plain-text SMILES with its documented `smiles2stock` command:
```bash
smiles2stock --files zinc_building_blocks.smi --output zinc.hdf5
```

## Forward Validation with Molecular Transformer

AiZynthFinder predicts retrosynthesis (target -> precursors), while a forward model predicts products from reactants. For each proposed reaction step, serialize reactants and reagents in the format expected by the installed forward model, request its ranked product predictions, and compare standardized product structures with the planned product. The Molecular Transformer literature does not define a universal `molecular_transformer.predict_forward` Python function, so use the documented interface of the chosen implementation. Report the observed top-k round-trip match rate for the model, reaction representation, and dataset; no universal 30–50% pass rate is established by the AiZynthFinder 4.0 paper.

## Template-Free with Chemformer

Chemformer uses a BART-style transformer trained on USPTO reactions for SMILES-to-SMILES prediction. Its official repository is archived and does not expose the `Chemformer.load_pretrained(...).predict(...)` convenience API sometimes shown in informal examples. To reproduce the published model, use the inference entry point, Hydra configuration, tokenizer, and checkpoint bundled with one pinned archived commit, and record that environment. For new work, prefer a maintained template-free implementation with a documented inference interface and benchmark it on the intended reaction domain.

**Trade-off:** A template-free model can propose disconnections outside a fixed template library, but its outputs require syntax checks, atom/reaction consistency checks, route-level review, and prospective validation. Treat it as a comparison or complementary hypothesis generator rather than assuming that merging its routes with AiZynthFinder is always superior.

## Disconnection-Aware Design (DAD)

Modify generative design to also score retrosynthetic feasibility with AiZynthFinder batch mode.

**Goal:** Add retrosynthetic feasibility scoring to generative design pipelines for hundreds-to-thousands of candidate molecules.

**Approach:** Batch-process generated SMILES through `aizynthcli`, use the reported solved state and number of reactions for each route, and feed a documented feasibility definition back into the generative scoring function.

```bash
aizynthcli --smiles compounds.smi --output routes.json \
           --config config.yaml --policy uspto --stocks zinc
```

For each compound, returns top-K routes. Score-feasibility for generative design:
- "Solved" = at least one extracted route has all leaves in the configured stock
- "Short solved route" = a solved route whose reaction count is below a project-defined threshold
- "Unsolved" = no extracted route is solved under the search budget and stock; this does not prove that the target is unsynthesizable

## Cost-Aware Synthesis

Add building-block pricing as objective:

```python
from rdkit import Chem

def route_cost(route, price_db):
    total = 0
    for leaf in route.leafs():
        smi = Chem.CanonSmiles(leaf.smiles)
        if smi not in price_db:
            raise KeyError(f'No observed building-block price for {smi}')
        total += price_db[smi]
    return total
```

Combine observed building-block prices with a project-specific reaction-cost model that documents labor, scale, yield, purification, and vendor assumptions. Do not apply a universal per-step cost.

## Per-Tool Failure Modes

### AiZynthFinder -- template coverage gap

**Trigger:** Target molecule uses bond formation not in training reactions.

**Mechanism:** USPTO templates are biased toward common transformations; novel chemistry (organometallics, exotic heterocycles) missing.

**Symptom:** No solved route or route uses unsuitable simplifications.

**Fix:** Use appropriately licensed additional templates, compare a maintained or exactly version-pinned template-free model, and perform manual review.

### Chemformer -- non-canonical SMILES output

**Trigger:** Default Chemformer output.

**Mechanism:** Transformer can produce non-canonical SMILES variants.

**Symptom:** SMILES round-trip fails; validation tools confused.

**Fix:** Canonicalize Chemformer output via RDKit before comparing.

### Route uses non-stock building block

**Trigger:** Leaf node not in stock database.

**Mechanism:** AiZynthFinder tree may end on non-purchasable molecules.

**Symptom:** Route "complete" but route has non-stock leaves.

**Fix:** Select routes whose `ReactionTree.is_solved` value is true, or explicitly require `tree.in_stock(leaf)` for every leaf. Expand the stock only when the additional availability definition is justified and versioned.

### MCTS iteration limit too low

**Trigger:** Complex target requiring deep tree search.

**Mechanism:** MCTS may not find route in default 100 iterations.

**Symptom:** No routes returned despite plausible target.

**Fix:** Increase and record the iteration or time budget in a controlled sensitivity analysis, inspect policy coverage and stock termination, and stop when additional search no longer changes the route conclusions. No fixed budget is universally adequate.

### Forward validation fails

**Trigger:** Retro route uses chemistry that doesn't actually work in forward.

**Mechanism:** Template-based retro lacks reaction conditions / catalysts; forward prediction more conservative.

**Symptom:** Forward predicts different product than target.

**Fix:** Use as confidence signal, not rejection; many routes don't round-trip but are still valid synthesis-wise.

### Building block stock obsolete

**Trigger:** Old ZINC catalog used; building blocks no longer purchasable.

**Mechanism:** Commercial catalogs and regional availability change over time.

**Symptom:** Routes recommend unavailable building blocks.

**Fix:** Refresh and date the selected stock snapshot, and verify vendor availability before synthesis.

## Reconciliation: AiZynthFinder vs Chemformer

| Aspect | AiZynthFinder | Chemformer |
|--------|---------------|------------|
| Approach | Templates + MCTS | Transformer encoder-decoder |
| Speed | Fast for shallow trees | Single-pass per target |
| Interpretability | High (template + atom mapping) | Low (black box) |
| Novel disconnections | Limited by selected templates | Can emit hypotheses outside a fixed template library, with no guarantee of validity |
| Production maturity | Maintained open-source package | Official repository archived; reproduce only with a pinned environment |
| Cost | CPU | GPU recommended |

If comparing both, standardize and validate their outputs independently before combining route hypotheses.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `tree_search()` returns no routes | Search budget, policy coverage, or stock criterion | Inspect each factor; compare a maintained or exactly version-pinned alternative |
| All extracted routes contain many reactions | Complex target or unsuitable disconnections | Compare scorer values and alternatives; review manually |
| Route appears solved but stock status is unclear | Reading node attributes instead of the route API, or stale stock provenance | Check `tree.is_solved` and `tree.in_stock(leaf)`; record the stock snapshot |
| Building block price not found | Compound not in pricing DB | Use Enamine quote or vendor inquiry |
| Chemformer truncates SMILES | Token limit | Increase max_length |
| Forward prediction wrong | Out-of-distribution reaction | Use as confidence signal only |
| MCTS slow on simple target | Default config | Reduce time_limit; use smaller template set |

## References

- Saigiridharan L, Hassen AK, Lai J, Torren-Peraire P, Engkvist O, Genheden S. "AiZynthFinder 4.0: developments based on learnings from 3 years of industrial application." *J. Cheminform.* 16:57 (2024). DOI: 10.1186/s13321-024-00860-x.
- Irwin R et al. "Chemformer: a pre-trained transformer for computational chemistry." *Mach. Learn.: Sci. Technol.* 3:015022 (2022). DOI: 10.1088/2632-2153/ac3ffb.
- Schwaller P et al. "Molecular Transformer: A Model for Uncertainty-Calibrated Chemical Reaction Prediction." *ACS Cent. Sci.* 5:1572–1583 (2019). DOI: 10.1021/acscentsci.9b00576.
- Tu Z et al. "ASKCOS: Open-Source, Data-Driven Synthesis Planning." *Acc. Chem. Res.* 58:1764–1775 (2025). DOI: 10.1021/acs.accounts.5c00155.
- Granqvist E, Mercado R, Genheden S. "Retrosynformer: planning multi-step chemical synthesis routes via a decision transformer." *Digital Discovery* 5:348–362 (2026). DOI: 10.1039/D5DD00153F.
- AiZynthFinder 4.4 Python interface and policy/stock selection: https://molecularai.github.io/aizynthfinder/python_interface.html.
- AiZynthFinder configuration documentation: https://molecularai.github.io/aizynthfinder/configuration.html.
- Chemformer official archived repository: https://github.com/MolecularAI/Chemformer.

## Related Skills

- chemoinformatics/molecular-io - Parse target and route SMILES
- chemoinformatics/molecular-standardization - Standardize before retrosynthesis
- chemoinformatics/generative-design - Add synthetic feasibility to scoring
- chemoinformatics/reaction-enumeration - Forward direction (template enumeration)
- chemoinformatics/admet-prediction - Filter targets before retrosynthesis
