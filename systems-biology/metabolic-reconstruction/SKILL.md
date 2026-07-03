---
name: bio-systems-biology-metabolic-reconstruction
description: Builds draft genome-scale metabolic models from an annotated genome using CarveMe (top-down carving of a BiGG universal model) or gapseq (bottom-up pathway-evidence reconstruction), then loads and sanity-checks the draft in COBRApy. Use when creating a model for an organism without one, choosing between CarveMe and gapseq, gap-filling to a target medium, understanding why a draft that grows is still only a hypothesis, handling BiGG-vs-ModelSEED namespace mismatch, or preparing a draft for curation and community modeling.
tool_type: cli
primary_tool: CarveMe
---

## Version Compatibility

Reference examples tested with: CarveMe 1.6+, gapseq 1.2+, COBRApy 0.29+, DIAMOND 2.1+, Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: CarveMe needs an LP solver (academic CPLEX/Gurobi; SCIP is a slow open-source fallback) and a DIAMOND install, and its universal model is BiGG-derived so the BiGG-model release matters. gapseq is cloned from GitHub (not pip-installable), emits ModelSEED-namespace models, and its reference DB version matters. Model predictions are only comparable within the same tool, DB, and namespace.

# Metabolic Reconstruction

**"Build a metabolic model from my organism's genome"** -> Map the annotated proteome/genome to reactions in a reference database, assemble a draft network with a biomass reaction, and gap-fill so it can grow on a chosen medium.
- CLI: `carve genome.faa -o model.xml` (CarveMe, top-down); `gapseq doall genome.fna` (gapseq, bottom-up)

## The governing principle: a draft is a hypothesis, and "it grows" is guaranteed by construction

Automated reconstruction produces a DRAFT, not a finished model. The single most misleading signal is growth: CarveMe and gapseq GAP-FILL specifically to force biomass production on a chosen medium, so a draft that grows proves nothing biological - it was made to grow. Consequences:

- Gap-filled reactions are the least-evidenced part of the model (added to close a hole, not because homology supports them), and the gap-fill medium determines what gets added. Gap-filling on the wrong medium bakes in the wrong reactions. Flag gap-filled reactions as low-confidence and record the medium.
- Draft quality is bounded by the annotation and the reference database. CarveMe can only ever include reactions in the BiGG universe (biased toward well-studied organisms); gapseq's homology thresholds and pathway logic set its floor. Peripheral/novel metabolism is systematically underrepresented.
- The draft is the START of curation, not the end. Different tools produce markedly different models from the same genome, and successive models of the same organism (iJR904 -> iJO1366 -> iML1515) give different predictions. Reconstruction feeds model-curation, never bypasses it.

## Decision: CarveMe vs gapseq vs ModelSEED

| Goal | Tool | Why / trade-off |
|------|------|-----------------|
| Fast draft(s) for well-studied bacteria; batch/community | CarveMe (`carve`) | top-down MILP carving of a curated BiGG universe; minutes; universe is simulation-ready but BiGG-centric; universal biomass; weak transporters |
| Non-model/environmental clade; carbon-source & fermentation phenotypes | gapseq | bottom-up homology + pathway-completeness; slower, more transparent; better SCFA/carbon-use recovery; ModelSEED namespace complicates merging |
| Fully-automated web pipeline (RAST annotation) | ModelSEED/KBase | template-based; convenient; template biomass and aggressive gap-fill can force implausible reactions |
| Eukaryotes / fungi / actinomycetes | RAVEN (MATLAB) | KEGG/MetaCyc-based, template or de novo; MATLAB license; the eukaryote-capable option |

Do NOT treat "CarveMe and gapseq do the same thing, pick the faster one" as true: different philosophies, namespaces (BiGG vs ModelSEED), and failure modes. The choice is scientific. No single tool dominates - which is why consensus/ensemble reconstruction exists.

## CarveMe (top-down)

```bash
pip install carveme          # also needs DIAMOND and an LP solver (CPLEX/Gurobi; SCIP fallback)

# Draft from a PROTEIN FASTA (default input). Raw/GenBank genomes are NOT accepted.
carve genome.faa -o model.xml

# Gram type and universe are VALUES of -u/--universe, NOT --grampos/--gramneg flags.
carve genome.faa -o model.xml -u grampos    # {bacteria (default), grampos, gramneg, archaea, cyanobacteria}

# Gap-fill to force growth on a medium (opt-in; records what was added for that medium).
carve genome.faa -o model.xml --gapfill M9
carve genome.faa -o model.xml -u gramneg --gapfill M9,LB   # multiple media

# Nucleotide input instead of protein, or download by accession:
carve genome.fna --dna -o model.xml
```

Community reconstruction uses a SEPARATE `merge_community` command (not `carve`); see systems-biology/community-metabolic-modeling.

## gapseq (bottom-up, pathway-evidence)

```bash
git clone https://github.com/jotech/gapseq && cd gapseq && ./gapseq test   # cloned, not pip; check deps

# One-shot: find + find-transport + draft + fill
./gapseq doall genome.fna

# Or the explicit steps (note find-transport is its OWN subcommand, not `find -t`):
./gapseq find -p all genome.fna          # -> genome-all-Reactions.tbl, genome-all-Pathways.tbl
./gapseq find-transport genome.fna       # -> genome-Transporter.tbl  (singular)
./gapseq draft -r genome-all-Reactions.tbl -t genome-Transporter.tbl \
               -p genome-all-Pathways.tbl -c genome.fna   # -> genome-draft.RDS, genome-rxnWeights.RDS
./gapseq fill -m genome-draft.RDS -n dat/media/M9.csv \
              -c genome-rxnWeights.RDS -g genome-rxnXgenes.RDS   # -> genome.xml / genome.RDS
```

## Load and Sanity-Check the Draft

**Goal:** Read the draft, confirm it grows on the gap-fill medium, and inventory the parts most likely to be wrong.

**Approach:** Load the SBML into COBRApy, report network size and gene coverage, test growth, and count orphan (gene-less) reactions and exchanges - the draft's soft spots before curation.

```python
import cobra

model = cobra.io.read_sbml_model('model.xml')
print(f'reactions={len(model.reactions)} metabolites={len(model.metabolites)} genes={len(model.genes)}')
print(f'grows on gap-fill medium: {model.slim_optimize() > 1e-3}')   # true by construction if gap-filled
orphans = [r for r in model.reactions if not r.genes]   # no GPR: gap-filled, spontaneous, or transport
print(f'orphan (gene-less) reactions: {len(orphans)}  exchanges: {len(model.exchanges)}')
# Typical bacterial draft: ~1000-2500 reactions. Far outside that range flags an annotation problem.
```

## Namespaces (the silent killer of model comparison)

```python
# Reaction/metabolite IDs come from the tool's reference DB: CarveMe = BiGG, gapseq/ModelSEED =
# ModelSEED (seed.*), RAVEN = KEGG/MetaCyc. Two models in different namespaces cannot be merged or
# compared directly. Reconcile through MetaNetX/MNXref (MNXM* metabolites, MNXR* reactions) BEFORE
# any cross-tool merge or community build. This BiGG-vs-ModelSEED split is exactly why community
# modeling of CarveMe + gapseq outputs breaks without reconciliation.
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `carve` errors on a genome file | GenBank/nucleotide passed where protein FASTA expected | supply a protein FASTA, or add `--dna` for nucleotide |
| `--grampos`/`--gramneg` not recognized | those are `-u/--universe` VALUES, not flags | `carve ... -u grampos` |
| Draft cannot grow at all | no gap-filling requested, or wrong medium | add `--gapfill <medium>`; confirm the medium supplies biomass precursors |
| Draft grows on everything / implausibly | gap-fill forced reactions for the chosen medium | flag gap-filled reactions low-confidence; re-gap-fill on the correct medium; curate |
| Two models will not merge / IDs mismatch | different namespaces (BiGG vs ModelSEED) | reconcile via MetaNetX/MNXref before merging |
| gapseq `find -t` fails | transport is the `find-transport` subcommand | use `./gapseq find-transport genome.fna` |
| Very few genes / tiny network | poor annotation or wrong input file | check the proteome/annotation; verify gene IDs |

## Related Skills

- systems-biology/model-curation - Curate, gap-fill deliberately, and validate the draft (the required next step)
- systems-biology/flux-balance-analysis - Predict growth/flux once the model is trustworthy
- systems-biology/community-metabolic-modeling - Combine reconstructions into a community model
- genome-annotation/prokaryotic-annotation - Produce the annotated protein FASTA CarveMe/gapseq consume
- database-access/ncbi-datasets-cli - Fetch genome/proteome inputs

## References

- Machado D, Andrejev S, Tramontano M, Patil KR. 2018. Fast automated reconstruction of genome-scale metabolic models for microbial species and communities. *Nucleic Acids Res* 46(15):7542-7553. (CarveMe)
- Zimmermann J, Kaleta C, Waschina S. 2021. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models. *Genome Biol* 22(1):81.
- Henry CS, DeJongh M, Best AA, et al. 2010. High-throughput generation, optimization and analysis of genome-scale metabolic models. *Nat Biotechnol* 28(9):977-982. (ModelSEED)
- Wang H, Marcisauskas S, Sanchez BJ, et al. 2018. RAVEN 2.0: a versatile toolbox for metabolic network reconstruction. *PLoS Comput Biol* 14(10):e1006541.
- Thiele I, Palsson BO. 2010. A protocol for generating a high-quality genome-scale metabolic reconstruction. *Nat Protoc* 5(1):93-121.
- Mendoza SN, Olivier BG, Molenaar D, Teusink B. 2019. A systematic assessment of current genome-scale metabolic reconstruction tools. *Genome Biol* 20(1):158.
- Moretti S, Tran VDT, Mehl F, et al. 2021. MetaNetX/MNXref: unified namespace for metabolites and biochemical reactions. *Nucleic Acids Res* 49(D1):D570-D574.
- Feist AM, Palsson BO. 2010. The biomass objective function. *Curr Opin Microbiol* 13(3):344-349.
- Monk JM, Lloyd CJ, Brunk E, et al. 2017. iML1515, a knowledgebase that computes Escherichia coli traits. *Nat Biotechnol* 35(10):904-908.
