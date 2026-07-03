---
name: bio-structural-biology-modern-structure-prediction
description: Predicts protein and complex structures with deep-learning models (ESMFold, AlphaFold2/ColabFold, AlphaFold3, Chai-1, Boltz-1/2) and reconciles them with confidence metrics. Use when choosing a predictor by input and question rather than novelty (ESMFold single-chain, no-MSA, fast, metagenomic-scale vs AlphaFold3/Chai-1/Boltz for complexes, ligands, nucleic acids, ions, PTMs); recognizing that MSA depth is the dominant accuracy determinant so ESMFold trades accuracy for speed and degrades on orphan proteins; gating a complex on ipTM plus inter-chain PAE, not per-chain pLDDT; reading pLDDT as local confidence, PAE as inter-domain/inter-chain positioning, pTM as global fold; knowing a single prediction is one dominant conformer not an ensemble (no apo/holo, allosteric, or fold-switch states), that these are not variant-effect/ddG/affinity engines, and that a confident prediction is a hypothesis, not an experiment. Keywords ESMFold, AlphaFold3, Chai-1, Boltz-1, ColabFold, ipTM, PAE, pLDDT, MSA depth.
tool_type: python
primary_tool: ESMFold
---

## Version Compatibility

Reference examples tested with: fair-esm 2.0+, biopython 1.83+, numpy 1.26+, requests 2.31+, chai_lab 0.6+, boltz 2.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Modern Structure Prediction

**"Predict the structure of my protein"** -> Map an amino-acid (and optionally ligand/nucleic-acid) sequence to a single 3D model plus per-residue and pairwise confidence.
- Python: ESMFold local via `esm.pretrained.esmfold_v1()` (no MSA); ColabFold/AlphaFold via MMseqs2 MSA; AlphaFold3/Chai-1/Boltz for complexes and ligands.

## Governing Principle: a prediction is a hypothesis, and the predictor is chosen by the input and the question

The trap is treating a prediction as an answer and picking a model by novelty. Two facts govern every decision here. First, MSA depth is the dominant accuracy determinant for the coevolution-based models (AlphaFold2/3, Chai-1, Boltz): quality tracks how well-represented the sequence's family is, not how hard the biology is, so these models excel on deep-MSA families and degrade on orphan, fast-evolving, viral, or de-novo-designed sequences (Jumper 2021 *Nature* 596:583; Lin 2023 *Science* 379:1123). ESMFold is single-sequence with no MSA, so it is fast enough for metagenomic scale but is lower-accuracy on average and degrades hardest exactly where evolutionary signal is thin. Second, a default prediction is ONE dominant conformer, not an ensemble: it does not give apo vs holo, allosteric states, or fold switches, and it carries no Boltzmann populations. MSA subsampling and AF-Cluster sample some alternate states but are unreliable, seed-sensitive hypotheses (Wayment-Steele 2024 *Nature* 625:832), a generality directly challenged by a Matters Arising (Schafer & Porter 2025 *Nature* 638:E8-E12).

Three category errors follow and must be avoided. (1) These are not variant-effect, ddG, or stability engines: a single point mutation barely changes a deep MSA, so wild-type vs mutant predictions come back near-identical with near-identical pLDDT, and the model is insensitive to the mutation by construction (Buel & Walters 2022 *Nat Struct Mol Biol* 29:1; Pak 2023 *PLoS ONE* 18:e0282689). Use AlphaMissense, FoldX/Rosetta ddG, or ESM/EVE variant scores instead. (2) No co-folder gives a trustworthy Kd from geometry: a plausible complex or ligand pose is not evidence of binding or affinity, and CASP16 assessors found co-fold affinity ranking essentially unreliable; Boltz-2's affinity module is a screening prior only, not a measured constant. (3) AF3-class diffusion models can hallucinate confident-looking order in genuinely disordered regions and have measurable chirality violations (~4.4% on PoseBusters) and atom clashes (Abramson 2024 *Nature* 630:493), so every ligand pose needs a physical-validity check. A confident prediction is a starting hypothesis with spatially varying reliability; validate it against experiment before believing any part of it.

## Decision: which predictor for the job

| Job | Preferred | Why / caveat |
|---|---|---|
| Single-domain monomer, MAX accuracy | AlphaFold2 (LocalColabFold) or AlphaFold3 | Deep MSA = best accuracy; AF2 mature and well-understood |
| Monomer, deep MSA, fast and free | ColabFold (MMseqs2 MSA) | 40-60x faster search, near-AF2 accuracy (Mirdita 2022) |
| Metagenomic / genome-scale / triage | ESMFold | Fastest (no MSA); lower accuracy, weak on large/low-family proteins |
| Single-sequence when no homologs exist | ESMFold or Chai-1 (single-seq mode) | Both skip MSA; expect reduced accuracy, sanity-check hard |
| Protein-protein COMPLEX | AF-Multimer / AF3 / Boltz / Chai-1 | Gate on ipTM + inter-chain PAE, NOT per-chain pLDDT |
| Complex WITH ligand/ion/nucleic acid/PTM | AF3, Boltz-1/2, or Chai-1 | Co-folders; validate the POSE (PoseBusters), NOT affinity |
| Binding-affinity PRIOR for screening | Boltz-2 (affinity module) | Screening prior only; a relative affinity score, not a trusted Kd |
| Commercial / on-prem deployment | Boltz-1/2 or Chai-1 (both Apache-2.0/MIT, commercial OK) | AF3 weights are non-commercial (Google terms) |
| Alternative conformational states | AF2 + MSA subsampling / AF-Cluster | UNRELIABLE; hypotheses only, not ensembles or populations |
| Variant effect / stability / pathogenicity | NOT these tools | Insensitive to point mutations; use AlphaMissense/FoldX/ESM |

Licenses drift; verify before deploying. AF2 code Apache-2.0, weights CC-BY-4.0 (permissive). AF3 code Apache-2.0, weights under the non-commercial "AlphaFold 3 Model Parameters Terms of Use", granted on request to non-commercial orgs and received directly from Google (open, not open-source). Boltz-1 and Boltz-2 are MIT (code + weights, commercial use permitted). Chai-1 was relicensed to Apache-2.0 for both code and weights in November 2024 (commercial use, including drug discovery, permitted; it launched Sept 2024 under a restrictive non-commercial license, so older notes may say otherwise - verify terms). ESMFold code/weights are MIT.

## Decision: which confidence metric answers which question

| Metric | Scope | Answers | Read it for |
|---|---|---|---|
| pLDDT (0-100) | Per-residue, LOCAL | How well-placed is this residue's local environment | Trimming; a long <50 stretch usually flags an intrinsically disordered region, not an error |
| PAE (Angstrom) | Residue-pair | Expected error at j when aligned on i | Domain packing, linker geometry, inter-chain arrangement |
| pTM (0-1) | Whole model, GLOBAL | Estimated TM-score of the overall fold | Is the topology plausible (>~0.5) |
| ipTM (0-1) | Interface | Accuracy of relative subunit positioning | Complex interface reliability (>~0.8 likely, <~0.6 unreliable) |

The load-bearing reads: high per-residue pLDDT with a high inter-domain PAE block means each domain is confident internally but their relative arrangement is unknown - do not trust the linker or domain-domain interface. For a complex, judge the interface on ipTM plus the inter-chain PAE block; a complex can have high pLDDT on both chains and still be a garbage interface. AF-Multimer ranks models by 0.8*ipTM + 0.2*pTM, deliberately weighting the interface. All these metrics are self-reported and can be confidently wrong together on out-of-distribution inputs.

## Predict a monomer with ESMFold (fast, no MSA)

**Goal:** Get a single-chain model in seconds without building an MSA, and read pLDDT off the B-factor column.

**Approach:** Run ESMFold locally with `esm.pretrained.esmfold_v1()`; the hosted esmatlas API is intermittently down (SSL/internal-server errors) so local is the reliable path. pLDDT rides in the B-factor column but is confidence, not a temperature factor.

```python
import torch
import esm

model = esm.pretrained.esmfold_v1().eval().to('cuda')  # needs ~16 GB GPU for typical proteins

sequence = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH'
with torch.no_grad():
    pdb_text = model.infer_pdb(sequence)  # returns a PDB string with pLDDT in B-factor column
with open('esmfold.pdb', 'w') as f:
    f.write(pdb_text)
```

Hosted-API fallback (only when no GPU and the endpoint is up):

```python
import requests

url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
resp = requests.post(url, data=sequence, timeout=300)  # 300 s: long sequences take minutes
resp.raise_for_status()
pdb_text = resp.text
```

## Read per-residue confidence

**Goal:** Summarize where a prediction is trustworthy so downstream use is restricted to confident cores.

**Approach:** pLDDT sits in the B-factor column of every prediction (ESMFold, AlphaFold, co-folders). Band it into the standard cutoffs; a contiguous very-low band usually marks a disordered region, not a failure.

```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('pred', 'esmfold.pdb')
plddt = {res.id[1]: res['CA'].get_bfactor() for res in structure[0].get_residues() if 'CA' in res}

# Bands from the AlphaFold/EBI convention: >90 very high, 70-90 confident, 50-70 low, <50 very low.
very_high = [r for r, s in plddt.items() if s > 90]
confident = [r for r, s in plddt.items() if 70 <= s <= 90]
very_low = [r for r, s in plddt.items() if s < 50]  # likely intrinsically disordered, not wrong
print(f'mean pLDDT {sum(plddt.values())/len(plddt):.1f}; {len(very_low)} very-low residues')
```

## Predict a complex and gate on the interface

**Goal:** Model a protein-protein or protein-ligand complex and decide whether to believe the interface.

**Approach:** Use a co-folder (Chai-1 or Boltz), then accept the interface only if ipTM and the inter-chain PAE block agree. Chai-1 and Boltz run from the CLI; both default to no MSA and can call an MSA server. Verify the exact CLI with `--help` since these packages evolve fast.

```python
import subprocess

# Chai-1: one FASTA with a header per chain; '--use-msa-server' fetches an MSA (improves accuracy).
# Reference invocation - confirm with `chai-lab fold --help`.
subprocess.run(['chai-lab', 'fold', '--use-msa-server', 'complex.fasta', 'chai_out/'], check=True)

# Boltz: FASTA or YAML input; YAML is required to request the Boltz-2 affinity module.
# Reference invocation - confirm with `boltz predict --help`.
subprocess.run(['boltz', 'predict', 'complex.fasta', '--use_msa_server'], check=True)
```

**Goal:** Gate the predicted interface before trusting any cross-chain distance.

**Approach:** Read ipTM and pTM from the confidence JSON the co-folder writes, and band the interface on ipTM. Below ~0.6 the interface is unreliable or the chains likely do not interact; 0.6-0.8 is uncertain and the inter-chain PAE block decides; a confident interface wants ipTM > ~0.8.

```python
import json

with open('chai_out/scores.model_idx_0.json') as f:  # exact filename varies by tool/version
    conf = json.load(f)
iptm = conf.get('iptm')
ptm = conf.get('ptm')
if iptm is None or iptm < 0.6:        # <0.6: unreliable or chains likely do not interact
    print(f'interface NOT reliable (ipTM={iptm}); inspect inter-chain PAE before any claim')
elif iptm < 0.8:                      # 0.6-0.8: uncertain - the inter-chain PAE block decides
    print(f'interface UNCERTAIN (ipTM={iptm:.2f}); gate on the inter-chain PAE block')
else:
    ptm_str = f'{ptm:.2f}' if ptm is not None else 'NA'  # some score files omit pTM
    print(f'interface confident (ipTM={iptm:.2f}, pTM={ptm_str}); still confirm with inter-chain PAE')
```

## Prepare an AlphaFold3 server job

**Goal:** Submit a monomer or complex to the AlphaFold Server without local weights.

**Approach:** The server takes a JSON job listing entities and seeds; multiple seeds sample the diffusion head, so request several and inspect the spread rather than trusting one sample.

```python
import json

def af3_job(sequences, name='prediction', seeds=(1, 2, 3)):
    entities = [{'proteinChain': {'sequence': s, 'count': 1}} for s in sequences]
    return json.dumps([{'name': name, 'modelSeeds': list(seeds), 'sequences': entities}], indent=2)

job_json = af3_job(['MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH'])
```

## Reconcile multiple predictions

**Goal:** Compare models from different predictors and locate the regions they agree on.

**Approach:** Superimpose on a fixed CA correspondence and report pairwise RMSD, but treat RMSD as fold-agreement only where the aligned selection is stated; prefer a length-normalized fold metric (TM-score) for cross-method fold claims. See geometric-analysis for TM-score and superposition caveats.

```python
from Bio.PDB import PDBParser, Superimposer

def ca_rmsd(pdb_a, pdb_b):
    parser = PDBParser(QUIET=True)
    a = [r['CA'] for r in parser.get_structure('a', pdb_a)[0].get_residues() if 'CA' in r]
    b = [r['CA'] for r in parser.get_structure('b', pdb_b)[0].get_residues() if 'CA' in r]
    n = min(len(a), len(b))  # Superimposer needs an equal-length ordered atom correspondence
    sup = Superimposer()
    sup.set_atoms(a[:n], b[:n])
    return sup.rms

print(f'ESMFold vs AF3 CA-RMSD: {ca_rmsd("esmfold.pdb", "af3.pdb"):.2f} Angstrom')
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Mutant and wild-type predictions look identical | A point mutation barely changes a deep MSA; the model is insensitive to it | Do not read structure/pLDDT deltas as variant effect; use AlphaMissense, FoldX, or ESM |
| Confident model but wrong in the lab | Prediction is one dominant conformer, not an ensemble; no apo/holo/allosteric states | Treat as a hypothesis; sample states cautiously (MSA subsampling) and validate experimentally |
| Complex accepted on high per-chain pLDDT | pLDDT is intra-chain local confidence, blind to the interface | Gate on ipTM + inter-chain PAE block; reject interface if ipTM < ~0.6 (0.6-0.8 uncertain) |
| Long low-pLDDT stretch treated as an error | Low pLDDT correlates with intrinsic disorder | Read <50 regions as likely IDRs (biologically real flexibility), not modeling failures |
| Two domains confident but arrangement wrong | High intra-domain pLDDT with high inter-domain PAE | Trust each domain, not the relative orientation or linker; split at high-PAE hinges |
| ESMFold much worse than AlphaFold on an orphan | ESMFold is single-sequence and degrades where evolutionary signal is thin | Use MSA-based ColabFold/AF for orphan/de-novo proteins; keep ESMFold for scale |
| Ligand pose has wrong chirality or clashes | AF3-class diffusion can violate stereochemistry (~4.4% chirality) | Run PoseBusters/validity checks on every pose; do not assume physical plausibility |
| Reported Kd from a co-fold pose | Co-folders give geometry, not affinity; CASP16 found affinity ranking unreliable | Use Boltz-2 affinity only as a screening prior; confirm with FEP or experiment |
| esmatlas API returns SSL / internal-server error | The hosted ESMFold endpoint is intermittently down | Run ESMFold locally via `esm.pretrained.esmfold_v1()` |
| Two predictions "disagree" but were run differently | Different MSA depth/source, recycles, seeds, or templates change the answer | Report the MSA pipeline and settings; two predictions are not comparable if these differ |
| RMSD between predictions looks huge for the same fold | Global all-atom RMSD is dominated by flexible loops and needs a stated selection | Superimpose on CA/core and report the selection; use TM-score for fold agreement |

## Related Skills

- alphafold-predictions - Retrieve precomputed AlphaFold DB models and read their pLDDT/PAE
- structure-io - Parse and write predicted PDB/mmCIF files
- geometric-analysis - RMSD, superposition, and TM-score caveats for comparing models
- structure-navigation - Walk chains/residues/atoms in a predicted structure
- structure-preparation - Trim, add hydrogens, and protonate a predicted model before docking or MD
- binding-site-detection - Detect pockets on a predicted model (inherits apo/rotamer uncertainty)
- alignment/structural-alignment - Structure-based alignment before comparing sequence-different models
- chemoinformatics/virtual-screening - Dock into a predicted pocket (inherits predicted rotamer/backbone error)
- chemoinformatics/ml-docking-rescoring - Rescore co-folded poses; co-fold geometry is not affinity

## References

Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. Nature 596:583-589 (2021). doi:10.1038/s41586-021-03819-2.
Abramson J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. Nature 630:493-500 (2024). doi:10.1038/s41586-024-07487-w.
Lin Z, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379:1123-1130 (2023). doi:10.1126/science.ade2574.
Evans R, et al. Protein complex prediction with AlphaFold-Multimer. bioRxiv 2021.10.04.463034 (2021, preprint). doi:10.1101/2021.10.04.463034.
Baek M, et al. Accurate prediction of protein structures and interactions using a three-track neural network. Science 373:871-876 (2021). doi:10.1126/science.abj8754.
Mirdita M, et al. ColabFold: making protein folding accessible to all. Nat Methods 19:679-682 (2022). doi:10.1038/s41592-022-01488-1.
Wayment-Steele HK, et al. Predicting multiple conformations via sequence clustering and AlphaFold2. Nature 625:832-839 (2024). doi:10.1038/s41586-023-06832-9.
Schafer JW, ..., Porter LL (2025) Sequence clustering confounds AlphaFold2 (Matters Arising). *Nature* 638:E8-E12. doi:10.1038/s41586-024-08267-2.
Buel GR, Walters KJ. Can AlphaFold2 predict the impact of missense mutations on structure? Nat Struct Mol Biol 29:1-2 (2022). doi:10.1038/s41594-021-00714-2.
Pak MA, et al. Using AlphaFold to predict the impact of single mutations on protein stability and function. PLoS ONE 18:e0282689 (2023). doi:10.1371/journal.pone.0282689.
Wohlwend J, et al. Boltz-1: democratizing biomolecular interaction modeling. bioRxiv 2024.11.19.624167 (2024, preprint). doi:10.1101/2024.11.19.624167.
Chai Discovery. Chai-1: decoding the molecular interactions of life. bioRxiv 2024.10.10.615955 (2024, preprint). doi:10.1101/2024.10.10.615955.
