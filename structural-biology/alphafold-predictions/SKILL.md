---
name: bio-structural-biology-alphafold-predictions
description: Retrieves and interprets AlphaFold Protein Structure Database (AFDB) models by UniProt accession, reading pLDDT and PAE confidence correctly. Use when treating pLDDT as PER-RESIDUE confidence (not global accuracy) and recognizing a long low-pLDDT stretch as an intrinsically disordered region rather than a modeling error; reading PAE to segment confident domains and judge inter-domain/relative-position confidence that high mean pLDDT cannot certify; recognizing a static AFDB model carries NO ligands, ions, cofactors, PTMs, quaternary assembly, or alternative conformations (pLDDT sits in the B-factor column with opposite polarity to thermal motion); and deciding an AFDB entry vs re-running prediction. Keywords AlphaFold DB, pLDDT, PAE, B-factor column, intrinsic disorder, UniProt, Foldseek.
tool_type: python
primary_tool: requests
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, numpy 1.26+, requests 2.31+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# AlphaFold Predictions

**"Get the AlphaFold model for my protein and tell me which parts to believe"** -> Fetch the precomputed AFDB entry by UniProt accession, then read its confidence files (pLDDT per-residue, PAE per-residue-pair) to decide which regions and which geometry are trustworthy.
- Python: `requests.get(f'https://alphafold.ebi.ac.uk/api/prediction/{accession}')` returns metadata with `cifUrl`/`pdbUrl`/`paeDocUrl` download links.

## Governing Principle

An AFDB entry is ONE AlphaFold2 prediction of a SINGLE UniProt sequence modeled as an isolated chain in vacuum - a per-chain hypothesis of a dominant fold, not an experimental structure and not a biological state. It carries NO ligands, ions, cofactors, metals, or substrates (the pocket is apo even when the protein only folds around a cofactor), NO post-translational modifications, NO quaternary structure or biological assembly (it is a monomer even for obligate oligomers), NO alternative conformations (one static snapshot - kinases render in one activation state, transporters in one gate state), and NO membrane context. The cardinal error is using an AFDB coordinate file as an experimental holo complex instead of as a scored guess whose own confidence files tell the reader which parts to believe.

Three confidence traps ride on top of this. First, pLDDT is written into the B-FACTOR COLUMN but is per-residue CONFIDENCE (0-100, higher=better) with OPPOSITE polarity to a real B-factor - any tool that reads that column as thermal motion inverts the meaning, and feeding raw pLDDT-as-B into crystallographic refinement mis-weights it. Second, a long low-pLDDT stretch is usually a genuine INTRINSICALLY DISORDERED REGION (pLDDT is competitive with dedicated IDR predictors; Akdel 2022 *Nat Struct Mol Biol* 29:1056; Piovesan 2022 *Protein Sci* 31:e4466), not a modeling failure - trimming it as "junk" discards real biology. Third, pLDDT is LOCAL and per-residue: high mean pLDDT does NOT certify inter-domain placement. Two rigid domains can each be 95 pLDDT yet float at an unknown relative orientation - reading that is PAE's job. "Confident" bounds the structural self-consistency of the prediction, NOT its biological correctness.

## Decision: pLDDT confidence bands

| pLDDT | Band | Operational meaning | Trust for |
|-------|------|---------------------|-----------|
| > 90 | Very high | Backbone AND well-oriented side chains | Side-chain detail, catalytic-geometry hypotheses, MR core |
| 70-90 | Confident | Backbone generally correct | Fold, domain topology, backbone-level MR; be wary of side chains |
| 50-70 | Low | Backbone uncertain, cautionary zone | Coarse topology at best; never trust details |
| < 50 | Very low | Ribbon is a placeholder; frequently an IDR | Disorder signal, NOT a conformation |

Band cutoffs 90/70/50 are the AFDB-defined thresholds (Jumper 2021 *Nature* 596:583). A very-low band is a disorder SIGNAL, not proof of error - a conditionally-folded binding region looks disordered in AFDB yet folds on binding a partner AFDB never sees.

## Decision: downstream-use suitability

| Use | AFDB monomer suitability | Key caveat |
|-----|--------------------------|------------|
| Remote-homology / fold search (Foldseek) | Excellent | Feed the confident core; a match is a hypothesis |
| Molecular replacement | Very good after processing | Trim + pLDDT->pseudo-B + PAE domain split first |
| Fold / domain-architecture analysis | Good | Segment by PAE, not the stitched cartoon |
| Disorder / IDR annotation | Good (pLDDT as predictor) | Low pLDDT = disorder signal, not error |
| Ligand docking / virtual screening | Poor to moderate | Apo pocket, unreliable rotamers, wrong gate/state, absent cofactor |
| Mechanism / catalytic geometry | Poor without holo | No ligands/metals/PTMs; wrong-state risk |
| Quaternary structure / interfaces | Not applicable | Monomer only - run AF3 / AF-Multimer |
| Conformational ensembles / allostery | Not applicable | Single static state |

Foldseek structure search over AFDB (van Kempen 2024 *Nat Biotechnol* 42:243) is the transformative use - it encodes backbone into the 3Di alphabet and finds structural homologs invisible to sequence search (see alignment/structural-alignment). Molecular replacement needs the model PROCESSED first: `phenix.process_predicted_model` (Oeffner 2022 *Acta Cryst D* 78:1303) reads pLDDT from the B column, converts pLDDT->pseudo-B, trims below ~0.7 fractional pLDDT, and splits into PAE-defined domains. Docking into an AFDB pocket gives confidently wrong poses when the backbone is confident but the rotamers, gate state, or cofactor are not.

## Decision: use the AFDB entry vs run a new prediction

| Situation | Use AFDB | Run a new prediction |
|-----------|----------|----------------------|
| Single well-covered UniProt monomer, fold-level question | Yes | No |
| Need a complex, assembly, or interface | No (monomer only) | Yes - AF3 / AF-Multimer |
| Need ligand / ion / cofactor / PTM context | No | AF3 or dock into experimental |
| Designed or mutant sequence not in UniProt | No | Yes |
| Want deeper / custom MSA depth | No (MSA is fixed) | Yes |
| Want a specific alternative conformation | No (single state) | Yes (subsampled MSA) - still hard |
| Very long non-human protein (> 2700 aa) | Often absent | Yes (domain-wise) or ESM Atlas |

Anything needing complexes, ligands, mutants, custom MSA depth, or a specific state points to modern-structure-prediction. AFDB is monomer-only and FIXED at deposition - a newer method or deeper MSA is not reflected.

## Fetch the AFDB entry via REST metadata

**Goal:** Download the coordinate file and PAE for a UniProt accession without hard-coding a version suffix.

**Approach:** Query the prediction metadata endpoint, which returns the current download URLs (`cifUrl`, `pdbUrl`, `paeDocUrl`); the version token drifts (v4 -> v6 as of 2025) so the URLs are discovered, never assembled by hand.

```python
import requests

def afdb_metadata(accession):
    url = f'https://alphafold.ebi.ac.uk/api/prediction/{accession}'
    r = requests.get(url)
    r.raise_for_status()
    return r.json()  # list; one object per fragment/isoform, empty if no model exists

def fetch_afdb(accession, out_dir='.'):
    entries = afdb_metadata(accession)
    if not entries:
        return None  # >2700-aa non-human proteins and non-UniProt sequences are often absent
    entry = entries[0]
    cif_text = requests.get(entry['cifUrl']).text
    pae_json = requests.get(entry['paeDocUrl']).json()
    cif_path = f"{out_dir}/AF-{accession}.cif"
    with open(cif_path, 'w') as f:
        f.write(cif_text)
    return cif_path, pae_json

result = fetch_afdb('P04637')  # human p53
```

Long proteins split into fragments `AF-{accession}-F1-...`, `-F2-...` (human > 2700 aa, ~1400-aa windows shifted by 200); `afdb_metadata` returns one entry per fragment. Relative placement ACROSS fragments is independent and must not be trusted.

## Read pLDDT from the B-factor column

**Goal:** Extract per-residue confidence and classify each residue into a band.

**Approach:** Parse the model, read the B-factor field of the CA atom (that is where AFDB stores pLDDT), and map the score through the 90/70/50 cutoffs.

```python
from Bio.PDB import MMCIFParser

def extract_plddt(cif_file):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('af', cif_file)
    plddt = {}
    for residue in structure[0].get_residues():
        if residue.id[0] == ' ' and 'CA' in residue:
            # pLDDT rides in the B-factor column but is CONFIDENCE (0-100, higher=better),
            # opposite polarity to a thermal B-factor - never read it as motion
            plddt[residue.id[1]] = residue['CA'].get_bfactor()
    return plddt

def plddt_band(score):
    if score > 90:                       # 90: side-chain-trustworthy core (AFDB very-high cut)
        return 'very_high'
    if score >= 70:                      # 70: backbone-reliable fold (AFDB confident cut)
        return 'confident'
    if score >= 50:                      # 50: coarse topology only below this
        return 'low'
    return 'very_low'                    # <50: usually an intrinsically disordered region

plddt = extract_plddt('AF-P04637.cif')
mean_plddt = sum(plddt.values()) / len(plddt)
disordered = [res for res, s in plddt.items() if s < 50]  # candidate IDR, not error
```

A high mean does not license inter-domain claims - a globally 90-pLDDT model can still have two domains at an unconstrained orientation. Check PAE before measuring any inter-domain distance.

## Read PAE to segment domains and judge relative orientation

**Goal:** Decide which residue pairs have a confident relative position and where to split the model into independent rigid bodies.

**Approach:** Load the compact PAE matrix; low off-diagonal blocks mark domains whose relative orientation is confident, bright (high) inter-block regions mark independently-placed domains to segment.

```python
import numpy as np

def load_pae(pae_json):
    entry = pae_json[0] if isinstance(pae_json, list) else pae_json
    # Compact format (2023+): 2D num_res x num_res array (values rounded to integer).
    # Legacy 1D 'distance'/'residue1'/'residue2' fields were removed - do not read them.
    return np.array(entry['predicted_aligned_error'])

def interdomain_confidence(pae, domain_a, domain_b):
    # PAE is asymmetric (aligning on i vs j differs); average both off-diagonal blocks.
    block = np.concatenate([pae[np.ix_(domain_a, domain_b)].ravel(),
                            pae[np.ix_(domain_b, domain_a)].ravel()])
    return block.mean()  # low (roughly < 5 A) = confident relative placement; high = a guess

pae = load_pae(fetch_afdb('P04637')[1])
```

Confident low-PAE squares along the diagonal define the confidently-predicted DOMAINS; a bright inter-block region means "these two domains are correctly folded individually but their relative arrangement is unconstrained - treat them as separate rigid bodies." This is exactly how AFDB defines predicted domains and how MR pipelines split a search model.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Model "colored by flexibility" looks inverted | Read the B-factor column as thermal motion | It is pLDDT (confidence, higher=better); color by pLDDT bands |
| Low-pLDDT tail deleted, then a known IDR/linker is missing | Treated low pLDDT as "wrong" | Low pLDDT usually = disorder; keep and annotate as IDR, cross-check a sequence disorder predictor |
| Confident domains, but inter-domain distance is nonsense | Trusted the stitched cartoon on high mean pLDDT | pLDDT is local; read PAE - high off-diagonal PAE = unconstrained relative placement |
| Docking scores look great but validate poorly | Docked into an apo AFDB pocket | Pocket lacks the ligand/cofactor and has unreliable rotamers; use a holo structure or flexible docking |
| Catalytic-geometry conclusion contradicts experiment | AFDB has no metals/ligands/PTMs and one static state | Do not read mechanism from a monomer apo model; get a holo structure |
| 404 / empty metadata list for a large protein | > 2700-aa non-human proteins are excluded; non-UniProt sequences absent | Run a new prediction (domain-wise) or query the ESM Metagenomic Atlas |
| Only residues 1-1400 returned for a long human protein | The model is fragmented (F1, F2 ...) | Iterate all metadata entries; never trust placement across fragments |
| `KeyError: 'distance'` loading PAE | Code written for the legacy 1D PAE JSON | Read the 2D `predicted_aligned_error` array from the compact format |
| Hard-coded `..._v4.cif` URL 404s | The version suffix advanced (v6 as of 2025) | Discover URLs from `/api/prediction/{accession}` metadata, do not assemble them |
| pLDDT looks fine but the biological state is wrong | Modeled the wrong assembly/conformation confidently | Confidence bounds self-consistency, not biological correctness; ask what context AFDB could not see |
| Foldseek hits are noisy or low-quality | Fed the low-pLDDT spaghetti into the search | 3Di is backbone geometry; search the confident core only |

## Related Skills

- structural-biology/modern-structure-prediction - run a new prediction when a complex, ligand, mutant, or specific state is needed
- structural-biology/structure-io - parse and convert the downloaded PDB/mmCIF
- structural-biology/geometric-analysis - RMSD, superposition, and per-residue deviation against an experimental structure
- structural-biology/structure-modification - trim low-pLDDT regions or write pLDDT into the B-factor column for coloring
- structural-biology/structure-preparation - add hydrogens and protonation states before docking or MD on the model
- structural-biology/binding-site-detection - find pockets on the predicted model (apo, unreliable rotamers)
- alignment/structural-alignment - Foldseek 3Di search over AFDB for remote-homology detection
- database-access/uniprot-access - resolve gene names and sequences to the UniProt accession AFDB is keyed on
- database-access/remote-homology - sequence-level homology search to complement structure search

## References

- Jumper J, et al. (2021) Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583-589.
- Tunyasuvunakool K, et al. (2021) Highly accurate protein structure prediction for the human proteome. *Nature* 596:590-596.
- Varadi M, et al. (2022) AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Res* 50:D439-D444.
- Varadi M, et al. (2024) AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Res* 52:D368-D375.
- Akdel M, et al. (2022) A structural biology community assessment of AlphaFold2 applications. *Nat Struct Mol Biol* 29:1056-1067.
- Piovesan D, Monzon AM, Tosatto SCE (2022) Intrinsic protein disorder and conditional folding in AlphaFoldDB. *Protein Sci* 31:e4466.
- van Kempen M, et al. (2024) Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* 42:243-246.
- Oeffner RD, et al. (2022) Putting AlphaFold models to work with phenix.process_predicted_model and ISOLDE. *Acta Cryst D* 78:1303-1314.
- Terwilliger TC, et al. (2024) AlphaFold predictions are valuable hypotheses and accelerate but do not replace experimental structure determination. *Nat Methods* 21:110-116.
