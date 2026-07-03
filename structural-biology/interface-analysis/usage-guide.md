# Interface Analysis

## Overview

This skill maps protein-protein and protein-ligand interfaces from a coordinate file and quantifies them. It covers contact-residue detection under an explicit distance cutoff, structural epitope and binding-site (ligand-contact) residues, and buried surface area (BSA) as the physical measure of interface size. It also frames the two judgment calls that dominate interface work: which contact definition to use (and always stating the cutoff), and whether an interface observed in the crystal is biological or a packing artifact. The primary tool is Biopython Bio.PDB (NeighborSearch and ShrakeRupley), with freesasa and PDBePISA as complements.

## Prerequisites

```bash
pip install biopython numpy freesasa
```

- A structure file (PDB or mmCIF), ideally the biological assembly, not just the asymmetric unit.
- Familiarity with SASA basics (see structural-biology/geometric-analysis).
- PDBePISA (https://www.ebi.ac.uk/pdbe/pisa/) is a web service, no install needed, for the biological-assembly call.

## Quick Start

Tell the agent what interface question you have:
- "List the residues where chain A contacts chain B in this complex"
- "Which residues line the ATP binding site?"
- "Map the structural epitope this antibody contacts on the antigen"
- "Compute the buried surface area of this dimer interface"
- "Is this crystal interface biological or just packing?"
- "Find the salt bridges across this interface"

## Example Prompts

### Contact and interface residues
> "Report every residue of chain B within 4.5A of chain A, and tell me which cutoff you used and why."
> "Give me the contact map between the two chains at an 8A CA-CA definition, and explain what that captures versus a heavy-atom cutoff."

### Ligand and epitope mapping
> "Identify the binding-site residues within 4A of the bound inhibitor (HETATM LIG)."
> "This is an antibody-antigen complex; map the structural epitope on the antigen chain and the paratope on the antibody."

### Buried surface area
> "Compute the buried surface area of the interface between chains A and B and report the per-partner area, keeping the SASA parameters identical across the three calculations."
> "Compare the interface size of these two dimers computed the same way."

### Biological vs crystal interface
> "This deposition has two chains in the asymmetric unit; is their interface biological or a crystal-packing contact? Weigh BSA, H-bonds, and what PISA says, and tell me how confident I can be."
> "Download the biological assembly for this PDB and compute the interface on that, not the deposited asymmetric unit."

## What the Agent Will Do

1. Confirm the oligomeric state to analyze, preferring the biological assembly over the asymmetric unit, and download it if needed.
2. Resolve provenance issues first: pick one altloc, decide on waters and crystallization additives, and note missing/disordered residues.
3. Choose and state a contact definition (heavy-atom 4-5A for physical contact, CA-CA 8A for topology) with its rationale.
4. Detect contact residues with NeighborSearch, keeping only cross-partner pairs, and report the interface residue lists per chain.
5. For a ligand or epitope, select the target group (HETATM or partner chain) and collect residues within the cutoff.
6. Quantify the interface as buried surface area, computing complex and isolated-part SASA with identical ShrakeRupley settings.
7. When biological relevance is in question, assemble the probabilistic signals (BSA magnitude, H-bonds/salt bridges, conservation, PISA CSS) and report a hypothesis with a confidence caveat, not a verdict.

## Tips

- A contact is defined by an arbitrary cutoff; the residue count changes with it, so always state the cutoff and the atom set (heavy-atom vs CA-CA vs including hydrogens).
- A contact list is not an interface. The physical measure is buried surface area; compute all three SASA terms with identical parameters or the subtraction is meaningless.
- Compute interfaces on the biological assembly. An interface in the asymmetric unit can be pure crystal packing.
- "Biological interface" is a hypothesis. PISA is ~80-90% accurate with known false positives; larger BSA, more H-bonds/salt bridges, complementarity, and conservation raise confidence but do not prove it.
- SASA depends on the 1.4A probe radius and the algorithm; an absolute BSA is only comparable to another BSA computed the same way. Use freesasa for Lee-Richards or full radii control.
- Salt-bridge and H-bond geometric definitions are loose and tool-dependent (distances quoted anywhere from 3.2 to 5.0A); state the exact criteria, and remember angles cannot be checked without modeled hydrogens.
- Double-counting alternate conformations invents impossible contacts; select one altloc before searching.
- Missing loops are disorder, not real chain gaps; reconcile against SEQRES before drawing interface conclusions near a gap.

## Related Skills

- geometric-analysis - SASA fundamentals, NeighborSearch, distances that this skill builds on
- structure-io - download the biological assembly before interface analysis
- structure-modification - resolve altlocs and strip waters before contact detection
- structure-navigation - select chains, residues, and HETATM ligands by identity
- structure-validation - check the region is well-fit before trusting an interface
- immunoinformatics/epitope-prediction - structural epitope mapping
- chemoinformatics/virtual-screening - binding-site definition for docking
- alignment/structural-alignment - superpose complexes before comparing interfaces
