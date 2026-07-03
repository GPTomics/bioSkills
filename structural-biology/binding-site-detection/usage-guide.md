# Binding Site Detection

## Overview

This skill finds putative ligand-binding pockets and druggable cavities de novo on a protein structure that has no bound ligand (an apo structure or a predicted model). It covers geometric pocket enumeration (fpocket alpha-spheres, CASTp surface topology), machine-learned ligandability ranking (P2Rank, DoGSiteScorer), and cryptic-pocket detection over a conformational ensemble (mdpocket). The governing idea is that a geometric cavity is a hypothesis, not automatically a functional or druggable site: pure geometry finds every concavity, including crystallization-additive sites and non-functional clefts, and druggability scores trained on holo/druggable reference sets systematically under-detect apo, shallow, allosteric, and cryptic pockets. When a ligand is already bound, mapping the residues around it (interface-analysis) is more reliable than any de-novo prediction.

## Prerequisites

- fpocket (external CLI, includes mdpocket) - install via conda (`conda install -c conda-forge fpocket`) or build from the Discngine/fpocket repository
- P2Rank (external CLI, `prank`) - download the standalone release from the rdk/p2rank repository (requires Java)
- CASTp and DoGSiteScorer are WEB SERVERS (no install): CASTp at http://sts.bioe.uic.edu/castp/ and DoGSiteScorer on ProteinsPlus at https://proteins.plus/
- Python for parsing outputs: `pip install biopython numpy`

## Quick Start

Tell the agent what you want to do:
- "Find the binding pockets on my apo structure and rank them by druggability"
- "Where could a ligand bind on this protein?"
- "Score the druggability of the pockets in this PDB file"
- "Detect cryptic pockets from my MD trajectory"
- "Is the top pocket on my AlphaFold model reliable enough to dock into?"
- "Run both fpocket and P2Rank and compare the top-ranked pocket"

## Example Prompts

### De-novo pocket detection
> "I have an apo crystal structure with no ligand bound. Detect all the surface pockets and tell me which one is most likely to bind a drug-like molecule."

> "Run fpocket on this PDB and give me the pockets ranked by druggability score, and flag any pocket that overlaps a crystallization additive."

### Ligandability scoring and method choice
> "Score the ligandability of the cavities in this structure with a machine-learning method that does not just use geometry, and compare its top pocket to fpocket's."

> "This is a protein-protein interaction target and the druggability scores all come out low. Does that mean there is no site, or is the scoring under-detecting a flat/apo pocket?"

### Cryptic and predicted-model cases
> "The known allosteric pocket is not open in my apo structure. How do I detect a cryptic pocket from a molecular dynamics ensemble?"

> "I only have an AlphaFold model. Can I trust the pocket it shows well enough to run docking, and what should I check first?"

## What the Agent Will Do

1. Establish the structure state (holo, apo, predicted model, or ensemble) and route accordingly: if a ligand is already bound, it maps the site with interface-analysis instead of detecting de novo.
2. For a predicted model, it flags the apo-like ground state and unreliable pocket-lining rotamers and recommends validation/preparation before detection.
3. Runs fpocket to enumerate pockets by alpha-sphere geometry and parses the `<stem>_out/<stem>_info.txt` descriptors, ranking by druggability score.
4. Runs P2Rank for a template-free ML ligandability ranking and parses the `_predictions.csv` probabilities, comparing the top pocket against fpocket.
5. Points to CASTp for pocket surface area/volume and DoGSiteScorer for a second SVM druggability opinion when a web-server assessment is wanted.
6. For a cryptic or transient site, runs mdpocket over an MD trajectory or ensemble and weights pockets by opening frequency.
7. Reports each pocket as a ranked hypothesis - stating the method and structure state - and hands the chosen site to virtual screening.

## Tips

- A geometric cavity is not a binding site. Cross-check every top pocket against the HETATM records in the input; a high-scoring cavity where a sulfate, glycerol, or PEG sat is a crystallization site, not a target.
- Druggability and ligandability scores were trained largely on holo, druggable sets. A low score on an apo, shallow, allosteric, or protein-protein-interface target does not mean there is no site - corroborate with a second tool and known biology.
- A cryptic pocket is closed in a single apo snapshot by definition. It only appears across a conformational ensemble, so use mdpocket over an MD or NMR ensemble, never one static file.
- On a predicted model, the pocket-lining side-chain rotamers are the least reliable atoms even when backbone confidence is high; validate and prepare the model first, and prefer a holo experimental structure of the target or a close homolog when one exists.
- Agreement between fpocket (geometry) and P2Rank (ML) on the same top pocket raises confidence; disagreement is a flag to inspect the site, not a tie to break by preference.
- Detect on the biological assembly, not the asymmetric unit, when a functional pocket may sit across an inter-chain interface.

## Related Skills

- structural-biology/interface-analysis - map residues around an already-bound ligand (holo); the complement of de-novo detection
- structural-biology/structure-validation - judge whether a structure or predicted model is trustworthy before detecting on it
- structural-biology/structure-preparation - clean, protonate, and fix the structure before detection (especially a predicted model)
- structural-biology/structure-modification - strip additives and resolve altlocs before detection
- structural-biology/alphafold-predictions - the apo-pocket and unreliable-rotamer caveats for predicted models
- structural-biology/structure-navigation - select chains and HETATM ligands to cross-check pockets against additives
- structural-biology/structure-io - download the biological assembly for inter-chain pockets
- chemoinformatics/virtual-screening - dock a library into a chosen pocket (presupposes a known site)
- chemoinformatics/ml-docking-rescoring - rescore docked poses in the detected site
