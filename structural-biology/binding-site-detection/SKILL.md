---
name: bio-structural-biology-binding-site-detection
description: Detects putative ligand-binding pockets and druggable cavities de novo on an apo protein structure with fpocket, P2Rank, CASTp, and DoGSiteScorer, ranking them by druggability/ligandability score. Use when detecting cavities on an apo structure with no bound ligand; choosing geometric pocket enumeration (fpocket alpha-spheres, CASTp) vs ML ligandability scoring (P2Rank, DoGSiteScorer); recognizing that a geometric cavity is a hypothesis, not automatically a functional or druggable site (may be a crystal-additive or non-functional cleft); knowing druggability scores were trained on holo sets and under-detect apo, shallow, and cryptic pockets; detecting cryptic or transient pockets over an MD or conformational ensemble (mdpocket); and detecting on a predicted model whose pocket-lining rotamers are the least reliable atoms. Keywords binding site, pocket, cavity, druggability, ligandability, fpocket, P2Rank, CASTp, DoGSiteScorer, apo, cryptic pocket, alpha sphere, mdpocket.
tool_type: mixed
primary_tool: fpocket
---

## Version Compatibility

Reference examples tested with: fpocket 4.1+, P2Rank 2.4+, biopython 1.83+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Binding Site Detection

**"Find the pockets / druggable cavities on my apo structure"** -> Enumerate surface concavities de novo and rank them by a ligandability or druggability score.
- CLI: `fpocket -f protein.pdb` (alpha-sphere / Voronoi) or `prank predict -f protein.pdb` (ML ligandability)

**"Where could a ligand bind, and how druggable is each site?"** -> Detect candidate pockets, then score each for the likelihood of binding a drug-like molecule.
- CLI: fpocket druggability score (per pocket) or P2Rank probability; CASTp / DoGSiteScorer web servers for surface area and SVM druggability

**"Detect a cryptic pocket that is not open in this single snapshot"** -> Track pocket opening across a conformational ensemble, not one static file.
- CLI: `mdpocket --trajectory_file traj.xtc --trajectory_format xtc -f topology.pdb` over an MD trajectory or NMR ensemble

## Governing Principle

A geometric cavity is NOT automatically a functional or druggable site. Pure geometry - fpocket alpha-spheres (Le Guilloux 2009 *BMC Bioinformatics* 10:168), CASTp surface topology (Tian 2018 *Nucleic Acids Res* 46:W363) - finds ALL concavities on the surface, including crystallization-additive sites (a sulfate, glycerol, or PEG cleft), non-functional surface pockets, and inter-domain crevices. A ranked pocket list is therefore a HYPOTHESIS, not a binding site. The single most reliable site is one that is already occupied: if the structure is holo (a ligand is bound), map the residues around that ligand directly - that is interface-analysis, not de-novo detection.

Druggability and ligandability SCORES (the fpocket druggability score, P2Rank probability, DoGSiteScorer SVM) rank pockets, but the ranking is not function and not a guarantee. These models were trained largely on holo, druggable reference sets (Schmidtke 2010 *J Med Chem* 53:5858; Volkamer 2012 *J Chem Inf Model* 52:360), so they systematically UNDER-detect and under-score apo, shallow, allosteric, and protein-protein-interface pockets, and they cannot see a CRYPTIC pocket at all - a cryptic site forms a pocket in the holo structure but is closed in the apo structure (Cimermancic 2016 *J Mol Biol* 428:709), so it only appears across an MD or conformational ensemble, never in one static apo snapshot.

Detecting on a PREDICTED model (AlphaFold, ESMFold) inherits that model's pocket-conformation unreliability. The model is usually an apo-like ground state that carries no ligand, cofactor, or metal that would shape the real site, and the side-chain rotamers lining a pocket are the LEAST reliable atoms even when the backbone pLDDT is high (see alphafold-predictions) - so a "confident" backbone can present a confidently wrong pocket. Validate and prepare the model first (structure-validation, structure-preparation). Whatever the input, state the method used and remember: ranking is not function, and the top pocket is a starting hypothesis to corroborate. The orthogonal corroboration signals are an independent detector, evolutionary conservation of the lining residues (a functional site is usually conserved, an additive/crystal cleft usually is not), hot-spot or mixed-solvent mapping (FTMap, mixed-solvent MD), a homolog holo structure, and known biology - not the geometric rank alone.

## Decision: method by question

| Question / situation | Method | Best when | Fails / misleads when |
|---|---|---|---|
| Enumerate every geometric concavity, fast | fpocket (alpha-sphere/Voronoi) or CASTp (surface topology) | Apo or holo, whole-surface scan, ranked candidate list | Reports additive sites, crystal clefts, and non-functional pockets as pockets; ranking != function |
| Rank pockets by ligandability / druggability | P2Rank (ML probability), fpocket druggability score, DoGSiteScorer (SVM) | Prioritizing among many cavities on a folded globular domain | Trained on holo/druggable sets -> under-scores apo, shallow, allosteric, and PPI sites |
| Detect a cryptic or transient pocket | mdpocket over an MD trajectory or conformational ensemble | Site opens only on binding / not present in a single apo snapshot | A single static structure - a cryptic site is invisible without conformational sampling |
| A ligand is already bound (holo structure) | interface-analysis (residues around the bound ligand) | The real binding site is known - map it, do not re-detect | Running de-novo detection when the answer is already in the file |
| Surface area / volume of a defined pocket | CASTp (analytic surface + volume) or fpocket descriptors | Quantifying a known pocket's size | Treating a large computed cavity as evidence of function |

The one-line rule: geometry enumerates concavities, scores rank them, an ensemble reveals cryptic ones, and a bound ligand settles it. Reach for the cheapest method that answers the actual question.

## Decision: structure state and what a pocket result means

| Structure state | What a detected pocket reveals | What it cannot reveal |
|---|---|---|
| Holo (ligand bound) | The real, occupied site - highest reliability; map it directly | Whether other detected pockets are functional (still hypotheses) |
| Apo experimental | Open, ligand-accessible pockets present in THIS conformation | Cryptic sites closed in this snapshot; induced-fit geometry |
| Predicted model (AlphaFold/ESMFold) | Candidate fold and gross concavities | Reliable pocket rotamers (least reliable atoms), cofactor/metal-shaped sites; validate first |
| MD / NMR ensemble | Transient and cryptic pockets, pocket dynamics/persistence | Which opening is biologically relevant (do not over-read rare openings) |

Apo detection answers "what could open here in this state", holo answers "what is actually bound", a predicted model answers "what might a pocket look like in the ground state", and an ensemble answers "what transiently opens". Match the claim to the state.

## Detect and rank pockets with fpocket

**Goal:** Enumerate candidate pockets on an apo structure and rank them by druggability, keeping the result explicitly a hypothesis list.

**Approach:** Run fpocket, which writes a `<stem>_out/` directory beside the input containing a `<stem>_info.txt` descriptor file and per-pocket `pockets/pocketN_atm.pdb` files. Parse the info file into per-pocket descriptors and sort by the druggability score. The 0.5 druggability cutoff is a rule-of-thumb triage line, not a law (Schmidtke 2010): above it a drug-like molecule is plausible, below it the pocket is likely undruggable by conventional small molecules, but apo and cryptic sites routinely fall below it while still being real.

```python
import subprocess
from pathlib import Path

def run_fpocket(pdb_path):
    subprocess.run(['fpocket', '-f', pdb_path], check=True)
    stem = Path(pdb_path).stem
    return Path(pdb_path).with_name(f'{stem}_out')  # fpocket writes <stem>_out/ beside the input

def parse_fpocket_info(out_dir):
    info = next(Path(out_dir).glob('*_info.txt'))
    pockets, current = [], None
    for raw in info.read_text().splitlines():
        line = raw.strip()
        if line.startswith('Pocket'):
            current = {'pocket': int(line.split()[1])}
            pockets.append(current)
        elif ':' in line and current is not None:
            key, val = (part.strip() for part in line.split(':', 1))
            current[key] = float(val)
    return pockets

out_dir = run_fpocket('protein.pdb')
pockets = sorted(parse_fpocket_info(out_dir),
                 key=lambda p: p['Druggability Score'], reverse=True)
druggable = [p for p in pockets if p['Druggability Score'] >= 0.5]  # 0.5 = triage, not a law
for p in pockets[:5]:
    print(f"pocket {p['pocket']}: drug={p['Druggability Score']:.2f} score={p['Score']:.2f}")
```

The pocket `Score` ranks how likely the cavity binds a small molecule; the `Druggability Score` (0-1) estimates likelihood of binding a DRUG-LIKE molecule specifically. A high-scoring pocket that coincides with a bound additive in the deposited file is a crystallization site, not a target - cross-check the pocket atoms against any HETATM in the input (structure-navigation).

## Score ligandability with P2Rank

**Goal:** Get a template-free, machine-learned ligandability ranking that does not depend on alpha-sphere geometry alone.

**Approach:** Run P2Rank, which scores points on the solvent-accessible surface and clusters them into ranked pockets, writing `<inputname>_predictions.csv`. Parse the calibrated `probability` per pocket. P2Rank is ML-based and template-free (Krivak 2018 *J Cheminform* 10:39), so it complements fpocket's pure geometry - agreement between the two on the same top pocket raises confidence; disagreement flags a site to scrutinize.

```python
import csv
import subprocess
from pathlib import Path

def run_p2rank(pdb_path, out_dir='p2rank_out'):
    subprocess.run(['prank', 'predict', '-f', pdb_path, '-o', out_dir], check=True)
    return Path(out_dir) / f'{Path(pdb_path).name}_predictions.csv'  # keeps the input extension

def parse_p2rank(csv_path):
    with open(csv_path) as fh:
        reader = csv.DictReader(fh, skipinitialspace=True)  # P2Rank pads columns with spaces
        rows = [{k.strip(): v.strip() for k, v in row.items()} for row in reader]
    return [{'rank': int(r['rank']), 'score': float(r['score']),
             'probability': float(r['probability'])} for r in rows]

predictions = parse_p2rank(run_p2rank('protein.pdb'))
top = predictions[0]  # rows are already ordered by rank
print(f"top pocket: probability={top['probability']:.2f} score={top['score']:.2f}")
```

P2Rank ranks pockets even when every probability is modest; a low top probability on an apo or PPI target is expected, not a failure to find a site.

## Surface topology and SVM druggability servers

CASTp 3.0 (Tian 2018 *Nucleic Acids Res* 46:W363, http://sts.bioe.uic.edu/castp/) is a web server that delineates pockets, interior cavities, and channels analytically and reports each one's molecular surface AREA and VOLUME - use it when the question is pocket geometry (size, mouth, buried volume) rather than a druggability call. DoGSiteScorer, on the ProteinsPlus server (Volkamer 2012 *J Chem Inf Model* 52:360, https://proteins.plus/), predicts pockets by a difference-of-Gaussians grid and returns an SVM druggability score trained on a druggable/undruggable set - a second, independent druggability opinion to compare against fpocket. Both are servers with no Python binding; submit a structure and parse the returned table. Treat a server druggability number the same way as fpocket's: trained on holo/druggable data, so it under-scores apo and cryptic sites.

## Cryptic pockets over an ensemble

**Goal:** Find a site that is closed in the single apo structure but opens on binding.

**Approach:** A cryptic pocket does not exist in one static snapshot, so run pocket detection over a conformational ENSEMBLE - an MD trajectory, an NMR ensemble, or frames from enhanced sampling that is designed to open cryptic sites (mixed-solvent MD / SWISH, accelerated MD) - and track pocket persistence and opening frequency. Generating that ensemble is the real upstream step; mdpocket (Schmidtke 2011 *Bioinformatics* 27:3276) then analyzes it, building a pocket-frequency grid across frames, and the `-S` flag adds a per-snapshot drug score so a druggable pocket that only appears transiently is visible.

```python
import subprocess

def run_mdpocket(topology_pdb, trajectory, fmt='xtc'):
    # -S adds the per-snapshot drug score so a transiently druggable pocket is flagged
    subprocess.run(['mdpocket', '--trajectory_file', trajectory,
                    '--trajectory_format', fmt, '-f', topology_pdb, '-S'], check=True)
```

Do not over-read a pocket that opens in a handful of frames - a rare, short-lived opening is a weaker hypothesis than a persistent one, and enhanced sampling can manufacture openings that are not physiologically populated. Once a cryptic site is chosen, the bound-state geometry (not the closed apo structure) is what a docking run needs.

## Once a site is chosen

A detected, ranked pocket is the INPUT to structure-based drug discovery, not the endpoint. Hand the chosen pocket (its center and lining residues) to chemoinformatics/virtual-screening to dock a library into it, and to chemoinformatics/ml-docking-rescoring to rescore poses - both of those presuppose a known site, which is exactly what this skill produces. Carry the caveat forward: docking into an apo or predicted pocket with unreliable rotamers gives confidently wrong poses unless the pocket is prepared and, ideally, cross-checked against a holo conformation.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Top-ranked pocket is a surface additive site | Geometry finds all concavities, including where a sulfate/glycerol/PEG sat | Cross-check pocket atoms against HETATM/additives in the input; strip additives first (structure-modification) |
| A known binding site is missed entirely | Site is cryptic - closed in the single apo snapshot | Run mdpocket over an MD/NMR ensemble, not one static structure |
| Druggability scores all low on a real target | Score trained on holo/druggable sets under-scores apo, shallow, allosteric, PPI sites | Low score is not "no site"; corroborate with an independent method (P2Rank/DoGSiteScorer), evolutionary conservation of the lining residues (ConSurf), hot-spot / mixed-solvent mapping (FTMap), homolog holo structures, and known biology |
| fpocket finds no `_info.txt` | Parsed the wrong directory or fpocket failed silently | fpocket writes `<stem>_out/` beside the input; check the run completed and the stem matches |
| P2Rank parse gives KeyError on column | The CSV has a space after each comma so header keys carry a leading space (` probability`) | Use `csv.DictReader(skipinitialspace=True)` and strip keys before indexing |
| Pockets differ between apo and holo of the same protein | Induced fit - the pocket reshapes on binding | Expected; detect on the state that matches the question, do not average them |
| "Confident" AlphaFold pocket docks poorly | Pocket-lining rotamers are the least reliable atoms; model is apo-like | Validate/prepare the model (structure-validation); prefer a holo experimental structure |
| Same pocket, different rank from two tools | fpocket geometry vs P2Rank ML weight features differently | Agreement raises confidence; treat disagreement as a flag to inspect, not a tie-break |
| Reported "N druggable pockets" as a finding | Ranking treated as function | State it as a ranked hypothesis list; a score is not a validated site |
| mdpocket pocket appears in few frames trusted as real | Over-reading a rare, short-lived opening | Weight by opening frequency/persistence; rare openings are weak hypotheses |
| Detection on the asymmetric unit misses an interface pocket | Functional pocket sits across a symmetry-generated interface | Detect on the biological assembly (structure-io) when the site may be inter-chain |

## Related Skills

- structural-biology/interface-analysis - map residues around an ALREADY-BOUND ligand (holo); the complement of de-novo detection
- structural-biology/structure-validation - judge whether a structure or predicted model is trustworthy before detecting on it
- structural-biology/structure-preparation - clean, protonate, and fix the structure before pocket detection (especially a predicted model)
- structural-biology/structure-modification - strip crystallization additives and resolve altlocs before pocket detection
- structural-biology/alphafold-predictions - the apo-pocket and unreliable-rotamer caveats for detecting on a predicted model
- structural-biology/structure-navigation - select chains and HETATM ligands to cross-check pockets against additives
- structural-biology/structure-io - download the biological assembly (not just the asymmetric unit) for inter-chain pockets
- chemoinformatics/virtual-screening - dock a library into a chosen pocket (presupposes a known site)
- chemoinformatics/ml-docking-rescoring - rescore docked poses in the detected site

## References

- Le Guilloux V, Schmidtke P, Tuffery P (2009) Fpocket: an open source platform for ligand pocket detection. *BMC Bioinformatics* 10:168. (fpocket; Voronoi tessellation and alpha spheres)
- Schmidtke P, Barril X (2010) Understanding and predicting druggability. A high-throughput method for detection of drug binding sites. *J Med Chem* 53(15):5858-5867. (fpocket druggability score; trained on druggable/undruggable cavities)
- Krivak R, Hoksza D (2018) P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. *J Cheminform* 10:39. (P2Rank; template-free ML ligandability)
- Tian W, Chen C, Lei X, Zhao J, Liang J (2018) CASTp 3.0: computed atlas of surface topography of proteins. *Nucleic Acids Res* 46(W1):W363-W367. (surface pockets, cavities, channels; area and volume)
- Volkamer A, Kuhn D, Grombacher T, Rippmann F, Rarey M (2012) Combining global and local measures for structure-based druggability predictions. *J Chem Inf Model* 52(2):360-372. (DoGSiteScorer; SVM druggability)
- Cimermancic P, et al. (2016) CryptoSite: expanding the druggable proteome by characterization and prediction of cryptic binding sites. *J Mol Biol* 428(4):709-719. (cryptic site = pocket in holo but not apo)
- Schmidtke P, Bidon-Chanal A, Luque FJ, Barril X (2011) MDpocket: open-source cavity detection and characterization on molecular dynamics trajectories. *Bioinformatics* 27(23):3276-3285. (pocket tracking over an ensemble; cryptic/transient pockets)
