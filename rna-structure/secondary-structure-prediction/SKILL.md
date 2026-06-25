---
name: bio-rna-structure-secondary-structure-prediction
description: Predicts RNA secondary structure with ViennaRNA, treating the Boltzmann ensemble (partition function, base-pair probabilities, centroid, MEA, stochastic samples) as the object rather than a single MFE fold. Covers consensus folding from alignments (RNAalifold), SHAPE-constrained folding, RNA-RNA interaction (RNAcofold/RNAduplex/RNAup), local and linear-time methods for long RNA, and pseudoknot-aware tools. Use when folding an RNA and choosing between MFE, centroid, MEA, or ensemble sampling; judging whether a single structure is well-defined; folding long RNAs where a global MFE is meaningless; handling suspected pseudoknots; or weighing thermodynamic versus comparative versus deep-learning prediction.
tool_type: cli
primary_tool: ViennaRNA
---

## Version Compatibility

Reference examples tested with: ViennaRNA 2.6+, matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Secondary Structure Prediction

**"Predict the secondary structure of my RNA sequence"** -> Compute base pairs under a nearest-neighbor thermodynamic model, but report the Boltzmann ENSEMBLE (partition function, base-pair probabilities, centroid/MEA, per-base confidence), not a single fold.
- CLI: `RNAfold -p` for single-sequence ensemble folding
- CLI: `RNAalifold` for consensus structure from an alignment
- CLI: `RNAcofold` / `RNAduplex` / `RNAup` for RNA-RNA interaction
- Python: `import RNA` (`fold_compound`) for scripted ensemble analysis

## The governing principle: the MFE is one sample from an ensemble, not "the structure"

The minimum free energy (MFE) structure is the single lowest-energy fold, but every possible structure has probability proportional to exp(-G/RT). The partition function (McCaskill 1990) sums over that whole Boltzmann ensemble and yields the probability of each base pair, not just one fold. The MFE is frequently NOT the biologically relevant structure: riboswitches sit at a poised switch between two folds, many RNAs are genuinely dynamic, and the functional fold is often slightly suboptimal. Report the MFE alone and the uncertainty is hidden.

Three consequences shape every decision below:
- Trust the ENSEMBLE, not the point estimate. Use partition-function quantities (base-pair probability > 0.9, low positional entropy, low ensemble diversity) as the per-pair and per-base confidence. A low MFE with a diffuse base-pair-probability matrix means the single structure is not trustworthy.
- Accuracy is modest and length-dependent. Single-sequence thermodynamic folding recovers ~70-73% of base pairs for short, well-behaved RNAs (tRNA, 5S rRNA) and degrades sharply beyond ~700 nt. A single global MFE of a multi-kilobase mRNA or viral genome is close to meaningless -> use local (RNAplfold), linear-time (LinearFold/LinearPartition), comparative, or probing-restrained methods.
- ViennaRNA folds only nested canonical + wobble (G-U) pairs. Pseudoknots are silently excluded (general pseudoknot prediction is NP-hard); non-canonical pairs (sheared G-A, base triples, the Leontis-Westhof geometric families) are invisible. A confidently wrong nested fold is the failure mode for frameshift elements, riboswitch aptamers, and viral UTRs.

## Which ViennaRNA program for which question

| Question / input | Tool | Why |
|---|---|---|
| Fold one sequence, want structure + confidence | RNAfold -p | MFE + partition function (centroid/MEA/diversity) |
| Ensemble free energy only, no dot plot (speed) | RNAfold -p0 | skips base-pair probabilities, ~50% faster |
| Consensus structure from an alignment of homologs | RNAalifold | thermodynamics + covariation |
| Two RNAs that dimerize (full model + concentrations) | RNAcofold (-c) | intra+inter pairs, equilibrium species |
| Fast two-strand hybridization screen | RNAduplex | inter-molecular pairs only, no internal structure |
| sRNA/miRNA-target where site accessibility matters | RNAup | opening energy + hybridization (correct for buried sites) |
| Sample alternative conformations | RNAsubopt -p N / `fc.pbacktrack(n)` | Boltzmann sampling |
| Long mRNA: local pairing / target accessibility | RNAplfold | windowed pair + unpaired probabilities |
| Long RNA/genome: find local structured elements | RNALfold | locally stable structures, bounded span |
| Sequence longer than a few kb | LinearFold / LinearPartition | O(n), avoids a meaningless global O(n^3) fold |
| Pseudoknot suspected | IPknot / ProbKnot / Knotty | nested folders structurally cannot |
| Have SHAPE/DMS reactivities | RNAfold --shape (Deigan) | restrain folding with experimental data |

## Which "answer" to report

| Goal | Answer | Note |
|---|---|---|
| Quick single estimate | MFE | over-calls weak pairs; not for long/low-complexity RNA |
| Conservative, ensemble-representative structure | centroid | minimum total base-pair distance to ensemble; fewer false pairs, can under-pair |
| Best single structure (esp. with probing data) | MEA (tune gamma) | high gamma -> more pairs/recall, low gamma -> fewer/precision |
| Per-pair / per-base confidence | base-pair probability (>0.9) + positional entropy | structure-agnostic confidence track |
| Conformational switching / multiple states | stochastic sampling | cluster the samples into populations |
| Is one structure well-defined? | ensemble diversity (low) + ensemble defect | length-relative, not an absolute cutoff |

## RNAfold: single-sequence ensemble folding

```bash
# MFE only
echo "GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA" | RNAfold --noPS

# Ensemble: partition function + base-pair probabilities + centroid + MEA + ensemble diversity
# --noPS suppresses the *_dp.ps / *_ss.ps PostScript files RNAfold writes to the CWD by default
echo ">myRNA" > rna.fa && echo "GGGCUAUUAGCUCAGUUGGUUAGAGCGCACC" >> rna.fa
RNAfold -p --MEA --noLP --noPS < rna.fa
```

Key flags (verified against the current manpage):

| Option | Effect |
|--------|--------|
| `-p` | partition function + base-pair-probability matrix (unlocks centroid/MEA/diversity) |
| `-p0` | ensemble free energy ONLY, no base-pair probabilities (faster) |
| `--MEA[=gamma]` | maximum-expected-accuracy structure (default gamma 1.0); `--MEA` implies `-p` |
| `-d2` | dangling-end model (default); use `-d0` for comparative/alignment folding to avoid dangle artifacts |
| `-d3` | also allow coaxial stacking of adjacent helices in multiloops (MFE folding only; the partition function `-p` ignores `-d3` and falls back to `-d2`, so ensemble quantities do not reflect it) |
| `--noLP` | forbid lonely (isolated) base pairs; standard for well-folded RNA |
| `--maxBPspan N` | cap base-pair span; crude knob for long sequences |
| `-T 37` | folding temperature in Celsius (default 37) |
| `--shape FILE` / `--shapeMethod` | SHAPE-directed folding (see structure-probing) |
| `-g` | allow G-quadruplex formation (default off; turn on for G-rich sequences where G4s compete with canonical pairing) |
| `--noPS` | suppress PostScript drawings (always set in scripts to avoid CWD clutter) |

## Centroid, MEA, sampling, and per-base confidence (Python)

`fc.pf()` MUST be called before `bpp()`, `centroid()`, `MEA()`, `pbacktrack()`, `positional_entropy()`, `ensemble_defect()`, or `mean_bp_distance()` -- they all read the partition-function matrices and silently return empty/garbage otherwise.

```python
import RNA

seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
fc = RNA.fold_compound(seq)

mfe_struct, mfe = fc.mfe()
_, ensemble_g = fc.pf()                 # partition function; ensemble G <= MFE always

centroid, _ = fc.centroid()             # conservative, fewer false-positive pairs
mea_struct, _ = fc.MEA()                # MEA(gamma); gamma>1 favors pairing (recall), gamma<1 precision
diversity = fc.mean_bp_distance()       # ensemble diversity: low = well-defined (read RELATIVE to length)
defect = fc.ensemble_defect(mfe_struct) # expected wrongly-paired positions of this structure vs ensemble
entropy = fc.positional_entropy()       # per-base Shannon entropy: low = confidently paired-or-unpaired

```

Sample alternative conformations from the Boltzmann ensemble (riboswitches, bistable RNA) with the CLI `RNAsubopt -p N` (N stochastic samples) or the Python `fc.pbacktrack(N)` after `fc.pf()`; cluster the samples to find conformational populations. `fc.pbacktrack` requires a ViennaRNA build with stochastic backtracking enabled -- if it returns nothing, use `RNAsubopt -p N`.

Decision rule: one well-defined structure -> centroid or MEA; report stability/confidence -> partition-function quantities (base-pair probabilities + positional entropy); conformational switching -> stochastic sampling; quick single estimate -> MFE (with the caveats above).

## Constrained and SHAPE-directed folding

```python
import RNA

seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'

# Hard constraints: force positions unpaired or paired
fc = RNA.fold_compound(seq)
fc.hc_add_up(35, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)   # 1-indexed; force position 35 unpaired
fc.hc_add_bp(1, 72, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
constrained, c_mfe = fc.mfe()

# Soft SHAPE pseudo-energy restraint (Deigan model). The reactivity vector is 1-INDEXED:
# prepend -999 as a placeholder for index 0; -999 elsewhere means "no data" (NOT zero reactivity).
# m=1.8, b=-0.6 is the standard SHAPE pair (Hajdin 2013, the ViennaRNA default), not Deigan 2009's own m=2.6/b=-0.8.
fc2 = RNA.fold_compound(seq)
reactivities = [-999.0] + [0.1, 0.05, 0.8, 0.9] + [-999.0] * (len(seq) - 4)
fc2.sc_add_SHAPE_deigan(reactivities, 1.8, -0.6)
shape_struct, shape_mfe = fc2.mfe()   # this energy INCLUDES the SHAPE pseudo-energy; do NOT compare it to the unrestrained MFE
```

The energy returned after a SHAPE restraint folds in the pseudo-energy bonus, so it is not on the same scale as an unconstrained MFE -- compare the STRUCTURES (base-pair distance, SHAPE agreement), not the two energy numbers. See structure-probing for obtaining reactivities and for the SHAPE-vs-DMS parameter choice.

## Comparative (consensus) folding: homologs beat a single sequence

Evolution conserves structure while sequence drifts, so a compensatory substitution (an A-U in one species becoming G-C at the same two columns) is direct evidence of a real pair that thermodynamics alone cannot see. RNAalifold folds an alignment with a combined thermodynamic + covariation score.

```bash
# Consensus structure; format (Stockholm/Clustal/FASTA) is auto-detected, alignment is positional.
# --ribosum_scoring improves covariation detection; -d0 avoids dangle artifacts at gapped columns.
RNAalifold --ribosum_scoring -d0 -p --noPS alignment.sto
```

| Option | Effect |
|--------|--------|
| `--cfactor` | covariation weight (default 1.0; lower leans on thermodynamics) |
| `--nfactor` | penalty for sequences that cannot form the consensus pair (default 1.0) |
| `--ribosum_scoring` | use RIBOSUM covariation matrices (recommended) |
| `-p` | consensus partition function + base-pair probabilities |

RNAalifold accuracy depends entirely on alignment quality and real covariation: near-identical sequences carry no covariation signal and it degrades toward noisy single-sequence folding. RNAalifold assumes a FIXED, correct alignment; when homologs cannot be aligned reliably, TurboFold II co-estimates alignment AND structure across the sequences jointly (Tan et al. 2017) and is the better choice. A predicted consensus structure is a HYPOTHESIS until covariation is statistically validated -- test it with R-scape (see covariation-analysis), which found no significant covariation support for the proposed HOTAIR/Xist/SRA lncRNA structures.

## Long RNA: do not fold one global structure

For a multi-kilobase mRNA, lncRNA, or viral genome, a single global O(n^3) MFE both over-pairs across long ranges and is slow, and folding is co-transcriptional and local in reality.

```bash
# Windowed local pairing + per-position UNPAIRED (accessibility) probabilities
RNAplfold -W 200 -L 150 -u 30 < long_rna.fa      # -W window, -L max base-pair span, -u accessibility region length

# Scan for locally stable structured elements with a bounded span
RNALfold -L 150 < long_rna.fa

# Linear-time approximate MFE (LinearFold) and partition function (LinearPartition), if installed
echo "GGGAAACCC..." | linearfold
echo "GGGAAACCC..." | linearpartition
```

LinearFold's 5'->3' beam search can match or improve accuracy versus the exact cubic algorithm on long RNAs (the exact global model is not more correct when a single structure is not meaningful), besides being far faster.

## Pseudoknots: ViennaRNA cannot, by construction

The standard dynamic programming forbids crossing pairs, and general pseudoknot prediction is NP-hard (Lyngso & Pedersen 2000) -- RNAfold/RNAalifold silently return the best NESTED structure. Suspect a pseudoknot for tmRNA, telomerase RNA, RNase P, many riboswitch aptamers (SAM-II, preQ1), -1 ribosomal frameshift elements, IRES, and group I/II intron cores.

| Tool | Class / method | Note |
|------|----------------|------|
| IPknot | integer programming over base-pair probabilities | fast, broad class, the pragmatic default |
| ProbKnot | MEA assembly from McCaskill probabilities (RNAstructure) | any topology, fastest/most scalable |
| Knotty | MFE over the broad CCJ class | more complex crossing topologies |
| pknotsRG | MFE over restricted simple recursive pseudoknots, O(n^4) | narrower class |

Pseudoknot prediction is substantially less accurate and more expensive than nested folding -- treat any predicted pseudoknot as a hypothesis to corroborate with a second tool, covariation, or probing.

## RNA-RNA interaction: pick by the binding question

| Tool | Models | Use when |
|------|--------|----------|
| RNAcofold | both intramolecular AND intermolecular pairs; `-c` gives equilibrium concentrations | full dimerization model |
| RNAduplex | inter-molecular pairs only (no internal structure), fast | first-pass target screen |
| RNAup | opening (accessibility) energy + hybridization energy | sRNA/miRNA-target where the site may be buried in structure (the physically correct choice) |

The two strands are concatenated with `&` (RNAcofold/RNAduplex); RNAup takes the two sequences on separate lines.

```bash
# Full dimer model (intra + inter pairs); -p for the heterodimer partition function
echo "GCGCGCAUAU&AUAUGCGCGC" | RNAcofold -p --noPS
# With -c, RNAcofold reads the two monomer concentrations and reports equilibrium fractions of the
# five species (AB, AA, BB, A, B) -- use it to ask how much dimer actually forms, not just whether it is favorable.

# Fast inter-molecular-only hybridization screen (no internal structure)
echo "GCGCGCAUAU&AUAUGCGCGC" | RNAduplex

# Accessibility-corrected sRNA/miRNA-target binding: opening energy + hybridization (-b includes both)
RNAup -b < two_sequences.fa
```

## Thermodynamics vs deep learning: DL is not a default for novel RNA

Deep-learning predictors (SPOT-RNA, UFold, E2Efold) report high accuracy ON FAMILIES SEEN IN TRAINING, but under family-fold cross-validation that removes train/test homology their accuracy collapses to at or below the thermodynamic baseline (Szikszai et al. 2022); the apparent gains are intra-family memorization, and benchmark sets are ~55% rRNA / >90% rRNA+tRNA (Flamm et al. 2022). For a genuinely novel RNA (unseen Rfam family), no single-sequence method (DL or thermodynamic) is reliable -- the robust evidence is covariation (R-scape) and experimental probing. If using DL, prefer the thermodynamics-integrated hybrid MXfold2 over end-to-end nets; never cite intra-family accuracy as proof of de-novo performance.

"Is this more structured than random?" (z-score vs shuffled controls, RNAz, randfold): shuffles MUST preserve DINUCLEOTIDE composition (Altschul-Erikson), because MFE is dominated by GC content and base stacking, a dinucleotide property -- a mononucleotide shuffle inflates significance and makes almost anything look stable. For an alignment, RNAz combines a dinucleotide-controlled z-score with a structure conservation index (SCI = consensus MFE / mean single-sequence MFE; ~1 = a conserved structure), but reads Clustal/MAF, not Stockholm. A negative z-score means "more stable than random," NOT "this structure is correct"; covariation is the stronger evidence standard.

## Method-class selection (the big fork)

| Situation | Recommended | Avoid as default |
|---|---|---|
| Single novel RNA, no homologs | thermodynamic ensemble (RNAfold -p / RNAstructure) | pure end-to-end DL |
| Aligned homologs with covariation | RNAalifold + R-scape validation | single-sequence MFE |
| Homologs but no trusted alignment | TurboFold II (joint alignment + structure) | align-then-RNAalifold on a poor alignment |
| Have SHAPE/DMS reactivities | probing-restrained folding (Deigan/Zarringhalam) | unrestrained MFE |
| Long mRNA / transcriptome scale | RNAplfold / LinearFold+LinearPartition | global O(n^3) MFE |
| Pseudoknot biology | IPknot/ProbKnot/Knotty + cross-check | RNAfold (cannot) |
| Willing to use DL | MXfold2 (thermodynamics-integrated) | E2Efold/UFold on unseen families |

When competing methods or parameters are in play, verify current behavior against the installed tool's `--help` and the latest docs before trusting a number.

## Getting the structure out and drawing it

Dot-bracket is the default text form; CT and BPSEQ are the interchange formats downstream tools (ProbKnot, IPknot, RNAstructure) read and write. To DRAW a structure, use forna (web), R2DT (template-based standard layouts for known families), or VARNA; RNAfold's own `*_ss.ps` PostScript drawing is exactly what `--noPS` suppresses, so drop `--noPS` when the built-in diagram is wanted. The example renders the base-pair probability dot plot with matplotlib.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `AttributeError: module 'RNA' has no attribute 'sequence_shuffle'` | no such ViennaRNA function | use a dinucleotide-preserving shuffle (`ushuffle`, `esl-shuffle -d`) for z-scores |
| `bpp()`/`centroid()`/`MEA()` return empty or garbage | `fc.pf()` not called first | call `fc.pf()` before any ensemble quantity |
| `*_dp.ps` / `*_ss.ps` files appearing in the working directory | RNAfold/RNAalifold write PostScript by default | pass `--noPS` (and run in a scratch dir) |
| Long mRNA gives one improbable global fold | global MFE is meaningless past ~700 nt | use RNAplfold / LinearFold / LinearPartition |
| Predicted structure has a pseudoknot the tool "missed" | RNAfold cannot represent crossing pairs | use IPknot / ProbKnot / Knotty |
| Consensus structure looks confident but is wrong | RNAalifold trusts the alignment; no real covariation | validate with R-scape; check alignment quality |
| `RNAalifold --aln alignment.sto` treated as input flag | `--aln` is an OUTPUT (annotated PostScript) flag | pass the alignment positionally: `RNAalifold alignment.sto` |
| SHAPE-constrained fold barely changes | reactivity vector mis-indexed or zeros where data is missing | vector is 1-indexed (prepend -999); use -999 for no-data, not 0 |

## Related Skills

- structure-probing - Obtain SHAPE/DMS reactivities to constrain folding
- ncrna-search - Classify structured RNAs by family with Infernal/Rfam
- covariation-analysis - Statistically validate a predicted conserved structure with R-scape
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation
- small-rna-seq/target-prediction - miRNA-target prediction using accessibility
- sequence-manipulation/sequence-properties - Sequence composition and GC content
- data-visualization/heatmaps-clustering - Rendering the base-pair probability matrix (dot plot)

## References

- McCaskill JS. 1990. The equilibrium partition function and base pair binding probabilities for RNA secondary structure. Biopolymers 29(6-7):1105-1119. doi:10.1002/bip.360290621
- Mathews DH, Sabina J, Zuker M, Turner DH. 1999. Expanded sequence dependence of thermodynamic parameters improves prediction of RNA secondary structure. J Mol Biol 288(5):911-940. doi:10.1006/jmbi.1999.2700
- Ding Y, Chan CY, Lawrence CE. 2005. RNA secondary structure prediction by centroids in a Boltzmann weighted ensemble. RNA 11(8):1157-1166. doi:10.1261/rna.2500605
- Lu ZJ, Gloor JW, Mathews DH. 2009. Improved RNA secondary structure prediction by maximizing expected pair accuracy. RNA 15(10):1805-1813. doi:10.1261/rna.1643609
- Lorenz R, Bernhart SH, Honer zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. 2011. ViennaRNA Package 2.0. Algorithms Mol Biol 6:26. doi:10.1186/1748-7188-6-26
- Lyngso RB, Pedersen CNS. 2000. RNA pseudoknot prediction in energy-based models. J Comput Biol 7(3-4):409-427. doi:10.1089/106652700750050862
- Sato K, Kato Y, Hamada M, Akutsu T, Asai K. 2011. IPknot: fast and accurate prediction of RNA secondary structures with pseudoknots using integer programming. Bioinformatics 27(13):i85-i93. doi:10.1093/bioinformatics/btr215
- Bellaousov S, Mathews DH. 2010. ProbKnot: fast prediction of RNA secondary structure including pseudoknots. RNA 16(10):1870-1880. doi:10.1261/rna.2125310
- Jabbari H, Wark I, Montemagno C, Will S. 2018. Knotty: efficient and accurate prediction of complex RNA pseudoknot structures. Bioinformatics 34(22):3849-3856. doi:10.1093/bioinformatics/bty420
- Tan Z, Fu Y, Sharma G, Mathews DH. 2017. TurboFold II: RNA structural alignment and secondary structure prediction informed by multiple homologs. Nucleic Acids Res 45(20):11570-11581. doi:10.1093/nar/gkx815
- Huang L, Zhang H, Deng D, Zhao K, Liu K, Hendrix DA, Mathews DH. 2019. LinearFold: linear-time approximate RNA folding by 5'-to-3' dynamic programming and beam search. Bioinformatics 35(14):i295-i304. doi:10.1093/bioinformatics/btz375
- Zhang H, Zhang L, Mathews DH, Huang L. 2020. LinearPartition: linear-time approximation of RNA folding partition function and base-pairing probabilities. Bioinformatics 36(Suppl_1):i258-i267. doi:10.1093/bioinformatics/btaa460
- Hajdin CE, Bellaousov S, Huggins W, Leonard CW, Mathews DH, Weeks KM. 2013. Accurate SHAPE-directed RNA secondary structure modeling, including pseudoknots. Proc Natl Acad Sci USA 110(14):5498-5503. doi:10.1073/pnas.1219988110
- Sato K, Akiyama M, Sakakibara Y. 2021. RNA secondary structure prediction using deep learning with thermodynamic integration (MXfold2). Nat Commun 12:941. doi:10.1038/s41467-021-21194-4
- Szikszai M, Wise M, Datta A, Ward M, Mathews DH. 2022. Deep learning models for RNA secondary structure prediction (probably) do not generalise across families. Bioinformatics 38(16):3892-3899. doi:10.1093/bioinformatics/btac415
- Flamm C, Wielach J, Wolfinger MT, Badelt S, Lorenz R, Hofacker IL. 2022. Caveats to deep learning approaches to RNA secondary structure prediction. Front Bioinform 2:835422. doi:10.3389/fbinf.2022.835422
- Rivas E, Clements J, Eddy SR. 2017. A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs. Nat Methods 14(1):45-48. doi:10.1038/nmeth.4066
