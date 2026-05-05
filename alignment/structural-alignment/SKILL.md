---
name: bio-alignment-structural
description: Align protein structures using Foldseek 3Di, TM-align, US-align, DALI, or Foldmason for structural MSA. Predict, score, and superpose backbone coordinates when sequence identity is below the twilight zone or remote-homology detection is required. Use when sequence MSA fails (<25% identity), when the dark proteome is the target, when AlphaFoldDB / ESM Atlas search is needed, or when structural superposition is the goal.
tool_type: mixed
primary_tool: Foldseek
---

## Version Compatibility

Reference examples tested with: Foldseek 8+, TM-align 20220412+, US-align 20231222+, Foldmason 1+, BioPython 1.83+, pymol-open-source 3.0+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `foldseek --version`, `TMalign`, `USalign`, `foldmason --version`
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Structural Alignment

**"Align two protein structures"** → Compute backbone-aware superposition and a fold-similarity score (TM-score, RMSD, or LDDT).
- CLI pairwise: `TMalign A.pdb B.pdb`, `USalign A.pdb B.pdb`
- CLI search at scale: `foldseek easy-search query/ AFDB result.m8 tmp/`
- CLI structural MSA: `foldmason easy-msa structures/*.pdb out tmp/`
- Python pairwise: `Bio.PDB.Superimposer`, or `subprocess` wrapping `TMalign` / `USalign` (see `examples/tm_align_pairwise.py`)
- GUI / scripted molecular-graphics superposition: ChimeraX `matchmaker`, PyMOL `super`/`cealign`

**"Find structural homologs of an AlphaFold model"** → Search a structure database by 3Di-encoded structural alphabet (Foldseek) or by full TM-align rotation (DALI, US-align).

## When to Use Structural Alignment

| Sequence identity | Recommended approach |
|-------------------|---------------------|
| >= 40% | Sequence DP (Bio.Align, BLASTP) is sufficient |
| 25-40% | Sensitive sequence (MMseqs2, jackhmmer, profile-profile HHsearch) |
| 15-25% | Profile-profile (HHsearch) OR Foldseek if structures available |
| < 15% (dark proteome / twilight zone) | Foldseek (3Di), TM-align, US-align, pLM aligners |

Sequence alignment below 15% identity is statistically indistinguishable from random pairings. The exact twilight-zone cutoff is length-dependent: Rost 1999 (Prot Eng) showed the curve drops to 25% at length 80, 20% at length 250 -- short alignments need higher identity for the same statistical signal, so a 15-25% rule of thumb is shorthand for "twilight zone for proteins of typical domain size (~150-300 residues)". If reasonable structural models exist (PDB, AlphaFoldDB, ESMFold), structural alignment is far more reliable in this regime.

### Twilight-Zone Threshold Exceptions

For families with strong functional or structural constraint (ribosomal proteins, class I aminoacyl-tRNA synthetases with HIGH/KMSKS motifs, histone fold, cytochrome c with CXXCH, or any long alignment with well-distributed conservation), sequence MSA may remain reliable down to ~12-15% identity. Verify with a profile-profile method (HHsearch) before committing to structural alignment.

## Pairwise Structural Alignment Tool Selection

| Tool | Reference | Score | Best for |
|------|-----------|-------|----------|
| TM-align | Zhang & Skolnick 2005 NAR | TM-score | Single-chain pairwise; standard fold-similarity benchmark |
| US-align | Zhang et al 2022 Nat Methods | TM-score | Multi-chain protein, RNA, DNA, complexes; original successor to TM-align |
| Foldseek-Multimer | Kim et al 2025 Nat Methods | TM-score | Database-scale multi-chain complex search; 10-100x faster than US-align in pairwise mode, 1000-10000x in DB search; preferred at AFDB-multimer scale |
| FATCAT | Ye & Godzik 2003 Bioinf | FATCAT score | Flexible alignment with allowed twists/breaks |
| CE | Shindyalov & Bourne 1998 Prot Eng | CE score | Combinatorial extension of fragments; PDB legacy |
| DALI | Holm 2022 NAR | Z-score | Distance-matrix alignment; superior at TM 0.3-0.5 (twilight fold); operates on AFDB at scale |
| Bio.PDB.Superimposer | Cock et al 2009 Bioinf | RMSD | Pure-Python superposition with known atom correspondence |

TM-align and US-align report two TM-scores by default: one normalised by chain 1 length and one by chain 2 length. The standard fold-similarity metric for asymmetric reporting is the larger of the two (effectively normalising by the shorter chain), since TM-score is length-asymmetric. The `-a T` flag adds a third score normalised by the average of the two chain lengths (useful for symmetric clustering). **TM > 0.5** indicates the same fold; **TM > 0.8** indicates equivalent topology / close structural relationship; **TM < 0.2** indicates random structural similarity. RMSD alone is misleading: RMSD scales with length, depends on outliers, and an "optimization RMSD" (used to fit) is not the same as an "evaluation RMSD" (used to compare). Always report TM-score alongside RMSD and the alignment length.

| Metric | Threshold | Source / Interpretation |
|--------|-----------|-------------------------|
| TM-score | > 0.5 | Same fold (Zhang & Skolnick 2004 Proteins; reported by TM-align, US-align, Foldseek, DALI) |
| TM-score | > 0.8 | Equivalent topology (homologous) |
| TM-score | < 0.2 | Random structural similarity |
| DALI Z-score | > 20 | Definitely homologous (Holm 2020 Prot Sci) |
| DALI Z-score | 8 - 19 | Probable homology |
| DALI Z-score | 2 - 8 | Candidate; verify with TM-score or biology |
| DALI Z-score | < 2 | Not significant (random) |
| GDT-TS | > 50 | Correct fold (CASP/LGA scoring; reported by LGA, MaxCluster, OpenStructure -- not by TM-align/Foldseek) |
| LDDT | > 0.6 | Conventional cutoff for "correctly modelled" residue (Mariani et al 2013 Bioinf shows the score is fold-architecture-dependent and does not define a universal hard threshold; >= 0.6 is widely used in CAMEO and Foldseek output as a working cutoff) |
| RMSD | < 2 A over >100 residues | Strong superposition; below 1.5 A is excellent |

### Foldseek vs DALI: Modern Comparison

Both target structural homolog search, but with different strengths:

| Question | Foldseek | DALI |
|----------|----------|------|
| Find any same-fold homolog quickly | YES (3Di indexed; 1000-1M seq/s) | Slow (full distance matrix) |
| Sensitivity at TM > 0.5 (same-fold) | Comparable to TM-align | Slightly higher recall |
| Sensitivity at TM 0.3-0.5 (twilight fold) | Recall drops sharply | DALI Z-score retains signal (Holm 2022) |
| Alignment quality for indel-rich pairs | Local 3Di+AA; can miss large indels | Distance-matrix optimal handles indels well |
| AFDB-scale all-vs-all | Tractable | Now tractable post-2022 (Holm DALI server update) |

Practical workflow for remote homology: (1) Foldseek `easy-search` to retrieve same-fold candidates fast, (2) for hits with TM < 0.5 or with structural-indel rich alignments, re-align with DALI for higher-quality residue equivalences. The two tools are complementary; Foldseek-only misses some twilight-fold relationships, DALI-only misses the AFDB-scale opportunity.

### TM-score Threshold Caveats

Apply the 0.5 same-fold rule only to globular domains of ~100-300 residues; for chains <60 residues TM > 0.5 occurs by chance ~3-5%, so use Xu & Zhang 2010's length-aware Gumbel p-value instead. Always report both length-normalised TM-scores (or use `-a T` for the symmetrised average) and pair RMSD with the atom count (`(RMSD, n_atoms_used)`); raw RMSD without length normalisation is misleading.

### TM-align Pairwise Run

**Goal:** Compute TM-score and RMSD between two structures and write a superposed PDB.

**Approach:** Invoke TMalign or USalign with output flags; parse the structured output for downstream filtering.

```bash
TMalign chainA.pdb chainB.pdb -o superposed.sup
TMalign chainA.pdb chainB.pdb -outfmt 2

USalign chainA.pdb chainB.pdb -mol prot -outfmt 2
USalign complex_A.pdb complex_B.pdb -mm 1 -ter 0
```

`-outfmt 2` returns tabular output (one line per pair). For multi-chain complexes, US-align with `-mm 1 -ter 0` aligns full assemblies; TM-align is single-chain only.

**US-align multi-chain score interpretation.** US-align `-mm 1 -ter 0` produces multiple TM-scores per complex: normalised by query length, normalised by template length, and (for symmetric homo-oligomers) the best-matching chain permutation. For different stoichiometries (query hetero-A2B2 vs template hetero-A4B4), use `-mm 4` and `-byresi 0`. The "best-of" TM-score is the appropriate cross-complex homology metric; per-chain TM-scores serve as a sanity check (low per-chain + high complex TM = topology match without local fold conservation, suggests promiscuous interaction not deep homology). See Zhang et al 2022 Nat Methods for the full mode taxonomy.

### Foldseek-Multimer for Database-Scale Complex Search

For multi-chain complex search at AFDB-multimer scale (~214k entries) or against PDB complexes, US-align is too slow; Foldseek-Multimer (Kim et al 2025 Nat Methods) is the modern default. Two modes share the chain-pairing prefilter:

| Mode | Algorithm | Use when |
|------|-----------|----------|
| `Foldseek-MM` (default) | 3Di+AA Gotoh per chain pair | Fast database search; default speed-sensitivity tradeoff |
| `Foldseek-MM-TM` | TM-align per chain pair after prefilter | Top-hit refinement with full TM-score |

```bash
foldseek easy-multimersearch query_complex.pdb afdb_multimer result tmp/
foldseek easy-multimersearch query_complex.pdb afdb_multimer result tmp/ --multimer-tm-threshold 0.5

foldseek easy-multimercluster *.pdb cluster_result tmp/ --multimer-tm-threshold 0.65
```

Decision guide: pairwise complex pair on a few hundred targets -> US-align with `-mm 1 -ter 0`. Database search across thousands to millions of complex entries -> Foldseek-Multimer. For AFDB-Multimer or PDB100-Multimer scale, US-align is computationally infeasible; Foldseek-Multimer matches US-align chain-pairing in >99% of cases at 10-100x speedup pairwise and 10^3-10^4x in database mode (Kim et al 2025 benchmark).

**Reporting convention:** Always quote BOTH the multimer TM-score (whole-complex) AND the worst per-chain TM-score; high multimer TM with one low per-chain score signals topology-matched but locally divergent chains (often promiscuous binders or paralog swaps), which is biologically distinct from a uniformly-conserved complex.

### Bio.PDB.Superimposer

**Goal:** Superpose a known atom correspondence and compute RMSD without running an external aligner.

**Approach:** Use when residue equivalence is already established (e.g. same sequence, different conformations). For unknown correspondence, prefer TM-align / US-align.

```python
from Bio.PDB import PDBParser, Superimposer

parser = PDBParser(QUIET=True)
mobile_atoms = list(parser.get_structure('mobile', 'mobile.pdb').get_atoms())
reference_atoms = list(parser.get_structure('ref', 'reference.pdb').get_atoms())

ca_mobile = [a for a in mobile_atoms if a.get_id() == 'CA']
ca_reference = [a for a in reference_atoms if a.get_id() == 'CA']
n = min(len(ca_mobile), len(ca_reference))

sup = Superimposer()
sup.set_atoms(ca_reference[:n], ca_mobile[:n])
sup.apply(mobile_atoms)
print(f'RMSD: {sup.rms:.3f} A over {n} CA atoms')
```

## Structural Search at Scale: Foldseek

Run Foldseek for AlphaFoldDB-scale structural search; it indexes a 20-letter 3Di alphabet for thousand- to million-fold speedup over TM-align at comparable same-fold sensitivity.

```bash
# Search query structures against AFDB (default: --alignment-type 2)
foldseek easy-search query.pdb afdb_database result.m8 tmp/

# Refine top hits with full TM-align rotation (slower but global TM-score)
foldseek easy-search query.pdb afdb_database result.m8 tmp/ --alignment-type 1

# All-versus-all clustering at TM > 0.5
foldseek easy-cluster structures/*.pdb cluster_result tmp/ --tmscore-threshold 0.5

# Custom output columns
foldseek easy-search query.pdb afdb_database result.m8 tmp/ \
    --format-output query,target,evalue,alntmscore,qtmscore,ttmscore,lddt,bits
```

| Foldseek `--alignment-type` | Algorithm | Use when |
|-----------------------------|-----------|----------|
| 0 | 3Di Gotoh-Smith-Waterman (local) | Not recommended; 3Di alone, no amino-acid signal |
| 1 | TMalign (global) | Refine top hits with full TM-score; slow |
| 2 | 3Di+AA Gotoh-Smith-Waterman (local) | Default; best speed-sensitivity tradeoff |

### pLDDT-Filtering AlphaFold Structures Before Foldseek

Mask residues with pLDDT < 70 before Foldseek indexing or search; their backbone coordinates encode as random 3Di letters and contaminate hits below TM ~ 0.4. AFDB clusters from Barrio-Hernandez et al 2023 are pre-filtered; ESMFold predictions need the same step.

```bash
# Mask low-pLDDT residues during database creation (B-factor column stores pLDDT)
# The *.pdb glob requires shell expansion; for Python subprocess calls use glob.glob()
# to expand the file list before passing as argv.
foldseek createdb --mask-bfactor-threshold 70.0 *.pdb afdb_masked
```

## Structural Multiple Sequence Alignment

| Tool | Reference | Best for |
|------|-----------|----------|
| Foldmason | Gilchrist et al 2024 (preprint) | Billion-protein-scale structural MSA on AFDB; the structural counterpart to MAFFT |
| 3D-Coffee / Expresso | Notredame group | <100 chains with mixed PDB and sequence input |
| mTM-align | Dong et al 2018 Bioinf | Multiple structure alignment; output suitable for tree building. Server: `yanglab.qd.sdu.edu.cn/mTM-align/` (Yang Lab moved from Nankai to Shandong University in 2021; verify URL before use) |
| PROMALS3D | Pei, Kim & Grishin 2008 NAR | Hybrid sequence-structure profile MSA; uses PSI-BLAST + DALI/SAP templates |
| MUSTANG | Konagurthu et al 2006 Proteins | Pure structural MSA via residue-residue equivalences |

### Foldmason easy-msa

```bash
foldmason easy-msa structures/*.pdb result tmp/ \
    --refine-iters 100 \
    --report-mode 1

# Outputs: result_aa.fa (amino-acid MSA), result_3di.fa (3Di MSA), result.nw (guide tree), result.html (LDDT report)
```

`--refine-iters` controls iterative MSA refinement; `--report-mode 1` produces an HTML report with per-column LDDT confidence. The 3Di MSA can be used directly for evolutionary analyses where structure rather than sequence is the appropriate signal.

**Per-column LDDT extraction:** Run with `--report-mode 2` to produce machine-readable JSON output:

```bash
foldmason easy-msa structures/*.pdb result tmp/ --report-mode 2
# Produces result.json with per-column LDDT in machine-readable form
```

```python
import json
with open('result.json') as f:
    report = json.load(f)
lddt_per_column = report.get('per_column_lddt')
```

Verify the exact JSON schema with `foldmason easy-msa --help` and inspect a sample `result.json`; the field name may evolve across releases.

## Hybrid Sequence-Structure Approaches

When a small dataset has known structures or templates, hybrid tools dramatically improve MSA accuracy:

- **T-Coffee Expresso** (Notredame et al) -- PSI-BLAST searches PDB for templates, runs SAP or TM-align pairwise structural library, then T-Coffee consistency over the resulting library. Reported in Notredame-group benchmarks to substantially improve correct-column rate over sequence-only T-Coffee. Requires internet access for PSI-BLAST and PDB template fetching; for offline / reproducible runs use 3D-Coffee with locally curated templates.
- **3D-Coffee** -- user-supplied PDB templates with `-template_file templates.txt`; uses SAP or TM-align as the structural method.
- **PROMALS3D** (Pei, Kim & Grishin 2008 NAR) -- PSI-BLAST profiles + secondary-structure prediction + DALI/SAP templates; used heavily in the Grishin lab's superfamily annotations. PROMALS3D server availability has been intermittent; if unreachable, T-Coffee Expresso (mode `expresso`) is a current alternative for structure-informed sequence MSA.

```bash
t_coffee input.fasta -mode expresso -output fasta_aln -outfile aligned.fasta

t_coffee input.fasta -template_file templates.txt -mode 3dcoffee
```

## AlphaFold Integration

The "predict-then-align" workflow has become standard for remote-homology problems:

1. Search sequence with MMseqs2 / jackhmmer to seed an MSA.
2. Predict structure with AlphaFold2 (ColabFold), ESMFold, or AlphaFold3 -- see `structural-biology/modern-structure-prediction` for prediction workflows.
3. Search the predicted structure against AFDB or PDB with Foldseek.
4. Use Foldmason or PROMALS3D to derive a structure-aware MSA from the hits.
5. Refine sequence MSA using the structure-derived column equivalences.

Search AFDB clusters (~2.3 M Foldseek-derived clusters from Barrio-Hernandez et al 2023) as the canonical remote-homology entry point; supplement with the ESM Atlas (~600 M ESMFold metagenomic structures) when AFDB recall is insufficient. For pLDDT semantics and curated AlphaFoldDB entry handling, see `structural-biology/alphafold-predictions`.

## pLM-Based Sequence Aligners

Run a pLM aligner (TM-Vec, vcMSA, DEDAL, pLM-BLAST) when no structure is available but identity is in the twilight zone (<15-25%). They recover signal that DP aligners miss but do NOT replace TM-align / US-align when structures exist; use them as a complementary first pass.

| Tool | Reference | Embedding source |
|------|-----------|------------------|
| vcMSA | McWhite, Armour-Garb & Singh 2023 Genome Res | ProtT5 vector clustering for MSA (alpha-stage tool per its repo README; last release Oct 2023; test on representative inputs before pipeline use) |
| DEDAL | Llinares-Lopez et al 2023 Bioinf | Differentiable end-to-end alignment with pLM features |
| TM-Vec | Hamamsy et al 2024 Nat Biotech | TM-score prediction from ProtT5 embeddings (active fork: `valentynbez/tmvec`; original `tymor22/tm-vec` is in limited maintenance) |
| pLM-BLAST | Kaminski et al 2023 Bioinf | BLAST-style hits via pLM cosine similarity |

These run on raw sequence (no structure required) and recover signal at <15% identity that DP aligners miss entirely. They do NOT replace structural alignment when structures exist; treat them as complementary.

**pLM aligner accuracy ceiling.** TM-Vec, vcMSA, and DEDAL predict alignment-equivalent quantities (TM-score, alignment columns) from sequence-only embeddings. Hamamsy et al 2024 report TM-Vec TM-score prediction RMSE ~0.07 on SCOP -- sufficient for coarse filtering ("is this hit at all structurally similar?") but NOT for fine-grained ranking (e.g. choosing among hits with TM 0.5-0.7). For final structural-similarity scoring, run TM-align or US-align on predicted structures (ColabFold + USalign) rather than relying on the pLM-predicted score.

## Visualisation and Inspection

Structural alignments are read by viewing the superposed structures, not the sequence text:

```bash
# PyMOL headless: -cq runs without GUI/quietly; -d passes commands directly
pymol -cq -d "load reference.pdb; load mobile.pdb; super mobile, reference; ray 800,600; png fig.png"

# ChimeraX (modern, replaces Chimera): MatchMaker uses Needleman-Wunsch + iterative refinement
ChimeraX --cmd "open reference.pdb mobile.pdb; matchmaker #2 to #1; save fig.png"
```

PyMOL `super` performs cycle-fitting for distantly related structures (better than `align` when sequence identity is low); `cealign` runs CE; `tmalign` (PyMOL plugin) wraps TM-align. ChimeraX MatchMaker is the modern replacement and uses Needleman-Wunsch + iterative refinement with the option to invoke external tools.

## Decision Tree by Goal

| Goal | First-line tool |
|------|-----------------|
| Two structures, known correspondence | `Bio.PDB.Superimposer` |
| Two structures, unknown correspondence | TMalign or USalign |
| Multi-chain complex, pairwise | USalign with `-mm 1 -ter 0` |
| Multi-chain complex, database search | `foldseek easy-multimersearch` (Foldseek-Multimer) |
| Twilight-fold homology (TM 0.3-0.5) | DALI (Z-score ranks low-similarity hits better than Foldseek) |
| Structural homolog search at AFDB scale (single chain) | `foldseek easy-search` |
| All-vs-all clustering of structures | `foldseek easy-cluster --tmscore-threshold 0.5` |
| Multiple structure alignment, < 100 chains | T-Coffee Expresso or mTM-align |
| Multiple structure alignment, > 1000 chains | Foldmason `easy-msa` |
| Hybrid sequence-structure MSA | PROMALS3D or T-Coffee Expresso |
| pLM aligner for dark proteome | TM-Vec, vcMSA, DEDAL |
| Distant homology with no structures | Predict with ColabFold/ESMFold first, then Foldseek |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| TM-align "atoms not enough" | Single-residue chains or ligand-only PDB | Filter to canonical amino acids before running |
| Foldseek "no hits" | Wrong database format | Run `foldseek databases` to confirm AFDB / PDB100 download |
| Bio.PDB Superimposer "different number of atoms" | Atom selection mismatch | Filter both lists to common atom names (CA only, or N/CA/C/O) |
| TM-score normalised by length 1 | Used `-L 1` | Drop `-L` flag; default normalisation is correct |
| Foldmason output empty | Input structures not in same chain ID | Renumber and rename chains uniformly before running |

## Related Skills

- alignment/multiple-alignment - Sequence MSA when identity > 25%; complementary for hybrid tools
- alignment/pairwise-alignment - Sequence pairwise; use first to filter before structural alignment
- alignment/msa-parsing - Parse and analyze structural MSA output for downstream metrics
- alignment/msa-statistics - Apply per-column conservation to structurally-derived MSAs
- alignment/alignment-io - Read and convert structural MSA files for downstream tools
- alignment/alignment-trimming - Trim structural MSAs the same way as sequence MSAs
- structural-biology/modern-structure-prediction - Predict structures (AlphaFold2/3, ESMFold) used as input
- structural-biology/alphafold-predictions - Curated AFDB / pLDDT handling
- structural-biology/structure-navigation - Map alignment columns to PDB residues
- phylogenetics/modern-tree-inference - Trees from structural MSAs (Foldmason output works directly)

## References

- van Kempen M et al. 2024. Fast and accurate protein structure search with Foldseek. Nat Biotech 42:243-246.
- Kim W, Mirdita M, Levy Karin E, Gilchrist CLM, Schweke H, Soding J, Levy E, Steinegger M. 2025. Rapid and sensitive protein complex alignment with Foldseek-Multimer. Nat Methods (in press; Foldseek-Multimer release).
- Zhang Y, Skolnick J. 2005. TM-align: a protein structure alignment algorithm based on the TM-score. NAR 33:2302-2309.
- Zhang C, Shine M, Pyle AM, Zhang Y. 2022. US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes. Nat Methods 19:1109-1115.
- Holm L. 2020. Using Dali for protein structure comparison. Methods Mol Biol 2112:29-42 (Z-score interpretation thresholds).
- Holm L. 2022. Dali server: structural unification of protein families. NAR 50:W210-W215.
- Gilchrist CLM et al. 2024. Foldmason: multiple protein structure alignment at scale with 3Di. Preprint.
- Barrio-Hernandez I et al. 2023. Clustering predicted structures at the scale of the known protein universe. Nature 622:637-645.
- Hamamsy T et al. 2024. Protein remote homology detection and structural alignment using deep learning. Nat Biotech 42:975-985.
- Xu J, Zhang Y. 2010. How significant is a protein structure similarity with TM-score = 0.5? Bioinf 26:889-895.
