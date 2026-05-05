---
name: bio-alignment-multiple
description: Perform multiple sequence alignment using MAFFT, MUSCLE5, ClustalOmega, or T-Coffee. Guides tool and algorithm selection based on dataset size, sequence divergence, and downstream application. Use when aligning three or more homologous sequences for phylogenetics, conservation analysis, or evolutionary studies.
tool_type: mixed
primary_tool: MAFFT
---

## Version Compatibility

Reference examples tested with: MAFFT 7.520+, MUSCLE 5.1+, ClustalOmega 1.2.4+, T-Coffee 13+, PAL2NAL 14+, BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mafft --version`, `muscle -version`, `clustalo --version`
- Python: `pip show biopython` then `help(module.function)` to check signatures

If code throws errors, introspect the installed tool and adapt the example to match the actual CLI flags rather than retrying.

# Multiple Sequence Alignment

**"Align multiple sequences"** → Compute an optimal alignment of three or more homologous sequences using progressive, iterative, or consistency-based methods.
- CLI: `mafft` (most versatile), `muscle` (highest accuracy), `clustalo` (scales well), `t_coffee` (consistency-based)
- Python: `subprocess.run()` wrapping CLI tools; BioPython `Bio.Align.Applications` was removed in BioPython 1.86 (verify with `pip show biopython`); use `subprocess` directly

## MSA Algorithm Taxonomy

When a tool is failing on a dataset, switch to a tool from a different algorithmic family rather than tuning flags. The six families and their characteristic failure modes:

| Family | Representative tools | Best at | Fails when |
|--------|---------------------|---------|-----------|
| Progressive | ClustalW, MAFFT FFT-NS-2 | Fast, large datasets, similar lengths | Early-stage gap errors propagate; no recovery |
| Iterative refinement | MAFFT L-INS-i, MUSCLE3, PRRN | Recovers from progressive errors at <2000 seqs | Slow on >2000; still guide-tree dependent |
| Consistency-based | T-Coffee, ProbCons | Highest accuracy <100 seqs; integrates evidence | O(N^2 to N^4) scaling; heavy compute |
| HMM-based | HMMER hmmalign, ClustalOmega (HHalign), UPP, WITCH | Adding sequences to a curated profile; fragmentary input | Needs an existing high-quality profile or backbone |
| Divide-and-conquer | PASTA, MAGUS, MUSCLE5 super5 | Heterogeneous large datasets (>10k seqs) | Sub-alignment merges can introduce artefacts |
| Structure or pLM-informed | Foldmason, PROMALS3D, vcMSA, 3D-Coffee | Dark proteome, <15% identity, dataset has structures | Requires structures or a working pLM |

## Tool Selection

Pick by dataset size and divergence; the default recommendations follow the table.

| Tool | Best For | Max Sequences | Accuracy | Speed |
|------|----------|---------------|----------|-------|
| MAFFT L-INS-i | Highest accuracy, <200 seqs | ~200 | Highest | Slow |
| MAFFT FFT-NS-2 | Large datasets, good balance | ~50,000 | Good | Fast |
| MAFFT E-INS-i | Sequences with long unalignable internal regions | ~200 | High | Slow |
| MUSCLE5 (PPP `-align`) | Benchmarked highest accuracy on Balifam-10000 | ~1000 | Highest | Medium |
| MUSCLE5 (`-super5`) | Large datasets via mBed clustering | ~100,000+ | Good | Medium |
| ClustalOmega | Very large datasets, HMM-based profiles | ~190,000 in published benchmark (Sievers et al 2011 Mol Syst Biol) | Good | Fast |
| T-Coffee (default) | Small datasets needing maximum accuracy | ~200 | Highest | Slowest |

**Default recommendation**: MAFFT L-INS-i for <200 sequences; MAFFT FFT-NS-2 or MUSCLE5 super5 for thousands; MUSCLE5 ensemble (`-stratified`) when alignment confidence estimates are needed.

### Cross-Aligner Sensitivity Check

When downstream analysis depends on a specific column (a candidate selection site, a contact-prediction position), run BOTH MAFFT L-INS-i and MUSCLE5 `-align` and verify that column is stable across the two outputs. If unstable, flag as low-confidence regardless of GUIDANCE2/TCS scores. The MUSCLE5 ensemble (`-stratified`/`-diversified`) accomplishes the same check directly within one tool and is preferred when ensemble output is acceptable downstream.

### Beyond MAFFT and MUSCLE: Scale and Domain

Some workloads exceed what the four main tools handle gracefully. Use the scale-and-domain table below to escape default-tool blind spots.

| Scenario | Recommended tool | Why |
|----------|------------------|-----|
| 100k - millions of sequences (UniRef cluster reps) | FAMSA v2 (Deorowicz et al 2016 SREP, v2 2024) | Million-scale MSA in hours with Pareto-optimal accuracy |
| Heterogeneous large dataset (variable length, divergence) | PASTA (Mirarab et al 2015 J Comp Biol) or MAGUS (Smirnov & Warnow 2021 Bioinf) | Divide-and-conquer plus iterative re-alignment |
| Fragmentary input (metagenomics, eDNA, ancient DNA) | UPP (Nguyen et al 2015 Genome Biol) or WITCH (Shen et al 2022 J Comp Biol) | HMM backbone tolerates partial sequences |
| Dark proteome / <15% identity | vcMSA (McWhite, Armour-Garb & Singh 2023 Genome Res; alpha-stage tool per its repo README; last release Oct 2023; test on representative inputs before pipeline use) or Foldmason (Gilchrist et al 2024) | pLM-embedding or structural alignment where sequence fails |
| RNA with secondary structure | Infernal `cmalign` (Nawrocki & Eddy 2013), R-Coffee, MAFFT-Q-INS-i | Sequence + base-pair consistency |
| Strand-unknown (Sanger, Nanopore raw) | `mafft --adjustdirection --globalpair` | Auto-detects and reverse-complements as needed |
| Adding many sequences to a curated profile | HMMER `hmmalign` (Eddy 2011 PLOS CB) | Profile-driven; better than `mafft --add` for Pfam-style use |

**vcMSA limitation:** Pre-filter input to within ~2x mean length; vcMSA degrades sharply on mixed full-length/fragment input because ProtT5 embeddings encode positional context. For mixed-length sets, segment long sequences via HHsearch domain decomposition before alignment, or use Foldmason on predicted structures.

## Critical Concepts

### Mitigate Guide-Tree Dependency

All major MSA tools build a guide tree, then align progressively along it; once a gap is inserted in the progressive phase, it is never removed, so early errors propagate. To mitigate: prefer iterative-refinement modes (MAFFT `-i`, MUSCLE5), use consistency scoring (T-Coffee) for small datasets, and quantify uncertainty with GUIDANCE2 bootstrapping or the MUSCLE5 ensemble before publishing column-specific conclusions.

### Joint MSA-Phylogeny Co-estimation (Small Datasets Only)

The theoretically correct answer to guide-tree dependency is to estimate the alignment and tree jointly under a statistical evolutionary model rather than treating MSA as a fixed input to phylogenetic inference. BAli-Phy version 3 (Redelings 2021 Bioinf 37:3032) does this via MCMC, producing a posterior distribution over alignments and trees with insertion/deletion rates as model parameters. Version 3 is O(n) instead of O(n^2) per likelihood evaluation, but the practical ceiling remains ~70-200 sequences before runtime becomes prohibitive (weeks for >100 sequences). When the dataset fits the cap, BAli-Phy gives the most defensible alignment+tree pair for publication; when it does not, the practical alternative is MUSCLE5 ensemble + IQ-TREE per-replicate (Edgar 2022) to approximate posterior alignment uncertainty without joint estimation.

| Dataset size | Recommended approach |
|--------------|----------------------|
| < 70 sequences | BAli-Phy v3 joint MSA+tree posterior; gold standard |
| 70 - 200 sequences | BAli-Phy v3 if compute allows (weeks); else MUSCLE5 ensemble + IQ-TREE per replicate |
| > 200 sequences | MUSCLE5 ensemble (`-stratified`) + IQ-TREE per replicate; BAli-Phy not feasible |
| > 1000 sequences | Single MAFFT/MUSCLE5 + standard bootstrap; ensemble methods become intractable |

### Sequence Divergence Thresholds

| Protein Identity | Signal Level | Recommendation |
|-----------------|--------------|----------------|
| >40% | Strong | Any MSA tool produces reliable alignment |
| 25-40% | Moderate (twilight zone begins) | Use iterative methods (L-INS-i, MUSCLE5); validate with GUIDANCE2 |
| 20-25% | Weak | Profile-profile methods (HHpred); consider structural alignment |
| <15-20% (length-dependent twilight) | Noise dominates signal | Sequence MSA is unreliable; switch to structural alignment (Foldseek, TM-align) or pLM aligners -- see `alignment/structural-alignment` |

## Running MAFFT

### Algorithm Selection

MAFFT offers multiple algorithms with explicit accuracy/speed tradeoffs. Selecting the right mode is critical; the difference between L-INS-i and FFT-NS-1 can be the difference between a correct and incorrect downstream phylogeny.

| Algorithm | Flag | Strategy | Best For |
|-----------|------|----------|----------|
| FFT-NS-1 | `--retree 1` | Progressive only | Quick look, >10,000 seqs |
| FFT-NS-2 | `--retree 2` | Progressive + guide tree rebuild | Default balance, 200-10,000 seqs |
| FFT-NS-i | `--maxiterate 1000` | Iterative refinement | Moderate improvement, 200-2,000 seqs |
| G-INS-i | `--globalpair --maxiterate 1000` | Global pairwise + iterative | Sequences alignable over full length, <200 |
| L-INS-i | `--localpair --maxiterate 1000` | Local pairwise + iterative | Single alignable domain amid divergent flanks, <200 |
| E-INS-i | `--genafpair --maxiterate 1000` | Local with generalized affine gaps | Multiple conserved motifs separated by unalignable regions, <200 |
| Auto | `--auto` | Auto-selects based on dataset size | When unsure |

**Decision guide**: If sequences share a single conserved domain (most common case), use L-INS-i. If sequences are globally similar (e.g., ortholog set of similar length), use G-INS-i. If sequences have multiple conserved blocks separated by highly variable linker regions (e.g., multi-domain proteins with variable interdomain regions), use E-INS-i.

**G-INS-i failure mode:** When sequences have long divergent N/C-terminal extensions (signal peptides, intrinsically-disordered regions, isoform-specific tails), the global pairwise distance G-INS-i uses for guide-tree construction is dominated by the extension's mismatch content. The guide tree then mis-clusters sequences by extension similarity rather than core homology. Symptom: the resulting alignment has core conserved domains poorly aligned despite high-identity flanks. Switch to L-INS-i (local pairwise; tolerant of divergent flanks) or pre-process to remove signal peptides/disordered regions before alignment.

### What `--auto` Picks (and Why to Specify Explicitly)

`mafft --auto` silently downgrades the algorithm based on dataset size. The decision tree is roughly:

| Sequences | `--auto` selects | Equivalent flags |
|-----------|------------------|------------------|
| < 200 | L-INS-i | `--localpair --maxiterate 1000` |
| 200 - 500 | FFT-NS-i | `--retree 2 --maxiterate 2` |
| 500 - 2000 | FFT-NS-2 | `--retree 2 --maxiterate 0` |
| 2000 - 50000 | FFT-NS-2 (one-pass) | `--retree 1 --maxiterate 0` |
| > 50000 | PartTree | `--parttree --retree 1 --maxiterate 0` |

The transition at 200 sequences flips the alignment from "iterative-refined accurate" to "single-pass progressive". Note that `--auto` invokes FFT-NS-i with only `--maxiterate 2` in the 200-500 range (a truncated form of the full FFT-NS-i which uses `--maxiterate 1000`); for best accuracy in this range, specify `--retree 2 --maxiterate 1000` explicitly. For publication-quality phylogenetics, specify the algorithm explicitly so reproducibility audits do not rely on internal threshold heuristics.

### Basic Usage

**Goal:** Run MAFFT on a FASTA file with appropriate algorithm selection.

**Approach:** Invoke MAFFT via command line or subprocess, selecting the algorithm based on dataset characteristics.

```bash
# Highest accuracy for <200 sequences (local pairwise iterative)
mafft --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Good balance for medium datasets
mafft --retree 2 input.fasta > aligned.fasta

# Auto-select algorithm based on dataset size
mafft --auto input.fasta > aligned.fasta

# Protein alignment with specific matrix (default BLOSUM62)
mafft --amino --localpair --maxiterate 1000 input.fasta > aligned.fasta

# DNA alignment (auto-detected, but can be explicit)
mafft --nuc --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Adjust gap penalties (op=gap open, ep=gap extension)
mafft --op 1.53 --ep 0.123 --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Multithreaded
mafft --thread 8 --localpair --maxiterate 1000 input.fasta > aligned.fasta
```

```python
import subprocess

def run_mafft(input_fasta, output_fasta, algorithm='linsi', threads=4):
    algo_flags = {
        'linsi': ['--localpair', '--maxiterate', '1000'],
        'ginsi': ['--globalpair', '--maxiterate', '1000'],
        'einsi': ['--genafpair', '--maxiterate', '1000'],
        'fftns2': ['--retree', '2'],
        'auto': ['--auto'],
    }
    cmd = ['mafft', '--thread', str(threads)] + algo_flags[algorithm] + [input_fasta]
    with open(output_fasta, 'w') as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'MAFFT failed (exit {result.returncode}):\n{result.stderr}')

run_mafft('sequences.fasta', 'aligned.fasta', algorithm='linsi')
```

MAFFT writes its progress log and errors to stderr, not stdout. Capturing stderr and surfacing the failure message is essential when MAFFT exits non-zero (e.g. encountering ambiguous characters, oversized input, missing libraries) -- otherwise `check=True` raises a CalledProcessError without showing the actionable message.

### Adding Sequences to an Existing Alignment

**Goal:** Add new sequences to an existing MSA without realigning the entire dataset.

**Approach:** Use MAFFT's `--add` for full-length sequences, `--addfragments` for partial / surveillance / metagenomic reads. The two flags are NOT interchangeable.

| Flag | Use when | Behaviour |
|------|----------|-----------|
| `--add` | New sequences are full-length homologues | Each new sequence is profile-aligned end-to-end |
| `--addfragments` | New sequences are partial (reads, contigs, surveillance amplicons) | Free terminal gaps; new sequences may align to a sub-region of the profile |
| `--addprofile` | Adding an entire pre-aligned profile | Profile-profile alignment |
| `--keeplength` | Output must have the same column count as input MSA | Insertions in new sequences are dropped (key for HMMER/Pfam-style use) |

```bash
# Full-length additions; may extend alignment with new columns
mafft --add new_seqs.fasta existing_alignment.fasta > updated.fasta

# Partial reads or fragments; common for SARS-CoV-2 surveillance
mafft --addfragments new_reads.fasta --keeplength reference_msa.fasta > updated.fasta

# Profile-profile merge
mafft --addprofile second_msa.fasta first_msa.fasta > merged.fasta
```

For HMM-curated families (Pfam, Rfam), `hmmalign --trim --outformat afa profile.hmm new_seqs.fa > out.fasta` is the canonical alternative; it constrains insertions to lowercase columns rather than introducing new alignment columns and is the format expected by downstream HMMER/HHsuite tooling.

## Running MUSCLE5

Pick `-align` (PPP) for peak-accuracy runs up to ~1000 sequences, or `-super5` (mBed clustering + chunked alignment) for thousands to millions. `-super5` is not a "lower quality" mode; both share the same HMM-perturbation ensemble machinery.

| Command | Algorithm | Designed for | Output |
|---------|-----------|--------------|--------|
| `-align` (PPP) | Posterior probability progressive (HMM) | <= ~1000 seqs, peak accuracy | Single MSA or .efa ensemble |
| `-super5` | mBed-clustering + chunked alignment | Thousands to millions of seqs | Single MSA or .efa ensemble |

### Basic Usage

```bash
# Peak accuracy PPP algorithm, <1000 sequences
muscle -align input.fasta -output aligned.fasta -threads 8

# super5 divide-and-conquer for thousands of sequences
muscle -super5 input.fasta -output aligned.fasta -threads 8
```

### Ensemble Mode for Alignment Confidence

**Goal:** Quantify alignment uncertainty by generating multiple HMM-perturbed alignments and measuring column consistency.

**Approach:** MUSCLE5 (Edgar 2022 Nat Comm) ships two ensemble modes: `-stratified` (4 replicates by default) and `-diversified` (100 replicates by default). Both write an Ensemble FASTA (.efa) file containing all replicates; column-level confidence is the fraction of replicates that place a given residue pair in the same column. `-perturb SEED` is a separate flag that sets the HMM-perturbation random seed, not an ensemble selector.

```bash
# Stratified ensemble: 4 replicates spanning HMM perturbation strata
muscle -super5 input.fasta -stratified -output ensemble.efa

# Diversified ensemble: 100 replicates exploring guide-tree and HMM space
muscle -super5 input.fasta -diversified -output ensemble.efa

# Optional: change replicate count and HMM-perturbation seed
muscle -super5 input.fasta -stratified -replicates 8 -perturb 42 -output ensemble.efa
```

The .efa output is consumed downstream to derive confidence-weighted bootstrap support: each replicate is fed to a tree builder and the resulting trees combined (Edgar 2022 supplement). Columns consistently aligned across replicates are reliable; high-divergence regions diverge between replicates and should be flagged before phylogenetic inference.

## Running ClustalOmega

Use ClustalOmega when datasets reach hundreds of thousands of sequences (the mBed guide tree scales O(N log N)) or for HMM-profile-driven alignment. Prefer MAFFT or MUSCLE5 below that scale.

```bash
# Basic alignment
clustalo -i input.fasta -o aligned.fasta --auto

# Force overwrite output
clustalo -i input.fasta -o aligned.fasta --force

# Specify output format
clustalo -i input.fasta -o aligned.phy --outfmt=phylip

# Use more iterations for better accuracy
clustalo -i input.fasta -o aligned.fasta --iter=5

# Profile-profile alignment (align two existing MSAs)
clustalo --p1 profile1.fasta --p2 profile2.fasta -o merged.fasta

# Add sequences to existing alignment
clustalo -i new_seqs.fasta --profile1 existing.fasta -o updated.fasta

# Multithreaded
clustalo -i input.fasta -o aligned.fasta --threads=8
```

## Running T-Coffee

Use T-Coffee for small datasets (<50 sequences) where maximum accuracy matters or where structural templates exist (Expresso, 3D-Coffee). It is slower than progressive aligners but integrates diverse evidence via consistency-based library scoring.

| Mode | Flag | What it does |
|------|------|-------------|
| Default | (none) | T-Coffee + Lalign pairwise library |
| M-Coffee | `-mode mcoffee` | Combines libraries from MAFFT, MUSCLE, ClustalW, ProbCons, T-Coffee, etc. |
| Expresso | `-mode expresso` | PSI-BLAST searches PDB for structural templates, runs SAP structural alignment (requires internet for PSI-BLAST and PDB lookups; offline alternative: 3D-Coffee with user-supplied templates) |
| 3D-Coffee | `-mode 3dcoffee -template_file templates.txt` | User-supplied PDB templates; SAP / TM-align pairwise structural library |
| R-Coffee | `-mode rcoffee` | RNA: combines sequence alignment with consensus secondary structure (RNAplfold) |
| Pro-Coffee | `-mode procoffee` | Promoter regions: enforces position-specific TF binding-site alignment |
| Reliability | `-evaluate -output score_ascii` | TCS column reliability score for an existing alignment |

```bash
t_coffee input.fasta -output fasta_aln -outfile aligned.fasta

t_coffee input.fasta -mode mcoffee -output fasta_aln -outfile aligned.fasta

t_coffee input.fasta -mode expresso -output fasta_aln -outfile aligned.fasta

t_coffee -infile aligned.fasta -evaluate -output score_ascii > tcs_scores.ascii
```

**When to use T-Coffee**: Small datasets (<50 sequences) where maximum accuracy matters, especially when structural information (PDB templates) is available. Expresso (Armougom et al 2006 NAR) and 3D-Coffee modes (Poirot et al 2004; O'Sullivan et al 2004 JMB) substantially improve correct-column rate over sequence-only T-Coffee when structures exist; verify the latest benchmark numbers in the project documentation. Expresso requires internet access for PSI-BLAST + PDB lookups. The TCS reliability score (Chang et al 2014 MBE) flags individual columns as reliable/unreliable for downstream filtering before phylogenetics.

## Codon-Aware Alignment

### When Codon Alignment Is Required

Coding sequences destined for selection analysis (dN/dS with PAML/codeml, HyPhy BUSTED/MEME/aBSREL) **must** be aligned respecting codon boundaries. Standard nucleotide MSA tools do not preserve reading frames and produce systematically incorrect dN/dS estimates -- Fletcher & Yang (2010 MBE) showed conventional aligners cause false-positive selection signals under PAML M8 even on clean simulated data.

### Codon-Alignment Tool Decision Tree

| Input cleanliness | Recommended tool | Notes |
|------------------|------------------|-------|
| Clean orthologs (no frameshifts, no internal stops) | MAFFT-protein + PAL2NAL (Suyama, Torrents & Bork 2006 NAR) | Fastest; standard PAML pipeline input |
| Recently duplicated paralogs (indel-rich) | PRANK +F codon (Loytynoja & Goldman 2008 Science) | Phylogeny-aware indel model; fewest false-positive selection calls (Fletcher & Yang 2010) |
| Frameshifts, pseudogenes, error-prone assemblies | MACSE v2 `alignSequences -fs <cost>` (Ranwez et al 2018 MBE) | Frameshift-tolerant; preserves reading frame across sequencing errors |
| Mixed dataset with some bad genes | OMM_MACSE pipeline (Scornavacca, Belkhir, Lopez et al 2019; Ranwez group) | HMMER homology check + MAFFT prealignment + MACSE refinement |
| HyPhy-grade quality (BUSTED, MEME, aBSREL input) | HyPhy `pre-msa.bf` / `post-msa.bf` (Pond lab; Kosakovsky Pond et al) | Strips stop codons, runs MSA at protein level, threads back, validates frames |

### PAL2NAL: Protein-Guided Codon Alignment

**Goal:** Thread a nucleotide coding sequence alignment onto a protein alignment to preserve reading frame.

**Approach:** Align protein sequences first (higher sensitivity), then use PAL2NAL to map the protein alignment back to codons. Suyama et al 2006 NAR; standard codeml input pipeline.

```bash
mafft --localpair --maxiterate 1000 proteins.fasta > proteins_aligned.fasta
pal2nal.pl proteins_aligned.fasta codons.fasta -output fasta > codons_aligned.fasta
pal2nal.pl proteins_aligned.fasta codons.fasta -output paml > codons_aligned.phy
```

**Non-standard genetic codes:** PAL2NAL by default uses the standard code (NCBI table 1). For mitochondrial (vertebrate=2, yeast=3, invertebrate=5), ciliate macronuclear (=6, UAA/UAG=Gln), or other non-standard codes, the protein-to-codon mapping mismatches and PAL2NAL silently produces wrong codon assignments. Specify explicitly with `-codontable N`:

```bash
pal2nal.pl proteins_aligned.fasta codons.fasta -output paml -codontable 2 > codons_mt.phy
```

For datasets mixing genetic codes (e.g. nuclear + mitochondrial CDS in one tree), translate each lineage with its own code BEFORE protein alignment, never apply a single code globally. See NCBI Translation Tables for the full numbering.

### PRANK +F Codon Mode

**Goal:** Align coding sequences under a phylogeny-aware indel model that does not over-collapse insertions.

**Approach:** PRANK +F (Loytynoja & Goldman 2008 Science) treats indels as evolutionary events on a tree, distinguishing insertions from deletions. Slower than MAFFT but recommended for selection-analysis prep.

```bash
prank -d=codons.fasta -o=prank_aligned -codon -F
```

The `+F` flag enforces "fewer false insertions"; without it PRANK behaves more like a conventional aligner. For dN/dS analysis under PAML M2a/M8, PRANK +F is the most conservative input choice.

### MACSE v2: Frameshift-Tolerant Codon Alignment

**Goal:** Align coding sequences directly at the codon level while tolerating frameshifts and internal stops.

**Approach:** MACSE scores based on amino acid translation while operating on DNA. The v2 toolbox (Ranwez et al 2018 MBE) ships a sub-program selector via `-prog`; the eight most-used sub-programs are listed below (run `macse -help` for the full set in the installed release):

| Sub-program | Purpose |
|-------------|---------|
| `alignSequences` | Align coding sequences (`-fs <cost>` sets frameshift penalty) |
| `enrichAlignment` | Add new sequences to existing codon alignment |
| `refineAlignment` | Refine an alignment with MACSE scoring |
| `splitAlignment` | Split into sub-alignments (e.g. by guide tree clades) |
| `exportAlignment` | Convert between MACSE-internal and standard formats |
| `trimAlignment` | Position-aware trimming preserving codon structure |
| `trimNonHomologousFragments` | Detect and remove non-homologous regions per sequence |
| `reportGapsAA2NT` | Map protein-level gap removals back to nucleotide alignment |

```bash
java -jar macse_v2.jar -prog alignSequences -seq coding_seqs.fasta \
    -out_NT aligned_nt.fasta -out_AA aligned_aa.fasta -fs 30

java -jar macse_v2.jar -prog enrichAlignment -align existing.fasta \
    -seq new_seqs.fasta -out_NT updated.fasta
```

### OMM_MACSE: Recommended Pipeline Wrapper

OMM_MACSE (Ranwez group, used in OrthoMaM v10+) chains: HMMER homology screening to drop non-homologous sequences, MAFFT pre-alignment for guide-tree, then MACSE refinement. This is the recommended pipeline for genome-scale ortholog datasets where some genes will contain frameshifts.

```bash
OMM_MACSE_v12.02.sif --in_seq_file orthogroup.fasta --out_dir omm_out --out_file_prefix orthogroup --genetic_code_number 1
```

### HyPhy pre-msa.bf / post-msa.bf

For HyPhy-grade dN/dS analyses (BUSTED, MEME, aBSREL, RELAX), Pond lab's standard workflow is:

```bash
hyphy pre-msa.bf --input cds.fasta
mafft --auto cds.fasta_protein.fas > cds.fasta_protein.msa
hyphy post-msa.bf --protein-msa cds.fasta_protein.msa --nucleotide-sequences cds.fasta_nuc.fas --output cds.codon.msa
```

`pre-msa.bf` strips internal stop codons and translates; `post-msa.bf` validates that all sequences remain in-frame after threading. Without this validation step, a single mis-aligned codon can cascade into a spurious episodic-selection call.

### Confidence Assessment

For publication-grade phylogenetics or selection analysis, alignment uncertainty must be quantified per column and unreliable columns excluded BEFORE downstream inference -- not after.

| Method | Tool | Output |
|--------|------|--------|
| Guide-tree perturbation | GUIDANCE2 (Sela et al 2015 NAR) | Per-column and per-residue reliability score (0-1); default cutoff 0.93. No GUIDANCE3 has been released; deep-learning successors are exploratory only. |
| Library-consistency | T-Coffee TCS (Chang et al 2014 MBE) | Per-column reliability score via `-evaluate -output score_ascii` |
| HMM ensemble | MUSCLE5 `-stratified` / `-diversified` | EFA file; column confidence = fraction of replicates supporting it |
| Co-optimal alignments | HoT (Landan & Graur 2007 MBE) | Heads-or-tails alternate optimal alignment for each column |

Mask columns below the reliability threshold before phylogenetic inference. Trees built from filtered alignments are markedly more stable to method choice.

**GUIDANCE2 thresholds are tool-specific.** The 0.93 default cutoff is calibrated against MAFFT-LINSI in Sela et al 2015. Running GUIDANCE2 on top of MUSCLE5 or PRANK uses different perturbation distributions and produces scores on a slightly different scale. For non-MAFFT base aligners, calibrate the threshold empirically using a benchmark with known-correct columns or use the tool-native ensemble scoring (MUSCLE5 `-stratified`) instead.

## Post-Alignment Validation Checklist

Before proceeding to downstream analysis, verify alignment quality:

1. **Visual inspection**: Scan for columns of mostly gaps with scattered residues (hallmark of misalignment)
2. **Gap distribution**: High gap fraction (>50% of columns with gaps) suggests problematic regions or inclusion of non-homologous sequences
3. **Sequence identity**: If average pairwise identity is <25% for proteins, alignment reliability is questionable
4. **Outlier sequences**: Sequences with excessive gaps relative to others may be non-homologous or fragments; consider removing and re-aligning
5. **Conservation pattern**: Functional domains should show clear conservation; absence of expected conserved motifs suggests alignment error or non-homology
6. **Run GUIDANCE2 or MUSCLE5 ensemble**: Quantify alignment confidence per column before phylogenetic inference

## When NOT to Run MSA

- **Non-homologous sequences**: MSA tools always produce an alignment, even for unrelated sequences; verify homology first (e.g., BLAST E-value < 1e-5)
- **Sequences below the twilight zone**: Below ~20% protein identity, sequence signal is lost in noise; structural alignment is needed
- **Different domain architectures**: Globally aligning multi-domain proteins with different domain orders produces meaningless results; align individual domains separately
- **Very different lengths without shared homology**: Aligning a 50-residue fragment against 1000-residue proteins globally forces biologically meaningless gaps; use local alignment or fragment-aware modes (E-INS-i)
- **Highly repetitive sequences**: Tandem repeats cause alignment ambiguity; specialized tools (e.g., TRUST for tandem repeats) may be needed

## Quick Reference

| Task | Command |
|------|---------|
| Best accuracy (<200 seqs) | `mafft --localpair --maxiterate 1000 in.fa > out.fa` |
| Large dataset | `mafft --retree 2 in.fa > out.fa` or `clustalo -i in.fa -o out.fa` |
| Uncertainty estimation | `muscle -super5 in.fa -stratified -output out.efa` |
| Codon-aware | Align protein first, then `pal2nal.pl prot.fa cds.fa -output fasta` |
| Add to existing MSA | `mafft --add new.fa existing.fa > updated.fa` |
| Profile merge | `clustalo --p1 msa1.fa --p2 msa2.fa -o merged.fa` |
| With structure info | `t_coffee in.fa -mode expresso` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| MAFFT runs out of memory | L-INS-i on too many sequences | Switch to FFT-NS-2 or auto mode |
| Alignment has many all-gap columns | Non-homologous sequences included | Filter input by BLAST hits first |
| ClustalOmega crashes on large input | Memory limits | Increase `--MAC` RAM or use MAFFT |
| PAL2NAL "length mismatch" | Protein/DNA sequences not from same genes | Verify correspondence with sequence IDs |
| MUSCLE5 slow on large datasets | Default algorithm not designed for >10K seqs | Use `-super5` mode |
| Poor alignment despite high identity | Wrong sequence type detection | Explicitly specify `--amino` or `--nuc` |

## Related Skills

- alignment/pairwise-alignment - Compare two sequences using PairwiseAligner
- alignment/alignment-io - Read/write MSA files in various formats
- alignment/msa-parsing - Parse, filter, trim, and assess MSA quality
- alignment/msa-statistics - Calculate identity, conservation, entropy metrics
- alignment/structural-alignment - Foldseek, TM-align, Foldmason for the dark proteome
- alignment/alignment-trimming - ClipKIT, trimAl, BMGE post-MSA column filtering
- phylogenetics/modern-tree-inference - Build phylogenetic trees from MSAs
- sequence-io/read-sequences - Read input sequences for alignment

## References

- Edgar RC. 2022. Muscle5: high-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. Nat Comm 13:6968.
- Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7. MBE 30:772-780.
- Sievers F et al. 2011. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol 7:539.
- Notredame C, Higgins DG, Heringa J. 2000. T-Coffee: a novel method for fast and accurate multiple sequence alignment. JMB 302:205-217.
- Loytynoja A, Goldman N. 2008. Phylogeny-aware gap placement prevents errors in sequence alignment and evolutionary analysis. Science 320:1632-1635.
- Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. 2018. MACSE v2: toolkit for the alignment of coding sequences accounting for frameshifts and stop codons. MBE 35:2582-2584.
- Suyama M, Torrents D, Bork P. 2006. PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. NAR 34:W609-W612.
- Sela I, Ashkenazy H, Katoh K, Pupko T. 2015. GUIDANCE2: accurate detection of unreliable alignment regions accounting for the uncertainty of multiple parameters. NAR 43:W7-W14.
- Chang JM, Di Tommaso P, Notredame C. 2014. TCS: a new multiple sequence alignment reliability measure to estimate alignment accuracy and improve phylogenetic tree reconstruction. MBE 31:1625-1637.
- Fletcher W, Yang Z. 2010. The effect of insertions, deletions, and alignment errors on the branch-site test of positive selection. MBE 27:2257-2267.
- Armougom F, Moretti S, Poirot O, Audic S, Dumas P, Schaeli B, Keduas V, Notredame C. 2006. Expresso: automatic incorporation of structural information in multiple sequence alignments using 3D-Coffee. NAR 34:W604-W608.
- McWhite CD, Armour-Garb I, Singh M. 2023. Leveraging protein language models for accurate multiple sequence alignments. Genome Res 33:1145-1153.
- Redelings BD. 2021. BAli-Phy version 3: model-based co-estimation of alignment and phylogeny. Bioinf 37:3032-3034.
