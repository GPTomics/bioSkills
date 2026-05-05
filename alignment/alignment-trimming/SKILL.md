---
name: bio-alignment-trimming
description: Trim multiple sequence alignments using ClipKIT, trimAl, BMGE, Divvier, or HMMcleaner with mode selection guidance per downstream goal. Use when removing unreliable columns or contaminating residues before phylogenetic inference, HMM building, or selection analysis.
tool_type: mixed
primary_tool: ClipKIT
---

## Version Compatibility

Reference examples tested with: ClipKIT 2.1+, trimAl 1.4+, BMGE 1.12+, Divvier 1.01+, HMMcleaner (current CPAN release of `Bio::MUST::Apps::HmmCleaner`), BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `clipkit --version`, `trimal --version`, `BMGE --help`, `Divvier --help`
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Alignment Trimming

**"Remove unreliable columns from this MSA"** → Filter or split columns based on gap fraction, conservation, entropy, or per-residue quality.
- CLI: `clipkit`, `trimal`, `BMGE`, `Divvier`, `HMMcleaner`
- Python: post-process via Bio.AlignIO with custom column masks

**"Make this alignment publication-grade for phylogenetics"** → Apply ClipKIT's `kpic-smart-gap` mode, or trimAl `-automated1`, then verify via tree-stability comparison before vs after trimming.

Tool choice and aggressiveness matter more than trimming vs not-trimming. Pick a mode by dataset character (table below), and always run a sensitivity check by building trees on trimmed and untrimmed alignments.

### Pick a Trimming Mode by Dataset Character

| Dataset character | Trimming effect | Recommended approach |
|-------------------|-----------------|----------------------|
| Deep-divergence orthologs (>500 Ma), saturated 3rd codons | Aggressive trimming HURTS (Tan-style result) | No trim or `kpic-smart-gap` only; report sensitivity to trimming choice |
| Mid-depth eukaryotic (animal phyla, fungal classes) | ClipKIT `kpic-smart-gap` HELPS | Steenwyk-style result |
| Shallow (within-genus) | All trimmers ~equivalent | Choose for downstream-tool compatibility |
| Concatenated supermatrix with very long alignments (>10 kb) | Trimming reduces phylogenetic noise | `kpic-smart-gap` or BMGE `-h 0.5` |
| Single short genes (<200 bp aligned) | Trimming amplifies stochastic error | Skip column trimming; use sequence-level outlier filtering |

**The 20%/40% rule.** Tan et al 2015 (Syst Biol) and Steenwyk et al 2020 (PLOS Bio) appear to disagree but their conclusions are reconciled by trimming aggressiveness: light trimming (<20% of columns removed) has minimal impact on tree accuracy regardless of method; heavy trimming (>40%) removes phylogenetic signal alongside noise and degrades tree accuracy on most empirical datasets. Steenwyk's ClipKIT `kpic-smart-gap` improves trees because it stays in the light-trim regime; the older Gblocks defaults fail because they over-trim. **Operational rule:** if your trimmer removes >40% of columns, the mode is too aggressive for the dataset; switch to a less aggressive mode or skip trimming.

Always run a sensitivity analysis: build the tree on trimmed AND untrimmed alignments. If topology and support are stable across trimming choices, the conclusion is robust; if unstable, report this and pick the result better supported by independent evidence (gene-tree concordance, biological priors).

## Goal-Driven Tool Selection

| Downstream goal | First-line tool | Rationale |
|----------------|-----------------|-----------|
| Phylogenetic-tree input (concatenated genes) | ClipKIT `kpic-smart-gap` (Steenwyk et al 2020 PLOS Bio) | Retains parsimony-informative + constant sites; consistently produces better trees |
| Phylogenetic-tree input (single gene) | ClipKIT `smart-gap` | Default mode; dynamic threshold determination |
| HMM profile building (HMMER, HHsuite) | trimAl `-gappyout` (Capella-Gutierrez et al 2009 Bioinf) | Aggressive gap removal acceptable; profile quality benefits |
| Selection / dN/dS analysis (PAML, HyPhy) | TCS column masking or GUIDANCE2 (NOT aggressive trimming) | Removing columns causes false-positive selection signals (Fletcher & Yang 2010 MBE) |
| Deep prokaryotic phylogenomics | BMGE (Criscuolo & Gribaldo 2010 BMC Evol Biol) | Entropy-based with BLOSUM62 context; standard in GToTree pipeline |
| Cross-contaminated sequences | HMMcleaner (Di Franco et al 2019 BMC Evol Biol) | Per-residue cleaning; targets contamination not column quality |
| Preserve phylogenetic signal in indels | Divvier (Ali, Bogusz & Whelan 2019 MBE) | Splits ambiguous columns rather than removing them |
| Column-mapping retention for site analysis | trimAl `-colnumbering` | Outputs original-column indices for downstream cross-reference |
| Codon-aware trimming | MACSE `trimAlignment` (Ranwez et al 2018 MBE) | Preserves codon boundaries; pairs with MACSE codon MSA |

## ClipKIT Modes

Pick a `-m` mode by what the downstream task needs the alignment to KEEP. Run `clipkit --help` for the full list of 15 modes.

| Pick this mode | When |
|----------------|------|
| `kpic-smart-gap` | Publication-grade phylogenomics on concatenated supermatrices (recommended default) |
| `smart-gap` | Single-gene tree input or unbalanced datasets where `kpic` over-trims outgroups |
| `gappy` / `gappyout` | Need an explicit gap threshold or trimAl-equivalent behaviour |
| `kpi-smart-gap` | Maximum-parsimony tree input (drops constant sites) |
| `c3` / `cst` / `entropy` / `block-gappy` / `composition-bias` / `heterotachy` | Specialist needs: third-codon stripping, manual masks, entropy thresholding, conserved-block preservation, long-branch mitigation, or heterotachy filtering -- consult `clipkit --help` |

**Goal:** Run ClipKIT on an MSA with the mode appropriate to the downstream goal.

**Approach:** Pick a `-m` mode from the table above, optionally enable `--log` for reproducibility, and write the trimmed alignment in the desired format.

```bash
clipkit input.fasta -m kpic-smart-gap -o trimmed.fasta
clipkit input.fasta -m gappy -g 0.9 -o trimmed.fasta
clipkit input.fasta -m kpic-smart-gap --output-format phylip -o trimmed.phy
clipkit input.fasta -m kpic-smart-gap --log -o trimmed.fasta
```

The `--log` flag produces a per-column kept/removed log -- essential for reproducibility audits. ClipKIT also has a Python API (`clipkit.api.clipkit`) for in-memory trimming.

### ClipKIT kpic Failure Modes

`kpic` defines parsimony-informative as "at least 2 sequences with each of at least 2 different residues". On taxonomically unbalanced datasets (e.g. 95 close relatives + 5 outgroups), columns informative WITHIN the dominant clade look informative GLOBALLY, while columns that distinguish the outgroup from the dominant clade may have a single residue in the outgroup and fail kpic. Result: the trimmed alignment is biased toward the dominant clade's signal. Mitigations:
- Use `kpic-smart-gap` AND verify outgroup branch length is preserved before vs after trimming
- For unbalanced datasets, consider `smart-gap` (no kpic constraint) or sequence-weight-aware trimming (no production tool exists; manual subsampling to balance clades is the practical fix)
- ClipKIT does NOT use Henikoff-style sequence weighting; the count-based informativeness is unweighted (Steenwyk et al 2020 supplement)

For bias-aware trimming, BMGE's entropy thresholding with `-w` window weighting (default 3) is partially weighted but not by clade structure. ClipKIT GitHub issues #71 and #88 document the unbalanced-dataset failure.

## trimAl Modes

Pick a trimAl mode by downstream tool. `-gappyout` for HMM profile builds, `-strictplus` for phylogenetic tree input, `-automated1` only when audit-grade reproducibility is not required (heuristic has changed across point releases).

| Mode | Flag | Description |
|------|------|-------------|
| Automated heuristic | `-automated1` | Picks `gappyout`, `strict`, or `strictplus` based on alignment characteristics |
| Gap-only | `-gappyout` | Automatic gap-fraction threshold from gap distribution |
| Strict | `-strict` | Combined gap + similarity criterion; recommended for phylogenetics |
| Strict-plus | `-strictplus` | Stricter; ~50% more aggressive than `-strict` |
| Manual gap | `-gt 0.5` | Remove columns with > 50% gaps (set fraction explicitly) |
| Manual similarity | `-st 0.5` | Remove columns with similarity below threshold |
| Manual conservation | `-cons 60` | Keep at least 60% of original columns regardless of other criteria |
| Sequence-quality | `-resoverlap 0.8 -seqoverlap 75` | Remove sequences with poor residue/sequence overlap |
| HTML report | `-htmlout report.html` | Visual diff of kept/removed columns |
| Column mapping | `-colnumbering` | Output original column indices preserved |

**Goal:** Trim an MSA with trimAl using a heuristic or manual threshold suited to the downstream tool.

**Approach:** Select an automated mode (`-automated1`, `-gappyout`, `-strictplus`) for typical use, or compose `-gt`, `-st`, `-cons` thresholds for manual control; export an HTML diff or column index mapping when audit-grade reproducibility is required.

```bash
trimal -in input.fasta -out trimmed.fasta -automated1
trimal -in input.fasta -out trimmed.fasta -gappyout
trimal -in input.fasta -out trimmed.fasta -strictplus -htmlout report.html
trimal -in input.fasta -out trimmed.fasta -gt 0.3 -st 0.5 -cons 60 -colnumbering > columns.txt
```

`-gappyout` is the recommended choice for HMM profile building (HMMER `hmmbuild` benefits from aggressive gap removal). `-strictplus` is recommended for phylogenetic tree input. `-automated1` is the safe default when characterising the dataset is impractical.

**Reproducibility note:** `trimAl -automated1` selects between `gappyout`, `strict`, and `strictplus` via internal heuristics that have changed between point releases (1.4 -> 1.4.1 -> 2.0-rc). For audit-grade reproducibility, do NOT use `-automated1`; specify the underlying mode explicitly. Record `trimal --version` in pipeline manifests. trimAl 2.0 (in development) is expected to change defaults; until release, pin to 1.4.1 in production environments.

## BMGE: Block Mapping and Gathering with Entropy

Use BMGE for deep prokaryotic phylogenomics (the GToTree default) or whenever entropy-based, matrix-aware column filtering is preferred over gap-only heuristics.

**Goal:** Filter MSA columns by matrix-aware entropy and gap-rate thresholds, particularly for deep prokaryotic phylogenomics.

**Approach:** Invoke BMGE.jar with `-t` (sequence type) and tune `-h` (entropy) and `-g` (gap rate) against the divergence depth of the dataset; recalibrate `-h` if the substitution matrix is changed via `-m`.

```bash
java -jar BMGE.jar -i input.fasta -t AA -of trimmed.fasta -h 0.5 -g 0.2

java -jar BMGE.jar -i input.fasta -t DNA -of trimmed.fasta -m DNAPAM100:2 -h 0.5
```

| Flag | Default | Meaning |
|------|---------|---------|
| `-h` | 0.5 | Entropy threshold (lower = more aggressive trimming) |
| `-g` | 0.2 | Gap rate threshold |
| `-m` | BLOSUM62 (AA) | Substitution matrix for entropy weighting |
| `-t` | -- | Required: AA, CODON, or DNA |
| `-b` | 5 | Minimum block size |

`-h 0.4` is the recommended setting for deep phylogenomics where conservation is low; `-h 0.6` retains more sites for shallower phylogenies.

**Matrix-aware threshold:** BMGE computes entropy weighted by the chosen substitution-matrix probabilities. The `-h 0.5` and `-h 0.4` thresholds are calibrated for BLOSUM62 (default). Switching to `-m BLOSUM30` (deep phylogenomics) makes the same numeric `-h` more permissive; `-m BLOSUM90` (close orthologs) makes it more aggressive. When using a non-default matrix, recalibrate by running on a known-good benchmark before applying.

## Divvier: Column Splitting Instead of Removal

Use Divvier when phylogenetic signal in indels matters and pure removal would discard block-boundary information; it splits ambiguous columns rather than dropping them.

```bash
# Divvier ships as a compiled binary; invoke as ./divvier or `divvier` if on PATH
./divvier -divvy input.fasta

./divvier -partial -mincol 4 -divvygap input.fasta
```

`-divvy` (default) outputs `input.fasta.divvy.fas` (full divvying); `-partial` only filters individual ambiguous characters. `-mincol N` sets the minimum number of confident characters required to keep a split column (integer count, default 2). `-divvygap` writes gaps instead of asterisks for phylogenetics-tool compatibility. Divvier is particularly useful when input alignments have many short conserved blocks separated by ambiguous regions -- removing the regions discards information about block boundaries that splitting preserves.

Maintenance note: Divvier's last release was 2019 (`simonwhelan/Divvier`); the tool is no longer actively maintained but results remain reproducible with the v2019 binary.

## HMMcleaner: Per-Residue Cleaning

Run HMMcleaner before column trimming when cross-contamination or annotation errors are suspected; it masks suspect residues with `X` rather than removing columns. Apply only to alignments with >= 15 sequences (false-positive rate climbs to ~10% at 10 sequences).

```bash
HmmCleaner.pl input.fasta
HmmCleaner.pl input.fasta --no-large-remove
```

Output is `input_hmm.fasta` with low-confidence residues masked. Used in OrthoMaM and Compositional Heterogeneity-aware phylogenomic pipelines where cross-contamination from sequencing or annotation errors is suspected.

**Sample-size sensitivity:** HMMcleaner builds a per-residue HMM from the OTHER sequences at each position, so on small alignments (<15 sequences) the HMM has insufficient signal and over-flags real divergent residues as contamination. Di Franco et al 2019 supplementary results show false-positive rate climbs to ~10% at 10 sequences and stabilises below 2% at >50 sequences. OrthoMaM v11 and PhyloHerb pipelines apply HMMcleaner only to alignments with >= 15 sequences. For smaller alignments, manual inspection or BLAST-based contamination screening is more reliable.

## Gblocks: The Legacy Default

Use Gblocks only when matching a legacy pipeline; prefer ClipKIT or trimAl for new work. Default parameters are too aggressive -- always relax with `-b1=50% -b2=85% -b3=10 -b4=5 -b5=h`.

```bash
Gblocks input.fasta -t=p -b1=50% -b2=85% -b3=10 -b4=5 -b5=h
```

## PhyIN: Phylogenetic Incompatibility Trimming

PhyIN (Maddison 2024 PeerJ) is a complementary trimmer that flags neighbouring columns whose split patterns are phylogenetically incompatible (i.e. cannot share a tree). It targets a different failure mode than gap-based or entropy-based trimmers: ClipKIT, trimAl, and BMGE all evaluate columns independently, so a region with consistent gap structure passes their filters even if its character patterns are pairwise tree-incompatible (alignment artefact). PhyIN catches that case.

```bash
phyin input.fasta -o trimmed.fasta -w 5 -t 0.5
```

Use PhyIN as a SECOND-PASS trimmer after ClipKIT/trimAl when alignment artefact is the suspected source of incongruence in a phylogenomic dataset; not a replacement for first-pass column filtering. PhyIN's incompatibility test is signal-direction agnostic, so it does not distinguish "the alignment is wrong here" from "true incongruent locus" (e.g. introgression, ILS); only gene-tree comparison after tree-building can disambiguate. Reference: PeerJ 2024 paper.

## Decision Tree by Downstream Goal

```
What is the next step?
+- Phylogenetic ML tree (RAxML, IQ-TREE)
|  +- Concatenated supermatrix? -> ClipKIT kpic-smart-gap
|  +- Single gene? -> ClipKIT smart-gap or trimAl -automated1
|  +- Deep prokaryotic? -> BMGE -h 0.4 -g 0.2
|  +- Suspected cross-contamination? -> HMMcleaner first, then ClipKIT
|
+- Bayesian inference (MrBayes, BEAST)
|  +- Same as ML; trimming is acceptable but log column mapping
|
+- HMM profile (HMMER hmmbuild, HHsuite)
|  +- trimAl -gappyout (aggressive gap removal helps profile quality)
|
+- Selection analysis (PAML codeml, HyPhy)
|  +- DO NOT aggressively trim
|  +- Use TCS column masking (T-Coffee -evaluate) or GUIDANCE2
|  +- See alignment/multiple-alignment confidence assessment
|
+- Sequence logo / motif scan
   +- Light gap-only trim (trimAl -gt 0.5)
   +- Preserve all variable positions
```

## TCS Column Masking for Selection Analysis

TCS (Chang et al 2014 MBE) scores each column on a 0-9 scale by consistency across a pairwise-alignment library. The recommended dN/dS workflow:

**Goal:** Mask unreliable columns before selection analysis to prevent alignment errors from inflating false-positive dN/dS calls.

**Approach:** Build a library-rich consistency score via M-Coffee, evaluate the alignment against that library to obtain per-column TCS scores, then keep columns at or above a chosen confidence threshold using `seq_reformat`.

```bash
# 1. Generate library-rich consistency scores via M-Coffee (combines libraries from MAFFT, MUSCLE, ClustalW, ProbCons)
t_coffee input.fasta -mode mcoffee -output fasta_aln -outfile aligned.fasta

# 2. Compute per-column TCS scores against the same library
t_coffee -infile aligned.fasta -evaluate -output score_ascii > tcs_scores.ascii

# 3. Filter columns at a chosen threshold via seq_reformat (5-9 = retain columns scoring 5 or higher)
t_coffee -other_pg seq_reformat -in aligned.fasta -struc_in tcs_scores.ascii \
    -struc_in_f number_aln -action +keep '[5-9]' -output fasta_aln > aligned_tcs5.fasta
```

Threshold guidance: TCS >= 5 retains columns confidently aligned by the majority of library methods; TCS >= 7 is conservative and used by Chang et al 2014 for HoT-validated benchmarks. TCS-from-mcoffee gives substantially better column-confidence ranking than TCS-from-default-tcoffee because the library is more diverse. For PAML branch-site test, Fletcher & Yang 2010 found TCS-5 masking eliminates false positives nearly completely on simulated data.

## MACSE Frameshift Markers Need Post-Processing

MACSE encodes detected frameshifts as `!` (within-codon insertion) and `*` (premature stop) in the nucleotide output. PAML codeml interprets `!` as missing data; HyPhy `BUSTED`/`MEME` interpret `!` as `N` but `*` triggers a parse error mid-sequence. Standard cleanup before downstream analysis:

**Goal:** Convert MACSE frameshift and stop markers into formats that PAML and HyPhy will parse correctly.

**Approach:** Replace `!` with gaps for PAML, and use MACSE's own `exportAlignment` sub-program with `-codonForInternalStop NNN` for HyPhy-safe output.

```bash
# Convert frameshift markers to gaps; replace internal stops with NNN
sed -e 's/!/-/g' aligned_nt.fasta > aligned_paml.fasta
java -jar macse_v2.jar -prog exportAlignment -align aligned_nt.fasta \
    -codonForFinalStop --- -codonForInternalStop NNN \
    -out_NT aligned_hyphy.fasta -out_AA aligned_hyphy_aa.fasta
```

`-codonForInternalStop NNN` replaces internal stops with `NNN` (HyPhy-safe). Without this step, the dN/dS run silently produces results that do not correspond to the alignment shown.

## Aggressiveness Cap

Compute the retained-column fraction and warn if it drops below 0.7. Trimming more than 20-30% of columns rapidly degrades ML tree accuracy on empirical data.

```python
def trimming_fraction(input_alignment, trimmed_alignment):
    return trimmed_alignment.get_alignment_length() / input_alignment.get_alignment_length()
```

If retention drops below 0.7, switch to a less aggressive mode or accept the original alignment unfiltered.

## Column Mapping for Reproducibility

When trimming for phylogenetic input, retain the column-index mapping so downstream site-specific analyses (selection per site, structure-mapped residues, etc.) can be back-traced:

```bash
trimal -in input.fasta -out trimmed.fasta -automated1 -colnumbering > kept_columns.txt

clipkit input.fasta -m kpic-smart-gap --log -o trimmed.fasta
```

These two flags emit different formats. `trimal -colnumbering` writes a comma-separated list of the original column indices that survived (e.g. `0, 1, 2, 4, ...`) on stdout. `clipkit --log` writes a `<output>.log` file with one row per original column reporting `position\tkeep_or_trim\tsite_classification\tgap_proportion`. To compare runs across the two tools, parse each format separately and reconcile to a common "kept-column index list" before applying the index map for site-specific cross-reference.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| ClipKIT empty output | All columns failed `kpic` filter | Switch to `smart-gap` mode (less aggressive); check input for non-amino-acid characters |
| trimAl "all sequences are identical" | Input has duplicates with no variation | Remove duplicate sequences first |
| BMGE Java OutOfMemoryError | Default JVM heap too small | `java -Xmx16g -jar BMGE.jar` |
| Divvier crashes on large input | Memory limits | Run on per-gene alignments rather than concatenated supermatrices |
| HMMcleaner reports "no contamination" | Input alignment too small (< 5 sequences) | Pool more sequences or skip per-residue cleaning |
| Trimming removes 80%+ of columns | Mode too aggressive for divergent input | Switch to `kpic-gappy` or relaxed-parameter Gblocks |

## Related Skills

- alignment/multiple-alignment - Generate the input MSA before trimming; confidence assessment (GUIDANCE2, TCS, MUSCLE5 ensemble) for selection-analysis input
- alignment/msa-parsing - Parse the trimmed alignment; sequence weighting; coordinate mapping
- alignment/msa-statistics - Compute conservation and gap statistics to inform mode selection
- alignment/structural-alignment - Trim structural MSAs the same way as sequence MSAs
- phylogenetics/modern-tree-inference - Build trees from trimmed alignments

## References

- Steenwyk JL, Buida TJ, Li Y, Shen XX, Rokas A. 2020. ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. PLOS Bio 18:e3001007.
- Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T. 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinf 25:1972-1973.
- Criscuolo A, Gribaldo S. 2010. BMGE: a new software for selection of phylogenetic informative regions from multiple sequence alignments. BMC Evol Biol 10:210.
- Maddison WP. 2024. PhyIN: trimming alignments by phylogenetic incompatibilities among neighbouring sites. PeerJ 12:e18504.
- Tan G, Muffato M, Ledergerber C, Herrero J, Goldman N, Gil M, Dessimoz C. 2015. Current methods for automated filtering of multiple sequence alignments frequently worsen single-gene phylogenetic inference. Syst Biol 64:778-791.
- Fletcher W, Yang Z. 2010. The effect of insertions, deletions, and alignment errors on the branch-site test of positive selection. MBE 27:2257-2267.
- Chang JM, Di Tommaso P, Notredame C. 2014. TCS: a new multiple sequence alignment reliability measure. MBE 31:1625-1637.
