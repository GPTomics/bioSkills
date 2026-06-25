---
name: bio-rna-structure-covariation-analysis
description: Tests whether a proposed or predicted RNA secondary structure is supported by evolutionary covariation using R-scape, which scores compensatory substitutions against a phylogeny-aware null and estimates the statistical power of the alignment. Use when validating a conserved-structure claim before trusting it (the test that found no support for HOTAIR/Xist/SRA lncRNA structures); separating real covariation from phylogenetic correlation; deciding whether an alignment even has the power to test structure; or building a covariation-supported consensus (CaCoFold) to seed a covariance model or folding.
tool_type: cli
primary_tool: R-scape
---

## Version Compatibility

Reference examples tested with: R-scape 2.0+, Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Covariation Analysis

**"Is my RNA's conserved secondary structure real?"** -> Measure whether the base pairs covary across an alignment more than phylogeny and base composition alone would produce, and whether the alignment has the power to detect such covariation.
- CLI: `R-scape -s alignment.sto` to test a given consensus structure
- CLI: `R-scape --cacofold alignment.sto` to build a covariation-supported structure de novo (CaCoFold)

## The governing principle: covariation is the gold standard, but a negative needs POWER

A compensatory (covarying) mutation is the strongest possible evidence for a base pair: if two columns change together across evolution so as to PRESERVE pairing (a G-C in one species becoming A-U in another at the same two positions), that directly evidences the pair, and a structure conserved by covariation across a deep alignment beats any single-sequence thermodynamic prediction. R-scape (Rivas, Clements & Eddy 2017) tests whether observed pairwise covariation EXCEEDS a phylogeny-aware null, separating real structural covariation from the apparent covariation that phylogenetic correlation and biased composition produce on their own.

The load-bearing nuance is that "no significant covariation" is NOT automatically "no structure" -- it can mean the alignment lacks the POWER to detect covariation (too few sequences, or sequences too similar, so there is not enough variation to observe compensatory changes). R-scape estimates, for each pair, the probability it would be called significant if it were a true pair (its power), so the result is a THREE-way verdict, not pass/fail:

| Verdict | Covariation | Power | Meaning |
|---------|-------------|-------|---------|
| Supports a conserved structure | significant pairs found | -- | the structure has evolutionary evidence |
| Rejects a conserved structure | none significant | adequate power | enough variation to detect covariation, yet none -> structure not supported (HOTAIR/Xist/SRA) |
| Cannot infer | none significant | low power | too few/too-similar sequences -> the alignment cannot test structure; gather more diverse homologs |

Reporting only "R-scape found 0 significant pairs" without the power context is the central misuse: a low-power negative says nothing about the structure. R-scape draws the low- vs high-power line at an explicitly arbitrary 10% mean alignment power (the sum of per-pair power over the number of base pairs; Rivas et al. 2020): below ~10%, treat a negative as "cannot infer."

## How R-scape decides

R-scape computes a per-pair covariation statistic (the G-test by default, with average-product correction to remove background phylogenetic signal), builds a null distribution by simulating alignments under the inferred phylogeny and base composition, and assigns each pair an E-value. A pair is significantly covarying when its E-value is at or below the target (default 0.05). It reports the number of expected covarying pairs found, their positions, the inferred substitutions at each, and the per-pair power. Significance is judged against the phylogenetic null, so a raw "positive covariation score" is not enough -- only covariation ABOVE the null counts.

## Test a given consensus structure

The input is a Stockholm alignment with a `#=GC SS_cons` line (the structure to test) -- e.g. an Rfam SEED, an RNAalifold consensus, or a hand-curated structure.

```bash
# -s evaluates the pairs in the alignment's SS_cons; -E sets the E-value target (default 0.05).
# --outdir keeps R-scape's outputs (.cov, .power, .sorted.cov, R2R .svg) out of the CWD.
R-scape -s -E 0.05 --outdir rscape_out alignment.sto
```

With `-s`, R-scape runs TWO tests: one on the pairs in the proposed SS_cons, and a separate one on all OTHER possible pairs -- so a significantly covarying pair OUTSIDE the proposed structure is evidence for a better or alternative fold, not just a yes/no on the given one. (A bare `R-scape alignment.sto` without `-s` tests all possible pairs as one set; `-s` is what scopes the primary test to the proposed structure.)

Outputs include `<msa>.cov` (covarying pairs: positions, score, E-value, substitutions, power), `<msa>.power` (power analysis), `<msa>.sorted.cov`, and an R2R `.svg`/`.pdf` diagram. Read the diagram by its legend: R-scape marks significantly covarying pairs distinctly from pairs that are merely structurally compatible and from pairs inconsistent with the covariation, so the highlighted pairs are the ones with evolutionary support. The header reports nseq, alignment length, average identity, and number of base pairs.

## Build a covariation-supported structure de novo (CaCoFold)

When there is no trusted structure to test, let covariation drive the fold. CaCoFold (`--cacofold`, also accepted as `--fold`) maximizes the support from significantly covarying pairs and can include pseudoknots as additional structure layers.

```bash
# Predict a structure from the alignment's covariation; writes a CaCoFold .sto with a new SS_cons.
R-scape --cacofold -E 0.05 --outdir rscape_out alignment.sto
```

The CaCoFold structure is grounded in evolutionary evidence rather than thermodynamics alone, which makes it a strong consensus to seed a covariance model (ncrna-search) or to compare against a thermodynamic fold (secondary-structure-prediction).

## The lncRNA cautionary tale

R-scape found NO statistically significant covariation support for the proposed secondary structures of the lncRNAs HOTAIR, SRA, and Xist (Rivas et al. 2017), despite their being thermodynamically plausible and widely cited. The lesson: a thermodynamically reasonable, even phylogenetically suggestive, structure is NOT established until covariation is statistically demonstrated, and for many lncRNAs the structural conservation is simply not there. Always run this test before asserting a conserved structure, and always report whether a negative is a power-limited "cannot infer" or a powered "rejects."

## Practical requirements

- Use a DEEP, DIVERSE alignment: covariation needs sequences that actually vary at paired columns while preserving the pair. Near-identical sequences carry no covariation signal (low power); a handful of sequences cannot test structure.
- Power is driven by the number of independent SUBSTITUTIONS at a pair, not the raw sequence count: a few near-identical sequences have essentially zero power, and meaningful power usually needs dozens of homologs spanning a broad identity range (well below ~90-95% average pairwise identity, ideally down toward ~60%). Size the alignment by R-scape's `.power` output, not by a fixed sequence count.
- Building the alignment is the real bottleneck for poorly conserved RNAs (lncRNAs especially): gather diverged homologs (e.g. synteny-anchored orthologs, an Infernal/cmsearch sweep, or RNAcentral) and align with a structure-aware aligner -- a bad alignment both destroys real covariation and manufactures spurious signal.
- Alignment quality matters: misaligned columns destroy real covariation and can manufacture spurious signal. Validate the alignment before trusting either a positive or a negative.
- Covariation tests CONSERVATION of structure, which is a different question from whether the RNA is a real, expressed transcript -- "is it real" in the transcription/processing sense needs expression and functional evidence, not R-scape.
- R-scape tests pairs, not whole helices alone; helix-level aggregation is available in recent versions for noisier alignments -- check `R-scape --help` for the current options.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "0 significant pairs" reported as "no structure" | ignoring power | check the `.power` output; a low-power negative is "cannot infer", not "rejects" |
| R-scape exits with no pairs tested | alignment has no `#=GC SS_cons` and `-s` was used | add a consensus structure, or use `--cacofold` to predict one |
| Outputs (.cov, .svg) dumped into the working directory | no output directory set | pass `--outdir <dir>` |
| Spurious covariation across the whole alignment | misaligned columns or strong phylogenetic correlation | improve the alignment; R-scape's null already corrects phylogeny, but bad alignments still mislead |
| A positive covariation score assumed to validate a pair | score is not significance | require E-value <= target (0.05) against the phylogenetic null, not a raw positive score |

## Related Skills

- secondary-structure-prediction - Predict the structure whose conservation is then tested
- ncrna-search - Validate a custom CM's SS_cons here before building the covariance model
- structure-probing - Experimental evidence complementary to evolutionary covariation
- alignment/msa-statistics - Assess the alignment depth and diversity covariation needs
- phylogenetics/tree-io - The phylogeny underlying the covariation null

## References

- Rivas E, Clements J, Eddy SR. 2017. A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs. Nat Methods 14(1):45-48. doi:10.1038/nmeth.4066
- Rivas E, Clements J, Eddy SR. 2020. Estimating the power of sequence covariation for detecting conserved RNA structure. Bioinformatics 36(10):3072-3076. doi:10.1093/bioinformatics/btaa080
- Rivas E. 2020. RNA structure prediction using positive and negative evolutionary information. PLoS Comput Biol 16(10):e1008387. doi:10.1371/journal.pcbi.1008387
- Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29(22):2933-2935. doi:10.1093/bioinformatics/btt509
