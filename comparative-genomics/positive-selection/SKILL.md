---
name: bio-comparative-genomics-positive-selection
description: Detect positive (diversifying / episodic / pervasive) selection using codon dN/dS frameworks. Implements PAML codeml site models (M0/M1a/M2a/M7/M8/M8a), branch models, branch-site model A (Zhang 2005), and HyPhy methods (BUSTED, BUSTED-S, BUSTED-MH, BUSTED-PH, MEME, FEL, FUBAR, aBSREL, SLAC, RELAX, GARD, FUBAR-MH). Includes McDonald-Kreitman framework (asymptotic alpha, impMKT, polyDFE, DFE-alpha, GRAPES) for within-species + divergence inference, RERconverge for trait-correlated rate shifts, CSUBST for convergent substitution, and PhyloAcc for accelerated noncoding evolution. Use when testing adaptive evolution at codons, branches, or full gene; running GARD recombination pre-screen; controlling alignment-error and gBGC false positives; reconciling PAML vs HyPhy results; or performing genome-scale selection scans.
tool_type: mixed
primary_tool: PAML
---

## Version Compatibility

Reference examples tested with: PAML 4.10.7+, HyPhy 2.5.62+ (BUSTED-MH from Lucaci 2023 bioRxiv; FUBAR-MH from same), datamonkey.org 2024+ for web jobs, IQ-TREE 2.3.6+, MACSE V2.07+, PRANK 170427+, MAFFT 7.526+, PREQUAL 1.02+, HmmCleaner 0.243+, GARD (HyPhy bundled), RDP5 5.59+, ete4 4.1.0+, BioPython 1.84+, scipy 1.13+, polyDFE 2.0+, DFE-alpha 2.16+, GRAPES 1.1.1+, RERconverge 0.3.0+, CSUBST 1.6.0+, PhyloAcc 2.4.0+. Quest-for-Selection benchmark refreshed annually.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `codeml` (PAML; check by `codeml /dev/null` -- prints version banner), `hyphy --version`, `gard --help`
- Python: `pip show pyhyphy`; introspect ete4 API for tree-labeling
- R: `packageVersion('RERconverge')`; `?correlateWithBinaryPhenotype`

If code throws `branch-site test LRT non-positive`, `omega2 hit upper bound 999`, `MEME ML mixed gradient`, the most common cause is alignment error or saturated dS -- inspect alignment with TCS / Guidance2 and dS-vs-divergence-time. PAML 4.10 changed several control-file keywords from 4.9 (`getSE = 1` syntax tightened).

# Positive Selection Analysis

**"Is this gene / branch / site under positive selection?"** -> dN/dS (omega = nonsynonymous-to-synonymous substitution rate ratio) framework with explicit choice of WHICH question is being asked (gene-wide / branch-specific / site-specific / episodic) and WHICH null is being rejected. The "test failed because of selection" claim has more known confounders than any other comparative-genomics inference; **mandatory pre-screens are: recombination (GARD), alignment errors (PREQUAL or HmmCleaner), saturation (dS distribution), and gBGC (W->S substitution bias)**. Skipping any one inflates Type-I error to ~20-50% (Anisimova & Yang 2007 MBE 24:1219; Pond 2006 Mol Biol Evol 23:1891).

- CLI: `codeml` PAML site, branch, branch-site models
- CLI: `hyphy busted` `hyphy meme` `hyphy fel` `hyphy fubar` `hyphy absrel` `hyphy relax` `hyphy gard`
- Web: datamonkey.org for HyPhy jobs without local install
- R: `RERconverge::correlateWithBinaryPhenotype()` for trait-rate associations
- CLI: `csubst analyze` for convergent substitution
- R/CLI: `phyloacc` for noncoding accelerated evolution

## Algorithmic Taxonomy

| Method | Question | Null model | Strength | Fails when |
|--------|----------|------------|----------|------------|
| PAML codeml M0 (Yang 1997 CABIOS 13:555) | Gene-wide single-omega estimate | -- (point estimate) | Standard reference omega; baseline test | Site heterogeneity (use M3+) |
| codeml M1a vs M2a (Yang 2000 Genetics 155:431) | Any site under selection? | Nearly neutral, 2-category | Conservative; LRT df=2 | Low power for episodic selection |
| codeml M7 vs M8 | More-sensitive site test | Beta(0,1) | Higher power than M1a/M2a | Higher false-positive rate; relaxed-constraint mimics selection |
| codeml M8 vs M8a (Swanson 2003 Genetics 165:725) | Conservative site test (omega2 = 1 null) | Beta + omega2=1 | Cleanest LRT df=1; preferred site test | Lower power than M7 vs M8 |
| codeml branch-site mod A (Zhang 2005 MBE 22:2472) | Selection on pre-specified foreground branch | A1 (omega2=1 fixed) | Most powerful for episodic per-branch selection | Foreground specified post hoc -> Type-I inflation |
| codeml clade model (Bielawski & Yang 2004 J Mol Evol 59:121) | Different omega between named clades | M3 with shared categories | Tests for shifted selection regime | Requires clade pre-specification |
| codeml free-ratio | Per-branch omega estimates (exploratory) | M0 | Visualizes branch-wise variation | Unidentifiable for short branches; no formal LRT |
| HyPhy BUSTED (Murrell 2015 MBE 32:1365) | Any episodic selection on any branch site? | No omega+ class | Site + branch joint; foreground assignable | Sensitive to alignment errors |
| HyPhy BUSTED-S (Wisotsky 2020 MBE 37:2430) | BUSTED with synonymous-rate variation | -- | Corrects for SRV; reduces false positives | Slightly less power than BUSTED |
| HyPhy BUSTED-MH (Lucaci 2023 bioRxiv 2022.12.02) | BUSTED with multi-nucleotide substitutions | -- | Captures complex (multi-hit) substitutions; reduces false positives from MNMs | Newer; limited benchmarking |
| HyPhy BUSTED-PH | Two phenotypes; selection on one not other | -- | Tests phenotype-specific selection | Requires phenotype branch label |
| HyPhy MEME (Murrell 2012 PLoS Genet 8:e1002764) | Per-site episodic selection | FEL | Detects sites under episodic positive selection | Higher false-positive rate at p threshold |
| HyPhy FEL (Kosakovsky Pond 2005 MBE 22:1208) | Per-site pervasive selection | -- | Fast; counts substitutions per site | No episodic detection |
| HyPhy FUBAR (Murrell 2013 MBE 30:1196) | Bayesian per-site pervasive selection | -- | Scales to 1000s of sequences; posterior probability | No episodic detection |
| HyPhy SLAC | Counting-based fast estimator | -- | Very fast; rough estimate | Lower power; no statistical model |
| HyPhy aBSREL (Smith 2015 MBE 32:1342) | Branch-specific selection without pre-specification | -- | Adaptive per-branch omega categories; corrects multiple testing | Multiple-testing burden across many branches |
| HyPhy RELAX (Wertheim 2015 MBE 32:820) | Selection relaxation (k<1) or intensification (k>1) | -- | Detects RELAXED selection; cannot be done by other tests | Not designed for adaptive evolution per se |
| HyPhy GARD (Pond 2006 MBE 23:1891) | Recombination breakpoint detection | No recombination | MANDATORY pre-screen for any selection test | Computationally heavy; > 50 sequences slow |
| McDonald-Kreitman (McDonald & Kreitman 1991 Nature 351:652) | Adaptive substitution rate alpha from poly + div data | Neutral mutation accumulation | Per-gene alpha; population genetics native | Slightly deleterious bias (downward); fixed by asymptotic alpha |
| Asymptotic alpha (Messer & Petrov 2013 PNAS 110:8615) | MK with slightly deleterious correction | -- | Unbiased alpha; works at low MAF SFS | Requires SFS data |
| impMKT (Murga-Moreno 2022 G3 12:jkac206) | MK with conservative imputation | -- | Gene-level evidence; faster than alpha asymptotic | Less unbiased than asymptotic alpha |
| polyDFE (Tataru & Bataillon 2019 Bioinformatics 35:2868) | Full DFE + alpha jointly | -- | Quantifies the distribution of fitness effects | Computational cost; requires polymorphism data |
| DFE-alpha (Eyre-Walker & Keightley 2009 MBE 26:2097) | Faster DFE method | -- | Standard DFE inference; many simulated DFEs | Requires demographic correction |
| GRAPES (Galtier 2016 PLoS Genet 12:e1005774) | DFE on neutral + selected sites | -- | Joint demography + alpha; robust | Genome-scale dataset required |
| RERconverge (Kowalczyk 2019 Bioinformatics 35:4815; Redlich 2024 MBE 41:msae210) | Relative-rate shifts correlated with categorical phenotype | -- | Phylogenome-wide trait associations | Inherits all dN/dS confounders |
| CSUBST (Fukushima & Pollock 2023 Nat Eco Evo 7:155) | Convergent substitutions across independent lineages | -- | Combinatorial-substitution omegaC ratio; null-corrected | Requires multi-clade dataset |
| PhyloAcc (Hu 2019 MBE 36:1086; Thomas 2024) | Bayesian convergent accelerated noncoding rate | -- | For noncoding elements (CNEs); convergent rate shifts | CDS analyses prefer codon-based methods |
| phyloP (Pollard 2010 GR 20:110) | Per-site noncoding rate test | -- | Simple; widely used for noncoding | No convergence; site-by-site |
| PRANK + codeml pipeline | Codon-aware MSA + codeml | -- | Standard publication-grade workflow | Slow for large datasets |

Methodology evolves; verify the latest HyPhy / PAML manuals and the Álvarez-Carretero "Beginner's Guide" (Álvarez-Carretero et al 2023 MBE 40:msad041) before locking on a single method. The BUSTED-MH and FUBAR-MH (multi-hit) extensions specifically address known Type-I inflation from multi-nucleotide substitutions and are now recommended over basic BUSTED / FUBAR.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Single gene, mammalian (~60 Myr), pre-specified foreground branch | codeml branch-site mod A AND HyPhy aBSREL on foreground | Mutual validation; mod A LRT df=1 + aBSREL adaptive site classes |
| Single gene, deep eukaryote (~500+ Myr), no foreground hypothesis | GARD pre-screen -> BUSTED-MH gene-wide -> MEME for sites | Episodic-selection-only methods; saturation-aware (HyPhy under MG94 codon model) |
| Genome-wide scan, vertebrates | codeml M7 vs M8 OR HyPhy FUBAR-MH per gene; FDR-correct | Pervasive-selection sites; multi-hit correction critical at scale |
| Episodic selection scan | HyPhy MEME genome-wide (per gene); FDR-correct | Site-level episodic detection |
| Branch-specific selection on unspecified branches | HyPhy aBSREL | Adaptive per-branch test with built-in multiple-testing |
| Comparing selection regimes between two phenotypes | HyPhy BUSTED-PH or RELAX | Phenotype-specific or relaxation-detection |
| Recently diverged species (low divergence) | MK / asymptotic alpha (population genetics) | Codon dN/dS unreliable at low divergence; SFS-based instead |
| Within-species, dense polymorphism + divergence | polyDFE / GRAPES / asymptotic-MK | Full DFE + alpha jointly; preferred for adaptive-substitution rate |
| Coding selection genome-wide, with SFS available | grapes -m AUTO_ALL | Demography-aware alpha; standard population-genetics-aware adaptive-substitution scan |
| Noncoding accelerated evolution (CNEs / ECRs) | PhyloAcc, phyloP-acc | Codon-based unsuitable; PhyloAcc Bayesian convergence |
| Convergent substitutions across independent lineages | CSUBST | Combinatorial-substitution omegaC; null-corrected |
| Trait-correlated rate shifts genome-wide | RERconverge | Categorical / binary phenotype; correlates RERs across thousands of genes |
| Suspected positive selection but dS > 2 | Use protein-level method or reduce taxon sampling | Codon-based methods unreliable at saturation; protein-only ASR can still work |
| Recombination expected (immune genes, viral genomes) | GARD pre-screen mandatory | Recombination + tree-based selection -> false positives (Anisimova 2003) |
| Convergent codon substitution at specific sites | TDG09 (Tamuri 2009) or PCOC (Rey 2018) | Detects convergence at specific sites associated with traits |
| Drug-target evolution screen | aBSREL on candidate genes; cross-validate with MEME | Recent positive selection at drug-target loci |
| Pathogen / immune-evasion gene with high dS variation | BUSTED-S (synonymous rate variation aware) | dS variation across sites violates basic BUSTED assumptions |
| Plasmodium / Trypanosoma / Plasmid analysis | BUSTED-MH (multi-hit aware) | Multi-nucleotide substitutions common in these; basic BUSTED inflates false positives |

## Per-Method Failure Modes

### Recombination producing false positive selection

**Trigger:** Running codeml or BUSTED on a gene with recombination breakpoints (viral genes, immune genes, paralog families).

**Mechanism:** All single-tree codon models assume one phylogeny across all sites. Recombination produces different trees for different segments; treating them as one tree forces the model to invent rate variation that mimics positive selection (Anisimova et al 2003 Genetics 164:1229).

**Symptom:** PAML M8 strongly rejects M7 (LRT > 50), with omega2 = 999 (PAML upper bound) at several "selected sites"; HyPhy BUSTED highly significant; sites clustered in specific gene regions.

**Fix:** **MANDATORY: run GARD before any positive selection test.** If GARD detects breakpoints (p < 0.05), partition the alignment at breakpoints and analyze each segment separately, or use the recombination-aware MEME with the partitioned tree set. RDP5 (Martin 2021 Virus Evol 7:veaa087) is an alternative for viral genomes. GARD output `.json` lists breakpoint positions and posterior support.

### Alignment errors producing false positives

**Trigger:** Using default MAFFT or MUSCLE alignment on divergent CDS sequences; skipping codon-aware aligner.

**Mechanism:** Frame-shifted or misaligned codons introduce apparent non-synonymous substitutions at every position; codon-aware tools see these as positive selection (Schneider 2009 GBE 1:114; Markova-Raina & Petrov 2011 GR 21:863).

**Symptom:** "Selected sites" cluster in alignment regions with > 30% gaps; per-site posteriors in BEB / FUBAR concentrate in ambiguous columns; PREQUAL or Guidance2 marks these regions as poorly aligned; protein alignment shows obvious mismatches.

**Fix:** Use codon-aware aligner: **PRANK** (Loytynoja 2014 Methods Mol Biol 1079:155) is the standard for selection analysis (correctly models insertions); MACSE V2 (Ranwez 2018 MBE 35:2582) handles frameshifts and pseudogenes natively; OMM_MACSE wrapper combines them. After alignment, filter with PREQUAL (segment-level) or HmmCleaner (Di Franco 2019 BMC Eco Evo 19:21); do NOT use block-filtering (Gblocks, trimAl) which removes informative sites. Segment-level filtering preferred for selection (Di Franco 2019).

### Saturated synonymous sites

**Trigger:** Comparing distantly related taxa (deep eukaryotic divergence, > 100 Myr); dS > 3 across most pairs.

**Mechanism:** Synonymous sites have undergone multiple substitutions; the observed dS underestimates true dS. The model can't recover the true rate; omega = dN/dS becomes unstable at the upper bound or low (depending on which direction the bias goes).

**Symptom:** PAML M0 omega = 999 or near-zero; per-branch dS variance huge; sites with omega > 1 in M8 BEB are at conserved residues (paradox).

**Fix:** Reduce taxon sampling to species with dS < 2 on internal branches. For deep selection inference on conserved residues, use protein-level methods (BUSTED with `--model GTR` AA codon translation; aBSREL with protein model option) or restrict to subclade with reasonable saturation. Yang 2007 PAML manual recommends dS < 1.5 per branch.

### gBGC inflating apparent positive selection

**Trigger:** Mammalian / vertebrate gene with W->S substitution bias on a fast-evolving lineage.

**Mechanism:** GC-biased gene conversion fixes A/T -> G/C alleles preferentially in regions of high recombination, independent of selection (Galtier & Duret 2007 Trends Genet 23:273; Capra 2013 PLoS Genet 9:e1003684). Standard codon models attribute this to positive selection because nonsynonymous substitutions are unequally distributed across codon positions.

**Symptom:** Branch with apparent positive selection sits in high-recombination region; W->S / S->W substitution ratio > 1.5; selected sites concentrate at non-degenerate codon positions; HyPhy MEME-MH and BUSTED-MH attribute signal to multi-hit rather than positive selection.

**Fix:** Test for gBGC: W->S substitution rates on selected branch / S->W rates; report ratio. Re-run selection analysis with HyPhy BUSTED-MH (multi-hit aware); if signal vanishes, the original "selection" was gBGC + multi-hit substitutions. For genome-wide scans, mask sub-telomeric / high-recombination regions.

### Branch-site test foreground specification

**Trigger:** Running codeml branch-site mod A after looking at the data to choose foreground branch.

**Mechanism:** The branch-site test is designed for a single a priori foreground; post hoc specification inflates Type-I by ~5x because the choice was informed by the data.

**Symptom:** Branch-site test highly significant for the "interesting" branch; aBSREL on same data shows no significant branch (aBSREL has built-in multiple-testing correction).

**Fix:** Specify foreground branches in registered protocol before looking at data. For exploratory branch-wise analysis, use aBSREL (Smith 2015 MBE 32:1342) which adaptively assigns branch-specific omega classes with multiple-testing built in. If branch-site test was post hoc, apply Bonferroni correction across all branches tested + report explicitly.

### LRT critical value confusion

**Trigger:** Computing branch-site test p-value using standard chi-square df=2.

**Mechanism:** The branch-site test compares mod A (4 omega classes) against mod A1 (omega2 fixed at 1). The LRT statistic distribution is a 50:50 mixture of point-mass-at-0 and chi-square(df=1), not chi-square(df=2) (Self & Liang 1987 JASA 82:605; Zhang 2005 MBE 22:2472; Wong 2004 Genetics 168:1041). Using df=2 makes the test conservative; using df=1 standard makes it anticonservative.

**Symptom:** Branch-site p-values incorrectly inflated or deflated; users report finding selection at very stringent thresholds.

**Fix:** Use the 50:50 mixture critical value: 2.71 at p=0.05 (NOT 3.84). PAML's `chi2 1 LRT` command applies the mixture. Many published applications use chi-square df=2 conservatively, which loses power but doesn't inflate; chi-square df=1 directly is wrong and inflates Type-I.

### omega2 hitting upper bound (999)

**Trigger:** PAML codeml output shows omega2 = 999 for an "under selection" site class.

**Mechanism:** PAML codeml uses an internal upper bound of 999 (= "infinity" in single precision). Hitting it indicates numerical issue: extremely few synonymous sites in the selected class, dS underestimation, or numerical optimization failure.

**Symptom:** Sites flagged as positive selection have omega2 = 999; BEB posteriors for those sites are weirdly distributed.

**Fix:** Re-run with multiple starting values of omega (`fix_omega=0`, vary `omega = 0.1, 0.5, 1.0, 2.0, 5.0` across runs); check that all converge to same omega. Inspect alignment at flagged sites for unusual residue conservation. If omega = 999 persists, the gene may have rare-substitution patterns; switch to BUSTED-MH which accounts for multi-hit substitutions.

### Multiple-testing burden in genome scans

**Trigger:** Running selection tests across thousands of genes without correction.

**Mechanism:** With ~5000 protein-coding genes in a typical analysis, 250 will be significant at p=0.05 under H0. The false-discovery rate without correction is 50%.

**Symptom:** Implausibly large gene lists "under selection"; functional categories enriched are non-specific (e.g. all immune genes by FDR).

**Fix:** Apply FDR correction (Benjamini-Hochberg). Genes in syntenic regions are non-independent; use Benjamini-Yekutieli for stronger control under dependence. For HyPhy site-level methods, the per-site p < 0.1 default is a starting point; multiple-test correction within a gene is typically not applied (sites within a gene are dependent), but cross-gene correction is necessary. Holm-Bonferroni for strict Type-I.

### Convergent substitution misinterpreted as positive selection

**Trigger:** Lineage-specific selection found at a residue that has independently changed in multiple unrelated lineages.

**Mechanism:** Convergent substitutions at the same site in independent lineages produce signals in branch-site and other tests; this is convergence, not adaptive evolution per se (though convergent residues often ARE adaptive).

**Symptom:** Same residue flagged in multiple unrelated lineages by branch-site test; alignment shows convergent substitutions.

**Fix:** Switch from selection test to convergence test: CSUBST (Fukushima & Pollock 2023 Nat Eco Evo 7:155) for combinatorial substitution analysis; RERconverge (Redlich 2024 MBE 41:msae210) for relative-rate-vs-phenotype across categorical traits; PCOC (Rey 2018) for biophysical convergence. Report both convergence test and selection test results.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| dN/dS interpretation | omega < 1 purifying; omega = 1 neutral; omega > 1 positive (per site, branch, or gene depending on model) | Yang & Bielawski 2000 TREE 15:496; foundational |
| Branch-site test LRT critical value | 2.71 at p=0.05 (50:50 mixture chi^2) | Self-Liang 1987 JASA 82:605; Zhang 2005 MBE 22:2472 |
| Site-level p-value default | p <= 0.1 (FEL, MEME, FUBAR); FUBAR posterior >= 0.9 | Murrell 2012/2013; Datamonkey conventions |
| BEB posterior probability | >= 0.95 significant; >= 0.99 highly significant | Yang & Bielawski 2000 |
| dS upper limit for reliability | dS < 1.5 per branch; dS < 3 overall | Yang 2007 PAML manual |
| Minimum sequences for codeml | >= 8 with sufficient divergence | Anisimova et al 2001 MBE 18:1585 |
| Branch-site test minimum lineages | >= 20 in tree; >= 4 background branches | Yang 2007 |
| GARD breakpoint significance | p < 0.05 to partition alignment | Pond 2006; mandatory pre-screen |
| MK alpha threshold | alpha > 0 indicates adaptive substitutions; report 95% CI | Smith & Eyre-Walker 2002 |
| Asymptotic alpha minimum SFS density | >= 50 sites per frequency bin | Messer & Petrov 2013 |
| FDR genome-wide selection scan | q < 0.05 Benjamini-Hochberg | Standard |
| MEME minimum site-level support | p < 0.1; +/-3 sequences with substitutions | Murrell 2012 |
| aBSREL p-value | p < 0.05 (corrected by Holm-Bonferroni internally) | Smith 2015 |
| RELAX k interpretation | k < 1 relaxed; k > 1 intensified | Wertheim 2015 |
| Codon usage bias ENC | ENC < 35 high bias; consider effect on dS | Wright 1990 Gene 87:23 |
| W->S substitution ratio for gBGC | > 1.5 suggests gBGC | Galtier 2013 |
| BUSTED-MH multi-hit threshold | omega_DH > 1 indicates multi-hit pattern | Lucaci 2023 |
| HyPhy SRV (Synonymous Rate Variation) | Use BUSTED-S when dS varies across sites > 2x | Wisotsky 2020 |

## Selection Scan Standard Pipeline

**Goal:** Test all coding genes in a clade for evidence of positive selection, with full quality control.

**Approach:** Align with PRANK -> filter with PREQUAL -> pre-screen with GARD -> run BUSTED-MH (gene-wide) + MEME (sites) + aBSREL (branches); FDR-correct across genes; verify top candidates pass alignment / saturation / gBGC checks.

```bash
# Per-gene pipeline (parallelizable)
for og in orthogroups/*.fa; do
    base=$(basename $og .fa)

    # 1. Codon-aware MSA
    prank -d=$og -o=msa/$base.prank -codon -F

    # 2. Filter alignment errors (segment-level)
    PREQUAL -i msa/$base.prank.best.fas -o msa_filt/$base

    # 3. Recombination pre-screen
    hyphy gard --alignment msa_filt/$base.filtered --output gard/$base.json
    # If breakpoints found: partition and treat per-segment

    # 4. Gene-wide test (multi-hit aware)
    hyphy busted --alignment msa_filt/$base.filtered \
        --tree species_tree.nwk --output busted_mh/$base.json \
        --srv Yes --multiple-hits Double+Triple

    # 5. Site-level
    hyphy meme --alignment msa_filt/$base.filtered \
        --tree species_tree.nwk --output meme/$base.json

    # 6. Branch-level
    hyphy absrel --alignment msa_filt/$base.filtered \
        --tree species_tree.nwk --output absrel/$base.json
done

# 7. Aggregate and FDR
python aggregate_selection_scan.py busted_mh/ meme/ absrel/ > selection_results.tsv
```

```python
'''Aggregate genome-wide selection scan results; FDR-correct.'''
import json, glob, pandas as pd
from scipy.stats import false_discovery_control

def parse_busted(p):
    d = json.load(open(p))
    return {'p_value': d.get('test results', {}).get('p-value'),
            'LRT': d.get('test results', {}).get('LRT'),
            'omega_DH': d.get('fits', {}).get('Unconstrained model', {}).get('omega3')}

def count_meme_sig(p, alpha=0.1):
    d = json.load(open(p))
    mle = d.get('MLE', {}).get('content', {}).get('0', {})
    headers = [h[0] for h in d.get('MLE', {}).get('headers', [[]])]
    pi = headers.index('p-value') if 'p-value' in headers else -1
    return sum(1 for v in mle.values() if pi >= 0 and v[pi] < alpha)

rows = []
for path in glob.glob('busted_mh/*.json'):
    gene = path.split('/')[-1].replace('.json', '')
    rows.append({'gene': gene, **parse_busted(path),
                 'meme_sig_sites': count_meme_sig(f'meme/{gene}.json')})
df = pd.DataFrame(rows)
df['busted_fdr'] = false_discovery_control(df['p_value'].fillna(1.0), method='bh')
df['adaptive'] = (df['busted_fdr'] < 0.05) & (df['meme_sig_sites'] > 0)
df.sort_values('busted_fdr').to_csv('selection_results.tsv', sep='\t', index=False)
```

## PAML Branch-Site Test (Operational)

**Goal:** Test for episodic positive selection on a pre-specified foreground branch.

**Approach:** Mark foreground in newick (`#1`) -> codeml branch-site mod A vs A1 -> LRT against 50:50 mixture chi^2(0):chi^2(1).

```bash
# Mark foreground branch: use ete4 or manually
python -c "
from ete4 import Tree
t = Tree('species_tree.nwk', format=1)
target = t.search_nodes(name='target_species')[0]
target.name = target.name + ' #1'
print(t.write(format=1))
" > foreground.nwk

# Branch-site mod A (alternative)
cat > codeml_modA.ctl << 'EOF'
seqfile = alignment.phy
treefile = foreground.nwk
outfile = mod_A.mlc
runmode = 0
seqtype = 1
CodonFreq = 2
model = 2
NSsites = 2
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
RateAncestor = 1
cleandata = 0
EOF
codeml codeml_modA.ctl

# Null model A1 (omega_2 = 1)
cp codeml_modA.ctl codeml_modA1.ctl
sed -i 's/^omega = 0.4/omega = 1/' codeml_modA1.ctl
sed -i 's/^fix_omega = 0/fix_omega = 1/' codeml_modA1.ctl
sed -i 's/outfile = mod_A.mlc/outfile = mod_A1.mlc/' codeml_modA1.ctl
codeml codeml_modA1.ctl
```

```python
'''Branch-site test LRT with 50:50 mixture critical value.'''
from scipy.stats import chi2

def branch_site_lrt(lnL_alt, lnL_null):
    lrt = 2 * (lnL_alt - lnL_null)
    if lrt <= 0:
        return {'LRT': lrt, 'p_value': 0.5}
    # 50:50 mixture of chi^2(0) and chi^2(1)
    p = 0.5 * (1 - chi2.cdf(lrt, df=1))
    return {'LRT': lrt, 'p_value': p}
```

Foreground branch must be specified before viewing data; for genome-wide screens with no a priori branch, use aBSREL instead. Bayes Empirical Bayes (BEB) sites with posterior > 0.95 on positive-selection class are the per-site call.

## McDonald-Kreitman with Asymptotic Alpha

**Goal:** Estimate adaptive substitution rate alpha = 1 - (Ds Pn) / (Dn Ps), corrected for slightly deleterious bias.

**Approach:** Compute counts of synonymous and nonsynonymous polymorphisms (P) and divergences (D); fit asymptotic alpha by binning by minor-allele frequency and extrapolating.

```r
# Standard MK
mk_alpha <- function(Dn, Ds, Pn, Ps) {
    1 - (Ds * Pn) / (Dn * Ps)
}

# Asymptotic alpha via the Messer-Petrov 2013 web tool
# (https://benhaller.com/messerlab/asymptoticMK.html) or the impMKT R package
# (Murga-Moreno 2022 G3 12:jkac206) which wraps the asymptotic computation.
library(impMKT)
# Inputs: per-frequency-bin (Pn, Ps) plus genome-wide (Dn, Ds)
freq_bins <- seq(0.01, 0.5, 0.01)
pn_by_freq <- c(...)  # nonsyn polymorphism count per bin
ps_by_freq <- c(...)  # syn polymorphism count per bin
fit <- asymptoticMK(
    Dn = total_dn, Ds = total_ds,
    Pn = pn_by_freq, Ps = ps_by_freq,
    x = freq_bins
)
fit$alpha_asymptotic       # adaptive substitution rate, corrected
fit$alpha_original         # original MK (biased)
```

For full DFE inference (alpha + distribution of fitness effects), use polyDFE (Tataru-Bataillon 2019):

```bash
polyDFE -d data.txt -m C -i estimates.init -o output_basename
```

DFE-alpha (Eyre-Walker 2009) and GRAPES (Galtier 2016) are alternatives; GRAPES is most robust for genome-wide adaptive-substitution scans.

## RERconverge for Trait-Correlated Rate Shifts

**Goal:** Identify genes whose evolutionary rate correlates with a binary or categorical phenotype across the species tree.

**Approach:** Compute per-gene relative evolutionary rates -> correlate against phenotype -> Bonferroni or FDR-correct across genes.

```r
library(RERconverge)

# Read alignments and tree
trees <- readTrees('orthogroup_trees.txt', minSpecies = 10)
rer <- getAllResiduals(trees, useSpecies = species_names, transform = 'sqrt',
                       weighted = TRUE, scale = TRUE)

# Define binary phenotype (e.g., echolocation in mammals)
phen_paths <- foreground2Paths(c('Bat1', 'Bat2', 'Dolphin'), trees, clade = 'terminal')
phen_vec <- foreground2Tree(c('Bat1', 'Bat2', 'Dolphin'), trees, clade = 'terminal')

# Correlate
cors <- correlateWithBinaryPhenotype(rer, phen_paths, min.sp = 10, min.pos = 2,
                                      weighted = 'auto')
top_genes <- cors[order(cors$P), ][1:50, ]
```

For categorical traits (more than binary), Redlich 2024 MBE 41:msae210 extends RERconverge.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| codeml M8 significant, BUSTED null | M8 vs M7 inflated by relaxed constraint mimicking selection | Trust BUSTED; check M8a vs M8 instead (stricter null) |
| codeml branch-site significant, aBSREL null | Branch-site test foreground post hoc | aBSREL with built-in multiple-testing is correct; downgrade claim |
| BUSTED significant, MEME no sites | Episodic at sites BUSTED can't pinpoint; or basic BUSTED detected SRV not selection | Run BUSTED-S; if signal vanishes, was SRV; if persists, gene-wide episodic |
| MEME positive, FEL null | Episodic selection (MEME-specific) | Trust MEME for episodic; FEL only detects pervasive |
| Multiple tests positive at same site | High-confidence site under selection | Report; consider experimental validation |
| Test positive but PREQUAL flagged 20% of alignment | Alignment artifact | Re-filter (HmmCleaner); re-test; downgrade if positive site is in filtered region |
| Test positive but in high-recombination region | gBGC | W->S substitution test; if gBGC-attributable, downgrade |
| BUSTED-MH null where BUSTED significant | Multi-hit substitutions misattributed | Trust BUSTED-MH; original positive was multi-hit pattern |
| RELAX k > 1 with branch-site test null | Selection regime intensification (more purifying) | RELAX captures regime shift; branch-site missed because foreground different |
| Branch-site significant on Drosophila branch but no signal in mammals | Lineage-specific adaptation; or dS saturation in mammals | Inspect dS distribution; if mammals dS < 0.5 across branch, signal is real; if dS > 2, saturation explanation |
| asymptotic alpha < 0 | DFE has high deleterious load; or demographic violation | Check polyDFE / GRAPES with demographic correction |

**Operational rule for publication:** GARD pre-screen documented as negative + PREQUAL/HmmCleaner filtering applied + dS < 1.5 per branch + W->S ratio not elevated + BUSTED-MH significant (gene-wide) + MEME flags sites + aBSREL flags branches with consistent direction = publication-ready evidence. Single-method significance (especially M8 vs M7 alone) should be downgraded.

## Cohort Gotchas

- **Immune / MHC loci:** intra-genic recombination is high; GARD pre-screen mandatory; high apparent positive selection often reflects gene conversion between alleles, not adaptive change
- **Viral genomes:** rapid evolution + recombination + multi-hit substitutions common; BUSTED-MH and FUBAR-MH essential; use RDP5 for recombination detection
- **Plasmodium / Trypanosoma:** high codon-usage bias and multi-hit substitutions; use BUSTED-MH and BUSTED-S
- **Mammalian X-chromosome:** higher dS than autosomes (male-driven evolution); gBGC asymmetry by chromosome; reduce dS threshold for X-linked genes
- **Recent human / population genetics:** dS dramatically underestimated at recent divergence; use SFS-based methods (asymptotic alpha, polyDFE)
- **Convergent evolution traits (echolocation, marine):** RERconverge / CSUBST / PhyloAcc-noncoding designed for these; codon-based methods alone miss the convergent signal

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "GARD pre-screen?" | Yes; no breakpoints (or partitioned at p < 0.05 breakpoints); per-segment results consistent |
| "Alignment filtering?" | PRANK codon-aware MSA; PREQUAL segment filter applied; Guidance2 scores reported |
| "Saturation?" | dS distribution shown; max per-branch dS < 1.5; analysis restricted to subclades meeting this |
| "Branch-site test foreground post hoc?" | Foreground pre-registered OR exploratory analysis acknowledged + aBSREL with built-in multiple-testing used |
| "Multiple-testing correction?" | FDR (Benjamini-Hochberg) across genes; per-site within gene not corrected (dependence) |
| "gBGC?" | W->S substitution ratio not elevated; non-sub-telomeric; BUSTED-MH null in candidates rules out |
| "Multi-hit?" | BUSTED-MH used; if signal persists, robust to multi-hit confounder |
| "Why this LRT df?" | Branch-site test uses 50:50 mixture (Self-Liang 1987; Zhang 2005); critical value 2.71 at p=0.05 |
| "Sensitivity to model choice?" | Cross-validated PAML vs HyPhy; consistent across both; reported both p-values |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| codeml runs but rst file empty | `RateAncestor = 0` or path not writable | Set `RateAncestor = 1`; check output directory |
| codeml omega2 = 999 | Numerical pathology / saturated dS | Vary starting omega; reduce taxon sampling to dS < 2 |
| codeml LRT negative (-0.001) | Numerical noise at convergence | Round; treat as no signal; rerun with different starting values |
| HyPhy "tree branches don't match alignment" | Mismatched taxa names | Use exact same labels in tree and alignment |
| HyPhy MEME returns "no sites significant" | High alignment uncertainty; or no episodic selection | Re-filter alignment; try BUSTED-S for gene-wide signal |
| GARD takes forever | > 50 sequences | Reduce to representative subset; or use RDP5 for viral data |
| MK alpha negative | Demographic issue or DFE has many slightly deleterious | Use polyDFE / GRAPES with demography correction |
| RERconverge "too few species per gene" | Stringent default | Reduce `min.sp = 5`; document |
| CSUBST omega_C unstable | Few combinations; small clade | Need >= 5 clades for stable convergence estimate |
| PhyloAcc convergence failure | Insufficient lineages | Re-run with relaxed prior; check input MAF distribution |

## Tool Installation Notes

```bash
conda install -c bioconda paml hyphy gard prank prequal hmmcleaner
# RDP5: http://web.cbio.uct.ac.za/~darren/rdp.html
# MACSE V2: wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.07.jar
pip install ete4 pyhyphy csubst
Rscript -e "install.packages(c('asymptoticMK', 'polyDFE'))"
Rscript -e "remotes::install_github('nclark-lab/RERconverge')"
# polyDFE / GRAPES / DFE-alpha source binaries at respective github / bioconda channels
```

For genome-wide scans (> 5000 genes), parallelize per-gene analyses with Snakemake / Nextflow.

## References

- Yang Z 1997 CABIOS 13:555 (PAML codeml)
- Yang Z et al 2000 Genetics 155:431 (codon models M0-M8)
- Yang Z & Bielawski JP 2000 TREE 15:496 (codon model framework)
- Zhang J et al 2005 MBE 22:2472 (branch-site mod A); Wong WSW et al 2004 Genetics 168:1041 (LRT mixture); Self SG & Liang K-Y 1987 JASA 82:605 (LRT boundary)
- Swanson WJ et al 2003 Genetics 165:725 (M8a null); Bielawski JP & Yang Z 2004 J Mol Evol 59:121 (clade models)
- Anisimova M & Yang Z 2007 MBE 24:1219 (PAML comparison); Anisimova M et al 2003 Genetics 164:1229 (recombination FP); Anisimova M, Bielawski JP & Yang Z 2001 MBE 18:1585 (LRT power)
- Pond SLK et al 2006 MBE 23:1891 (GARD); Martin DP et al 2021 Virus Evol 7:veaa087 (RDP5)
- Kosakovsky Pond SL & Frost SDW 2005 MBE 22:1208 (FEL); Murrell B et al 2012 PLoS Genet 8:e1002764 (MEME); Murrell B et al 2013 MBE 30:1196 (FUBAR)
- Murrell B et al 2015 MBE 32:1365 (BUSTED); Wisotsky SR et al 2020 MBE 37:2430 (BUSTED-S); Lucaci AG et al 2023 bioRxiv 2022.12.02.518889 (BUSTED-MH)
- Smith MD et al 2015 MBE 32:1342 (aBSREL); Wertheim JO et al 2015 MBE 32:820 (RELAX)
- McDonald JH & Kreitman M 1991 Nature 351:652 (MK); Smith NGC & Eyre-Walker A 2002 Nature 415:1022 (alpha); Messer PW & Petrov DA 2013 PNAS 110:8615 (asymptotic alpha)
- Murga-Moreno J et al 2022 G3 12:jkac206 (impMKT); Tataru P & Bataillon T 2019 Bioinformatics 35:2868 (polyDFE); Eyre-Walker A & Keightley PD 2009 MBE 26:2097 (DFE-alpha); Galtier N 2016 PLoS Genet 12:e1005774 (GRAPES)
- Galtier N & Duret L 2007 Trends Genet 23:273 (gBGC); Capra JA et al 2013 PLoS Genet 9:e1003684 (gBGC genome-scale)
- Schneider A et al 2009 GBE 1:114 + Markova-Raina P & Petrov D 2011 GR 21:863 (alignment-error FP)
- Loytynoja A 2014 Methods Mol Biol 1079:155 (PRANK); Ranwez V et al 2018 MBE 35:2582 (MACSE V2); Whelan S et al 2018 Bioinformatics 34:3929 (PREQUAL); Di Franco A et al 2019 BMC Eco Evo 19:21 (HmmCleaner)
- Yang Z 2007 PAML manual; Álvarez-Carretero S et al 2023 MBE 40:msad041 (Beginner's Guide PAML)
- Kowalczyk A et al 2019 Bioinformatics 35:4815 + Redlich R et al 2024 MBE 41:msae210 (RERconverge)
- Fukushima K & Pollock DD 2023 Nat Eco Evo 7:155 (CSUBST); Hu Z et al 2019 MBE 36:1086 (PhyloAcc); Pollard KS et al 2010 GR 20:110 (phyloP); Rey C et al 2018 MBE 35:2296 (PCOC)

## Related Skills

- comparative-genomics/ortholog-inference - Single-copy ortholog alignments as input
- comparative-genomics/ancestral-reconstruction - Branch-specific ancestral sequence inference
- comparative-genomics/gene-tree-species-tree-reconciliation - Reconciled gene trees as PAML input
- alignment/multiple-alignment - PRANK / MACSE codon-aware MSA
- alignment/alignment-trimming - PREQUAL / HmmCleaner segment filtering
- phylogenetics/modern-tree-inference - Tree inference required for codeml
- population-genetics/selection-statistics - SFS-based alpha + DFE methods
- causal-genomics/heritability-partitioning - LDSC partition includes positive-selection annotations
- variant-calling/variant-annotation - Functional annotation of selected sites
