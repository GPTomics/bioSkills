---
name: bio-comparative-genomics-ancestral-reconstruction
description: Reconstruct ancestral states at internal phylogenetic nodes for sequences (PAML codeml, IQ-TREE --ancestral, GRASP, FastML), discrete traits (corHMM hidden-rate Markov, ape::ace, phytools::make.simmap stochastic mapping, BayesTraits), and continuous traits (phytools::fastAnc, geiger Brownian/OU, RPANDA). Use when designing constructs for ancestral protein resurrection, tracing trait evolution along a tree, performing stochastic character mapping, testing models of trait evolution (BM vs OU vs EB), inferring ancestral genome content via Dollo or DTL reconciliation, or quantifying ancestral-state uncertainty for downstream comparative analyses.
tool_type: mixed
primary_tool: PAML
---

## Version Compatibility

Reference examples tested with: PAML 4.10.7+, IQ-TREE 2.3.6+, GRASP 2024+ (web/CLI), FastML 3.11+, RevBayes 1.2.4+, BayesTraits V4.1+, R 4.4+, ape 5.8+, phytools 2.3+, corHMM 2.9+, geiger 2.0.11+, phangorn 2.12+, RERconverge 0.3.0+, BioPython 1.84+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('corHMM')` then `?ancRECON`, `?make.simmap`, `?ace`
- CLI: `codeml` (no `--version`; check `paml -h` or examine `Phylip.tre` example), `iqtree2 --version`, `revbayes --version`
- Python: `pip show biopython`; check `Bio.Phylo.PAML.codeml` API

If code throws `AttributeError`, `ImportError`, missing slot errors on R S4 objects, or PAML `mlc` parsing failures, introspect the installed package (`?` in R, `help()` in Python) and adapt the example rather than retrying. PAML output formats are stable across 4.9 -> 4.10; IQ-TREE's `--ancestral` flag replaced `-asr` in v2.0+.

# Ancestral State Reconstruction

**"What did this gene / trait / genome look like at an internal node?"** -> Choose the reconstruction framework that matches the data class (sequence / discrete trait / continuous trait / gene content) and the inference question (point estimate vs full posterior; marginal vs joint vs scaled-conditional). The single most common mistake is reconstructing under a site-independent or trait-stationary model when the underlying biology demands a hidden-rate or epistatic model -- the resulting "ancestor" is mathematically optimal under the wrong model and is silently wrong (Beaulieu & O'Meara 2016 Syst Biol 65:583; Boyko & Beaulieu 2021 MEE 12:468).

- Sequence ASR (protein resurrection): PAML codeml `RateAncestor=1`; IQ-TREE2 `--ancestral`; GRASP (graph-based, handles indels); FastML (Bayesian)
- Discrete traits: R `ape::ace(type='discrete')`; `corHMM::corHMM()` (rate categories); `phytools::make.simmap()` stochastic mapping; BayesTraits MultiState/Discrete
- Continuous traits: `phytools::fastAnc()`; `phytools::contMap()`; `geiger::fitContinuous(model='BM'|'OU'|'EB')`; RPANDA `fit_t_env()`
- Ancestral gene content (presence/absence): `ape::ace(type='discrete', model='ARD')`; Dollo parsimony in phangorn; ALE/GeneRax for full DTL (see [[gene-tree-species-tree-reconciliation]])

## Algorithmic Taxonomy

| Framework | Data class | Inference | Strength | Fails when |
|-----------|------------|-----------|----------|------------|
| ML marginal (codeml RateAncestor=1; IQ-TREE --ancestral) | Sequence / discrete | Site-by-site MAP + posterior | Per-site uncertainty; fast (Yang 1995 Genetics 141:1641; Pupko 2000 MBE 17:890) | Strong epistasis; site-independent assumption violated |
| ML joint (Pupko 2000 MBE 17:890; IQ-TREE marginal+joint output) | Sequence / discrete | Single most-likely joint history across all nodes | Internally consistent ancestral sequence | Loses per-site uncertainty; epistasis hidden |
| Stochastic mapping (Nielsen 2002 Syst Biol 51:729; Huelsenbeck 2003 Syst Biol 52:131; phytools::make.simmap) | Discrete | Full posterior over character histories along branches | Quantifies transition timing and rates per branch; supports posterior arithmetic | Long trees (mixing slow); rare-state biases |
| Bayesian MCMC (RevBayes, BayesTraits, MrBayes) | Sequence / discrete / continuous | Full posterior; supports model averaging | Honest uncertainty; rate-variable; hierarchical | Slow; convergence diagnostics required (ESS > 200) |
| Parsimony (Fitch 1971 Syst Zool 20:406; Dollo) | Discrete | MP states at nodes | Fast; assumption-light | Felsenstein-zone LBA artifact; biased toward fast change (Felsenstein 1978 Syst Zool 27:401) |
| Hidden Markov / hidden rates (corHMM; HiSSE; HMM) | Discrete | State + rate-class jointly | Captures rate heterogeneity across the tree; non-stationarity (Beaulieu 2013 Syst Biol 62:725) | Requires enough state changes to identify hidden rates |
| Threshold model (Felsenstein 2012 Am Nat 179:145; phytools::threshBayes) | Binary on continuous liability | MCMC on latent liabilities | Models polygenic / underlying-quantitative discrete traits | Slow MCMC; complex liability covariance |
| BM / OU / EB on continuous (geiger::fitContinuous) | Continuous | Phylogenetic regression on BM, OU mean-reverting, EB time-decay | Standard for body-size / niche-shape continuous traits | Model adequacy ignored (Boettiger 2012 Evolution 66:2240; Cooper 2016 Biol J Linn Soc 118:64) |
| Multi-rate BM / OUwie (Beaulieu 2012 Evolution 66:2369) | Continuous | Rate / optimum varies by clade or discrete regime | Models regime shifts; integrates with discrete trait history | Regime mismapping cascades to spurious rate differences |
| Phylogenetic generalized least squares (PGLS) | Continuous (multivariate) | Mean expected under BM; covariance from tree | Tests for correlation while controlling shared ancestry (Felsenstein 1985 Am Nat 125:1) | Strong evolutionary rate heterogeneity; non-BM trait |
| DTL reconciliation for gene content (ALE, GeneRax) | Gene tree / orthogroup | Ancestral gene presence + duplications/transfers/losses | Joint sequence + gene-content posterior | See [[gene-tree-species-tree-reconciliation]] |
| Indel-aware ASR (GRASP, FastML) | Sequence | Treats gaps as a separate process | Handles indel evolution explicitly; supports protein engineering | Slower; limited model families |

Methodology evolves; verify the latest `corHMM` / `phytools` vignettes before locking on a single approach. For continuous-trait macroevolution, consult Cooper 2016 model-adequacy reviews.

## Decision Tree by Experimental Scenario

| Scenario | Recommended method | Why |
|----------|---------------------|-----|
| Protein resurrection (~50-500 Myr divergences) | IQ-TREE2 `--ancestral` + GRASP indel reconstruction | Per-site marginal probabilities for alt-construct design; GRASP fixes indel ambiguity that PAML treats as missing data |
| Codon-level sequence ASR with selection inference | PAML codeml `RateAncestor=1`, model M0 (single omega), `seqtype=1` | Codon model native; integrates with branch reconstruction; produces `rst` with BEB-style site probs |
| Deep eukaryote / archaeal ASR (> 1 Bya) | Bayesian (RevBayes / PhyloBayes-MPI CAT-GTR) | Site-heterogeneous CAT model corrects compositional LBA (Szánthó 2023 Syst Biol 72:767); ML site-homogeneous models fail at this depth |
| Binary discrete trait with 5-30 taxa | `ape::ace(type='discrete', model='ARD')` + bootstrap | Standard for simple binary; ER/SYM/ARD model comparison via AIC |
| Binary discrete trait with 30+ taxa, suspected rate variation | `corHMM(rate.cat=2)` HMM | Hidden rates capture rate heterogeneity; mandatory if Beaulieu 2013 sensitivity test fails |
| State-dependent diversification (correlation with speciation/extinction) | HiSSE (Beaulieu & O'Meara 2016) NOT BiSSE | BiSSE has catastrophic Type-I rate when rate heterogeneity is misattributed (Rabosky & Goldberg 2015 Syst Biol 64:340); HiSSE is the required null |
| Multi-state with phylogenetic uncertainty | `phytools::make.simmap(nsim=1000)` over a tree distribution | Marginalize over tree + state uncertainty; report posterior probabilities |
| Continuous trait, single regime | `phytools::fastAnc()` + `contMap` | Fast BM ML reconstruction; visual continuous reconstruction along branches |
| Continuous trait, suspected regime shifts | OUwie or `bayou` (Uyeda 2014 Syst Biol 63:902) | Multi-optimum OU models infer optimum shifts and their tree positions |
| Binary trait expected to be polygenic underlying | Threshold model `phytools::threshBayes` | Models latent liability properly; binary -> continuous bridge (Felsenstein 2012) |
| Ancestral gene family content | DTL reconciliation (ALE / GeneRax) | See [[gene-tree-species-tree-reconciliation]]; full posterior over D/T/L events |
| Ancestral genome architecture (gene order) | AGORA (Muffato 2023 Nat Eco Evo 7:355); DeCoSTAR | Joint reconciliation + adjacency posterior |
| Convergent rate shifts in noncoding | PhyloAcc (Hu 2019 MBE 36:1086); Thomas 2024 update | Bayesian Markov model on conserved noncoding elements |
| Convergent amino-acid substitutions | CSUBST (Fukushima & Pollock 2023 Nat Eco Evo 7:155) | Combinatorial substitution ratio omega_C; null-corrected |
| Categorical trait correlated rates | RERconverge (Kowalczyk 2019; Redlich 2024 MBE 41:msae210) | Relative evolutionary rates linked to a binary or categorical phenotype |

## Per-Tool Failure Modes

### Long-branch attraction (LBA) at deep nodes

**Trigger:** Tree with two long terminal branches separated by a short internal branch; mixed amino-acid compositions across taxa.

**Mechanism:** Site-homogeneous models (LG, WAG, JTT) assume constant amino-acid equilibrium frequencies across the tree. When real compositions differ (e.g. thermophiles vs mesophiles), models underestimate the probability of convergent substitutions at compositionally-constrained sites, producing apparent shared derived states between long branches (Szánthó 2023 Syst Biol 72:767). The reconstructed ancestor at the deep node is biased toward whichever long-branch composition the model favors.

**Symptom:** Bootstrap support at the contested node remains high under site-homogeneous models but collapses under CAT-GTR / CAT-PMSF. Posterior predictive checks for compositional homogeneity reject the model (Foster 2004 Syst Biol 53:485).

**Fix:** Move to PhyloBayes-MPI with CAT-GTR or IQ-TREE2 with CAT-PMSF (`-m LG+C60+F+R` then `-ft <tree>` for posterior mean site frequencies). For ASR specifically, use ancestral reconstruction only when the model passes compositional adequacy. Slow-fast site removal (recoded amino acids; Susko & Roger 2007 MBE 24:2139) is an alternative.

### Epistasis breaking site-independent ASR

**Trigger:** Multiple sites in the same protein evolve under coupled constraints (compensatory pairs in RNA secondary structure; buried-residue covariance; allosteric networks).

**Mechanism:** ML/Bayesian ASR assumes sites are independent given the tree; the joint ancestral sequence is the product of per-site posteriors. Real proteins evolve through compensatory substitutions where a destabilizing mutation at site i is compensated by a mutation at site j (Pollock 2012 PNAS 109:E1352; Shah 2015 PNAS 112:E3226). The ML ancestral sequence can contain a never-tested combination of states.

**Symptom:** The reconstructed protein fails to fold or is non-functional when expressed; positions flagged ambiguous (P < 0.9) are non-random and cluster in 3D space when mapped to structure.

**Fix:** Use GRASP indel-aware reconstruction; design 4-8 alternative constructs varying ambiguous positions (P < 0.9), prioritizing residues that are structurally coupled to high-confidence ML states; experimentally test each construct; report the range of functional reconstructions, not the single ML sequence. Hochberg & Thornton 2017 Annu Rev Biophys 46:247 review epistasis strategies.

### Root placement error cascading to deep ancestors

**Trigger:** Trees rooted by midpoint, outgroup with very long branch, or `--prefix` auto-root.

**Mechanism:** Marginal ASR posteriors at internal nodes depend on the root's position because the root defines the time direction of substitution. A wrong root flips state inferences for deep nodes (especially when ancestral state is asymmetric, e.g. presence -> absence is more common than reverse).

**Symptom:** Re-rooting the tree changes the inferred ancestral state at the deepest node by > 0.2 posterior probability; STRIDE rooting (Emms 2017 MBE 34:3267) disagrees with outgroup rooting.

**Fix:** Run ASR over a set of candidate roots; report robust nodes (state invariant) and root-sensitive nodes separately. For phylogenomic-scale data, use STRIDE / MAD rooting (Tria 2017 Nat Eco Evo 1:0193) or ALE-rooting (Williams 2017) and document the rooting strategy.

### BiSSE false-positive in state-dependent diversification

**Trigger:** Testing whether a discrete trait influences speciation/extinction using BiSSE (Maddison 2007 Syst Biol 56:701).

**Mechanism:** BiSSE attributes ALL rate variation to the focal trait. When real diversification heterogeneity is caused by a hidden character correlated with the focal trait, BiSSE reports a spurious significant association (Rabosky-Goldberg 2015 Syst Biol 64:340 — ~40% Type-I rate at moderate trees).

**Symptom:** BiSSE LRT highly significant but biological mechanism unclear; HiSSE rejects BiSSE in favor of a hidden-state model with the focal trait neutral.

**Fix:** Run HiSSE as the required null model (Beaulieu & O'Meara 2016 Syst Biol 65:583). Report BiSSE only if HiSSE-null is rejected. For traits with deep clade structure, use FiSSE (Rabosky & Goldberg 2017 Evolution 71:1432) which is robust to model misspecification by design.

### Parsimony vs ML on asymmetric rates (Felsenstein bias)

**Trigger:** Trait has strongly asymmetric forward vs reverse rates (e.g. gene loss > gene gain).

**Mechanism:** Parsimony minimizes total changes, implicitly assuming symmetric rates. ML/Bayesian methods estimate the rate matrix from the data and reconstruct accordingly. Under strong asymmetry, parsimony over-reconstructs the rarer state at ancestors (Cunningham 1999 Syst Biol 48:665).

**Symptom:** Parsimony and ML reconstructions disagree at > 30% of nodes; ML rates fit AIC-better with ARD (all-rates-different) than ER (equal-rates).

**Fix:** Always run ER vs SYM vs ARD model comparison via `ape::ace(model='...')` AIC; use ARD when asymmetry is supported. For gene-content evolution, Dollo parsimony (gain rare, loss common) is often the better prior than equal-rates ML.

### Continuous-trait BM-only model with non-BM evolution

**Trigger:** Fitting `phytools::fastAnc()` (which assumes BM) to a trait with strong directional or stabilizing selection.

**Mechanism:** fastAnc returns the BM-MLE ancestral state, which is a weighted mean of descendant values with weights from the BM covariance matrix. If the trait evolved under OU (stabilizing), real ancestor values were closer to the optimum than fastAnc returns; if under EB (early burst), real ancestors were more variable than fastAnc returns.

**Symptom:** BM model fits with `geiger::fitContinuous(model='BM')` give AIC > 4 above OU or EB; phylogenetic signal Pagel's lambda < 0.5; Blomberg's K significantly < 1 (Blomberg 2003 Evolution 57:717).

**Fix:** Run `fitContinuous` with multiple models (BM/OU/EB/lambda/kappa/delta); use the best AIC model's ancestral reconstruction. For OU, use `OUwie::ace()`; for regime shifts, `bayou::bayou.mcmc`. Always report Pagel's lambda alongside ancestral estimates as a phylogenetic-signal indicator. Boettiger 2012 Evolution 66:2240 and Cooper 2016 Biol J Linn Soc 118:64 detail model-adequacy testing.

### Alignment error propagating to ASR

**Trigger:** Using `MUSCLE` or default MAFFT alignment on highly diverged sequences (< 30% identity).

**Mechanism:** Misaligned columns place non-homologous residues into the same column. ML ASR treats those residues as states of the same character, producing impossible ancestral inferences (Vialle 2018 MBE 35:1783).

**Symptom:** Ambiguous regions of the alignment correspond to low-confidence ASR sites; gappy columns dominate the low-confidence set; alignment scoring (TCS, Guidance2) marks the same regions as poorly aligned.

**Fix:** Filter alignment with HmmCleaner (Di Franco 2019 BMC Eco Evo 19:21) or PREQUAL (Whelan et al 2018 Bioinformatics 34:3929) before ASR. Segment-level filtering outperforms block-level filtering (Gblocks, trimAl) for downstream evolutionary inference. For ASR specifically, mask ambiguous columns (treat as missing) rather than removing them, to preserve coordinates.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| ASR site high confidence | posterior >= 0.95 | Standard convention (Yang 1995); above this treat state as fixed |
| ASR site moderate confidence | 0.80 <= posterior < 0.95 | Worth alternative-construct testing in resurrection studies |
| ASR site ambiguous | posterior < 0.80 | Design alternative constructs; cluster against structure |
| Pagel's lambda interpretation | lambda > 0.7 strong signal; 0.3-0.7 moderate; < 0.3 weak (ad-hoc operational convention; Pagel 1999 introduced lambda but did not prescribe these cutoffs) | Pagel 1999 Nature 401:877 (method); community convention (thresholds) |
| Blomberg K interpretation | K > 1 conserved; K = 1 BM; K < 1 weak signal | Blomberg 2003 Evolution 57:717 |
| Bootstrap support for ancestral clade | >= 70% before trusting the state at that node | Standard; below this, root-sensitivity tests required |
| MCMC ESS for Bayesian ASR | ESS >= 200 per parameter; ASRV at least 200 | RevBayes / Tracer convention; Lakner 2008 Syst Biol 57:86 |
| Stochastic mapping nsim | >= 1000 simulations per tree | Bollback 2006 BMC Bioinf 7:88; for asymmetric rates raise to >= 5000 |
| Tree depth limit for protein ASR | dS / branch length < 1.0 at deepest node | Above this, signal saturated; Yang 2007 PAML manual |
| Codon ASR minimum sequences | >= 8 with sufficient divergence (~0.5 substitutions/site total) | Operational convention; below this, codon-model parameters poorly constrained |
| Minimum taxa for binary trait ER vs ARD AIC | >= 20 tips; below this, rates often unidentifiable | Beaulieu 2016 |
| GRASP indel posterior threshold | >= 0.8 to call indel present at node | GRASP documentation; below this, both states tested experimentally |
| OUwie regime requires | >= 10 tips per regime | Beaulieu 2012 Evolution 66:2369; below this, optima unidentifiable |
| HMM rate categories | start with rate.cat=2; AIC compare against 1 | Boyko & Beaulieu 2021 MEE 12:468 |

## PAML codeml Ancestral Reconstruction

**Goal:** Reconstruct ancestral codon or protein states at internal nodes using ML under a stationary codon/protein model.

**Approach:** Build a codon-aware alignment (PRANK or MACSE) -> infer rooted ML tree with the same model intended for ASR -> create codeml control file with `RateAncestor=1` -> run codeml -> parse `rst` for per-site posteriors and ancestral sequences.

```python
'''PAML codeml ancestral reconstruction with site-confidence parsing'''

import subprocess
import re
import os
from collections import defaultdict


def write_codeml_ctl(alignment, tree, out_dir, seqtype='codon', model='M0'):
    '''Write codeml control file for ancestral reconstruction.

    seqtype: 1=codon, 2=aa, 3=codon translated
    For protein resurrection use seqtype=2 with model=3 (empirical+gamma).
    RateAncestor=1 produces rst file with per-site posteriors.
    '''
    if seqtype == 'codon':
        ctl = f'''
      seqfile = {alignment}
     treefile = {tree}
      outfile = {out_dir}/mlc
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
 RateAncestor = 1
    cleandata = 0
   getSE = 1
'''
    else:
        ctl = f'''
      seqfile = {alignment}
     treefile = {tree}
      outfile = {out_dir}/mlc
      runmode = 0
      seqtype = 2
        model = 3
   aaRatefile = lg.dat
       alpha = 0.5
        ncatG = 4
 RateAncestor = 1
    cleandata = 0
'''
    ctl_path = os.path.join(out_dir, 'codeml.ctl')
    open(ctl_path, 'w').write(ctl)
    return ctl_path


def parse_rst_posteriors(rst_file):
    '''Parse rst file. Returns per-node, per-site (state, prob) records.
    PAML rst section headers: "Prob distribution at node X" and "List of extant and reconstructed sequences".
    '''
    text = open(rst_file).read()
    node_posteriors = defaultdict(list)

    for block in re.finditer(r'Prob distribution at node (\d+)(.*?)(?=Prob distribution at node|\Z)', text, re.DOTALL):
        node = int(block.group(1))
        for ln in block.group(2).splitlines():
            m = re.match(r'\s*(\d+)\s+\S+:\s*(\w)\((\d\.\d+)\)', ln)
            if m:
                site, state, prob = int(m.group(1)), m.group(2), float(m.group(3))
                node_posteriors[node].append({'site': site, 'state': state, 'prob': prob})
    return node_posteriors


def summarize_node_confidence(node_posteriors, ambig_cutoff=0.8):
    '''Per-node confidence summary for protein-resurrection construct design.'''
    summary = {}
    for node, sites in node_posteriors.items():
        probs = [s['prob'] for s in sites]
        n_high = sum(p >= 0.95 for p in probs)
        n_amb = sum(p < ambig_cutoff for p in probs)
        summary[node] = {
            'mean_post': sum(probs) / len(probs),
            'frac_high': n_high / len(probs),
            'n_ambiguous': n_amb,
            'ambiguous_sites': [s['site'] for s in sites if s['prob'] < ambig_cutoff]
        }
    return summary
```

## IQ-TREE2 Marginal Reconstruction

**Goal:** Faster ASR with native model selection, for protein/DNA alignments where PAML is too slow.

**Approach:** Run IQ-TREE with `-m TEST` to pick model -> `--ancestral` produces `.state` table of per-site posteriors; output rooted by `-o <outgroup>`.

```bash
iqtree2 -s alignment.fasta -m MFP --ancestral -o outgroup_taxon -B 1000 -nt AUTO --prefix asr_iqtree
# .state columns: Node | Site | State | p_A | p_C | p_G | p_T  (DNA)
# or            : Node | Site | State | p_A p_R p_N ...        (protein)
```

```python
import pandas as pd

def load_iqtree_state(state_file):
    '''Returns dict[node] -> DataFrame[site, state, max_post, all_post].'''
    df = pd.read_csv(state_file, sep='\t', comment='#')
    df.columns = [c.strip() for c in df.columns]
    state_cols = [c for c in df.columns if c.startswith('p_')]
    df['max_post'] = df[state_cols].max(axis=1)
    return df.groupby('Node')
```

## Stochastic Mapping (phytools::make.simmap)

**Goal:** Sample full character histories along branches for a discrete trait; quantify ancestral state probabilities with proper uncertainty.

**Approach:** Fit Mk rate matrix -> sample `nsim` simulated maps under the posterior -> summarize state probabilities per node and transitions per branch.

```r
library(phytools)
library(corHMM)

tree <- read.tree('species_tree.nwk')
traits <- read.csv('traits.csv', row.names = 1)
x <- setNames(traits$state, rownames(traits))

# Fit Mk model and pick ARD vs SYM vs ER by AIC
fit_er <- fitMk(tree, x, model = 'ER')
fit_sym <- fitMk(tree, x, model = 'SYM')
fit_ard <- fitMk(tree, x, model = 'ARD')
aic <- sapply(list(fit_er, fit_sym, fit_ard), AIC)
best_model <- c('ER', 'SYM', 'ARD')[which.min(aic)]

# Stochastic mapping under best model, nsim >= 1000 (raise for asymmetric rates)
smaps <- make.simmap(tree, x, model = best_model, nsim = 1000, pi = 'fitzjohn')
node_pp <- summary(smaps)$ace          # posterior probabilities per state per node

# Hidden-rate alternative for clade-rate heterogeneity
hmm_fit <- corHMM(phy = tree, data = data.frame(species = names(x), trait = x),
                   rate.cat = 2, model = 'ARD', node.states = 'marginal')
hmm_fit$states                          # marginal ancestral state matrix (rows: nodes)
```

`pi='fitzjohn'` uses the Fitzjohn 2009 Syst Biol 58:595 root prior (root state proportional to its equilibrium probability), preferable to `pi='estimated'` which can over-fit, or `pi='equal'` which can bias toward the rarer state when data are asymmetric.

## Continuous-Trait ASR with Model Adequacy

**Goal:** Reconstruct ancestral values for a continuous trait while honestly reporting which model class the data support.

**Approach:** Fit BM/OU/EB/lambda models -> AIC compare -> reconstruct under best model -> report Pagel's lambda and Blomberg's K as signal quality.

```r
library(phytools); library(geiger)

trait <- setNames(traits$mass, rownames(traits))
tree_p <- multi2di(tree)               # phytools/geiger require fully bifurcating

models <- c('BM', 'OU', 'EB', 'lambda', 'kappa', 'delta')
fits <- lapply(models, function(m) fitContinuous(tree_p, trait, model = m))
aic_tbl <- data.frame(model = models, AIC = sapply(fits, function(f) f$opt$aic))
best <- aic_tbl$model[which.min(aic_tbl$AIC)]

# Phylogenetic signal
lambda_fit <- phylosig(tree_p, trait, method = 'lambda', test = TRUE)
K_fit <- phylosig(tree_p, trait, method = 'K', test = TRUE)

if (best == 'BM') {
    asr <- fastAnc(tree_p, trait, CI = TRUE)        # BM ML ancestral with 95% CI
} else if (best == 'OU') {
    # Use OUwie for proper OU ancestral with regime
    library(OUwie)
    df <- data.frame(species = names(trait), regime = rep(1, length(trait)), trait = trait)
    asr_ou <- OUwie.anc(OUwie(tree_p, df, model = 'OU1'), data = df)
} else {
    # phytools::contMap can use lambda/EB rescaled tree
    asr_contmap <- contMap(tree_p, trait, method = 'fastAnc', plot = FALSE)
}
```

Report Pagel's lambda alongside ancestral values: lambda > 0.7 indicates BM-like, signal supports tree-based reconstruction; lambda < 0.3 indicates ecology rather than phylogeny drives variance and ASR is poorly constrained.

## GRASP Indel-Aware Sequence ASR (Protein Resurrection Workflow)

**Goal:** Reconstruct ancestral protein sequence INCLUDING indel states, for experimental resurrection.

**Approach:** GRASP (Foley 2022 PLoS Comp Biol 18:e1010633) uses a partial-order alignment graph to model indels probabilistically; outputs alternative reconstructions and explicit indel uncertainty.

```bash
grasp -aln alignment.fasta -tree species.nwk -out grasp_run --inference joint --threads 8
# Produces:
#   grasp_run/asr.fasta      ancestral sequences at every internal node
#   grasp_run/asr_posterior  per-site, per-state posterior
#   grasp_run/asr_indels     indel posterior per node per gap-block
```

**Construct-design protocol:** (1) Take ML sequence as primary construct. (2) For each site with max posterior < 0.8 and a second state with posterior > 0.2, build a single-mutant alternative. (3) For each indel block with posterior < 0.8, build a present and an absent alternative. (4) Order constructs by structural compactness (avoid surface loops first). Typical batch: 4-8 constructs per ancestral node. Hochberg & Thornton 2017 Annu Rev Biophys 46:247 review the operational pipeline.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Parsimony and ML disagree at > 30% nodes | Strongly asymmetric rates | Trust ML with ARD model after AIC support; check Felsenstein 1978 LBA zone |
| BiSSE highly significant, HiSSE neutral | Hidden trait drives diversification | Report HiSSE as primary (Rabosky-Goldberg 2015) |
| Codeml marginal and joint disagree at deep node | Strong epistasis or model misspecification | Test alternative constructs; consider GRASP and BAli-Phy joint inference |
| Site-homogeneous LG agrees with PhyloBayes CAT-GTR at all nodes | No deep compositional heterogeneity | Site-homogeneous is fine; report both as confirmation |
| Site-homogeneous LG and CAT-GTR disagree at deep node | Compositional LBA | CAT-GTR is correct; report CAT-PMSF for fast follow-up |
| Marginal and stochastic-map states disagree | Asymmetric rate matrix; root prior matters | Use `pi='fitzjohn'`; trust stochastic map for asymmetric cases |
| BM ASR vs OU ASR disagree on direction of trait change at root | Trait under stabilizing selection | OU model is correct if AIC supports; reconstruct toward the optimum, not the BM weighted mean |
| Reconstruction at deepest node flips under re-rooting | Insufficient outgroup support; LBA | Run STRIDE / MAD rooting; report deep-node state with uncertainty |
| GRASP indel and PAML "missing-data" reconstructions disagree | PAML ignores indels | GRASP is correct for protein resurrection; PAML codon-only is fine for selection analysis |

**Operational rule for publication:** Reconstruct under multiple models (at minimum: ER/SYM/ARD for discrete; BM/OU/EB for continuous; site-homogeneous + CAT-PMSF or PhyloBayes for sequence ASR at deep nodes); report only nodes whose state is invariant across models OR explicitly flag model-sensitive nodes. Single-model claims should be downgraded.

## Cohort Gotchas

- **Polyploid species in continuous-trait analyses:** body size, genome size, gene count are confounded with ploidy; assign subgenomes (see [[whole-genome-duplication]]) and treat as separate tips, or use multilabel-tree methods.
- **Hybrid taxa break the bifurcating-tree assumption:** ape and phytools assume strict bifurcations; for hybrids, use phylogenetic networks (phangorn, RevBayes admixture) or remove hybrid tips before ASR.
- **Tip-dated trees from molecular clock require careful root prior:** RevBayes / BEAST2 with calibrated tip dates produce trees in absolute time; PAML/IQ-TREE work in relative substitutions. Match the time unit when integrating across tools.
- **OrthoFinder species trees from gene-tree summary (STAG):** branch lengths are coalescent-units when used for ASR; convert to substitutions via concatenated alignment if downstream tools require it.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why marginal not joint reconstruction?" | Marginal exposes per-site uncertainty necessary for resurrection construct design; joint is internally consistent but hides ambiguity (Pupko 2000 MBE 17:890) |
| "How was epistasis handled?" | Designed and tested N alternative constructs at ambiguous (P < 0.8) sites; report functional range, not just ML sequence (Hochberg & Thornton 2017) |
| "Why these models?" | AIC compared ER/SYM/ARD for discrete; BM/OU/EB/lambda for continuous; site-homogeneous + CAT-PMSF for deep sequence ASR; reported ancestral state only at model-invariant nodes |
| "Phylogenetic signal?" | Pagel's lambda = X; Blomberg's K = Y; signal supports tree-based ASR (or: signal weak, ASR exploratory only) |
| "Effect of rooting?" | Reconstructed under multiple rootings; state at root invariant across STRIDE / MAD / outgroup, or explicit caveat for root-sensitive nodes |
| "Multiple-testing across nodes?" | Ancestral states are estimated quantities, not tested hypotheses; report posteriors per node, no multiplicity correction needed |
| "Why not BiSSE for trait-diversification correlation?" | HiSSE replaces BiSSE per Rabosky-Goldberg 2015; BiSSE Type-I rate ~40% on simulated data |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `codeml` rst file missing posteriors section | `RateAncestor` set to 0 or not set | Set `RateAncestor=1` in control file |
| IQ-TREE `--ancestral` empty `.state` file | Tree not rooted | Specify outgroup with `-o`; IQ-TREE requires rooting for ancestral output |
| `make.simmap` chains never mix | Asymmetric rate matrix with sparse data | Increase nsim to 5000+; switch root prior to `pi='fitzjohn'` |
| `ape::ace` fails with `NA/NaN/Inf in foreign function call` | Polytomies in tree | `multi2di()` to resolve to bifurcating; or use `ape::ace(method='ML')` with phytools |
| GRASP runs but produces empty asr.fasta | Alignment includes stop codons or X characters | Strip non-canonical residues; check that the alignment passes Bio.SeqIO validation |
| `fitContinuous(model='OU')` returns lambda = 0 | OU collapsed to white noise (no signal) | Trait variance unexplained by tree; reconsider whether continuous-trait ASR is meaningful |
| Stochastic mapping returns all-or-nothing at deep nodes | Insufficient signal; pi='equal' default | Use `pi='fitzjohn'`; check that taxa span both states |
| HiSSE convergence failures | `bound.par` too restrictive; hessian singular | Use `starting.vals=NULL`, restart with `output.liks=TRUE`; switch to `BiSSE-ness` (Magnuson-Ford & Otto 2012 Am Nat 180:225) as fallback |
| codeml omega = 0 or 999 at branch | Saturation or alignment artifact | Increase model complexity (M0 -> M3); check dS at branch; if dS > 3, ASR unreliable on that branch |
| Ancestral genome content reconstruction inflated | Assembly fragmentation produces false absences | Use BUSCO-completeness-corrected presence/absence; or run [[gene-tree-species-tree-reconciliation]] which handles loss-vs-missing |

## Tool Installation Notes

```bash
# CLI
conda install -c bioconda paml iqtree
# RevBayes via source build (https://revbayes.github.io/download) or homebrew on macOS
# GRASP from https://github.com/bodenlab/GRASP (Java)
# FastML web at fastml.tau.ac.il or CLI binary

# R
install.packages(c('ape', 'phytools', 'geiger', 'corHMM', 'phangorn', 'OUwie', 'bayou', 'RERconverge'))
remotes::install_github('thej022214/hisse')

# Python
pip install biopython ete4
```

For protein resurrection workflows, GRASP is the modern standard; for selection-context codon ASR, PAML remains the reference; for trait macroevolution, phytools + corHMM + OUwie are the working tier; for Bayesian rigor with model uncertainty, RevBayes is required.

## References

- Yang Z et al 1995 Genetics 141:1641 (marginal ASR likelihood framework)
- Pupko T et al 2000 MBE 17:890 (joint ASR efficient algorithm)
- Nielsen R 2002 Syst Biol 51:729 (stochastic mapping)
- Huelsenbeck JP et al 2003 Syst Biol 52:131 (Bayesian stochastic mapping)
- Felsenstein J 1985 Am Nat 125:1 (phylogenetic independent contrasts)
- Felsenstein J 2012 Am Nat 179:145 (threshold model)
- Pagel M 1999 Nature 401:877 (lambda phylogenetic signal)
- Blomberg SP et al 2003 Evolution 57:717 (K statistic)
- Beaulieu JM et al 2013 Syst Biol 62:725 (hidden Markov state-rate decoupling)
- Beaulieu JM & O'Meara BC 2016 Syst Biol 65:583 (HiSSE)
- Rabosky DL & Goldberg EE 2015 Syst Biol 64:340 (BiSSE Type-I rates)
- Boyko JD & Beaulieu JM 2021 MEE 12:468 (generalized HMM corHMM)
- Bollback JP 2006 BMC Bioinf 7:88 (SIMMAP)
- Pollock DD et al 2012 PNAS 109:E1352 (compensatory epistasis)
- Shah P et al 2015 PNAS 112:E3226 (contingency and entrenchment epistasis)
- Hochberg GKA & Thornton JW 2017 Annu Rev Biophys 46:247 (ASR for protein resurrection)
- Foley G et al 2022 PLoS Comp Biol 18:e1010633 (GRASP indel-aware ASR)
- Szánthó LL et al 2023 Syst Biol 72(4):767–780 (compositional LBA; CAT-PMSF) — DOI 10.1093/sysbio/syad013
- Boettiger C et al 2012 Evolution 66:2240 (model adequacy for continuous-trait macroevolution)
- Cooper N et al 2016 Biol J Linn Soc 118:64 (cautionary note OU)
- Cunningham CW 1999 Syst Biol 48:665 (asymmetric rate parsimony bias)
- Felsenstein J 1978 Syst Zool 27:401 (long branch attraction)
- Maddison WP et al 2007 Syst Biol 56:701 (BiSSE)
- Beaulieu JM et al 2012 Evolution 66:2369 (OUwie)
- Uyeda JC & Harmon LJ 2014 Syst Biol 63:902 (bayou)
- FitzJohn RG 2009 Syst Biol 58:595 (root prior)
- Whelan S et al 2018 Bioinformatics 34:3929 (PREQUAL)
- Di Franco A et al 2019 BMC Eco Evo 19:21 (HmmCleaner)
- Emms DM & Kelly S 2017 MBE 34:3267 (STRIDE rooting)
- Tria FDK et al 2017 Nat Eco Evo 1:0193 (MAD rooting)
- Hu Z et al 2019 MBE 36:1086 (PhyloAcc)
- Fukushima K & Pollock DD 2023 Nat Eco Evo 7:155 (CSUBST)
- Muffato M et al 2023 Nat Eco Evo 7:355 (AGORA)

## Related Skills

- comparative-genomics/positive-selection - Branch- and site-level selection inference on ancestral branches
- comparative-genomics/ortholog-inference - Define orthogroups whose alignments feed ASR
- comparative-genomics/gene-tree-species-tree-reconciliation - DTL-aware ancestral gene-content inference; root inference via ALE
- comparative-genomics/whole-genome-duplication - Ks-dating provides time scale for ancestral state inference
- comparative-genomics/comparative-annotation-projection - Project ancestral CDS to descendants for validation
- phylogenetics/modern-tree-inference - Generate rooted ML/Bayesian trees as ASR scaffold
- phylogenetics/bayesian-inference - RevBayes / MrBayes priors for Bayesian ASR
- phylogenetics/divergence-dating - Time-calibrated trees as input for absolute-time ASR
- alignment/multiple-alignment - PRANK / MACSE indel-aware alignment before sequence ASR
- alignment/alignment-trimming - PREQUAL / HmmCleaner filtering before ASR
