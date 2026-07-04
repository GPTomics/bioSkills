---
name: bio-comparative-genomics-gene-tree-species-tree-reconciliation
description: Reconcile gene trees against a species tree under probabilistic models of duplication, transfer, and loss (DTL) using ALE (Szöllősi 2013 amalgamated likelihood), GeneRax (Morel 2020 ML reconciliation), AleRax (Morel 2024 co-estimation), Whale.jl (Bayesian DTL+WGD), RANGER-DTL 2 parsimony, NOTUNG, ecceTERA, and Treerecs. Use when inferring ancestral gene-family content, distinguishing duplication from horizontal transfer from differential loss, rooting deep species trees from gene-content signals (STRIDE / Williams 2017 ALE-rooting), counting DTL events per branch, refining noisy gene trees against a species tree, modeling WGD events jointly with DTL, or producing publication-grade gene-family histories for phylogenomic / comparative analyses.
tool_type: cli
primary_tool: ALE
---

## Version Compatibility

Reference examples tested with: ALE 1.0+ (ssolo/ALE github), GeneRax 2.1.3+ (BenoitMorel/GeneRax), AleRax 1.2.0+ (BenoitMorel/AleRax; Morel 2024 Bioinformatics 40:btae162), Whale.jl 2.0+ (arzwa/Whale.jl), RANGER-DTL 2.0+ (Bansal lab; Bansal 2018 Bioinformatics 34:3214), NOTUNG 2.9.1.5+ (Stolzer 2012; Chen 2000), ecceTERA 1.2.5+, Treerecs 1.2+, IQ-TREE 2.3.6+, MrBayes 3.2.7+, BUSCO 5.7+, ete4 4.1.0+, BioPython 1.84+. Open Tree of Life and NCBI Taxonomy reference databases at 2024-Q3 minimum for species-tree-aware inference.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `ALEml_undated --help`, `ALEml --help` (dated), `generax --help`, `alerax --help`
- Julia: `using Whale; Whale.WhaleProblem`; `]status` for package versions
- Python: `pip show ete4`; `ete4 --help`

If code throws `species tree mismatch`, `gene tree taxa not in species tree`, or `MPI process pool failure`, these reconciliation tools share strict label-consistency requirements: species labels must match exactly across the species tree and gene trees (case-sensitive, no whitespace), and gene IDs typically encode species via prefix (`species|gene_id` separator convention). Use `sed` / `awk` normalization scripts before reconciliation.

# Gene Tree Species Tree Reconciliation

**"Where did this gene family come from, and what events shaped its history?"** -> Reconcile gene trees against species trees under explicit probabilistic models of duplication (D), horizontal transfer (T), and loss (L). The reconciliation framework converts gene-tree-species-tree discordance into a quantitative history of evolutionary events. Modern probabilistic methods (ALE, GeneRax, AleRax) **distinguish gene-tree-error-driven discordance from biological discordance** by integrating over gene-tree uncertainty -- a critical advance over parsimony reconciliation (NOTUNG, RANGER) which treats input gene trees as fixed and inflates duplication/loss counts from gene-tree noise (Boussau 2013 Genome Res 23:323; Morel 2020 MBE 37:2763).

- CLI: `ALEobserve` + `ALEml_undated` -- Bayesian DTL on a sample of gene trees (amalgamated likelihood)
- CLI: `generax` -- ML reconciliation; refines gene trees jointly with reconciliation
- CLI: `alerax` -- co-estimation of gene and species trees + DTL rates (Morel 2024)
- Julia: `using Whale` -- Bayesian DTL + WGD modeling
- CLI: `ranger-dtl` -- parsimony DTL with cost weights
- CLI: `notung` -- duplication-loss only (DL); user-friendly GUI; legacy

## Algorithmic Taxonomy

| Tool | Approach | Events modeled | Inference | Strength | Fails when |
|------|----------|----------------|-----------|----------|------------|
| ALE undated (Szöllősi 2013 Syst Biol 62:901) | Amalgamated likelihood over gene-tree distribution; species-tree-aware | D, T, L | Bayesian | Posterior over D/T/L events at every species-tree branch; integrates over gene-tree uncertainty | Requires gene-tree posterior sample (>= 100 bootstrap/UFBoot trees); slow for many families |
| ALE dated | Same as undated but uses time-calibrated species tree | D, T, L | Bayesian | Time-aware; better donor inference | Requires dated species tree (BEAST2 / RevBayes calibration) |
| GeneRax (Morel 2020 MBE 37:2763) | ML reconciliation + joint gene-tree refinement | D, T, L | ML | Faster than ALE; refines noisy gene trees; species-tree-aware | Less uncertainty quantification than ALE |
| AleRax (Morel 2024 Bioinformatics 40:btae162) | Co-estimation of gene tree, species tree, and DTL rates | D, T, L | Bayesian / ML hybrid | Gold standard 2024; corrects gene-tree-error feedback into species tree | Computationally heaviest; needs >= 20 species |
| Whale.jl (Zwaenepoel & Van de Peer 2019 MBE 36:1384) | Bayesian DTL + WGD via amalgamated likelihood | D, T, L, WGD | Bayesian (Turing.jl) | Native WGD modeling; modern Bayesian framework | Julia ecosystem dependency |
| RANGER-DTL 2.0 (Bansal 2018 Bioinformatics 34:3214) | Parsimony DTL with user cost weights (D-cost, T-cost, L-cost) | D, T, L | Parsimony | Fast; deterministic; many gene families per minute | Cost weights are user choices; results sensitive to costs |
| NOTUNG (Chen 2000 JCB 7:429; Stolzer 2012 Bioinformatics 28:i409) | Parsimony DL; HGT extension | D, L (optional T) | Parsimony | User-friendly GUI; widely used | DL-only by default; HGT extension less rigorous than ALE |
| ecceTERA (Jacox 2016 Bioinformatics 32:2056) | DTL on input set of trees; ILS extension | D, T, L (+ ILS) | Parsimony / DP | Handles ILS jointly with DTL | Less popular than ALE; smaller community |
| Treerecs (Comte 2020 Bioinformatics 36:4822) | Joint species tree + gene tree refinement | D, L | ML | Refines gene trees by species-tree constraint | No HGT; eukaryote-focused |
| DLCpar (Wu 2014 GR 24:475) | DLC parsimony for DL + coalescence (ILS) | D, L, C | Parsimony | Models ILS explicitly | No HGT; older |
| GraphDTL (Tofigh 2011) | Graph algorithm for DTL | D, T, L | Parsimony | Fast on small instances | Less used today |
| Phyldog (Boussau 2013 GR 23:323) | Joint species-tree-gene-tree DL with site-rate variation | D, L | Bayesian | Joint inference; refines gene trees | Bacteria-unfriendly; eukaryote-only |

Methodology evolves; verify the AleRax / ALE documentation before locking on a single approach. The probabilistic ALE / GeneRax / AleRax tools have largely superseded parsimony reconciliation for serious phylogenomic work; parsimony is fine for screening but not for publication-grade DTL inference.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Bacterial / archaeal phylogenomics, 50-500 genomes | GeneRax (refinement) -> ALE undated (posterior) | Two-stage: GeneRax refines, ALE provides posterior |
| Eukaryote DL inference, no HGT expected | NOTUNG (legacy) or Treerecs | DL is the dominant signal; HGT rare in animals |
| Mixed prokaryote/eukaryote with HGT | ALE undated | Probabilistic D/T/L; ALE-rooting (Williams 2017) for deep questions |
| Plant comparative genomics with WGD | Whale.jl | Native WGD modeling; Bayesian |
| Need uncertainty quantification | ALE or AleRax | Posteriors on every branch; ML methods give point estimates only |
| Need fastest possible per-family analysis | RANGER-DTL parsimony | Deterministic; multi-gene parallel |
| Co-estimate species tree from many gene families | AleRax | Modern gold-standard; corrects gene-tree-error |
| Root a deep species tree from DTL signal | ALE undated rooting (Williams 2017 method) | Root inference from D/T/L event distribution |
| Detect ancient HGT in archaea / bacteria | ALE undated | Probabilistic T detection at each branch; donor inferred |
| Identify ancestral gene family content | ALE; report origination events per branch | Posterior over presence/absence at internal nodes |
| Test specific HGT hypothesis (e.g. plant -> nematode) | ALE on filtered OG set; manual gene tree inspection | Quantitative T posterior |
| Distinguish HGT from differential gene loss | ALE event posteriors (T vs L on candidate branch) | Probabilistic ratio between alternatives |
| WGD detection alongside DTL | Whale.jl explicit WGD modeling | Joint inference; replaces post hoc Ks plotting |
| Refine noisy gene trees against species tree | GeneRax `--strategy SPR` | Species-tree-aware gene-tree refinement |
| Gene family birth-death modeling | See [[gene-family-evolution]] (CAFE5) | Reconciliation is per-family; CAFE5 is across families |
| Single gene of interest, single species | Manual gene-tree placement; ALE not needed | Reconciliation framework is genome-scale |

## Per-Tool Failure Modes

### Gene-tree-error feedback inflating duplications

**Trigger:** Using GeneRax or ALE with poorly-supported gene trees (low bootstrap, short alignments).

**Mechanism:** Noisy gene trees show spurious topology that, when reconciled, produces apparent duplications-followed-by-losses or transfers. The reconciliation framework cannot distinguish gene-tree noise from real DTL events; the output is biased toward more events (Boussau 2013).

**Symptom:** D + T + L event counts exceed reasonable rates (e.g. > 5 events per gene per Myr in eukaryotes); per-branch event posteriors are diffuse; ALE convergence (in `_uTs` files) is slow.

**Fix:** Use ALE (which integrates over gene-tree posterior sample) rather than GeneRax (which uses a single ML tree). For GeneRax users, provide UFBoot trees with high (`-B 1000`) bootstraps, and refine via `--strategy SPR`. AleRax (Morel 2024) co-estimates gene trees + species tree + DTL rates, addressing this feedback directly.

### Species labels and gene IDs mismatch

**Trigger:** Different naming conventions across gene trees, species tree, and orthology files.

**Mechanism:** Reconciliation tools require species labels in the species tree to match a defined prefix or suffix in each gene ID. Inconsistencies cause silent failures or partial reconciliation.

**Symptom:** ALE fails with "taxon not in species tree"; or runs but produces zero reconciliation events; or reconciliation matrix is sparse.

**Fix:** Strict normalization before reconciliation. ALE convention: gene IDs as `species|gene_id` (pipe separator); species labels match the species tree leaf names exactly. Pre-process all gene trees with:
```bash
for tree in gene_trees/*.nwk; do
    sed -i 's/_gene_/|/g' "$tree"   # adjust separator
done
```
Verify with `nw_labels -I species_tree.nwk` vs `nw_labels -I gene_trees/OG0000001.nwk` -- both species sets must be subsets of the species tree leaves.

### Cost-weight sensitivity in parsimony reconciliation (RANGER)

**Trigger:** Running RANGER-DTL with default costs (D=2, T=3, L=1).

**Mechanism:** Parsimony reconciliation minimizes total cost; the inferred D / T / L event count depends linearly on these costs. Default costs are biased toward favoring duplication-loss explanations over transfer.

**Symptom:** RANGER reports fewer transfers than ALE/GeneRax on the same data; sensitivity-analysis varies event counts dramatically with cost changes.

**Fix:** Run RANGER with cost sensitivity sweep: D in {1, 2, 3, 4}, T in {1, 2, 3, 4, 5}, L = 1. Report consensus events appearing in all cost combinations. Or move to probabilistic ALE/GeneRax which infers rates from data, not user costs. ecceTERA also allows cost-sweep mode.

### Root sensitivity in ALE undated

**Trigger:** ALE on a species tree with poorly-supported root.

**Mechanism:** ALE undated treats the species tree root as fixed; event posteriors at deepest branches depend on the root location. Wrong root flips D vs T inference at deep nodes.

**Symptom:** Running ALE under multiple candidate rootings (STRIDE, MAD, outgroup) produces qualitatively different event histories at deep branches.

**Fix:** Run ALE under multiple rootings; report only robust events. For deep phylogenomic questions, use ALE-rooting (Williams 2017 PNAS 114:E4602): run ALE under all possible rootings, choose the root maximizing the joint likelihood. AleRax can co-estimate the root, removing this issue.

### WGD events misattributed as duplications

**Trigger:** Reconciliation on a clade with known WGD (vertebrates 2R, fish 3R, salmonids Ss4R, plant lineages).

**Mechanism:** Standard DTL models (ALE, GeneRax) treat WGD as a series of individual duplications; the posterior at the WGD branch is dominated by D events but the joint event is whole-genome.

**Symptom:** D events at the WGD branch are 10-100x higher than other branches; many "duplications" cluster temporally.

**Fix:** Use Whale.jl which natively models WGD as a single event with a flexible retention parameter (Zwaenepoel 2019 MBE 36:1384). Otherwise, post hoc identify clusters of synchronized duplications and label as WGD.

### Saturation at deep timescales

**Trigger:** Reconciliation on extremely deep clades (>1 Gyr divergence in bacteria; >500 Myr in eukaryotes).

**Mechanism:** Gene families have undergone many cycles of D, T, L; the observed pattern is consistent with many DTL histories. Parameter identifiability is lost.

**Symptom:** ALE convergence requires hundreds of cycles; posterior on D/T/L rates is uninformative; branch event posteriors are diffuse.

**Fix:** Restrict to subclades with more recent divergence for quantitative DTL claims; for deep questions, qualitative event-class identification only (e.g. "transfers occurred along this branch" without precise count). Williams 2017 PNAS 114:E4602 demonstrates how ALE on deep archaeal phylogeny still resolves event class.

### MPI parallelization failures in GeneRax

**Trigger:** Running GeneRax on cluster with many gene families.

**Mechanism:** GeneRax uses MPI to parallelize per-family analysis; misconfigured MPI environment (wrong `srun`/`mpirun`/`mpiexec`, wrong allocator) silently runs serial or hangs.

**Symptom:** GeneRax progresses through few families per hour; cluster CPU usage shows only 1 core per node.

**Fix:** Verify MPI: `mpirun -n 8 hostname` should show 8 different hostnames or threads. Use `--per-family-rates` and proper MPI launch: `mpirun -n $SLURM_NTASKS generax ...`. Reduce `--threads` to 1 per family (let MPI handle parallelism). AleRax uses the same convention.

### ILS misattributed as transfers

**Trigger:** Reconciliation on rapidly radiated clade (incomplete lineage sorting expected).

**Mechanism:** ILS produces gene-tree-species-tree discordance indistinguishable from HGT at short internodes. ALE / GeneRax cannot separate them.

**Symptom:** Many "transfers" at short internal branches; transfer rate per branch correlates with branch length (more transfers at short branches).

**Fix:** Use ecceTERA or DLCpar which model coalescence (ILS) jointly with DTL; or pre-screen for ILS-likely loci via Dsuite ABBA-BABA (see [[introgression-detection]]); restrict ALE to gene families where ILS unlikely (long internodes). For phylogenomic-scale ILS, consider Phyldog (Boussau 2013) or the ASTRAL-Pro2 coalescent species tree.

### Multifurcations in the species tree

**Trigger:** Using a species tree with polytomies (unresolved nodes).

**Mechanism:** ALE / GeneRax / AleRax assume strictly bifurcating species trees; multifurcations break the inference.

**Symptom:** Tool fails with "polytomy detected" or runs but produces nonsensical results at multifurcating nodes.

**Fix:** Resolve polytomies before reconciliation via `ape::multi2di()` (random resolution) or with an outgroup-informed resolution (RAxML / ASTRAL-Pro2 on more data). Document the resolution.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| ALE gene-tree sample size | >= 100 bootstrap or UFBoot trees per family | ALE documentation; below this, posterior poorly resolved |
| Bacterial transfer rate (per gene per branch) | typically 0.001-0.05 in ALE inferences | Szöllősi 2013; varies clade |
| Eukaryote duplication rate | typically 0.0001-0.005 | eukaryote-specific convention; calibrate per clade |
| Branch-wise DTL event posterior | > 0.5 for "called" event | ssolo/ALE convention |
| Minimum gene families for AleRax | >= 100 families; >= 20 species | Morel 2024 |
| Minimum species for species-tree rooting via ALE | >= 30 species across the clade | Williams 2017 |
| GeneRax `--strategy` choices | EVAL only (no refinement), SPR (refinement), HYBRID (random + SPR) | GeneRax docs |
| RANGER cost weights default | D=2, T=3, L=1 (Bansal 2018) | Sensitivity sweep recommended |
| Whale.jl MCMC burn-in | >= 1000 samples; convergence by ESS >= 200 | Whale.jl docs |
| Reasonable runtime per family (ALE) | 1-30 minutes | Modern CPU; varies with gene-tree count |
| Reasonable runtime per family (GeneRax) | 0.1-5 minutes | Modern CPU; SPR is slower than EVAL |
| Maximum families per AleRax run | < 5000 (computational) | Above this, partition |
| Species labels case-sensitivity | strict; case mismatch = silent failure | Universal convention |
| Gene-ID separator | ALE / GeneRax accept a single-character separator via `separators="X"` (commonly `|` or `_`); Whale.jl uses a mapping file. Verify per-tool. | Tool-specific |
| Transfer branch detection minimum support | branch posterior > 0.5 + AU test on alternative placement | Conservative publication-grade |
| Distinguish T from L | T posterior - L posterior on candidate branch | Subtraction approach; ALE outputs both |

## ALE Standard Workflow

**Goal:** Quantify D/T/L events per branch of a species tree, integrating gene-tree uncertainty.

**Approach:** Build UFBoot gene trees per orthogroup -> `ALEobserve` to encode tree samples -> `ALEml_undated` for reconciliation -> aggregate per-branch event posteriors.

```bash
# 1. Build UFBoot gene trees per orthogroup
mkdir -p gene_trees
for og in orthogroups/*.fa; do
    base=$(basename $og .fa)
    iqtree2 -s $og -m TEST -B 1000 -nt 2 --prefix gene_trees/$base
done

# 2. Encode for ALE
for ufb in gene_trees/*.ufboot; do
    ALEobserve $ufb
done
# Produces gene_trees/*.ale files

# 3. Reconcile against species tree
mkdir -p reconciled
for ale in gene_trees/*.ale; do
    ALEml_undated species_tree.nwk $ale \
        separators="|" \
        sample=100 \
        output_format=newick
    mv ${ale%.*}_*.uml ${ale%.*}.uml
    mv ${ale%.*}_*.uTs ${ale%.*}.uTs
done

# 4. Aggregate branch-wise events
python aggregate_ale_events.py reconciled/ > branch_events.tsv
```

```python
'''Aggregate ALE outputs per species-tree branch.'''
import glob
from collections import defaultdict
import pandas as pd


def parse_uts(path):
    '''ALE _uTs format: branch  duplications  transfers  losses  originations  speciations.'''
    rows = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            parts = ln.split()
            if len(parts) >= 5:
                rows.append({
                    'branch_id': parts[0],
                    'duplications': float(parts[1]),
                    'transfers': float(parts[2]),
                    'losses': float(parts[3]),
                    'originations': float(parts[4]),
                })
    return pd.DataFrame(rows)


def aggregate(reconciled_dir):
    agg = defaultdict(lambda: defaultdict(float))
    for uts in glob.glob(f'{reconciled_dir}/*.uTs'):
        df = parse_uts(uts)
        family = uts.split('/')[-1].split('.')[0]
        for _, row in df.iterrows():
            for event in ('duplications', 'transfers', 'losses', 'originations'):
                agg[row['branch_id']][event] += row[event]
    return pd.DataFrame(agg).T.fillna(0)


branch_events = aggregate('reconciled')
print(branch_events.sort_values('transfers', ascending=False).head(20))
```

## GeneRax for ML Reconciliation with Refinement

**Goal:** Reconcile gene trees against a species tree while jointly refining noisy gene trees.

**Approach:** Provide gene trees + alignments + species tree -> GeneRax SPR strategy refines each gene tree -> reconciliation produces D/T/L history.

```bash
# Prepare families file (one line per family with paths)
cat > families.txt << 'EOF'
[FAMILIES]
- OG0000001
starting_gene_tree = gene_trees/OG0000001.nwk
alignment = alignments/OG0000001.fa
mapping = mapping.txt
subst_model = GTR+G
- OG0000002
starting_gene_tree = gene_trees/OG0000002.nwk
alignment = alignments/OG0000002.fa
mapping = mapping.txt
subst_model = GTR+G
EOF

# Mapping file: gene_id  species_label
# One line per gene
generate_mapping.py orthogroups.tsv > mapping.txt

# Run GeneRax with SPR refinement
mpirun -n 16 generax \
    --families families.txt \
    --species-tree species_tree.nwk \
    --rec-model UndatedDTL \
    --strategy SPR \
    --prefix generax_run \
    --per-family-rates \
    --max-spr-radius 5

# Output:
#   generax_run/results/<family>/inferredGeneTree.newick  refined gene tree
#   generax_run/results/<family>/reconciliation.nhx       reconciliation in NHX format
#   generax_run/species_trees/species_tree.newick         species tree (with optional re-rooting)
#   generax_run/families.txt                              refined families
```

## AleRax Co-Estimation

**Goal:** Co-estimate gene trees, species tree, and DTL rates from gene-family alignments.

**Approach:** Provide alignments + initial species tree -> AleRax runs joint Bayesian / ML co-estimation.

```bash
# Prepare AleRax families file (similar to GeneRax but with bootstrap trees)
cat > families.txt << 'EOF'
[FAMILIES]
- OG0000001
starting_gene_tree = gene_trees/OG0000001.ufboot
alignment = alignments/OG0000001.fa
mapping = mapping.txt
subst_model = GTR+G
EOF

mpirun -n 32 alerax \
    --families families.txt \
    --species-tree initial_species_tree.nwk \
    --rec-model UndatedDTL \
    --output alerax_run \
    --gene-tree-samples 100 \
    --species-tree-samples 1
```

AleRax outputs include `alerax_run/inferred_species_tree.nwk` (possibly re-rooted; D/T/L-aware) and per-family reconciled gene trees.

## Whale.jl Bayesian DTL + WGD

**Goal:** Bayesian DTL inference with explicit WGD modeling.

**Approach:** Define species tree with WGD nodes -> Whale.jl posterior over D/T/L + WGD retention.

```julia
# IMPORTANT: Whale.jl public API evolves; verify against current docs at
# https://github.com/arzwa/Whale.jl before scripting. The sketch below is
# illustrative; introspect via `?WhaleModel`, `?WhaleProblem`, `?read_ale`.
using Whale, NewickTree, Distributions, DynamicHMC

# Read species tree (Whale uses readnw or SlicedTree)
species_tree = readnw(read("species_with_wgd.nwk", String))

# Define rates model and WGD retention parameters
# Whale parameterization: λ (duplication), μ (loss), q (WGD retention per WGD node), ρ (sampling)
# The exact constructor signature is package-version-dependent; see ?WhaleModel
rates = ConstantDLWGD(λ=0.1, μ=0.1, q=Dict(1=>0.3), η=0.66)

# Read amalgamated gene-tree distributions (a directory of .ale files)
ale_dir = "families/"
ale_data = read_ale(ale_dir, species_tree)

# Build problem and sample with DynamicHMC (verify with current Whale.jl examples)
problem = WhaleProblem(ale_data, species_tree, rates)
results = mcmc(problem, n=1000)
```

**Operational note:** the Whale.jl 2.x DSL has moved between releases; do NOT copy-paste this without verifying with `?` on the actual installed version. Reference examples live in the upstream `examples/` directory.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ALE high transfer posterior, RANGER low | Cost weights in RANGER bias against T | Trust ALE; sensitivity-sweep RANGER costs |
| GeneRax fewer events than ALE | GeneRax uses ML on single tree; ALE integrates over tree distribution | Trust ALE; GeneRax may have underestimated due to gene-tree uncertainty |
| AleRax revises species tree from initial guess | Gene-tree-error feedback was biasing initial estimate | Trust AleRax (co-estimation corrects feedback) |
| Whale.jl WGD posterior high; ALE shows duplication burst | Same event; Whale models as WGD with retention parameter | Whale.jl interpretation is more biological for clades with known WGD |
| NOTUNG DL reconciliation contradicts ALE | NOTUNG ignores T; ALE includes T | Trust ALE for HGT-affected clades |
| Same family: D in one tool, T in another | Boundary case; both events possible | Report both with explicit caveats; or use AleRax for joint co-estimation |
| Many "transfers" at short internal branches | ILS confounded with T | Switch to ecceTERA / DLCpar / Phyldog for ILS-aware inference |
| ALE posterior diffuse across branches | Saturated; very ancient family or short tree | Restrict to subclade; or report event class only |
| GeneRax fails on 1% of families | Specific orthology / alignment issue | Inspect failed families manually; often related to species-label mismatch |

**Operational rule for publication:** ALE with 100+ bootstrap gene trees + Bayesian event posteriors > 0.5 + multiple-rooting robustness + biological corroboration (e.g. HGT predictions cross-checked with composition / [[hgt-detection]]) = publication-grade DTL inference. Single parsimony reconciliation (RANGER, NOTUNG) is appropriate for screening but should be backed by ALE for published claims.

## Cohort Gotchas

- **Bacterial clades with rampant HGT (e.g. Enterobacteriaceae, Streptomyces):** ALE transfer posteriors will dominate; calibrate expected per-branch transfer rate against published values (Szöllősi 2013; Williams 2017)
- **Endosymbionts with genome reduction (Buchnera, Wolbachia):** rampant gene loss; loss rates may exceed all other events; calibrate against known biology
- **Salmonid 4R Ss4R WGD:** add WGD node to species tree before reconciliation; use Whale.jl with native WGD parameter; treating Ss4R as duplication burst inflates DL counts
- **Plant tetraploids:** subgenome assignment first ([[whole-genome-duplication]]); reconcile each subgenome lineage separately
- **Rapid radiations (cichlids, Drosophila species groups, hominoids):** ILS confounded with transfers; use ASTRAL-Pro2 coalescent species tree as input; or apply ecceTERA with ILS extension
- **Polyploid lineages with WGD-derived "extra" genes:** appear as duplications in DTL methods; AleRax with WGD node specification handles this; Whale.jl native

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Gene-tree uncertainty?" | ALE integrates over 100+ UFBoot trees per family; AleRax co-estimates |
| "Species tree rooting?" | ALE-rooting (Williams 2017) under multiple candidate rooting hypotheses; events robust across rootings reported |
| "ILS vs HGT?" | Tested via ABBA-BABA / Dsuite; or used ecceTERA / DLCpar for joint inference |
| "WGD vs DL?" | Whale.jl native WGD modeling for known-WGD lineages; otherwise WGD node added explicitly |
| "Why ALE over RANGER?" | RANGER is cost-sensitive parsimony; ALE provides posterior probabilities |
| "Why AleRax over ALE?" | AleRax co-estimates gene tree + species tree + DTL rates, correcting gene-tree-error feedback; preferred for publication-grade work since 2024 |
| "ILS at rapid radiation?" | ILS-aware reconciliation (ecceTERA, DLCpar); or restrict to gene families with long internodes |
| "Cost weight sensitivity (RANGER)?" | If RANGER used, sensitivity sweep across cost weights performed; consensus events reported |
| "Contamination?" | Pre-filtered with FCS-GX / BlobTools (cross-ref [[hgt-detection]]) |
| "Gene-tree quality?" | UFBoot bootstrap >= 95 average; PREQUAL / HmmCleaner alignment filter applied |
| "Species sampling?" | At least 30 species for ALE-rooting; clade-balanced sampling reported |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| ALEml "taxon not found" | Species label mismatch | Normalize labels exactly; verify with `nw_labels` |
| GeneRax MPI hangs | Wrong MPI launcher; per-thread vs per-process confusion | Use `mpirun -n $NTASKS generax --threads 1`; check `mpirun -n 4 hostname` |
| ALE convergence slow / never | Insufficient gene-tree samples | Increase UFBoot samples to 1000; verify each `.ale` file has >= 100 trees |
| Whale.jl Turing sampling fails | NUTS step-size; tree dimension high | Try `HMC()` with manual tuning; or scale down number of WGD nodes |
| RANGER outputs single most-parsimonious history | One MPR returned by default | Run with `--explore-mpr` for multiple equally parsimonious histories |
| AleRax species tree differs from initial | Co-estimation revised it (intended) | Use AleRax's tree; report initial vs refined |
| Reconciliation gives 1000 events on a 5-gene family | Numerical instability or bad input | Inspect family for paralog confusion or sequence error |
| Per-branch event posterior > 1.0 | Multiple events on same branch | Expected for high-DTL-rate branches; sum over event classes |
| NOTUNG DL counts vastly exceed ALE D + L counts | NOTUNG attributes T events to L (no T model) | Use ALE for HGT-affected clades |
| ecceTERA "ILS detected at every branch" | ILS extension over-applied | Restrict ILS-aware to known short-internode regions only |

## Tool Installation Notes

```bash
# ALE (C++)
git clone https://github.com/ssolo/ALE && cd ALE && mkdir build && cd build && cmake .. && make
# Or via bioconda
conda install -c bioconda ale

# GeneRax (with MPI)
git clone --recursive https://github.com/BenoitMorel/GeneRax && cd GeneRax && ./install.sh

# AleRax
git clone --recursive https://github.com/BenoitMorel/AleRax && cd AleRax && ./install.sh

# Whale.jl (Julia)
julia -e 'using Pkg; Pkg.add("Whale")'

# RANGER-DTL
wget https://compbio.engr.uconn.edu/software/RANGER-DTL/RANGER-DTL-Linux.tar.gz
tar xf RANGER-DTL-Linux.tar.gz

# NOTUNG
wget https://www.cs.cmu.edu/~durand/Notung/download/Notung-2.9.1.5.tar.gz

# ecceTERA
git clone https://github.com/cchauve/ecceTERA && cd ecceTERA && make

# Treerecs
conda install -c bioconda treerecs

# Newick utilities (label inspection)
conda install -c bioconda newick_utils
```

For cluster deployment, AleRax / GeneRax require MPI; verify `mpicc --version` and run a small `hostname` test before deploying to genome-scale data.

## References

- Szöllősi GJ et al 2013 Syst Biol 62:901 (ALE undated)
- Szöllősi GJ et al 2015 Syst Biol 64:e42 (ALE rooting concept)
- Morel B et al 2020 MBE 37:2763 (GeneRax)
- Morel B et al 2024 Bioinformatics 40:btae162 (AleRax co-estimation)
- Williams TA et al 2017 PNAS 114:E4602 (ALE-rooting of the archaeal tree of life)
- Zwaenepoel A & Van de Peer Y 2019 MBE 36:1384 (Whale.jl Bayesian DTL+WGD)
- Bansal MS et al 2018 Bioinformatics 34:3214 (RANGER-DTL 2.0)
- Chen K et al 2000 J Comp Biol 7:429 (NOTUNG)
- Stolzer M et al 2012 Bioinformatics 28:i409 (NOTUNG-HGT extension)
- Jacox E et al 2016 Bioinformatics 32:2056 (ecceTERA)
- Comte N et al 2020 Bioinformatics 36:4822 (Treerecs)
- Wu Y-C et al 2014 GR 24:475 (DLCpar)
- Boussau B et al 2013 Genome Res 23:323 (Phyldog joint inference)
- Sjöstrand J et al 2012 Bioinformatics 28:2994 (PrIME-DLRS Bayesian)
- Tofigh A, Hallett M & Lagergren J 2011 IEEE/ACM TCBB 8:517 (DTL graph algorithm)
- Maddison WP 1997 Syst Biol 46:523 (gene tree discordance causes)
- Rabier C-E et al 2014 MBE 31:750 (rate inference under WGD)
- Tria FDK et al 2017 Nat Eco Evo 1:0193 (MAD rooting alternative)
- Emms DM & Kelly S 2017 MBE 34:3267 (STRIDE rooting)

## Related Skills

- comparative-genomics/hgt-detection - DTL reconciliation underlies probabilistic HGT inference
- comparative-genomics/ortholog-inference - Orthogroups feed reconciliation pipeline
- comparative-genomics/gene-family-evolution - CAFE5 birth-death across families (complementary to per-family reconciliation)
- comparative-genomics/whole-genome-duplication - WGD modeling in Whale.jl
- comparative-genomics/ancestral-reconstruction - DTL informs ancestral gene-content inference
- phylogenetics/modern-tree-inference - UFBoot bootstrap gene trees for ALE input
- phylogenetics/bayesian-inference - MrBayes / RevBayes alternative gene-tree posteriors
- phylogenetics/species-trees - ASTRAL-Pro2 coalescent species tree as ALE/AleRax input
- alignment/multiple-alignment - High-quality MSA precedes gene-tree inference
- alignment/alignment-trimming - PREQUAL / HmmCleaner for clean gene trees
