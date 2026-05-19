---
name: bio-comparative-genomics-introgression-detection
description: Detect introgression and admixture between species or populations using Dsuite (Malinsky 2021 fast D-statistics), Patterson's D / ABBA-BABA test (Green 2010; Durand 2011), f4-ratio and f-branch statistic (Malinsky 2018), TreeMix (Pickrell & Pritchard 2012), HyDe (Blischak 2018), QuIBL (Edelman 2019), sprime (Browning 2018), Twisst (Martin 2017), PhyloNet (Solis-Lemus 2017) for explicit phylogenetic networks, and qpAdm / qpGraph (Patterson 2012; Lipson 2013). Distinguish introgression from incomplete lineage sorting (ILS), ancestral structure, ghost-lineage admixture, and rate variation. Use when testing inter-species gene flow, dating admixture events, identifying introgressed segments, building phylogenetic networks for reticulate evolution, or applying the ABBAclustering (Koppetsch-Malinsky-Matschiner 2024) framework for divergent-species gene flow.
tool_type: cli
primary_tool: Dsuite
---

## Version Compatibility

Reference examples tested with: Dsuite 0.5+ (millanek/Dsuite; Malinsky 2021 Mol Ecol Res 21:584; ABBAclustering option from Koppetsch-Malinsky-Matschiner 2024 Syst Biol), HyDe 0.4.3+ (Blischak 2018 Syst Biol 67:821), QuIBL (Edelman 2019 Science 366:594), TreeMix 1.13+ (Pickrell & Pritchard 2012 PLoS Genet 8:e1002967), sprime (Browning 2018 Cell 173:53), Twisst (Martin & Van Belleghem 2017 Genetics 206:429), PhyloNet 3.8.2+ (NakhlehLab/PhyloNet; Than-Ruths-Nakhleh 2008 BMC Bioinf 9:322) and PhyloNetworks 0.16+ (JuliaPhylo/PhyloNetworks; Solis-Lemus & Ane 2017 PLoS Comp Biol 13:e1005485), qpAdm / qpGraph (AdmixTools v2.0+; Maier 2023), ADMIXTOOLS2 R wrapper (Maier 2023 eLife 12:e85492), MaCS-like simulators (msprime 1.3+ for testing), bcftools 1.21+, samtools 1.21+, vcftools 0.1.16+, R 4.4+. See upstream Dsuite docs for visualization helpers.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `Dsuite --version`; `treemix --help`; `qpDstat --help` (AdmixTools)
- Python: `pip show msprime hyde`
- R: `packageVersion('admixtools')`

If code throws `Dsuite SETS file format error`, `TreeMix matrix singular`, `qpAdm rotation failed`, these tools share strict input requirements: Dsuite needs a SETS file mapping samples to populations + an outgroup; TreeMix needs allele-frequency matrix; qpAdm needs ind / snp / geno trio (EIGENSTRAT format).

# Introgression and Admixture Detection

**"Has there been gene flow between these species / populations?"** -> Tests for inter-population admixture span site-frequency (ABBA-BABA, f4), tree-topology (Dsuite f-branch, QuIBL, Twisst), explicit-network (PhyloNet), and gene-trajectory (sprime, hapne) approaches. The fundamental confounder is **incomplete lineage sorting (ILS)**: under a symmetric tree with no gene flow, the ABBA and BABA patterns occur with equal frequency from ancestral polymorphism. **A significant D-statistic indicates EITHER (a) introgression, OR (b) ancestral structure, OR (c) sampling from a ghost lineage** -- additional evidence is required to distinguish (Green 2010 Science 328:710; Durand 2011 MBE 28:2239; Eriksson & Manica 2012 PNAS 109:13956). For high-confidence claims, combine D-statistic with **f-branch** mapping (assigns admixture to specific branches), Twisst / QuIBL topology-weighting, and TreeMix migration edges.

- CLI: `Dsuite Dtrios` -- standard D-statistic for trios of populations
- CLI: `Dsuite Fbranch` -- assign admixture signal to specific tree branches
- CLI: `treemix -i freq.gz -k 1000 -m 0 -o out` -- migration edges
- CLI: `qpAdm` (AdmixTools) -- rotation-based admixture proportion testing
- R: `admixtools::qpadm()` modern wrapper
- CLI: `HyDe -i sequences -t taxa.txt --hyptest` -- hybridization detection at site level
- CLI: `python QuIBL.py inputfile.txt` (config file with treefile path + parameters) -- topology weighting for ILS vs introgression

## Algorithmic Taxonomy

| Tool / Statistic | Approach | Output | Strength | Fails when |
|------------------|----------|--------|----------|------------|
| Patterson's D / ABBA-BABA (Green 2010 Science 328:710; Durand 2011 MBE 28:2239) | Counts ABBA vs BABA site patterns in fixed tree (((P1, P2), P3), Out) | D = (ABBA - BABA) / (ABBA + BABA) | Standard introgression test; tractable | D > 0 means EITHER admixture OR ancestral structure OR ghost lineage |
| Dsuite Dtrios (Malinsky 2021 Mol Ecol Res 21:584) | Fast D for all population trios | D + jackknife SE + p-value | Scales to many populations; integrated with phylogenetic tree | Inherits ABBA-BABA assumptions; symmetric tree only |
| f4-ratio statistic (Patterson 2012 Genetics 192:1065) | Quartet-based admixture proportion | Admixture proportion alpha | Cleaner than D for inferring admixture amount | Requires correct phylogeny; ILS confound |
| Dsuite Fbranch (Malinsky 2018 Nat Eco Evo 2:1940) | Tree-aware f4-ratio mapping | Branch-specific admixture signal | Best at distinguishing admixture among related lineages | Tree must be reliable; ghost branches not detectable |
| ABBAclustering (Koppetsch-Malinsky-Matschiner 2024 Syst Biol) | Dsuite option for divergent-species gene flow | Cluster-based admixture | Designed for cases where standard ABBA-BABA fails | Newer; specific use case |
| TreeMix (Pickrell & Pritchard 2012 PLoS Genet 8:e1002967) | ML tree with migration edges from drift covariance | Tree + migration edges with weights | Visualizes complex admixture | Migration-edge selection subjective; allele-frequency-based |
| qpAdm (Patterson 2012; AdmixTools v2.0; Maier 2023 eLife 12:e85492) | Tests if target derives from sources via rotation | Pass/fail per source set + admixture proportion | Robust statistical framework; rotation control | Source population sampling critical |
| qpGraph (Patterson 2012) | ML phylogenetic graph with admixture | Best-fit graph with admixture nodes | Quantitative; explicit tree+admixture model | Computationally heavy; manual graph topology design |
| HyDe (Blischak 2018 Syst Biol 67:821) | Site-level hybridization detection | Per-site / per-locus hybrid evidence | Sensitive to recent hybridization | Less specific to ancient admixture |
| QuIBL (Edelman 2019 Science 366:594) | Topology weighting on gene trees | Per-tree-weight relative to alternative | Distinguishes ILS from introgression at locus level | Per-locus inference; cross-validation needed |
| Twisst (Martin & Van Belleghem 2017 Genetics 206:429) | Topology weighting on phylogenetic trees | Tree topology weights per genomic window | Visualize topology variation across genome | Computational cost; tree-window decisions |
| PhyloNetworks (Solis-Lemus & Ane 2017 PLoS Comp Biol 13:e1005485; Julia, JuliaPhylo/PhyloNetworks) | SNaQ pseudo-likelihood network inference from gene trees / concordance factors | Explicit reticulation network | Modern Julia ecosystem; SNaQ scales well | Different software than the Java "PhyloNet" |
| PhyloNet (Than-Ruths-Nakhleh 2008 BMC Bioinf 9:322; Java, NakhlehLab/PhyloNet) | ML / MP / Bayesian network inference | Explicit reticulation network | Mature Java tool; many inference modes | Computationally heavy; Maven build |
| sprime (Browning 2018 Cell 173:53) | Per-individual archaic introgression detection | Introgressed haplotype tracts | Designed for archaic-human-like cases | Specific to closely related introgression source |
| Relate (Speidel 2019 Nat Genet 51:1321) + hapne (Fournier 2023 Nat Commun 14:7517) | Phasing-aware genealogy + introgression haplotype dating | Introgressed segments + ages from tract-length distribution | Modern haplotype-aware methods | Require phased data; Relate genealogy precedes hapne dating |
| F3 statistic (Patterson 2012) | Three-population test for admixture | F3 with SE | Specifically detects admixture (vs negative drift) | F3 < 0 indicates admixture; mixed signals |
| F4-statistic (Patterson 2012) | Four-population test | F4 + p-value | Generalizes D; allows variable outgroup distance | Symmetric assumption |

Methodology evolves; verify the Dsuite documentation and Malinsky 2024 review (eLife) before locking on a single approach. The combination of (1) Dsuite D + Fbranch + (2) Twisst / QuIBL + (3) network method (PhyloNet) is the modern best practice for publication-grade introgression claims.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Test for introgression between two species, outgroup available | Dsuite Dtrios + Fbranch | Standard ABBA-BABA + tree-aware branch mapping |
| Multi-population complex admixture inference | TreeMix + qpAdm | Migration edges + rotation tests |
| Distinguish introgression from ILS | QuIBL + Twisst across many loci | Topology weighting at locus level |
| Date the admixture event | Relate (Speidel 2019) + hapne (Fournier 2023) + Twisst window analysis | Genealogy reconstruction + haplotype-length-distribution dating + topology |
| Quantify admixture proportion alpha | f4-ratio (Patterson 2012) or qpAdm | Standard framework |
| Identify introgressed genomic regions | sprime + Twisst | Per-individual + per-window |
| Explicit phylogenetic network | PhyloNet with model selection (DLRS-NL or InferNetwork_MP) | Quantitative reticulation |
| Test for ghost-lineage admixture | qpGraph manual topology + AdmixTools rotation | Compare alternative graphs |
| Recent hybridization | HyDe + per-individual analysis | Site-level hybrid evidence |
| Ancient introgression | Dsuite + qpAdm; Patterson f-statistics for archaic | Tree-aware methods preferred |
| Plant-plant hybridization with parental polyploidy | HyDe + sprime + chromosome-level analysis | Polyploid context |
| Human-archaic introgression (Neanderthal, Denisovan) | sprime + qpAdm + ChromoPainter | Sample-specific haplotype analysis |
| Cichlid radiation | Dsuite + TreeMix + QuIBL (Malinsky 2018 Nat Eco Evo 2:1940) | Established workflow |
| Drosophila species group | Dsuite + Twisst + qpAdm | Standard workflow |
| Genomic regions of introgression vs vertical inheritance | Twisst per-window + chromosome painting | Visualization of variation across genome |
| Test for introgression direction (P1 -> P2 vs P2 -> P1) | qpAdm rotation; Twisst directionality from phasing | qpAdm handles direction better than D |
| Symmetric phylogeny violation (no clear outgroup) | f3 statistic + qpGraph topology | Avoid D-statistic if outgroup unclear |
| Suspect ABBA-BABA ILS-confounded | Switch to QuIBL or Twisst for locus-level discrimination | These distinguish ILS from introgression |

## Per-Method Failure Modes

### ILS confounded with introgression in D-statistic

**Trigger:** D-statistic significantly different from zero for a population trio.

**Mechanism:** Under symmetric phylogeny with no gene flow, ABBA and BABA patterns occur with equal frequency from ancestral polymorphism (ILS). Random allele sorting in descendant lineages produces the same site patterns as introgression. D > 0 is therefore consistent with introgression OR with asymmetric ILS OR with ancestral structure (Eriksson & Manica 2012 PNAS 109:13956).

**Symptom:** Many tested trios show D > 0 with low effect size; Twisst / QuIBL on the same loci shows mixed topologies; phylogenetic distance to outgroup is large.

**Fix:** Use Dsuite f-branch (Malinsky 2018) to assign admixture to specific branches; verify with QuIBL or Twisst topology weighting at the locus level. ABBAclustering (Koppetsch-Malinsky-Matschiner 2024) specifically handles divergent-species gene flow.

### Ghost lineage admixture mimicking signal

**Trigger:** Inferring P1 -> P2 admixture without considering that the admixed lineage could come from an unsampled extinct sister to P3.

**Mechanism:** D > 0 with P3 as donor is identical to D > 0 with a ghost sister of P3 as donor. The data cannot distinguish without explicit modeling.

**Symptom:** Hypothesized donor population doesn't match expected geography or ecology; qpGraph cannot fit a parsimonious admixture model without invoking ghost lineages.

**Fix:** Build qpGraph with potential ghost lineages explicitly; report multiple compatible scenarios; consult known paleontology / biogeography to constrain possibilities. The Neanderthal-Denisovan-anatomically-modern-human framework is heavily ghost-lineage-aware (Slon 2018; Mafessoni 2020).

### Ancestral structure (population subdivision before admixture)

**Trigger:** Significant D-statistic in a recently diverged radiation.

**Mechanism:** Ancestral subdivision in the ancestor of (P1, P2, P3) produces D > 0 from biased ancestral allele inheritance, not from gene flow between descendants.

**Symptom:** D positive across multiple trios but pattern not consistent with simple admixture; effect size differs from f-statistic predictions; rapid radiation context.

**Fix:** Examine F-statistics framework (f3, f4); use qpGraph with structured ancestry; consider PhyloNet for explicit reticulate model. Eriksson & Manica 2012 PNAS 109:13956 documents this confound.

### Outgroup-distance effect on D-statistic

**Trigger:** Using a distantly related outgroup (>200 Myr) for the D-test.

**Mechanism:** Long branches to the outgroup accumulate convergent substitutions; some sites called ABBA / BABA are actually convergent across branches, not true ABBA / BABA patterns.

**Symptom:** D-statistic at convergent-evolution-prone sites is elevated; D varies with outgroup choice.

**Fix:** Use a closer outgroup (< 100 Myr); avoid extremely deep outgroups. Mask sites with extreme conservation across the four taxa. Standard practice in animal species groups uses Drosophila simulans / yakuba for D. melanogaster trios.

### Sample-size bias in D-statistic

**Trigger:** P1 and P2 have very different sample sizes (e.g. 1 vs 50).

**Mechanism:** Allele frequencies in small samples are more extreme; counting ABBA/BABA across genome is dominated by the population with larger sample.

**Symptom:** D-statistic dominated by the larger-sample population; jackknife SE underestimated.

**Fix:** Use one individual per population (Dsuite default for Dtrios); or downsample to consistent N; report effective per-population sample size.

### TreeMix migration-edge selection

**Trigger:** Choosing number of migration edges via likelihood ratio or AIC alone.

**Mechanism:** TreeMix likelihood increases with each migration edge added (more parameters); the "right" number requires statistical justification, but published criteria vary.

**Symptom:** Different selection criteria yield different optimal migration edge counts; results sensitive to k-block size (in jackknife).

**Fix:** Use OptM (Fitak 2021 Biol Methods Protoc 6:bpab017) for principled migration-edge selection; report cross-validation across k-block sizes. Combine with qpAdm rotation tests as independent corroboration.

### qpAdm rotation failures

**Trigger:** qpAdm reporting "rotation rejected" for proposed source population set.

**Mechanism:** qpAdm tests whether target population can be modeled as admixture of source populations; rotation tests check that swapping sources for the outgroup set doesn't break the model. Rotation failure means the proposed sources don't capture the actual admixture.

**Symptom:** qpAdm fits with feasible-sounding sources, but rotation indicates poor fit; admixture proportions sum to weird values.

**Fix:** Reconsider source populations; include "right" populations representing ancestral structures; document rotation results. AdmixTools v2.0 (Maier 2023) provides better diagnostics.

### Multiple-testing correction across population trios

**Trigger:** Running Dsuite Dtrios across all population trios genome-wide and reporting significant ones.

**Mechanism:** With N populations, there are N choose 4 quartets to test; with 50 populations, ~1.5M tests. Naive 5% Type-I gives 75,000 false positives.

**Symptom:** Many "significant" trios reported but with low effect sizes; biological interpretation impossible.

**Fix:** Apply FDR correction across trios; or restrict to a priori hypothesized trios. Dsuite computes Bonferroni-corrected p-values via `-c` flag.

### Genomic-window choice in Twisst / QuIBL

**Trigger:** Running Twisst on 50 kb windows vs 100 kb windows giving different results.

**Mechanism:** Smaller windows have more topology stochasticity; larger windows average across recombination events. The "correct" window depends on local recombination rate.

**Symptom:** Twisst conclusions reverse with window size; spatial pattern of introgression unstable.

**Fix:** Choose windows matching local LD decay (typically 10-100 kb for vertebrates; 5-20 kb for Drosophila). Report sensitivity to window choice. Combine with Twisst-explore.py for spatial visualization.

### Phylogenetic network not unique

**Trigger:** PhyloNet inferring different reticulation networks with similar likelihood.

**Mechanism:** Phylogenetic networks are non-unique; multiple networks can have similar likelihoods. PhyloNet's MP / ML criterion doesn't disambiguate beyond a tolerance.

**Symptom:** PhyloNet returns several networks within delta-LRT = 5; biological interpretation differs.

**Fix:** Report top-N networks; combine with f-statistics (qpAdm / qpGraph) as independent confirmation; use biological knowledge (geography, ecology) to constrain.

### sprime contamination from non-archaic introgression source

**Trigger:** sprime applied to populations with both ancient and recent admixture sources.

**Mechanism:** sprime detects archaic-like haplotypes by comparing to outgroup; haplotype tracts from any non-outgroup-resembling source produce false positives.

**Symptom:** sprime calls many haplotype tracts in populations with known recent introgression that's not from the "archaic" target.

**Fix:** Use sprime only when sources are well-characterized; cross-validate with paired analysis using known archaic vs known recent sources separately.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Patterson's D significance | jackknife z-score |z| > 3 (Bonferroni-corrected) | Green 2010; standard |
| D-statistic interpretation | |D| > 0 with sig z is admixture OR ILS OR ghost | Green 2010 caveat |
| f4-ratio admixture proportion | 0 <= alpha <= 1 | Patterson 2012 |
| Dsuite Fbranch significance | jackknife z > 3 per branch | Malinsky 2018 |
| TreeMix migration weight | typically 0.05-0.30 for biologically relevant | Pickrell 2012 |
| TreeMix k-block size | 500-1000 SNPs default | Empirical |
| qpAdm sources passing rotation | rotation-corrected p > 0.05 | Patterson 2012; Maier 2023 |
| qpAdm admixture proportion | 0 <= alpha <= 1 with SE | Standard |
| QuIBL trio test | best-trio topology weight > 0.5 + significant z | Edelman 2019 |
| Twisst window size | 10-100 kb (vertebrate); 5-20 kb (Drosophila) | Empirical |
| Per-window Twisst topology weight | dominant topology > 50% | Standard |
| HyDe per-individual significance | LRT p < 0.05 vs no-hybridization | Blischak 2018 |
| sprime archaic haplotype length | > 50 kb typical archaic; < 20 kb questionable | Browning 2018 |
| Sample size per population (Dsuite) | 1+ (defaults to one per pop); 5+ preferred | Empirical |
| Outgroup distance for D-stat | < 100 Myr (close); 50-200 Myr (typical) | Empirical |
| FDR threshold across trios | q < 0.05 (BH) | Standard |
| Bonferroni for N choose 4 trios | adjust per number tested | Standard |
| f3 statistic for admixture | f3 < 0 indicates admixture | Patterson 2012 |
| ABBAclustering threshold | depends on clustering algorithm | Koppetsch 2024 |
| Phylogenetic network tolerance | delta-LRT < 2-5 for "similar" | Standard |

## Dsuite Standard Workflow

**Goal:** Compute Patterson's D and Fbranch across all population trios from a VCF.

**Approach:** Prepare SETS file mapping samples to populations -> Dsuite Dtrios -> Dsuite Fbranch -> visualize.

```bash
# 1. Prepare SETS file (samples to populations)
cat > SETS.tsv << 'EOF'
sample_id    population
A1           Population_A
A2           Population_A
B1           Population_B
B2           Population_B
C1           Population_C
D1           Outgroup
EOF

# 2. Run Dsuite Dtrios
# Real CLI flags (verify with `Dsuite Dtrios --help` against installed version):
#   -t/--tree=FILE      species tree (Newick)
#   -o/--out-prefix=    output prefix
#   -k/--no-f4-ratio    skip f4-ratio computation
#   -c/--no-combine     do not write the _combine.txt file
# Fbranch is its own subcommand (`Dsuite Fbranch`); no Dtrios flag toggles it.
Dsuite Dtrios \
    --tree=species_tree.nwk \
    -o trios_run \
    population_genotypes.vcf.gz \
    SETS.tsv

# Output: trios_run_BBAA.txt, trios_run_Dmin.txt, trios_run_tree.txt, trios_run_combine.txt

# 3. Run Fbranch for tree-aware admixture mapping
Dsuite Fbranch species_tree.nwk trios_run_tree.txt > trios_run_fbranch.txt

# 4. ABBAclustering test (Koppetsch-Malinsky-Matschiner 2024).
# Implementation specifics depend on the Dsuite branch / fork distributing the ABBAclustering option;
# consult the Koppetsch 2024 supplement and `Dsuite --help` for current invocation. The option may
# be exposed as a subcommand or per-trio flag rather than a global Dtrios switch.
```

```python
'''Parse Dsuite output for per-trio admixture signal.'''
import pandas as pd


def load_dsuite_trios(path):
    '''SETS.tsv_BBAA.txt columns: P1, P2, P3, D, p-value, Z, BBAA, ABBA, BABA, f_d, f_dM, df, ...'''
    df = pd.read_csv(path, sep='\t')
    df['admixture_signal'] = (df['Z'] > 3) & (df['p_value'] < 0.05)
    return df


def filter_top_admixture(df, by='Z', top=20):
    return df.sort_values(by, ascending=False).head(top)
```

## TreeMix for Population Tree with Migration Edges

**Goal:** Build a population tree with migration edges showing admixture.

**Approach:** Convert VCF to allele-frequency matrix -> run TreeMix -> select migration count via OptM.

```bash
# 1. Convert VCF to TreeMix input
python tree_mix_input.py --vcf population.vcf.gz --popmap popmap.tsv \
    --output input.frq.gz

# 2. Run TreeMix
treemix -i input.frq.gz -m 0 -bootstrap -k 1000 -o output_m0
for m in 1 2 3 4 5; do
    treemix -i input.frq.gz -m $m -bootstrap -k 1000 -o output_m${m}
done

# 3. Select optimal m via OptM (R)
Rscript -e "
library(OptM)
# optM(folder, method, ...): folder is the TreeMix output directory (not the input frq.gz).
# Accepts method='Evanno' (delta-m), 'linear', or 'SiZer'. See ?OptM::optM for details.
optM_result <- optM(folder='.', method='Evanno', tsv=NULL)
"
```

## QuIBL for Locus-Level ILS vs Introgression

**Goal:** Distinguish ILS from introgression at the locus level via topology weighting.

**Approach:** Build per-locus phylogenetic trees -> QuIBL on tree list.

```bash
# Build per-locus gene trees
for region in regions/*.fa; do
    iqtree2 -s $region -m GTR+G -B 1000 -nt 2 --prefix gene_trees/$(basename $region .fa)
done

# Combine into tree list
cat gene_trees/*.treefile > combined_trees.txt

# Run QuIBL: input is a config file pointing at the tree list + parameters
cat > quibl_input.txt << 'EOF'
treefile: combined_trees.txt
outputfile: quibl_results.txt
likelihoodthresh: 0.01
gradascentscalar: 0.5
totaloutgroup: outgroup_taxon
overallnumlambda: 2
EOF
python QuIBL.py quibl_input.txt
```

## qpAdm Rotation Test

**Goal:** Test if target population is admixed from candidate sources.

**Approach:** Define source set + outgroups; AdmixTools qpAdm rotation tests source population set.

```bash
# AdmixTools formatted input: ind / snp / geno (EIGENSTRAT)
# Or use admixtools R wrapper

# In R:
Rscript -e "
library(admixtools)
# Setup left (sources) + right (outgroups)
left <- c('SourceA', 'SourceB')
right <- c('Outgroup1', 'Outgroup2', 'Outgroup3')
target <- 'AdmixedPop'

results <- qpadm(prefix='eigenstrat_data', target=target,
                 left=left, right=right)
results$rankdrop  # rotation test rejection
results$weights   # admixture proportions
"
```

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| D-statistic significant, Twisst dominant topology supports species tree | ILS confound | Trust Twisst; D was ILS-driven |
| D significant, Fbranch attributes to specific branch | True introgression | Confirmed admixture |
| qpAdm fails rotation; TreeMix doesn't show migration edge | Inconsistent admixture scenarios; insufficient data | Reconsider sources; report exploratory |
| QuIBL supports introgression, D-statistic null | D-statistic underpowered or sample-size biased | Trust QuIBL locus-level |
| TreeMix migration edge with TreeMix likelihood gain only with high m | Overfitting | Use OptM; restrict to robust edges |
| PhyloNet network is non-unique | Multiple valid reticulation models | Report all; use prior knowledge to constrain |
| sprime calls many tracts in non-archaic-source population | False positive from non-archaic | Restrict sprime to populations with known archaic source |
| HyDe site-level signal at conserved genes | Convergent evolution masquerading as hybridization | Filter conserved sites; verify with non-coding regions |
| D-statistic stable; Fbranch shows nothing | All admixture attributable to a single ghost branch | Ghost-lineage candidate |

**Operational rule for publication:** Patterson D + Fbranch + Twisst/QuIBL window-level + at least one network method (PhyloNet or qpGraph) agreeing on admixture; multiple-testing correction across trios; ILS quantified via Twisst or coalescent simulation; explicit acknowledgement of ghost-lineage / ancestral-structure alternatives.

## Cohort Gotchas

- **Recent radiations:** ILS confounded with introgression; require additional evidence beyond D-statistic
- **Hybridization vs introgression:** hybridization is contemporary; introgression is ancestral admixture
- **Plant polyploids:** subgenome assignment first; cross-subgenome "introgression" is typically homeologous
- **Reduced-genome organisms:** low marker density; ABBA-BABA underpowered
- **Sex chromosomes:** non-recombining; restrict introgression tests to autosomal data
- **Distantly related introgression source (>200 Myr):** sample more carefully; convergent substitutions confounding
- **Cytonuclear discordance:** mtDNA introgression may differ from autosomal; report separately
- **Migration vs introgression:** D-statistic detects allele introgression; demographic models needed for migration history
- **Rapid radiation example:** cichlids (Malinsky 2018 Nat Eco Evo 2:1940); standard workflow documented

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "ILS vs introgression?" | Multi-method approach: D + Fbranch + Twisst + QuIBL; ILS-confounded loci identified |
| "Ghost lineage?" | qpGraph models tested with potential ghost lineages; biological context constrains; report alternatives |
| "Ancestral structure?" | qpAdm rotation rules out simple admixture from sampled sources; Eriksson & Manica 2012 alternative explicitly considered |
| "Outgroup choice?" | Outgroup distance reported; sensitivity tested with multiple outgroups |
| "Multiple-testing correction?" | FDR across population trios; or restricted to a priori hypothesized |
| "TreeMix migration count?" | OptM Evanno method; cross-validated across k-block sizes |
| "qpAdm sources?" | Rotation tests passed; multiple alternative source sets evaluated |
| "Twisst window size?" | Matched local LD decay; sensitivity tested |
| "Direction of introgression?" | qpAdm directional + Twisst phasing data |
| "Date of admixture?" | hapne + tract-length-distribution-based dating; CIs reported |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Dsuite SETS file empty | Missing samples in VCF | Verify VCF + SETS sample names |
| TreeMix matrix singular | Insufficient SNPs or perfect collinearity | Increase k-block size; remove duplicates |
| qpDstat AdmixTools fails | EIGENSTRAT format issue | Re-convert with convertf |
| HyDe per-individual no result | Default thresholds too strict | Relax --signal threshold |
| QuIBL slow | Many trees | Parallelize; reduce per-tree size |
| Twisst output empty | Window-tree mismatch | Verify gene trees for each window |
| PhyloNet OOM | Many taxa; network complexity | Reduce taxon set; restrict to candidate reticulations |
| sprime no haplotypes | Phasing or outgroup issue | Verify input phased; check outgroup |
| D-statistic sign confused | Allele orientation issue | Use derived alleles consistently |
| qpAdm convergence failure | Source overlap or rare alleles | Filter MAF > 0.05; check source overlap |

## Tool Installation Notes

```bash
# Dsuite
git clone https://github.com/millanek/Dsuite && cd Dsuite && make
# Or: conda install -c bioconda dsuite

# TreeMix
conda install -c bioconda treemix

# AdmixTools (qpDstat, qpAdm, qpGraph)
git clone https://github.com/DReichLab/AdmixTools && cd AdmixTools && make

# Modern AdmixTools v2 (R)
remotes::install_github('uqrmaie1/admixtools')

# HyDe
git clone https://github.com/pblischak/HyDe && cd HyDe && pip install -e .

# QuIBL
git clone https://github.com/miriamtnzr/QuIBL

# Twisst
git clone https://github.com/simonhmartin/twisst && pip install -e .

# PhyloNet
git clone https://github.com/NakhlehLab/PhyloNet && cd PhyloNet && mvn package

# sprime
conda install -c bioconda sprime
```

For population-genetic analyses, the Dsuite + AdmixTools v2 (R) + TreeMix combination is standard. PhyloNet is heavier and requires Maven build.

## References

- Green RE et al 2010 Science 328:710 (Patterson D / ABBA-BABA, Neanderthal admixture)
- Durand EY et al 2011 MBE 28:2239 (D-statistic theoretical framework)
- Patterson N et al 2012 Genetics 192:1065 (f-statistics framework)
- Lipson M et al 2013 PLoS Genet 9:e1004149 (qpAdm)
- Maier R et al 2023 eLife 12:e85492 (AdmixTools v2 / admixtools R)
- Pickrell JK & Pritchard JK 2012 PLoS Genet 8:e1002967 (TreeMix)
- Malinsky M et al 2021 Mol Ecol Resources 21:584 (Dsuite)
- Malinsky M et al 2018 Nat Eco Evo 2:1940 (f-branch statistic; Lake Malawi cichlid radiation)
- Koppetsch T et al 2024 Syst Biol (ABBAclustering)
- Edelman NB et al 2019 Science 366:594 (QuIBL)
- Martin SH & Van Belleghem SM 2017 Genetics 206:429 (Twisst)
- Blischak PD, Chifman J, Wolfe AD & Kubatko LS 2018 Syst Biol 67:821 (HyDe)
- Solis-Lemus C & Ane C 2017 PLoS Comp Biol 13:e1005485 (PhyloNetworks / SNaQ)
- Than C, Ruths D & Nakhleh L 2008 BMC Bioinf 9:322 (PhyloNet)
- Browning SR et al 2018 Cell 173:53 (sprime archaic introgression)
- Eriksson A & Manica A 2012 PNAS 109:13956 (ancestral structure produces D-stat without admixture); Soraggi S et al 2018 G3 8:551 (D-statistic with low-coverage data)
- Slon V et al 2018 Nature 561:113 (Denisovan-Neanderthal hybrid)
- Mafessoni F et al 2020 PNAS 117:15132 (ghost-lineage modeling)
- Speidel L et al 2019 Nat Genet 51:1321 (Relate genealogy)
- Fournier R et al 2023 Nat Commun 14:7517 (hapne; haplotype-length introgression dating)
- Fitak RR 2021 Biol Methods Protoc 6:bpab017 (OptM)
- Skoglund P et al 2012 Genetics 191:295 (admixture inference review)
- Lawson DJ et al 2018 Nat Comm 9:3258 (haplotype-based admixture)
- Patin E et al 2017 Science 356:543 (sub-Saharan admixture example)
- Frantz LAF et al 2020 PNAS 117:17231 (animal domestication admixture)

## Related Skills

- comparative-genomics/hgt-detection - Distinguish HGT from hybridization in microbes
- comparative-genomics/gene-tree-species-tree-reconciliation - Reconciliation alternative for introgression
- comparative-genomics/synteny-analysis - Per-window topology analysis
- population-genetics/population-structure - PCA / ADMIXTURE precedes introgression testing
- population-genetics/selection-statistics - Selection on introgressed regions
- population-genetics/linkage-disequilibrium - LD haplotype context for sprime
- phylogenetics/modern-tree-inference - Gene-tree inference for QuIBL / Twisst
- phylogenetics/species-trees - Coalescent species tree under ILS
- causal-genomics/heritability-partitioning - Inheritance partitioning by genomic region
- variant-calling/joint-calling - Multi-sample VCF for introgression analysis
- read-alignment/bwa-alignment - Alignment underlies variant calling for D-statistic
