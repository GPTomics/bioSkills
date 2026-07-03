---
name: bio-comparative-genomics-whole-genome-duplication
description: Detect, date, and contextualize whole-genome duplication (WGD / paleopolyploidy) events using wgd v2 (Chen & Zwaenepoel 2024), KsRates (Sensalari 2022 substitution-rate-corrected Ks dating), DupGen_finder (Qiao 2019), MAPS (Li 2018 phylogenomic), POInT (Conant 2008 ordered-block), SLEDGe (2024 ML-based), Whale.jl (Bayesian DTL+WGD), and synteny-anchored paranome construction. Use when identifying ancient polyploidy from Ks distributions and synteny block analysis, positioning WGD events relative to speciation, distinguishing tandem from segmental from WGD duplications, dating the 2R/3R vertebrate / fish / salmonid WGDs, building paranome and Ks-age mixture models, applying KsRates substitution-rate correction across lineages, or testing alternative biased-fractionation / dosage-balance models post-WGD.
tool_type: mixed
primary_tool: wgd
---

## Version Compatibility

Reference examples tested with: wgd v2.0.31+ (heche-psb/wgd; Chen & Zwaenepoel 2024 Bioinformatics 40:btae272), KsRates 1.1.3+ (VIB-PSB/ksrates; Sensalari 2022 Bioinformatics 38:530), DupGen_finder (Qiao 2019 Genome Biol 20:38), MAPS 1.0 (Li 2018), POInT (Conant lab), SLEDGe (bioRxiv 2024.01.17.574559), Whale.jl 2.0+, ksrates pip 1.1+, MCScanX 1.0+, PAML 4.10+ (yn00/codeml for Ks), BLAT 36+, DIAMOND 2.1+, R 4.4+, mclust 6.1+ (for mixture models). Python 3.10+ required for wgd v2.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `wgd --version`; `ksrates --version`; `wgd ksd --help`
- Python: `pip show wgd ksrates`
- R: `packageVersion('mclust')`

If code throws `wgd ksd: cannot find PAML output`, `KsRates: insufficient sister species`, `MAPS: missing tree`, these tools have specific input expectations: wgd needs codon-aware MAFFT/MUSCLE alignment; KsRates needs configured `config_ksrates.txt`; MAPS needs nucleotide tree. The deprecated arzwa/wgd v1 is replaced by heche-psb/wgd v2.

# Whole Genome Duplication Analysis

**"Are there WGD events in this lineage and when did they occur?"** -> WGD detection combines **Ks distributions** (synonymous-substitution rates between gene paralog pairs, showing peaks at past polyploidy events) and **synteny block analysis** (parallel collinear blocks within a genome). Modern best practice uses **wgd v2** (Chen & Zwaenepoel 2024 Bioinformatics 40:btae272) as an integrated pipeline. **KsRates** (Sensalari 2022 Bioinformatics 38:530) is mandatory for cross-lineage comparison because substitution rates vary across the tree -- ignoring this places WGDs incorrectly relative to speciation events. The fundamental tradeoff: Ks plot peaks are visually obvious but biologically ambiguous between (1) small-scale tandem duplications, (2) segmental duplications, and (3) true WGD; combining Ks with synteny anchors disambiguates these.

- CLI: `wgd dmd` and `wgd ksd` -- paranome construction + Ks distribution
- CLI: `wgd syn` -- synteny-anchored WGD signal extraction
- CLI: `ksrates init` then `ksrates wgd-paralogs ortho` -- substitution-rate-corrected positioning
- CLI: `MCScanX -h` then `dupgen_finder` -- duplication-class assignment
- CLI: `mapsR` -- gene-tree-based WGD phylogenetic placement

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| wgd v2 (Chen & Zwaenepoel 2024 Bioinformatics 40:btae272) | Integrated paranome + Ks + synteny + WGD dating | Ks distributions, collinearity plots, GMM/ELMM mixture fits, dating | Standard 2024 pipeline; replaces deprecated arzwa/wgd v1 | Saturation at Ks > 2; single-lineage substitution rate variation |
| KsRates (Sensalari 2022 Bioinformatics 38:530) | Substitution-rate correction via outgroup pairs | Adjusted Ks ages of focal-species paralogs vs orthologs | MANDATORY when comparing WGDs across lineages with different rates | Requires at least 2 outgroups for rate calibration |
| DupGen_finder (Qiao 2019 Genome Biol 20:38) | Classifies duplications by genomic context | Per-gene class: tandem, proximal, dispersed, segmental, WGD | Disambiguates duplication type | Class assignment depends on intervening-gene-count windows |
| MAPS (Li 2018) | Phylogenomic placement of WGD via gene-tree topology mapping | WGD position on species tree | Detects WGDs from gene-tree-species-tree discordance | Computationally heavy; requires many gene trees |
| POInT (Conant & Wolfe 2008 Genetics 179:1681) | Order-aware reconstruction of WGD chromosomes | Reconstructed ancestral WGD genome | Strong inference for syntenic-block ages | Lineage-specific tuning required |
| SLEDGe (bioRxiv 2024.01.17.574559) | ML classifier on Ks plot features | WGD vs no-WGD binary call + confidence | Reduces visual-peak-fitting subjectivity | Newer; less validated |
| wgd v1 (arzwa/wgd; DEPRECATED) | Predecessor of v2 | -- | Historical | Use v2 (heche-psb/wgd) |
| dupHMM (Zwaenepoel 2019 different MBE paper) | HMM on Ks for WGD/SSD attribution | Per-pair WGD/SSD classification | Probabilistic disambiguation | Less integrated than wgd v2 |
| Whale.jl (Zwaenepoel & Van de Peer 2019 MBE 36:1384) | Bayesian DTL + WGD reconciliation | WGD posterior at species-tree nodes | Native WGD modeling; integrates with [[gene-tree-species-tree-reconciliation]] | Julia ecosystem |
| WGDexploreR (legacy) | Visualize Ks plots | Plots only | Visualization aid | Not for inference |
| McLuster-WGD (custom workflows) | Mixture-model fitting on Ks | GMM / ELMM components | For custom peak fitting | Not a standard tool |
| ksrates web (sensalari 2022) | Web interface to ksrates | Same as CLI | User-friendly | Manual config; not scriptable for genome-wide |
| wgs2pep (Vandepoele/wgs) | Pep-to-WGD ortholog identification | Per-paralog WGD assignment | Plant-focused | Less general |

Methodology evolves; verify the current wgd v2 manual and the Chen & Zwaenepoel 2023 review chapter (in *Polyploidy: Methods and Protocols*) before locking on a single workflow. The 2R vertebrate / 3R teleost / Ss4R salmonid WGDs are well-established; novel WGD claims require concordance across Ks, synteny, and phylogenomic placement.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Plant comparative genomics; suspect WGD | wgd v2 + KsRates | Modern standard pipeline; plants are WGD-prone |
| Vertebrate 2R or 3R WGD analysis | wgd v2 + MAPS phylogenomic placement | Vertebrate-specific; combines Ks + tree topology |
| Salmonid Ss4R WGD | Whale.jl with WGD node in species tree | Native WGD modeling in Bayesian framework |
| Distinguish WGD from sequential small duplications | Ks GMM fit + synteny block analysis (wgd syn) | Ks peak alone insufficient; synteny confirms |
| Date WGD relative to a speciation event | KsRates with outgroup speciation calibration | Required when substitution rates differ |
| Identify duplicates retained vs lost after WGD | DupGen_finder + Ks distribution | Classification + age estimation |
| Detect novel WGD in non-model organism | wgd dmd -> wgd ksd -> visual GMM/ELMM mixture | Standard discovery workflow |
| Reconstruct ancestral WGD genome architecture | POInT | Order-aware ancestral reconstruction |
| Test post-WGD bias in fractionation (dosage-balance) | DupGen_finder + per-gene functional annotation | Compare retained vs lost genes |
| Distinguish recent vs ancient WGD | KsRates + Ks peak location after correction | Ks peak < 0.1 recent; 0.5-1.5 ancient |
| WGD in clade with rapid evolution (e.g. fish) | KsRates mandatory | Without rate correction, Ks peaks displaced |
| ML classification of WGD signal | SLEDGe | Newer ML-based alternative to manual fitting |
| Integrate WGD with DTL inference | Whale.jl OR ALE with WGD branch | See [[gene-tree-species-tree-reconciliation]] |
| Pan-clade WGD survey (e.g. all angiosperms) | wgd v2 batch + ksrates aggregated | Standard for sustained phylogenomic surveys |
| Polyploid genome with subgenomes | DupGen_finder for tandem/segmental + AnchorWave proali for subgenome-aware synteny | Subgenome assignment first |
| Recently diverged species pair, possible recent WGD | minimap2 -x asm5 + SyRI + Ks distribution | High-resolution recent-WGD detection |

## Per-Tool Failure Modes

### Saturation at Ks > 2

**Trigger:** Computing Ks for ancient WGD candidates from distantly related taxa.

**Mechanism:** Synonymous substitutions saturate after Ks ~2; each site has undergone multiple substitutions. Observed Ks underestimates true Ks; the relationship between Ks and time becomes non-monotonic above 1.5. WGD peaks at age 200+ Myr in vertebrates / 100+ Myr in plants are at Ks > 2 and unreliable (Vanneste 2013 MBE 30:177).

**Symptom:** wgd ksd output shows peak at Ks > 1.5 with broad distribution; KsRates rate-corrected Ks even more uncertain; visual fitting yields ambiguous components.

**Fix:** Restrict Ks-based WGD inference to Ks < 1.5; for older WGDs, use phylogenomic methods (MAPS, POInT, Whale.jl) which don't rely on Ks alone. Report Ks-based dating with explicit saturation caveat. For 2R vertebrate WGD (~500 Myr), MAPS / Whale.jl is required because Ks is saturated.

### Substitution-rate variation across lineages

**Trigger:** Comparing WGD position in two lineages with different molecular evolutionary rates.

**Mechanism:** If lineage A evolves twice as fast as lineage B, the same biological time is at twice the Ks in lineage A. A WGD at Ks = 0.5 in lineage A and Ks = 0.25 in lineage B might be the same event biologically.

**Symptom:** WGD claimed at "different times" in different lineages; speciation-vs-WGD relative timing flips depending on which lineage is the focal species.

**Fix:** Always use KsRates (Sensalari 2022) when comparing WGDs across lineages. KsRates uses outgroup-species speciation events to calibrate the relative substitution rates, then rescales Ks to a common scale. Single-lineage Ks dating is unreliable for inter-lineage comparison.

### Tandem duplication masquerading as WGD peak

**Trigger:** Ks distribution shows a peak at low Ks (< 0.5); user concludes "recent WGD."

**Mechanism:** Tandem duplications create paralog pairs with low Ks; these accumulate in clusters (especially in NLR genes in plants, OR genes in mammals). The Ks peak reflects clustered tandem origins, not WGD.

**Symptom:** "WGD" peak Ks-distribution is dominated by paralogs from tandem clusters; synteny analysis shows few cross-chromosome parallel blocks.

**Fix:** Use DupGen_finder to classify each paralog pair as tandem / proximal / dispersed / segmental / WGD. The true WGD signal comes from segmental + WGD classes. Recompute Ks distribution restricted to WGD-classified pairs (or anchor-pair restricted to synteny blocks). wgd v2 integrates this filtering.

### Synteny block age inconsistency

**Trigger:** Different synteny blocks from the same suspected WGD have different Ks distributions.

**Mechanism:** Real WGD blocks should share a synchronized Ks distribution centered at the WGD age; if blocks have very different ages, the inferred "WGD" was probably a series of segmental duplications.

**Symptom:** Per-block Ks median varies widely; some blocks show Ks ~ 0.3 and others Ks ~ 1.0; not consistent with single WGD.

**Fix:** Compute per-block Ks distribution; require synchronization (interquartile range overlapping across blocks). wgd v2 reports per-block age dispersion. If dispersion is high, downgrade to "potentially WGD-like" or attribute to ancient segmental duplication history.

### Mixture model under/overfitting components

**Trigger:** Fitting GMM (gaussian mixture model) or ELMM (exponential-lognormal mixture model) to Ks distribution.

**Mechanism:** Mixture models with too few components miss real WGD signal; too many components find spurious peaks. BIC-based model selection is standard but sensitive to bin choice in histograms.

**Symptom:** Different runs with different bin counts give different component numbers; visual peaks don't match component means.

**Fix:** Use BIC-based component selection with 5-fold cross-validation; report uncertainty in component count. wgd v2 fits both GMM and ELMM and reports BIC for each model up to 5 components. Visually validate against synteny block analysis.

### Comparing wgd v1 (deprecated) vs v2 outputs

**Trigger:** Mixing scripts written for arzwa/wgd v1 with heche-psb/wgd v2.

**Mechanism:** wgd v1 was deprecated 2023; v2 has different command-line interface, default parameters, and output file structure. Scripts from v1 fail or produce different results in v2.

**Symptom:** Older scripts using `wgd ksd` v1 syntax don't work in v2; output paths differ.

**Fix:** Update all scripts to wgd v2 syntax. The v1 wgd is no longer maintained; v2 is the current standard.

### KsRates failing with insufficient species

**Trigger:** Running KsRates with one focal species and one outgroup.

**Mechanism:** KsRates calibrates substitution-rate variation via at least two outgroup speciation events; one outgroup species provides only one rate point and cannot resolve rate heterogeneity.

**Symptom:** KsRates errors or produces unreliable rate-corrected Ks.

**Fix:** Include >= 2 outgroup species at different phylogenetic distances; for plants, a model angiosperm + gymnosperm pair; for vertebrates, lamprey + invertebrate outgroup. Document outgroup sampling.

### Polyploid genome confusion

**Trigger:** Running wgd v2 on an unassigned-subgenome polyploid (e.g. wheat hexaploid without A/B/D subgenome labels).

**Mechanism:** WGD signal is muddled when subgenomes are not separated; orthologs across subgenomes appear as paralogs, intermixing tandem / segmental / WGD / homeolog signal.

**Symptom:** Ks distribution shows multiple overlapping peaks at varying intensities; component fits are unstable.

**Fix:** Assign subgenomes first using k-mer methods (KMC2 + SubPhaser; Jia 2022 New Phytol 235:801) or synteny-based assignment (GENESPACE; Lovell 2022). Then run wgd v2 on each subgenome separately. Document subgenome assignment.

### Confusion of homeologous (WGD) vs orthologous (speciation) pairs

**Trigger:** Computing Ks across recent WGD species without explicit homeolog labeling.

**Mechanism:** In recent WGD lineages (e.g. salmonid Ss4R), the "ortholog" between species and the "homeolog within species" can have similar Ks; standard ortholog detection (OrthoFinder) lumps both.

**Symptom:** wgd v2 reports Ks distributions but doesn't distinguish homeologs from orthologs; downstream dating confused.

**Fix:** Use synteny-anchored ortholog detection (GENESPACE, ProteinOrtho-synteny) to separate cross-species orthologs from within-species homeologs. Run wgd v2 separately on each class. Report each Ks distribution explicitly.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Ks saturation upper limit | Ks < 1.5 for reliable inference; >= 2 saturated | Vanneste 2013 MBE 30:177 |
| Recent WGD Ks range | Ks 0.1-0.5 | Standard convention |
| Ancient WGD Ks range | Ks 0.5-1.5 | Standard convention |
| WGD synteny block minimum | >= 5 anchors per block | wgd v2 / GENESPACE default |
| Tandem duplicate window for DupGen_finder | 5 genes default; species-tunable | Qiao 2019 |
| Proximal duplicate window | 5-25 genes | Standard |
| Number of mixture components to consider | 1-5 with BIC selection | wgd v2 default |
| Cross-block age synchrony for WGD | per-block Ks IQR overlap | Visual + statistical |
| KsRates minimum outgroups | >= 2 species at different distances | Sensalari 2022 |
| KsRates substitution rate correction valid range | Ks < 1.5 in original; correction can extend slightly | KsRates docs |
| MAPS minimum gene trees | >= 1000 single-copy ortholog trees | Li 2018 |
| Whale.jl WGD detection power | depends on retention rate; > 0.3 retention typically detectable | Zwaenepoel 2019 |
| Vertebrate 2R Ks | saturated; estimated 500-700 Myr | Dehal 2005 |
| Teleost 3R Ks | saturated; estimated 250-350 Myr | Glasauer 2014 |
| Salmonid Ss4R Ks | 80-100 Myr; Ks ~0.1-0.2 | Lien 2016 |
| Plant 1R / 2R | varies clade; angiosperm zeta ~120 Myr | Soltis 2009 |
| Synonymous codon site count | dS reliable when >= 30 synonymous sites per pair | Yang 2007 PAML |
| Per-pair gene length minimum | >= 300 bp CDS | wgd v2 default |
| Codon-aware MSA aligner | MAFFT --auto or MUSCLE; for distant pairs PRANK | wgd v2 supports all |
| Bootstrap support for MAPS | >= 80 on gene-tree branches | Li 2018 |

## wgd v2 Standard Workflow

**Goal:** Detect and date WGD events from a focal-species proteome and CDS, with synteny-anchored Ks distribution.

**Approach:** Build paranome -> compute Ks distribution -> identify synteny anchors -> mixture model fit -> visualize.

```bash
# Install (Python 3.10+ required)
pip install wgd

# 1. Build paranome (all-vs-all paralog identification)
wgd dmd cds.fasta -o output/paranome.tsv -t 16

# 2. Compute Ks distribution
# Verify exact flags with `wgd ksd --help` (the wgd v2 CLI evolves; --pairwise / --ks-method spelling differs across versions).
wgd ksd output/paranome.tsv cds.fasta -o output/ksd \
    --aligner mafft --n-threads 16

# 3. Synteny anchors (intra-genome)
wgd syn output/paranome.tsv gff.bed cds.fasta -o output/syn \
    --feature gene --gene-attribute Name --min-block-size 5

# 4. Mixture model on synteny-anchored Ks
wgd mix output/ksd/ks_distributions.tsv -o output/mix \
    --model both --components 1 2 3 4 5

# 5. Visualize
wgd viz output/mix/mix_results.tsv -o output/viz \
    --type histogram --bins 50
```

```python
'''Parse wgd v2 outputs and identify WGD peak from Ks distribution.'''
import pandas as pd
from scipy import stats


def load_ks_distribution(ksd_file):
    '''Load wgd ksd output: pair  Ks  Ka  ...'''
    df = pd.read_csv(ksd_file, sep='\t')
    # Filter saturated
    df_filt = df[(df['Ks'] > 0) & (df['Ks'] < 2)]
    return df_filt


def fit_mixture(ks_values, n_components=3):
    '''Fit GMM (gaussian mixture model) on Ks via sklearn.'''
    from sklearn.mixture import GaussianMixture
    gmm = GaussianMixture(n_components=n_components, random_state=42)
    gmm.fit(ks_values.reshape(-1, 1))
    return {
        'means': gmm.means_.flatten(),
        'variances': gmm.covariances_.flatten(),
        'weights': gmm.weights_,
        'bic': gmm.bic(ks_values.reshape(-1, 1))
    }


def identify_wgd_peaks(ks_values, min_components=1, max_components=5):
    '''BIC-based model selection.'''
    results = {}
    for n in range(min_components, max_components + 1):
        results[n] = fit_mixture(ks_values, n)
    best_n = min(results.keys(), key=lambda k: results[k]['bic'])
    return best_n, results[best_n]
```

## KsRates Substitution-Rate-Corrected Workflow

**Goal:** Position a focal-species WGD relative to speciation events with substitution-rate correction.

**Approach:** Define focal species + 2+ outgroups -> compute orthologous Ks (focal vs outgroup) -> compute paralogous Ks (focal-internal) -> rescale via outgroup calibration.

```bash
# KsRates is best driven through its Snakemake pipeline (which orchestrates ortholog Ks,
# paralog Ks, rate correction, and plotting). Subcommand naming differs across releases;
# always verify with `ksrates --help` against the installed version.

# 1. Generate a working config (subcommand spelling varies by release; e.g. `init` in some,
#    `generate-config` in others). Inspect `ksrates --help`.
ksrates init config_ksrates.txt        # OR: ksrates generate-config config_ksrates.txt

# 2. Edit config_ksrates.txt to set focal_species, outgroups, FASTA + GFF paths, tree.

# 3. Run the full pipeline (Snakemake-driven):
ksrates --config config_ksrates.txt --n-threads 16
# Or invoke individual stages (subcommand names: see `ksrates --help`).
```

The output plot shows the rate-corrected focal-species paralog Ks distribution with vertical lines indicating outgroup-species speciation events; WGD peaks before vs after speciation events can be distinguished.

## DupGen_finder for Duplication Class Assignment

**Goal:** Classify each paralog pair as tandem / proximal / dispersed / segmental / WGD by genomic context.

**Approach:** MCScanX collinearity + intervening-gene count -> classify each duplicate by class.

```bash
# Input: MCScanX collinearity file + GFF
git clone https://github.com/qiao-xin/DupGen_finder
cd DupGen_finder

# Run with intervening-gene count
./DupGen_finder.pl \
    -i input_collinearity_file \
    -t species_name \
    -c gene_count_file \
    -o output

# Output: tandem, proximal, dispersed, segmental, wgd duplicate lists
```

```python
'''Aggregate DupGen_finder output for downstream Ks distribution per class.'''
def load_dupgen(class_dir):
    classes = {}
    for cls in ('tandem', 'proximal', 'dispersed', 'segmental', 'wgd'):
        path = f'{class_dir}/{cls}.pairs'
        try:
            classes[cls] = pd.read_csv(path, sep='\t', header=None,
                                       names=['gene1', 'gene2'])
        except FileNotFoundError:
            classes[cls] = pd.DataFrame(columns=['gene1', 'gene2'])
    return classes
```

## MAPS Phylogenomic WGD Placement

**Goal:** Place a WGD event on a species tree using gene-tree-species-tree mapping.

**Approach:** Build many single-copy ortholog gene trees -> map each tree's topology to the species tree -> identify branches with topology consistent with a WGD event.

MAPS is heavyweight; requires CRG database setup and significant compute. See https://bitbucket.org/barker-lab/maps/src for the current pipeline.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| wgd v2 peak at Ks 0.5, MAPS supports WGD on different branch | Rate variation across lineages | KsRates with outgroups; trust rate-corrected position |
| Ks peak prominent, synteny anchors sparse | Tandem-driven peak | DupGen_finder classification; verify WGD class |
| Synteny blocks present, Ks peak diffuse | Old WGD; saturated Ks | Restrict to recent paralogs; expect Ks 0.5-1.5 |
| One outgroup's KsRates says WGD before speciation, other after | Outgroup substitution rate variation | Use multiple outgroups; report ranges; explicit caveat |
| GMM fits 2 components, ELMM fits 3 | Model class differs | Visual inspection; report both; choose by BIC |
| DupGen_finder calls "WGD" but Ks > 2 | Saturation | Ks unreliable; verify via synteny + phylogenomic placement |
| Two recently sequenced sister species disagree on WGD date | Different annotation pipelines | Re-annotate consistently; verify orthology |
| wgd v2 detects WGD; Whale.jl posterior doesn't | Different inference frameworks | Whale.jl is Bayesian and gene-tree-based; trust if MAPS / phylogeny supports |
| Salmonid Ss4R detectable by wgd v2 in young salmonids | Recent WGD; clearer signal | Confirmation, not contradiction |

**Operational rule for publication:** Concordance across Ks distribution (wgd v2), synteny-anchored block analysis (wgd syn / GENESPACE), KsRates rate-corrected positioning, and at least one phylogenomic method (MAPS / Whale.jl / ALE with WGD node) = publication-grade WGD claim. Single-Ks-peak claims should be downgraded to "Ks evidence for possible WGD."

## Cohort Gotchas

- **Plant WGD legacy:** All angiosperms share at least one ancestral WGD; the zeta WGD ~120 Myr (Soltis 2009 Am J Bot 96:336); confirm consensus against published plant lineages
- **Vertebrate 2R:** ~500-600 Myr; Ks saturated; use MAPS phylogenomic placement (Dehal & Boore 2005)
- **Teleost 3R:** ~320 Myr; in addition to 2R; doubles ohnologs in fish
- **Salmonid Ss4R:** ~80-100 Myr; recent enough that Ks is informative; ohnologs identifiable
- **Plant polyploids (wheat, Brassica):** recent; subgenomes assignable; subgenome-stratify before wgd v2
- **Allopolyploids vs autopolyploids:** allopolyploids have higher inter-subgenome divergence (easier to assign subgenomes); autopolyploids have similar subgenomes (chimeric assembly risk)
- **Lineage-specific extra duplications post-WGD:** make Ks peak broader; differentiate from sequential WGDs

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Saturation?" | Ks < 1.5 reported; for older WGDs, MAPS or Whale.jl used |
| "Rate variation?" | KsRates applied with 2+ outgroups; rate-corrected positioning |
| "Tandem vs WGD?" | DupGen_finder classification; WGD-class subset reported |
| "Synteny confirmation?" | wgd syn or GENESPACE; per-block Ks synchrony reported |
| "Phylogenomic support?" | MAPS, Whale.jl, or ALE-with-WGD agrees |
| "Mixture model uncertainty?" | BIC-based model selection; 1-5 components tested |
| "Subgenome assignment (polyploid)?" | k-mer-based or synteny-based subgenome assignment first |
| "Why wgd v2 over wgd v1?" | v1 deprecated 2023; v2 (heche-psb/wgd) is current |
| "Comparison to published WGD dates?" | Match published vertebrate 2R / teleost 3R / salmonid Ss4R; cite primary literature |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `wgd ksd` reports few Ks values | CDS sequences not in-frame or short | Verify CDS quality; filter short sequences |
| KsRates "no orthologs found" | Outgroup species names mismatched | Normalize species labels exactly |
| GMM fits with 5 components, all similar means | Overfitting | Use BIC; restrict max components to 3-4 |
| Per-block Ks IQR > 0.5 | Block age dispersion; likely segmental not WGD | Trust dispersion; reclassify as segmental |
| wgd v2 ELMM converges to single component | Insufficient pairs | Reduce filtering; include more paralog pairs |
| MAPS slow / OOM | Many gene trees; large dataset | Cluster; reduce to representative single-copy orthologs |
| Whale.jl Turing fails to mix | NUTS step size; tree dimension | HMC manual tuning; or reduce species tree |
| DupGen_finder produces empty tandem set | Default window inappropriate | Adjust window for genome density |
| Ks distribution has secondary peak at 1.8 | Saturation artifact (not WGD) | Recompute with PRANK or PAML codeml; restrict Ks < 1.5 |
| Polyploid "WGD" appears at Ks 0.05 | Recent polyploidy is ohnologous WGD, expected | Verify with subgenome assignment |
| Plant analysis lacks reference outgroup gymnosperm | KsRates needs gymnosperm + angiosperm | Add Selaginella or moss outgroup |

## Tool Installation Notes

```bash
# wgd v2
pip install wgd
# Or: git clone https://github.com/heche-psb/wgd && cd wgd && pip install .

# KsRates
pip install ksrates

# DupGen_finder (Perl + MCScanX dependency)
git clone https://github.com/qiao-xin/DupGen_finder
# Install MCScanX (Wang 2012): git clone https://github.com/wyp1125/MCScanX && cd MCScanX && make

# MAPS (Python; requires CRG, ete3)
git clone https://bitbucket.org/barker-lab/maps

# Whale.jl (Julia)
julia -e 'using Pkg; Pkg.add("Whale")'

# SLEDGe
git clone https://github.com/SLEDGe-team/SLEDGe

# POInT (C++)
git clone https://github.com/gconant0/PoInT
cd PoInT && make

# PAML for Ks computation (yn00 method)
conda install -c bioconda paml

# DIAMOND / BLAST for paranome
conda install -c bioconda diamond blast

# R packages
install.packages(c('mclust', 'rmixmod'))
```

For new analyses, default to wgd v2 + KsRates as the primary pipeline; MAPS or Whale.jl for confirmation on phylogenomic placement.

## References

- Chen H & Zwaenepoel A 2024 Bioinformatics 40:btae272 (wgd v2)
- Sensalari C et al 2022 Bioinformatics 38:530 (KsRates)
- Qiao X et al 2019 Genome Biol 20:38 (DupGen_finder)
- Li Z et al 2018 PNAS 115:4713 (MAPS phylogenomic WGD placement)
- Conant GC & Wolfe KH 2008 Genetics 179:1681 (POInT)
- Zwaenepoel A & Van de Peer Y 2019 MBE 36:1384 (Whale.jl)
- SLEDGe team 2024 bioRxiv 2024.01.17.574559 (SLEDGe ML classifier)
- Vanneste K et al 2013 MBE 30:177 (Ks saturation)
- Holland PWH et al 1994 Development Suppl:125 (2R hypothesis)
- Dehal P & Boore JL 2005 PLoS Biol 3:e314 (vertebrate 2R confirmation)
- Glasauer SMK & Neuhauss SCF 2014 Mol Genet Genomics 289:1045 (teleost 3R)
- Lien S et al 2016 Nature 533:200 (salmonid Ss4R)
- Soltis DE et al 2009 Am J Bot 96:336 (polyploidy and angiosperm diversification)
- Birchler JA & Veitia RA 2007 Plant Cell 19:395 (gene balance hypothesis)
- Force A et al 1999 Genetics 151:1531 (subfunctionalization)
- Maere S et al 2005 PNAS 102:5454 (post-WGD retention bias)
- Tang H et al 2008 GR 18:1944 (synteny / MCScan)
- Lovell JT et al 2022 eLife 11:78526 (GENESPACE)
- Jia K-H et al 2022 New Phytol 235:801 (SubPhaser subgenome assignment)
- Yang Z 2007 PAML manual (yn00 codon Ks)
- Chen H & Zwaenepoel A 2023 Methods Mol Biol 2545:3 (inference of ancient polyploidy from genomic data)

## Related Skills

- comparative-genomics/synteny-analysis - Synteny-anchored Ks (wgd syn / GENESPACE)
- comparative-genomics/ortholog-inference - Ortholog detection feeds DupGen / wgd
- comparative-genomics/gene-tree-species-tree-reconciliation - Whale.jl native WGD modeling; ALE with WGD branch
- comparative-genomics/gene-family-evolution - CAFE5 birth-death often shows post-WGD retention bias
- comparative-genomics/positive-selection - Selection on retained ohnologs (post-WGD diversification)
- comparative-genomics/ancestral-reconstruction - Ancestral pre-WGD gene state
- phylogenetics/divergence-dating - Time-calibrated WGD positioning
- alignment/multiple-alignment - Codon-aware MSA for Ks computation
- alignment/pairwise-alignment - Pairwise Ks via PAML yn00 / codeml
- genome-assembly/assembly-qc - BUSCO / Compleasm before WGD analysis
