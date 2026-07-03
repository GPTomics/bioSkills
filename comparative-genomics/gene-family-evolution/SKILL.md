---
name: bio-comparative-genomics-gene-family-evolution
description: Model gene-family birth-death dynamics across a species tree using CAFE5 (Mendes et al 2020 Bioinformatics 36:5516 gamma-distributed rate categories), CAFE5-error (annotation-error-aware), Count (Csurös 2010 ancestral state reconstruction), BadiRate (Librado 2012 likelihood + parsimony), DupliPHY-Family, and ALE/AleRax (for per-family DTL; see [[gene-tree-species-tree-reconciliation]]). Test lineage-specific gene-family expansions and contractions, distinguish biological dynamics from annotation artifacts, account for assembly fragmentation, identify functional enrichment in expanded / contracted families. Use when correlating gene-family changes with phenotype evolution, ranking lineages by adaptive gene-family-rate shifts, post-WGD dosage-balance analysis, or building Birth-death models from OrthoFinder presence/absence matrices.
tool_type: cli
primary_tool: CAFE5
---

## Version Compatibility

Reference examples tested with: CAFE5 5.1.0+ (Mendes et al 2020 Bioinformatics 36(22-23):5516-5518), Count 11.0319+ (Csurös 2010 Bioinformatics 26:1910), BadiRate 1.35+ (Librado 2012 Bioinformatics 28:279), DupliPHY-Family (Ames et al 2012), CAFExp (legacy CAFE 4.2 -- DEPRECATED; use CAFE5), OrthoFinder 3.0+ for HOG input, R 4.4+, mclust 6.1+, phytools 2.3+, ETE4 4.1.0+ for tree manipulation. ALE/GeneRax/AleRax in companion skill [[gene-tree-species-tree-reconciliation]].

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `cafe5 --help`; `Count.exe` (Java); `badirate --help`
- R: `packageVersion('phytools')`
- Python: `pip show ete4`

If code throws `CAFE5: lambda did not converge`, `Count negative branch length`, `BadiRate gamma not initialized`, the most common causes are: (1) annotation heterogeneity inflating family sizes, (2) saturated families (CAFE5 needs reasonable rate variation), (3) negative branch lengths in input tree (Count requires ultrametric). Pre-process: filter OG matrix to families present in >= 50% of species; resolve polytomies; ultrametricize tree.

# Gene Family Evolution

**"Which gene families expanded or contracted in which lineages?"** -> Birth-death models on phylogeny (Hahn 2007; Csurös 2010) treat each orthogroup's per-species count as evolving under a stochastic birth-death process; lineage-specific rate shifts are detected as departures from a global rate. **Annotation heterogeneity is the single largest confounder**: different annotation pipelines predict different numbers of genes per family, producing apparent lineage-specific expansions that are artifacts of annotation choice (Tonkin-Hill 2020 demonstrated this for bacterial pangenomes). Consistent annotation + BUSCO/Compleasm completeness filtering are mandatory before any birth-death model interpretation. CAFE5 (Mendes et al 2020 Bioinformatics 36:5516) replaces older CAFE versions with gamma-distributed rate categories for more biologically realistic modeling.

- CLI: `cafe5 -i orthogroup_counts.tsv -t species_tree.nwk -p` -- main CAFE5 workflow
- CLI: `cafe5 -e` -- error-aware mode for annotation uncertainty
- CLI: `Count` -- Java GUI / CLI for parsimony + likelihood ASR
- CLI: `badirate -t tree.nwk -d counts.tsv` -- likelihood birth-death + branch parsimony

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| CAFE5 (Mendes et al 2020 Bioinformatics 36(22-23):5516-5518) | Birth-death with gamma rate categories | Global / per-family lambda + significant rate shifts | Modern standard; handles rate heterogeneity; explicit Type-I control | Annotation heterogeneity confounds; needs > 100 families |
| CAFE5-error | Annotation-error-aware extension | Same plus error estimates | Critical for noisy annotations | Manual error-rate specification or estimation |
| Count (Csurös 2010 Bioinformatics 26:1910) | Both ML and parsimony ASR | Branch event counts (D, L) per family | Comprehensive output; GUI | Slower than CAFE5; less modern UX |
| BadiRate (Librado 2012 Bioinformatics 28:279) | Likelihood birth-death + branch parsimony | Lineage-specific rate shifts | Combines stochastic + parsimony | Less commonly used; older |
| DupliPHY-Family (Ames et al 2012) | Per-family birth-death | Ancestral counts per family | Family-level granularity | Older; less integrated with modern OrthoFinder |
| ALE / GeneRax / AleRax (Szöllősi 2013; Morel 2024) | Per-family DTL reconciliation | Per-family D/T/L event counts | Direct integration with [[gene-tree-species-tree-reconciliation]] | Slower; per-family rather than across-family |
| CAFExp / CAFE 4.2 (DEPRECATED) | Earlier CAFE | -- | Historical | Use CAFE5 |
| Whale.jl with WGD (Zwaenepoel 2019) | Bayesian DTL+WGD | WGD-aware family dynamics | Native WGD integration | Julia ecosystem |
| Functional enrichment downstream | clusterProfiler / topGO on expanded/contracted | GO/KEGG enrichment | Standard | Multiple testing across families |

Methodology evolves; CAFE5 is the modern standard; verify the current CAFE5 manual (hahnlab/CAFE5) and Hahn lab papers before locking on a single approach.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Standard CAFE-style birth-death analysis | CAFE5 with gamma rate categories | Modern standard; handles rate variation |
| Annotation-pipeline-heterogeneity | CAFE5-error mode | Explicit error modeling |
| Post-WGD retention bias | CAFE5 + DupGen_finder classification + functional enrichment | Combine birth-death with WGD-specific analysis |
| Lineage-specific gene-family-rate shifts correlated with phenotype | CAFE5 with binary phenotype | Standard CAFE workflow |
| Ancestral gene-family counts at internal nodes | Count ASR | Per-node count posteriors |
| Test for "fast-evolving" family on specific lineage | CAFE5 lambda-tree (per-clade lambda) | Compares lambdas across clades |
| Functional enrichment in expanded families | clusterProfiler / topGO on expansion lists | Standard |
| HGT-affected families (prokaryotes) | ALE / GeneRax per-family DTL (see [[gene-tree-species-tree-reconciliation]]) | DTL framework explicit |
| Test if all families share single lambda | CAFE5 global lambda hypothesis | Restricted model |
| Specific family analysis (e.g. immune gene family) | ALE / AleRax per-family | Per-family detail |
| Distinguish gain from loss as primary driver | Count separate D and L counts | Standard parsimony |
| Convergent gene-family-rate shifts | RERconverge on family-count vectors | Trait-correlated rate shifts |
| Identify ancestral pan-clade family complement | CAFE5 ASR at MRCA | Pre-radiation family complement |
| Sub-clade-specific expansions in plant genomes | CAFE5 with clade-specific lambda | Compare angiosperm to gymnosperm |

## Per-Tool Failure Modes

### Annotation heterogeneity inflating expansions / contractions

**Trigger:** Running CAFE5 on counts from genomes annotated by different pipelines (Augustus, MAKER, BRAKER, NCBI).

**Mechanism:** Different annotation tools predict different numbers of genes per family; the same biological gene family may be annotated with 5 genes in BRAKER and 8 in MAKER. CAFE5 sees the difference as an "expansion" in the MAKER-annotated species (Tonkin-Hill 2020 documented this for bacterial pangenomes; same principle in eukaryotes).

**Symptom:** "Most expanded families" cluster in species annotated by a single pipeline (often more permissive tool); per-species "expansion rate" correlates with annotation pipeline rather than biology.

**Fix:** Re-annotate all genomes with a single pipeline (currently BRAKER3 or Funannotate for eukaryotes; Bakta for bacteria) before CAFE5. Alternatively use CAFE5-error mode with explicit error rates per species. Document annotation pipeline + version per species.

### Assembly fragmentation creating false contractions

**Trigger:** Including draft assemblies with low N50 in CAFE5 analysis.

**Mechanism:** Fragmented assemblies miss genes; the same family appears with fewer genes than expected. CAFE5 reports this as a contraction in the affected species.

**Symptom:** "Contracted" families in fragmented assemblies; correlation between BUSCO completeness and CAFE5 "contraction" rate; per-species missing-gene count varies 5-10x.

**Fix:** Require >= 90% BUSCO/Compleasm completeness for inclusion. Exclude species with > 5% lower BUSCO than median. Document N50 + BUSCO per assembly. For unavoidably fragmented assemblies, exclude from CAFE5 or use CAFE5-error with empirical error estimates.

### CAFE5 lambda non-convergence

**Trigger:** CAFE5 reports "lambda did not converge"; lambda jumps between values across runs.

**Mechanism:** Insufficient data (< 100 families); strong rate heterogeneity not captured; or input tree non-ultrametric.

**Symptom:** lambda estimate unstable; AIC of model selection variable.

**Fix:** Require >= 100 orthogroups (preferably > 1000) in input. Ensure tree is ultrametric (`ape::chronos()` or `treePL`); CAFE5 expects time-scaled tree. Use gamma rate categories (CAFE5 default `-k 4`) for rate heterogeneity. If still non-convergent, restrict to single-copy or small families; check tree branch lengths for negative values.

### Gamma rate-category misinterpretation

**Trigger:** Reporting "gamma rate categories" as biological gene-family clusters.

**Mechanism:** CAFE5 gamma categories are statistical buckets representing rate heterogeneity across families; they're not "fast-evolving family clusters" with biological meaning per se.

**Symptom:** Confusion in interpretation; "category-1 families are special."

**Fix:** Treat gamma categories as a statistical device. Report per-family lambda (estimated under per-family rate model) or per-clade lambda. Functional interpretation comes from family-level statistical tests, not category assignment.

### Multiple-testing across many families

**Trigger:** Reporting "significant" expansions / contractions without correction.

**Mechanism:** With ~10,000 families tested, ~500 will be significant at p=0.05 under H0. Without correction, false-discovery rate is high.

**Symptom:** Long lists of "expanded" families; functional enrichment dominated by chance hits.

**Fix:** Apply FDR (Benjamini-Hochberg) across families; or restrict to a priori hypothesized families. CAFE5 reports per-family p-values; FDR-correct downstream.

### Tree non-ultrametric / negative branches

**Trigger:** Using ML tree directly (substitutions per site) as input to CAFE5.

**Mechanism:** CAFE5 expects ultrametric (time-calibrated) tree; substitution-based branch lengths are not time-calibrated and may even have negative branches after rate variation.

**Symptom:** CAFE5 errors with "negative branch length"; or produces unreliable lambda estimates.

**Fix:** Time-calibrate the tree: use `ape::chronos()`, treePL (Smith & O'Meara 2012 Bioinformatics 28:2689), or LSD2 (To et al 2016 Syst Biol 65:82; gascuel-lab/LSD2) for fast NP-like calibration. Or use existing time-calibrated tree from TimeTree database.

### Outlier-family-driven lambda estimate

**Trigger:** Including extreme-size families (e.g., NLR resistance gene clusters in plants with 500+ members) in CAFE5.

**Mechanism:** Birth-death model's likelihood is dominated by large families; one or two outlier families can overwhelm the global lambda estimate.

**Symptom:** Excluding the top-5 largest families changes lambda by > 30%; lambda confidence intervals huge.

**Fix:** Robust analysis: report lambda with and without outliers; consider per-family lambda for largest families. Functional annotation of outliers reveals if they're biologically expected expansions (rapidly evolving gene families).

### Convergent gene-family-rate shifts not captured

**Trigger:** Phylogenomic question is whether multiple independent lineages show similar gene-family-rate shifts.

**Mechanism:** CAFE5 doesn't natively test for convergent rate shifts; per-clade lambda is the closest analog.

**Symptom:** Manual inspection of expansions across independent lineages shows pattern; CAFE5 doesn't formalize it.

**Fix:** Combine CAFE5 per-family lambdas with RERconverge (Redlich et al 2024 MBE 41:msae210) for trait-correlated rate shifts; use CSUBST (Fukushima 2023 Nat Eco Evo 7:155) for convergent substitution patterns.

### CAFE5 vs ALE for HGT-affected family

**Trigger:** Applying CAFE5 to bacterial families where transfer is common.

**Mechanism:** CAFE5 models birth-death (D, L) only; ignores transfer (T). For HGT-affected families, the dynamics include T events that CAFE5 cannot capture.

**Symptom:** CAFE5 lambda for bacterial families is implausibly high; family counts don't fit birth-death model.

**Fix:** Use ALE / GeneRax / AleRax for HGT-affected families ([[gene-tree-species-tree-reconciliation]]); CAFE5 is appropriate for vertical-inheritance-dominated families. Combine: CAFE5 for the bulk; ALE for HGT-confirmed families.

### Family classification (single-copy / multi-copy) affecting interpretation

**Trigger:** Treating single-copy orthogroups (always = 1 per species) and multi-copy uniformly.

**Mechanism:** Single-copy OGs have no count variation; including them dilutes the analysis. Multi-copy families show meaningful variation.

**Symptom:** Including single-copy OGs lowers lambda; per-family results dominated by them.

**Fix:** Restrict CAFE5 input to multi-copy orthogroups (max count >= 2 in any species). Filter via OrthoFinder HOG-classification single-copy column.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Minimum orthogroups for CAFE5 | >= 100; preferably > 1000 | Mendes et al 2020 Bioinformatics 36:5516 |
| Maximum gene-family count | depends on tree depth; typical < 1000 per family | Practical |
| BUSCO/Compleasm completeness | >= 90% for inclusion | Standard |
| FDR for expansions / contractions | q < 0.05 (Benjamini-Hochberg) | Standard |
| Significant lambda shift | LRT p < 0.05 / corrected for clade tests | CAFE5 LRT |
| Gamma categories | 4 default (`-k 4`) | CAFE5 docs |
| Tree ultrametricity | required; calibrate with ape::chronos or treePL | CAFE5 requirement |
| Per-clade lambda significance | LRT vs global lambda | CAFE5 docs |
| Annotation pipeline | one pipeline across all species; document version | Best practice |
| Functional enrichment q-value | q < 0.05 (BH-corrected) | Standard |
| Multi-copy family filter | max count across species >= 2 | Practical |
| Tree time-calibration source | TimeTree, treePL, LSD2 | Standard alternatives |
| Excluded singletons | yes, in some workflows | Standard |
| Outlier family exclusion threshold | top 5% by max count, exclude or per-family lambda | Robust |
| Convergent shift detection | RERconverge or CSUBST | Companion methods |

## CAFE5 Standard Workflow

**Goal:** Identify gene families with lineage-specific rate shifts (expansions / contractions).

**Approach:** Prepare ortholog count matrix -> ultrametricize tree -> run CAFE5 with gamma categories -> FDR-correct -> annotate.

```bash
# 1. Generate ortholog matrix from OrthoFinder v3 HOG file (inline Python; no external script needed)
python3 <<'PY'
import pandas as pd
hog = pd.read_csv('Phylogenetic_Hierarchical_Orthogroups/N0.tsv', sep='\t')
sp_cols = [c for c in hog.columns if c not in ('HOG', 'OG', 'Gene Tree Parent Clade')]
counts = pd.DataFrame({sp: hog[sp].fillna('').apply(
    lambda x: len(str(x).split(',')) if x and ',' in str(x) else (1 if str(x).strip() else 0)
) for sp in sp_cols})
counts.insert(0, 'family_id', hog['HOG'])
counts.insert(0, 'Description', '(null)')
counts.to_csv('cafe_input.tsv', sep='\t', index=False)
PY

# Output format: Description  family_id  Species1  Species2  ... (counts per species)

# 2. Ultrametricize tree
Rscript -e "
library(ape)
tree <- read.tree('SpeciesTree_rooted.txt')
tree_ultra <- chronos(tree)
write.tree(tree_ultra, 'tree_ultrametric.nwk')
"

# 3. Run CAFE5 with gamma categories
cafe5 \
    -i cafe_input.tsv \
    -t tree_ultrametric.nwk \
    -p \
    -k 4 \
    -e \
    -o cafe_output
# `-e` (no argument) lets CAFE5 estimate a global error model. To supply a pre-built
# error model file, use `-eerror.txt` (CAFE5 concatenates the flag and argument).
# Verify with `cafe5 --help`.

# Output:
#   cafe_output/Base_results.txt        Per-family results
#   cafe_output/Base_clade_results.txt   Per-clade lambda
#   cafe_output/Base_asr.tre            ASR tree
#   cafe_output/Base_clade_results.txt   LRT against null

# 4. FDR-correct across families (inline Python; see filter_significant_families() below)
```

```python
'''Filter expanded / contracted families with FDR-correction.'''
import pandas as pd
from statsmodels.stats.multitest import multipletests


def filter_significant_families(cafe_results_path, change_table_path, fdr_threshold=0.05):
    '''Identify families with significant rate shifts and their per-branch expansion/contraction direction.

    CAFE5 outputs (Base_family_results.txt or similar) include a per-family p-value (vs the null
    of a single global lambda). The per-branch expansion/contraction comes from a separate
    `Base_change.tab` file (per-family x per-branch count change). Per-family lambda is NOT
    a column in Base_family_results.txt; that's a global parameter.
    '''
    df = pd.read_csv(cafe_results_path, sep='\t')
    pcol = 'p-value' if 'p-value' in df.columns else 'pvalue'
    df['fdr'] = multipletests(df[pcol].fillna(1.0), method='fdr_bh')[1]
    significant = df[df['fdr'] < fdr_threshold]

    # Pair with per-branch change table; positive cell = expansion on that branch
    changes = pd.read_csv(change_table_path, sep='\t')
    return {'significant_families': significant, 'per_branch_changes': changes}
```

## CAFE5-error for Annotation Heterogeneity

```bash
# Prepare per-species error rates (from BUSCO completeness or empirical)
cat > error_rates.tsv << 'EOF'
species_id    error_rate
Species_A     0.05
Species_B     0.08
Species_C     0.03
EOF

# Run with error-aware mode (supply pre-built error model file).
# CAFE5 concatenates the -e flag with its argument (no space): -e<error_model_file>
cafe5 \
    -i cafe_input.tsv \
    -t tree_ultrametric.nwk \
    -p \
    -k 4 \
    -eerror_model.txt \
    -o cafe_error_output
```

## Count for Per-Branch Ancestral State

```bash
# Count requires Java
java -jar Count.jar -i count_format.tsv -t species_tree.nwk \
    -o count_output --ancestral_counts

# Per-branch D and L counts in count_output/branch_events.tsv
```

## Functional Enrichment Downstream

```r
library(clusterProfiler)
library(org.Hs.eg.db)  # or appropriate species DB

# Enrichment of expanded families
expanded_genes <- read.csv('expanded_genes.tsv', stringsAsFactors = FALSE)$gene_id

go_enrich <- enrichGO(gene = expanded_genes,
                      OrgDb = org.Hs.eg.db,
                      ont = 'BP',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05)

kegg_enrich <- enrichKEGG(gene = expanded_genes,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)
```

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| CAFE5 significant; ALE shows no DTL signal | CAFE5 detects count change; ALE detects events | Both can be true; CAFE5 is count-level, ALE is event-level |
| CAFE5 significant in bacterial family | Possibly HGT-driven; CAFE5 ignores T | Re-run with ALE for HGT-affected families |
| Per-clade lambda differs in CAFE5 vs uniform lambda | Real rate heterogeneity | Trust per-clade |
| Annotation heterogeneity hypothesis | Re-annotation eliminates "expansion" | Confirm annotation artifact; report as such |
| CAFE5 + RERconverge agree on expansion + trait shift | Convergent biological mechanism | Strong evidence |
| Outlier-family-driven global lambda | Top-5 families distort estimate | Report robust lambda; manually flag outliers |
| CAFE5 expansion in fragmented species | False; assembly fragmentation | Re-assess BUSCO; exclude or correct |
| Single-copy family flagged as significant | Statistical artifact; no variation | Exclude single-copy OGs; restrict to multi-copy |
| Whale.jl WGD branch shows D burst, CAFE5 shows expansion | Same event; WGD modeling preferred | Whale.jl is more biological for known WGD |
| Count parsimony vs CAFE5 likelihood disagree | Parsimony underestimates losses | Trust CAFE5 likelihood |

**Operational rule for publication:** CAFE5 with gamma rate categories + annotation pipeline normalized + BUSCO completeness > 90% + FDR-corrected significance + functional enrichment of expanded/contracted + (for bacteria) ALE complement = publication-grade gene-family evolution analysis.

## Cohort Gotchas

- **WGD lineages:** post-WGD retention bias; gene balance hypothesis (Birchler & Veitia 2007); analyze with [[whole-genome-duplication]] context
- **Plant gene families:** NLR clusters (resistance) and ribosomal proteins are inherently large; expect lineage variation
- **Mammalian gene families:** olfactory receptors are highly variable; expect lineage-specific changes
- **Bacterial gene families:** HGT-driven dynamics; use ALE/AleRax instead of CAFE5
- **Polyploid species:** subgenome assignment first; analyze each subgenome separately
- **Rapidly evolving lineages:** higher branch-specific rates; per-clade lambda model
- **Conserved species (e.g., extant cyanobacteria):** lambda may be very low; few changes to detect
- **Recent radiations:** insufficient time for divergence; CAFE5 may have low power
- **Highly fragmented MAGs:** include only high-quality MAGs

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Annotation pipeline?" | Single pipeline (BRAKER3 / Funannotate / Bakta) across all species; version pinned |
| "BUSCO completeness?" | >= 90% required; per-species reported |
| "Assembly fragmentation?" | N50 reported; species with > 5% lower BUSCO than median excluded |
| "Multiple testing?" | FDR (Benjamini-Hochberg) applied across families |
| "Lambda non-convergence?" | CAFE5 with -k 4 gamma categories; tree ultrametricized |
| "Outlier families?" | Robust analysis with and without; outliers individually annotated |
| "HGT in bacteria?" | ALE / GeneRax cross-checked for HGT-affected families |
| "Functional enrichment?" | clusterProfiler / topGO with FDR; pathway-level interpretation |
| "Tree time-calibration?" | TimeTree-based or treePL/LSD2; documented |
| "WGD effects?" | DupGen_finder + WGD-specific analysis; subgenome-aware |
| "Convergent shifts?" | RERconverge complementary analysis |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| CAFE5 "lambda did not converge" | Insufficient data or non-ultrametric tree | Ultrametricize; increase families; check tree |
| CAFE5 "negative branch length" | ML tree input | Time-calibrate with chronos / treePL / LSD2 |
| CAFE5 errors on input format | Wrong column order | Check OrthoFinder HOG format; species in column 2 onward |
| Count GUI hangs | Java memory | Increase Java heap: `java -Xmx16G -jar Count.jar` |
| BadiRate gamma errors | Initialization issue | Use default `-g 1` or specify per family |
| Per-clade lambda implausibly high | Outlier family or tree issue | Exclude top families; re-run |
| All families significant | No FDR correction | Apply BH correction |
| Annotation heterogeneity not addressed | Mixed pipelines | Re-annotate consistently |
| Family with 0 count for all species | Annotation issue | Filter rows with all zeros |
| Highly variable family count (e.g. 500-5000) | Real biological variation or annotation | Annotate manually; consider exclusion |

## Tool Installation Notes

```bash
# CAFE5
conda install -c bioconda cafe
# Or: git clone https://github.com/hahnlab/CAFE5

# Count
wget http://www.iro.umontreal.ca/~csuros/gene_content/count.tar.gz
tar xf count.tar.gz

# BadiRate
git clone https://github.com/PauloRoldan/badirate

# DupliPHY-Family
# Web only; no public CLI

# Whale.jl (Julia) — see [[gene-tree-species-tree-reconciliation]]
julia -e 'using Pkg; Pkg.add("Whale")'

# R packages
install.packages(c('ape', 'phytools', 'clusterProfiler', 'org.Hs.eg.db'))

# Time calibration
conda install -c bioconda treepl
# Or use ape::chronos (R)

# For OrthoFinder input
conda install -c bioconda orthofinder
```

For Funannotate / BRAKER3 reannotation (essential pre-CAFE5):
```bash
conda install -c bioconda funannotate braker3
```

## References

- Hahn MW et al 2005 Genome Res 15:1153 (CAFE original framework)
- Mendes FK et al 2020 Bioinformatics 36:5516 (CAFE5)
- Csurös M 2010 Bioinformatics 26:1910 (Count)
- Librado P et al 2012 Bioinformatics 28:279 (BadiRate)
- Ames RM et al 2012 Bioinformatics 28:48 (DupliPHY-Family)
- Tonkin-Hill G et al 2020 Genome Biol 21:180 (Panaroo; annotation heterogeneity)
- Smith SA & Dunn CW 2008 Bioinformatics 24:715 (Phyutility)
- Smith SA & O'Meara BC 2012 Bioinformatics 28:2689 (treePL)
- To T-H, Jung M, Lycett S & Gascuel O 2016 Syst Biol 65:82 (LSD2)
- Birchler JA & Veitia RA 2007 Plant Cell 19:395 (gene balance)
- Redlich R et al 2024 MBE 41:msae210 (RERconverge categorical)
- Fukushima K & Pollock DD 2023 Nat Eco Evo 7:155 (CSUBST)
- Lynch M & Conery JS 2000 Science 290:1151 (gene duplication mechanism)
- Force A et al 1999 Genetics 151:1531 (subfunctionalization)
- De Bie T et al 2006 Bioinformatics 22:1269 (CAFE 2)
- Han MV et al 2013 MBE 30:1987 (CAFE 3)
- Otto SP & Whitton J 2000 Annu Rev Genet 34:401 (polyploidy mechanisms)
- TimeTree (database, http://www.timetree.org)

## Related Skills

- comparative-genomics/ortholog-inference - OrthoFinder HOG matrix is CAFE5 input
- comparative-genomics/gene-tree-species-tree-reconciliation - ALE per-family DTL, complement to CAFE5
- comparative-genomics/whole-genome-duplication - Post-WGD retention bias context
- comparative-genomics/positive-selection - Selection within expanded families
- comparative-genomics/ancestral-reconstruction - Ancestral count reconstruction
- phylogenetics/divergence-dating - Time-calibrated tree for CAFE5
- phylogenetics/modern-tree-inference - Species tree input
- pathway-analysis/go-enrichment - Functional enrichment of expanded families
- pathway-analysis/gsea - GSEA on family expansions
- single-cell/cell-annotation - Cell-type-specific gene-family expansions
