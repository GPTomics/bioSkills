---
name: bio-workflows-outbreak-pipeline
description: Orchestrates genomic-epidemiology outbreak investigation from pathogen isolates to transmission networks, forking bacterial (snippy -> Gubbins recombination-masking -> IQ-TREE -> TreeTime -> TransPhylo) vs viral (Nextstrain/augur), with parallel MLST typing (cgMLST delegated to epidemiological-genomics/pathogen-typing) and AMR surveillance. Use when committing ONE reference genome for SNP calling (every isolate and distance inherits its coordinates), applying MANDATORY Gubbins recombination-masking on core.full.aln before the tree for recombining bacteria (skipping it inflates the clock 2-5x), gating time-scaling on a temporal-signal test (TempEst R2 >= 0.3), using a pathogen- AND population-specific cluster threshold rather than a universal SNP cutoff, or pinning pangolin-data/Nextclade/Freyja versions for the viral route. Hands mechanism to the epidemiological-genomics component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: mlst
workflow: true
depends_on:
  - epidemiological-genomics/pathogen-typing
  - epidemiological-genomics/amr-surveillance
  - epidemiological-genomics/phylodynamics
  - epidemiological-genomics/transmission-inference
  - epidemiological-genomics/variant-surveillance
qc_checkpoints:
  - after_typing: "Valid ST assigned (7-locus mlst); cgMLST via chewBBACA is deferred to epidemiological-genomics/pathogen-typing"
  - after_amr: "AMR genes identified; AMRFinderPlus applies its curated per-gene thresholds (it does not use a flat 90% identity cutoff)"
  - after_phylodynamics: "Root-to-tip R2 >=0.3, clock rate plausible"
  - after_transmission: "Transmission pairs consistent with epi data"
---

## Version Compatibility

Reference examples tested with: ncbi-amrfinderplus 4.0+, hamronization 1.1+, tb-profiler 6.2+, mlst 2.23+, chewBBACA 3.3+, pangolin 4.3+ (pangolin-data 1.30+), nextclade 3.8+, snippy 4.6+, gubbins 3.3+, clonalframeml 1.13+, IQ-TREE 2.3.6+, TreeTime 0.11+, BEAST 2.7.6+ (BDSKY 1.5+, MASCOT 3.0+, BICEPS), TransPhylo 1.4+ (R), outbreaker2 1.2+ (R), bactdating 1.1+ (R), freyja 1.4+, mob_suite 3.1+, BioPython 1.84+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>`; CLI: `<tool> --version` then `--help`
- R: `packageVersion('<pkg>')` then `?function_name`
- Pangolin: `pangolin --all-versions` (records pangolin + pangolin-data + scorpio + constellations)
- Nextclade: `nextclade dataset list --tag latest sars-cov-2`
- TB-Profiler: `tb-profiler list_db` (verify WHO catalogue edition)

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Outbreak Pipeline

**"Characterize a pathogen outbreak from my isolate sequences"** -> Orchestrate MLST typing, SNP phylogeny, TreeTime time-scaled tree construction, TransPhylo transmission inference, AMR profiling, and variant surveillance for genomic epidemiology.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

A transmission tree is the MAP estimate among many equally probable trees, and its trustworthiness is decided at four seams.

1. **One reference genome for SNP calling, committed before any isolate.** snippy calls SNPs against ONE reference; every isolate, the core alignment, and the cluster distances inherit its coordinates. A build/coordinate mismatch (isolates called against different references, or an AMR panel keyed to a different build) fabricates or hides SNPs. Commit the reference (GenBank format for snippy) up front.
2. **Recombination-masking is MANDATORY for recombining bacteria, on core.full.aln, before the tree.** Gubbins masking is not optional for S. pneumoniae, N. gonorrhoeae, E. coli, Klebsiella, Campylobacter, H. pylori — skipping it inflates the clock rate 2-5x and biases R_e. Input MUST be `core.full.aln` (full positions incl. invariant), NOT `core.aln` (variable-only). Skip masking ONLY for documented-clonal Mtb.
3. **Time-scaling requires a passing temporal-signal test, committed before trusting any dated tree.** TempEst root-to-tip regression (Rambaut 2016) with R2 >= 0.3 as a field convention (the paper sets no threshold and cautions R2 is not a valid significance test); the date-randomisation test is a secondary check (can pass with narrow sampling windows), NOT a substitute. The outbreak-scale clock is lineage/host-specific, not a universal constant.
4. **The cluster threshold is pathogen- AND population-specific — never a universal cutoff.** From the literature per pathogen (Mtb <=12/<=5 SNP; S. aureus <=15; K. pneumoniae KPC <=21; C. difficile <=2 masked; Salmonella <=5 cgMLST alleles). Walker's 5-SNP Mtb threshold was calibrated in low-transmission UK and inflates apparent recent transmission 2-5x in high-burden settings. A genomic distance is not an epidemiological distance without a time-scaled prior.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Reference genome + build (bacterial) | Every isolate's SNPs, the core alignment, cluster distances; a mismatch fabricates/hides SNPs |
| Recombination-masking scheme (Gubbins on core.full.aln) | The clock rate and R_e; skipping inflates the clock 2-5x for recombining taxa |
| Clock + temporal-signal (TempEst R2 >= 0.3) | Whether the dated tree is supported at all |
| Cluster threshold (pathogen + population-specific) | Who is "linked"; a universal SNP cutoff over-clusters in high-burden settings |
| Viral tool versions (pangolin-data/Nextclade/Freyja) | The lineage call; same genome, different call across versions — pin them |

## Workflow Overview

```
Pathogen Isolate Genomes (FASTA/FASTQ) + collection dates + (optional) contact data
        |
        v
   +---------+---------+
   |                   |
   v                   v
[1a. MLST + serotyping    [1b. AMR + species mode:
     + Pangolin/UShER          AMRFinderPlus --organism,
     for SARS-CoV-2;           TB-Profiler for Mtb,
     cgMLST ->                 hAMRonization across tools]
     pathogen-typing]
   |                   |
   +--------+----------+
            |
            v
[2. snippy + snippy-core (bacteria) -> Gubbins on core.full.aln to mask recombination
    (mandatory for bacteria; skip only for clonal Mtb)]
            |
            v
[3. IQ-TREE on recombination-masked alignment + TempEst R^2 >= 0.3 + date-randomisation;
    TreeTime --coalescent skyline --clock-filter 4 OR BactDating;
    BEAST2 BDSKY (origin > rootHeight, multi-chain) for posterior R_e]
            |
            v
[4. Transmission inference: outbreaker2 (dense + contact data) OR TransPhylo (sparse,
    from dated tree) OR transcluster (pair-level probability); pathogen-specific SNP
    threshold for cluster definition -- NEVER a universal cutoff]
            |
            v
Transmission tree posterior + R_e(t) + lineage / clone context + AMR phenotype
```

## Prerequisites

```bash
conda install -c bioconda mlst chewbbaca ncbi-amrfinderplus hamronization tb-profiler \
    snippy snp-dists gubbins clonalframeml iqtree treetime pangolin nextclade freyja \
    mob_suite plasmidfinder sistr_cmd seqsero2 kleborate kaptive seroba

conda install -c bioconda beast2
packagemanager -add BDSKY BEASTLabs feast ORC MASCOT BICEPS

Rscript -e "install.packages(c('TransPhylo', 'outbreaker2', 'BactDating', 'bdskytools', 'coda', 'ape'))"

amrfinder -u
tb-profiler update_tbdb
```

## Primary Path: Bacterial Outbreak Investigation

### Step 1a: MLST Typing (Parallel)

**Goal:** Assign 7-locus PubMLST sequence types to all isolates for clonal-context interpretation.

**Approach:** Run Seemann's `mlst` per assembly; auto-detect scheme; concatenate the per-isolate output into a cohort TSV.

```bash
#!/bin/bash
ISOLATES="isolate1.fasta isolate2.fasta isolate3.fasta"
OUTDIR="outbreak_results"
mkdir -p ${OUTDIR}/{mlst,amr,alignment,phylo,transmission}

# Run MLST on all isolates
echo "=== MLST Typing ==="
for fasta in $ISOLATES; do
    sample=$(basename $fasta .fasta)
    mlst $fasta > ${OUTDIR}/mlst/${sample}.mlst.txt
done

# Combine results
cat ${OUTDIR}/mlst/*.mlst.txt > ${OUTDIR}/mlst/all_mlst.tsv
echo "MLST complete: ${OUTDIR}/mlst/all_mlst.tsv"
```

### Step 1b: AMR Detection (Parallel) -- AMRFinderPlus with species mode

**Goal:** Produce per-isolate AMR calls with species-specific point-mutation panel activated, then harmonise across the cohort to the PHA4GE schema for cross-lab comparison.

**Approach:** AMRFinderPlus `-n` for nucleotide assembly with `--organism` and `--plus`; pipe each per-isolate TSV through `hamronize amrfinderplus` with mandatory PHA4GE metadata; `hamronize summarize` merges to a cohort table. For *M. tuberculosis*, switch to TB-Profiler -- AMRFinderPlus has no Mtb organism mode.

```bash
echo "=== AMR Detection ==="
SPECIES="Klebsiella_pneumoniae"
for fasta in $ISOLATES; do
    sample=$(basename $fasta .fasta)
    amrfinder -n $fasta --organism $SPECIES --plus --threads 8 \
        -o ${OUTDIR}/amr/${sample}.amrfinder.tsv

    hamronize amrfinderplus \
        --analysis_software_version $(amrfinder -V | awk '/Software/{print $NF}') \
        --reference_database_version $(amrfinder -V | awk '/Database/{print $NF}') \
        --input_file_name ${sample} \
        ${OUTDIR}/amr/${sample}.amrfinder.tsv > ${OUTDIR}/amr/${sample}.hamr.tsv
done

hamronize summarize -t tsv -o ${OUTDIR}/amr/cohort.hamr.tsv ${OUTDIR}/amr/*.hamr.tsv
echo "AMR summary: ${OUTDIR}/amr/cohort.hamr.tsv"
```

For *M. tuberculosis*, route to TB-Profiler instead -- AMRFinderPlus has no Mtb organism mode. For colistin / mcr surveillance and any plasmid-mobility claim, follow with MOB-suite (`mob_recon` + `mob_typer`) to determine plasmid context. See epidemiological-genomics/amr-surveillance for the full decision tree.

### Step 2: Core Genome Alignment + Recombination Masking (Bacteria)

**Goal:** Build a recombination-aware core-genome alignment that is safe for downstream clock inference.

**Approach:** Snippy per isolate against the reference; snippy-core to merge into the core alignment; Gubbins on `core.full.aln` (NOT `core.aln`) to mask recombinant tracts. Skipping recombination masking inflates the clock rate 2-5x for recombining bacteria (S. pneumoniae, N. gonorrhoeae, E. coli, Klebsiella, Campylobacter, H. pylori); the date-randomisation test is NOT a sufficient guard.

```bash
echo "=== Core Genome Alignment ==="
REFERENCE="reference.gbk"  # Reference genome in GenBank format

# Run snippy for each isolate
for fasta in $ISOLATES; do
    sample=$(basename $fasta .fasta)
    snippy --outdir ${OUTDIR}/alignment/snippy_${sample} \
           --ref $REFERENCE \
           --ctgs $fasta \
           --cpus 8
done

# Core SNP alignment
snippy-core --ref $REFERENCE --prefix core ${OUTDIR}/alignment/snippy_*

# Mandatory for recombining bacteria (S. pneumoniae, N. gonorrhoeae, E. coli, Klebsiella,
# Campylobacter, H. pylori). Skip ONLY for clonal Mtb cross-lineage analyses where
# recombination is documented to be rare; even then a recombination check is defensible.
# Input MUST be core.full.aln (full positions including invariant); core.aln (variable-only)
# gives wrong recombination calls because Gubbins cannot estimate background SNP density.
run_gubbins.py --prefix gubbins core.full.aln

mv core.* gubbins.* ${OUTDIR}/alignment/
echo "Recombination-masked alignment: ${OUTDIR}/alignment/gubbins.filtered_polymorphic_sites.fasta"
```

### Step 3: Phylodynamics with TreeTime

**Goal:** Time-scale the recombination-masked phylogeny with a global clock-rate estimate, gated by temporal-signal QC.

**Approach:** IQ-TREE on the recombination-masked alignment with `+ASC` ascertainment correction; TreeTime with coalescent skyline prior and `--clock-filter 4`; inspect `root_to_tip_regression.pdf` BEFORE trusting downstream output (R^2 >= 0.3 minimum as a field convention; TempEst sets no threshold).

```python
import subprocess
from Bio import Phylo, AlignIO
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

outdir = Path('outbreak_results')

# Build ML tree on the recombination-masked alignment with ascertainment-bias correction
# +ASC is required because the input contains variable positions only post-Gubbins
subprocess.run([
    'iqtree2', '-s', str(outdir / 'alignment/gubbins.filtered_polymorphic_sites.fasta'),
    '-m', 'GTR+G+ASC', '-B', '1000', '-bnni', '-T', 'AUTO',
    '--prefix', str(outdir / 'phylo/outbreak')
], check=True)

# Prepare metadata with dates
# Format: name\tdate (YYYY-MM-DD or decimal year)
metadata = pd.DataFrame({
    'name': ['isolate1', 'isolate2', 'isolate3', 'isolate4', 'isolate5'],
    'date': ['2024-01-15', '2024-01-22', '2024-02-01', '2024-02-10', '2024-02-15']
})
metadata.to_csv(outdir / 'phylo/metadata.tsv', sep='\t', index=False)

# Run TreeTime
subprocess.run([
    'treetime',
    '--tree', str(outdir / 'phylo/outbreak.treefile'),
    '--aln', str(outdir / 'alignment/gubbins.filtered_polymorphic_sites.fasta'),
    '--dates', str(outdir / 'phylo/metadata.tsv'),
    '--outdir', str(outdir / 'phylo/treetime_output'),
    '--coalescent', 'skyline',
    '--clock-filter', '4',  # SD multiplier for TreeTime's clock filter
    '--confidence',
    '--reroot', 'best'
], check=True)

# Temporal-signal QC: inspect root_to_tip_regression.pdf BEFORE trusting any downstream output.
# R^2 >= 0.3 minimum (field convention, NOT from Rambaut 2016 -- TempEst sets no threshold and
# states R^2 is an informal dispersion measure, not a significance test). If R^2 < 0.3, time-scaling is not
# supported -- report uncertainty and consider extending the sampling window. The
# date-randomisation test is a secondary check; it can pass with narrow sampling windows
# (false negative).
print('TreeTime output:', outdir / 'phylo/treetime_output')
```

### Step 4: Transmission Inference with TransPhylo

**Goal:** Reconstruct the posterior who-infected-whom transmission tree and R_e from the dated phylogeny.

**Approach:** Convert the TreeTime dated tree to TransPhylo `ptree`; supply pathogen-tuned generation-time and sampling-time Gamma priors; run MCMC at >=1e5 iterations (10k is smoke-test only); summarise via medoid transmission tree and per-pair WIWS probability. For dense outbreaks with contact-tracing data, outbreaker2 with `ctd` is preferred over TransPhylo (genomic-only).

```r
library(TransPhylo)
library(ape)

# Load dated tree from TreeTime
tree <- read.nexus("outbreak_results/phylo/treetime_output/timetree.nexus")

# Set parameters
# dateT: date when sampling stopped
# w.shape, w.scale: generation time distribution (Gamma)
# For many bacteria: mean ~14 days, shape=2, scale=7
dateT <- 2024.2  # Decimal year when sampling ended (end of observation)
w_shape <- 2     # Generation time shape (Gamma)
w_scale <- 7/365 # Gamma SCALE = 7 days; mean generation time = shape*scale = 2*7 = ~14 days

# TransPhylo operates on a `ptree` (dated phylogeny + last-sample date), NOT a raw ape phylo;
# convert first or inferTTree errors on a NULL ptree$ptree/$nam.
ptree <- ptreeFromPhylo(tree, dateLastSample = dateT)

# Run TransPhylo with enough iterations for posterior convergence; 10k is a smoke-test only.
# For publication, run >=1e5 (small outbreaks) to >=1e6+ iterations and inspect trace plots.
res <- inferTTree(ptree, dateT = dateT,
                   w.shape = w_shape, w.scale = w_scale,
                   mcmcIterations = 1e5,
                   startNeg = 1, startPi = 0.5)

# medTTree returns a coloured transmission tree (ctree); plot it with plotCTree
med_ctree <- medTTree(res)

# Plot transmission tree
pdf("outbreak_results/transmission/transmission_tree.pdf", width=10, height=8)
plotCTree(med_ctree)
dev.off()

# Who infected whom matrix (same 0.5 burn-in as the R_e estimate below, so both discard pre-convergence)
wiw <- computeMatWIW(res, burnin = 0.5)
write.csv(wiw, "outbreak_results/transmission/who_infected_whom.csv")

# R_e estimate (effective reproduction number under current immunity / interventions).
# This is NOT R_0 (basic reproduction number in a fully susceptible population);
# the phylodynamics literature is explicit about this distinction.
# getOffspringDist(record, burnin, k) gives the per-case offspring distribution; average
# across sampled hosts for a cohort R_e (or use BEAST2 BDSKY for a posterior Re(t)).
# Host names come from the ptree (res has no $ttree$nam field).
offspring <- sapply(ptree$nam, function(k) mean(getOffspringDist(res, k = k, burnin = 0.5)))
# The interval is the 2.5-97.5% spread of per-host mean offspring (across-host dispersion), NOT a
# posterior credible interval (per-host posteriors were collapsed by mean() first).
cat("R_e estimate:", mean(offspring), "(across-host 2.5-97.5%:", quantile(offspring, 0.025), "-", quantile(offspring, 0.975), ")\n")
```

### Python Alternative: TransPhylo via rpy2

**Goal:** Drive the same TransPhylo workflow from Python pipelines that prefer not to fork into R.

**Approach:** rpy2 bridges into the TransPhylo R package with named-argument passing; same priors and MCMC iteration discipline apply.

```python
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import pandas as pd
from pathlib import Path

pandas2ri.activate()

transphylo = importr('TransPhylo')
ape = importr('ape')

outdir = Path('outbreak_results')

tree = ape.read_nexus(str(outdir / 'phylo/treetime_output/timetree.nexus'))

date_t = 2024.2
w_shape = 2
w_scale = 7/365

# Convert to a TransPhylo ptree before inference (inferTTree needs ptree, not a raw phylo).
ptree = transphylo.ptreeFromPhylo(tree, dateLastSample=date_t)
res = transphylo.inferTTree(ptree, dateT=date_t, w_shape=w_shape, w_scale=w_scale,
                             mcmcIterations=10000, startNeg=1, startPi=0.5)

# medTTree returns a ctree; hand it to R's global env and plot with plotCTree.
med_ctree = transphylo.medTTree(res)
ro.globalenv['med_ctree'] = med_ctree
ro.r(f'''
pdf("{outdir}/transmission/transmission_tree.pdf", width=10, height=8)
plotCTree(med_ctree)
dev.off()
''')

print(f'Transmission tree saved to {outdir}/transmission/')
```

## Visualization: Outbreak Timeline

**Goal:** Plot isolates over time coloured by sequence type to communicate cluster expansion and clonal context.

**Approach:** Merge collection-date metadata with MLST output; plot per-isolate timestamps as a strip chart with per-ST colour.

```python
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

metadata = pd.read_csv('outbreak_results/phylo/metadata.tsv', sep='\t')
metadata['date'] = pd.to_datetime(metadata['date'])

mlst = pd.read_csv('outbreak_results/mlst/all_mlst.tsv', sep='\t', header=None,
                    names=['file', 'scheme', 'ST'] + [f'locus{i}' for i in range(7)])
mlst['sample'] = mlst['file'].apply(lambda x: x.split('/')[-1].replace('.fasta', ''))

# Merge data
combined = metadata.merge(mlst[['sample', 'ST']], left_on='name', right_on='sample')

fig, ax = plt.subplots(figsize=(12, 6))

colors = {'ST11': 'red', 'ST258': 'blue', 'ST307': 'green'}
for st in combined['ST'].unique():
    subset = combined[combined['ST'] == st]
    ax.scatter(subset['date'], [1]*len(subset), label=f'ST{st}',
               s=100, c=colors.get(f'ST{st}', 'gray'), alpha=0.7)

ax.set_xlabel('Date')
ax.set_ylabel('')
ax.set_title('Outbreak Timeline by Sequence Type')
ax.legend()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('outbreak_results/outbreak_timeline.pdf')
```

## Parameter Recommendations

| Step | Parameter | Value | Rationale |
|------|-----------|-------|-----------|
| snippy | --mincov | 10 | Minimum coverage for variant call |
| Gubbins | input | core.full.aln | Full positions required to estimate background SNP density; core.aln is wrong |
| IQ-TREE | -m | GTR+G+ASC | +ASC ascertainment correction for SNP-only post-Gubbins input |
| TreeTime | --clock-filter | 4 | SD multiplier on root-to-tip residual; TreeTime convention |
| TreeTime | R^2 minimum | 0.3 | Below this, temporal signal treated as insufficient (field convention; TempEst itself sets no cutoff) |
| TransPhylo | w.shape, w.scale | 2, 7/365 | Gamma scale 7 days x shape 2 = ~14-day mean generation time; cite the pathogen-specific literature |
| TransPhylo | mcmcIterations | 1e5-1e6+ | 10k is a smoke-test only; inspect trace and ESS before reporting |
| BEAST2 BDSKY | origin | > rootHeight | Initialise to ~(tMRCA + 0.1*tMRCA); origin == rootHeight biases R_e upward (Stadler 2013) |
| BEAST2 chains | independent runs | >=3-4 | Single-chain ESS >=200 is necessary but not sufficient; combine after marginal overlap |
| Pangolin | --analysis-mode | usher | pangoLEARN deprecated mid-2023; UShER default since v4 (de Bernardi Schneider 2024, Virus Evol 10:vead085) |

## Pathogen-Specific SNP / cgMLST Cluster Thresholds

Cluster definition is pathogen- AND population-specific. NEVER apply a universal cutoff. See epidemiological-genomics/transmission-inference for full table with citations.

| Pathogen | Cluster threshold | Source |
|----------|-------------------|--------|
| *M. tuberculosis* (core SNP) | <=12 SNPs (likely transmission); <=5 (recent) | Walker 2013 *Lancet Infect Dis* 13:137 (UK low-transmission setting -- inflates 2-5x in high-burden) |
| *S. aureus* (core SNP) | <=15 SNPs (within hospital) | Coll 2020 *Lancet Microbe* 1:e328 |
| *K. pneumoniae* (KPC outbreak) | <=21 SNPs | Field convention; no threshold source |
| *C. difficile* (recombination-masked core SNP) | <=2 SNPs (likely direct) | Eyre 2013 *NEJM* 369:1195 |
| *Salmonella* (cgMLST, EnteroBase) | <=5 alleles | EnteroBase / EFSA convention |
| *Listeria* (PulseNet cgMLST) | <=4 alleles | PulseNet protocol |
| SARS-CoV-2 | NOT defined by SNP alone | 0-2 SNPs + epi link + sampling window |
| HIV-1 subtype B | 1.5% TN93 distance | HIV-TRACE US-CDC default (re-tune for non-B subtypes) |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Fabricated or hidden SNPs / wrong distances | Build/coordinate mismatch (isolates called against different references, or AMR panel on another build) | One committed reference for every isolate + the AMR panel; verify contig/seqid consistency |
| Clock inflated 2-5x, false transmission links | Recombination masking skipped | Gubbins on core.full.aln BEFORE the tree; the date-randomisation test is NOT a sufficient guard |
| Over-clustered outbreak in a high-burden setting | Universal SNP cutoff | Pathogen- AND population-specific threshold from the literature; a genomic distance is not an epidemiological distance |
| "Who infected whom" overclaimed | Transmission read from single-isolate SNP distances | "Transmission consistent with genomics"; single-isolate-per-host trees are under-identified (MAP among many) |
| Same genome, different lineage across labs/dates (viral) | pangolin-data/Nextclade/Freyja version churn | Pin the versions in metadata; re-run all samples against ONE version before comparing |
| Poor temporal signal | Insufficient sampling / recombination | Mask recombination (Gubbins); check dates; do not time-scale below TempEst R2 0.3 |
| Missing AMR genes | Database mismatch | Try multiple databases (ncbi/card/resfinder); report allele identity, not just family |

## Output Files

| File | Description |
|------|-------------|
| `mlst/all_mlst.tsv` | Sequence types for all isolates |
| `amr/cohort.hamr.tsv` | AMR gene presence/absence matrix |
| `alignment/core.aln` | Core genome SNP alignment |
| `phylo/outbreak.treefile` | ML phylogenetic tree |
| `phylo/treetime_output/` | Dated tree and molecular clock |
| `transmission/transmission_tree.pdf` | Inferred transmission network |
| `transmission/who_infected_whom.csv` | Transmission probability matrix |

## Related Skills

- database-access/sra-data - Download outbreak FASTQ from SRA / ENA
- database-access/ncbi-datasets-cli - Bulk-pull pathogen reference genomes (e.g. `datasets download virus`)
- epidemiological-genomics/pathogen-typing - MLST and cgMLST details
- epidemiological-genomics/amr-surveillance - AMRFinderPlus, ResFinder
- epidemiological-genomics/phylodynamics - TreeTime, BEAST2 parameters
- epidemiological-genomics/transmission-inference - TransPhylo configuration
- epidemiological-genomics/variant-surveillance - Nextclade for viral outbreaks
- phylogenetics/modern-tree-inference - IQ-TREE2 model selection

## References

- Croucher NJ, Page AJ, Connor TR, et al (2015) Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. *Nucleic Acids Research* 43:e15. DOI 10.1093/nar/gku1196. (recombination masking.)
- Didelot X, Fraser C, Gardy J, Colijn C (2017) Genomic infectious disease epidemiology in partially sampled and ongoing outbreaks (TransPhylo). *Molecular Biology and Evolution* 34:997-1007. DOI 10.1093/molbev/msw275. (transmission inference.)
- Walker TM, Ip CLC, Harrell RH, et al (2013) Whole-genome sequencing to delineate Mycobacterium tuberculosis outbreaks: a retrospective observational study. *Lancet Infectious Diseases* 13:137-146. DOI 10.1016/S1473-3099(12)70277-3. (5-SNP threshold, low-transmission calibration.)
- Sagulenko P, Puller V, Neher RA (2018) TreeTime: maximum-likelihood phylodynamic analysis. *Virus Evolution* 4:vex042. DOI 10.1093/ve/vex042.
