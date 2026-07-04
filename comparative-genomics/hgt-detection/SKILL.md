---
name: bio-comparative-genomics-hgt-detection
description: Detect horizontal gene transfer (HGT / LGT) using compositional methods (GC%, codon usage, tetranucleotide z-scores via SIGI-HMM, AlienHunter, IslandViewer 4, IslandPath-DIMOB), phylogenetic-incongruence methods (AvP, HGTphyloDetect, ALE / GeneRax / AleRax reconciliation, RANGER-DTL), and BLAST-distribution methods (HGTector v2, DarkHorse, Alien Index). Use when screening prokaryote genomes for genomic islands and HGT events, distinguishing HGT from incomplete lineage sorting / differential gene loss / hybridization, mapping donor lineages via phylogenetic placement, separating eukaryotic HGT from contamination, ruling out gBGC as a false signal, or quantifying DTL rates with ALE/GeneRax on bacterial trees.
tool_type: mixed
primary_tool: HGTector
---

## Version Compatibility

Reference examples tested with: HGTector 2.0b3+, AvP 1.0.4+, HGTphyloDetect 1.0+, ALE 1.0+ (ssolo/ALE github), GeneRax 2.1.3+, AleRax 1.2.0+ (Morel 2024), RANGER-DTL 2.0+, IslandViewer 4 (web), mobileOG-db 1.0+, MetaCHIP 1.10+, IQ-TREE 2.3.6+, BioPython 1.84+, DIAMOND 2.1.10+. Open Tree of Life and NCBI Taxonomy reference databases updated 2024-Q3 minimum for HGTector/AvP.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show hgtector` then `hgtector search --help`
- CLI: `ALEml_undated --help`, `generax --help`, `alerax --help`
- DB: `hgtector database --check` for taxonomy version

If code throws `Taxonomy ID not found`, `database version mismatch`, or `KeyError` on NCBI taxids, refresh the local taxonomy dump (NCBI updates monthly). ALE/GeneRax expect newick gene trees with bootstraps; AleRax expects gene-tree distributions (uniform bootstrap samples or UFBoot trees).

# Horizontal Gene Transfer Detection

**"Are these genes horizontally acquired, and from where?"** -> HGT signal lives in three orthogonal signal classes: composition (recent transfers carry donor codon usage; erodes by Lawrence-Ochman 1998 amelioration in ~50-200 Myr), phylogeny (gene tree nests within distant clade), and phyletic distribution (patchy taxonomic presence). No single class proves HGT; **claims require concordance across at least two classes** plus mandatory exclusion of contamination and differential gene loss (DGL). The most consequential failure mode in eukaryotic HGT detection is contamination passing all three classes silently (Boothby 2015 tardigrade retraction; Crisp 2015 human "145 HGTs" refuted by Salzberg 2017 GB 18:85).

- Python: `hgtector search` -> `hgtector analyze` for BLAST-distribution screen
- Python: `AvP` (Koutsovoulos 2022 PLoS Comp Biol 18:e1010686) for eukaryotic phylogenetic HGT with automated tree workflow
- CLI: `ALEml_undated` (Szöllősi 2013 Syst Biol 62:901), `generax` (Morel 2020 MBE 37:2763), `alerax` (Morel 2024 Bioinformatics 40:btae162) for prokaryote DTL reconciliation
- Web: IslandViewer 4 (Bertelli 2017 NAR 45:W30) for bacterial genomic islands
- CLI: `metachip` (Song 2019 Microbiome 7:36) for metagenomic HGT inference

## Algorithmic Taxonomy

| Method class | Tool | Signal | Strength | Fails when |
|--------------|------|--------|----------|------------|
| Composition (parametric) | SIGI-HMM (Waack 2006 BMC Bioinf 7:142) | HMM on codon-usage anomaly | Recent transfers (<50 Myr); single-locus resolution | Amelioration eroded composition; native composition heterogeneous (Streptomyces, Borrelia) |
| Composition (parametric) | AlienHunter (Vernikos & Parkhill 2006 Bioinformatics 22:2196) | Variable-window tetranucleotide IVOM | Recent island detection; window-size aware | Old transfers; high-variance native composition |
| Composition (parametric) | IslandPath-DIMOB (Bertelli 2017) | Dinucleotide bias + mobility genes | Combines composition + mobile-element signature | Mobile-element-free transfers |
| Composition aggregator | IslandViewer 4 web (Bertelli 2017 NAR 45:W30) | Consensus of IslandPath / SIGI / IslandPick | Best single-genome bacterial screen; curated benchmark | Recent radiations; close-relative donors |
| BLAST-distribution | HGTector v2 (Zhu 2014 BMC Genomics 15:717) | Close vs distal BLAST hit ratio against full taxonomy | No tree required; scales to thousands of genomes | Close-relative donors (similar lineage hits); incomplete taxonomic sampling |
| BLAST-distribution | DarkHorse (Podell & Gaasterland 2007 GB 8:R16) | Lineage Probability Index (LPI) | Quantifies taxonomic prior of best hits | Same lineage-coverage caveat as HGTector |
| BLAST-distribution | Alien Index (Gladyshev 2008 Science 320:1210; recent: AI tools) | Best hit metazoan vs non-metazoan score | Standard for eukaryote HGT screen | False positives from rapid evolution / contamination |
| Phylogenetic incongruence | AvP (Koutsovoulos 2022) | Automated tree-building + AU test on candidate HGTs | End-to-end eukaryote pipeline; orthogroup-aware | Poor taxon sampling at putative donor lineage |
| Phylogenetic incongruence | HGTphyloDetect (Yuan 2023 Brief Bioinform 24:bbad035) | Web-friendly; Bayesian incongruence test | Modern eukaryote-friendly | Computationally heavy at genome scale |
| Probabilistic DTL reconciliation | ALE (Szöllősi 2013) | Amalgamated likelihood over gene-tree sample | Bayesian-posterior over D/T/L events; explicit donor inference | Requires gene tree distribution (bootstrap or UFBoot) |
| Probabilistic DTL reconciliation | GeneRax (Morel 2020) | ML reconciliation; species-tree-aware | Faster than ALE; refines noisy gene trees | Less uncertainty quantification than ALE |
| Probabilistic DTL reconciliation | AleRax (Morel 2024) | Co-estimates gene + species trees + DTL rates | Gold standard 2024; corrects gene-tree-error feedback | Computationally heavy; needs >= 20 species |
| Parsimony reconciliation | RANGER-DTL 2.0 (Bansal 2018 Bioinformatics 34:3214) | Min-cost DTL parsimony | Fast; deterministic; many trees | No likelihood; sensitive to event-cost choice |
| Metagenome HGT | MetaCHIP (Song 2019 Microbiome 7:36) | BLAST + phylogeny on MAG-pairs | Designed for metagenome-assembled genomes | Requires high-quality MAGs (>= 90% complete) |
| HGT-aware species tree | ASTRAL-Pro2 (Zhang & Mirarab 2022 Bioinformatics 38:4949) | Quartet-coalescent with paralog handling | Robust to HGT-inflated gene-tree discord | Still assumes ILS-coalescent; not explicit HGT model |

Methodology evolves; verify the current AleRax / AvP documentation before committing to a single approach. ALE/GeneRax/AleRax have largely superseded older parsimony reconciliation for bacterial phylogenomics.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Single bacterial genome, recent HGT screen | IslandViewer 4 (web) + HGTector v2 | Composition + BLAST distribution; standard 1-genome workflow |
| 5-200 bacterial genomes, ancient HGT focus | ALE or AleRax on orthogroup trees | Probabilistic DTL; recovers ancient transfers obscured by amelioration |
| 200+ bacterial genomes, phylogenomic HGT rates | GeneRax (faster) -> ALE (validation on top candidates) | Scales; ALE only for the candidates needing posterior support |
| Eukaryote genome, suspected HGT | **Contamination check first** (see below); then AvP or HGTphyloDetect | Eukaryote HGT field is dominated by contamination false positives |
| Metagenome / MAG analysis | MetaCHIP (Song 2019) | Designed for MAGs; tolerates fragmentation |
| Gene-family-level HGT rate inference | ALE/AleRax with site-rate variation | Posterior over D/T/L rates; quantitative inference |
| Donor inference (which lineage was the source) | ALE branch-wise D/T/L map; AvP donor placement | Both attach donor branch posterior |
| Plant-plant HGT (parasitic-host) | AvP with broader plant taxa; manual gene-tree inspection | Standard methods often miss plant-plant transfers |
| Putative HGT correlated with antibiotic resistance | IslandViewer + AMRFinderPlus + mobileOG-db | Combine HGT detection with mobile-element + resistance annotation |
| Phage-mediated transfer screen | PHASTER / PhageBoost + IslandViewer | Phage detection orthogonal to general HGT |
| Endosymbiotic gene transfer (organelle -> nucleus) | Custom: BLAST nuclear proteome against mitochondrial/plastid proteome; tree per hit | Standard HGT tools miss EGT context; expect ~5-15% of nuclear plant proteome is plastid-derived |
| Putative HGT shows incongruence but no other evidence | Test against ILS, DGL, hybridization | Phylogenetic discordance has biological alternatives (Maddison 1997 Syst Biol 46:523) |

## Per-Tool Failure Modes

### Contamination masquerading as eukaryotic HGT (THE critical failure)

**Trigger:** Eukaryote genome assembly, especially from non-axenic culture, microbiome-associated organism, or low-coverage shotgun.

**Mechanism:** Bacterial DNA in the sample is assembled as contigs separate from the eukaryote nuclear contigs but is reported as part of the assembly. Genes on contaminant contigs phylogenetically nest within bacteria, compositionally differ from the eukaryote, and have patchy phyletic distribution -- triggering all three HGT signal classes simultaneously.

**Symptom:** "HGT" genes are concentrated on short, low-coverage contigs; have GC% dramatically different from the bulk genome; show codon usage indistinguishable from bacteria; cluster on contigs lacking eukaryotic gene order; the genome assembly's BUSCO completeness anomaly indicates contamination (e.g. tardigrade Hypsibius dujardini: Boothby 2015 PNAS 112:15976 -> Koutsovoulos 2016 PNAS 113:5053 retraction).

**Fix:** MANDATORY contamination filter before any eukaryotic HGT analysis. Use BlobTools2 (Challis 2020 G3) coverage-vs-GC visualization to identify contaminant blobs; Kraken2 (Wood 2019 GB 20:257) or Conterminator (Steinegger 2020 GB 21:115) to taxonomically classify contigs; FCS-GX (Astashyn 2024 GB 25:60) is the NCBI tool now required for GenBank submission. Apply these BEFORE running AvP / HGTector / Alien Index. After cleaning, re-screen. Crisp 2015 Genome Biol 16:50 "145 HGT in humans" was reduced to ~17 after Salzberg 2017 GB 18:85 contamination-aware reanalysis.

### Amelioration eroding compositional signal

**Trigger:** Compositional-only HGT calls on bacterial genomes diverged > 50 Myr from donor.

**Mechanism:** After transfer, point mutations occur under host mutation bias (e.g. AT bias in obligate symbionts), driving codon usage and GC% toward host values (Lawrence & Ochman 1997 J Mol Evol 44:383; 1998 PNAS 95:9413). Half-life of compositional signal in bacteria is roughly 50-200 Myr depending on mutation rate and selection. Ancient transfers are compositionally indistinguishable from native genes.

**Symptom:** Phylogenetic methods detect transfer but composition-based tools (SIGI, IslandPath) miss it; GC anomaly is weak (|z| < 2); codon adaptation index resembles native genes.

**Fix:** For HGT older than ~50 Myr, rely on phylogenetic methods (ALE / AvP); composition is informative only as supporting evidence for recent transfers. Report transfer age estimate from ALE branch posterior alongside composition.

### Differential gene loss mimicking HGT

**Trigger:** Gene present in distantly related taxa but absent in immediate sister lineages.

**Mechanism:** An ancestral gene is independently lost in multiple intermediate lineages; the surviving taxa appear "incompatible" with the species tree. Phylogenetic incongruence and patchy taxonomic distribution are identical to HGT signatures (Maddison 1997 Syst Biol 46:523).

**Symptom:** Gene-tree topology actually matches species-tree topology of the surviving taxa, just with intermediate taxa missing. Composition matches the recipient genome (no HGT amelioration to explain). ALE infers high loss rate at branches and may favor a "transfer + loss" or "duplication + losses" event class depending on cost weights.

**Fix:** Inspect ALE event-class posteriors (D vs T vs L); when loss rate inference is non-negligible at sister branches, prefer the loss-explanation. Check independently sequenced relatives if available. Quantify Dollo-parsimony loss vs ML-DTL event likelihoods. RANGER-DTL with loss cost = 1, transfer cost = 3 favors loss; rerun with transfer cost = 1.5 to see sensitivity.

### Hybridization / reticulate evolution

**Trigger:** Sister species suspected to have hybridized or share recent introgression; rapid radiations.

**Mechanism:** Gene flow between closely related lineages produces gene-tree species-tree discordance indistinguishable from HGT at short divergence times. Phylogenetic networks (not trees) are required.

**Symptom:** Multiple genes show identical pattern of incongruence (same direction, same sister); D-statistic (see [[introgression-detection]]) significantly different from zero between candidate hybrid and its inferred sister.

**Fix:** Run ABBA-BABA / Dsuite for the candidate (see [[introgression-detection]]); if introgression signal is uniform across the genome, prefer hybridization explanation. ALE on bacterial trees can handle this via transfer-with-replacement; for eukaryotes, phylogenetic networks (e.g. PhyloNetworks / SNaQ via Solis-Lemus & Ane 2016 PLoS Genet 12:e1005896, or NakhlehLab PhyloNet) properly model reticulation.

### gBGC-driven AT->GC substitution bias

**Trigger:** Mammalian / vertebrate gene with apparent unusual nucleotide composition triggering compositional HGT call.

**Mechanism:** GC-biased gene conversion (gBGC) in regions of high recombination drives Weak->Strong (AT->GC) fixation independently of selection (Galtier & Duret 2007 Trends Genet 23:273). Sub-telomeric and recombination-hot regions show elevated GC indistinguishable from a GC-rich HGT donor at the composition level.

**Symptom:** "HGT" candidates cluster in sub-telomeric regions or high-recombination zones; W->S substitution bias is elevated; phylogenetic placement is consistent with host species; selection scans (BUSTED) show no signal.

**Fix:** Test for gBGC explicitly via W->S vs S->W substitution ratios (Galtier 2025 Genetics 230:iyaf111; Capra 2013 PLoS Genet 9:e1003684); require non-zero phylogenetic incongruence in addition to composition; exclude regions with recombination rate > 90th percentile.

### Close-relative donor making BLAST methods fail

**Trigger:** Transfer between species in the same genus or family (genus-level transfer).

**Mechanism:** HGTector / DarkHorse / Alien Index compare close vs distal BLAST hits; when donor is closely related, hits cluster in "close" and the method reports no anomaly.

**Symptom:** HGTector hgt_score near zero, but the gene shows clear phylogenetic placement within a sister taxon rather than expected vertical inheritance.

**Fix:** Use phylogenetic methods (ALE on orthogroup trees) for genus-level HGT; BLAST-distribution methods are designed for inter-phylum transfer. Restrict HGTector taxonomic comparison level (`--rank order` or higher).

### Poor sampling of donor lineage

**Trigger:** Putative HGT phylogenetically places "near" a poorly sampled taxon (e.g. Archaea or candidate phyla).

**Mechanism:** Long branches in the donor clade attract the query gene by LBA (Felsenstein 1978 Syst Zool 27:401); poor sampling means real intermediate relatives are absent, and the gene appears nested within the wrong clade.

**Symptom:** AU/SH topology test rejects alternative placements only marginally (p ~ 0.05); adding any newly available genome from the candidate donor lineage changes placement substantially.

**Fix:** Use site-heterogeneous models (CAT-PMSF, PhyloBayes-MPI) to mitigate LBA at deep nodes; add taxa from undersampled lineages where possible; report donor inference with explicit support-level caveats. Cross-check with ALE branch posterior, which integrates over alternative donor lineages probabilistically.

### Tetranucleotide z-score thresholding in heterogeneous genomes

**Trigger:** Applying generic |z| > 2 cutoff to genomes with naturally high composition variance (Streptomyces, Burkholderia, Halophiles).

**Mechanism:** rRNA operons, prophage remnants, and CRISPR arrays have native composition anomalies; generic z-score thresholds flag them as HGT.

**Symptom:** "HGT" calls concentrate in rRNA-flanking regions, ribosomal protein operons, or annotated prophages.

**Fix:** Mask known native-anomaly features (rRNA, tRNA, CRISPR) before composition analysis; use genome-specific thresholds calibrated against confirmed native genes; raise |z| threshold to 3 for high-variance genomes.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| HGT call requires concordance | >= 2 of 3 signal classes (composition, phylogeny, distribution) | Standard convention (Ravenhall 2015 PLoS Comp Biol 11:e1004095) |
| HGTector hgt_score | > 0.5 moderate; > 1.0 strong | HGTector documentation; Zhu 2014 |
| Alien Index | AI > 45 (P_metazoan / P_non-metazoan) | Gladyshev 2008 Science 320:1210; modern eukaryote default |
| GC z-score generic threshold | |z| > 2 suggestive; > 3 strong | But raise to >= 3 for high-variance genomes |
| Genomic island minimum size | >= 5 contiguous genes / >= 8 kb | IslandViewer 4 default |
| ALE D/T/L event posterior | branch-wise event > 0.5 to call transfer | ssolo/ALE convention |
| AU test for tree-topology rejection | p < 0.05 to reject vertical inheritance topology | Shimodaira 2002 Syst Biol 51:492 |
| Tree-based HGT minimum sequences | >= 8 taxa per gene tree; >= 4 from putative recipient clade | Below this, donor inference unreliable |
| Required contamination-screen completeness | FCS-GX or BlobTools applied; documented in methods | NCBI GenBank now requires FCS-GX for new submissions |
| ALE undated vs dated | Use undated for ancient comparisons; dated for time-calibrated trees | Szöllősi 2013; AleRax recommends dated when sub-clade times are available |
| Bayesian gene-tree sample size for ALE | >= 100 bootstrap or UFBoot trees | ALE convention |
| Donor lineage minimum sampling | >= 5 genomes from candidate donor order | Below this, donor inference is exploratory |
| Eukaryotic HGT post-cleanup threshold | Re-screen after FCS-GX / BlobTools; expect 10-50x reduction | Boothby->Koutsovoulos 2016 (~80% reduction); Crisp->Salzberg 2017 (~90% reduction) |
| MetaCHIP MAG quality | >= 90% complete, < 5% contamination (CheckM2) | Song 2019 |

## HGTector v2 Workflow

**Goal:** Identify HGT candidates in a prokaryote genome via taxonomic BLAST hit distribution.

**Approach:** Build / download taxonomy database -> `hgtector search` to run DIAMOND against reference proteomes -> `hgtector analyze` to score close-vs-distal taxonomic distribution -> rank candidates by score.

```python
'''HGTector v2 wrapper with quality filters'''

import subprocess
import pandas as pd


def run_hgtector(proteome_faa, out_dir, db_dir, threads=8, rank='order'):
    '''HGTector v2 search + analyze.

    rank: NCBI taxonomic rank below which hits are "close"; "order" is default.
    Use "family" to detect family-level transfers, "phylum" only for ancient inter-phylum HGT.
    '''
    subprocess.run([
        'hgtector', 'search', '-i', proteome_faa, '-o', f'{out_dir}/search',
        '-m', 'diamond', '-d', f'{db_dir}/db.dmnd', '-t', str(threads),
        '--queries-per-batch', '200'
    ], check=True)
    subprocess.run([
        'hgtector', 'analyze', '-i', f'{out_dir}/search', '-o', f'{out_dir}/analyze',
        '-t', f'{db_dir}/taxdump', '--rank-self', rank,
        '--rank-close', rank, '--bandwidth', 'auto'
    ], check=True)
    return f'{out_dir}/analyze'


def parse_hgtector(analyze_dir, min_score=0.5):
    '''HGTector outputs scores.tsv with: gene, hits, close, distal, ...
    distal > close indicates likely foreign origin.
    '''
    df = pd.read_csv(f'{analyze_dir}/scores.tsv', sep='\t')
    df['hgt_score'] = df['distal'] - df['close']
    df['confidence'] = pd.cut(df['hgt_score'], bins=[-1e9, 0.5, 1.0, 1e9],
                              labels=['low', 'moderate', 'strong'])
    return df[df['hgt_score'] > min_score].sort_values('hgt_score', ascending=False)
```

## ALE Probabilistic DTL Reconciliation

**Goal:** Quantify D / T / L events with branch-wise posterior on a bacterial / archaeal species tree.

**Approach:** Infer per-orthogroup gene-tree bootstrap distributions (UFBoot from IQ-TREE) -> observe (encode as `.ale` files) -> run `ALEml_undated` to reconcile against species tree -> extract per-branch transfer count posteriors.

```bash
# 1. Build per-orthogroup UFBoot trees
for og in orthogroups/*.fa; do
    iqtree2 -s $og -m TEST -B 1000 -nt 2 --prefix ${og%.fa}_ufb
done

# 2. Observe and reconcile
for tree in orthogroups/*_ufb.ufboot; do
    ALEobserve $tree
    ALEml_undated species_tree.nwk ${tree}.ale separators="|" sample=100
done

# 3. Aggregate D/T/L counts per species-tree branch
python aggregate_ale.py orthogroups/*_uTs > dtl_branchwise.tsv
```

```python
'''Aggregate ALE outputs into branch-wise event tables.'''
import re
from collections import defaultdict
import pandas as pd

def parse_ale_uts(uts_file):
    '''ALE _uTs format: per-branch D/T/L counts (mean over sampled reconciliations).'''
    branch_events = {}
    with open(uts_file) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 5:
                branch_id = parts[0]
                # columns: branch  duplications  transfers  losses  originations
                branch_events[branch_id] = {
                    'duplications': float(parts[1]),
                    'transfers': float(parts[2]),
                    'losses': float(parts[3]),
                    'originations': float(parts[4])
                }
    return branch_events
```

Bacterial datasets characteristically show transfers >> duplications, while opisthokonts (eukaryotes) show the opposite (Szöllősi 2013; Williams 2017 PNAS 114:E4602). ALE-rooting of species trees (Williams 2017) is now the standard rooting approach for deep bacterial phylogenomics.

## AvP Eukaryotic HGT Phylogenetic Workflow

**Goal:** Detect HGT in eukaryotic proteomes with automated phylogenetic-tree validation.

**Approach:** `avp prepare` builds taxonomic database from candidate hits -> `avp detect` builds gene trees per candidate -> `avp classify` runs AU/SH tests on alternative placements.

```bash
avp prepare -t taxonomy.tsv -d nr.fasta -o database
avp detect -i query_proteome.fasta -d database -o avp_results -threads 8
avp classify -i avp_results -d database --au-test --bootstrap 1000
```

The output `avp_results/classification.tsv` reports per-gene: HGT / putative-HGT / unknown / non-HGT, donor lineage assignment, gene-tree support, and AU-test p-value. **Apply only after FCS-GX or BlobTools contamination filter.**

## IslandViewer 4 + mobileOG-db Integration

For single-genome bacterial analysis, the combined workflow is: IslandViewer 4 (web; submit assembly, returns islands via IslandPath-DIMOB / SIGI-HMM / IslandPick consensus) -> annotate genes with mobileOG-db (Brown 2022 Appl Environ Microbiol 88:e0099122) for mobile element families -> cross-check candidate islands for integrase / transposase / phage signatures -> if AMR is annotated, cross-reference with AMRFinderPlus or CARD.

The five most reliable island indicators in order of strength: (1) integrase / transposase / phage protein co-located; (2) tRNA gene at one boundary (integration hotspot); (3) GC% deviation > 2 SD from genome mean; (4) absence in close relatives; (5) DNA recognition sites for known integrases (e.g. attL/attR).

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| HGTector positive, ALE rejects transfer at any branch | Sequencing artifact, BLAST distribution skewed by uneven taxonomic sampling | Trust ALE (phylogenetic, posterior-based); investigate HGTector's "distal" hit composition |
| Composition positive, phylogeny negative | Native compositional outlier (rRNA region, prophage native to lineage) | Mask known features; treat as native unless mobile-element signature corroborates |
| Composition negative, phylogeny positive (ancient HGT) | Amelioration completed | Trust phylogeny; report transfer age from ALE branch posterior |
| ALE high transfer posterior, RANGER-DTL low | Cost-weight sensitivity in RANGER | Trust ALE (likelihood-based); RANGER costs are user choice |
| Multiple methods positive, but candidate is on low-coverage contig | Contamination | Apply FCS-GX; re-screen post-cleanup |
| AvP positive in eukaryote, but candidate is single-exon, no introns | Contamination (bacterial-origin) or recent bacterial gene without intron acquisition | Inspect contig genome architecture; eukaryote-genomic flanking required |
| HGT signal in tardigrade / human / sponge / Nematoda | Suspect contamination | Re-examine assembly quality; FCS-GX mandatory (Boothby 2015 / Crisp 2015 lessons) |
| Putative HGT in a clade with known hybridization | Reticulate evolution, not HGT | Switch to phylogenetic networks (PhyloNet); D-stat across the clade |
| Donor placement weak (AU p ~ 0.05) | Poor donor sampling or LBA | Add available genomes from candidate donor; switch to CAT-PMSF model |

**Operational rule for publication:** Concordant evidence from 2+ signal classes + contamination-cleared assembly + donor lineage with > 5 sampled relatives + ALE/AvP phylogenetic support with AU p < 0.05 + DGL ruled out by sister-lineage inspection. Single-method (especially composition-only) claims should be downgraded to "candidate" status.

## Cohort Gotchas

- **Genome MAGs from metagenomes:** completeness and contamination matter more than for isolates. Apply CheckM2 first; HGT detection on incomplete MAGs is unreliable.
- **Streptomyces, Borrelia, Mollicutes:** genomes with extreme GC variance or unusual replication; composition methods produce many false positives. Raise z-threshold to 3 and rely more on phylogenetic methods.
- **Endosymbionts (Buchnera, Wolbachia):** strong AT bias from genome reduction; gene from external donor with normal composition appears as HGT but may be ancestral. Use AleRax with reduced-effective-population-size aware priors.
- **Plant genomes:** plastid-to-nucleus EGT is widespread (~5-18% of nuclear plant proteome plastid-derived per Martin 2002 PNAS 99:12246); standard HGT tools detect these but they are not "horizontal" in the prokaryote sense.
- **Animal microbiome-associated organisms:** consortium contamination is common in non-axenic cultures; pre-2010 published "HGTs" in such organisms require FCS-GX re-screening.
- **Single-cell amplified genomes (SAGs):** chimeric MDA artifacts mimic HGT; require co-assembly validation before HGT calls.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Contamination?" | FCS-GX applied; BlobTools coverage-vs-GC plot showed no contaminant blobs; HGT candidates re-screened post-cleanup |
| "Why is this transfer, not differential gene loss?" | ALE branch-wise event posteriors favor transfer (posterior > 0.5); sister taxa show no trace of gene; loss-only model fits AIC-worse |
| "Why not just composition?" | Amelioration window is ~50-200 Myr (Lawrence-Ochman 1998); for older transfers, phylogenetic methods are required; we report composition only as supporting evidence |
| "How was donor inferred?" | ALE branch transfer posterior > X on donor branch; AU test rejects alternative placements (p < 0.05); >= 5 genomes sampled from donor order |
| "ILS not explanation?" | Tested via Dsuite D-statistic (introgression-detection); not significant; gene-tree incongruence cannot be explained by recent ILS at the divergence depth involved |
| "gBGC?" | W->S substitution bias not elevated; gene is not in high-recombination region; selection scan (positive-selection) is null |
| "Hybridization?" | D-statistic and f4 between candidate donor / recipient pair tested; non-significant; or, if significant, reported as introgression not HGT |
| "Mobile element signature?" | Integrase / transposase co-located OR adjacent tRNA gene OR attL/attR; or, if none, explicitly flagged as composition+phylogeny-only evidence |
| "Why ALE over RANGER?" | ALE provides posterior probabilities; RANGER is parsimony with arbitrary cost weights |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| HGTector "Taxonomy ID X not found" | Outdated `taxdump` | `hgtector database --update`; ensure taxdump is < 6 months old |
| AvP runs but no candidates | Strict default cutoffs | Relax `--alien-index-cutoff` (default 45 -> try 30); inspect intermediate files |
| ALE rejects all gene trees | Gene-tree taxon names don't match species tree | Use orthogroup-consistent species labels; check separator ' \t' not '_' |
| GeneRax MPI hangs | One CPU per tree; thread-count too high | Reduce `--threads` to match cores actually allocated by scheduler |
| AleRax "no parameter convergence" | Too few gene trees or species | Need >= 20 species and >= 100 gene families; convergence is hard for small datasets |
| FCS-GX returns 0 contaminants when external evidence suggests contamination | Reference DB outdated or strain-specific | Update FCS-GX reference; cross-check with BlobTools and Kraken2 |
| IslandViewer returns nothing | Genome < 1 Mb or fragmented | Use stand-alone IslandPath-DIMOB or SIGI on full assembly; manual annotation of small genomes |
| Codon usage CAI all near 1.0 | Reference codon usage from same gene set | Use external highly-expressed-gene reference (ribosomal genes) for CAI scaling |
| RANGER-DTL "transfer cost too low gives transfers everywhere" | Default cost balance | Run with cost-sensitivity sweep (T cost 1-5, L cost 1, D cost 1-2); report robust events only |
| Eukaryotic HGT call from contig < 50 kb | Likely contaminant | Filter for HGT calls on contigs > 100 kb with >= 5 native eukaryote genes flanking |

## Tool Installation Notes

```bash
# HGTector v2
pip install hgtector
hgtector database --output ~/hgtector_db --reference refseq --rank species

# AvP
git clone https://github.com/GDKO/AvP && cd AvP && pip install .
# Requires DIAMOND, IQ-TREE, mafft on PATH

# ALE + GeneRax + AleRax
conda install -c bioconda ale generax
# AleRax: build from source per https://github.com/BenoitMorel/AleRax
git clone --recursive https://github.com/BenoitMorel/AleRax && cd AleRax && ./install.sh

# Contamination tools (REQUIRED for eukaryote HGT)
conda install -c bioconda fcs-gx blobtoolkit kraken2
pip install conterminator

# Mobile elements
git clone https://github.com/clb21565/mobileOG-db
conda install -c bioconda amrfinder
```

NCBI now mandates FCS-GX screening for new eukaryote genome submissions to GenBank (Astashyn 2024 GB 25:60); apply it to any pre-2023 published eukaryote assembly before HGT analysis.

## References

- Lawrence JG & Ochman H 1997 J Mol Evol 44:383 (amelioration concept)
- Lawrence JG & Ochman H 1998 PNAS 95:9413 (amelioration quantification)
- Maddison WP 1997 Syst Biol 46:523 (gene tree species tree discordance causes)
- Felsenstein J 1978 Syst Zool 27:401 (LBA)
- Galtier N & Duret L 2007 Trends Genet 23:273 (gBGC review)
- Galtier N 2025 Genetics 230:iyaf111 (gBGC selection model)
- Vernikos GS & Parkhill J 2006 Bioinformatics 22:2196 (AlienHunter)
- Waack S et al 2006 BMC Bioinf 7:142 (SIGI-HMM)
- Podell S & Gaasterland T 2007 Genome Biol 8:R16 (DarkHorse)
- Gladyshev EA et al 2008 Science 320:1210 (Alien Index; bdelloid rotifer HGT)
- Szöllősi GJ et al 2013 Syst Biol 62:901 (ALE undated)
- Bansal MS et al 2018 Bioinformatics 34:3214 (RANGER-DTL 2.0)
- Bertelli C et al 2017 NAR 45:W30 (IslandViewer 4)
- Williams TA et al 2017 PNAS 114:E4602 (ALE rooting of the archaeal tree of life)
- Boothby TC et al 2015 PNAS 112:15976 (tardigrade HGT claim)
- Koutsovoulos G et al 2016 PNAS 113:5053 (tardigrade contamination refutation)
- Crisp A et al 2015 Genome Biol 16:50 (human HGT claim)
- Salzberg SL 2017 Genome Biol 18:85 (human HGT contamination refutation)
- Zhu Q, Kosoy M & Dittmar K 2014 BMC Genomics 15:717 (HGTector)
- Koutsovoulos G et al 2022 PLoS Comp Biol 18:e1010686 (AvP)
- Yuan L et al 2023 Brief Bioinform 24:bbad035 (HGTphyloDetect)
- Morel B et al 2020 MBE 37:2763 (GeneRax)
- Morel B et al 2024 Bioinformatics 40:btae162 (AleRax co-estimation)
- Brown CL et al 2022 Appl Environ Microbiol 88:e0099122 (mobileOG-db)
- Song W et al 2019 Microbiome 7:36 (MetaCHIP)
- Astashyn A et al 2024 GB 25:60 (FCS-GX)
- Wood DE et al 2019 GB 20:257 (Kraken2)
- Steinegger M et al 2020 GB 21:115 (Conterminator)
- Challis R et al 2020 G3 10:1361 (BlobTools2)
- Ravenhall M et al 2015 PLoS Comp Biol 11:e1004095 (HGT inference review)
- Martin W et al 2002 PNAS 99:12246 (plastid EGT)
- Capra JA et al 2013 PLoS Genet 9:e1003684 (gBGC quantification)
- Solis-Lemus C & Ane C 2016 PLoS Genet 12:e1005896 (SNaQ / PhyloNetworks reticulation)
- Shimodaira H 2002 Syst Biol 51:492 (AU test)

## Related Skills

- comparative-genomics/gene-tree-species-tree-reconciliation - Deep DTL methodology (ALE / GeneRax / AleRax)
- comparative-genomics/ortholog-inference - Orthogroup construction for ALE input
- comparative-genomics/introgression-detection - Distinguish HGT from hybridization via D-statistic
- comparative-genomics/genome-distance-and-species-delineation - ANI / dDDH for donor-recipient distance
- comparative-genomics/pangenome-analysis - Bacterial pangenome assembly for accessory-gene HGT analysis
- phylogenetics/modern-tree-inference - Gene-tree inference with UFBoot for ALE
- phylogenetics/bayesian-inference - PhyloBayes-MPI CAT-GTR for deep-node LBA mitigation
- alignment/multiple-alignment - Codon-aware alignment for compositional analyses
- metagenomics/amr-detection - Mobile element / resistance gene context for HGT candidates
- genome-annotation/prokaryotic-annotation - Mobile-element annotation alongside gene annotation
