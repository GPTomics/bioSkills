---
name: bio-microbiome-functional-prediction
description: Predicts community functional POTENTIAL from 16S/ITS amplicon ASVs with PICRUSt2 (or q2-picrust2) by phylogenetic interpolation of reference-genome gene content - EPA-ng placement, gappa, castor hidden-state prediction of KO/EC/Pfam copy number, 16S copy-number normalization, and MinPath MetaCyc/KEGG pathways - gated by the NSTI quality index. Covers why predicted function is taxonomy re-encoded (never measured gene content and never activity), the mandatory NSTI report (--max_nsti 2 silently drops novel ASVs), why accuracy IS reference coverage (gut Spearman ~0.8, soil/marine collapse), the circularity trap, and Tax4Fun2/FAPROTAX/BugBase alternatives. Use when inferring KO/EC/MetaCyc potential from an ASV table, gating on NSTI, or choosing a prediction method. For MEASURED shotgun function see metagenomics/functional-profiling; for enrichment of KO lists see pathway-analysis/go-enrichment; for DA of predicted tables see differential-abundance.
tool_type: cli
primary_tool: PICRUSt2
---

## Version Compatibility

Reference examples tested with: PICRUSt2 2.5+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The bundled REFERENCE is the version that matters. PICRUSt2 ships its own reference tree, alignment, and per-genome trait tables (16S/KO/EC/Pfam/COG/TIGRFAM) for ~20,000 genomes - there is no separate multi-GB download as in HUMAnN. That fixed reference is the entire ceiling on accuracy, and NSTI is computed against THAT reference. The reference is versioned with the PICRUSt2 release; record the release so a prediction is reproducible.

# Functional Prediction with PICRUSt2

**"Predict the functional pathways from my 16S data"** -> Place ASVs on a reference genome tree and report the gene content of their nearest sequenced relatives - because PICRUSt2 never sees a functional gene from the sample, it infers POTENTIAL from who-is-there.
- CLI: `picrust2_pipeline.py -s asv_seqs.fna -i asv_table.biom -o picrust2_out -p 8 --max_nsti 2`

Scope: PREDICTED gene-content potential from amplicon ASVs. MEASURED shotgun function (HUMAnN) -> metagenomics/functional-profiling. ASV table + rep-seqs come from amplicon-processing + taxonomy-assignment. The QIIME2 q2-picrust2 path -> qiime2-workflow. DA of predicted tables -> differential-abundance (compositional, same theory as metagenomics/abundance-estimation). Reading KO/pathway lists -> pathway-analysis/go-enrichment.

## The Single Most Important Modern Insight -- Predicted Function Is Taxonomy Re-Encoded, Not a Measurement

PICRUSt2 never sequences a functional gene from the sample. It places each ASV on a tree of ~20,000 reference genomes, asks "what genes do this ASV's nearest sequenced relatives carry?", and reports that guess as the community's function. Every output number is the gene content of OTHER organisms, weighted by how much of the 16S resembles them. Three corollaries each common misuse violates:

1. **Potential, never activity.** The honest claim is "increased butyrate-production capacity," never "increased butyrate production / upregulated / more active." The inference chain is 16S abundance -> who-is-there -> relatives' genomes -> predicted gene presence -> predicted copy number -> potential. Transcription and flux are further measurement layers away.
2. **Accuracy IS reference coverage.** Predicted-vs-shotgun Spearman tops out around 0.8 in the densely-referenced human gut (Douglas 2020); in soil, marine, sediment, and novel environments the references are sparse, NSTI rises, and the prediction collapses toward a restatement of taxonomy.
3. **Circularity.** Predicted function is a DETERMINISTIC function of the ASV table - same input, same KO table, every time. So a "functional difference" between groups is the TAXONOMIC difference projected through a fixed lookup, not independent corroboration. Predicted-function analysis is hypothesis-generating; a functional conclusion needs shotgun or metatranscriptomics.

Organize the analysis around defending these three (report NSTI, forbid activity verbs, do not double-count taxonomy as a second finding), not around listing flags.

## The Predicted < Measured Ladder

State this explicitly when reporting. Each rung is a real measurement layer above the one below:

predicted potential (PICRUSt2, this skill) < measured gene carriage (HUMAnN / shotgun, metagenomics/functional-profiling) < measured transcription (metatranscriptomics) < measured flux (fluxomics/metabolomics).

PICRUSt2 is the bottom rung. A PICRUSt2 `path_abun_unstrat` table looks identical to a HUMAnN `pathabundance` table - same MetaCyc IDs, same shape, same MinPath logic - but HUMAnN counts reads that actually aligned to genes in the sample, while PICRUSt2 reports relatives' genomes. Never merge or compare the two as interchangeable.

## How PICRUSt2 Builds the Table (the five-stage interpolation)

Every stage is an inference layered on the previous, so error compounds (Douglas 2020 *Nat Biotechnol* 38:685; original Langille 2013 *Nat Biotechnol* 31:814). PICRUSt2's headline advance over PICRUSt1 is that it places arbitrary de-novo ASVs (not closed-reference Greengenes OTUs):

1. **Placement.** Align ASVs to the reference alignment (hmmalign), place into the reference tree with EPA-ng (Barbera 2019 *Syst Biol* 68:365), resolve placements with gappa (Czech 2020 *Bioinformatics* 36:3263). Default `--placement_tool epa-ng` (alternative `sepp`).
2. **Hidden-state prediction (HSP).** For each placed ASV - whose genome was never sequenced - interpolate the copy number of every gene family from neighboring reference genomes using castor (Louca & Doebeli 2018 *Bioinformatics* 34:1053). Default `--hsp_method mp` (maximum parsimony; `pic` is faster but `mp` is the recommended default).
3. **16S copy-number normalization.** Divide each ASV's abundance by its PREDICTED 16S copy number (a genome with seven 16S copies is over-counted 7x in an amplicon survey). The correction is itself an HSP output, not a measurement.
4. **Metagenome inference.** Per-sample gene abundance = sum over ASVs of (ASV abundance / 16S copies) x (predicted gene copy number). Emits KO/EC/Pfam tables.
5. **Pathway inference.** Map genes to MetaCyc reactions and call pathways with MinPath (Ye & Doak 2009 *PLoS Comput Biol* 5:e1000465). Pathway COVERAGE (`--coverage`, presence confidence) and pathway ABUNDANCE are different questions - do not conflate.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| PICRUSt2 | Douglas 2020 *Nat Biotechnol* 38:685 | EPA-ng placement + castor HSP into ~20k-genome tree -> KO/EC/Pfam, MetaCyc | de-novo ASVs, broad function; strongest in human gut; the default |
| Tax4Fun2 | Wemheuer 2020 *Environ Microbiome* 15:11 | nearest-BLAST to reference 16S + habitat-specific reference + functional-redundancy index | habitat-specific reference available; want a functional-redundancy metric |
| FAPROTAX | Louca 2016 *Science* 353:1272 | curated literature lookup: taxon -> biogeochemical group (nitrification, methanogenesis, sulfate reduction) | ENVIRONMENTAL/biogeochemistry, cultured-taxon-dominated marine/soil; NOT gene-content prediction |
| BugBase | Ward 2017 *bioRxiv* 133462 | organism-level PHENOTYPE prediction (PICRUSt-style gene content) | coarse community phenotypes (aerobic/anaerobic, Gram, oxidative stress); PREPRINT-only, never journal-published |
| PanFP | Jun 2015 *BMC Res Notes* 8:479 | per-lineage pangenome profile weighted by OTU abundance | lineage-level functional summary; no de-novo placement |
| Piphillin | Iwai 2016 *PLoS ONE* 11:e0166104 | direct nearest-neighbor BLAST of ASVs to genome DB (no tree, no HSP) | avoids the phylogenetic-interpolation assumption |

FAPROTAX is conceptually different: a curated taxon-to-function lookup answering "is this community nitrifying?", saying nothing about uncharacterized taxa (they map to nothing). Use FAPROTAX for biogeochemical-cycle questions, PICRUSt2 for a predicted KO/pathway profile - different questions.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Human gut / well-referenced host, broad KO/pathway profile | PICRUSt2 + report NSTI | dense references -> low NSTI -> Spearman ~0.8 vs shotgun; credible broad strokes |
| Soil / marine / sediment / novel environment, biogeochemical question | FAPROTAX (broad groups) | reference tree is gut/host-biased; high NSTI makes PICRUSt2 KOs mostly extrapolation |
| Need REAL (measured) functional gene content | -> metagenomics/functional-profiling | only shotgun reads measure genes; PICRUSt2 predicts, HUMAnN measures |
| Coarse community phenotype (aerobe/anaerobe, Gram, pathogenic potential) | BugBase (flag: preprint-only) | organism-level phenotype, not pathways; cite its unpublished status |
| Already in QIIME2 with .qza artifacts | -> qiime2-workflow (q2-picrust2 plugin) | same engine; FeatureTable in, FeatureTable out, into qiime diversity/composition |
| DA between groups on predicted table | -> differential-abundance, run >=2 CoDA tools | compositional + prediction error; frame as hypothesis-generating (circularity) |
| Reading the resulting KO/pathway list for biology | -> pathway-analysis/go-enrichment | predicted abundance is summed gene content, not an over-representation test |

## Run the Pipeline

```bash
# Full pipeline: place -> HSP -> 16S-normalize -> metagenome -> MetaCyc pathways
picrust2_pipeline.py \
    -s asv_seqs.fna \              # representative ASV sequences (FASTA)
    -i asv_table.biom \            # ASV abundance table (BIOM or TSV; samples as columns)
    -o picrust2_out \
    -p 8 \
    --hsp_method mp \              # maximum parsimony (recommended default); pic is faster but not recommended
    --max_nsti 2 \                 # ASVs ABOVE 2 are DROPPED before inference; report how many (see below)
    --verbose
# Key outputs (gzipped):
#   KO_metagenome_out/pred_metagenome_unstrat.tsv.gz   KEGG ortholog abundances
#   EC_metagenome_out/pred_metagenome_unstrat.tsv.gz   EC-number abundances
#   pathways_out/path_abun_unstrat.tsv.gz              MetaCyc pathway abundances
#   marker_predicted_and_nsti.tsv.gz                   per-ASV 16S copies + metadata_NSTI (the quality file)
```

`--stratified` additionally emits per-ASV contribution tables (large, much slower); `--per_sequence_contrib` is only meaningful with it. `--coverage` adds pathway coverage (a different question from abundance). The unrolled per-step scripts are `place_seqs.py` -> `hsp.py -i {16S,KO,EC} -m mp [-n]` -> `metagenome_pipeline.py --max_nsti 2` -> `pathway_pipeline.py`; `--max_nsti` filtering and 16S normalization happen in `metagenome_pipeline.py`. `add_descriptions.py -m METACYC` attaches human-readable names.

## Report NSTI (mandatory)

**Goal:** Quantify how much of the prediction is extrapolation and how much of the sampled community the NSTI gate discarded, so the result is interpretable.

**Approach:** Read `marker_predicted_and_nsti.tsv.gz`, summarize the `metadata_NSTI` distribution, and report the number of ASVs AND the fraction of READS dropped at `--max_nsti 2` (a study that loses 40% of reads predicted function for a different community than it sampled).

```python
import pandas as pd

nsti = pd.read_csv('picrust2_out/marker_predicted_and_nsti.tsv.gz', sep='\t')   # cols: sequence, metadata_NSTI
asv_counts = pd.read_csv('asv_table.tsv', sep='\t', index_col=0)                # ASVs x samples
nsti = nsti.set_index('sequence')
reads_per_asv = asv_counts.sum(axis=1)

max_nsti = 2.0   # PICRUSt2 default; ASVs above this are dropped before metagenome inference
dropped = nsti.index[nsti['metadata_NSTI'] > max_nsti]
reads_dropped_frac = reads_per_asv.reindex(dropped).sum() / reads_per_asv.sum()
print(f'mean NSTI {nsti.metadata_NSTI.mean():.3f}  median {nsti.metadata_NSTI.median():.3f}')
print(f'ASVs dropped at NSTI>{max_nsti}: {len(dropped)}/{len(nsti)}  reads dropped: {reads_dropped_frac:.1%}')
```

## Per-Method Failure Modes

### Claiming activity or expression
**Trigger:** reporting "increased butyrate production" / "upregulated" / "more metabolically active." **Mechanism:** PICRUSt2 measured no genes and no transcripts - only inferred gene presence from relatives' genomes. **Symptom:** a results sentence with an activity verb on a predicted pathway. **Fix:** restrict every claim to "potential" / "predicted capacity"; for activity, cite metatranscriptomics, not this skill.

### NSTI ignored or under-reported
**Trigger:** accepting the default `--max_nsti 2` filter without reporting the distribution or the dropped fraction. **Mechanism:** the filter silently deletes the most novel/under-referenced ASVs - exactly the organisms an environmental study cares about. **Symptom:** a predicted-function result with no NSTI numbers in the methods. **Fix:** report mean/median NSTI, the distribution, and the ASV AND read fraction dropped (the helper above); treat high mean NSTI as a red flag that the result is mostly extrapolation.

### Wrong environment (reference coverage too sparse)
**Trigger:** running PICRUSt2 on soil/marine/sediment/plant/novel hosts and reporting fine-grained KO differences. **Mechanism:** the ~20k-genome reference tree is gut/host-biased; sparse references mean high NSTI and predictions interpolated from distant relatives. **Symptom:** high mean NSTI yet confident KO/pathway tables. **Fix:** report the environment and NSTI; prefer FAPROTAX for the broad biogeochemical question, or do real shotgun. "Relatively better than other predictors" (per the paper) is not "trustworthy in absolute terms."

### Predicted treated as measured
**Trigger:** describing `path_abun_unstrat` as "the pathways present in the community" or merging it with a HUMAnN table. **Mechanism:** the two tables share IDs and shape but PICRUSt2 reports relatives' genomes, HUMAnN reports reads that aligned to genes in the sample. **Symptom:** predicted and measured tables combined or compared as one. **Fix:** label predictions as potential; never merge with shotgun; keep coverage and abundance as separate questions.

### Circularity (taxonomy double-counted as a second finding)
**Trigger:** reporting "groups differed taxonomically AND functionally" as two lines of evidence. **Mechanism:** predicted function is a deterministic function of the ASV table, so the functional difference IS the taxonomic difference re-encoded. **Symptom:** a predicted-function DA result presented as orthogonal corroboration of a taxonomic result. **Fix:** present predicted function as a hypothesis-generating summary of the taxonomic signal; for orthogonal functional evidence use shotgun/metatranscriptomics.

### DA without compositional correction
**Trigger:** uncorrected Wilcoxon/t-test on relative abundances of the predicted table. **Mechanism:** the table is compositional, depth-confounded, and zero-inflated, on top of prediction error. **Symptom:** a long list of "significant" pathways that do not replicate across methods. **Fix:** use >=2 CoDA tools (ALDEx2, ANCOM-BC2, MaAsLin2, LinDA) and report the intersection (Nearing 2022 *Nat Commun* 13:342); ALDEx2 wants count-like features-as-rows, NOT relab-normalized output - see differential-abundance.

### Strain-level function invisible
**Trigger:** inferring strain-specific function (toxin, resistance, pathogenicity island) from a 16S-based prediction. **Mechanism:** 16S resolves to roughly genus/species; accessory genome, HGT, plasmids, and prophage-borne genes vary within a species and are assigned the reference neighbors' core content. **Symptom:** a strain-level functional claim from amplicon data. **Fix:** state that the species core is the ceiling regardless of NSTI; strain function needs isolate genomes or shotgun.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `--max_nsti` 2.0 (default) | Douglas 2020 *Nat Biotechnol* 38:685 | ASVs whose nearest sequenced genome is >2 substitutions/site away are too extrapolated to trust; dropped before inference - always report the dropped read fraction |
| `--hsp_method mp` (default) | PICRUSt2 manual | maximum parsimony is the recommended HSP method; `pic` is faster but not recommended (the legacy skill wrongly defaulted to `pic`) |
| `--placement_tool epa-ng` (default) | Barbera 2019 *Syst Biol* 68:365 | EPA-ng + gappa is the default placement path; `sepp` is the alternative |
| `--min_align` 0.8 (default) | PICRUSt2 manual | an ASV must align over >=80% of its length to be placed; poorly aligning ASVs are excluded |
| Predicted-vs-shotgun Spearman ~0.79-0.88 (gut) | Douglas 2020 *Nat Biotechnol* 38:685 | the empirical ceiling in the BEST case (dense human-gut references); lower elsewhere - the anchor for every accuracy caveat |
| NSTI ~0.5 / ~0.5-1 / ->2 (heuristic bands) | community practice | rough well-characterized / moderate / weak guide; no NSTI value converts predicted potential into measured function |
| DA: >=2 CoDA tools, report intersection | Nearing 2022 *Nat Commun* 13:342 | tool choice changes the predicted-pathway hit list; consensus beats any single tool |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `metadata_NSTI` / NSTI file not found | reading `marker_nsti_predicted.tsv` (does not exist) | the real file is `marker_predicted_and_nsti.tsv.gz`; the column is `metadata_NSTI` |
| Near-empty output / most ASVs dropped | high NSTI (wrong environment) | check the NSTI distribution; consider FAPROTAX or shotgun; do not just lower `--max_nsti` to keep them |
| ALDEx2 gives implausible results on predicted table | fed relab-normalized output, features-as-columns | pass raw (count-like) predicted abundances with features as rows - see differential-abundance |
| Predicted and shotgun pathway tables disagree | comparing PICRUSt2 to HUMAnN as if interchangeable | they are different objects (predicted vs measured); do not merge |
| `--per_sequence_contrib` produces nothing | used without `--stratified` | it is only meaningful with `--stratified` |
| Pathway "presence" and "abundance" conflated | reading `path_abun` as coverage | abundance and `--coverage` are different questions |

## References

- Douglas GM, Maffei VJ, Zaneveld JR, Yurgel SN, Brown JR, Taylor CM, Huttenhower C, Langille MGI. 2020. PICRUSt2 for prediction of metagenome functions. *Nat Biotechnol* 38:685-688.
- Langille MGI, Zaneveld J, Caporaso JG, et al. 2013. Predictive functional profiling of microbial communities using 16S rRNA marker gene sequences. *Nat Biotechnol* 31:814-821.
- Barbera P, Kozlov AM, Czech L, Morel B, Darriba D, Flouris T, Stamatakis A. 2019. EPA-ng: massively parallel evolutionary placement of genetic sequences. *Syst Biol* 68:365-369.
- Czech L, Barbera P, Stamatakis A. 2020. Genesis and Gappa: processing, analyzing and visualizing phylogenetic (placement) data. *Bioinformatics* 36:3263-3265.
- Louca S, Doebeli M. 2018. Efficient comparative phylogenetics on large trees. *Bioinformatics* 34:1053-1055.
- Ye Y, Doak TG. 2009. A parsimony approach to biological pathway reconstruction/inference for genomes and metagenomes. *PLoS Comput Biol* 5:e1000465.
- Wemheuer F, Taylor JA, Daniel R, Johnston E, Meinicke P, Thomas T, Wemheuer B. 2020. Tax4Fun2: prediction of habitat-specific functional profiles and functional redundancy based on 16S rRNA gene sequences. *Environ Microbiome* 15:11.
- Louca S, Parfrey LW, Doebeli M. 2016. Decoupling function and taxonomy in the global ocean microbiome. *Science* 353:1272-1277.
- Ward T, Larson J, Meulemans J, et al. 2017. BugBase predicts organism-level microbiome phenotypes. *bioRxiv* 133462 (preprint; not peer-reviewed).
- Jun SR, Robeson MS, Hauser LJ, Schadt CW, Gorin AA. 2015. PanFP: pangenome-based functional profiles for microbial communities. *BMC Res Notes* 8:479.
- Iwai S, Weinmaier T, Schmidt BL, et al. 2016. Piphillin: improved prediction of metagenomic content by direct inference from human microbiomes. *PLoS ONE* 11:e0166104.
- Nearing JT, Douglas GM, Hayes MG, et al. 2022. Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun* 13:342.

## Related Skills

- amplicon-processing - Generate the ASV table and representative sequences consumed here
- taxonomy-assignment - Taxonomic labels for the same ASVs (predicted function tracks these)
- differential-abundance - Compositional DA of the predicted KO/pathway table
- qiime2-workflow - The q2-picrust2 plugin path inside QIIME2
- metagenomics/functional-profiling - MEASURED shotgun function (HUMAnN); the predicted-vs-measured wall
- pathway-analysis/go-enrichment - Reading/enriching the predicted KO/MetaCyc lists
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
