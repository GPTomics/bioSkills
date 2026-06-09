---
name: bio-metagenomics-functional-profiling
description: Profiles the functional potential of shotgun metagenomes with HUMAnN 3's tiered search (MetaPhlAn prescreen, Bowtie2 pangenome, translated DIAMOND vs UniRef), giving gene-family (RPK) and MetaCyc pathway abundances stratified by species. Covers why a metagenome measures potential not activity, why dropping UNMAPPED/UNINTEGRATED biases everything, why stratification is an estimate, coverage-vs-abundance and MinPath/gap-fill, UniRef90-vs-50 and biome database bias, and the assembly/eggNOG/dbCAN/antiSMASH alternatives. Use when obtaining pathway or gene-family abundances, regrouping to KO/EC/GO, normalizing functional tables, or choosing read-based vs assembly-based functional profiling. For AMR genes see amr-detection; for host-gene enrichment see pathway-analysis.
tool_type: cli
primary_tool: HUMAnN
---

## Version Compatibility

Reference examples tested with: HUMAnN 3.6+, MetaPhlAn 4.1+, DIAMOND 2.1+, pandas 2.2+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `humann --version` then `humann --help` to confirm flags and defaults
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Output is driven by the reference databases: the ChocoPhlAn nucleotide pangenome, the UniRef90/50 protein database, and the MetaPhlAn database version HUMAnN calls. Record all three. Match the MetaPhlAn database version to the HUMAnN version when supplying `--taxonomic-profile` (MetaPhlAn 3 and 4 databases differ), and pin UniRef90 vs UniRef50, which changes both sensitivity and the UNMAPPED fraction.

# Functional Profiling

**"What can my community do?"** -> Quantify gene families and pathways with a tiered search that only translates the reads the fast steps could not place - measuring functional POTENTIAL the community encodes, never what it is expressing.
- CLI: `humann --input reads.fastq.gz --output out/ --taxonomic-profile sample_metaphlan.tsv --threads 8`

Scope: read-based community function (HUMAnN) and the assembly/specialized-database alternatives. Read classification -> kraken-classification, metaphlan-profiling. AMR gene quantification -> amr-detection. Host-gene over-representation (GO/KEGG/GSEA) -> the pathway-analysis category. Assembly/ORF mechanics -> genome-assembly/metagenome-assembly. Host depletion and trimming -> contamination-controls, read-qc.

## The Single Most Important Modern Insight -- A Metagenome Measures Potential, Not Activity

A gene family or a "complete pathway" in HUMAnN output is a CAPABILITY the community encodes - never a rate, a flux, or proof of expression. The DNA says the cell could ferment pyruvate; only RNA (metatranscriptome), protein, or metabolite data says it is. RNA functional profiles decouple from gene carriage (Franzosa 2014 *PNAS* 111:E2329), so narrating a pathabundance table as "the disease microbiome upregulates X" is the cardinal sin - it carries the gene more abundantly, nothing more. If the question is about activity, pair with metatranscriptomics: run the same HUMAnN on RNA and divide RNA-CPM by matched DNA-CPM per feature to get expression per gene copy. Memorable form: HUMAnN reports what the community CAN do, never what it IS doing. And every cell of the table is a model-dependent artifact of a reference database plus a tiered search at chosen thresholds - a hypothesis conditioned on the reference, not a measurement.

## The Tiered Search Is the Algorithm

HUMAnN does not brute-force-translate every read. Three tiers each filter the input to the next, so the slow translated step only sees what the fast steps could not place:

1. **Taxonomic prescreen (MetaPhlAn).** Detect species, then build a sample-specific ChocoPhlAn pangenome from only species above `--prescreen-threshold` (default 0.01%). Reuse via `--taxonomic-profile` to skip re-running MetaPhlAn.
2. **Nucleotide pangenome search (Bowtie2).** Reads mapping to that pangenome get a UniRef90 family WITH a high-confidence species label - this is where confident stratification comes from.
3. **Translated protein search (DIAMOND).** Reads that failed tier 2 are 6-frame translated and aligned to full UniRef90/50; they get a family but the species is INFERRED or `|unclassified` - lower-confidence stratification. Reads matching nothing become UNMAPPED.

This model is why `--bypass-*` flags and `--prescreen-threshold` change results, and why a stratified contribution from the translated tier is an estimate, not a measurement.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| HUMAnN 3 | Beghini 2021 *eLife* 10:e65088 | tiered read-based gene-family + MetaCyc pathway abundance | quantitative community function across samples |
| eggNOG-mapper v2 | Cantalapiedra 2021 *Mol Biol Evol* 38:5825 | orthology annotation (KO/GO/EC/CAZy/COG) of predicted ORFs | assembly route; flat functional catalogue |
| DIAMOND | Buchfink 2021 *Nat Methods* 18:366 | fast sensitive protein search (blastx) | custom read-vs-protein-DB profiling; backend of many tools |
| dbCAN3 | Zheng 2023 *Nucleic Acids Res* 51:W115 | CAZyme family/subfamily + substrate | carbohydrate-active enzymes (UniRef under-resolves these) |
| antiSMASH 7 | Blin 2023 *Nucleic Acids Res* 51:W46 | biosynthetic gene cluster detection | secondary-metabolite BGCs; CONTIGS only |
| MinPath | Ye & Doak 2009 *PLoS Comput Biol* 5:e1000465 | parsimony pathway calling inside HUMAnN | suppresses naive any-gene-implies-pathway over-calling |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Quantitative community function across samples | HUMAnN 3 (read-based) | counts every read; comprehensive; DB-bounded |
| Gut/host-associated, fine resolution | HUMAnN + UniRef90 | well-covered biome; specific families |
| Soil/marine/novel, big UNMAPPED | HUMAnN + UniRef50, or assembly route | UniRef90 cannot map divergent homologs |
| Need gene-to-organism/operon context or novel function | assembly + Prodigal + eggNOG-mapper | genomic context; accept loss of the unassembled majority |
| Carbohydrate-active enzymes | dbCAN3 | CAZy families/substrate beat generic UniRef |
| Biosynthetic gene clusters | antiSMASH (-> assembly first) | clusters span kb; require contigs |
| Activity, not capability | metatranscriptome (HUMAnN RNA mode) / RNA-DNA ratio | DNA cannot report expression |
| AMR gene quantification | -> amr-detection | dedicated ARG databases are the standard |
| Host-gene pathway enrichment | -> pathway-analysis | community pathway abundance is not GSEA |

## Run HUMAnN and Build a Functional Table

```bash
# Pre-QC first: adapter/quality trim AND host-deplete (e.g. KneadData). Host reads inflate UNMAPPED
# and waste DIAMOND time. Paired-end has no native pairing - concatenate R1+R2 into one file.
cat sample_R1.fq.gz sample_R2.fq.gz > sample.fq.gz
humann --input sample.fq.gz --output out/ \
    --taxonomic-profile sample_metaphlan.tsv \   # reuse the MetaPhlAn profile; do NOT --remove-temp-output and lose it
    --threads 8                                   # defaults: prescreen 0.01, translated-id 80 (uniref90), gap-fill on, minpath on

# Normalize PER SAMPLE before cross-sample stats (RPK is depth-dependent), then join and split.
humann_renorm_table -i out/sample_genefamilies.tsv -o out/sample_cpm.tsv -u cpm   # cpm preferred for models
humann_join_tables -i out -o merged_pathabundance.tsv --file_name pathabundance
humann_regroup_table -i merged_pathabundance.tsv -g uniref90_ko -o merged_ko.tsv  # adds an UNGROUPED row - keep it
humann_split_stratified_table -i merged_pathabundance.tsv -o .                     # run stats on the UNSTRATIFIED file
```

## Differential Abundance Without Biasing the Denominator

**Goal:** Test pathways between conditions without inventing abundance by discarding the unmapped fraction.

**Approach:** Keep UNMAPPED/UNINTEGRATED through normalization, check they do not differ by group (they often track the phenotype), then test the unstratified community totals with a compositional method (MaAsLin2/ANCOM-BC), not a bare Mann-Whitney on proportions.

```python
import pandas as pd

df = pd.read_csv('merged_pathabundance_unstratified.tsv', sep='\t', index_col=0)
meta = pd.read_csv('metadata.tsv', sep='\t', index_col=0)
unmapped = df.loc[['UNMAPPED', 'UNINTEGRATED']]            # the denominator - inspect, do not drop
g1 = meta.index[meta['condition'] == 'healthy']
g2 = meta.index[meta['condition'] == 'disease']
# If UNMAPPED differs by group, an assigned-feature comparison is confounded - report it.
unmapped_shift = unmapped[g2].mean(axis=1) - unmapped[g1].mean(axis=1)
print('UNMAPPED/UNINTEGRATED shift (disease - healthy):')
print(unmapped_shift.round(1))
# Then hand the unstratified table (UNMAPPED retained) to MaAsLin2/ANCOM-BC, which model
# compositionality and zero-inflation - not a raw t-test/Mann-Whitney on relative abundance.
```

## Per-Method Failure Modes

### Dropping UNMAPPED/UNINTEGRATED then renormalizing
**Trigger:** `.drop(['UNMAPPED','UNINTEGRATED'])` before relab. **Mechanism:** rescales assigned features to sum to 1, inventing abundance and erasing the database-coverage signal. **Symptom:** two samples with 20% vs 60% UNMAPPED look identical; a hit appears or vanishes. **Fix:** keep them through normalization; report them; if comparing assigned features only, confirm UNMAPPED does not differ by group.

### Stratification read as ground truth
**Trigger:** "species X contributes Y% of pathway Z." **Mechanism:** tier-2 species labels are confident, tier-3 (translated) are inferred or `|unclassified`. **Symptom:** over-confident organism-of-origin claims; ignored unclassified mass. **Fix:** treat contributions as estimates; flag the unclassified fraction; for confident gene-to-organism linkage use the assembly+binning route.

### Gut-centric QC thresholds applied to other biomes
**Trigger:** flagging a soil run as failed because UNMAPPED > 50%. **Mechanism:** UniRef is biased toward well-studied microbes; environmental biomes have large genuine dark function. **Symptom:** healthy environmental runs labeled failures. **Fix:** in novel biomes a big UNMAPPED is the environment, not a failure - drop to UniRef50 for sensitivity or switch to the assembly route.

### Coverage vs abundance confusion; gap-fill/MinPath artifacts
**Trigger:** reading `pathcoverage` as abundance or trusting a "complete" pathway. **Mechanism:** gap-fill (default on) scores pathways with missing reactions; MinPath parsimony prunes redundant pathways so absence is not biological absence; coverage is de-emphasized in recent docs. **Symptom:** fabricated completeness or missing real-but-redundant pathways. **Fix:** use abundance for stats; treat coverage as a soft presence prior; know gap-fill and MinPath are on by default.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `--prescreen-threshold` 0.01% | HUMAnN docs | species below this are excluded from the tier-2 pangenome |
| `--translated-identity-threshold` 80 (UniRef90) / 50 (UniRef50) | HUMAnN docs | the sensitivity knob; divergent homologs fail the 80% cutoff |
| nucleotide/translated coverage 90/50 | HUMAnN docs | query 90%, subject 50% coverage filters |
| `--evalue` 1.0 | HUMAnN docs | permissive DIAMOND e-value; tier design controls specificity |
| gap-fill on, MinPath on | HUMAnN docs | sensitivity/parsimony defaults that shape pathway presence |
| Normalize RPK -> CPM per sample before stats | HUMAnN docs | RPK is depth-dependent; CPM preferred for linear/log models |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Huge UNMAPPED on a gut sample | host reads not removed | host-deplete + trim before HUMAnN |
| MetaPhlAn profile gone after run | `--remove-temp-output` deleted `_humann_temp/` | omit it; keep the profile for reuse and taxonomy |
| `--taxonomic-profile` rejected | MetaPhlAn DB version mismatch with HUMAnN | match the MetaPhlAn database to the HUMAnN version |
| Pathways look over-called | MinPath off / naive any-gene mapping | keep MinPath on (default) |
| Stratified DA is all zeros/noise | bare t-test on zero-inflated strata | run stats on the unstratified table with MaAsLin2/ANCOM-BC |
| CAZymes under-annotated by UniRef | generic database under-resolves CAZy | use dbCAN3 on predicted ORFs |

## References

- Beghini F, McIver LJ, Blanco-Miguez A, et al. 2021. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife* 10:e65088.
- Franzosa EA, Morgan XC, Segata N, et al. 2014. Relating the metatranscriptome and metagenome of the human gut. *PNAS* 111:E2329-E2338.
- Ye Y, Doak TG. 2009. A parsimony approach to biological pathway reconstruction/inference for genomes and metagenomes. *PLoS Comput Biol* 5:e1000465.
- Buchfink B, Reuter K, Drost HG. 2021. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods* 18:366-368.
- Cantalapiedra CP, Hernandez-Plaza A, Letunic I, Bork P, Huerta-Cepas J. 2021. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol* 38:5825-5829.
- Hyatt D, Chen GL, LoCascio PF, et al. 2010. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11:119.
- Zheng J, Hu B, Zhang X, et al. 2023. dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Res* 51:W115-W121.
- Blin K, Shaw S, Augustijn HE, et al. 2023. antiSMASH 7.0: new and improved predictions. *Nucleic Acids Res* 51:W46-W50.

## Related Skills

- metaphlan-profiling - The taxonomic prescreen HUMAnN reuses via --taxonomic-profile
- kraken-classification - Alternative taxonomic input
- abundance-estimation - Compositional normalization shared with functional tables
- amr-detection - Dedicated ARG quantification (HUMAnN can surface AMR families but is not standard)
- metagenome-visualization - Plot and test functional tables
- contamination-controls - Host depletion before HUMAnN
- genome-assembly/metagenome-assembly - Assembly route for contextualized/novel function
- pathway-analysis/kegg-pathways - Organism-centric pathway interpretation (not community abundance)
