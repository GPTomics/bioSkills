---
name: bio-epidemiological-genomics-pathogen-typing
description: Assigns isolate identity at the right resolution for the question -- ANI/Mash species triage, 7-locus MLST historical comparability, cgMLST/wgMLST outbreak resolution (chewBBACA, BIGSdb, Ridom SeqSphere, EnteroBase HierCC), in-silico serotyping (SISTR/SeqSero2 Salmonella, SerotypeFinder E. coli, Kaptive Klebsiella, SeroBA pneumococcus, spa+SCCmec S. aureus), and lineage callers (TB-Profiler/Mykrobe barcode for MTBC, Pangolin + Nextclade for SARS-CoV-2, PopPUNK GPSC for S. pneumoniae). Use when typing bacterial isolates for surveillance or outbreak investigation, choosing between cgMLST allele distance and core-SNP distance for cluster definition, harmonising calls across schemas/database versions, assigning MTBC lineage with the Napier 90-SNP barcode, calling Salmonella serovar via SISTR with monophasic Typhimurium awareness, running Pangolin UShER mode with pangolin-data version pinning, or selecting a typing resolution to match the surveillance question.
tool_type: mixed
primary_tool: chewBBACA
---

## Version Compatibility

Reference examples tested with: mlst 2.23+, chewBBACA 3.3+, SISTR 1.1+, SeqSero2 1.3+, SerotypeFinder 2.0+, Kleborate 3.0+, Kaptive 3.0+, SeroBA 1.0+, PopPUNK 2.7+, pangolin 4.3+ (pangolin-data 1.30+), nextclade 3.8+, tb-profiler 6.2+, mykrobe 0.13+, mash 2.3+, skani 0.2+, snippy 4.6+, snp-dists 0.8+, pandas 2.2+, BioPython 1.84+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Pangolin: `pangolin --all-versions` records pangolin + pangolin-data + scorpio + constellations
- Nextclade dataset: `nextclade dataset list --tag latest sars-cov-2`
- chewBBACA schema: schemas are versioned independently; record the schema source AND the schema-fetch date for any published call
- TB-Profiler bundled DB: `tb-profiler list_db` for the WHO catalogue edition currently bundled

If a tool reports an unexpected serovar / lineage / ST, the schema or barcode version is the first thing to check; introspect the installed package and database date before re-running.

# Pathogen Typing

**"What is this isolate, and is it the same as that one?"** -> Pick a typing resolution that matches the epidemiological question, then run the appropriate caller and report the call with explicit schema / database / lineage-barcode versions. Resolution mismatch is the single most common typing error -- using 7-locus MLST to investigate a 20-isolate outbreak (multiple unrelated isolates share an ST) and using cgMLST to track decade-long lineage trends (allele drift obscures the lineage signal) are equally wrong, in opposite directions.

- CLI: `mlst assembly.fa` -- 7-locus MLST via Seemann's mlst with PubMLST schemas
- CLI: `chewBBACA.py AlleleCall -i assemblies/ -g schema/ -o alleles/ --cpu 8` -- cgMLST allele profile
- CLI: `pangolin sequences.fasta --analysis-mode usher --all-versions` -- SARS-CoV-2 Pango lineage with version provenance
- CLI: `tb-profiler profile -1 r1.fq.gz -2 r2.fq.gz -p sample` -- MTBC lineage (Coll/Napier barcode) plus DR call
- Python: `pandas` + `snp-dists` for cluster definition with pathogen-specific thresholds

## The Single Most Important Modern Insight -- cgMLST and core-SNP distance answer different questions

cgMLST counts allele changes (one wobble within a 1-kb locus = 1 allelic difference, regardless of how many SNPs sit inside that allele) and is robust to small assembly artifacts; core-SNP counts every position and is sensitive to alignment / mapping artifacts but recovers higher resolution. The EFSA harmonised foodborne approach uses cgMLST for *Salmonella* / *Listeria* / *E. coli* (Salm cgMLST <=5 alleles = cluster; PulseNet *Listeria* <=4 alleles); UK / EU TB outbreak literature uses core-SNP with recombination masking (Walker 2013: <=12 SNPs = likely transmission; <=5 = recent). Mixing the two yields incompatible cluster definitions: the same 50-isolate outbreak under cgMLST may cluster as 1 group at threshold 5, under core-SNP as 3 groups at threshold 12. The output of pathogen-typing is not a single distance -- it is a distance metric, a threshold derived from a specific population, and a schema or reference version. All three must travel with the call.

## Algorithmic Taxonomy

| Method | Mechanism | Resolution | Strength | Fails when |
|--------|-----------|-----------|----------|------------|
| ANI / fastANI / skani | Pairwise nucleotide identity over orthologous regions | Species (>=95% ANI = same species) | Species-level QC; cross-genus is non-metric for Mash | Species below ~80% ANI; Mash distance violates triangle inequality |
| Mash (Ondov 2016 *Genome Biol* 17:132) | MinHash sketches; fast pairwise distances | Species triage | Seconds per pair; the rapid screening standard | Distance non-metric at low identity; ANI <80% unreliable |
| 7-locus MLST (Maiden 2013 *Nat Rev Microbiol* 11:728; Seemann mlst tool) | PubMLST allele lookup at 7 housekeeping loci | Sequence type | Historical comparability across decades | Insufficient resolution for outbreaks; multiple unrelated isolates may share ST |
| cgMLST (chewBBACA; BIGSdb; Ridom SeqSphere) | Allele lookup across ~1000-3000 core loci | Outbreak / multi-country surveillance | Allele-distance is robust to small mapping errors | Schema-version-dependent; cross-schema NON-comparable; missing-locus handling matters |
| wgMLST | Allele lookup across all genes including accessory | Highest typing resolution | Maximum resolution within schema | Even more schema-dependent than cgMLST |
| Core-SNP typing (snippy + snp-dists; Parsnp; Lyve-SET) | Reference-based SNP calling on core genome | Single-SNP resolution | Highest resolution; underlies most Mtb work | Reference-dependent; recombination must be masked for bacteria |
| HierCC (Zhou 2021 *Bioinformatics* 37:3645) | Hierarchical clustering on cgMLST at multiple thresholds (HC5, HC10, HC50, etc.) | Stable nomenclature across schema updates | Cross-schema-update stability for EnteroBase pathogens | Limited to EnteroBase organisms |
| PopPUNK GPSC (Lees 2019 *Genome Res* 29:304) | k-mer-based clustering with Gaussian mixture on core+accessory distance | Population cluster | Stable IDs across additions; scales to >100k genomes | Cluster membership of an individual isolate can shift as model is refined |
| Pangolin (O'Toole 2021 *Virus Evol* 7:veab064) | Phylogenetic placement (UShER) or ML classifier (pangoLEARN) | SARS-CoV-2 Pango lineage | Curated dynamic nomenclature; recombinant X-prefix designations | pangoLEARN deprecated mid-2023; lab-to-lab version skew silently flips lineage calls |
| Nextclade (Aksamentov 2021 *JOSS* 6:3773) | Reference-tree placement + clade calling + QC | Clade nomenclature + mutations + QC | Mutation reports + QC integrated; multi-pathogen datasets | Dataset version drift changes which mutations count as "lineage-defining" |
| TB-Profiler + Coll/Napier barcode (Phelan 2019; Coll 2014; Napier 2020) | Reference-based SNP call against H37Rv + barcode SNP set | MTBC lineage 1-9 + drug resistance | Integrated lineage + DST; supports Napier 90-SNP barcode covering lineages 7-9 | Pre-Napier 2020 barcodes miscall lineage 7-9 isolates |
| Mykrobe (Hunt 2019 *Wellcome Open Res* 4:191) | k-mer presence/absence panels | Species + AMR for TB / S. aureus / Salmonella | Fast; cross-check for TB-Profiler | Panel may lag WHO catalogue |
| SISTR (Yoshida 2016 *PLoS ONE* 11:e0147101) | cgMLST + ribosomal MLST + serovar inference | Salmonella serovar + antigenic formula | ~94% concordance with traditional sero-typing; monophasic-aware | Novel/rare serovars without reference panel; antigen-cluster regulatory mutations are silent |
| SeqSero2 (Zhang 2019 *AEM* 85:e01746-19) | k-mer + targeted-assembly for Salmonella | Salmonella serovar | Designed for low-coverage / fragmented data | Slightly different output schema than SISTR |
| SerotypeFinder (Joensen 2015 *J Clin Microbiol* 53:2410) | BLAST against O- and H-antigen biosynthesis genes | E. coli O:H | Standard for E. coli serotyping | Misses novel O / H types; fimH typing is separate |
| Kleborate + Kaptive (Lam 2021 *Nat Commun* 12:4188) | Integrated MLST + K/O typing + virulence (ICEKp / iuc / ybt / clb / iro) + AMR | Klebsiella surveillance | Hypervirulence vs classical distinction; K/O loci typed by Kaptive | Kaptive K-locus DB versioned; KL calls can flip between Kaptive v1 / v2 / v3 |
| spa + SCCmec (Harmsen 2003 *J Clin Microbiol* 41:5442; Kaya 2018 *mSphere* 3:e00612-17) | spa repeat-region typing + SCCmec cassette typing | S. aureus typing | Historical comparability; clinical surveillance | spa repeat array fragments at borderline read length; assembler-dependent |
| SeroBA (Epping 2018 *Microb Genom* 4:e000186) | k-mer-based serotyping from raw reads | S. pneumoniae serotype | 98% concordance; no assembly needed; runs on >=15x coverage | Vaccine replacement makes serotype the load-bearing surveillance unit |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| "Is this strain even what we think it is?" | Mash / skani ANI triage; ANI >=95% confirms species | 7-locus MLST cannot answer species; cgMLST schema may not load on wrong species |
| Routine *Salmonella* serotyping | SISTR (assembled) or SeqSero2 (low coverage / fragmented); flag monophasic 1,4,[5],12:i:- explicitly | Single tool without monophasic awareness; reporting Typhimurium when the antigen cluster is deleted |
| *Klebsiella* surveillance | Kleborate (integrates MLST + K/O via Kaptive + ICEKp virulence + AMR via AMRFinderPlus) | MLST + Kaptive + AMR run separately and not integrated -- missing hypervirulence call |
| *S. pneumoniae* surveillance | SeroBA serotype + PopPUNK GPSC; vaccine-replacement is the central post-PCV story | Reporting MLST ST without serotype (serotype IS the vaccine-actionable surveillance unit) |
| *S. aureus* typing | spa (Ridom) + SCCmec (SCCmecFinder) + MLST + clonal complex; flag CC8 USA300 / CC22 EMRSA-15 / CC30 EMRSA-16 / CC398 livestock | Just spa without CC; spa repeat assembly artifacts unflagged |
| *M. tuberculosis* lineage | TB-Profiler primary (Coll/Napier barcode) + Mykrobe cross-check; verify Napier 2020 90-SNP barcode (covers lineages 7-9) | Pre-Napier barcodes miscall L7-9; MIRU-VNTR is obsolete for new surveillance |
| SARS-CoV-2 lineage | Pangolin with `--analysis-mode usher` (UShER default since v4; pangoLEARN deprecated mid-2023) + Nextclade cross-check; document pangolin-data version | pangoLEARN alone; reporting lineage without pangolin-data version pin |
| Define outbreak cluster | cgMLST allele-distance with pathogen-tuned threshold (Salm <=5; *Listeria* PulseNet <=4; TB <=12 SNPs core; *C. difficile* <=2 SNPs core); cite the threshold's source population | Universal SNP threshold across pathogens (10x variation across taxa); applying Walker 2013 UK thresholds in high-transmission settings |
| Reproducible nomenclature across years | PopPUNK GPSC (S. pneumoniae) / Pangolin (SARS-CoV-2) / HierCC (EnteroBase pathogens) / SISTR serovar -- curated stable IDs | Raw cgMLST allele profile as the surveillance unit -- not stable across schema updates |
| Multi-lab cross-comparison | Same schema source AND same schema version AND same tool version; Pathogenwatch / EnteroBase shared platform | Locally computed cgMLST profiles compared across labs without schema versioning |

Methodology evolves; before any high-stakes typing report, verify Pangolin's current default analysis-mode and the Napier barcode currently bundled in TB-Profiler.

## Running 7-Locus MLST and cgMLST

**Goal:** Produce reproducible ST and cgMLST allele-distance calls for a cohort of bacterial isolates, with schema versioning preserved for cross-lab comparison.

**Approach:** Run Seemann's `mlst` for the 7-locus baseline (uses bundled PubMLST schemas); for cgMLST, fetch the schema explicitly with `chewBBACA.py DownloadSchema` (records the source and date), then `AlleleCall`, then `ExtractCgMLST` with the per-organism missing-locus threshold; compute allele distances on the pairwise-complete intersection of called loci.

```bash
mlst --threads 8 assemblies/*.fa > cohort.mlst.tsv

chewBBACA.py DownloadSchema \
    -sp "Salmonella enterica" \
    -sc cgMLST \
    -o schema_dir
SCHEMA_DATE=$(date -u +%Y-%m-%d)

chewBBACA.py AlleleCall \
    -i assemblies/ \
    -g schema_dir/cgMLST \
    -o alleles_out/ \
    --cpu 8

chewBBACA.py ExtractCgMLST \
    -i alleles_out/results_alleles.tsv \
    -o cgmlst_profile.tsv \
    --threshold 0.95

echo "schema_source: chewBBACA SalmonellaEnterica cgMLST" > cgmlst_profile.metadata
echo "schema_fetched: ${SCHEMA_DATE}" >> cgmlst_profile.metadata
```

`AlleleCall` output uses special codes: `LNF` (locus not found), `PLOT` (truncated), `NIPH` (non-informative paralog), `ASM` (allele small/short), `ALM` (allele large), and asterisk (new allele). These are MISSING DATA. Counting "0" or "-" against another "0" or "-" as 1 allelic difference is a quiet correctness bug that affects most cgMLST cluster definitions outside specialist labs.

## SNP-Based Outbreak Cluster Definition

**Goal:** Identify isolates within an outbreak threshold using core-SNP distances on a recombination-aware alignment, with the pathogen-specific threshold cited from its source population.

**Approach:** Snippy against a high-quality reference; `snippy-core` to extract core SNPs; for bacteria, Gubbins to mask recombinant tracts before counting (otherwise apparent SNP distance is inflated by recombination); `snp-dists` for pairwise distance matrix; apply the published pathogen-specific threshold (citing Walker 2013 for TB, Eyre 2013 for *C. difficile*, EFSA convention for *Salmonella*).

```bash
for r1 in reads/*_R1.fq.gz; do
    sample=$(basename "${r1}" _R1.fq.gz)
    r2="reads/${sample}_R2.fq.gz"
    snippy --outdir snippy_out/${sample} --R1 "${r1}" --R2 "${r2}" --reference reference.fa --cpus 8
done

snippy-core --ref reference.fa --prefix core snippy_out/*

run_gubbins.py --prefix gubbins core.full.aln

snp-dists -c gubbins.filtered_polymorphic_sites.fasta > cohort.snp_dists.csv
```

`run_gubbins.py` input MUST be `core.full.aln` (full-position alignment with reference). Passing `core.aln` (variable positions only) produces wrong recombination calls because Gubbins cannot estimate background SNP density without invariant positions.

## Lineage Calling for Mtb and SARS-CoV-2

**Goal:** Assign MTBC lineage via the Napier 2020 barcode or SARS-CoV-2 Pango lineage via UShER placement, with explicit version pinning so the call is reproducible.

**Approach:** TB-Profiler bundles the Coll 2014 + Napier 2020 barcode and reports lineage 1-9 (Napier 2020 added lineages 7-9 that older 62-SNP barcodes miss); Pangolin with `--analysis-mode usher` performs phylogenetic placement on the daily-updated UShER tree and is preferred over pangoLEARN since v4 (Pongmoragot 2024 *Virus Evol* 10:vead085); always record `pangolin --all-versions` output alongside the lineage call.

```bash
tb-profiler profile -1 reads_R1.fq.gz -2 reads_R2.fq.gz -p sample --dir tbp_out
mykrobe predict --sample sample --species tb --output sample.mykrobe.json --format json reads_R1.fq.gz reads_R2.fq.gz

pangolin sequences.fasta --analysis-mode usher --outfile lineage_report.csv
pangolin --all-versions > pangolin_versions.txt

nextclade dataset get --name sars-cov-2 --output-dir nc_dataset/sars-cov-2
NC_DATASET_TAG=$(jq -r '.tag' nc_dataset/sars-cov-2/pathogen.json)

nextclade run \
    --input-dataset nc_dataset/sars-cov-2 \
    --output-tsv nextclade.tsv \
    --output-json nextclade.json \
    sequences.fasta

echo "nextclade_dataset_tag: ${NC_DATASET_TAG}" > nextclade.metadata
```

## Per-Method Failure Modes

### chewBBACA missing-locus codes counted as allelic differences

**Trigger:** A naive cgMLST distance computation that treats LNF / PLOT / NIPH / ASM / ALM / "0" / "-" as integer alleles and counts them against other missing codes.

**Mechanism:** chewBBACA output uses special codes for missing or paralogous loci. The pairwise allele distance must be computed on the intersection of called loci, not the union. Treating two LNFs as "same allele" (=0 difference) or as "different alleles" (=1 difference) are both wrong; the locus must be excluded from the comparison.

**Symptom:** Outbreak cluster definition flips depending on completeness of the assemblies; samples with more missing data appear artificially close to each other (if missing-vs-missing counts 0) or far from everyone (if missing-vs-allele counts 1).

**Fix:** Use chewBBACA's `ExtractCgMLST` with a `--threshold` (typically 0.95) to drop loci called in fewer than that fraction of samples; for the pairwise distance, restrict to loci called in BOTH samples (pairwise-complete) and report the number of loci compared alongside the distance.

### Pangolin lineage call flips between lab A (older pangolin-data) and lab B (current)

**Trigger:** Two labs submit the same consensus genome to Pangolin with different pangolin-data versions; the lineage call differs (e.g., BA.2 vs BA.2.86 vs JN.1).

**Mechanism:** Lineage designation happens through pango-designation GitHub issues -- community-driven, often days-to-weeks before pangolin-data releases include the lineage. During this window, the same genome is callable as the parent lineage (older pangolin-data) or the child (current pangolin-data). pangolin-data is updated weekly; lab-to-lab version skew is routine.

**Symptom:** Cross-lab lineage prevalence comparisons over time show implausible jumps that coincide with pangolin-data release dates rather than biology.

**Fix:** Pin pangolin-data version explicitly with `pangolin --all-versions` recorded alongside every call. For published or regulatory output, re-run the WHOLE archive against a single pangolin-data version before reporting. Comparing today's BA.2.86 call to last month's "Unassigned" call is invalid.

### MTBC lineage 7-9 miscalled by pre-Napier barcode

**Trigger:** TB-Profiler / Mykrobe running an older bundled barcode (Coll 2014 62-SNP, pre-2020 builds); isolate is from Ethiopia, Rwanda, or East Africa.

**Mechanism:** The Coll 2014 62-SNP barcode covers lineages 1-7 but predates the formal designation of lineages 8 (Rwanda) and 9. The Napier 2020 90-SNP barcode (Napier *Genome Med* 12:114) adds these and refines L4 sublineages. An older bundled barcode silently maps lineage 8 / 9 isolates to "unknown" or to a nearest-barcode-match L4 sub-lineage based on partial SNPs.

**Symptom:** Cross-lab Mtb surveillance dataset shows Ethiopia / Rwanda isolates as a mix of "unknown" and unexpected L4 sublineages.

**Fix:** Verify the bundled barcode version (`tb-profiler list_db`); update to the Napier 2020 barcode or later before any lineage-stratified analysis. Re-run historical Mtb data against the current barcode whenever the barcode is updated.

### Beijing-lineage Mtb literature based on spoligotype, not WGS sublineage

**Trigger:** Citing pre-2014 "Beijing family" findings (hypervirulence claims, vaccine-escape claims, faster-evolution claims) as if they apply to a single modern WGS sublineage.

**Mechanism:** Pre-WGS, the *Mtb* "Beijing" family was defined by spoligotype pattern (absence of DR spacers 1-34, presence of 35-43). WGS revealed Beijing is a paraphyletic grouping containing multiple sublineages with distinct phenotypes (modern Beijing = lineage 2.2.1.1; ancestral Beijing = 2.2.1.2; proto-Beijing = 2.1). Much of the pre-2014 Beijing literature was based on the spoligotype-defined paraphyletic group.

**Symptom:** Beijing-as-a-monolith conclusions inappropriately applied to modern WGS-typed isolates that may fall in different sublineages of L2.

**Fix:** Always specify the WGS sublineage (e.g., "2.2.1.1 modern Beijing", "2.2.1.2 ancestral Beijing") rather than "Beijing". Treat pre-2014 Beijing claims as hypotheses to be re-validated on WGS-typed cohorts.

### Mash distance used to cluster genomes below 80% ANI

**Trigger:** Mash-distance hierarchical clustering of distantly related genomes (multi-genus or low-identity comparisons).

**Mechanism:** Mash distance (Ondov 2016 *Genome Biol* 17:132) is a Jaccard-derived estimator of mutation rate, validated for ANI >=80% (the Mash docs state 95% confidence at ANI >=90%). Below that, the Mash distance becomes non-metric -- A-B + B-C can be < A-C -- and clustering algorithms that assume metric distances produce undefined output.

**Symptom:** Cross-genus or low-identity comparison clusters do not match phylogenetic expectations; the same cluster definition is unstable to addition of new genomes.

**Fix:** Use Mash only within the validated ANI range (>=80%, ideally >=90%). For cross-genus or low-identity comparisons use AAI or ANI-from-alignment (`skani`, `pyani`). Document the ANI range of the cohort before computing any Mash-based distance.

### Cross-schema cgMLST distances pooled as if comparable

**Trigger:** Combining cgMLST distances from a chewBBACA local schema with distances from a Ridom SeqSphere schema in the same outbreak comparison.

**Mechanism:** Ridom SeqSphere, chewBBACA-built schemas, and EnteroBase schemas are independently curated and use different locus sets and different allele numbering. A "cgMLST distance of 5" between two isolates does NOT equal a "cgMLST distance of 5" under a different schema.

**Symptom:** Multi-country outbreak comparison reports irreconcilable cluster definitions; the same isolate pair appears in-cluster in one lab and out-of-cluster in another.

**Fix:** Document schema source AND schema-fetch date for every cgMLST profile; do not pool distances across schema sources. For multi-country collaboration, agree on a single schema (typically EnteroBase HierCC for *Salmonella* / *E. coli* / *Listeria*) at the outset.

### K-locus typing flips between Kaptive versions

**Trigger:** Comparing *K. pneumoniae* K-locus prevalence between studies that used different Kaptive database versions.

**Mechanism:** The Kaptive K-locus DB has been versioned multiple times since release; K-locus calls have flipped (KL1 / KL2 boundary refinements; KL149+ additions) between Kaptive v1 and v2. Kaptive v3 (2024) added O-locus typing and an L-locus scheme.

**Symptom:** Longitudinal K-locus prevalence trend has implausible jumps coinciding with Kaptive release dates rather than biology.

**Fix:** Document Kleborate AND Kaptive version with every K/O call. For longitudinal trend analysis, re-run all historical assemblies against a single Kaptive version.

## Reconciliation: When Typing Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Pangolin "BA.2.86", Nextclade clade "23I" | Equivalent at different resolutions -- BA.2.86 is within 23I | Report both; Pango lineage for sub-clade resolution |
| Pangolin "BA.5.2", Nextclade "Unassigned" | Nextclade dataset older than pangolin-data; OR Nextclade QC failed | Update Nextclade dataset; re-run; inspect QC fields |
| Pangolin UShER and pangoLEARN disagree | pangoLEARN is the deprecated decision-tree classifier (deprecated mid-2023) | Trust UShER call |
| SISTR "monophasic Typhimurium 1,4,[5],12:i:-", slide agglutination "Typhimurium" | fljB gene deleted in monophasic variant; slide agglutination cannot detect | Trust genome call; flag for confirmatory testing if surveillance requires |
| cgMLST cluster definition flips between two schemas | Different locus sets; non-comparable | Pick one schema for the entire analysis; document |
| TB lineage call differs between TB-Profiler and Mykrobe | Different bundled barcode versions; one may predate Napier 2020 | Update both; if disagreement persists, manual barcode SNP inspection |
| Two consecutive pangolin-data releases call the same consensus differently | Lineage definitions revised between releases | Pin pangolin-data; re-run whole archive on dataset update |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| Mash / ANI species boundary | >=95% ANI | ANI species-delineation convention |
| Mash distance validity range | >=80% ANI (>=90% for 95% confidence) | Ondov 2016 *Genome Biol* 17:132 |
| cgMLST cluster -- *Salmonella* (chewBBACA Salm scheme) | <=5 allelic differences | EFSA harmonised approach |
| cgMLST cluster -- *Listeria monocytogenes* (PulseNet) | <=4 allelic differences | PulseNet protocol convention |
| cgMLST cluster -- *E. coli* (EnteroBase) | <=10 allelic differences (STEC outbreak) | EnteroBase convention |
| Core SNP -- *M. tuberculosis* | <=5 SNPs (recent transmission); <=12 SNPs (likely transmission) | Walker 2013 *Lancet Infect Dis* 13:137 (UK low-transmission setting) |
| Core SNP -- *Staphylococcus aureus* | <=15 SNPs (within hospital outbreak); <=40 SNPs (broader temporal cluster) | Coll 2017 *Clin Infect Dis* 65:1781 |
| Core SNP -- *Klebsiella pneumoniae* (KPC outbreak) | <=21 SNPs | Snitkin 2012 *Sci Transl Med* 4:148ra116 |
| Core SNP -- *Clostridioides difficile* (recombination-masked) | <=2 SNPs (likely direct transmission); <=10 (plausible within 6 months) | Eyre 2013 *NEJM* 369:1195 |
| Core SNP -- *Neisseria gonorrhoeae* | <=25 SNPs (transmission) | UKHSA STI framework |
| SARS-CoV-2 cluster definition | NOT defined by SNP alone; combine 0-2 SNPs + epi link + sampling window | SARS-CoV-2 within-host diversity literature |
| chewBBACA cgMLST extraction completeness | 0.95 (drop loci called in <95% of samples) | chewBBACA convention |
| SeroBA minimum coverage | >=15x | Epping 2018 *Microb Genom* 4:e000186 |

CRITICAL: a number from one pathogen does NOT transfer to another. Substitution rate, recombination, host range, generation interval, and within-host diversity vary by 100x across pathogens. Always cite the threshold's source population, especially for Walker 2013 (UK low-transmission setting) which routinely inflates apparent recent-transmission rates by 2-5x in high-burden settings.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Outbreak split into spurious clusters | cgMLST missing-locus codes counted as allelic differences | Restrict to pairwise-complete loci |
| Lineage prevalence shows implausible jump | pangolin-data version drift | Pin pangolin-data; re-run archive |
| Mtb lineage call "unknown" for East African isolate | Pre-Napier barcode | Update to Napier 2020 |
| Salmonella sample called "Typhimurium" when antigen cluster deleted | Tool did not flag monophasic variant | Use SISTR or SeqSero2 with explicit monophasic check |
| `nextclade run --input-dataset` rejected | v2 syntax in v3 | v3 uses `--input-dataset` for pre-downloaded folder; verify `nextclade --version` |
| `pangolin --inference usher` not recognised | Flag is `--analysis-mode usher` | Use `--analysis-mode usher` |
| chewBBACA fails on assembly with broken loci | Schema BSR cutoff too strict; assembly too fragmented | Document quality; consider re-assembly with long reads |
| GPSC cluster differs between studies | PopPUNK model updated; individual cluster membership can shift | Re-run PopPUNK against current reference DB before comparison |
| Mash hierarchical clustering unstable | Cohort spans wide ANI range; below validity threshold | Use ANI-from-alignment (skani / pyANI) for cross-genus comparisons |
| MOB-suite plasmid context conflicts with cgMLST cluster | cgMLST is chromosome-focused; plasmid distance is independent | Report both with framing |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why cgMLST and not core-SNP?" | EFSA harmonised foodborne uses cgMLST; cross-lab comparability via shared schema. For outbreak-internal who-infected-whom, core-SNP supplements |
| "What threshold was used, and on what population?" | Cite Walker 2013 / Eyre 2013 / EFSA per pathogen; caveat the population for non-universal thresholds |
| "How were missing cgMLST loci handled?" | Pairwise-complete distance; locus must be called in BOTH samples to count |
| "Pangolin version?" | `pangolin --all-versions` recorded; re-run on dataset update |
| "Why TB-Profiler over Mykrobe?" | TB-Profiler is primary (WHO catalogue integration + Napier barcode); Mykrobe is cross-check on R/XDR calls |
| "Was the Napier 2020 barcode checked?" | Verified `tb-profiler list_db`; lineage 7-9 callable |
| "Was the Beijing-as-monolith literature cited?" | Disaggregated to WGS sublineage (2.2.1.1 modern / 2.2.1.2 ancestral); pre-2014 Beijing claims treated as hypotheses |
| "Multi-country outbreak: harmonised schema?" | EnteroBase HierCC for *Salmonella* / *E. coli* / *Listeria*; documented schema source + date |

## References

- Maiden MCJ, van Rensburg MJJ, Bray JE et al (2013) MLST revisited: the gene-by-gene approach to bacterial genomics. *Nat Rev Microbiol* 11(10):728-736. doi:10.1038/nrmicro3093
- Silva M, Machado MP, Silva DN et al (2018) chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. *Microb Genom* 4(3):e000166. doi:10.1099/mgen.0.000166
- Zhou Z, Charlesworth J, Achtman M (2021) HierCC: a multi-level clustering scheme for population assignments based on core genome MLST. *Bioinformatics* 37(20):3645-3646. doi:10.1093/bioinformatics/btab234
- Yoshida CE, Kruczkiewicz P, Laing CR et al (2016) The Salmonella In Silico Typing Resource (SISTR). *PLoS ONE* 11(1):e0147101. doi:10.1371/journal.pone.0147101
- Zhang S, den Bakker HC, Li S et al (2019) SeqSero2: rapid and improved Salmonella serotype determination using whole-genome sequencing data. *Appl Environ Microbiol* 85(23):e01746-19. doi:10.1128/AEM.01746-19
- Joensen KG, Tetzschner AMM, Iguchi A et al (2015) Rapid and easy in silico serotyping of Escherichia coli isolates by use of whole-genome sequencing data. *J Clin Microbiol* 53(8):2410-2426. doi:10.1128/JCM.00008-15
- Lam MMC, Wick RR, Watts SC et al (2021) A genomic surveillance framework and genotyping tool for Klebsiella pneumoniae and its related species complex. *Nat Commun* 12:4188. doi:10.1038/s41467-021-24448-3
- Lees JA, Harris SR, Tonkin-Hill G et al (2019) Fast and flexible bacterial genomic epidemiology with PopPUNK. *Genome Res* 29(2):304-316. doi:10.1101/gr.241455.118
- Epping L, van Tonder AJ, Gladstone RA et al (2018) SeroBA: rapid high-throughput serotyping of Streptococcus pneumoniae from whole genome sequence data. *Microb Genom* 4(7):e000186. doi:10.1099/mgen.0.000186
- Harmsen D, Claus H, Witte W et al (2003) Typing of methicillin-resistant Staphylococcus aureus in a university hospital setting by using novel software for spa repeat determination and database management. *J Clin Microbiol* 41(12):5442-5448. doi:10.1128/JCM.41.12.5442-5448.2003
- Kaya H, Hasman H, Larsen J et al (2018) SCCmecFinder, a web-based tool for typing of staphylococcal cassette chromosome mec in Staphylococcus aureus using whole-genome sequence data. *mSphere* 3(1):e00612-17. doi:10.1128/mSphere.00612-17
- Coll F, McNerney R, Guerra-Assunção JA et al (2014) A robust SNP barcode for typing Mycobacterium tuberculosis complex strains. *Nat Commun* 5:4812. doi:10.1038/ncomms5812
- Coll F, Harrison EM, Toleman MS et al (2017) Longitudinal genomic surveillance of MRSA in the UK reveals transmission patterns in hospitals and the community. *Clin Infect Dis* 65(11):1781-1789. doi:10.1093/cid/cix645
- Snitkin ES, Zelazny AM, Thomas PJ et al (2012) Tracking a hospital outbreak of carbapenem-resistant Klebsiella pneumoniae with whole-genome sequencing. *Sci Transl Med* 4(148):148ra116. doi:10.1126/scitranslmed.3004129
- Napier G, Campino S, Merid Y et al (2020) Robust barcoding and identification of Mycobacterium tuberculosis lineages for epidemiological and clinical studies. *Genome Med* 12(1):114. doi:10.1186/s13073-020-00817-3
- Phelan JE, O'Sullivan DM, Machado D et al (2019) Integrating informatics tools and portable sequencing technology for rapid detection of resistance to anti-tuberculous drugs. *Genome Med* 11:41. doi:10.1186/s13073-019-0650-x
- Hunt M, Bradley P, Lapierre SG et al (2019) Antibiotic resistance prediction for Mycobacterium tuberculosis from genome sequence data with Mykrobe. *Wellcome Open Res* 4:191. doi:10.12688/wellcomeopenres.15603.1
- Walker TM, Ip CLC, Harrell RH et al (2013) Whole-genome sequencing to delineate Mycobacterium tuberculosis outbreaks: a retrospective observational study. *Lancet Infect Dis* 13(2):137-146. doi:10.1016/S1473-3099(12)70277-3
- Eyre DW, Cule ML, Wilson DJ et al (2013) Diverse sources of C. difficile infection identified on whole-genome sequencing. *N Engl J Med* 369(13):1195-1205. doi:10.1056/NEJMoa1216064
- Ondov BD, Treangen TJ, Melsted P et al (2016) Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol* 17:132. doi:10.1186/s13059-016-0997-x
- O'Toole Á, Scher E, Underwood A et al (2021) Assignment of epidemiological lineages in an emerging pandemic using the pangolin tool. *Virus Evol* 7(2):veab064. doi:10.1093/ve/veab064
- Aksamentov I, Roemer C, Hodcroft EB, Neher RA (2021) Nextclade: clade assignment, mutation calling and quality control for viral genomes. *J Open Source Softw* 6(67):3773. doi:10.21105/joss.03773
- Pongmoragot J, Pearson C, Borg ML et al (2024) Comparison of UShER-based and pangoLEARN-based Pangolin lineage assignments for SARS-CoV-2 sequences. *Virus Evol* 10(1):vead085. doi:10.1093/ve/vead085

## Related Skills

- amr-surveillance - Strain context complements AMR; Kleborate integrates typing + AMR for Klebsiella
- transmission-inference - SNP-cluster definition from cgMLST or core-SNP feeds outbreak transmission inference
- phylodynamics - Time-scaled tree from typed isolates for R_e estimation
- variant-surveillance - SARS-CoV-2 Pango / Nextclade lineage assignment overlaps; this skill owns the typing call, variant-surveillance owns longitudinal frequency tracking
- comparative-genomics/pangenome-analysis - Core / accessory genome partitioning underlies cgMLST schema design
- comparative-genomics/whole-genome-alignment - Core-genome alignment for SNP-typing
- variant-calling/vcf-basics - Per-isolate VCF for SNP-typing
- variant-calling/variant-calling - Per-isolate variant calling that feeds cgMLST and SNP-typing
- read-alignment/bwa-alignment - Read mapping upstream of variant calling and snippy
- alignment/multiple-alignment - Multiple sequence alignment for core SNP extraction
- database-access/entrez-fetch - Reference genome retrieval for snippy / Snippy-core
- metagenomics/strain-tracking - Community strain tracking via Kraken2 / StrainPhlAn (NOT isolate-focused)
