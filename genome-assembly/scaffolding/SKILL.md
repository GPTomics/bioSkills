---
name: bio-genome-assembly-scaffolding
description: Orders and orients assembled contigs into chromosome-scale scaffolds from long-range linking data, inserting N-gap spacers (adds no sequence). Covers Hi-C/Omni-C scaffolding (YaHS, SALSA2, 3D-DNA/Juicer), Hi-C read-mapping prerequisites (map each end separately, no mate rescue, dedup, enzyme-aware), reading the contact map for misjoins/inversions/false-duplications, manual curation in Juicebox/PretextView (the VGP/DToL standard), reference-guided scaffolding (RagTag) and its karyotype-erasure hazard, genetic-map (ALLMAPS) and Bionano optical-map integration, chimera-breaking before scaffolding, gap-filling, and telomere/contig-vs-scaffold-N50 QC (tidk). Use when turning contigs into chromosomes with Hi-C, integrating a linkage map or optical map, choosing a scaffolder by available linking data, or judging whether a chromosome-scale assembly is trustworthy.
tool_type: cli
primary_tool: YaHS
---

## Version Compatibility

Reference examples tested with: YaHS 1.2+, SALSA2 2.3+, juicer_tools 1.22/2.0+, bwa 0.7.17+, chromap 0.2+, samtools 1.19+, RagTag 2.1+, tidk 0.2.3+, seqkit 2.6+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

hifiasm output filenames and YaHS/SALSA2 enzyme and resolution defaults have changed across releases; confirm output names and `-e`/`-r` behaviour against the installed version. The Hi-C-mapping flag set (`bwa mem -5SP` vs the Arima/VGP per-end pipeline) varies by vendor (Arima, Dovetail/Omni-C, Phase Genomics) - check the current mapping guide before pasting flags. If a command errors, introspect the tool and adapt rather than retrying.

# Genome Scaffolding

**"Turn my contigs into chromosomes"** -> Order and orient contigs by long-range linking signal, insert N-gap spacers, then curate the contact map - scaffolding adds no sequence, so the output is a draft chromosome structure, not finished sequence.
- CLI: `yahs contigs.fa hic_to_contigs.bam` (Hi-C default), `run_pipeline.py -a contigs.fa -b sorted.bed -m yes` (SALSA2), `ragtag.py scaffold ref.fa contigs.fa` (reference-guided), `allmaps path` (linkage maps)

## The Single Most Important Modern Insight -- Automated Scaffolding Produces a Draft; the Genome Is Trustworthy Only After Manual Curation of the Contact Map

YaHS finishes in minutes and emits a file literally named `_scaffolds_final.agp`. It is **not final**. Hi-C scaffolding is a contact-frequency *inference* - it orders/orients contigs by the polymer-physics signal that contact frequency decays with 1D distance - and it is error-prone exactly at the joins. Every reference-grade pipeline (VGP, Darwin Tree of Life, Earth BioGenome) treats the automated AGP as a *first pass* that a human curator then corrects in PretextView or Juicebox: breaking misjoins, flipping inversions, removing false duplications, assigning chromosomes by eye. Howe 2021 quantified the hidden labor across 111 VGP/DToL assemblies: on average **221 interventions per Gb (67 breaks, 105 joins, 49 false-duplication removals)**. Three load-bearing consequences:

1. **The contact map is the QC, not a decoration.** A correct chromosome shows one bright diagonal with smooth off-diagonal decay; off-diagonal blocks, anti-diagonal "bowties," and bleeding between chromosomes are errors to inspect and break (see Reading the Contact Map). An assembly paper whose methods say "scaffolded with [tool]" but never mention curation/Juicebox/PretextView is shipping a draft - downgrade any claim that depends on large-scale order (synteny, fusions, structural variants).

2. **Scaffold N50 is not contig N50, and conflating them is the classic reviewer catch.** Scaffolding inserts runs of **N** (estimated, often arbitrary length) and adds zero sequence. A genome can post "scaffold N50 = 60 Mb, chromosome-scale!" while its contig N50 is 200 kb - the contiguity is borne by gaps, not finished sequence. Always report contig N50, scaffold N50, gap count, and total N bases separately.

3. **Reordering cannot rescue a bad contig.** Scaffolders take contigs as given; a chimeric contig drags the wrong sequence into the wrong chromosome. Chimeras must be **broken before scaffolding** (see Chimera-Breaking).

## Linking-Data Modality Taxonomy

| Modality | How it links | Scale | Status |
|----------|--------------|-------|--------|
| Hi-C / Omni-C | proximity ligation; contact frequency ~ 1/(1D distance) | chromosome-scale | Dominant modern method. Omni-C is enzyme-free/sequence-agnostic; Arima uses fixed enzyme motifs |
| Bionano optical maps (DLS) | labeled motif patterns aligned to in-silico contig maps | megabase | Resolves large SVs and bridges complex repeats; separate instrument + library (cost decision) |
| 10x linked reads | barcoded short reads from one long molecule | ~10-100 kb | DISCONTINUED by 10x (2020). Legacy data only; Tigmint/ARCS still run, no new data |
| Genetic / linkage maps | marker recombination order = chromosome order | chromosome-scale, low resolution | ALLMAPS integrates multiple maps; orthogonal validation of Hi-C |
| Reference (homology) | align contigs to a related genome, copy its order | as good as synteny | RagTag. Fast, but imposes the reference's karyotype - see hazard below |
| Long reads as linkers | k-mer pairs spanning junctions | read length | LINKS; mostly superseded by assembling with the long reads directly |
| Mate-pair / jumping | large-insert paired short reads | 2-40 kb | Legacy (SSPACE); chimeric-insert artifacts; never recommend today |

## Decision Tree by Available Linking Data

| Available linking data | Use | Why |
|-------------|-----|-----|
| Hi-C/Omni-C, want fast contiguous default | YaHS | community default; fast, high N90, AGP + Juicebox outputs in one command |
| Hi-C, want graph-aware conservative joins + input-error correction | SALSA2 `-m yes` (`-g graph.gfa`) | iterative; breaks chimeric input; uses assembly graph to avoid orientation errors |
| Hi-C, interactive curation central to workflow | 3D-DNA + Juicer + Juicebox (JBAT) | the Aiden-lab `.assembly` <-> Juicebox round-trip is built around hand-editing |
| Hi-C, vertebrate/reference-grade | YaHS or SALSA2 -> PretextView/JBAT curation | VGP/DToL standard: automate then curate |
| Diploid, Hi-C, want both haplotypes as chromosomes | hifiasm `--h1/--h2` to PHASE first, then YaHS to SCAFFOLD each haplotype | same Hi-C, two jobs - phasing picks the homolog, scaffolding orders along it (see Phasing vs Scaffolding) |
| Closely related reference, karyotype NOT a question | RagTag `correct` then `scaffold` | homology ordering in minutes - but never if karyotype is the biology |
| One or more genetic/linkage maps | ALLMAPS | integrates multiple maps; robust to marker errors; validates/anchors Hi-C |
| Bionano optical maps + sequence assembly | Bionano Solve hybrid scaffold (before Hi-C) | megabase maps bridge complex repeats and large SVs |
| 10x linked reads (legacy data) | Tigmint (break) -> ARCS/ARKS (link) | break misassemblies first; no new 10x data exists |
| Contigs not yet QC'd / haplotigs present | -> assembly-qc, purge_dups | scaffolding haplotigs strings them up as fake chromosomes |
| Hi-C contact map for TADs/loops, not chromosomes | -> hi-c-analysis/matrix-operations | different use of the same assay |

**Modern reference-grade recipe (vertebrate/eukaryote):** HiFi -> hifiasm contigs -> (optional Bionano hybrid scaffold) -> Hi-C scaffold (YaHS/SALSA2) -> **manual curation in PretextView/JBAT** -> gap-fill (TGS-GapCloser) -> QC (tidk telomeres, contact-map diagonal, contig-vs-scaffold N50).

## Hi-C Read Mapping (the step people botch)

Hi-C reads are **not** a normal paired-end library: the two ends come from *different* genomic loci ligated together, so standard PE proper-pair/insert-size logic mis-handles them.

```bash
# Map each end SEPARATELY with no mate rescue/pairing; -5 reports the 5' (junction) portion as primary.
bwa index contigs.fa
bwa mem -5SP -T0 -t 16 contigs.fa hic_R1.fq.gz hic_R2.fq.gz | \
    samtools view -@ 8 -b - > aligned.bam
# Mandatory: mark/remove PCR + optical duplicates (Hi-C is duplicate-rich) before scaffolding.
samtools sort -@ 8 -n aligned.bam | samtools fixmate -m - - | \
    samtools sort -@ 8 - | samtools markdup -@ 8 - hic_to_contigs.bam
samtools index hic_to_contigs.bam
```

- **Each end separately, no mate rescue** (`-5SP`): `-5` = 5' portion of a chimeric junction read as primary; `-S`/`-P` skip mate rescue and pairing. The Arima/VGP per-end pipeline (`bwa mem` on R1 and R2 independently, `filter_five_end.pl`, `two_read_bam_combiner.pl`) is the higher-fidelity equivalent.
- **MAPQ filter** at scaffolding (YaHS `-q`, SALSA2 reads MAPQ) to drop repeat-ambiguous reads.
- **Enzyme awareness:** Arima/Dovetail-DpnII cut at fixed motifs (`GATC`); the scaffolder uses `-e` to model legitimate read starts. **Omni-C / DNase Hi-C is enzyme-free** -> omit `-e` in YaHS, use `-e DNASE` in SALSA2.
- **chromap** (Zhang 2021 *Nat Commun* 12:6566) is the fast modern alternative: `chromap --preset hic -r contigs.fa -x index -1 R1 -2 R2 ...` aligns + dedups Hi-C ~10x faster; YaHS accepts its output.

## YaHS (the current default)

```bash
samtools faidx contigs.fa
# Enzyme: -e GATC (DpnII/Dovetail), -e GATC,GANTC (Arima 2-enzyme); OMIT -e for Omni-C/enzyme-free.
yahs -e GATC -o out contigs.fa hic_to_contigs.bam
# Outputs: out_scaffolds_final.agp + out_scaffolds_final.fa, intermediate out_rNN.agp, out.bin
# Contig error-correction is ON by default (--no-contig-ec to disable; usually a mistake to disable).

# Juicebox prep for curation (.hic + .assembly to load in JBAT):
# NOTE: `juicer` here is the small utility BUNDLED with YaHS (operates on the .bin), NOT Aiden-lab Juicer;
# `juicer_tools.jar` below IS the separate Aiden-lab jar. Do not confuse the two.
juicer pre -a -o out_JBAT out.bin out_scaffolds_final.agp contigs.fa.fai 2>tmp_assembly.log
java -Xmx48G -jar juicer_tools.jar pre out_JBAT.txt out_JBAT.hic <(cat tmp_assembly.log | grep PRE_C_SIZE | awk '{print $2" "$3}')
# After hand-editing in Juicebox, export out_JBAT.review.assembly, then:
juicer post -o out_curated out_JBAT.review.assembly out_JBAT.liftover.agp contigs.fa
```

YaHS also accepts a BED (`bamToBed`) or PA5/pairs input. The `_scaffolds_final.agp` is the editable, reviewable object - curation edits the AGP, then regenerates FASTA. `--telo-motif` lets YaHS use telomere signal during scaffolding.

## SALSA2 (graph-aware, iterative, breaks chimeras)

```bash
bamToBed -i aligned.bam > alignment.bed
sort -k4 alignment.bed > sorted.bed                 # MUST sort by read name (column 4)
python run_pipeline.py -a contigs.fa -l contigs.fa.fai -b sorted.bed \
       -e GATC -o salsa_out -m yes -g contigs_graph.gfa -p yes
# -e DNASE for Omni-C/enzyme-free; -m yes enables input-error (chimera) correction; -p yes writes AGP+FASTA per iteration
```

The signature `-m yes` step breaks input contigs where the Hi-C signal contradicts the assembler's join *before* scaffolding. Assembly headers must not contain `:`. `-i` sets iterations (default 3), `-c` min contig length (default 1000).

## Reference-Guided Scaffolding (RagTag) -- Power and Peril

```bash
ragtag.py correct ref.fa contigs.fa                 # break query at putative misassemblies vs reference FIRST
ragtag.py scaffold ref.fa ragtag_output/ragtag.correct.fasta -t 16   # order/orient by homology (minimap2 default)
# Outputs ragtag.scaffold.agp + .fasta; inserts 100 bp N gaps by default; -C collapses unplaced into chr0.
```

**The peril (load-bearing):** RagTag orders contigs *to match the reference*, so by construction the output looks like the reference. Every real biological rearrangement that distinguishes the target organism - a chromosome fusion, fission, translocation, or inversion polymorphism - is silently **erased** and replaced by the reference's structure. **Never reference-scaffold a genome whose karyotype or large-scale structure is a biological question** - the pipeline "discovers" the reference's karyotype because it was imposed (this has happened in the literature). Acceptable uses: a structurally conserved close relative needing only a coordinate system for gene content; a scaffold *hypothesis* to then test/correct against Hi-C (`ragtag.py merge -b hic.bam` lets Hi-C arbitrate). On a genome with interesting structure the honest move is Hi-C + curation, full stop.

## Genetic/Linkage Maps (ALLMAPS) and Optical Maps (Bionano)

ALLMAPS (Tang 2015 *Genome Biol* 16:3) integrates one or more genetic maps to order/orient scaffolds: `allmaps merge map1.csv map2.csv -o maps.bed` then `allmaps path maps.bed contigs.fa`. It is the gold standard for *validating* chromosome assignment and for species where Hi-C is hard; weight multiple maps by quality. Bionano Solve hybrid scaffolding (no journal paper; cite the *Bionano Solve Theory of Operation: Hybrid Scaffold* technical document) aligns DLS optical maps to in-silico contig maps, resolves conflicts, and merges into hybrid scaffolds - run it *before* Hi-C on large/repetitive genomes where megabase maps bridge repeat arrays Hi-C fumbles. The decision is economic (separate instrument + library), not purely technical.

## Reading the Contact Map (the skill nobody writes down)

The map *displays its own errors* to a trained eye - render with PretextMap -> PretextView (interactive, emits a curated AGP), PretextSnapshot (static PNG), or Juicebox/JBAT.

- **Correct chromosome:** one bright diagonal, contacts decaying smoothly off-diagonal; inter-chromosome space dim and uniform.
- **Misjoin:** the diagonal breaks - signal jumps off-diagonal into a separate block at the junction (two diagonals stitched at a corner). Break it.
- **Inversion:** an **anti-diagonal "bowtie"/butterfly** - the gradient runs backwards through the segment. The most common curation fix once recognized; flip it.
- **Translocation / wrong assignment:** off-diagonal bright blocks bleeding between what should be separate chromosomes.
- **False duplication (leaked haplotig):** a faint off-diagonal stripe parallel to the diagonal contacting the same neighborhood. Remove it (49 of the 221 interventions/Gb).
- **Hardest judgment:** the map is noisy near the diagonal even when correct. A *real* weak join is faint but coherent with the right decay shape; *mapping noise* is structureless speckle. Treat ambiguous near-diagonal signal as a flag to inspect, not a number to trust.

## Chimera-Breaking, Gap-Filling, and Phasing-vs-Scaffolding

- **Break, then scaffold:** SALSA2 `-m yes` and YaHS default contig-EC break chimeric input where Hi-C contradicts the join; Tigmint (Jackman 2018 *BMC Bioinformatics* 19:393) breaks before ARCS for linked reads; RagTag `correct` for the reference-based version. Disabling these to "preserve the assembly" is a common self-inflicted wound.
- **Gap-filling is a SEPARATE downstream step** that replaces N-runs with real sequence by bridging long reads across the gap (TGS-GapCloser; LR_Gapcloser). Gaps sit at repeats (that is why the assembler stopped there), so a read can bridge into the *wrong* repeat copy and insert locally-plausible but globally-wrong sequence - worse than an honest N. Gap-fill *after* curation, sanity-check filled lengths against the gap estimate, and report which gaps were closed vs left.
- **Hi-C-for-scaffolding vs Hi-C-for-phasing:** the same library, two jobs. Scaffolding asks *where on the chromosome* (orders collapsed contigs); phasing asks *which homolog* (hifiasm `--h1/--h2` during assembly). When handed "Hi-C, make chromosomes," the first question is *haploid scaffolding or diploid phasing+scaffolding?* - the tools and failure modes differ entirely.

## Per-Method Failure Modes

### Shipping the automated AGP as final
**Trigger:** publishing `_scaffolds_final` without curation. **Mechanism:** automated joins are inferences error-prone at junctions (~221 edits/Gb in VGP/DToL). **Symptom:** misjoins/inversions in the contact map; later-discovered wrong synteny/fusions. **Fix:** curate in PretextView/JBAT before any large-scale-order claim.

### Conflating scaffold N50 with contig N50
**Trigger:** reporting only scaffold N50. **Mechanism:** scaffolding adds N-gaps, not sequence. **Symptom:** spectacular scaffold N50 over a small contig N50 - contiguity is nitrogen. **Fix:** report both N50s, gap count, total N bases separately.

### Reference-scaffolding a structurally interesting genome
**Trigger:** RagTag on a non-model genome whose karyotype is studied. **Mechanism:** ordering to the reference imposes the reference's structure. **Symptom:** fusions/translocations/inversions silently erased; circular "discovery" of the reference karyotype. **Fix:** Hi-C + curation; use RagTag only as a hypothesis to test.

### Scaffolding before breaking chimeras
**Trigger:** scaffolding contigs as given, contig-EC off. **Mechanism:** a chimeric contig drags region B into region A's chromosome; every join to it is wrong. **Symptom:** misjoins the contact map blames on the scaffolder. **Fix:** keep SALSA2 `-m yes` / YaHS contig-EC on; Tigmint for linked reads.

### Treating Hi-C reads as a normal PE library
**Trigger:** `bwa mem` without `-5SP`, no dedup, or proper-pair logic. **Mechanism:** the two ends are different loci; mate rescue and insert-size filtering corrupt placement. **Symptom:** sparse/noisy contact map, weak joins. **Fix:** map each end separately (`-5SP` or per-end pipeline), dedup, MAPQ filter, enzyme-aware.

### Gap-filling across a repeat
**Trigger:** running a gap-closer over the whole assembly to "finish" it. **Mechanism:** a long read bridges into a different copy of the repeat the gap sits in. **Symptom:** filled length far from the gap estimate; locally clean, globally mis-assembled. **Fix:** gap-fill after curation, sanity-check lengths, report closed-vs-left.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| ~221 curation interventions/Gb (67 breaks, 105 joins, 49 false-dup) | Howe 2021 *GigaScience* | expect the automated draft to be substantially editable, not finished |
| Hi-C ~50-100M valid pairs (vertebrate); rule-of-thumb ~1x coverage per ~1 Mb | VGP practice (~approx) | too little -> weak/missing joins, sparse noisy map; depends on genome size/repeats |
| Scaffold N50 / contig N50 > ~10x | curator reflex | contiguity is gap-borne; flag genome as contiguous-on-paper but unfinished |
| T2T chromosome = telomere at BOTH ends + zero internal gaps | T2T convention | tidk both-end recovery with internal Ns = chromosome-scale but not finished |
| RagTag default gap = 100 bp N | Alonge 2022 | placeholder length is arbitrary; never a real distance estimate |
| MAPQ filter on Hi-C reads (e.g. YaHS `-q 10`) | mapping practice | drop repeat-ambiguous reads that create spurious joins |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `_scaffolds_final` looks chromosome-scale but synteny is wrong | shipped uncurated draft | curate the contact map in PretextView/JBAT |
| Huge scaffold N50, small contig N50 | contiguity borne by N-gaps | report both; do not conflate |
| Reference-scaffolded genome "has" the reference's karyotype | RagTag imposed the structure | re-scaffold with Hi-C + curation |
| Sparse/noisy contact map, few joins | Hi-C mapped as normal PE / under-sequenced | remap with `-5SP`+dedup; add Hi-C coverage |
| YaHS makes spurious joins | enzyme/MAPQ misconfigured; chimeric contigs | set correct `-e` (omit for Omni-C), raise `-q`, keep contig-EC on |
| SALSA2 errors on input | `:` in FASTA headers or BED not name-sorted | rename headers; `sort -k4` the BED |
| Filled gaps far from estimated length | gap-filler bridged wrong repeat copy | gap-fill post-curation; sanity-check lengths; report closed-vs-left |

## References

- Zhou C, McCarthy SA, Durbin R. 2023. YaHS: yet another Hi-C scaffolding tool. *Bioinformatics* 39:btac808.
- Ghurye J, et al. 2019. Integrating Hi-C links with assembly graphs for chromosome-scale assembly (SALSA2). *PLoS Comput Biol* 15:e1007273.
- Dudchenko O, et al. 2017. De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds (3D-DNA). *Science* 356:92-95.
- Durand NC, et al. 2016. Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. *Cell Syst* 3:95-98.
- Alonge M, et al. 2022. Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. *Genome Biol* 23:258.
- Tang H, et al. 2015. ALLMAPS: robust scaffold ordering based on multiple maps. *Genome Biol* 16:3.
- Yeo S, et al. 2018. ARCS: scaffolding genome drafts with linked reads. *Bioinformatics* 34:725-731.
- Jackman SD, et al. 2018. Tigmint: correcting assembly errors using linked reads from large molecules. *BMC Bioinformatics* 19:393.
- Zhang H, et al. 2021. Fast alignment and preprocessing of chromatin profiles with chromap. *Nat Commun* 12:6566.
- Rhie A, et al. 2021. Towards complete and error-free genome assemblies of all vertebrate species (VGP). *Nature* 592:737-746.
- Howe K, et al. 2021. Significantly improving the quality of genome assemblies through curation. *GigaScience* 10:giaa153.
- Brown MR, Gonzalez de la Rosa PM, Blaxter M. 2025. tidk: a toolkit to rapidly identify telomeric repeats from genomic datasets. *Bioinformatics* 41:btaf049.
- Bionano Genomics. Bionano Solve Theory of Operation: Hybrid Scaffold (technical document) - hybrid scaffolding has no journal paper.

## Related Skills

- long-read-assembly - Produces the contigs this skill orders into chromosomes
- hifi-assembly - hifiasm phases haplotypes (Hi-C `--h1/--h2`) before each is scaffolded
- assembly-polishing - Polish contigs before scaffolding; gaps remain N until gap-filling
- assembly-qc - Contig-vs-scaffold N50, BUSCO, and Merqury QV on the scaffolded result
- hi-c-analysis/hic-data-io - Hi-C read/pairs handling feeding the scaffolder
- comparative-genomics/synteny-analysis - Validate scaffold order/orientation against a related genome
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> curate -> QC
