---
name: bio-genome-assembly-metagenome-assembly
description: Assembles microbial-community sequencing into metagenome-assembled genomes (MAGs) with metaFlye (ONT), metaSPAdes/MEGAHIT (Illumina), and hifiasm-meta/metaMDBG (PacBio HiFi), then recovers genomes via multi-binner consolidation (MetaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB -> DAS_Tool) and QCs them against MIMAG with CheckM2, GUNC, and GTDB-Tk. Covers why a metagenome is not a genome (uneven coverage, micro-diversity, strain collapse to consensus), differential-coverage binning, co-assembly vs per-sample, the rRNA-operon collapse that fails short-read MAGs, and strain resolution with inStrain. Use when reconstructing genomes from a microbiome, soil, ocean, or gut community, recovering MAGs, or resolving strain-level variation.
tool_type: cli
primary_tool: metaFlye
---

## Version Compatibility

Reference examples tested with: Flye 2.9+, SPAdes 3.15+ (metaSPAdes), MEGAHIT 1.2+, hifiasm-meta 0.3+, metaMDBG 1.0+, MetaBAT2 2.15+, MaxBin 2.2.7+, CONCOCT 1.1+, SemiBin 2.0+, VAMB 4.1+, DAS_Tool 1.1.6+, CheckM2 1.0+, GUNC 1.0+, GTDB-Tk 2.4+, inStrain 1.7+, minimap2 2.26+, samtools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

GTDB-Tk results track the reference-package RELEASE (e.g. R214 vs R220); the DB release MUST match the GTDB-Tk binary or classification silently fails. CheckM2 and GUNC each download their own DIAMOND DB. SemiBin2's pretrained `--environment` models are versioned. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Metagenome Assembly

**"Assemble genomes from my metagenome"** -> Co-assemble a community at uneven, strain-mixed coverage, then bin the contigs into a set of consensus population genomes (MAGs) and QC each against MIMAG. The deliverable is MAGs, not a single assembly.
- CLI: `flye --meta --nano-hq reads.fq` (ONT), `spades.py --meta -1 R1.fq -2 R2.fq` or `megahit -1 R1.fq -2 R2.fq` (Illumina), `hifiasm_meta`/`metaMDBG` (HiFi); then binners -> `DAS_Tool` -> `checkm2 predict` + `gunc run` + `gtdbtk classify_wf`

## The Single Most Important Modern Insight -- A Metagenome Is Not a Genome; the Assembler Cannot Assume Uniform Coverage

Every isolate assembler is built on the premise that the true sequence sits at roughly one depth, so a coverage drop or spike signals a repeat or an error. In a community that premise is false by construction: an abundant species at 500x and a rare one at 3x are both real. Running plain SPAdes/Unicycler or single-genome Flye on a community treats the abundance spread and strain bubbles as errors to "fix" and produces garbage. Use `--meta` modes. Three consequences cascade:

1. **A MAG is a population consensus, not an organism's genome.** Co-occurring strains differing by <1% ANI become bubbles the assembler collapses into one consensus path -- a sequence that may match no actual cell in the sample. A 99%-complete circular MAG is still the consensus of the dominant strain; minority-strain accessory genome is averaged away. Treat every per-strain, per-allele, or pangenome claim from a single consensus MAG as suspect until read-level microdiversity (inStrain) or strain-aware assembly confirms it.
2. **The deliverable is a community of MAGs, not one assembly -- a community has no N50.** N50 is dominated by whichever few abundant genomes assembled well and says nothing about the community; a "better N50" assembly can have recovered fewer genomes. Report MAG count split by MIMAG tier (HQ/medium/low) and community fraction binned. Bigger total assembly size is not better -- it can mean more chimeras and strain-fragmentation.
3. **Modern practice is multi-binner -> consolidate -> CheckM2 + GUNC -> GTDB-Tk.** Never trust one binner; run several (each weights composition vs coverage differently and recovers a partially-different genome set), reconcile with DAS_Tool, then QC every bin with CheckM2 AND GUNC (completeness lies about chimeras) before classifying against the MIMAG 90%/5% bar. HiFi/long reads are the single biggest quality jump: they span the conserved rRNA operons and strain bubbles short reads shred, yielding complete, circular, genuinely-HQ MAGs.

## Assembler Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| metaFlye | Kolmogorov 2020 *Nat Methods* | repeat-graph long-read meta-assembler (`--meta` mandatory) | ONT/CLR community de novo; polish after |
| metaSPAdes | Nurk 2017 *Genome Res* | strain-aware multi-k de Bruijn (`spades.py --meta`) | Illumina, contiguity priority; ONE paired library only |
| MEGAHIT | Li 2015 *Bioinformatics* | succinct de Bruijn, low-memory | huge/complex Illumina co-assemblies, soil; lower contiguity |
| hifiasm-meta | Feng 2022 *Nat Methods* | strain-resolved HiFi string graph | PacBio HiFi communities; keeps strains apart |
| metaMDBG | Benoit 2024 *Nat Biotechnol* | minimizer de Bruijn for HiFi | HiFi; often ~2x the HQ circular MAGs, low RAM |
| OPERA-MS | Bertrand 2019 *Nat Biotechnol* | short-read meta scaffolded by long reads | hybrid short + long |

## Binner Taxonomy

| Binner | Citation | Signal | When |
|--------|----------|--------|------|
| MetaBAT2 | Kang 2019 *PeerJ* | tetranucleotide freq (TNF) + coverage, parameter-free | fast default workhorse (`-m 1500`) |
| MaxBin2 | Wu 2016 *Bioinformatics* | EM over TNF + marker genes + coverage | single/few samples |
| CONCOCT | Alneberg 2014 *Nat Methods* | GMM on composition+coverage of cut-up contigs | many samples; 4-step pipeline, not one command |
| SemiBin2 | Pan 2023 *Bioinformatics* | self-supervised contrastive deep learning | current SOTA; short + long; pretrained env models |
| VAMB | Nissen 2021 *Nat Biotechnol* | variational autoencoder of coabundance + k-mer | multi-sample; separates close strains |
| DAS_Tool | Sieber 2018 *Nat Microbiol* | dereplicate-aggregate-score across binners (consolidation) | ALWAYS run; non-redundant set beats any single binner |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Illumina, complex/huge/soil, low RAM | MEGAHIT `--presets meta-sensitive` | succinct dBG fits in memory; multi-library |
| Illumina, contiguity priority, tractable size | metaSPAdes (`spades.py --meta`) | strain-aware repeat resolution; merge libraries first (one paired lib only) |
| ONT-only community | metaFlye `--meta --nano-hq` -> polish | repeat-graph meta mode; ONT needs polishing -> assembly-polishing |
| PacBio HiFi community | hifiasm-meta or metaMDBG | complete circular strain-resolved MAGs; fixes rRNA collapse |
| Hybrid short + long | OPERA-MS | short-read meta scaffolded with long reads |
| Recover MAGs (any assembly) | >=2-3 binners -> DAS_Tool -> CheckM2 + GUNC -> GTDB-Tk | ensemble beats one binner; chimera + taxonomy gates |
| Only ONE sample | composition-only binning, expect weak bins | differential coverage needs multiple samples; add samples, not tuning |
| Multiple samples available | map ALL samples to each assembly for binning depth | differential-coverage is the strongest binning signal |
| Strain-level question | -> inStrain on reads mapped to MAGs | consensus MAGs blur strains; needs read-level microdiversity |
| Read-based taxonomy / rare biosphere | -> metagenomics/kraken-classification | assembly is blind below the abundance-detection limit |
| Reads not QC'd / host-contaminated | -> long-read-sequencing/long-read-qc | remove host reads vs a T2T reference before assembly |
| MAG contamination forensics | -> contamination-detection | detailed CheckM2/GUNC interpretation |

## metaFlye (ONT / Long Reads)

```bash
flye --meta --nano-hq ont.fastq.gz --out-dir flye_out -t 32
#   --meta        uneven-coverage metagenome mode (REQUIRED for communities)
#   read-type flag (mutually exclusive): --nano-hq (Guppy5+/Q20) | --nano-raw (older) |
#       --pacbio-hifi | --pacbio-raw (CLR)
# outputs: assembly.fasta, assembly_graph.gfa, assembly_info.txt (circularity flag in col 'circ.')
```
ONT contigs are contiguous but error-prone (indels in homopolymers); polish before downstream use (-> assembly-polishing; medaka needs the matching basecaller model). HiFi usually needs no polishing.

## metaSPAdes / MEGAHIT (Illumina)

```bash
# metaSPAdes -- contiguity priority; exactly ONE paired library
spades.py --meta -1 R1.fastq.gz -2 R2.fastq.gz -o spades_out -t 32 -m 500
#   -m memory cap in GB (SPAdes aborts if exceeded); -k auto by default
# outputs: contigs.fasta, scaffolds.fasta

# MEGAHIT -- huge/low-RAM; accepts comma-separated multiple libraries
megahit -1 a1.fq.gz,b1.fq.gz -2 a2.fq.gz,b2.fq.gz -o megahit_out -t 32 \
        --presets meta-sensitive --min-contig-len 1000
#   --presets meta-sensitive | meta-large (huge complex); raise --min-contig-len to ~1000 for binning
```
metaSPAdes `--meta` supports exactly ONE paired-end library -- a real constraint people miss; concatenate libraries first or use MEGAHIT for many. metaSPAdes is heavier on RAM/time and chokes on soil-scale co-assembly; MEGAHIT assembled a 252 Gbp soil set on one node at the cost of somewhat more fragmentation.

## HiFi (the transformative case)

```bash
hifiasm_meta -t 32 -o asm hifi.fastq.gz
awk '/^S/{print ">"$2"\n"$3}' asm.p_ctg.gfa > asm.p_ctg.fa   # GFA -> FASTA

metaMDBG asm --out-dir mdbg_out --in-hifi hifi.fastq.gz --threads 32   # often ~2x HQ circular MAGs
```

## Coverage for Binning, then Bin

```bash
# Map reads back to the assembly -- per sample for differential coverage
minimap2 -ax map-ont -t 32 contigs.fa reads.fq.gz | samtools sort -@ 32 -o s1.sorted.bam -
samtools index s1.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth depth.txt s1.sorted.bam s2.sorted.bam   # all samples

metabat2 -i contigs.fa -a depth.txt -o metabat/bin -m 1500 -t 32   # -m 1500 = min contig (do not go below ~1000)
run_MaxBin.pl -contig contigs.fa -abund abund1.txt -out maxbin/bin -thread 32 -min_contig_length 1000
SemiBin2 single_easy_bin -i contigs.fa -b s1.sorted.bam -o semibin_out   # add --environment human_gut for a pretrained model
```
CONCOCT is a 4-step pipeline (`cut_up_fasta.py` -> `concoct_coverage_table.py` -> `concoct` -> `merge_cutup_clustering.py` -> `extract_fasta_bins.py`), not one command. Differential-coverage binning needs MULTIPLE samples with abundance variation; with one sample binning collapses to weak composition-only signal -- more samples, not more tuning.

## Consolidate (DAS_Tool), then QC

```bash
# Convert each binner's output to contig2bin tables, then aggregate-and-score
Fasta_to_Contig2Bin.sh -i metabat/ -e fa    > metabat.tsv
Fasta_to_Contig2Bin.sh -i maxbin/  -e fasta > maxbin.tsv               # MaxBin emits .fasta
gunzip -k semibin_out/output_bins/*.gz 2>/dev/null || true             # SemiBin2 bins are gzipped
Fasta_to_Contig2Bin.sh -i semibin_out/output_bins/ -e fa > semibin.tsv
DAS_Tool -i metabat.tsv,maxbin.tsv,semibin.tsv -l metabat,maxbin,semibin \
         -c contigs.fa -o dastool/DAS --write_bins -t 32
# DAS_Tool assumes the SAME contig set across binners -- feeding bins from different assemblies is a silent error

checkm2 predict --input dastool/DAS_DASTool_bins/ -x fa --output-directory checkm2_out -t 32
gunc run --input_dir dastool/DAS_DASTool_bins/ --file_suffix .fa --out_dir gunc_out --threads 32
gtdbtk classify_wf --genome_dir dastool/DAS_DASTool_bins/ -x fa --out_dir gtdbtk_out --cpus 32
```
Run CheckM2 AND GUNC: CheckM2 counts marker copy number (completeness/contamination), GUNC tests whether a genome's genes share one lineage (chimerism). A bin made of two half-genomes can score high completeness, low contamination, and still be a chimera -- GUNC is the orthogonal catch. Dereplicate MAGs across samples (dRep ~95% ANI species, ~99% strain) before reporting counts.

## Strain Resolution

```bash
# Consensus MAGs blur strains; recover read-level microdiversity
inStrain profile sample.sorted.bam mags.fa -o instrain_out -p 16 -g genes.fna
inStrain compare -i instrain_A instrain_B -o instrain_compare   # popANI: shared-strain detection across samples
```

## Per-Method Failure Modes

### Isolate assembler on a community
**Trigger:** plain `spades.py`/`flye` (no `--meta`) on community reads. **Mechanism:** uniform-coverage assumption deletes rare-taxon contigs and mis-resolves strain bubbles as errors. **Symptom:** few short contigs, missing abundant taxa. **Fix:** `--meta` mode for every meta-assembler.

### Single-binner pipeline
**Trigger:** "we used SemiBin2 because it's SOTA," one binner. **Mechanism:** each binner recovers a partially-different genome set. **Symptom:** real MAGs left on the table; lower count than peers. **Fix:** run >=2-3 binners -> DAS_Tool consolidation.

### Single-sample differential-coverage expectation
**Trigger:** MetaBAT2/CONCOCT on one sample, surprised bins are bad. **Mechanism:** one coverage value -> composition (TNF) only, which is weak (related genera share TNF). **Symptom:** poor, split bins. **Fix:** more samples with abundance variation, then map all back; do not retune.

### Short-read MAG reported HQ on 90/5 alone
**Trigger:** calling a 95%-complete/2%-contam short-read MAG "high-quality." **Mechanism:** conserved + multi-copy rRNA operon tangles the short-read graph and stays unbinned. **Symptom:** MAG fails MIMAG HQ for missing 16S/23S/5S despite great completeness. **Fix:** check the FULL MIMAG HQ definition (rRNA + tRNA); use HiFi/long reads to span the operon.

### High CheckM2 completeness read as quality
**Trigger:** trusting completeness/contamination alone. **Mechanism:** non-overlapping markers from two half-genomes score high completeness, low contamination. **Symptom:** "clean" MAG that is a chimera. **Fix:** always pair CheckM2 with GUNC -> contamination-detection.

### Strain claims from a consensus MAG
**Trigger:** per-allele/per-strain interpretation of one MAG. **Mechanism:** the assembler collapsed strains to consensus. **Symptom:** strain findings that read-level data contradict. **Fix:** inStrain popANI or strain-aware (hifiasm-meta) assembly.

### Reading MAG richness as community richness
**Trigger:** "200 MAGs cover 60% of reads" treated as the whole community. **Mechanism:** the rare biosphere never crosses the assembly detection limit. **Symptom:** species richness massively undercounted. **Fix:** complement with read-based profiling -> metagenomics/kraken-classification.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| MIMAG high-quality: completeness >90%, contamination <5%, AND 16S/23S/5S rRNA + >=18 tRNAs | Bowers 2017 *Nat Biotechnol* | rRNA criterion is the silent killer; short-read MAGs fail it |
| MIMAG medium-quality: completeness >=50%, contamination <10% | Bowers 2017 | most short-read MAGs land here |
| Min contig length for binning ~1000-1500 bp | MetaBAT2 default 1500 | shorter contigs have unreliable TNF/coverage -> noise/chimeras |
| Differential-coverage binning needs N >= ~3-5 samples with abundance variation | binning practice | one coverage value cannot separate co-abundant genomes |
| Assembly detection limit ~3-5x coverage | de Bruijn requirement | below this the rare biosphere does not assemble at all |
| MAG dereplication 95% ANI (species), 99% (strain) | dRep/ANI convention | collapse redundant per-sample MAGs before counting |
| metaSPAdes input: exactly ONE paired library | `--meta` constraint | merge libraries or use MEGAHIT for many |
| N50 / contiguity | not a community metric | report MAG count + MIMAG tiers instead |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Few short contigs, missing abundant taxa | isolate mode on a community | add `--meta` |
| metaSPAdes rejects multiple libraries | `--meta` supports one paired library | concatenate, or use MEGAHIT |
| metaSPAdes aborts / out of memory | RAM cap hit on a complex co-assembly | raise `-m`, or switch to MEGAHIT |
| Poor bins from one sample | no differential-coverage signal | add samples; map all back for depth |
| "HQ MAG" lacks rRNA | rRNA-operon collapse in short reads | full MIMAG check; HiFi/long reads |
| Clean CheckM2 but suspect bin | chimera invisible to marker counts | run GUNC alongside CheckM2 |
| GTDB-Tk classify_wf errors | DB release != binary version | match the GTDB reference-package release |
| Strain finding not reproducible | consensus MAG blurs strains | inStrain popANI / strain-aware assembly |

## References

- Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. 2017. metaSPAdes: a new versatile metagenomic assembler. *Genome Res* 27:824-834.
- Li D, Liu CM, Luo R, Sadakane K, Lam TW. 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 31:1674-1676.
- Kolmogorov M, et al. 2020. metaFlye: scalable long-read metagenome assembly using repeat graphs. *Nat Methods* 17:1103-1110.
- Feng X, Cheng H, Portik D, Li H. 2022. Metagenome assembly of high-fidelity long reads with hifiasm-meta. *Nat Methods* 19:671-674.
- Benoit G, et al. 2024. High-quality metagenome assembly from long accurate reads with metaMDBG. *Nat Biotechnol* 42:1378-1383.
- Bertrand D, et al. 2019. Hybrid metagenomic assembly enables high-resolution analysis of resistance determinants and mobile elements in human microbiomes (OPERA-MS). *Nat Biotechnol* 37:937-944.
- Kang DD, et al. 2019. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 7:e7359.
- Wu YW, Simmons BA, Singer SW. 2016. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics* 32:605-607.
- Alneberg J, et al. 2014. Binning metagenomic contigs by coverage and composition (CONCOCT). *Nat Methods* 11:1144-1146.
- Pan S, Zhao XM, Coelho LP. 2023. SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. *Bioinformatics* 39:i21-i29.
- Nissen JN, et al. 2021. Improved metagenome binning and assembly using deep variational autoencoders (VAMB). *Nat Biotechnol* 39:555-560.
- Sieber CMK, et al. 2018. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy (DAS_Tool). *Nat Microbiol* 3:836-843.
- Chklovski A, et al. 2023. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 20:1203-1212.
- Orakov A, et al. 2021. GUNC: detection of chimerism and contamination in prokaryotic genomes. *Genome Biol* 22:178.
- Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. 2020. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics* 36:1925-1927.
- Bowers RM, et al. 2017. Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nat Biotechnol* 35:725-731.
- Olm MR, et al. 2021. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. *Nat Biotechnol* 39:727-736.

## Related Skills

- contamination-detection - CheckM2/GUNC interpretation and chimerism forensics for MAGs
- assembly-qc - Isolate-assembly QC; the uniform-coverage assumption metagenomes abandon
- assembly-polishing - Polish ONT/CLR meta-contigs before binning (HiFi usually needs none)
- metagenomics/kraken-classification - Read-based taxonomy; recovers the rare biosphere assembly cannot
- metagenomics/abundance-estimation - Community abundance downstream of recovered MAGs
- metagenomics/functional-profiling - Functional potential complementary to genome recovery
- long-read-sequencing/long-read-qc - Read-level QC and host removal before assembly
