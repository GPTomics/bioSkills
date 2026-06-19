---
name: bio-read-qc-contamination-screening
description: Detects contamination in sequencing reads - cross-species (FastQ Screen, Kraken2), vector/PhiX/adapter, rRNA, and same-species cross-sample/index-hopping and sample swaps (SNP fingerprints via verifyBamID2/NGSCheckMate/somalier). Use when suspecting cross-contamination, PDX host reads, microbial carry-over, or sample swaps, and to decide whether to report, filter, or align to a combined reference. For deep taxonomic profiling use metagenomics/kraken-classification.
tool_type: cli
primary_tool: fastq_screen
---

## Version Compatibility

Reference examples tested with: FastQ Screen 0.15+, Bowtie2 2.5+, Kraken2 2.1+, BBTools 39.0+, MultiQC 1.21+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Contamination Screening -- a species screen cannot see a same-species swap

Screen reads against a genome panel (FastQ Screen / Kraken2) for foreign ORGANISMS, and against SNP fingerprints for foreign or wrong INDIVIDUALS.

**"Check my reads for contamination"** -> Map a subsample against multiple references and/or fingerprint sample identity to find foreign DNA and mislabels.
- CLI: `fastq_screen --conf fastq_screen.conf sample.fastq.gz` (cross-species)
- CLI: `verifyBamID2 --SVDPrefix resource --BamFile sample.bam` (same-species contamination)

Scope: this skill OWNS contamination detection and the report-vs-filter decision. Deep taxonomic profiling/abundance -> metagenomics/kraken-classification, metagenomics/metaphlan-profiling. rRNA depletion as an RNA prep metric -> read-qc/rnaseq-qc. OUT OF SCOPE: adapter removal (read-qc/adapter-trimming).

## The Single Most Important Modern Insight

1. **A SPECIES screen answers "what ORGANISMS are here?"; a SNP FINGERPRINT answers "WHOSE DNA is this, and is it a mixture?" -- and these are orthogonal.** A species screen is structurally BLIND to same-species cross-sample contamination and sample swaps: a human-A + human-B mixture, or a mislabeled human file, produces a perfectly clean single-species profile. Most pipelines run only a species screen and declare the data clean. Any human or single-species-cohort pipeline needs BOTH a taxonomic screen AND a SNP-fingerprint identity/contamination check (verifyBamID2, NGSCheckMate, somalier, conpair).

2. **Index hopping is the same-species contamination that lives inside one run, and unique dual indexing (UDI) is the only clean fix.** On patterned flowcells (HiSeq X/4000, NovaSeq) ExAmp chemistry lets a free index adapter tag a fragment from another sample, spreading ~0.1-2% of reads into incorrect samples. Irrelevant for germline common variants, CATASTROPHIC for low-VAF work (ctDNA, single-cell, somatic) where hopped reads look like phantom low-frequency variants. Combinatorial indexing cannot detect it; UDI (a unique i7 AND i5 per sample) lets the demultiplexer drop impossible index pairs. UDI is effectively mandatory for ctDNA/plasma, single-cell, and low-input somatic.

3. **Default to SCREEN-AND-REPORT, not filter -- removing reads biases composition.** A species screen is a QC gate; if a contaminant is low, aligning to the correct reference simply will not place the foreign reads. Filter only a specific, named contaminant (PhiX before assembly, adapters before alignment) with a precise k-mer remover (BBDuk `ref=phix`), never "remove anything that hits the screen" (that also discards conserved rRNA/mito reads that belong to the sample). For PDX, align to a COMBINED human+mouse reference and keep human-assigned reads, OR use a dedicated post-alignment classifier (XenofilteR benchmarks above Xenome for variant false-positive rate); both beat hard pre-filtering on one genome, which mis-assigns conserved-region reads.

Deeper trap: reference-genome contamination corrupts the screen itself. A "human" hit can be bacterial sequence mis-deposited inside the human assembly (Conterminator found >2M contaminated GenBank entries). No `--confidence` setting fixes a wrong database; trust a deconned/curated DB and treat surprising single-source hits as DB artifacts until ruled out.

## The Contamination Taxonomy -- five classes, five fixes

| Class | What it is | Detect with | Fix |
|-------|-----------|-------------|-----|
| Cross-species | Mouse in human PDX; bacteria in culture | FastQ Screen, Kraken2/Bracken, sourmash, Xenome/XenofilteR | Combined-reference alignment; k-mer bin (report, do not blindly remove) |
| Cross-sample / index hopping | Same-species reads on the wrong sample (INVISIBLE to species screens) | verifyBamID2, NGSCheckMate, somalier, conpair | UDI at prep; drop impossible index pairs at demux |
| Vector / PhiX / adapter | Spike-in, cloning vector, linkers | UniVec/VecScreen; FastQ Screen adapter DB; BBDuk `ref=phix`/`ref=adapters` | k-mer trim/filter; upstream library QC |
| rRNA over-representation | Library-prep failure, not contamination | SortMeRNA, ribodetector | Re-prep / better depletion; filter only to recover depth |
| Cell-line: mycoplasma / misID | Mollicutes infection; HeLa cross-contamination | Kraken2 for Mollicutes; STR profiling for line identity | Clear culture or discard; STR-authenticate |

A pipeline that runs FastQ Screen and stops has checked exactly one of five boxes.

## Tool Taxonomy

| Tool | Mechanism | When |
|------|-----------|------|
| FastQ Screen | Map a subsample to a genome panel; classify hit categories | Cross-species QC gate; the workhorse screen |
| Kraken2 + Bracken | Exact-k-mer minimizer LCA classification; Bracken re-estimates abundance | Read-level taxonomy; many possible contaminants |
| BBSplit / BBDuk | k-mer binning / named-contaminant k-mer removal | Decontamination of a NAMED contaminant (PhiX, adapters) |
| Xenome / XenofilteR | Classify reads human vs mouse (k-mer / dual-alignment) | PDX host-graft disambiguation |
| sourmash | MinHash/FracMinHash containment sketches | Fast low-memory "what is in here?" screen |
| verifyBamID2 | Per-sample within-species contamination from population SNP AFs | Same-species contamination level (FREEMIX) |
| NGSCheckMate / somalier | SNP-fingerprint identity / relatedness | Sample swaps, tumor-normal pairing, longitudinal identity |
| conpair | Tumor-normal concordance + independent contamination | Matched T/N pairs |

## Decision Tree by Scenario

| Question | Use | Why |
|----------|-----|-----|
| Is a foreign ORGANISM present? | FastQ Screen or Kraken2 | Maps reads to species references |
| Is this the right INDIVIDUAL / one person? | NGSCheckMate / somalier | SNP fingerprint, species-screen-blind |
| What is the contamination LEVEL (human)? | verifyBamID2 (FREEMIX) | Estimates mixture fraction from SNP AFs |
| Tumor-normal pair: matched and clean? | conpair | Concordance + per-sample contamination |
| PDX host vs graft | Combined reference or a benchmarked classifier | XenofilteR > Xenome for SNV FP rate; both beat hard pre-filtering |
| Remove a NAMED contaminant | BBDuk `ref=...` | Precise k-mer removal, not "hits the screen" |
| Strip HUMAN reads before public deposition (non-human library) | hostile / NCBI sra-human-scrubber (HRRT) | Deposition compliance, not a QC gate; a masked T2T reference avoids stripping conserved microbial regions |

Default when uncertain: FastQ Screen as the QC gate for organisms, PLUS a SNP-fingerprint check (somalier/NGSCheckMate) for any human cohort.

## FastQ Screen

Maps a SUBSAMPLE (`--subset`, default 100000) with bowtie2 reporting >1 alignment, then classifies each read across the panel. Read the bar chart, not just "% mapped": contamination concentrates in `One_hit_one_genome` of an UNEXPECTED genome; homology (rRNA, mito, conserved loci) spreads into the `*_multiple_genomes` categories; high `Hit_no_genomes` means adapter dimer, a missing reference, or a novel organism (a diagnostic, not a verdict).

```bash
# Config: aligner binary + DATABASE lines (bowtie2 index prefixes)
cat > fastq_screen.conf <<'EOF'
BOWTIE2  /usr/local/bin/bowtie2
THREADS  8
DATABASE  Human  /refs/GRCh38_bt2/GRCh38
DATABASE  Mouse  /refs/GRCm39_bt2/GRCm39
DATABASE  Ecoli  /refs/Ecoli_bt2/Ecoli
DATABASE  PhiX   /refs/phix_bt2/phix
DATABASE  rRNA   /refs/rRNA_bt2/rRNA
EOF

fastq_screen --conf fastq_screen.conf --threads 8 --outdir screen/ *.fastq.gz
multiqc screen/                      # MultiQC parses *_screen.txt across samples

# Tag every read with a per-genome status, then extract a subset by pattern
fastq_screen --conf fastq_screen.conf --tag --filter 10000 sample.fastq.gz   # maps only to genome 1
fastq_screen --conf fastq_screen.conf --nohits sample.fastq.gz               # reads hitting nothing
```

`--filter` digits (one per genome, config order): 0=no map, 1=unique, 2=multi, 3=maps, 4=pass 0 or 1, 5=pass 0 or 2, -=ignore. `--subset 0` screens the whole file; `--bisulfite` uses Bismark.

## Kraken2 + Bracken (read-level taxonomy)

Default `--confidence 0.0` over-reports a long tail of spurious low-abundance species (a few shared k-mers suffice); raise to 0.05-0.1 and keep `--minimum-hit-groups 2` (or 3 for custom DBs). Kraken2 gives CLASSIFICATION; Bracken redistributes higher-rank reads to species for ABUNDANCE.

```bash
kraken2 --db /db/k2_standard --threads 8 --confidence 0.1 --paired \
        --report sample.kreport --use-names R1.fq.gz R2.fq.gz > sample.kraken
bracken -d /db/k2_standard -i sample.kreport -o sample.bracken -r 150 -l S
```

## Same-species: SNP fingerprints and index hopping

```bash
# verifyBamID2: FREEMIX = contamination fraction (action threshold ~0.02);
# FREEMIX~0 with CHIPMIX~1 indicates a SWAP, not contamination.
# --SVDPrefix points to the panel resource that ships with verifyBamID2 (resource/1000g.phase3...).
verifyBamID2 --SVDPrefix /res/1000g.phase3.100k.b38.vcf.gz.dat --BamFile sample.bam --Reference ref.fa

# somalier: extract genome sketches, then relate to find swaps / identity across a cohort.
# --sites = somalier's released sites.<build>.vcf.gz (github releases), not a custom panel.
somalier extract -d sites/ --sites sites.hg38.vcf.gz -f ref.fa sample.bam
somalier relate sites/*.somalier            # off-diagonal identity flags swaps

# conpair (tumor-normal): concordance + independent per-sample contamination
```

Index hopping is mitigated at demultiplexing with UDI (drop impossible i7,i5 pairs); residual contamination is then quantified by the SNP-fingerprint tools above.

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| "Single species, data is clean" but a swap is suspected | Species screen is blind to same-species swaps | Run somalier / NGSCheckMate on SNP fingerprints |
| Phantom low-VAF variants in ctDNA/single-cell | Index hopping on patterned flowcell | Use UDI; quantify residual with verifyBamID2/conpair |
| Kraken2 reports dozens of trace species | Default confidence 0.0 over-reports | Raise `--confidence` to 0.05-0.1; raise hit-groups |
| A "human" Kraken hit on a microbial isolate | Reference/DB contamination | Use a deconned DB; treat as artifact until confirmed |
| Filtering "contaminant" reads skews composition | Removed conserved rRNA/mito too | Remove a NAMED contaminant with BBDuk, not "hits the screen" |
| PDX human counts look biased | Hard pre-filtering of ambiguous reads | Align to combined human+mouse reference instead |

## References

Wingett SW, Andrews S. 2018. FastQ Screen: a tool for multi-genome mapping and quality control. F1000Research 7:1338.
Wood DE, Lu J, Langmead B. 2019. Improved metagenomic analysis with Kraken 2. Genome Biology 20:257.
Lu J, Breitwieser FP, Thielen P, Salzberg SL. 2017. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3:e104.
Steinegger M, Salzberg SL. 2020. Terminating contamination: large-scale search identifies more than 2,000,000 contaminated entries in GenBank. Genome Biology 21:115.
Conway T, Wazny J, Bromage A, et al. 2012. Xenome - a tool for classifying reads from xenograft samples. Bioinformatics 28(12):i172-i178.
Costello M, Fleharty M, Abreu J, et al. 2018. Characterization and remediation of sample index swaps by non-redundant dual indexing. BMC Genomics 19:332.
Zhang F, Flickinger M, Taliun SAG, et al. 2020. Ancestry-agnostic estimation of DNA sample contamination from sequence reads. Genome Research 30(2):185-194.
Lee S, Lee S, Ouellette S, Park WY, Lee EA, Park PJ. 2017. NGSCheckMate: software for validating sample identity in next-generation sequencing studies within and across data types. Nucleic Acids Research 45(11):e103.
Pedersen BS, Bhetariya PJ, Brown J, et al. 2020. Somalier: rapid relatedness estimation for cancer and germline studies using efficient genome sketches. Genome Medicine 12:62.

## Related Skills

read-qc/quality-reports - Bimodal GC and overrepresented sequences flag contamination
read-qc/adapter-trimming - Remove adapter contamination
read-qc/rnaseq-qc - rRNA fraction as a prep-efficiency metric
metagenomics/kraken-classification - Deeper taxonomic classification and profiling
variant-calling/joint-calling - Where SNP-fingerprint sample swaps do the most damage
