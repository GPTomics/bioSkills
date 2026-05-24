---
name: bio-comparative-genomics-pangenome-analysis
description: Build and analyze pangenomes for prokaryotes (Panaroo, PPanGGOLiN, PEPPAN, GET_HOMOLOGUES, anvi'o pangenomics) and eukaryotes (Minigraph-Cactus, PGGB, vg pangenome graphs). Implement Tettelin core/accessory/cloud genome decomposition (Tettelin 2005), Heap's law open/closed pangenome modeling, gene presence/absence GWAS (Scoary, pyseer), pangenome graph variant calling (vg, PanGenie), and structural-variation graph indexing. Use when assembling species- or genus-level pan-gene catalogs, separating core from accessory/shell/cloud genes, testing gene-content associations with phenotypes, building pangenome graphs from haplotype-resolved assemblies, calling SVs from pangenome graphs, or selecting between bacterial-pangenome and eukaryotic-pangenome workflows.
tool_type: cli
primary_tool: Panaroo
---

## Version Compatibility

Reference examples tested with: Panaroo 1.5.1+ (Tonkin-Hill 2020 Genome Biol 21:180), PPanGGOLiN 2.2.0+ (Gautreau 2020 PLoS Comp Biol 16:e1007732), PEPPAN 1.0.5+ (Zhou 2020 GR 30:1667), GET_HOMOLOGUES 25102023+, anvi'o 8.0+ (Eren 2021 Nat Microbiol 6:3), Minigraph-Cactus (Hickey 2024 Nat Biotech 42:663; bundled with Cactus 2.5+), PGGB 0.7.5+ (Garrison 2024 Nat Methods 21:2008), vg 1.59.0+ (Sirén et al 2024 NAR Genom Bioinform 6:lqae001), PanGenie 3.1.0+ (Ebler 2022 Nat Genet 54:518), PGR-TK 0.3.6+ (Chin 2023 Nat Methods 20:1290; cschin/pgr-tk; repo archived April 2026 transitioning to PANGEA), PANGEA (in development by DGI / Diploid Genomics as PGR-TK's successor for pangenome graph exploration + analysis -- check https://github.com/cschin/pgr-tk for current repository pointer), Bakta 1.10.4+ (annotation for input), Roary 3.13.0+ (DEPRECATED; use Panaroo), Scoary 1.6.16+, pyseer 1.3.11+, BUSCO 5.7+, FastTree 2.1.11+, RAxML-NG 1.2+. Python 3.10+ required for Panaroo / PPanGGOLiN.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `panaroo --version`; `ppanggolin --version`; `peppan --help`; `cactus-pangenome --help`; `pggb --version`; `vg version`
- Python: `pip show panaroo ppanggolin`

If code throws `Bakta annotation incompatible`, `GFA file inconsistent`, `vg index version mismatch`, the bacterial pangenome ecosystem expects consistent annotation; re-annotate all input genomes with the same tool and version before pangenome analysis.

# Pangenome Analysis

**"What genes are universal vs accessory across this set of genomes?"** -> The pangenome is the union of all genes across a sampled group; the **Tettelin partition** (Tettelin 2005 PNAS 102:13950) splits it into core (universal), shell (in many but not all), cloud (rare), and species-specific (private) genes. The fundamental dichotomy is **bacterial pangenome** (clusters genes into orthogroups; Panaroo / PPanGGOLiN / PEPPAN for compact genomes) vs **eukaryotic pangenome** (graph-based; Minigraph-Cactus / PGGB / vg for haplotype-resolved sequences). The choice depends on what's being represented: bacterial pangenome captures gene-content variation in a species/genus; eukaryotic pangenome graph captures haplotype-level structural and sequence variation. **Roary (Page 2015) is now deprecated** in favor of Panaroo, which handles annotation-error noise that previously inflated bacterial pangenomes by 30-50%.

- CLI: `panaroo -i annotated_gffs/ -o panaroo_out --clean-mode strict --remove-invalid-genes` -- bacterial pangenome
- CLI: `ppanggolin workflow --fasta fasta_list.tsv --output ppanggolin_out` -- partitioned bacterial pangenome
- CLI: `peppan -i genomes.gff -o peppan_out` -- genus-scale bacterial pangenome
- CLI: `cactus-pangenome jobStore seqFile.txt --reference name --vcf --gfa` -- Minigraph-Cactus pangenome
- CLI: `pggb -i genomes.fa.gz -n 90 -t 32 -o pggb_out` -- PGGB pangenome graph
- CLI: `vg autoindex --workflow giraffe -r ref.fa -v variants.vcf.gz` -- vg pangenome indexing

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| Panaroo (Tonkin-Hill 2020 GB 21:180) | Graph-based ortholog clustering + annotation-error correction | Core, shell, accessory pangenome with cleaned annotations | Best for clonal bacteria (Mtb 413-genome benchmark; Tonkin-Hill 2020) | Slow for > 10000 genomes; assumes Prokka/Bakta input |
| PPanGGOLiN (Gautreau 2020 PLoS CB 16:e1007732) | Hidden Markov partition: persistent/shell/cloud | Partitioned pangenome with HMM-based class assignment | Scales to many genomes; interpretable partitions | Probabilistic class boundaries differ from strict Tettelin |
| PEPPAN (Zhou 2020 GR 30:1667) | Bacterial pangenome for diverse genera | Pan + core genomes from 1000s of genomes | Designed for high diversity (whole genus) | Slower than newer alternatives at small scales |
| Roary (Page 2015; DEPRECATED) | Original bacterial pangenome | Same outputs as Panaroo | Legacy; widely cited | Inflates accessory by ~30-50% due to annotation-error tolerance; use Panaroo |
| GET_HOMOLOGUES (Contreras-Moreira 2013) | Multi-algorithm consensus (OrthoMCL, BDBH, COG) | Consensus pangenome | Cross-validates across methods | Slower; multi-program output integration |
| anvi'o pangenomics (Eren 2021 Nat Microbiol 6:3) | Interactive pangenome with metadata | Visual pangenome browsing + integration | Standard for interactive microbial pangenome | Less automated; manual curation expected |
| Minigraph-Cactus (Hickey 2024 Nat Biotech 42:663) | Cactus base-level + minigraph SV-graph integration | Pangenome graph (GFA, VCF, GBZ) | Production-grade for HPRC-scale haplotypes | Requires reference; designed for intra-species pangenome |
| PGGB (Garrison 2024 Nat Methods 21:2008) | All-vs-all wfmash + seqwish | Pangenome graph (GFA) | Modern reference-free graph; HPRC-validated | Computationally heavy at > 100 genomes |
| vg pangenome (Sirén et al 2024 NAR Genom Bioinform 6:lqae001) | Pangenome graph indexing + Giraffe / GiraffeY mapping | Mapped reads to graph + variant calling | vg ecosystem standard for graph-based variant calling | Setup complex; learning curve |
| PanGenie (Ebler 2022 Nat Genet 54:518) | Pangenome-graph-based genotyping | SV genotype calls | Efficient genotyping from short reads via graph | Requires pre-built pangenome graph |
| PGR-TK (Chin 2023 Nat Methods 20:1290; cschin/pgr-tk) | Minimizer Anchored Pangenome (MAP) graph + principal bundle decomposition | Multiscale pangenome graph; bundle SVGs; AGC-backed sequence db | Designed for repetitive / clinically-relevant genes (MHC class II, DAZ1-4, OPN1LW/OPN1MW); decomposes tangled graph into interpretable bundles; complements Minigraph-Cactus by exposing fine-grained allele structure | Repo archived April 2026 -> PANGEA succession; pinned to Peregrine-assembler-derived workflow; not a drop-in for variant-calling pipelines |
| PANGEA (in development by DGI / Diploid Genomics; succeeds PGR-TK; check cschin/pgr-tk for current pointer) | Next-generation MAP-graph framework | Same conceptual outputs as PGR-TK with modernized API | Active development 2026+; expected to add tighter integration with HPRC / T2T workflows | API surface in flux; pin specific version when scripting |
| Heaps law / Tettelin (Tettelin 2005 PNAS 102:13950) | Statistical model of pangenome openness | Open / closed pangenome classification | Foundational framework | Class boundaries depend on sampling |
| Scoary (Brynildsrud 2016) | Pan-GWAS on gene presence/absence | Phenotype-gene associations | Standard bacterial pan-GWAS tool | Limited to binary phenotypes |
| pyseer (Lees 2018 Bioinformatics 34:4310) | Continuous + binary phenotype association on k-mers/genes | Pan-GWAS with k-mer / gene-content units | More flexible than Scoary | Computational cost |
| pirate (Bayliss 2019) | Bacterial pangenome from multiple methods | Cross-method consensus | Alternative to GET_HOMOLOGUES | Less popular now |

Methodology evolves; verify the current Panaroo and PPanGGOLiN manuals + the 2024-2025 microbial pangenome reviews (Tonkin-Hill 2024 Nat Comm). The HPRC draft pangenome (Liao 2023 Nature 617:312) sets the modern eukaryotic pangenome standard; for bacterial work, Panaroo + PPanGGOLiN is the standard combination.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Bacterial strain set (5-1000 genomes) of one species | Panaroo + PPanGGOLiN | Cross-validation; Panaroo's annotation-cleaning + PPanGGOLiN's partition |
| Bacterial genus-level pangenome (> 1000 genomes) | PEPPAN | Designed for high genus-level diversity |
| Mycobacterium tuberculosis (clonal) | Panaroo | Tonkin-Hill 2020 benchmark; clonal pangenomes |
| E. coli (highly diverse) | PPanGGOLiN or PEPPAN | High accessory diversity |
| Eukaryotic intra-species pangenome (e.g. human, soybean) | Minigraph-Cactus or PGGB | Graph-based; SV-aware |
| HPRC-style 90 haplotype graph | Minigraph-Cactus | Production-grade for HPRC scale |
| Reference-free eukaryotic pangenome | PGGB | All-vs-all alignment-free graph |
| Pangenome graph for variant calling | vg autoindex -> vg giraffe | Standard graph-aligner ecosystem |
| Bacterial pan-GWAS for phenotype | Panaroo + Scoary or pyseer | Pangene matrix from Panaroo; pan-GWAS tool |
| Visualize pangenome interactively | anvi'o pangenomics workflow | Standard for interactive analysis |
| Distinguish core / shell / cloud genes | PPanGGOLiN (HMM-partitioned) or Tettelin manual partition on Panaroo output | Standard Tettelin framework |
| Open vs closed pangenome (Heaps law) | wgd v2 statistical fit OR custom mclust on Panaroo output | Tettelin 2005 framework |
| Detect HGT-acquired accessory genes | Cross-reference with [[hgt-detection]] | Pangenome + phylogeny |
| Eukaryotic structural-variation indexing | Minigraph-Cactus -> vg + PanGenie | SV-aware genotyping pipeline |
| Bacterial functional pangenome | Panaroo + eggNOG-mapper + KEGG | Functional annotation |
| Pangenome-aware reference for read alignment | vg giraffe with pangenome graph | Reduces reference bias |
| Genome-graph-based fine-mapping | vg + GraphAligner or vg giraffe | SV-aware short-read alignment |
| Repetitive / clinically relevant gene (MHC class II, DAZ1-4, OPN1LW/OPN1MW) | PGR-TK MAP graph + principal bundle decomposition | Built for tangled repeat graphs; bundle decomposition reveals haplotype-allele structure that linear refs collapse |
| Next-gen pangenome graph exploration (2026+) | PANGEA (PGR-TK successor, in development by DGI / Diploid Genomics) | Modernized successor to PGR-TK; check cschin/pgr-tk pointer for current repo |
| HLA / KIR / immune-locus pangenome | PGR-TK + manual bundle inspection | Standard tools collapse repeat alleles; PGR-TK's bundle decomposition preserves them |

## Per-Tool Failure Modes

### Annotation heterogeneity inflating accessory genome

**Trigger:** Running Panaroo on Prokka- vs Bakta- vs RefSeq- annotated genomes mixed.

**Mechanism:** Different annotation tools predict different gene boundaries; the same gene is annotated slightly differently across tools, appearing as separate orthogroups. Roary's tolerance of these differences inflated bacterial pangenome accessory by 30-50% (Tonkin-Hill 2020 GB 21:180). Panaroo's graph-based correction reduces this but cannot eliminate it.

**Symptom:** Per-strain "accessory" gene count is inflated relative to known biology; comparison to a single-pipeline reference reveals 30-50% spurious gene-content differences.

**Fix:** Re-annotate ALL genomes with one pipeline (currently Bakta 1.10.4+ for bacteria; Bakta is GenBank-compliant and faster than Prokka). Use `panaroo --clean-mode strict --remove-invalid-genes` to apply graph cleaning. Document annotation pipeline + version in methods.

### Tettelin partition class boundary artifacts

**Trigger:** Reporting "core genome" vs "accessory" as percent-of-strains thresholds (e.g. 99% = core).

**Mechanism:** Tettelin 2005 used 100%-presence = core; pragmatic studies use 95%-99%. The class boundary is arbitrary; small variation in the threshold dramatically changes core/accessory ratio.

**Symptom:** Core genome size varies 10-30% depending on whether threshold is 95% or 99%.

**Fix:** Report core/shell/cloud at multiple thresholds; PPanGGOLiN's HMM partition is more principled but still has tunable parameters. Standard reporting: core = present in >=95% (or >=99%); shell = 15-95%; cloud = < 15%. Document threshold.

### Roary's annotation-error noise (DEPRECATED tool)

**Trigger:** Using Roary in new analyses.

**Mechanism:** Roary tolerates annotation errors (low-identity matches; protein-vs-DNA matches), producing thousands of artifactual accessory genes. Panaroo's introduction (Tonkin-Hill 2020) demonstrated this by re-analyzing 413 Mtb genomes and finding Panaroo's accessory genome was 30-50% smaller than Roary's.

**Symptom:** Roary output has many "lineage-specific genes" with poor evidence (single-strain hits, short proteins, no functional annotation).

**Fix:** Migrate to Panaroo. Panaroo can read Roary's input format; for legacy projects, re-run with Panaroo and compare. Panaroo's `--clean-mode strict` enforces strict graph-based quality control.

### Pangenome graph reference bias (eukaryote)

**Trigger:** Using Minigraph-Cactus with a single reference; calling variants against the "reference" path.

**Mechanism:** Minigraph-Cactus is reference-anchored; the chosen reference appears throughout the graph as a privileged path. Variants are called relative to the reference path; non-reference haplotypes are under-represented in the variant calling.

**Symptom:** Variant call density on non-reference haplotypes is lower than on reference; allele frequencies skewed toward reference.

**Fix:** Use PGGB (reference-free) for less reference-biased analysis. Alternatively, treat the reference choice as a methodological parameter and document. For HPRC, multiple references can be used and results pooled.

### Heaps law misapplied

**Trigger:** Concluding "open pangenome" from a single Heaps-law fit on insufficient data.

**Mechanism:** Heaps law parameter alpha distinguishes open (alpha < 1) from closed (alpha > 1) pangenome; estimation requires sampling many genomes. Few genomes give unstable estimates.

**Symptom:** Heaps-law alpha varies by > 0.2 across resampling; conclusions about pangenome openness flip.

**Fix:** Require >= 50 genomes (preferably > 100) for Heaps-law estimation; report 95% CI from resampling. Tettelin 2005 demonstrated open Streptococcus pyogenes; Vernikos 2015 reviews open vs closed across taxa.

### PGGB memory exhaustion at many genomes

**Trigger:** Running PGGB with > 30 large eukaryotic genomes on a single node.

**Mechanism:** PGGB's wfmash all-vs-all step has O(N^2) memory pattern; large eukaryotic genomes (> 1 Gb) make memory prohibitive for > 30 input genomes.

**Symptom:** PGGB OOMs at the wfmash stage; cluster job killed by OOM-killer.

**Fix:** Use Minigraph-Cactus for > 30 genomes (it scales better); or split PGGB into chromosomes/regions. PGGB recommendation is <= 20 large genomes per run.

### vg index version mismatch breaking giraffe

**Trigger:** Pre-built vg index used with a different vg version for read alignment.

**Mechanism:** vg index format evolved; pre-built indexes from one version may not be compatible with another.

**Symptom:** vg giraffe fails with "index version" error.

**Fix:** Rebuild vg index with current version; or pin vg version for an analysis. Future vg releases promise backward compatibility but verify.

### PanGenie genotype false positives in repetitive regions

**Trigger:** PanGenie on highly repetitive regions (centromeres, segmental duplications).

**Mechanism:** Pangenome-graph-based genotyping requires unique paths in the graph; highly repetitive regions create graph-spaghetti paths that are difficult to genotype reliably.

**Symptom:** PanGenie calls many heterozygous SVs in known-repetitive regions; quality scores low.

**Fix:** Restrict PanGenie to non-repetitive regions; combine with traditional read-based SV callers (DELLY, Manta) for repeat regions. The HPRC paper documents this limitation (Liao 2023).

### Bacterial pangenome with frequent gene-content recombination

**Trigger:** Building a pangenome of a species with extensive recombination (e.g. Neisseria, Streptococcus pneumoniae).

**Mechanism:** Frequent recombination breaks the "vertical inheritance" assumption underlying ortholog clustering; the same gene appears in many phylogenetic positions across strains, complicating orthology and inflating accessory genome.

**Symptom:** Phylogenetic trees from core genome are unstable; per-gene trees show extensive incongruence; pangenome accessory genome appears artificially large.

**Fix:** Use ClonalFrameML (Didelot 2015 PLoS Comp Biol 11:e1004041) to identify recombinant regions; mask them before pangenome analysis. Restrict core genome analysis to non-recombinant regions.

### Annotation density variation across genomes

**Trigger:** Mixing well-annotated reference genomes with newly assembled, draft-annotation genomes.

**Mechanism:** Draft annotations miss small genes, pseudogenes, and lineage-specific genes; well-annotated genomes have these. Comparing them inflates "accessory" in draft genomes.

**Symptom:** Draft genomes have 200-500 fewer accessory genes than expected; per-genome BUSCO completeness > annotation completeness.

**Fix:** Re-annotate all genomes consistently with Bakta + Prodigal; document BUSCO completeness for each. Exclude genomes with > 5% lower BUSCO than median.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Core genome threshold | >=95% (relaxed) to 100% (strict) of strains | Tettelin 2005; pragmatic |
| Shell genome | 15-95% (or 5-95% per PPanGGOLiN) | PPanGGOLiN docs |
| Cloud genome | < 15% of strains | Tettelin 2005 |
| Heaps law alpha (open) | < 1 | Tettelin 2005 |
| Heaps law alpha (closed) | > 1 | Tettelin 2005 |
| Minimum genomes for Heaps law fit | >= 50; >= 100 preferred | Vernikos 2015 |
| Panaroo gene cluster identity | >=70% (default); stricter for clonal | Panaroo defaults |
| PPanGGOLiN coverage | 80% gene-length coverage in clustering | Default |
| PEPPAN BLAT threshold | identity >= 70% | Zhou 2020 |
| Mycobacterium tuberculosis core | ~3500-3700 genes (Tonkin-Hill 2020 Mtb benchmark) | Bench results |
| E. coli pangenome (open) | core ~2400; pangenome >15000 (Touchon 2016) | Reference |
| Plasmodium falciparum core | ~5300 genes (eukaryotic prokaryote-like) | Otto 2018 |
| Minimum strains for bacterial pangenome | >= 5; >= 20 for shell/cloud meaningful | Empirical |
| HPRC pangenome size | 90 haplotypes, ~6.4M variants | Liao 2023 |
| PGGB recommended max genomes | <= 20 large eukaryotic; 100+ for compact | Garrison 2024 |
| PanGenie minimum k-mer | k = 31 default | Ebler 2022 |
| vg index Haplotype Sampling | --haplotype-sampling YES for multi-pop graph | Sirén et al 2024 NAR Genom Bioinform 6:lqae001 |
| Scoary pan-GWAS p-value threshold | Bonferroni-corrected p < 0.05 | Brynildsrud 2016 |
| pyseer continuous-trait power | requires > 1000 isolates for solid signal | Lees 2018 |
| anvi'o pangenome minimum | 5+ genomes for non-trivial visualization | Eren 2021 |

## Panaroo Bacterial Pangenome Workflow

**Goal:** Construct a high-quality bacterial pangenome with annotation-error correction.

**Approach:** Annotate genomes consistently with Bakta -> run Panaroo strict mode -> partition with PPanGGOLiN.

```bash
# 1. Annotate all genomes with Bakta (consistent annotation)
mkdir -p annotated
for fa in genomes/*.fa; do
    name=$(basename $fa .fa)
    bakta --db /path/to/bakta-db --threads 16 \
        --output annotated/${name} --prefix $name \
        --genus Escherichia --species coli \
        $fa
done

# 2. Run Panaroo
panaroo -i annotated/*.gff -o panaroo_out -t 16 \
    --clean-mode strict --remove-invalid-genes

# 3. Extract pangenome matrix
# panaroo_out/gene_presence_absence.csv      strains x genes matrix
# panaroo_out/core_gene_alignment.aln        core gene MSA for phylogeny
# panaroo_out/pan_genome_reference.fa        consensus pangenome sequence

# 4. Tettelin partition (custom)
python tettelin_partition.py \
    --presence panaroo_out/gene_presence_absence.csv \
    --core-threshold 0.99 --shell-threshold 0.15 \
    --output panaroo_out/tettelin_classification.tsv

# 5. PPanGGOLiN HMM partition (alternative)
# PPanGGOLiN expects a TSV index: `genome_name<TAB>path/to.gff3` per row
for f in annotated/*.gff3; do
    printf "%s\t%s\n" "$(basename "$f" .gff3)" "$(realpath "$f")"
done > gff_list.tsv
ppanggolin all --anno gff_list.tsv -o ppanggolin_out --threads 16
# Output: ppanggolin_out/pangenome.h5 (HDF5 with HMM-partitioned genes)
```

```python
'''Tettelin core/shell/cloud partition from Panaroo gene presence/absence matrix.'''
import pandas as pd


def tettelin_partition(presence_matrix, core_threshold=0.99,
                       shell_threshold=0.15):
    '''Returns DataFrame[gene_name] -> Tettelin class.'''
    # presence_matrix: rows = genes, cols = strains, 0/1 entries
    n_strains = presence_matrix.shape[1]
    fraction = presence_matrix.sum(axis=1) / n_strains
    classes = pd.cut(fraction,
                     bins=[-0.01, shell_threshold, core_threshold, 1.01],
                     labels=['cloud', 'shell', 'core'])
    return pd.DataFrame({'fraction': fraction, 'class': classes})


def heaps_law(presence_matrix, n_iters=100):
    '''Estimate Heaps law alpha from genome sampling order.'''
    import numpy as np
    n_strains = presence_matrix.shape[1]
    pan_sizes = []
    for _ in range(n_iters):
        order = np.random.permutation(n_strains)
        pan = set()
        sizes = []
        for i in order:
            genes_in_i = presence_matrix.iloc[:, i] == 1
            pan.update(genes_in_i.index[genes_in_i].tolist())
            sizes.append(len(pan))
        pan_sizes.append(sizes)
    pan_array = np.array(pan_sizes)  # n_iters x n_strains
    n_sampled = np.arange(1, n_strains + 1)
    # Fit log-log
    mean_pan = pan_array.mean(axis=0)
    log_n = np.log(n_sampled)
    log_pan = np.log(mean_pan)
    alpha = np.polyfit(log_n, log_pan, 1)[0]
    return alpha
```

## Minigraph-Cactus for Eukaryotic Pangenome

**Goal:** Build pangenome graph from haplotype-resolved assemblies.

**Approach:** Provide reference + haplotypes -> Minigraph-Cactus produces GFA + VCF + GBZ for downstream genotyping.

```bash
# Prepare seqFile (Cactus convention)
cat > pangenome_seqs.txt << 'EOF'
GRCh38      GRCh38.fa
HG002.hap1  HG002.hap1.fa
HG002.hap2  HG002.hap2.fa
HG003.hap1  HG003.hap1.fa
HG003.hap2  HG003.hap2.fa
EOF

cactus-pangenome jobStore_path pangenome_seqs.txt \
    --outDir hprc_pangenome \
    --outName hprc_pangenome \
    --reference GRCh38 \
    --vcf \
    --gfa \
    --gbz \
    --indexCores 32 \
    --mapCores 32

# Outputs:
#   hprc_pangenome/hprc_pangenome.full.hal      Full Cactus HAL
#   hprc_pangenome/hprc_pangenome.gfa.gz        Graph Fragment Assembly format
#   hprc_pangenome/hprc_pangenome.vcf.gz        Short variants relative to GRCh38
#   hprc_pangenome/hprc_pangenome.gbz           GBZ compressed graph
#   hprc_pangenome/hprc_pangenome.giraffe.gbz   Giraffe-indexed graph for mapping
```

## vg Pangenome Genotyping

**Goal:** Genotype short reads against a pre-built pangenome graph.

**Approach:** Pre-built pangenome -> vg autoindex -> vg giraffe (fast mapping) -> vg call variants.

```bash
# Pre-build index
vg autoindex --workflow giraffe --threads 16 \
    --ref-graph hprc_pangenome.gfa.gz \
    --output hprc_index

# Map reads
vg giraffe \
    --gbz-name hprc_index.giraffe.gbz \
    --dist-name hprc_index.dist \
    --minimizer-name hprc_index.min \
    --fastq-in sample.R1.fq.gz \
    --fastq-in sample.R2.fq.gz \
    --output-format GAM \
    --threads 16 \
    > sample.gam

# Pack alignment information
vg pack -x hprc_index.giraffe.gbz -g sample.gam -o sample.pack

# Call variants
vg call hprc_index.giraffe.gbz -k sample.pack -a > sample.vcf
```

## Pan-GWAS with Scoary

**Goal:** Identify gene presence/absence associated with a phenotype.

**Approach:** Panaroo presence/absence matrix + phenotype file -> Scoary -> phenotype-gene associations.

```bash
# Run Scoary
scoary \
    --gene-presence-absence panaroo_out/gene_presence_absence.csv \
    --traits phenotypes.tsv \
    --output scoary_out \
    --threads 16 \
    --upgma-tree

# Output:
#   scoary_out/*_results.csv   per-trait gene associations
```

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Panaroo accessory >> PPanGGOLiN accessory | Panaroo "accessory" includes singleton; PPanGGOLiN "cloud" is a separate class | Compare strict mode definitions; Panaroo + PPanGGOLiN cross-validate |
| Roary accessory >> Panaroo accessory | Roary annotation-error inflation | Trust Panaroo; Roary deprecated |
| PEPPAN core != Panaroo core | Different clustering thresholds | Panaroo for clonal; PEPPAN for genus-scale; consistent within method |
| Minigraph-Cactus VCF and Cactus pairwise differ | Pangenome integrates SV; pairwise is direct | Minigraph-Cactus for variant-aware pangenome |
| PGGB and Minigraph-Cactus disagree on graph topology | PGGB reference-free; MC reference-anchored | Both valid; report both for transparency |
| Heaps-law alpha differs across resampling | Stochasticity; insufficient sampling | Require > 100 genomes; report 95% CI |
| PanGenie genotype contradicts read-based SV caller | Repetitive region (graph-spaghetti) | Trust read-based for repetitive; PanGenie for unique regions |
| Bakta and Prokka annotation give different gene counts | Different gene-prediction defaults | Use Bakta (GenBank-compliant); document |
| Scoary and pyseer disagree on top genes | Different statistical assumptions | Cross-validate; trust consensus |

**Operational rule for publication:** Bacterial pangenome uses Panaroo + PPanGGOLiN cross-validation; eukaryotic pangenome uses Minigraph-Cactus (HPRC scale) or PGGB (reference-free). Document annotation pipeline + version; report core/shell/cloud at multiple thresholds; verify Heaps-law on > 100 genomes for openness claims.

## Cohort Gotchas

- **Endosymbiont genomes (Buchnera, Wolbachia, mitochondria):** core genome is dominant; accessory is minimal; pangenome analysis less informative
- **Hypothetical proteins:** unknown function genes dominate "accessory" in non-model species; functional analysis needed
- **Phage genes:** prophages contribute to accessory; mask prophages for "core species genes"
- **Mobile genetic elements:** plasmids/IS elements appear in accessory; tag with mobileOG-db
- **Highly recombinogenic species (Neisseria, S. pneumoniae):** core genome unstable; restrict to non-recombinant regions
- **Plasmids:** chromosome vs plasmid distinction matters for pangenome; classify with PlasmidFinder
- **Polyploid eukaryotic pangenome:** subgenomes must be assigned first (see [[whole-genome-duplication]])
- **Pangenome graph for SV detection:** PanGenie + graph aligners; reference-anchored variant callers miss SVs

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Annotation pipeline?" | Bakta 1.10+ on all genomes; consistent settings; BUSCO completeness reported per strain |
| "Why Panaroo over Roary?" | Panaroo's annotation-error correction reduces accessory inflation 30-50% (Tonkin-Hill 2020 Mtb benchmark) |
| "Tettelin thresholds?" | Reported at multiple thresholds (95%, 99%); PPanGGOLiN HMM partition as cross-validation |
| "Heaps law inference?" | >= 100 genomes; resampling 95% CI reported |
| "Recombination?" | ClonalFrameML applied; recombinant regions masked or analyzed separately |
| "Why Minigraph-Cactus?" | HPRC standard; production-grade for haplotype-resolved pangenomes |
| "Reference bias?" | Reference choice documented; PGGB cross-validation for reference-free comparison |
| "Pan-GWAS multiple testing?" | Bonferroni-corrected across genes; or pyseer with k-mer-based |
| "Open vs closed pangenome?" | Heaps-law alpha reported with CI; openness claim conditional on alpha < 1 |
| "Annotation density consistency?" | BUSCO completeness verified per strain; outliers excluded |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Panaroo "input format error" | GFF3 missing required fields | Verify GFF3 from Bakta has correct attribute fields |
| PPanGGOLiN HDF5 unreadable | Version mismatch | Pin PPanGGOLiN version; rebuild |
| Roary used (legacy script) | Roary deprecated | Migrate to Panaroo |
| PEPPAN OOM | Too many genomes | Reduce to representative subset; or use PEPPAN with chunking |
| Cactus pangenome unrelated to expectation | Wrong seqFile syntax | Tabs not spaces; correct file paths |
| PGGB wfmash hangs | Too many large genomes | Reduce to <= 20 eukaryotic genomes |
| vg autoindex memory error | Insufficient RAM | Increase to 200+ GB; or split by chromosome |
| PanGenie genotype empty | Index mismatch between graph + reads | Re-build vg index with same vg version |
| Scoary "no significant traits" | Few strains or low effect | Increase strains; verify phenotype variance |
| pyseer p-values uniform | Population structure inflation | Use --lmm flag with kinship matrix |
| anvi'o display error | Database version mismatch | Re-create with current anvi'o version |
| Minigraph-Cactus runs but VCF empty | All inputs identical | Verify inputs differ |
| Bakta annotation gives 0 genes | Reference data not configured | Set BAKTA_DB env or use `--db /path` |

## Tool Installation Notes

```bash
# Bacterial pangenome
conda install -c bioconda panaroo ppanggolin peppan get_homologues anvio

# Eukaryotic pangenome
conda install -c bioconda cactus pggb vg pangenie

# PGR-TK (repeat-rich / clinical gene focus)
conda install -c bioconda pgr-tk
# Or: cargo install pgrtk; or Docker quay.io/cschin/pgr-tk
# Note: cschin/pgr-tk archived April 2026; transitioning to PANGEA (developed by DGI / Diploid Genomics; check upstream repo for pointer)

# Annotation
conda install -c bioconda bakta prokka

# Pan-GWAS
conda install -c bioconda scoary pyseer

# Recombination
conda install -c bioconda clonalframeml

# Mobile elements
git clone https://github.com/clb21565/mobileOG-db

# QC
conda install -c bioconda busco compleasm
```

For HPRC-scale eukaryotic pangenome, use cluster with >= 500 GB RAM and HPC scheduler integration via Toil (see [[whole-genome-alignment]]).

## References

- Tettelin H et al 2005 PNAS 102:13950 (core/accessory pangenome framework)
- Tonkin-Hill G et al 2020 Genome Biol 21:180 (Panaroo)
- Gautreau G et al 2020 PLoS Comp Biol 16:e1007732 (PPanGGOLiN)
- Zhou Z et al 2020 Genome Res 30:1667 (PEPPAN)
- Page AJ et al 2015 Bioinformatics 31:3691 (Roary; DEPRECATED)
- Contreras-Moreira B & Vinuesa P 2013 BMC Bioinf 14:64 (GET_HOMOLOGUES)
- Eren AM et al 2021 Nat Microbiol 6:3 (anvi'o pangenomics)
- Hickey G et al 2024 Nat Biotech 42:663 (Minigraph-Cactus)
- Garrison E et al 2024 Nat Methods 21:2008 (PGGB)
- Sirén J et al 2024 NAR Genom Bioinform 6:lqae001 (vg pangenome update)
- Ebler J et al 2022 Nat Genet 54:518 (PanGenie)
- Liao W-W et al 2023 Nature 617:312 (HPRC draft pangenome)
- Brynildsrud O et al 2016 Genome Biol 17:238 (Scoary)
- Lees JA et al 2018 Bioinformatics 34:4310 (pyseer)
- Didelot X & Wilson DJ 2015 PLoS Comp Biol 11:e1004041 (ClonalFrameML)
- Schwengers O et al 2021 Microb Genom 7 (Bakta)
- Vernikos GS et al 2015 Microb Genom 1:e000034 (pangenome openness review)
- Touchon M et al 2016 PLoS Genet 12:e1006413 (E. coli pangenome)
- Otto TD et al 2018 BMC Genomics 19:319 (Plasmodium pangenome)
- Brown CL et al 2022 mSystems 7:e0099122 (mobileOG-db)
- Mikheenko A et al 2018 Bioinformatics 34:i142 (QUAST-LG; long-read assembly evaluation -- earlier "pangenome review" attribution was incorrect; QUAST-LG benchmarks large-scale assembly QC).
- Chin C-S et al 2023 Nat Methods 20:1290 (PGR-TK; MAP graphs + principal bundle decomposition for repeat-rich / clinical genes; MHC, DAZ1-4, OPN1LW/OPN1MW examples)
- cschin/pgr-tk GitHub (repo; archived April 2026, transitioning to PANGEA developed by DGI / Diploid Genomics; consult upstream README for PANGEA pointer)

## Related Skills

- comparative-genomics/whole-genome-alignment - Minigraph-Cactus builds on Cactus; PGGB underlies eukaryotic pangenome
- comparative-genomics/ortholog-inference - Pangenome clusters are bacterial orthologs at species/genus level
- comparative-genomics/hgt-detection - Accessory genes often HGT-derived; mobile-element annotation cross-references
- comparative-genomics/gene-family-evolution - CAFE5 modeling on Panaroo presence/absence matrix
- genome-annotation/prokaryotic-annotation - Bakta annotation is pangenome input
- genome-annotation/repeat-annotation - Repeat masking before eukaryotic pangenome
- variant-calling/structural-variant-calling - Pangenome graph SV calling complements read-based callers
- variant-calling/joint-calling - vg + graph variant calling integrates with traditional VCF
- metagenomics/amr-detection - AMR genes often in bacterial accessory
- metagenomics/strain-tracking - Strain-specific accessory genes for tracking
- population-genetics/association-testing - Pan-GWAS for phenotype-gene-content association
