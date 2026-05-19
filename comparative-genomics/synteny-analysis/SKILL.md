---
name: bio-comparative-genomics-synteny-analysis
description: Detect syntenic blocks and structural rearrangements between genomes using MCScanX (Wang 2012), JCVI/MCScan (Tang 2008 Python), GENESPACE (Lovell 2022) for orthology-anchored riparian visualization, SyRI for structural variation, AnchorWave for sequence-level synteny, i-ADHoRe 3.0 for highly diverged species, SynNet for synteny networks, and ntSynt for multi-genome macrosynteny. Use when identifying collinear gene blocks across species, distinguishing macrosynteny from microsynteny, detecting inversions/translocations/duplications, anchoring orthology in WGD lineages, producing publication riparian plots, computing synteny block age via Ks (cross-references whole-genome-duplication), or running synteny-aware ortholog inference in polyploids.
tool_type: mixed
primary_tool: MCScanX
---

## Version Compatibility

Reference examples tested with: MCScanX 1.0+ (wyp1125/MCScanX commit 2020+), JCVI 1.4.21+ (Python port of MCScan), GENESPACE 1.4.0+ (Lovell 2022 eLife 78526), SyRI 1.7.1+ (Goel 2019 Genome Biol 20:277), plotsr 1.1.1+, AnchorWave 1.2.5+ (Song 2022 PNAS 119:e2113075119), i-ADHoRe 3.0.01+, SynNet (Zhao 2017 NAR 45:e108), ntSynt 1.0.4+ (2023), minimap2 2.28+, MUMmer 4.0.0+, OrthoFinder 3.0+, R 4.4+. plotsr requires pysam 0.22+ and seaborn 0.13+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `MCScanX -h`; `syri --version`; `python -m jcvi.compara.catalog ortholog --help`
- R: `packageVersion('GENESPACE')`; `?run_genespace`
- Python: `pip show jcvi`

If code throws `MCScanX: argument bad format`, `syri: input alignment file missing required columns`, or `GENESPACE: GFF parse error`, these tools have brittle input parsing: MCScanX requires 4-column `species_chr  gene  start  end` BED (non-standard), JCVI expects 4-column simple BED, GENESPACE requires GFF3 with `gene` feature type. Pre-process with `jcvi.formats.gff bed` or custom AWK.

# Synteny Analysis

**"Compare genome architecture between these species"** -> Detect conserved gene order (synteny) and infer rearrangement history. Synteny is NOT the same as collinearity: synteny is "genes on same chromosome", collinearity is "same order on same chromosome" (Fitch 1976 J Mol Evol 7:271 distinction; modern usage often conflates them). The choice of tool depends on whether the question is **gene-level co-linearity** (MCScanX, JCVI), **whole-genome structural rearrangements** (SyRI, AnchorWave), **multi-genome macrosynteny** (GENESPACE, ntSynt), or **synteny-aware orthology** (GENESPACE, ProteinOrtho-synteny). Repeat-masking quality is the dominant determinant of result reliability -- unmasked TEs produce ~100x more false anchor pairs than real syntenic anchors.

- CLI: `MCScanX` for collinear gene blocks via dynamic programming
- CLI: `python -m jcvi.compara.catalog ortholog A B` for JCVI/MCScan Python pipeline
- R: `run_genespace()` (Lovell 2022) for orthology-anchored riparian plots + pan-gene tracks
- CLI: `syri` for inversion / translocation / duplication detection
- CLI: `anchorwave proali` for sequence-level WGD-aware synteny

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| MCScanX (Wang 2012 NAR 40:e49) | Dynamic programming on BLAST hits with collinearity scoring | `.collinearity` blocks; tandem duplications collapsed | Most widely used; 14 downstream tools; well-benchmarked | Repeat-derived false BLAST hits; brittle 4-column input |
| MCScanX-h | Within-genome variant for WGD detection | Self-comparison blocks | Standard for paranome construction; WGD-aware | Same input fragility |
| JCVI / MCScan Python (Tang 2008 GR 18:1944) | Re-implementation of MCScanX with anchor stringency control | Publication-grade dotplots, karyotype, riparian | Best plotting; `--cscore` reciprocal-hit control | Slower than MCScanX C version on huge genomes |
| GENESPACE (Lovell 2022 eLife 78526) | OrthoFinder + MCScanX orthogroup-constrained synteny | Riparian plots; pan-gene tracks; CNV across genomes | Modern standard for plant comparative; integrates orthology | Slow at > 30 genomes; less rigorous for non-plant clades |
| i-ADHoRe 3.0 (Proost 2012 NAR 40:e11) | Iterative ordered-gene-list detection | Ghost gene families; deeply diverged synteny | Catches ancient synteny obscured by rearrangements | Computationally heavy at genome scale |
| AnchorWave (Song 2022 PNAS 119:e2113075119) | CDS/exon anchors -> wave-front alignment | Whole-genome alignment with WGD awareness; SVs | Best for plant WGD analyses with known ploidy | CDS-anchored only; intergenic resolution limited |
| SyRI (Goel 2019 GB 20:277) | Pairwise WGA -> syntenic-path identification -> SV calls | INV / TRANS / DUP / SYN / INS / DEL annotated | Comprehensive SV detection; works on chromosome-level assemblies | Requires chromosome-level assemblies; pairwise only |
| plotsr (Goel 2022 Bioinformatics 38:2922) | Multi-genome SyRI visualization | Stacked synteny + SV maps across N genomes | Best for visualizing 3-10 genome rearrangement histories | Inherits SyRI's pairwise input limitation |
| ntSynt (2023) | Minimizer-based alignment-free synteny | Multi-genome macrosynteny blocks | Alignment-free; handles > 15% divergence | Macrosynteny only; misses microsynteny |
| SynNet (Zhao 2017 NAR 45:e108) | Synteny block adjacency graphs | Synteny networks across many genomes | Phylogenetic network from synteny; detects deep ancestry | Less standard than block-based methods |
| Satsuma / progressive Cactus | Reference-free WGA | Whole-genome alignment (HAL format) | Underlies large-scale orthology; sequence-level synteny | See [[whole-genome-alignment]] |
| LASTZ chain/net (UCSC; Schwartz 2003 GR 13:103) | Pairwise WGA with chains and nets | Chains + nets in UCSC genome browser | Reference-anchored synteny; standard for UCSC tracks | See [[whole-genome-alignment]] |
| nucmer + dnadiff (MUMmer4; Marçais 2018 PLoS Comp Biol 14:e1005944) | MUM-anchored pairwise alignment | Whole-genome alignment with SV summary | Fast pairwise WGA for closely related | Sensitive only above ~70% identity |
| MashMap (Jain 2018 Bioinformatics 34:i748) | Approximate mapping for fragment-fragment synteny | Pairwise mappings with identity | Scales to thousands of genomes | Coarse (window-based); no SV inference |

Methodology evolves; GENESPACE has emerged as the de facto standard for plant comparative genomics (2022-2026). For non-plant clades, MCScanX or JCVI remain the workhorses. Whole-genome alignment for synteny is increasingly delegated to Cactus / Minigraph-Cactus and HAL toolkit (see [[whole-genome-alignment]]).

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Plant comparative genomics, 2-20 species | GENESPACE | Integrated orthology + synteny + visualization; modern standard |
| Animal / fungal, 2 species comparison | JCVI/MCScan | Best plots; flexible; widely cited |
| Animal / fungal, 10+ species | OrthoFinder + MCScanX + custom plot OR JCVI multi-genome | Avoid plant-specific GENESPACE assumptions |
| Bacterial / prokaryote synteny | progressiveMauve OR MCScanX-bacterial | Tools designed for compact genomes |
| Chromosome-level rearrangement inventory | SyRI + plotsr | Comprehensive SV (INV/TRANS/DUP) with publication viz |
| Polyploid genome analysis | AnchorWave proali mode with ploidy | WGD-aware synteny; subgenome-aware |
| Closely related strains (>= 95% identity) | nucmer + dnadiff | Faster than chromosome-level WGA; appropriate sensitivity |
| Distantly related species (< 70% identity) | i-ADHoRe 3.0 or LASTZ chains/nets | Sequence-level approaches lose power; ordered-list methods retain it |
| Ancient WGD detection | See [[whole-genome-duplication]] | wgd v2, KsRates; Ks-based dating outside synteny scope |
| WGA at clade-level (10+ vertebrates) | Cactus / Minigraph-Cactus -> halSynteny | See [[whole-genome-alignment]]; multi-genome reference-free |
| Pan-genome of bacterial strains | See [[pangenome-analysis]] | Panaroo / PPanGGOLiN / PEPPAN; different problem |
| Synteny-aware ortholog disambiguation | ProteinOrtho `-synteny` or GENESPACE | Tandem duplicates collapsed; co-orthologs assigned to syntenic position |
| Reference-guided assembly bias check | DO NOT use synteny across reference-guided scaffolds | Reference-guided scaffolds propagate scaffold assumptions; circular reasoning |
| Microsynteny network analysis | SynNet | Network from many genomes; detects ancient micro-conserved clusters |
| Synteny across > 15% divergent genomes | ntSynt or i-ADHoRe | Alignment-free or ordered-list approaches |
| Identifying syntelogs (syntenic orthologs only) | GENESPACE `riparian()` output; or MCScanX filtered by chr-pair | Syntelog tracks only same-chromosome-pair orthologs |

## Per-Tool Failure Modes

### Repeat-derived false synteny

**Trigger:** Running BLAST input for MCScanX/JCVI without softmasking repeats; or weak softmasking.

**Mechanism:** Transposable elements occupy 30-80% of plant and animal genomes; unmasked TEs produce millions of paralogous BLAST hits forming "synteny" blocks that are TE-driven, not real ancestry. The dynamic programming sees long apparent collinear blocks of TE-encoded proteins.

**Symptom:** MCScanX/JCVI output dominated by short (5-10 anchor) blocks; blocks cluster in TE-rich regions (heterochromatin, centromeric); BLAST hit count > 10x what is expected for the divergence.

**Fix:** Softmask genomes with RepeatModeler2 (de novo TE library) + RepeatMasker before BLAST; report softmasking statistics. For MCScanX, set `-e 1e-10` (stricter) for close species, `-e 1e-5` for distant. JCVI `--cscore 0.95` enforces near-reciprocal-best-hit, dramatically reducing TE-driven blocks. Always check that block size distribution drops at 5-10 anchors (TE-driven) vs proper > 20 anchor blocks (real synteny).

### Reference-guided assembly creating circular synteny

**Trigger:** Comparing a de novo assembly with a reference-guided scaffold of the same / related species.

**Mechanism:** Reference-guided scaffolding orders contigs based on synteny to the reference. Comparing the resulting scaffold to that same reference produces "synteny" that is actually the imposed reference order, not biological gene order. The comparison is circular.

**Symptom:** Suspicious uniformity in synteny across most chromosomes; visualization shows reference and query as perfectly parallel; SyRI reports few SVs where divergence-time-appropriate counts would predict many.

**Fix:** Verify assembly method via metadata. For reference-guided assemblies, treat synteny inferences with the reference as suspect; use a sister species reference or pure de novo assembly. Hi-C-scaffolded assemblies are not "reference-guided" in this sense and are safe.

### Fragmented assemblies underestimating macrosynteny

**Trigger:** Running synteny on draft assemblies with N50 < 1 Mb.

**Mechanism:** Macrosynteny detection requires that syntenic blocks be on single contigs; if contigs are short, blocks are artificially truncated. JCVI reports show artifactual "synteny loss" between species when one assembly is fragmented.

**Symptom:** Synteny "loss" correlates with assembly N50; per-chromosome dotplots show dotted lines where solid diagonals expected; SyRI reports excessive INVs / TRANS where the divergence is actually low.

**Fix:** Require minimum N50 of 1 Mb for cross-species synteny; chromosome-level assemblies preferred for SyRI. For fragmented assemblies, restrict analysis to long-contig regions or use SyRI's contig-level mode with caveats.

### Tandem duplicate inflation

**Trigger:** Including tandemly duplicated genes in synteny anchor input.

**Mechanism:** Tandem duplicates produce many same-region BLAST hits between species; the dynamic programming includes them as multiple anchors per gene position, artificially inflating block size and confidence.

**Symptom:** "High-confidence" synteny blocks in regions known for tandem expansions (NLR clusters in plants, OR genes in mammals); blocks anchored largely on adjacent tandem duplicates of the same gene.

**Fix:** MCScanX automatically collapses tandem duplicates (within 5 genes by default; `-w` parameter). Verify in the output the tandem-collapse step ran. JCVI `--tandem_Nmax 10` is more aggressive. For NLR / OR genes specifically, manual tandem filtering before BLAST is preferred.

### Cross-genome chromosome name mismatch

**Trigger:** GFF and FASTA from different sources with different chromosome naming (`chr1` vs `Chr1` vs `1` vs `NC_001234.5`).

**Mechanism:** MCScanX requires chromosomes to be named consistently; mismatch causes empty output or partial blocks.

**Symptom:** MCScanX completes but `.collinearity` file is small or empty; warnings about unknown chromosome IDs.

**Fix:** Normalize chromosome names with `awk` or `samtools faidx --regions` to unify naming convention; verify GFF + FASTA share the same names. For NCBI assemblies, use `datasets` to download with consistent labels.

### SyRI complaining about inversions in absence of true inversions

**Trigger:** Running SyRI on alignments with many small-scale rearrangements; or when one genome has been re-ordered (different convention for which strand is "+").

**Mechanism:** SyRI's syntenic-path algorithm penalizes any rearrangement; high-divergence pairs accumulate small-scale rearrangements that are not "real" inversions but inheritance differences.

**Symptom:** SyRI reports thousands of small "inversions" (< 1 kb); blocks dominate the SV count; total INV length unrealistic.

**Fix:** Filter SyRI output by size: real biologically-relevant INVs are typically > 5 kb. Tighten minimap2 sensitivity (`-x asm5` for close species; `-x asm10` or `-x asm20` for more divergent). Reverse-complement one assembly's chromosomes if the strand convention differs.

### GENESPACE OrthoFinder version mismatch

**Trigger:** Running GENESPACE with an outdated bundled OrthoFinder.

**Mechanism:** GENESPACE bundles a specific OrthoFinder version; if user's installed OrthoFinder is newer (v3 vs v2.5), the HOG output layout differs and GENESPACE parsing fails.

**Symptom:** GENESPACE reports "Phylogenetic_Hierarchical_Orthogroups not found" or similar; orthogroup table is empty.

**Fix:** Use the OrthoFinder version GENESPACE expects (currently OrthoFinder 2.5.4 for GENESPACE 1.4.x). Pin via `conda env`. Future GENESPACE versions will adapt to OrthoFinder 3 layout (Lovell update expected).

### Microsynteny vs macrosynteny conflation

**Trigger:** Reporting "synteny" without specifying scale.

**Mechanism:** Macrosynteny (same chromosome) is detectable across hundreds of millions of years but decays as rearrangements accumulate; microsynteny (same order within local region) can be deeply conserved even when macrosynteny is lost. Conflating them produces misleading "synteny loss" claims.

**Symptom:** Reports of "synteny lost between species X and Y" when both share microsyntenic gene clusters; or vice versa.

**Fix:** Always specify scale. SynNet (Zhao 2017) explicitly separates the two. MCScanX `-s` parameter (minimum block anchors) at 5 = microsynteny; at 20+ = macrosynteny. Report block-size distribution.

### Polyploid / WGD-affected genome ambiguity

**Trigger:** Running MCScanX on a polyploid (e.g. allohexaploid wheat) without subgenome assignment.

**Mechanism:** WGD doubles all genes; synteny is between subgenomes within the polyploid AND across to outgroup. Without subgenome assignment, all paralogs and orthologs collapse into the same blocks.

**Symptom:** Many "1:many" or "many:many" synteny relationships; per-chromosome synteny counts indicate doubled or tripled blocks; Ks distribution multimodal.

**Fix:** Use AnchorWave proali with explicit ploidy specification. Alternatively, run synteny twice: within polyploid (subgenome-vs-subgenome) and across to outgroup. See [[whole-genome-duplication]] for subgenome assignment workflow.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Assembly N50 for synteny | >= 1 Mb | Below this, results tool-dependent and unreliable; chromosome-level strongly preferred |
| MCScanX -s (minimum anchors per block) | 5 default; 10 stringent; 3 sensitive | Standard configuration; 5 is balance |
| MCScanX -m (maximum gaps between anchors) | 25 default | Higher for distant species; lower for recent |
| MCScanX -k (match score per anchor) | 50 default | Higher rewards longer blocks |
| MCScanX -e (BLAST e-value) | 1e-5 default; 1e-10 for close species | Stricter for recent radiations |
| BLAST evalue threshold | 1e-5 to 1e-10 typical | Depends on divergence |
| JCVI --cscore for reciprocal best | 0.7 default; 0.95 near-RBH; 0.99 RBH-only | Higher = fewer false positives, fewer hits |
| Tandem duplicate window | 5-10 genes (MCScanX default 5; JCVI 10) | Wang 2012; species-specific tuning |
| SyRI INV minimum size for biological significance | >= 5 kb | Below this, alignment noise dominates |
| SyRI TRANS minimum size | >= 1 kb | Standard convention |
| GENESPACE minimum syntenic block | 5 orthogroups | Lovell 2022 default |
| Synteny block decay (macrosynteny half-life) | ~150 Myr in vertebrates | Naruse 2004; Murat 2010 |
| Microsynteny conservation | up to 1 Gyr for metabolic gene clusters | Slot & Rokas 2010 |
| minimap2 preset for synteny | -x asm5 for < 5% divergence; asm10 for ~10%; asm20 for ~20% | minimap2 docs |
| MUMmer nucmer maxmatch | --maxmatch for SyRI; --mum default | MUMmer4 manual |
| Ks for syntenic block age (cross-references WGD) | Ks 0.1-0.5 recent; 0.5-1.5 older; > 1.5 saturated | See [[whole-genome-duplication]] for Ks plot interpretation |
| Repeat masking minimum | >= 90% of known TE families masked (RepeatMasker `.tbl`) | Below this, expect spurious synteny |

## MCScanX Standard Pipeline

**Goal:** Detect collinear gene blocks between two genomes.

**Approach:** Prepare 4-column BED with species prefix -> all-vs-all BLASTP -> run MCScanX -> parse `.collinearity` -> classify blocks.

```bash
# 1. Prepare MCScanX-format BED (species_chr  gene  start  end)
python -m jcvi.formats.gff bed --type=gene --key=ID species_A.gff > A.gff.tmp
python -m jcvi.formats.gff bed --type=gene --key=ID species_B.gff > B.gff.tmp
awk 'BEGIN{OFS="\t"}{print "A"$1, $4, $2, $3}' A.gff.tmp > work/A.gff
awk 'BEGIN{OFS="\t"}{print "B"$1, $4, $2, $3}' B.gff.tmp > work/B.gff
cat work/A.gff work/B.gff > work/A_B.gff

# 2. All-vs-all DIAMOND BLAST (faster than BLAST+)
diamond makedb --in species_A.faa --db work/A.dmnd
diamond makedb --in species_B.faa --db work/B.dmnd
diamond blastp --db work/A.dmnd --query species_A.faa --threads 16 \
    --outfmt 6 --evalue 1e-10 --max-target-seqs 5 --out work/AA.blast
diamond blastp --db work/B.dmnd --query species_B.faa --threads 16 \
    --outfmt 6 --evalue 1e-10 --max-target-seqs 5 --out work/BB.blast
diamond blastp --db work/B.dmnd --query species_A.faa --threads 16 \
    --outfmt 6 --evalue 1e-10 --max-target-seqs 5 --out work/AB.blast
diamond blastp --db work/A.dmnd --query species_B.faa --threads 16 \
    --outfmt 6 --evalue 1e-10 --max-target-seqs 5 --out work/BA.blast
cat work/AA.blast work/BB.blast work/AB.blast work/BA.blast > work/A_B.blast

# 3. Run MCScanX
cd work && MCScanX -s 5 -m 25 -k 50 -e 1e-10 A_B
# Output: A_B.collinearity (block table), A_B.tandem (tandem clusters), A_B.html (browsing)
```

```python
'''Parse MCScanX .collinearity and classify synteny relationships.'''
import re
from collections import defaultdict


def parse_collinearity(path):
    '''Returns list of dicts: {block_id, n_anchors, e_value, score, gene_pairs: [(g1, g2), ...]}'''
    blocks = []
    current = None
    with open(path) as fh:
        for line in fh:
            if line.startswith('## Alignment'):
                if current:
                    blocks.append(current)
                m = re.match(r'## Alignment (\d+): score=([0-9.]+) e_value=([0-9.e\-]+) N=(\d+)', line)
                if m:
                    current = {
                        'block_id': int(m.group(1)),
                        'score': float(m.group(2)),
                        'e_value': float(m.group(3)),
                        'n_anchors': int(m.group(4)),
                        'gene_pairs': []
                    }
            elif current and ':' in line and '\t' in line:
                parts = line.strip().split()
                if len(parts) >= 3:
                    current['gene_pairs'].append((parts[1], parts[2]))
    if current:
        blocks.append(current)
    return blocks


def classify_chromosome_synteny(blocks, gene_to_chr):
    '''Classify syntenic chromosome relationships: 1-1, 1-many, many-many.'''
    a_partners = defaultdict(set)
    for blk in blocks:
        for g1, g2 in blk['gene_pairs']:
            c1, c2 = gene_to_chr.get(g1), gene_to_chr.get(g2)
            if c1 and c2:
                a_partners[c1].add(c2)
    result = {}
    for chr_a, partners in a_partners.items():
        n = len(partners)
        result[chr_a] = '1-1' if n == 1 else ('1-many' if n <= 3 else 'many-many')
    return result
```

## GENESPACE Plant-Focused Pipeline

**Goal:** Build pan-gene tracks across N plant genomes with orthology-anchored synteny visualization.

**Approach:** Format input: per-species GFF + protein FASTA -> initialize -> run -> riparian plot.

```r
library(GENESPACE)

dir.create('work_dir')
parsedPaths <- parse_annotations(
    rawGenomeRepo = 'raw_genomes/',
    genomeDirs = list.dirs('raw_genomes/', recursive = FALSE),
    headerEntryIndex = 1,
    gffString = 'gff',
    faString = 'protein.faa',
    genespaceWd = 'work_dir/'
)

gpar <- init_genespace(
    wd = 'work_dir/',
    nCores = 8,
    rawOrthofinderDir = NULL,  # GENESPACE bundles OrthoFinder 2.5.x
    onewayBlast = TRUE
)

out <- run_genespace(gsParam = gpar, overwrite = TRUE)

# Riparian plot for 5 genomes vs reference
ripDat <- plot_riparian(
    gsParam = out,
    refGenome = 'Reference_species',
    useOrder = FALSE,
    backgroundColor = 'white'
)
```

GENESPACE outputs include `results/syntenicHits.txt` (anchor pairs), `results/pangenes.txt` (pan-gene presence/absence matrix), and `results/riparian.pdf` (riparian visualization). Pan-gene tracks across species are particularly useful for identifying lineage-specific genes (cf. pangenome-analysis for bacterial scope).

## SyRI for Structural Rearrangement Inventory

**Goal:** Identify INVs, TRANSs, DUPs, INS, DEL between two chromosome-level genome assemblies.

**Approach:** Align with minimap2 -x asm5 -> filter to chromosome-level mappings -> run SyRI -> visualize with plotsr.

```bash
# 1. Whole-genome alignment with minimap2
minimap2 -ax asm5 --eqx -t 16 reference.fa query.fa | \
    samtools sort -@ 8 -O bam - > work/aln.bam
samtools index work/aln.bam

# 2. SyRI on the alignment
syri -c work/aln.bam -r reference.fa -q query.fa -F B --prefix work/syri_out -k

# Outputs:
#   syri_out.syri.out      structural variant table (SYN/INV/TRANS/DUP/INS/DEL)
#   syri_out.syri.vcf      VCF-format variants
#   syri_out.syri.summary  block-level summary

# 3. Visualize with plotsr
echo "ref_species  query_species" > work/genome_pairs.txt
plotsr --sr work/syri_out.syri.out --genomes work/genome_pairs.txt \
    --output work/syri_plot.png --markers gene_markers.bed
```

For multi-genome rearrangement: run SyRI pairwise for each adjacent species in a phylogeny, then plotsr stacks the panel.

## JCVI / MCScan Python for Publication Plots

**Goal:** Produce publication-grade dotplots, karyotype, and synteny visualization.

**Approach:** Use JCVI's catalog ortholog workflow to derive anchors -> generate dotplot or karyotype directly.

```bash
# 1. Detect orthologs and synteny
python -m jcvi.formats.gff bed --type=mRNA --key=Name species_A.gff > A.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name species_B.gff > B.bed
python -m jcvi.formats.fasta format species_A.faa A.fasta
python -m jcvi.formats.fasta format species_B.faa B.fasta

python -m jcvi.compara.catalog ortholog --no_strip_names A B
# Produces A.B.lifted.anchors and other files

# 2. Dotplot
python -m jcvi.graphics.dotplot A.B.anchors --dpi 300 --output dotplot.png

# 3. Karyotype (chromosome-level synteny)
# Layout file: each row = species  chr  start  end  reverse
cat > layout.csv << 'EOF'
A.chr1 species_A 1 50000000 1
B.chr2 species_B 1 60000000 1
EOF
python -m jcvi.graphics.karyotype seqids.txt layout.csv

# 4. Microsynteny block plot
python -m jcvi.graphics.synteny blocks.layout
```

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MCScanX many blocks, GENESPACE few | GENESPACE uses orthogroup-constrained anchors (more conservative) | GENESPACE excludes non-orthogroup hits; MCScanX includes any high BLAST hit; trust GENESPACE for orthology-anchored synteny |
| GENESPACE 1:1 syntelog, MCScanX 1:many | Tandem duplicates in MCScanX collapsed at different granularity | Trust GENESPACE for clear syntelogs; MCScanX 1:many often reflects tandem expansions |
| SyRI reports many INVs, plotsr shows mostly SYN | Small INVs below biological-relevance threshold | Filter SyRI INVs < 5 kb; SyRI was correct, just noisy |
| MCScanX block on chromosomes that JCVI doesn't pair | JCVI `--cscore 0.7` more restrictive | Lower JCVI cscore to 0.5 OR raise MCScanX -e to 1e-10 |
| AnchorWave finds collinearity in WGD region, MCScanX doesn't | AnchorWave WGD-aware; MCScanX treats all anchors equally | AnchorWave correct for WGD lineages |
| ntSynt macrosynteny across genus, MCScanX microsynteny only | Different scales | Both correct at their scale; report both with explicit scale labels |
| SyRI calls "translocation" but it's a known chromosome split | Centromere break / Robertsonian translocation | Manual cytogenetic context; SyRI's mechanical detection is correct but may not be biologically novel |
| Recent species pair has > 1000 SVs in SyRI | Either real chromosome instability OR assembly errors | Check assembly QC (telomere completeness; BUSCO); cross-check with PacBio long-read SVs |
| Synteny on reference-guided scaffold | Circular reasoning | Discard synteny calls between scaffold and its reference |

**Operational rule for publication:** Synteny block annotation requires (1) BUSCO/Compleasm > 90% complete on both genomes; (2) softmasked repeats verified; (3) at least one cross-validation tool (MCScanX vs JCVI or MCScanX vs GENESPACE); (4) microsynteny vs macrosynteny scale explicit; (5) for SV calls, plotsr / SyRI followed by manual review of large rearrangements with read-coverage support.

## Cohort Gotchas

- **Polyploid plants:** WGD events confound 1:1 synteny; use AnchorWave or subgenome-assigned analysis
- **Salmonids / fish 4R:** Ts3R + Ss4R WGD events; ohnologs (WGD paralogs) appear as syntenic but not orthologous; cross-reference [[whole-genome-duplication]]
- **Mammalian X chromosome:** Massive recombination suppression; "synteny" sometimes inverted in opposite sex; verify strand convention
- **Centromeric regions:** generally unalignable; synteny tools may produce false breaks near centromeres
- **Telomeres:** assembly often incomplete; synteny calls near telomere boundaries unreliable
- **Sex chromosomes:** rapidly evolving; lower synteny signal than autosomes
- **Plant B chromosomes:** supernumerary; exclude from synteny analysis

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Assembly completeness?" | BUSCO/Compleasm > 90%; N50 reported per assembly; chromosome-level (or scaffold N50 > 1 Mb) |
| "Repeat masking?" | RepeatModeler2 de novo TE library + RepeatMasker; > 90% of known TE families masked; softmasked, not hardmasked |
| "Reference-guided scaffolding?" | Verified de novo OR Hi-C scaffolded; no reference-guided steps |
| "Tandem duplicates?" | Collapsed by MCScanX automatic detection (window 5); JCVI tandem_Nmax 10 |
| "Cross-validation?" | MCScanX results validated against JCVI (or GENESPACE) at consistent stringency |
| "Microsynteny vs macrosynteny?" | Reported both scales; minimum block size 5 anchors for microsynteny, 20+ for macrosynteny |
| "WGD-aware?" | AnchorWave proali with ploidy specification (for polyploids); or restricted to non-WGD lineages |
| "SV calling sensitivity?" | SyRI INVs > 5 kb; cross-validated with PacBio long-read SV calls where available |
| "Synteny block age (Ks)?" | Ks plotted per block; saturation threshold (Ks > 2) applied; see [[whole-genome-duplication]] |
| "Bacterial pangenome reference?" | Not applicable here; see [[pangenome-analysis]] for prokaryotes |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| MCScanX produces empty .collinearity | 4-column BED format wrong (extra columns, wrong delimiter, missing prefix) | Check exact 4-column format: `species_prefix_chr  gene_id  start  end`; tabs only |
| MCScanX dumps "no matches" message | BLAST e-value too strict; or gene IDs don't match BED | Verify gene IDs in BED match BLAST input; relax e-value |
| GENESPACE "OrthoFinder version mismatch" | Bundled OrthoFinder version conflict | Pin OrthoFinder 2.5.x for GENESPACE 1.4.x; or update GENESPACE |
| JCVI dotplot empty | `--cscore` too strict; few anchors survive | Lower cscore to 0.5; or check BLAST hit count |
| SyRI "no chromosomes match" | Reference and query use different chromosome IDs | Normalize names: `samtools faidx --regions` to extract by name |
| SyRI INVs appear in all assemblies | Strand convention issue OR repeat-driven | Reverse-complement one assembly's chromosomes if convention differs; mask repeats |
| AnchorWave proali times out | Genome > 1 Gb with high TE density | Pre-mask repeats; consider GENESPACE / MCScanX alternative |
| GENESPACE riparian plot too crowded | > 10 genomes | Subset to representative species or use multi-page riparian |
| Per-chromosome BLAST hit count > 100,000 | Unmasked repeats | Softmask before BLAST; verify by counting TE-encoded protein hits |
| Synteny blocks on Y chromosome appear absent | Y not assembled fully; or sex-specific | Document assembly limitation; exclude Y from synteny analysis |

## Tool Installation Notes

```bash
# MCScanX
git clone https://github.com/wyp1125/MCScanX && cd MCScanX && make
# JCVI / MCScan Python
pip install jcvi
# GENESPACE
remotes::install_github('jtlovell/GENESPACE')
# SyRI + plotsr
conda install -c bioconda syri plotsr
# AnchorWave
conda install -c bioconda anchorwave
# i-ADHoRe
git clone https://github.com/VIB-PSB/i-ADHoRe && cd i-ADHoRe && mkdir build && cd build && cmake .. && make
# ntSynt
pip install ntsynt
# minimap2 + MUMmer4
conda install -c bioconda minimap2 mummer4

# Repeat masking
conda install -c bioconda repeatmodeler repeatmasker
```

For GENESPACE, OrthoFinder 2.5.x must be pinned; install via `conda install -c bioconda orthofinder=2.5.5`. JCVI is the most generally useful Python toolkit; install per-project alongside specific synteny tools.

## References

- Fitch WM 1976 J Mol Evol 7:271 (synteny vs collinearity)
- Wang Y et al 2012 NAR 40:e49 (MCScanX)
- Tang H et al 2008 GR 18:1944 (synteny / MCScan Python)
- Lovell JT et al 2022 eLife 11:78526 (GENESPACE)
- Proost S et al 2012 NAR 40:e11 (i-ADHoRe 3.0)
- Song B et al 2022 PNAS 119:e2113075119 (AnchorWave)
- Goel M et al 2019 Genome Biol 20:277 (SyRI)
- Goel M et al 2022 Bioinformatics 38:2922 (plotsr)
- Zhao T et al 2017 NAR 45:e108 (SynNet synteny network)
- Marçais G et al 2018 PLoS Comp Biol 14:e1005944 (MUMmer4)
- Schwartz S et al 2003 GR 13:103 (LASTZ chains and nets)
- Jain C et al 2018 Bioinformatics 34:i748 (MashMap)
- Li H 2018 Bioinformatics 34:3094 (minimap2)
- Slot JC & Rokas A 2010 GBE 2:362 (microsynteny conservation in metabolic clusters)
- Naruse K et al 2004 GR 14:820 (synteny block decay)
- Murat F et al 2010 GR 20:1545 (vertebrate macrosynteny)
- Holland PWH et al 1994 Development Suppl:125 (2R hypothesis)
- Vanneste K et al 2013 GR 23:1304 (Ks saturation)
- Freeling M 2007 PNAS 104:8723 (gene balance)
- Force A et al 1999 Genetics 151:1531 (subfunctionalization)
- Zhao T & Schranz ME 2017 NAR 45:e108 (synteny network for phylogeny)
- Smith MR & Hahn MW 2021 PNAS 118:e2103725118 (gene-tree-aware synteny)

## Related Skills

- comparative-genomics/whole-genome-duplication - Ks-based WGD detection, paranome construction, KsRates
- comparative-genomics/whole-genome-alignment - Cactus / Minigraph-Cactus / LASTZ chains-and-nets for sequence-level synteny
- comparative-genomics/ortholog-inference - GENESPACE depends on OrthoFinder; synteny verifies orthology in WGD lineages
- comparative-genomics/pangenome-analysis - Bacterial pangenome (different problem; for prokaryotes)
- comparative-genomics/comparative-annotation-projection - TOGA + CESAR project genes through WGA-derived synteny
- comparative-genomics/positive-selection - Selection on syntenic gene pairs (1:1 syntelogs)
- phylogenetics/modern-tree-inference - Phylogenetic context for synteny block dating
- genome-assembly/assembly-qc - BUSCO / Compleasm assembly completeness check before synteny
- alignment/structural-alignment - Sequence-level alignment underlying SyRI
- variant-calling/structural-variant-calling - SV detection from short reads (orthogonal to SyRI WGA-based)
