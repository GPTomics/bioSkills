---
name: bio-comparative-genomics-synteny-analysis
description: Analyze genome collinearity and syntenic blocks using MCScanX, SyRI, and JCVI for comparative genomics. Detect conserved gene order, chromosomal rearrangements, and whole-genome duplications. Use when comparing genome structure between species or identifying conserved genomic regions.
tool_type: mixed
primary_tool: MCScanX
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, JCVI 1.3+, PAML 4.10+, matplotlib 3.8+, minimap2 2.26+, numpy 1.26+, pandas 2.2+, scipy 1.12+, SyRI 1.6+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Synteny Analysis

**"Compare genome structure between my species"** → Detect conserved gene order (syntenic blocks), chromosomal rearrangements, and whole-genome duplications by aligning genomic collinearity.
- CLI: `MCScanX` for collinear block detection from BLAST results
- Python: `jcvi.compara.synteny` for synteny visualization (dot plots, macro/micro)

## Assembly Quality Requirements

Synteny analysis is highly sensitive to assembly fragmentation. Key thresholds:

- **Minimum N50 of 1 Mb** for robust synteny detection; below this, results are tool-dependent and unreliable
- Fragmented assemblies can **underestimate synteny by up to 40%** (different tools disagree extensively)
- Gene-dense genomes tolerate more fragmentation than gene-sparse ones
- Chromosome-level assemblies are strongly preferred for detecting inversions and translocations
- **Reference-guided scaffolding creates false synteny** -- scaffolding against a related species propagates gene order assumptions, producing circular reasoning in synteny comparisons

Always verify assembly contiguity (N50, BUSCO completeness) before interpreting synteny results. Report assembly quality alongside synteny findings.

## Macrosynteny vs Microsynteny

- **Macrosynteny**: Chromosome-level conservation of gene content (same chromosome, not necessarily same order). Detectable across hundreds of millions of years but decays via accumulation of rearrangements.
- **Microsynteny**: Local conservation of gene order at sub-chromosomal level (few to dozens of genes). Can be deeply conserved, especially for metabolic gene clusters. Used in synteny network analysis for phylogenetic inference.

MCScanX and JCVI detect both scales; SyRI focuses on structural rearrangements at the macrosynteny level.

## Tool Selection

| Tool | Best For | Key Characteristic |
|---|---|---|
| MCScanX | General synteny, WGD detection, downstream analyses (14 tools) | Dynamic programming on BLAST hits; most widely used |
| JCVI/MCScan (Python) | Visualization, publication figures, multi-genome | Superior plotting; `--cscore 0.99` for reciprocal best hits |
| i-ADHoRe 3.0 | Ancient WGD, highly diverged species | Ordered gene lists (no sequence needed); iterative detection |
| AnchorWave | Plant genomes with WGD, polyploidy | CDS/exon anchors; WGD-aware alignment for known ploidy levels |
| SyRI | Structural rearrangements between assemblies | Inversions, translocations, duplications from whole-genome alignment |
| ntSynt | Multi-genome macrosynteny at scale | Alignment-free, minimizer graphs; handles >15% divergence |

## Repeat Masking

Softmask genomes before synteny analysis. Unmasked repeats produce millions of spurious BLAST alignments from transposable elements, drowning real syntenic signal. Use RepeatMasker with a species-specific library (or build one with RepeatModeler2 for non-model organisms). Always softmask (lowercase), never hardmask (N's), to preserve information for downstream tools.

## MCScanX Workflow

**Goal:** Detect conserved gene order (syntenic blocks) between two genomes.

**Approach:** Prepare GFF and all-vs-all BLASTP input files, run MCScanX to identify collinear gene blocks, parse the collinearity output, and classify syntenic relationships by coverage ratios.

```python
'''Synteny analysis with MCScanX and visualization'''

import subprocess
import pandas as pd
from collections import defaultdict

def prepare_mcscanx_input(gff_file, fasta_file, species_prefix):
    '''Prepare input files for MCScanX

    MCScanX requires:
    1. .gff file: gene positions (sp  gene  chr  start  end)
    2. .blast file: all-vs-all BLASTP results
    '''
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'gene':
                chrom = parts[0]
                start, end = int(parts[3]), int(parts[4])
                gene_id = parts[8].split('ID=')[1].split(';')[0]
                genes.append(f'{species_prefix}\t{gene_id}\t{chrom}\t{start}\t{end}')

    with open(f'{species_prefix}.gff', 'w') as f:
        f.write('\n'.join(genes))

    return f'{species_prefix}.gff'


def run_mcscanx(gff1, gff2, blast_file, output_prefix):
    '''Run MCScanX for synteny detection

    Key parameters:
    -k 50: Match score per collinear gene pair
    -g -1: Gap penalty per intervening gene
    -s 5: Minimum block size (5 genes default; use 10+ for stringent analysis)
    -m 25: Maximum gaps between anchor pairs (25 default)
    -e 1e-5: BLAST E-value threshold (use 1e-10 for close species)

    Default parameters: ~0.59 sensitivity / 0.81 precision on benchmarks.
    Tandem duplications are automatically collapsed before detection.
    '''
    # Combine GFF files
    subprocess.run(f'cat {gff1} {gff2} > {output_prefix}.gff', shell=True)

    # Copy BLAST file
    subprocess.run(f'cp {blast_file} {output_prefix}.blast', shell=True)

    # Run MCScanX
    # -s 5: Min genes per block (raise to 10 for stringent, lower to 3 for sensitive)
    # -m 25: Max gaps (raise for degraded synteny, lower for recent comparisons)
    cmd = f'MCScanX -s 5 -m 25 {output_prefix}'
    subprocess.run(cmd, shell=True)

    return f'{output_prefix}.collinearity'


def parse_collinearity(collinearity_file):
    '''Parse MCScanX collinearity output

    Output format:
    ## Alignment X: score=N e_value=X N genes
    X-Y: gene1  gene2
    '''
    blocks = []
    current_block = None

    with open(collinearity_file) as f:
        for line in f:
            if line.startswith('## Alignment'):
                if current_block:
                    blocks.append(current_block)
                parts = line.strip().split()
                score = int(parts[3].split('=')[1])
                e_value = float(parts[4].split('=')[1])
                n_genes = int(parts[5])
                current_block = {
                    'score': score,
                    'e_value': e_value,
                    'n_genes': n_genes,
                    'gene_pairs': []
                }
            elif current_block and '-' in line and ':' in line:
                parts = line.strip().split()
                if len(parts) >= 3:
                    gene1, gene2 = parts[1], parts[2]
                    current_block['gene_pairs'].append((gene1, gene2))

    if current_block:
        blocks.append(current_block)

    return blocks


def classify_synteny_type(blocks, species1_chroms, species2_chroms):
    '''Classify syntenic relationships

    Types:
    - 1:1: Direct orthology (conserved)
    - 1:many: Lineage-specific duplication
    - many:many: Ancient WGD or complex rearrangement
    '''
    sp1_coverage = defaultdict(list)
    sp2_coverage = defaultdict(list)

    for block in blocks:
        for gene1, gene2 in block['gene_pairs']:
            chr1 = species1_chroms.get(gene1)
            chr2 = species2_chroms.get(gene2)
            if chr1 and chr2:
                sp1_coverage[chr1].append(chr2)
                sp2_coverage[chr2].append(chr1)

    classifications = []
    for chr1, partners in sp1_coverage.items():
        unique_partners = len(set(partners))
        if unique_partners == 1:
            classifications.append(('1:1', chr1, partners[0]))
        else:
            classifications.append(('1:many', chr1, set(partners)))

    return classifications
```

## SyRI for Structural Variants

**Goal:** Identify structural rearrangements (inversions, translocations, duplications) between two genome assemblies.

**Approach:** Align genomes with minimap2, run SyRI on the alignment to detect syntenic regions and structural variants, and parse the output into a variant table.

```python
def run_syri(ref_genome, query_genome, alignment_file, output_prefix):
    '''Run SyRI for structural rearrangement identification

    SyRI detects:
    - Syntenic regions (SYN)
    - Inversions (INV)
    - Translocations (TRANS)
    - Duplications (DUP)
    - Insertions/Deletions (INS/DEL)

    Requires whole-genome alignment (minimap2 or MUMmer)
    '''
    # Align genomes with minimap2
    align_cmd = f'minimap2 -ax asm5 {ref_genome} {query_genome} > {output_prefix}.sam'
    subprocess.run(align_cmd, shell=True)

    # Run SyRI
    syri_cmd = f'syri -c {output_prefix}.sam -r {ref_genome} -q {query_genome} -F S --prefix {output_prefix}'
    subprocess.run(syri_cmd, shell=True)

    return f'{output_prefix}syri.out'


def parse_syri_output(syri_file):
    '''Parse SyRI structural variant output'''
    variants = []

    with open(syri_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                var_type = parts[10]
                ref_chr, ref_start, ref_end = parts[0], int(parts[1]), int(parts[2])
                qry_chr, qry_start, qry_end = parts[5], int(parts[6]), int(parts[7])
                variants.append({
                    'type': var_type,
                    'ref_chr': ref_chr,
                    'ref_start': ref_start,
                    'ref_end': ref_end,
                    'qry_chr': qry_chr,
                    'qry_start': qry_start,
                    'qry_end': qry_end,
                    'size': ref_end - ref_start
                })

    return pd.DataFrame(variants)
```

## JCVI Visualization

```python
def create_synteny_plot(blocks_file, layout_file, output_file):
    '''Create synteny dot plot with JCVI

    JCVI (python-jcvi) provides publication-ready figures.

    Key JCVI parameters for pairwise synteny:
    --cscore 0.70: Default hit stringency (ratio of best to 2nd-best hit)
    --cscore 0.99: Reciprocal best hits only (strictest, fewest false positives)
    --minspan 20: Minimum syntenic block size in genes
    --tandem_Nmax 10: Filter tandem duplicates within this gene distance
    '''
    from jcvi.graphics.dotplot import dotplot
    from jcvi.graphics.karyotype import karyotype

    # Dotplot shows collinear blocks as diagonal lines
    # Good for detecting WGD (parallel diagonals)
    dotplot_cmd = f'python -m jcvi.graphics.dotplot {blocks_file}'
    subprocess.run(dotplot_cmd, shell=True)

    # Karyotype view for chromosome-level comparison
    karyotype_cmd = f'python -m jcvi.graphics.karyotype {layout_file}'
    subprocess.run(karyotype_cmd, shell=True)


def detect_wgd(blocks, min_parallel_blocks=3):
    '''Detect whole-genome duplication signatures

    WGD indicators:
    - Multiple parallel syntenic blocks in self-comparison
    - 2:1 or 4:1 chromosome ratios between genomes
    - Ks distribution peak at consistent value across syntenic pairs
    - MCScanX "duplication depth" of 1 = one WGD round

    Pitfalls:
    - Very recent WGDs (low Ks) can be masked by continuous small-scale duplications
    - Ks peaks from different WGDs can merge if timing is similar
    - Ks > 2 is saturated and unreliable for dating ancient WGDs
    - Use wgd v2 tool for mixture model fitting and phylogenetic dating

    min_parallel_blocks: Minimum parallel blocks to call WGD (3 minimum)
    '''
    chrom_ratios = defaultdict(list)

    for block in blocks:
        if block['n_genes'] >= 10:  # Focus on substantial blocks
            pairs = block['gene_pairs']
            # Extract chromosome info from gene names
            # Implementation depends on naming convention
            pass

    return chrom_ratios
```

## Ks Analysis for Dating

**Goal:** Date gene duplication and speciation events using synonymous substitution rates (Ks).

**Approach:** Calculate Ks for syntenic gene pairs using PAML yn00/codeml, plot the Ks distribution, and identify peaks corresponding to whole-genome duplication events.

```python
def calculate_ks_for_pairs(cds_file1, cds_file2, gene_pairs):
    '''Calculate synonymous substitution rate (Ks) for gene pairs

    Ks interpretation:
    - Ks < 0.1: Very recent divergence or gene conversion
    - Ks 0.1-0.5: Within-species duplication (recent WGD)
    - Ks 0.5-1.5: Between closely related species (older WGD)
    - Ks 1.5-2.0: Approaching saturation, interpret with caution
    - Ks > 2.0: Saturated -- synonymous sites have undergone multiple substitutions

    Ks peaks indicate WGD events. Use mixture models (e.g., wgd v2) to
    formally fit components rather than visually identifying peaks.
    Polyploids show multiple peaks (one per WGD event).
    '''
    from Bio import SeqIO
    from Bio.Seq import Seq

    # Load CDS sequences
    cds1 = SeqIO.to_dict(SeqIO.parse(cds_file1, 'fasta'))
    cds2 = SeqIO.to_dict(SeqIO.parse(cds_file2, 'fasta'))

    ks_values = []
    for gene1, gene2 in gene_pairs:
        if gene1 in cds1 and gene2 in cds2:
            # Run yn00 or codeml for Ks calculation
            # Simplified - actual implementation uses PAML
            pass

    return ks_values


def plot_ks_distribution(ks_values, output_file):
    '''Plot Ks distribution to identify WGD peaks'''
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import gaussian_kde

    # Filter saturated values
    # Ks > 2 is typically saturated and unreliable
    ks_filtered = [k for k in ks_values if 0 < k < 2]

    fig, ax = plt.subplots(figsize=(10, 6))

    # Histogram
    ax.hist(ks_filtered, bins=50, density=True, alpha=0.7, label='Ks distribution')

    # KDE for peak detection
    if len(ks_filtered) > 10:
        kde = gaussian_kde(ks_filtered)
        x = np.linspace(0, 2, 200)
        ax.plot(x, kde(x), 'r-', linewidth=2, label='KDE')

    ax.set_xlabel('Ks (synonymous substitution rate)')
    ax.set_ylabel('Density')
    ax.legend()
    plt.savefig(output_file)

    return fig
```

## Polyploidy Considerations

For polyploid genomes (common in plants, ~35% are recent polyploids):
- **Autopolyploids**: Subgenome assignment is extremely challenging due to high homeolog similarity; risk of chimeric assembly
- **Allopolyploids**: Easier to resolve subgenomes (greater parental divergence); assign subgenomes before comparative analysis
- AnchorWave handles WGD-aware alignment when ploidy level is specified (`proali` mode with `-R` for ploidy)
- OrthoFinder and similar tools cannot natively handle polyploid genomes -- each subgenome appears as "extra" genes; assign subgenomes first or use multi-labeled tree methods
- Ks plots in polyploids show multiple peaks (one per WGD event); peaks from different events can merge if timing overlaps

## Related Skills

- comparative-genomics/positive-selection - dN/dS analysis on syntenic gene pairs
- comparative-genomics/ortholog-inference - Identify orthologs for synteny context
- phylogenetics/modern-tree-inference - Phylogenetic context for synteny dating
- alignment/pairwise-alignment - Sequence alignment for Ks calculation
- genome-annotation/annotation-transfer - Transfer annotations using syntenic context
