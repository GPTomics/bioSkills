---
name: bio-comparative-genomics-ortholog-inference
description: Infer orthologous gene groups across species using OrthoFinder and ProteinOrtho. Identify orthologs, paralogs, and co-orthologs for comparative genomics and functional annotation transfer. Use when identifying gene orthologs across species or building orthogroups for evolutionary analysis.
tool_type: cli
primary_tool: OrthoFinder
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, BUSCO 5.5+, NCBI BLAST+ 2.15+, OrthoFinder 2.5+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Ortholog Inference

**"Find orthologs across my species"** → Identify orthologous gene groups, paralogs, and co-orthologs across multiple species using sequence similarity clustering and gene tree reconciliation.
- CLI: `orthofinder -f proteomes/` for all-vs-all orthogroup inference

## Method Selection

| Method | Approach | Best For | Tradeoff |
|---|---|---|---|
| OrthoFinder | Tree-based (gene tree reconciliation with species tree) | Accuracy, evolutionary analysis, gene duplication events | Slower, needs sufficient species |
| ProteinOrtho | Graph-based (reciprocal best hits + connectivity) | Speed, many genomes, quick surveys | Less accurate for complex gene families |
| OMA/FastOMA | Graph-based (strict pairwise, hierarchical groups) | Precision-critical applications, large-scale (1000+ genomes) | Lowest recall (misses distant orthologs) |
| SonicParanoid2 | Graph-based (ML predictor + protein language model) | Fast + accurate graph-based | Newer, less community testing |

**Tree-based methods** (OrthoFinder) build gene trees and reconcile with the species tree to distinguish speciation (orthology) from duplication (paralogy). More accurate but computationally expensive.

**Graph-based methods** (ProteinOrtho, OMA, SonicParanoid) use sequence similarity with clustering. Faster but can confuse paralogs with orthologs when evolutionary rates vary.

Default recommendation: OrthoFinder for most analyses. ProteinOrtho for quick surveys or 50+ genomes. OMA/FastOMA when precision is paramount.

## Input Quality

Annotation quality directly affects orthology inference. Heterogeneous annotations across species spuriously inflate lineage-specific gene counts, creating false gene family expansions/contractions in downstream CAFE analysis.

- Use consistent annotation pipelines across species when possible
- Verify proteome completeness with BUSCO/Compleasm before running orthology
- Remove isoforms (keep longest per gene) to avoid inflating copy numbers
- Incomplete gene models produce truncated proteins that split true orthogroups

## Orthology Subtypes

- **One-to-one orthologs**: single gene in each species, ideal for phylogenomics
- **One-to-many / many-to-many**: lineage-specific duplications after speciation
- **In-paralogs**: paralogs from duplication AFTER the speciation event of reference
- **Out-paralogs**: paralogs from duplication BEFORE the speciation event
- **Co-orthologs**: in-paralogous genes collectively orthologous to a gene in the outgroup

## OrthoFinder Workflow

**Goal:** Infer orthologous gene groups across multiple species from their proteomes.

**Approach:** Run OrthoFinder on a directory of per-species FASTA files to perform all-vs-all DIAMOND search, gene/species tree inference, and ortholog/paralog classification, then parse the resulting orthogroups and classify by copy number pattern.

```python
'''Ortholog inference with OrthoFinder'''

import subprocess
import pandas as pd
import os


def run_orthofinder(proteome_dir, output_dir=None, threads=4):
    '''Run OrthoFinder on directory of proteomes

    Input: Directory with one FASTA file per species
    File naming: Species name derived from filename

    OrthoFinder pipeline:
    1. All-vs-all DIAMOND/BLAST search
    2. Gene tree inference per orthogroup
    3. Species tree inference (STAG/STRIDE)
    4. Gene tree rooting and reconciliation
    5. Ortholog/paralog classification via DLC model

    Key options:
    -M msa: Use MSA-based gene trees (more accurate, slower; recommended for <20 species)
    -M dendroblast: Distance-based trees (default, faster; sufficient for >20 species)
    -S diamond: Fast search (default)
    -S blast: More sensitive (use for divergent species or small proteomes)
    '''
    cmd = f'orthofinder -f {proteome_dir} -t {threads}'

    if output_dir:
        cmd += f' -o {output_dir}'

    # Add -M msa for MSA-based gene trees (more accurate for evolutionary analysis)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Output location
    if output_dir:
        results_dir = output_dir
    else:
        # OrthoFinder creates Results_MonDD in proteome_dir
        results_dir = None
        for d in os.listdir(proteome_dir):
            if d.startswith('OrthoFinder/Results_'):
                results_dir = os.path.join(proteome_dir, d)
                break

    return results_dir


def parse_orthogroups(orthogroups_file):
    '''Parse OrthoFinder Orthogroups.tsv

    Columns: Orthogroup, Species1, Species2, ...
    Values: Gene IDs (comma-separated if multiple)

    Orthogroup types:
    - Single-copy: One gene per species (ideal for phylogenomics)
    - Multi-copy: Duplications in some lineages
    - Species-specific: Genes unique to one species
    '''
    df = pd.read_csv(orthogroups_file, sep='\t')
    df = df.set_index('Orthogroup')

    orthogroups = {}
    for og_id, row in df.iterrows():
        genes = {}
        for species in df.columns:
            cell = row[species]
            if pd.notna(cell) and cell:
                genes[species] = cell.split(', ')
            else:
                genes[species] = []
        orthogroups[og_id] = genes

    return orthogroups


def classify_orthogroups(orthogroups, species_list):
    '''Classify orthogroups by copy number pattern

    Categories:
    - single_copy: Exactly one gene per species (best for phylogenomics)
    - universal: Present in all species (possibly multicopy)
    - partial: Missing from some species
    - species_specific: Only in one species
    '''
    classification = {
        'single_copy': [],
        'universal': [],
        'partial': [],
        'species_specific': []
    }

    for og_id, genes in orthogroups.items():
        present_in = [sp for sp in species_list if genes.get(sp)]
        copy_counts = [len(genes.get(sp, [])) for sp in species_list]

        if len(present_in) == 1:
            classification['species_specific'].append(og_id)
        elif len(present_in) == len(species_list):
            if all(c == 1 for c in copy_counts):
                classification['single_copy'].append(og_id)
            else:
                classification['universal'].append(og_id)
        else:
            classification['partial'].append(og_id)

    return classification


def get_single_copy_orthologs(orthogroups_file):
    '''Extract single-copy orthologs for phylogenomics

    Single-copy orthologs are ideal because:
    - Clear 1:1 relationships
    - No paralogy complications
    - Suitable for concatenated alignments
    '''
    df = pd.read_csv(orthogroups_file, sep='\t')
    df = df.set_index('Orthogroup')

    single_copy = []
    for og_id, row in df.iterrows():
        is_single = True
        for species in df.columns:
            cell = row[species]
            if pd.isna(cell) or cell == '':
                is_single = False
                break
            if ',' in str(cell):
                is_single = False
                break
        if is_single:
            single_copy.append(og_id)

    return df.loc[single_copy]
```

## Gene Trees and Reconciliation

```python
def parse_gene_trees(gene_trees_dir):
    '''Load gene trees from OrthoFinder

    Gene trees show evolutionary history within orthogroups
    Duplication/loss events inferred by species tree reconciliation
    '''
    from Bio import Phylo
    import glob

    trees = {}
    for tree_file in glob.glob(f'{gene_trees_dir}/*.txt'):
        og_id = os.path.basename(tree_file).replace('_tree.txt', '')
        trees[og_id] = Phylo.read(tree_file, 'newick')

    return trees


def identify_paralogs(orthogroup, species):
    '''Identify in-paralogs within an orthogroup

    In-paralogs: Duplications AFTER speciation (within one lineage)
    Out-paralogs: Duplications BEFORE speciation (separate orthogroups)
    Multiple genes from same species in an orthogroup = in-paralogs

    Distinguishing in- vs out-paralogs requires the species tree context
    and depends on which speciation event is being considered.
    OrthoFinder resolves this via gene tree reconciliation.
    '''
    genes = orthogroup.get(species, [])
    if len(genes) > 1:
        return {
            'species': species,
            'paralogs': genes,
            'count': len(genes)
        }
    return None


def find_co_orthologs(orthogroups, gene_id, species):
    '''Find co-orthologs of a gene

    Co-orthologs: Multiple genes in one species that are
    all orthologous to a single gene in another species

    Result of gene duplication after speciation
    '''
    for og_id, genes in orthogroups.items():
        if gene_id in genes.get(species, []):
            co_orthologs = {}
            for sp, sp_genes in genes.items():
                if sp != species and sp_genes:
                    co_orthologs[sp] = sp_genes
            return {'orthogroup': og_id, 'co_orthologs': co_orthologs}

    return None
```

## ProteinOrtho Alternative

**Goal:** Detect orthologs using ProteinOrtho as a faster alternative for many-genome comparisons.

**Approach:** Run ProteinOrtho with DIAMOND backend on multiple proteome FASTA files and parse the output table for orthologous groups with connectivity scores.

```python
def run_proteinortho(proteome_files, output_prefix, threads=4):
    '''Run ProteinOrtho for ortholog detection

    Faster than OrthoFinder for many genomes
    Uses synteny information if available

    -p=blastp+: Use DIAMOND (faster)
    -conn: Connectivity threshold (default 0.1)
    '''
    files_str = ' '.join(proteome_files)
    cmd = f'proteinortho -cpus={threads} -project={output_prefix} {files_str}'

    subprocess.run(cmd, shell=True)

    return f'{output_prefix}.proteinortho.tsv'


def parse_proteinortho(ortho_file):
    '''Parse ProteinOrtho output

    Columns: # Species, Genes, Alg.-Conn., Species1, Species2, ...
    '''
    df = pd.read_csv(ortho_file, sep='\t')

    orthogroups = {}
    for i, row in df.iterrows():
        og_id = f'OG{i:06d}'
        n_species = row['# Species']
        conn = row['Alg.-Conn.']

        genes = {}
        for col in df.columns[3:]:
            val = row[col]
            if pd.notna(val) and val != '*':
                genes[col] = val.split(',')
            else:
                genes[col] = []

        orthogroups[og_id] = {
            'genes': genes,
            'n_species': n_species,
            'connectivity': conn
        }

    return orthogroups
```

## Functional Annotation Transfer

```python
def transfer_annotation(query_gene, orthologs, annotation_db):
    '''Transfer functional annotation via orthology

    Confidence hierarchy:
    - One-to-one orthologs: Highest confidence; direct functional equivalence
    - Co-orthologs: Transfer to all, but note potential sub/neofunctionalization
    - In-paralogs (recent duplicates): Transfer with caution; function may have diverged
    - Distant orthologs (dS > 2): Lowest confidence; verify with domain conservation

    GO evidence codes:
    - ISO: Inferred from Sequence Orthology (recommended for 1:1 orthologs)
    - IBA: Inferred from Biological Aspect of Ancestor (phylogenetic propagation)
    - IEA: Inferred from Electronic Annotation (automated, lower confidence)

    Synteny context (see synteny-analysis) increases transfer confidence
    for genes in conserved genomic neighborhoods.
    '''
    annotations = []

    for species, genes in orthologs.items():
        for gene in genes:
            if gene in annotation_db:
                ann = annotation_db[gene]
                annotations.append({
                    'source_gene': gene,
                    'source_species': species,
                    'annotation': ann,
                    'evidence': 'ISO'  # Sequence orthology
                })

    return annotations
```

## Completeness Assessment

Before orthology analysis, verify proteome completeness with BUSCO or Compleasm:

```bash
# BUSCO: standard benchmark against OrthoDB single-copy orthologs
busco -i proteome.fasta -m proteins -l <lineage> -o busco_out

# Compleasm: 14x faster alternative using miniprot
compleasm run -a genome.fasta -l <lineage> -o compleasm_out
```

BUSCO categories: Complete (single-copy + duplicated), Fragmented, Missing. Expect >90% complete for well-assembled genomes. High duplication rates may indicate assembly collapse or recent WGD. Choose the most specific available lineage for the clade being compared.

## Related Skills

- comparative-genomics/synteny-analysis - Synteny-based ortholog verification and context
- comparative-genomics/positive-selection - Selection analysis on ortholog alignments
- phylogenetics/modern-tree-inference - Build species trees from single-copy orthologs
- alignment/pairwise-alignment - Align orthogroup sequences
- genome-annotation/annotation-transfer - Transfer annotations via orthology
