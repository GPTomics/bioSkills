---
name: bio-genome-annotation-annotation-transfer
description: Transfer gene annotations between genome assemblies using Liftoff for same-species annotation liftover and MiniProt for cross-species protein-to-genome alignment. Enables rapid annotation of new assemblies using existing reference annotations. Use when annotating a new assembly of a species with an existing reference annotation or mapping annotations across related species.
tool_type: cli
primary_tool: Liftoff
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Annotation Transfer

**"Transfer annotations from a reference to my new assembly"** → Map gene models from a well-annotated reference genome onto a new assembly using coordinate liftover or protein-to-genome alignment.
- CLI: `liftoff -g reference.gff -o target.gff ref.fa target.fa` (same species), `miniprot ref.mpi target.fa` (cross-species)

Transfer gene annotations from a reference genome to a new assembly (same species with Liftoff) or across species (with MiniProt protein-to-genome alignment). Faster and more consistent than de novo prediction when a high-quality reference annotation exists.

## Liftoff (Same-Species Transfer)

Liftoff maps annotations from a reference genome to a target assembly using Minimap2 alignments. Ideal for transferring annotations between different assemblies of the same species.

### Basic Usage

```bash
# Transfer annotations from reference to target
liftoff \
    -g reference_annotation.gff3 \
    -o lifted_annotation.gff3 \
    -u unmapped_features.txt \
    -dir liftoff_intermediates \
    -p 16 \
    target_assembly.fasta \
    reference_genome.fasta
```

### Key Options

| Option | Description |
|--------|-------------|
| `-g` | Reference annotation (GFF3 or GTF) |
| `-o` | Output annotation file |
| `-u` | File listing unmapped features |
| `-dir` | Directory for intermediate files |
| `-p` | CPU threads |
| `-sc` | Coverage threshold (default: 0.5; fraction of ref feature aligned) |
| `-s` | Sequence identity threshold (default: 0.5) |
| `-a` | Alignment coverage cutoff (default: 0.5) |
| `-copies` | Look for extra gene copies in target |
| `-exclude_partial` | Exclude partially mapped genes |
| `-chroms` | Chromosome name mapping file (tab-separated: ref\ttarget) |

### Strict Parameters for High-Quality Transfer

```bash
# Stricter thresholds for closely related assemblies
# sc 0.95: 95% of reference feature must align
# s 0.90: 90% sequence identity required
liftoff \
    -g reference.gff3 \
    -o lifted.gff3 \
    -u unmapped.txt \
    -dir liftoff_tmp \
    -sc 0.95 \
    -s 0.90 \
    -exclude_partial \
    -p 16 \
    target.fasta \
    reference.fasta
```

### With Chromosome Name Mapping

```bash
# Create chromosome mapping file (tab-separated)
# ref_chr1    target_scaffold_1
# ref_chr2    target_scaffold_2
liftoff \
    -g reference.gff3 \
    -o lifted.gff3 \
    -chroms chrom_map.txt \
    -p 16 \
    target.fasta \
    reference.fasta
```

### Output

The output GFF3 contains transferred annotations with additional attributes:

| Attribute | Description |
|-----------|-------------|
| `coverage` | Fraction of reference feature aligned |
| `sequence_ID` | Sequence identity of alignment |
| `extra_copy_number` | Copy number if `-copies` used |
| `valid_ORF` | Whether transferred CDS has valid ORF |

## LiftOn (Newer Successor)

LiftOn improves on Liftoff by combining Liftoff liftover with MiniProt protein alignment to correct gene models that do not transfer cleanly.

```bash
# LiftOn combines Liftoff + MiniProt
lifton \
    -g reference.gff3 \
    -o lifton_annotation.gff3 \
    -ref reference.fasta \
    -p 16 \
    target.fasta
```

## MiniProt (Cross-Species Protein Alignment)

MiniProt aligns protein sequences to a genome with splicing awareness. Ideal for cross-species annotation transfer using proteins from related species.

### Basic Usage

```bash
# Index target genome
miniprot -t 16 -d target.mpi target_assembly.fasta

# Align proteins to genome
miniprot -t 16 --gff target.mpi reference_proteins.faa > miniprot_alignments.gff
```

### Key Options

| Option | Description |
|--------|-------------|
| `-t` | CPU threads |
| `-d` | Build index database |
| `--gff` | Output in GFF3 format |
| `--gtf` | Output in GTF format |
| `-G` | Max intron size (default: 200000) |
| `-S` | Output alignment score |
| `--outs` | Output secondary alignments (for paralogs) |
| `-C` | Min alignment coverage (0-1; default: 0.5) |
| `-k` | K-mer size for indexing |

### Cross-Species Transfer

```bash
# Use proteins from closely related species
# -G: Adjust max intron size based on target species
# Vertebrates: -G 500000; Insects: -G 50000; Fungi: -G 5000
miniprot -t 16 --gff -G 500000 target.mpi related_species_proteins.faa > cross_species.gff
```

### Convert MiniProt GFF to Gene Models

```python
import gffutils

def miniprot_gff_to_gene_models(miniprot_gff, output_gff):
    '''Convert MiniProt alignment GFF to standard gene models.

    MiniProt outputs mRNA features with CDS children.
    This adds gene-level parent features for compatibility.
    '''
    db = gffutils.create_db(miniprot_gff, ':memory:', merge_strategy='merge')

    gene_id = 0
    with open(output_gff, 'w') as out:
        out.write('##gff-version 3\n')
        for mrna in db.features_of_type('mRNA'):
            gene_id += 1
            gene_line = f'{mrna.seqid}\tMiniProt\tgene\t{mrna.start}\t{mrna.end}\t{mrna.score}\t{mrna.strand}\t.\tID=mpgene_{gene_id}\n'
            mrna_line = f'{mrna.seqid}\tMiniProt\tmRNA\t{mrna.start}\t{mrna.end}\t{mrna.score}\t{mrna.strand}\t.\tID={mrna.id};Parent=mpgene_{gene_id}\n'
            out.write(gene_line)
            out.write(mrna_line)
            for child in db.children(mrna):
                child_line = f'{child.seqid}\tMiniProt\t{child.featuretype}\t{child.start}\t{child.end}\t{child.score}\t{child.strand}\t{child.frame}\tParent={mrna.id}\n'
                out.write(child_line)

    return output_gff
```

### Distinguish from Orthology-Based Transfer

MiniProt performs protein-to-genome alignment, which maps protein sequences to genomic coordinates with intron prediction. This is different from orthology-based transfer (see comparative-genomics/ortholog-inference), which identifies evolutionary relationships between gene families without genome alignment.

## Quality Assessment

**Goal:** Evaluate annotation transfer quality by comparing gene/transcript counts and validating that transferred CDSs have intact open reading frames.

**Approach:** Count genes and transcripts in both reference and transferred GFF files to compute a transfer rate, then extract each transferred CDS sequence from the target assembly and check for valid start codon, single stop codon, and correct frame.

```python
import gffutils
import pandas as pd

def compare_annotations(reference_gff, transferred_gff):
    '''Compare reference and transferred annotations for QC.'''
    ref_db = gffutils.create_db(reference_gff, ':memory:', merge_strategy='merge')
    tgt_db = gffutils.create_db(transferred_gff, ':memory:', merge_strategy='merge')

    ref_genes = list(ref_db.features_of_type('gene'))
    tgt_genes = list(tgt_db.features_of_type('gene'))

    ref_mrnas = list(ref_db.features_of_type(['mRNA', 'transcript']))
    tgt_mrnas = list(tgt_db.features_of_type(['mRNA', 'transcript']))

    stats = {
        'ref_genes': len(ref_genes),
        'transferred_genes': len(tgt_genes),
        'transfer_rate': len(tgt_genes) / len(ref_genes) if ref_genes else 0,
        'ref_transcripts': len(ref_mrnas),
        'transferred_transcripts': len(tgt_mrnas),
    }

    print('=== Annotation Transfer QC ===')
    print(f'Reference genes: {stats["ref_genes"]}')
    print(f'Transferred genes: {stats["transferred_genes"]}')
    print(f'Transfer rate: {stats["transfer_rate"]:.1%}')
    print(f'Reference transcripts: {stats["ref_transcripts"]}')
    print(f'Transferred transcripts: {stats["transferred_transcripts"]}')

    # Transfer rate > 95% is excellent for same-species liftover
    # Transfer rate > 80% is typical for closely related species
    # Transfer rate < 70% suggests distant species or assembly issues
    if stats['transfer_rate'] > 0.95:
        print('Quality: Excellent (>95% transfer rate)')
    elif stats['transfer_rate'] > 0.80:
        print('Quality: Good (>80% transfer rate)')
    else:
        print('Quality: Low transfer rate - check assembly quality or species distance')

    return stats

def check_transferred_orfs(transferred_gff, target_fasta):
    '''Check how many transferred CDSs have valid open reading frames.'''
    from Bio import SeqIO

    genome = SeqIO.to_dict(SeqIO.parse(target_fasta, 'fasta'))
    db = gffutils.create_db(transferred_gff, ':memory:', merge_strategy='merge')

    valid, invalid, total = 0, 0, 0
    for cds in db.features_of_type('CDS'):
        total += 1
        seq = genome[cds.seqid].seq[cds.start - 1:cds.end]
        if cds.strand == '-':
            seq = seq.reverse_complement()

        protein = seq.translate()
        if protein.startswith('M') and protein.endswith('*') and protein.count('*') == 1:
            valid += 1
        else:
            invalid += 1

    print(f'\n=== ORF Validation ===')
    print(f'Total CDSs: {total}')
    print(f'Valid ORFs: {valid} ({valid/total:.1%})')
    print(f'Invalid ORFs: {invalid} ({invalid/total:.1%})')

    return valid, invalid, total
```

## Troubleshooting

### Many Unmapped Features with Liftoff
- Check assembly contiguity (fragmented assemblies lose features at contig boundaries)
- Relax thresholds: `-sc 0.5 -s 0.5`
- Verify chromosome naming consistency

### MiniProt Misses Short Genes
- Reduce minimum alignment coverage: `-C 0.3`
- Check that protein sequences include short ORFs

### Invalid ORFs After Transfer
- Assembly may have variants causing frameshifts
- Try LiftOn which combines Liftoff + MiniProt for correction
- Consider re-predicting genes de novo in problem regions

## Related Skills

- eukaryotic-gene-prediction - De novo prediction alternative
- comparative-genomics/ortholog-inference - Orthology-based functional transfer
- comparative-genomics/synteny-analysis - Synteny context for annotation transfer
- genome-intervals/gtf-gff-handling - Parse and manipulate transferred annotations
