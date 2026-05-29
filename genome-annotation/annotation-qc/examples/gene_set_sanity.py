#!/usr/bin/env python3
'''Gene-set sanity panel for a genome annotation.

Gene count is a vanity metric; the diagnostic signal is in the mono-exonic
fraction, the protein-length distribution, and the mRNA:gene ratio. Read these
against the nearest well-annotated relative and the species' ploidy.
'''
# Reference: gffutils 0.12+, matplotlib 3.8+ | Verify API if version differs

import sys
import gffutils
import matplotlib.pyplot as plt

MONOEXONIC_FLAG = 0.30   # >30% single-exon in a vertebrate suggests unmasked TEs/pseudogenes/fragments; fungi run higher, calibrate per clade
SHORT_PROTEIN_AA = 100   # a spike below ~100 aa flags spurious/fragmented predictions (real small peptides exist; a spike is pathology)


def sanity_panel(gff_file, genome_size=None):
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy='merge')
    genes = list(db.features_of_type('gene'))
    mrnas = list(db.features_of_type(['mRNA', 'transcript']))

    exon_counts, cds_lengths = [], []
    for tx in mrnas:
        exons = list(db.children(tx, featuretype='exon'))
        exon_counts.append(len(exons))
        cds_bp = sum(c.end - c.start + 1 for c in db.children(tx, featuretype='CDS'))
        cds_lengths.append(cds_bp // 3)

    mono_frac = sum(1 for e in exon_counts if e == 1) / len(exon_counts) if exon_counts else 0
    mrna_per_gene = len(mrnas) / len(genes) if genes else 0
    short_frac = sum(1 for p in cds_lengths if p < SHORT_PROTEIN_AA) / len(cds_lengths) if cds_lengths else 0

    print(f'Genes: {len(genes)}   transcripts: {len(mrnas)}   mRNA:gene ratio: {mrna_per_gene:.3f}')
    print(f'Mono-exonic fraction: {mono_frac:.1%}')
    print(f'Proteins < {SHORT_PROTEIN_AA} aa: {short_frac:.1%}')
    if genome_size:
        coding_bp = sum(c.end - c.start + 1 for c in db.features_of_type('CDS'))
        print(f'Coding density: {coding_bp / genome_size:.1%}')

    if mrna_per_gene <= 1.001:
        print('WARNING: mRNA:gene == 1.00 -- isoform/UTR-naive; alternative-splicing and 3-prime-tag scRNA-seq analyses untrustworthy')
    if mono_frac > MONOEXONIC_FLAG:
        print(f'WARNING: mono-exonic fraction {mono_frac:.1%} high -- check repeat masking and contamination')
    if short_frac > 0.20:
        print(f'WARNING: {short_frac:.1%} of proteins < {SHORT_PROTEIN_AA} aa -- spurious/fragmented predictions likely')

    if cds_lengths:
        plt.hist([p for p in cds_lengths if p < 2000], bins=60)
        plt.axvline(SHORT_PROTEIN_AA, color='red', linestyle='--')
        plt.xlabel('protein length (aa)')
        plt.ylabel('count')
        plt.title('Protein-length distribution (healthy = unimodal ~300-450 aa)')
        plt.savefig('protein_length_distribution.png', dpi=150, bbox_inches='tight')
        plt.close()

    return {'genes': len(genes), 'mrna_per_gene': mrna_per_gene,
            'mono_exonic_fraction': mono_frac, 'short_protein_fraction': short_frac}


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: gene_set_sanity.py <annotation.gff3> [genome_size_bp]')
        sys.exit(1)
    size = int(sys.argv[2]) if len(sys.argv) > 2 else None
    sanity_panel(sys.argv[1], size)
