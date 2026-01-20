# Protein Inference

## Overview

Protein inference determines which proteins are present in a sample based on identified peptides, accounting for shared peptides between homologous proteins.

## The Problem

- One peptide can match multiple proteins (paralogs, isoforms)
- Proteins share peptides due to sequence similarity
- Cannot always determine which protein is present

## Key Concepts

### Unique Peptides
Peptides mapping to only one protein - strongest evidence

### Razor Peptides
Shared peptides assigned to protein with most total peptides

### Protein Groups
Proteins with indistinguishable peptide evidence grouped together

## Inference Strategies

| Method | Description |
|--------|-------------|
| Parsimony | Minimum protein set explaining all peptides |
| Occam's Razor | Assign shared peptides to "winning" protein |
| Probabilistic | ProteinProphet, EPIFANY - probability-based |
| All peptides | Report all possible proteins (most inclusive) |

## Protein Group Reporting

MaxQuant format example:
```
Protein IDs: P12345;Q67890
Majority protein IDs: P12345
Gene names: GENE1;GENE2
Unique peptides: 3
Razor + unique peptides: 5
```

## FDR Control

Protein-level FDR is different from peptide-level:
- **Peptide FDR**: 1% means 1% of peptides are false
- **Protein FDR**: 1% means 1% of protein groups are false
- Protein FDR typically 1-5%

## Best Practices

1. Report protein groups, not just lead proteins
2. Use unique peptides for confident identification (>= 2)
3. Apply separate protein-level FDR
4. Be cautious with single-peptide identifications
5. Document inference method used

## Common Outputs

- **proteinGroups.txt**: MaxQuant protein-level output
- **Protein.tsv**: FragPipe protein output
- **mzIdentML**: Standard format with protein groups
