# TCR-Epitope Binding - Usage Guide

## Overview

Infer or annotate TCR antigen specificity with the honesty the field requires: unsupervised clustering plus database lookup is the defensible workhorse, while de-novo binding prediction for unseen epitopes essentially does not work and must be framed as a validation-bound hypothesis. The skill routes each request to the appropriate task.

## Prerequisites

```bash
pip install tcrdist3 pandas scipy   # clustering + meta-clonotypes
# GLIPH2, clusTCR, GIANA are separate repos; ERGO-II / NetTCR / pMTnet are separate
# repos with pretrained weights. VDJdb, IEDB, and McPAS-TCR are downloadable databases.
```

## Quick Start

Tell your AI agent what you want to do:
- "Cluster my TCR repertoire and annotate clusters with known specificities"
- "Look up these CDR3b sequences in VDJdb at high confidence"
- "Group these TCRs by likely shared antigen, keeping HLA as a covariate"
- "Can I predict what this neoantigen-specific TCR binds de novo?"

## Example Prompts

### Clustering and Annotation

> "Run tcrdist3 on this cohort, cluster, and propagate VDJdb labels to clusters with a known member"

> "Use GLIPH2 to find specificity groups and predict the restricting HLA"

### Database Lookup

> "Match these CDR3b against VDJdb and McPAS-TCR, filtering on confidence, and report the epitope and HLA"

### De-Novo Caveats

> "This is a patient-specific neoantigen with no database neighbors - what can I actually claim?"

> "A model reports AUC 0.95 for TCR-epitope binding - should I trust it?"

## What the Agent Will Do

1. Ask whether known specificities exist (tetramer sort / database hits) or the epitope is truly de-novo
2. For known specificities: cluster (tcrdist3/GLIPH2/clusTCR) within an HLA-coherent cohort and annotate by lookup
3. Propagate labels by guilt-by-association to clusters containing a known member, with confidence and HLA stated
4. For de-novo: rank candidates with pMTnet/PanPep explicitly as a hypothesis and name the validation experiment
5. Interrogate any supervised model's negative-sampling scheme and demand epitope-disjoint evaluation
6. Refuse to let a per-pair probability substitute for a tetramer/functional assay

## Tips

- **Clustering is honest; prediction is not** - cluster within a dataset; do not extrapolate to unseen epitopes
- **No true negative set** - manufactured negatives dominate reported AUC; read the negative-sampling sentence first
- **CDR3b-only is information-poor** - alpha chain and V/J carry signal; paired data is the real substrate
- **HLA is a confounder** - the same CDR3 on a different allele is a different specificity; keep cohorts HLA-coherent
- **Metrics that lie** - per-epitope AUC averaged over seen epitopes hides the novel-epitope collapse; demand disjoint splits
- **10x labels are noisy** - dextramer calls are thresholds, not gold; require replicate/donor concordance
- **Structure is for refinement** - TCRdock ranks/rationalizes candidates; it does not screen unseen pairs

## Related Skills

- immunoinformatics/mhc-binding-prediction - the pMHC context a TCR recognizes
- immunoinformatics/neoantigen-prediction - de-novo neoantigen TCRs are the unseen-epitope case where prediction fails
- tcr-bcr-analysis/mixcr-analysis - upstream TCR repertoire extraction
- single-cell/clustering - paired single-cell TCR data and embedding-based grouping
