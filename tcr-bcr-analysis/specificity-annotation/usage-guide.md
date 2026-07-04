# Specificity Annotation - Usage Guide

## Overview

Mapping a TCR or BCR sequence to the antigen it recognizes is largely unsolved in the general case, and this skill exists to keep that honesty in the analysis. A specificity method does one of three things and none of them proves specificity: it annotates against a biased database, it finds a sequence neighborhood enriched for shared specificity, or it predicts binding only for epitopes seen in training. So a VDJdb hit or a GLIPH2 cluster label is a hypothesis with a confidence level, never a statement that a receptor is specific for antigen Y. The problem is hard because a TCR recognizes a peptide-MHC complex (MHC restriction cannot be dropped), a single TCR is massively cross-reactive, bulk sequencing gives only the unpaired beta chain, and the reference data is dominated by a handful of immunodominant epitopes with heavy HLA-A*02 and CD8/MHC-I skew. The single biggest quality lever most analyses skip is a generation-probability null: high-Pgen sequences recur across donors and match databases by chance, so any "public"/convergent/shared claim needs an OLGA Pgen or SONIA Ppost baseline.

## Prerequisites

- Python 3.9+ with pandas and numpy
- tcrdist3 for sequence clustering and meta-clonotypes: `pip install tcrdist3`
- OLGA for generation-probability nulls: `pip install olga`
- A clonotype table (AIRR rearrangement TSV or an equivalent CSV) with CDR3 amino acid, V gene, and J gene in IMGT nomenclature; ideally UMI/error-corrected upstream (see mixcr-analysis)
- For database annotation: a VDJdb / McPAS-TCR export, and the donor HLA typing for the restriction filter
- Optional: GLIPH2, GIANA, clusTCR for additional clustering; IEDB TCRMatch for similarity to characterized receptors
- Optional: IGoR and SONIA/soNNia for learning custom generative and selection models

## Quick Start

Tell the AI agent what the analysis needs:
- "Annotate my beta CDR3s against VDJdb, keeping only high-confidence, V- and HLA-concordant matches"
- "Cluster my TCRs into shared-specificity neighborhoods with tcrdist3 and report the reference used"
- "Test whether these public clonotypes are convergent or just high-Pgen using OLGA"
- "Explain why I should not call this repertoire specific for CMV from a database match"

## Example Prompts

### Database annotation
> "I have a bulk TRB repertoire and the donor's HLA type. Annotate the CDR3s against VDJdb, require V-gene concordance and that the donor carries the restricting HLA, filter to confidence score >= 1, and report matches as hypotheses with counts, not as a specificity call."

### Sequence clustering
> "Cluster these TCR beta clonotypes into shared-specificity neighborhoods with tcrdist3, build meta-clonotypes as centroid plus radius, and corroborate with a second method; report the reference repertoire and flag that clusters are enrichment, not per-receptor antigen labels."

### Generation-probability null
> "Several clonotypes are shared across my donors. Compute OLGA Pgen for each and tell me which are just easy-to-generate public sequences versus candidates for convergent antigen-driven selection."

### Overclaiming guardrail
> "A collaborator says this repertoire is specific for influenza because CDR3s match GILGFVFTL entries in a database. Explain the base-rate false-positive problem and what evidence would be needed to support that claim."

### BCR analogue
> "We see a shared IGHV3-53 CDRH3 motif across several SARS-CoV-2 convalescent donors. Is that convergent selection, and what null and confirmation are required?"

## What the Agent Will Do

1. Clarify the question (annotate, cluster, or test sharing) and the receptor context (chain, TCR vs BCR, bulk vs single-cell, donor HLA available).
2. For annotation, filter the database to a confidence threshold, join on CDR3 plus V-gene rather than CDR3 alone, and drop matches whose restricting HLA the donor does not carry.
3. For clustering, compute TCRdist neighborhoods or run GLIPH2/GIANA/clusTCR, build meta-clonotypes, and corroborate with a second method while reporting the reference repertoire.
4. For any sharing/convergence claim, compute a generation-probability null (OLGA Pgen, or SONIA Ppost) and down-weight high-Pgen sequences.
5. Frame ML predictor output as reliable only for well-trained epitopes and unreliable on novel epitopes, deferring predictor detail to the immunoinformatics skills.
6. Report every result as an annotation or hypothesis with its confidence and required concordances, reserving "specific" for experimental confirmation.

## Tips

- A database match or a cluster label is a hypothesis, never a specificity call; keep the word "specific" for tetramer/dextramer or functional confirmation.
- Never annotate on a bare CDR3-beta match: the number of chance hits scales with repertoire size, database size, and match permissiveness. Require V-gene and HLA concordance and a confidence filter.
- Attach a Pgen null to every "public"/convergent/shared claim. Most publicity is high generation probability plus convergent recombination, not antigen selection.
- Beta-only specificity is low-confidence because the paired alpha chain carries much of the signal; prefer single-cell paired data where the claim matters.
- Run at least two clustering methods and report agreement; there is no gold standard and benchmarks disagree by dataset and epitope. GLIPH2 clusters are reference- and parameter-sensitive.
- State the biases whenever reporting coverage or accuracy: epitope imbalance, HLA-A*02 skew, and CD8/MHC-I dominance inflate aggregate metrics.
- For BCR, cluster clones (shared V, J, junction length plus within-partition distance) before any convergence analysis, and remember conformational epitopes limit sequence-only prediction.

## Related Skills

- mixcr-analysis - Produce clonotype tables to annotate
- scirpy-analysis - Paired single-cell clonotypes for specificity work
- vdjtools-analysis - Public-clonotype context and overlap
- immunoinformatics/tcr-epitope-binding - ML epitope-binding prediction details
- immunoinformatics/mhc-binding-prediction - Upstream pMHC restriction
- immunoinformatics/neoantigen-prediction - Neoantigen-directed specificity
