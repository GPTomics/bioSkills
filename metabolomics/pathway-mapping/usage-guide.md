# Pathway Mapping Usage Guide

## Overview

Pathway mapping places metabolomics results in biochemical context through over-representation (ORA), metabolite-set enrichment (MSEA/QEA), mummichog/PSEA on raw m/z features, and network-diffusion enrichment (FELLA). The central failure it guards against: enrichment destroys annotation uncertainty (a tentative ID becomes a confident p-value), the background set silently controls every result, topology "impact" is a hub artifact, and a steady-state pool size is not flux. The honest output is a hypothesis about network activity, conditional on the chosen annotations, background, database boundary, and ionization settings.

## Prerequisites

```bash
# R packages
# MetaboAnalystR (GitHub): see https://github.com/xia-lab/MetaboAnalystR
BiocManager::install("FELLA")
BiocManager::install("KEGGREST")
```

Conceptual prerequisites: know whether the metabolites are confidently identified (KEGG/HMDB IDs -> ORA/MSEA) or are raw m/z features with no IDs (-> mummichog/PSEA); know the MSI confidence level (Schymanski 5-level scale) of the driving compounds; be able to state the background set (the assay-coverage metabolome) in one sentence; know the ionization mode and ppm of the run.

## Quick Start

Tell your AI agent what you want to do:
- "Run ORA on my list of identified significant metabolites against KEGG human pathways with an assay-specific background"
- "Predict pathway activity from my untargeted LC-MS feature table using mummichog with the full feature table as background"
- "Find the enzymes and reactions linking my metabolites using FELLA network diffusion"
- "Tell me whether my pathway result is being driven by a single hub metabolite"

## Example Prompts

### Identified-Compound Enrichment
> "Test which KEGG pathways are over-represented in my 18 confidently identified metabolites, using only the ~300 compounds my assay can detect as the background."
> "Run quantitative MSEA on my ranked metabolite fold changes instead of a cutoff-based ORA."

### Raw-Feature Activity Prediction
> "I have an untargeted negative-mode LC-MS peak table with m/z, p-value, and t-score but no IDs; predict perturbed pathways with mummichog and make sure the full table is the background."
> "Run integrated mummichog + GSEA PSEA at 5 ppm and report the predicted-active pathways as activity, not metabolite identifications."

### Mechanism and Sanity-Checking
> "Use FELLA diffusion to return the intermediate enzymes and reactions linking my KEGG compounds, and list which compounds did not map."
> "Check whether my top pathway is significant only because of L-alanine centrality, and re-run the enrichment across KEGG and SMPDB to see if it is robust."

## What the Agent Will Do

1. Determine the input type (identified IDs vs raw m/z) and pick ORA/MSEA vs mummichog/PSEA accordingly.
2. Construct an explicit background: the assay-coverage metabolome for ORA, or the full feature table (R_all) for mummichog.
3. Map IDs / features, declare ionization mode and ppm, and report mapping coverage before trusting any p-value.
4. Run the enrichment, treating topology "impact" only as a secondary tiebreaker.
5. Check whether a single hub/cofactor drives the result and whether the call is robust across databases.
6. State the regime ("predicted activity" vs "measured enrichment") and the MSI level of the driving compounds, and refuse flux/activity language from concentrations.

## Tips

- The background set is the null hypothesis made concrete; with the correct assay-specific background, many "significant" pathways disappear after FDR.
- Mummichog's #1 user error is supplying only significant features; its permutation null must draw from the entire feature table.
- The mummichog query p-cutoff defaults looser than 0.05 (often ~0.2) so the query is large enough to score; document the value.
- "TCA cycle / amino-acid metabolism enriched" is closer to a null result than a finding -- it is what the database can map; always report coverage.
- Pathway granularity moves p-values more than multiple-testing correction does; prefer cross-database consensus over one library.
- Pool size is not flux; downgrade "pathway activated" to "members co-varied with phenotype" and reserve activity claims for stable-isotope tracing (SIRM/13C-MFA).

## Related Skills

- metabolomics/metabolite-annotation - Annotation confidence levels (MSI) that feed ORA/MSEA and set the interpretive ceiling
- metabolomics/statistical-analysis - Upstream differential testing that produces the significant compound or feature list
- pathway-analysis/go-enrichment - Gene-set over-representation concepts
- pathway-analysis/gsea - Ranked-list enrichment concepts
- multi-omics-integration/mofa-integration - Joint gene+metabolite integration and its coverage-asymmetry traps
