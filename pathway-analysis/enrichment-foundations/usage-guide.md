# Enrichment Foundations

The conceptual spine for the pathway-analysis category: how to choose the enrichment method before any tool runs, which null hypothesis each tool computes, and how to decide whether an enrichment result is trustworthy. The per-tool skills (go-enrichment, gsea, kegg-pathways, reactome-pathways, wikipathways, enrichment-visualization) defer to this skill for the method-selection and trustworthiness decisions.

## Overview

Functional enrichment asks whether a gene set or pathway is over-represented in, or coordinately shifted across, an experiment's genes. The result is not a discovery - it is a claim about a gene list versus a specific database version, computed under a null hypothesis the analyst usually did not choose on purpose. This skill makes that choice explicit. It maps the input shape to the enrichment generation (over-representation analysis for a list, gene set enrichment for a ranking, pathway topology for a signed graph), names which null each tool tests (competitive vs self-contained, gene-sampling vs subject-sampling), and runs the trustworthiness checklist (correct testable-gene universe, FDR, redundancy collapse, leading-edge multifunctionality check, version reporting). The agent should make the method-selection and trustworthiness call from the SKILL.md alone, then hand the chosen method to the relevant tool skill.

## Prerequisites

```r
BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
# The tool skills add: enrichplot, ReactomePA, rWikiPathways, goseq, msigdbr, limma (CAMERA/ROAST)
```

Conceptual prerequisites the analyst must settle first:

- The input determines the method. A pre-selected gene list with no ranking calls for over-representation analysis (ORA); a per-gene statistic for nearly all genes calls for gene set enrichment (GSEA); a ranked list plus a curated signed signaling graph calls for pathway topology (SPIA). The ranking statistic itself comes from differential-expression/de-results.
- The background universe is the analysis, not the gene list. For ORA the universe is the genes that passed the same filter that produced the list (the testable genes), not the whole genome. The clusterProfiler default universe is silently all-annotated-genes.
- GSEA needs a NAMED numeric vector sorted in DECREASING order, not a filtered list. Every clusterProfiler/fgsea GSEA is a preranked, gene-permutation, competitive null - a discovery screen, not a calibrated test.
- Gene IDs must match the method (OrgDb keyType for enrichGO, kegg-id/ENTREZ for enrichKEGG, ENTREZ for ReactomePA/enrichWP); KEGG and WikiPathways query live databases and are not reproducible across releases.

## Quick Start

Tell your AI agent what you want to do:

- "Should I run ORA or GSEA on this data?"
- "Help me build a defensible background universe for my gene list"
- "Which null hypothesis does clusterProfiler GSEA actually test?"
- "Is this enrichment result trustworthy, or is it an artifact?"
- "Which pathway database should I use for this question?"

## Example Prompts

### Choosing the method

> "I have a full DESeq2 result with the test statistic for every gene and no clear significance cutoff. Should I run over-representation analysis or GSEA, and why?"

> "My gene list is 180 CRISPR screen hits with no per-gene ranking. Which enrichment generation is valid here and which databases make sense?"

### Defending the result

> "I ran enrichGO with the defaults and got 60 significant biological-process terms. Walk me through why that number might be misleading and what to check before reporting it."

> "A reviewer says my GSEA p-values are inflated because I used gene permutation. Explain whether they are right and what I should do about it."

### Database and null choice

> "Explain the difference between a competitive and a self-contained null, and tell me which one clusterProfiler GSEA computes."

> "I want to compare metabolic pathways across three conditions. Which database and which method fit, and where does the multiple-testing burden differ from GO?"

## What the Agent Will Do

1. Inspect the input shape: is there a ranked statistic for nearly all genes (GSEA), a pre-selected list plus a defensible background (ORA), or a signed signaling topology with fold changes (pathway topology)?
2. Name the enrichment generation from the input, not from preference, and route the actual run to the relevant tool skill (go-enrichment, gsea, kegg-pathways, reactome-pathways, wikipathways).
3. For an ORA, construct the testable-gene universe from the same filter that produced the list and contrast it with the all-genome default to show the difference.
4. State which null the chosen tool computes, and flag that every clusterProfiler GSEA is the inter-gene-correlation-uncorrected competitive null - a discovery screen, with CAMERA or ROAST as the calibrated alternative.
5. Run the trustworthiness checklist: correct universe, FDR-adjusted p reported, redundancy collapsed before interpretation, leading-edge genes inspected for a multifunctionality artifact, and the database release and access date recorded.
6. Recommend confirming a headline finding with a second method class, since no single method is universally best.

## Tips

- The first decision is the generation, and the input sets it. A bare list goes to ORA, a full ranking goes to GSEA, a ranking plus a signed graph goes to SPIA. Do not pick by habit.
- The background decides ORA significance. If the background is the whole genome when only the expressed genes were testable, the analysis measures expression bias, not enrichment.
- Treat a preranked-GSEA or ORA p-value as a discovery screen, not a calibrated competitive test. If calibration matters, run CAMERA; if the question is "is this set involved at all," run ROAST or fry on the expression matrix.
- A list of 40 significant GO terms is usually about 3 biological stories told 40 times. Collapse redundancy with semantic similarity before interpreting, and report term clusters.
- Inspect the leading-edge or geneID core genes. If the same 3-5 multifunctional hub genes explain the top 20 sets, that is one finding, not twenty.
- "Cancer / immune / apoptosis pathways enriched" is often the field's null result, because genes are annotated in proportion to how studied they are, not how important.
- Record package versions, the database release and access date, the ranking metric, the p-adjust method, and the universe in a provenance block, or the result is unreproducible.
- No single method is universally best. A robust finding survives a second method class; do not claim one method beats another by a number.

## Related Skills

- go-enrichment - ORA mechanics: enrichGO, the GO DAG, the background universe, GOseq length bias
- gsea - FCS mechanics: ranking metric, running-sum ES, leading edge, CAMERA/ROAST, GSVA/ssGSEA
- kegg-pathways - KEGG pathway/module enrichment and SPIA topology
- reactome-pathways - Reactome curated-pathway ORA and GSEA
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - Dot/cnet/emap plots and redundancy collapse
- differential-expression/de-results - Source of the gene list and the ranking statistic
- experimental-design/multiple-testing - FDR / p.adjust method theory
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
