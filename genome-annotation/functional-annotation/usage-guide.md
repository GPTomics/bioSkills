# Functional Annotation - Usage Guide

## Overview

Assign GO terms, KEGG orthologs, Pfam/InterPro domains, EC numbers, and product names to predicted proteins using eggNOG-mapper (orthology), InterProScan (domain signatures), and KofamScan (KEGG), routing specialized functions to dbCAN/antiSMASH/AMRFinderPlus/SignalP. The orienting principle: a functional label is a transferred hypothesis whose confidence is capped by its weakest link, and the goal is not "maximally annotated" but "honestly tiered" - prefer curated/orthology donors, demote specificity as identity/coverage fall, separate "domain present" from "function known", and keep the honest hypothetical fraction rather than over-annotating under pressure.

## Prerequisites

```bash
# eggNOG-mapper (pin the version for reproducibility) + DB
conda install -c bioconda eggnog-mapper=2.1.15
download_eggnog_data.py --data_dir /path/to/eggnog_db -y   # ~44 GB (DIAMOND DB installed by default; -D skips it)

# InterProScan (Java 11+, tens of GB of member-DB data)
conda install -c bioconda interproscan

# KEGG (free route) + GFF surgery
conda install -c bioconda kofamscan agat

# Python utilities
pip install pandas biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Add functional annotations to my predicted proteins with eggNOG-mapper and InterProScan"
- "Annotate KEGG orthologs with KofamScan (KEGG is paywalled, so use the free route)"
- "How much of my genome is honestly hypothetical, and is my named fraction too high for this taxon?"
- "Annotate CAZymes / AMR genes / signal peptides correctly, not from generic Pfam hits"

## Example Prompts

### Orthology and Domains

> "Run eggNOG-mapper with tax_scope auto and InterProScan with GO terms, then merge keeping each tool's evidence in its own namespace and preferring the orthology product name."

> "My proteins are from a deep-branching uncultured lineage - annotate honestly and tell me the hypothetical fraction rather than forcing names."

### Specialized Functions (Route to the Right Tool)

> "Annotate CAZymes with dbCAN (consensus of >=2 tools) - don't read substrate off a plain GH domain."

> "Find AMR genes including resistance point mutations with AMRFinderPlus, not generic homology."

### Interpretation

> "I'm about to run GO enrichment on these IEA annotations - what's the circularity risk and how do I control for it?"

> "A distant homolog got a specific substrate name - should I demote it to the superfamily level?"

## What the Agent Will Do

1. Verify the input protein FASTA and record DB versions (eggNOG, InterPro/Pfam, KofamScan profiles)
2. Run eggNOG-mapper (orthology) with `--tax_scope auto` and InterProScan (domains + GO + pathways)
3. Reconcile, preserving provenance per source; prefer curated/orthology names; demote specificity as identity/coverage fall
4. Route CAZymes/BGCs/AMR/signal peptides to dbCAN/antiSMASH/AMRFinderPlus/SignalP
5. Report annotation coverage AND the honest hypothetical fraction, flagging if the named fraction is too high for the taxon
6. Caveat any downstream GO enrichment for IEA circularity and study bias

## Tips

- **Honestly tiered beats maximally annotated** - 40% hypothetical correctly tiered is better than 95% named with half wrong. Resist pressure to "annotate everything"; each forced name is a future percolating error.
- **Prefer curated/orthology donors** - eggNOG OG consensus and Swiss-Prot over TrEMBL/NR best-hits; the cascade lives in the bulk databases.
- **Demote specificity as similarity falls** - full EC -> partial `1.1.1.-`; specific name -> superfamily; whole-protein -> per-domain. There is no safe identity in mechanistically-diverse superfamilies (enolase, amidohydrolase, TIM-barrel).
- **Domain != function** - ~10% of kinase-domain proteins are dead pseudokinases; report architecture and activity as separate claims; weight GO Molecular Function over Biological Process.
- **Reason in bits-per-residue + bidirectional coverage**, not raw e-value (which scales with DB size). One-domain coverage justifies a domain-level claim only.
- **GO enrichment on IEA is partly circular** - it can measure the interpro2go mapping and study popularity, not biology. Match background to pipeline, pin the GO release, prefer non-IEA where it exists, and say it is annotation-derived.
- **KEGG is paywalled** - use KofamScan/KofamKOALA or eggNOG's `KEGG_ko`; never scrape the REST API for whole genomes. A "complete module" is a reconstruction, not proof of flux.
- **Record provenance** - method, donor, evidence code, identity/coverage/bitscore, DB version; it is the only thing that stops a transfer becoming the next genome's "fact."

## Related Skills

- genome-annotation/prokaryotic-annotation - Bakta/PGAP product names before functional layers
- genome-annotation/eukaryotic-gene-prediction - Produces proteins for annotation
- genome-annotation/annotation-qc - Annotation-coverage and hypothetical-fraction sanity
- pathway-analysis/go-enrichment - Downstream GO enrichment (mind IEA circularity)
- pathway-analysis/kegg-pathways - KEGG pathway mapping
