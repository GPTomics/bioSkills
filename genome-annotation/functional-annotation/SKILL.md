---
name: bio-genome-annotation-functional-annotation
description: Assigns GO terms, Pfam/InterPro domains, KEGG orthologs, EC numbers, and product names to predicted proteins using eggNOG-mapper (orthology), InterProScan (domain signatures), and KofamScan (KEGG), routing specialized functions to dbCAN/antiSMASH/AMRFinderPlus/SignalP. Covers the orthology-vs-domain-vs-homology paradigms, the annotation-error percolation cascade, domain-presence-is-not-function, GO IEA circularity in enrichment, evidence tiering, and bit-score/coverage thresholds. Use when adding functional annotation to predicted genes, choosing between eggNOG-mapper and InterProScan, or judging how much to trust a functional label.
tool_type: cli
primary_tool: eggNOG-mapper
---

## Version Compatibility

Reference examples tested with: eggNOG-mapper 2.1.15 (pin for reproducibility), InterProScan 5.66+, KofamScan 1.3+, pandas 2.2+, AGAT 1.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

Annotation content tracks **database release**: record the eggNOG DB version, InterPro/Pfam release, and KEGG/KofamScan profile date, and note whether InterProScan used the EBI precalculated lookup service. eggNOG-mapper v3 is under testing (not production) - pin v2.1.15. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Functional Annotation

**"Functionally annotate my predicted proteins"** -> Transfer GO/KEGG/Pfam/EC/product labels from characterized proteins by orthology and domain signatures, attaching a confidence tier and provenance to each.
- CLI: `emapper.py -i proteins.faa --itype proteins -m diamond` (eggNOG-mapper), `interproscan.sh -i proteins.faa -f TSV,GFF3 -goterms -pa` (InterProScan)

## The Single Most Important Modern Insight -- Annotation Is a Propagated Hypothesis, Not a Measurement

Almost every label on a new genome is *transferred* by homology/orthology/ML from a small island of experimentally characterized proteins. The transfer chain is lossy and self-reinforcing - it behaves like a **percolation cascade** (Gilks 2002 *Bioinformatics* 18:1641): an over-specific name assigned in year 0, deposited with **no record that it was transferred**, becomes the nearest hit for the next genome, whose label becomes evidence for the next. By the time a query reaches NR, "number of hits agreeing" measures *how far an error spread*, not correctness. Schnoes 2009 (*PLoS Comput Biol* 5:e1000605) found misannotation reaching ~80% in bulk databases (TrEMBL/NR) and near-zero in curated Swiss-Prot - the gap *is* the curation. Three load-bearing consequences:

1. **The goal is not "maximally annotated" - it is "honestly tiered."** A genome that is 40% "hypothetical protein" with the rest correctly tiered by evidence is a better scientific object than one 95% named with half the names wrong. Prefer curated/orthology donors (Swiss-Prot, eggNOG OG consensus) over best-hits, and **demote specificity as identity/coverage fall** (full EC -> partial `1.1.1.-`; specific name -> superfamily; whole-protein -> per-domain). PI/reviewer pressure to "annotate everything" manufactures the next genome's percolating error.
2. **"Domain present" and "function known" are different claims.** A Pfam hit reports architecture, not activity - ~10% of the human kinome are catalytically dead pseudokinases that carry a confident "protein kinase" domain. Moonlighting (GAPDH), promiscuity, and mechanistically-diverse superfamilies (enolase, amidohydrolase, HAD, TIM-barrel: shared fold, divergent substrate) make this a first-order effect. **Fold conservation != function conservation** any more than sequence does - so structure-based transfer (Foldseek) inherits the same trap with *higher* false confidence.
3. **Record provenance on every label** (method, donor, donor evidence code, identity/coverage/bitscore, DB version). That is the only thing that stops the provenance-amnesia step that turns a transfer into a "fact."

## Tool Taxonomy

| Paradigm | Tool | Mechanism | Failure mode |
|----------|------|-----------|--------------|
| Orthology | eggNOG-mapper | seed-ortholog -> orthologous group -> consensus transfer | tax-scope sensitive; HGT/xenologs break the orthology assumption |
| Domain/signature | InterProScan | profile HMMs/matrices -> integrated InterPro entries | a domain implies a capability, not the substrate; broad families uninformative |
| KEGG ortholog | KofamScan | per-KO HMMs + adaptive thresholds | KO assignment, not pathway proof |
| Homology best-hit | DIAMOND vs Swiss-Prot | top-hit similarity, transfer label | best-hit != ortholog; transitive error propagation |
| ML / structure | DeepGO, DeepFRI, Foldseek | learned sequence/structure -> GO | low precision; ontology terms not products; reaches twilight zone only |

**Default workhorse pair:** eggNOG-mapper + InterProScan (orthogonal evidence: orthology vs signatures), reconciled afterward. Add KofamScan if KEGG pathway reconstruction is the goal (its adaptive per-KO thresholds are stricter than eggNOG's `KEGG_ko`). DIAMOND-vs-Swiss-Prot is the cheap product-name layer; never use it alone for GO.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bacterial isolate | Bakta/PGAP product names + eggNOG-mapper + InterProScan | structural pipeline first, then orthology + domains |
| Eukaryotic proteome | InterProScan (domains+GO+pathways) + eggNOG-mapper | orthogonal evidence, reconcile |
| Metagenome / MAG | eggNOG-mapper `--itype metagenome` (+ KofamScan, dbCAN) | built-in gene calling; KEGG modules |
| Twilight-zone / ORFan (no homolog) | ML (DeepGOPlus) or structure (ESMFold -> Foldseek -> DeepFRI) | only handle on the homology-free fraction; low-confidence leads |
| CAZymes / BGCs / AMR / signal peptides | -> dbCAN / antiSMASH / AMRFinderPlus / SignalP6 | a generic Pfam hit gives no substrate/phenotype/cluster |
| GO enrichment downstream | -> pathway-analysis/go-enrichment (mind IEA circularity) | enrichment on IEA partly tests the pipeline against itself |

## eggNOG-mapper

```bash
download_eggnog_data.py --data_dir db/ -y               # ~44 GB (DIAMOND DB installed by default; -D skips it)
emapper.py -i proteins.faa --itype proteins -m diamond \
    --tax_scope auto --data_dir db/ --cpu 16 -o annot --output_dir out/
```

Three stages: (1) **seed-ortholog search** (DIAMOND/MMseqs2/HMMER) anchors the query - this is a best-hit and is *not* the annotation; (2) **orthology assignment** retrieves the seed's fine-grained orthologs within the chosen taxonomic scope; (3) **functional transfer** pools terms across the *set of orthologs* (which damps single-entry misannotation - this is why eggNOG-mapper beats raw DIAMOND-vs-NR). `--tax_scope` is **the single most consequential parameter**: too broad gathers distant orthologs and over-generalizes function; `auto` lets each seed take its most-informative phylogenetic ceiling. `--itype {proteins,CDS,genome,metagenome}` (genome/metagenome runs Prodigal first). Output `.emapper.annotations` columns include `seed_ortholog`, `eggNOG_OGs`, `COG_category`, `Description`, `Preferred_name`, `GOs`, `EC`, `KEGG_ko`, `PFAMs` (read the actual header; `-` = empty).

## InterProScan

```bash
interproscan.sh -i proteins.faa -f TSV,GFF3 -goterms -pa -cpu 16
```

Runs member-database scanners (Pfam, PANTHER, NCBIfam, SUPERFAMILY, CDD, SMART, Gene3D, Hamap, PROSITE, ...) and **integrates overlapping signatures into InterPro entries** (stable IPRxxxxxx, with a type: Family/Domain/Repeat/Site/Homologous Superfamily). **Report at the InterPro-entry level** - it is the consensus that survives one member DB being wrong. `-goterms` adds the interpro2go mapping (these GO are IEA/electronic); `-pa` maps Reactome/MetaCyc. By default it queries the EBI precalculated lookup service (fast, MD5-keyed); `-dp` forces local compute (novel/confidential sequences, reproducibility). Java 11+ and a tens-of-GB data bundle required; for millions of proteins, chunk the FASTA into array jobs.

## Reconciling Multi-Tool Output with Python

**Goal:** Merge eggNOG and InterProScan per protein while preserving provenance, so a curated name is never silently overwritten by a generic domain.

**Approach:** Parse each tool's table, keep source namespaces separate, union GO with source tags, and prefer the orthology `Preferred_name`/`Description` for the human-readable product.

```python
import pandas as pd

def parse_eggnog(path):
    df = pd.read_csv(path, sep='\t', comment='#', header=None)
    cols = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl',
            'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko']
    df.columns = (cols + [f'c{i}' for i in range(len(df.columns) - len(cols))])[:len(df.columns)]
    return df

def best_product_name(row):
    name = row.get('Preferred_name', '-')
    return name if name not in ('-', '', None) else 'hypothetical protein'   # honest default, not a forced guess
```

Use AGAT (`agat_sp_manage_functional_annotation.pl`) to graft BLAST/InterProScan results onto a GFF3 (it handles the spec edge cases). For GO deliverables use GAF (carries the evidence code); keep each tool in its own `Dbxref` namespace.

## Ontology Rigor and the IEA Circularity

- **GO MF vs BP transfer with different reliability.** Molecular Function ("DNA helicase activity") is local/chemical and transfers with the fold; Biological Process ("DNA replication") is systemic context and does *not* transfer reliably - CAFA confirmed BLAST beats naive baselines for MF but not BP (Radivojac 2013 *Nat Methods* 10:221). Weight MF over BP when evidence is limited; never let a transferred BP term drive a conclusion alone.
- **Essentially all genome-derived GO is IEA** (Inferred from Electronic Annotation, never curator-reviewed). Running GO enrichment on IEA against an IEA background **partly tests the pipeline against itself** - if interpro2go maps a common domain to a term, every genome with that domain looks "enriched." Compounded by annotation bias (58% of human GO covers 16% of genes; Haynes 2018 *Sci Rep* 8:1362) and True-Path-Rule inflation of shallow terms. Pin versions, match background to foreground pipeline, prefer non-IEA where it exists, and state that the result is annotation-derived.
- **EC numbers are not stable.** Deleted/transferred numbers are tombstones (never reused); a partial EC `1.1.1.-` is a valid statement of ignorance (the EC equivalent of "hypothetical"). Demand orthology or a curated rule before asserting a full four-level EC.
- **KEGG bulk access is paywalled.** Free routes: KofamScan/KofamKOALA (local HMMs + adaptive per-KO thresholds; an `*` marks above-threshold hits) or eggNOG's `KEGG_ko`. A "complete module" is a reconstruction (a gap can be non-orthologous gene displacement; a filled step can be a paralog doing something else), not proof of flux.

## Per-Method Failure Modes

### Best-hit-as-ortholog
**Trigger:** transferring a specific function from one DIAMOND/BLAST top hit (esp. TrEMBL/NR). **Mechanism:** best-hit != ortholog; the bulk-DB hit is likely itself an auto-annotation. **Symptom:** confident specific names with no provenance. **Fix:** orthology consensus + Swiss-Prot donors.

### Over-specific transfer
**Trigger:** copying the exact substrate/EC of a characterized homolog onto a distant relative. **Mechanism:** mechanistically-diverse superfamilies share fold, not substrate. **Symptom:** a "muconate cycloisomerase" that does something else. **Fix:** demote to superfamily / partial EC as identity and coverage fall.

### Wrong eggNOG tax_scope
**Trigger:** leaving scope too broad/narrow or unpinned. **Mechanism:** distant orthologs over-generalize, or no informative orthologs. **Symptom:** vague or missing function. **Fix:** `auto`, or pin the known clade.

### Circular GO enrichment
**Trigger:** enriching IEA annotations against a mismatched background. **Mechanism:** measures the mapping table and study popularity, not biology. **Symptom:** "enriched" for whatever well-studied genes are annotated for. **Fix:** non-IEA where possible; matched background; pin versions; caveat the result.

### Reading specialized function from a generic hit
**Trigger:** inferring CAZyme substrate / AMR phenotype / BGC product from a plain Pfam domain. **Mechanism:** the substrate/phenotype/cluster signal is not in a generic domain. **Symptom:** wrong substrate or phenotype call. **Fix:** route to dbCAN / AMRFinderPlus / antiSMASH.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Reason in bits-per-residue, not raw e-value | alignment statistics | e-value scales with DB size (a database-size artifact); bits/residue is density |
| Bidirectional coverage ≥50-70% query and subject | transfer practice | one-domain coverage justifies only a domain-level claim |
| ~40% identity over full length (well-behaved families only) | soft floor | no safe identity in mechanistically-diverse superfamilies; demote specificity instead |
| Named fraction "too high for the taxon" (>90% on a novel isolate) | over-annotation smell test | loose thresholds manufacturing names; expect 20-50% hypothetical |
| eggNOG `--tax_scope auto` | eggNOG-mapper | per-seed informative ceiling |
| KofamScan adaptive per-KO threshold (`*`) | Aramaki 2020 | a single global e-value misfires across KO families |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low annotation rate | fragmented ORFs / narrow scope | check protein quality; `--tax_scope auto`; run both tools and merge |
| Specific name on a distant homolog | over-specific transfer | demote to superfamily / partial EC; record identity |
| eggNOG DB errors | DB/version mismatch | re-download; pin emapper 2.1.15 |
| InterProScan memory/time | full proteome at once | chunk FASTA; keep lookup service on; drop PANTHER/Gene3D if not needed |
| Enrichment "too clean" | IEA circularity / study bias | matched background; pin GO release; caveat |
| Multidomain protein mislabeled | named by first/best domain | report all domains with coordinates |

## References

- Cantalapiedra CP, et al. 2021. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol* 38:5825-5829.
- Huerta-Cepas J, et al. 2019. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource. *Nucleic Acids Res* 47:D309-D314.
- Jones P, et al. 2014. InterProScan 5: genome-scale protein function classification. *Bioinformatics* 30:1236-1240.
- Blum M, et al. 2025. InterPro: the protein sequence classification resource in 2025. *Nucleic Acids Res* 53:D444-D456.
- Schnoes AM, et al. 2009. Annotation error in public databases: misannotation of molecular function in enzyme superfamilies. *PLoS Comput Biol* 5:e1000605.
- Gilks WR, et al. 2002. Modeling the percolation of annotation errors in a database of protein sequences. *Bioinformatics* 18:1641-1649.
- Aramaki T, et al. 2020. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics* 36:2251-2252.
- Zheng J, et al. 2023. dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Res* 51:W115-W121.
- Teufel F, et al. 2022. SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nat Biotechnol* 40:1023-1025.
- Feldgarden M, et al. 2021. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep* 11:12728.
- Blin K, et al. 2023. antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures and visualisation. *Nucleic Acids Res* 51:W46-W50.
- Radivojac P, et al. 2013. A large-scale evaluation of computational protein function prediction (CAFA). *Nat Methods* 10:221-227.
- Haynes WA, et al. 2018. Gene annotation bias impedes biomedical research. *Sci Rep* 8:1362.

## Related Skills

- prokaryotic-annotation - Bakta/PGAP product names + locus tags before functional layers
- eukaryotic-gene-prediction - Produces the protein FASTA to annotate
- annotation-qc - Annotation-coverage and hypothetical-fraction sanity
- pathway-analysis/go-enrichment - Enrichment using GO annotations (mind IEA circularity)
- pathway-analysis/kegg-pathways - Pathway mapping with KEGG orthologs
- epidemiological-genomics/amr-surveillance - AMRFinderPlus/CARD for resistance genes and point mutations
