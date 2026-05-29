---
name: bio-genome-annotation-annotation-qc
description: Assesses the quality and completeness of a genome annotation with BUSCO (conserved single-copy ortholog recovery), OMArk (proteome completeness, consistency, and contamination), CheckM2 (prokaryotic completeness/contamination), and a gene-set sanity panel (gene count, mono-exonic fraction, protein-length distribution, mRNA:gene ratio, coding density). Covers the assembly-BUSCO-vs-proteome-BUSCO diagnostic, what BUSCO-Duplicated really means, why gene count is a vanity metric, and the QC of transferred annotations. Use when judging whether an annotation is good enough to publish or submit, diagnosing a suspect annotation, or comparing annotation completeness across pipelines.
tool_type: cli
primary_tool: BUSCO
---

## Version Compatibility

Reference examples tested with: BUSCO 5.5+, OMArk 0.3+, CheckM2 1.0+, compleasm 0.2.6+, gffutils 0.12+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

BUSCO results depend on the **lineage dataset** (e.g. `vertebrata_odb10` vs the shallow `eukaryota_odb10`) and the OrthoDB version - record them; a 99% on a shallow ~255-gene set and a 99% on a deep ~5,500-gene clade set are very different claims. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Annotation QC

**"Is my genome annotation any good?"** -> Measure conserved-gene completeness, proteome consistency/contamination, and gene-set sanity, and decide whether the limiting factor is the assembly or the predictor.
- CLI: `busco -i proteins.faa -m proteins -l <lineage>_odb10` (proteome) and `busco -i genome.fa -m genome -l <lineage>_odb10` (assembly), `omark`, `checkm2 predict` (prokaryotes)

## The Single Most Important Modern Insight -- BUSCO Measures the Easy 10%; the Diagnostic Is Genome-vs-Proteome

BUSCO completeness measures recovery of the ~1,000-5,500 most conserved single-copy orthologs - the most conserved, highly expressed, intron-stable genes any half-broken pipeline will find. A 98%-complete BUSCO confirms the housekeeping core is intact; it says **nothing** about whether the other ~20,000 models are chimeric, fragmented, frame-shifted, fused, or hallucinated from TEs. BUSCO is a smoke detector in one room of a burning house. Three load-bearing moves:

1. **Run assembly-BUSCO vs proteome-BUSCO on the same assembly - the diagnostic almost nobody runs.** `-m genome` does its own gene-finding; `-m proteins` scores the delivered proteome. If assembly is 98% and proteome is 85%, the genes are physically present and the *predictor* missed them - a training/evidence/masking problem, fixable without touching the assembly. If both are low, the genes aren't in the assembly - stop annotating and fix the assembly. This single fork redirects more wasted effort than any other check. (Always run BUSCO on the *delivered proteome*, `-m proteins`, not genome mode reported as if it described the annotation.)
2. **Gene count is a vanity metric.** The *same* assembly annotated by two pipelines routinely differs 20-40% in gene count, and a higher count is as likely to mean spurious TE ORFs and split models as more real genes. The diagnostic signal is in the **mono-exonic fraction, protein-length distribution, mean exons/gene, and mRNA:gene ratio** - not the headline count.
3. **Complement BUSCO with OMArk.** BUSCO's single-copy-ortholog lens is blind to over-prediction, chimeras, and contamination; OMArk assesses the *whole* proteome for completeness AND consistency (are genes consistent with their homologs?) and detects contaminant species. Use both.

## Tool Taxonomy

| Tool | Citation | Measures | Domain |
|------|----------|----------|--------|
| BUSCO | Manni 2021 *Mol Biol Evol* | conserved single-copy ortholog recovery (C/D/F/M) | euk + prok + viral; proteome/genome/transcriptome modes |
| OMArk | Nevers 2025 *Nat Biotechnol* | proteome completeness + consistency + contamination | eukaryotic proteomes; catches over-prediction BUSCO misses |
| compleasm | Huang 2023 | faster miniprot-based BUSCO-compatible completeness | euk; genome mode |
| CheckM2 | Chklovski 2023 *Nat Methods* | completeness + contamination (ML, lineage-agnostic) | bacteria/archaea, isolates + MAGs |
| Gene-set sanity panel | (gffutils) | gene count, mono-exonic %, protein length, mRNA:gene, coding density | any; the panel BUSCO can't provide |
| LAI | Ou 2018 *NAR* | LTR-RT assembly resolution (repeat-space contiguity) | LTR-rich genomes; assembly-side QC |

## Decision Tree by Scenario

| Scenario | QC to run | Why |
|----------|-----------|-----|
| Prokaryotic genome/MAG | CheckM2 (completeness/contamination) + coding density + tRNA/rRNA counts | CheckM2 is the field-standard completeness call; BUSCO is conservative on prokaryotes |
| Eukaryotic annotation | BUSCO `-m proteins` AND `-m genome` + OMArk + gene-set sanity panel | the genome-vs-proteome fork + over-prediction/contamination |
| Suspect high gene count | mono-exonic fraction + protein-length + BUSCO-Duplicated | distinguish haplotigs / unmasked TEs / split models |
| High BUSCO-Duplicated | check synteny + clade ploidy | uncollapsed haplotigs (purge) vs real WGD (keep) |
| Transferred annotation | BUSCO on the lifted set vs reference + ORF integrity | quantify silently lost conserved genes |
| Low proteome-BUSCO, high genome-BUSCO | fix evidence/masking/training | the genes are present; the predictor missed them |
| Both BUSCO low | -> genome-assembly/assembly-qc (purge_dups, more data) | the limit is the assembly |

## BUSCO and How to Read It

```bash
busco -i proteins.faa -m proteins -l vertebrata_odb10 -o busco_prot -c 16   # the delivered proteome
busco -i genome.fa    -m genome   -l vertebrata_odb10 -o busco_genome -c 16 # the assembly
```

Reported as `C:[S,D],F,M`: **Complete** (full-length match), **Duplicated** (subset of Complete, found ≥2x), **Fragmented** (partial), **Missing**. Use the **deepest applicable clade dataset**, not the shallow `eukaryota_odb10` (a 99% on a 255-gene set is trivially easy and not comparable to a 99% on a 5,500-gene clade set). **High Duplicated is the most-misread signal** - three causes, opposite responses: uncollapsed haplotigs (genome-wide, heterozygosity-scaled, non-syntenic -> purge_dups *before* annotating), real recent WGD (syntenic, ploidy-consistent, Ks duplication peak -> keep), or split models from a fragmented assembly. A clean haploid annotation runs ~1-3% Duplicated; >5-8% with no known WGD screams purge_dups first. compleasm is a faster miniprot-based alternative that often reports higher completeness than BUSCO's metaeuk path.

## OMArk (Proteome Consistency and Contamination)

```bash
omamer search --db LUCA.h5 --query proteins.faa --out proteins.omamer
omark -f proteins.omamer -d LUCA.h5 -o omark_out
```

OMArk catches what BUSCO's single-copy lens misses: it classifies the *whole* proteome as consistent / inconsistent (a gene whose structure conflicts with its homologs - a chimera or fragment) / unknown, and detects **contaminant species** mixed into the proteome. A high "inconsistent" fraction flags over-prediction or fused/split models even when BUSCO is green.

## CheckM2 (Prokaryotic Completeness/Contamination)

```bash
checkm2 predict --input genome.fna --output-directory checkm2_out --threads 16
```

For bacteria/archaea, CheckM2's lineage-agnostic ML model is the standard completeness/contamination call (handles reduced/novel lineages where marker sets fail). Gate annotation on it: **contamination >5%** mixes two organisms' genes into a chimeric set; **completeness <90% with contamination >5-10%** makes gene count, coding density, and hypothetical fraction uninterpretable - fix the assembly/binning first.

## Gene-Set Sanity Panel with Python

**Goal:** Compute the triage panel that reveals annotation health where gene count and BUSCO cannot.

**Approach:** Load the GFF3 into gffutils; compute the mRNA:gene ratio (1.00 = isoform/UTR-naive), the mono-exonic fraction, and the protein-length distribution; flag clade-anomalous values.

```python
import gffutils

MONOEXONIC_FLAG = 0.30   # >30% single-exon in a vertebrate suggests unmasked TEs/pseudogenes/fragments (calibrate per clade; fungi run higher)

def sanity_panel(gff_file):
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy='merge')
    genes = list(db.features_of_type('gene'))
    mrnas = list(db.features_of_type(['mRNA', 'transcript']))
    exon_counts = [len(list(db.children(tx, featuretype='exon'))) for tx in mrnas]
    mono_frac = sum(1 for e in exon_counts if e == 1) / len(exon_counts) if exon_counts else 0
    mrna_per_gene = len(mrnas) / len(genes) if genes else 0
    if mrna_per_gene <= 1.001:
        print('WARNING: mRNA:gene == 1.00 -- isoform/UTR-naive; AS/3-prime-tag scRNA-seq analyses untrustworthy')
    if mono_frac > MONOEXONIC_FLAG:
        print(f'WARNING: mono-exonic fraction {mono_frac:.1%} high -- check masking/contamination')
    return {'genes': len(genes), 'mrna_per_gene': mrna_per_gene, 'mono_exonic_fraction': mono_frac}
```

Read gene count against the *nearest well-annotated relative* and the species' ploidy, never in isolation (1.5-2x with no WGD = haplotigs/over-prediction; ~0.5x = over-masking/under-training). A healthy protein-length distribution is unimodal near the clade-typical ~300-450 aa; a sub-100-aa spike = spurious/fragmented calls; a fat left tail = partials from a fragmented assembly.

## Per-Method Failure Modes

### Reporting genome-mode BUSCO as the annotation's score
**Trigger:** running BUSCO `-m genome` and citing it as annotation quality. **Mechanism:** genome mode does its own gene-finding, often better than the author's pipeline on the conserved core. **Symptom:** reported BUSCO higher than the real proteome BUSCO. **Fix:** run `-m proteins` on the delivered proteome; report both and compare.

### High Duplicated read as success
**Trigger:** treating high BUSCO-D as good. **Mechanism:** uncollapsed haplotigs vs real WGD vs split models. **Symptom:** inflated gene count. **Fix:** plot duplicated pairs against synteny + clade ploidy; if no WGD and D>5-8%, purge_dups the assembly first.

### Shallow lineage dataset
**Trigger:** `eukaryota_odb10` (~255 genes) instead of the deep clade set. **Mechanism:** a small, conserved set is trivially easy to hit 99%. **Symptom:** misleadingly high completeness. **Fix:** use the deepest applicable clade dataset; record it.

### BUSCO alone, no OMArk
**Trigger:** judging a proteome by BUSCO completeness only. **Mechanism:** the single-copy-ortholog lens misses over-prediction, chimeras, and contamination. **Symptom:** green BUSCO on a proteome full of fused/fragmented or contaminant models. **Fix:** add OMArk for consistency + contamination.

### Annotating before the QC gate (prokaryote)
**Trigger:** Bakta/Prokka before CheckM2. **Mechanism:** contamination mixes organisms; low completeness truncates. **Symptom:** chimeric/inflated gene set, uninterpretable coding density. **Fix:** CheckM2 first; contamination >5% -> decontaminate.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Assembly-BUSCO vs proteome-BUSCO gap | the diagnostic fork | large gap = predictor missed present genes; both low = fix assembly |
| BUSCO-Duplicated ~1-3% (clean haploid); >5-8% no WGD | assembly norm | purge_dups before annotating |
| Use deepest applicable clade `_odb10` | BUSCO guidance | shallow sets trivially hit 99% |
| Mono-exonic ~10-20% (vertebrate); >25-30% flag | clade norm | unmasked TEs/pseudogenes/fragments; fungi legitimately higher |
| Protein length unimodal ~300-450 aa | eukaryote norm | sub-100-aa spike = spurious; fat left tail = partials |
| mRNA:gene ratio == 1.00 | annotation structure | isoform/UTR-naive |
| CheckM2 contamination ≤5%, completeness ≥90% | MIMAG-aligned | above/below -> prokaryotic QC numbers uninterpretable |
| Prokaryotic coding density ~88-90% | bacterial norm | <85% wrong table/fragmentation; >93% over-calling |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| BUSCO high, annotation still bad | BUSCO measures only the conserved core | add OMArk + sanity panel; inspect mono-exonic models |
| Proteome-BUSCO < genome-BUSCO | predictor missed present genes | fix evidence/masking/training, not the assembly |
| High Duplicated | haplotigs vs WGD vs split models | synteny + ploidy check; purge_dups if no WGD |
| 99% completeness looks too good | shallow lineage dataset | rerun with the deep clade set |
| Cross-pipeline gene counts disagree | gene count is not a quality metric | compare sanity panels, not counts |
| CheckM2 high contamination | mixed organisms / poor binning | decontaminate before annotation |

## References

- Manni M, et al. 2021. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. *Mol Biol Evol* 38:4647-4654.
- Simão FA, et al. 2015. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. *Bioinformatics* 31:3210-3212.
- Nevers Y, et al. 2025. Quality assessment of gene repertoire annotations with OMArk. *Nat Biotechnol* 43:124-133.
- Chklovski A, et al. 2023. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 20:1203-1212.
- Huang N, Li H. 2023. compleasm: a faster and more accurate reimplementation of BUSCO. *Bioinformatics* 39:btad595.
- Guan D, et al. 2020. Identifying and removing haplotypic duplication in primary genome assemblies (purge_dups). *Bioinformatics* 36:2896-2898.
- Ou S, Chen J, Jiang N. 2018. Assessing genome assembly quality using the LTR Assembly Index (LAI). *Nucleic Acids Res* 46:e126.

## Related Skills

- eukaryotic-gene-prediction - The annotation whose proteome/gene-set this QC evaluates
- prokaryotic-annotation - CheckM2 gate and coding-density sanity for prokaryotes
- annotation-transfer - BUSCO on a lifted set quantifies silently lost conserved genes
- repeat-annotation - LAI and over-masking diagnosis in the assembly-to-annotation handoff
- genome-assembly/assembly-qc - Assembly-side completeness; purge haplotigs before annotating
