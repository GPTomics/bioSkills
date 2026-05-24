---
name: bio-crispr-screens-in-vivo-screens
description: Designs and analyzes in vivo CRISPR screens in animal tumor models, organoids, and immune-cell adoptive transfers. Covers bottleneck math (250x cells/sgRNA requires ~25M cells implanted; impossible for most syngeneic models, forcing focused libraries), focused library design (Manguso 2017 Nature 547:413 immune screen; Chen 2015 tumor screens), CRISPR-StAR intrinsic-control screening (Uijttewaal 2025 Nat Biotechnol 43:1848), clonal-dynamics-limited detection, tumor-explant DNA recovery, syngeneic vs xenograft vs PDX considerations, and the relationship to downstream MAGeCK / drugZ analysis. Use when designing in vivo CRISPR screens for tumor / immune / metastasis biology, choosing focused vs genome-wide for animal models, addressing bottleneck-induced clonal collapse, picking the syngeneic / xenograft / PDX model, integrating in vivo with in vitro results, or applying CRISPR-StAR for animal experiments.
tool_type: mixed
primary_tool: MAGeCK
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5.9+, MAGeCK-VISPR 0.5.6+, pandas 2.2+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version`
- Reference focused libraries: Manguso 2017, Chen 2015, public Addgene aliquots

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## In Vivo CRISPR Screen Analysis

**"Design or analyze an in vivo CRISPR screen"** -> Account for the dramatic bottleneck during animal implantation and tumor growth; use focused libraries; recover DNA from tumor explants; analyze with bottleneck-adjusted hit calling.

- CLI: `mageck count` + `mageck test` for standard analysis
- Special handling: bottleneck-adjusted coverage thresholds; per-tissue per-animal replicate structure

## The In Vivo Bottleneck Problem

**Why in vivo screens differ from in vitro:**

| Constraint | In vitro | In vivo |
|------------|----------|---------|
| Cells per condition | 10M-100M (unlimited) | Limited by injection volume (1-5M cells typical) |
| Implant -> early tumor cell count | N/A | 10-100x drop typical (~94% library complexity may survive in CD45+ TILs) |
| Late tumor cell count | N/A | Further 5-10x reduction; final ~3.93 sgRNAs/gene |
| Bottleneck per animal | None | Tens of millions of cells fail to engraft |
| Library coverage achievable | 500-1000x | Often 50-100x effective at endpoint |
| sgRNAs survivable | Full library | 80-94% in early tumors; 60-80% in late tumors |

**Math:** A 70,000-sgRNA library at 500x coverage requires 35M cells in pool. Most syngeneic models can implant 1-5M cells. Result: real coverage is 70x at best; effective coverage at endpoint is even lower after bottleneck.

**Solution:** Use focused libraries (500-3,000 genes; ~3,000-15,000 sgRNAs) to maintain reasonable coverage despite the bottleneck.

## Focused Library Design for In Vivo

**Manguso et al 2017 *Nature* 547:413** established the canonical in vivo CRISPR screen methodology with a focused library:

- 2,368 genes (kinases + cell surface proteins + immune factors)
- 4 sgRNAs per gene (~9,500 sgRNAs total)
- Targeted at immune-evasion biology in syngeneic mouse melanoma
- Identified PTPN2, NLRC5, CD47 as immune-evasion genes

**Standard focused-library principles:**

1. **Gene selection:** Define the biology you want to test (e.g., immune evasion, metastasis); restrict library to genes plausibly involved (kinases, surface proteins, regulators)
2. **Library size:** 500-3,000 genes; 3-6 sgRNAs/gene; total 3,000-18,000 sgRNAs
3. **Coverage achievable:** With 5M cells implanted, 100-300x coverage is achievable

**Public focused libraries:**
- Manguso 2017 immune library (Addgene)
- Joung 2017 pooled screen reference
- DepMap focused panels for specific pathways
- Custom: order from Twist via CRISPick or CRISPOR

## CRISPR-StAR (Temporal Activation; Uijttewaal 2025)

**Uijttewaal et al 2025 *Nat Biotechnol* 43(11):1848** (published online Dec 2024) introduced CRISPR-StAR (Staggered Activation Reporter), which uses inducible Cas9 expression to delay gene editing until after cells have engrafted in the animal.

**How it works:**
1. Library is delivered into cells with inducible (Dox-controlled) Cas9
2. Cells are implanted in animal at MOI 0.3
3. Cells engraft into tumor (no editing yet); library complexity preserved
4. Days 5-7 post-implantation: induce Cas9 with doxycycline
5. Editing begins; screen proceeds for additional weeks
6. Tumor harvest, DNA extraction, sequencing

**Quantified gain:** CRISPR-StAR enables genome-scale in vivo screens (vs focused libraries) by generating intrinsic per-clone controls; outperforms conventional in vivo screens in therapy-resistant mouse melanoma models (Uijttewaal 2025).

## Syngeneic vs Xenograft vs PDX

| Model | Immune system | Use case |
|-------|---------------|----------|
| Syngeneic (e.g., B16 melanoma in C57BL/6) | Intact mouse immunity | Tumor-immune interaction; checkpoint biology |
| Xenograft (human cancer line in NSG) | Absent / impaired | Tumor cell-intrinsic biology; drug response in human cells |
| PDX (patient-derived xenograft) | Absent / impaired | Patient-specific biology; therapy testing |
| Humanized mouse | Reconstituted human immunity | Tumor-immune in human context (limited) |
| Organoid in vivo | None (in vitro) | Tumor cell-intrinsic in 3D structure |

**Decision rule:** For immune-targeting drug screens, use syngeneic. For human-cancer cell-intrinsic biology, use xenograft. For patient-specific drug screens, use PDX. Each requires different cell numbers and bottleneck planning.

## Tumor DNA Extraction and Sequencing

**Goal:** Recover sufficient sgRNA-containing DNA from tumor explants for sequencing.

**Approach:** Dissect tumor; lyse with proteinase K; extract genomic DNA; amplify the sgRNA locus by PCR; sequence on MiSeq / NextSeq / NovaSeq.

```bash
# Typical PCR + sequencing parameters for in vivo screens
# Per-tumor DNA: 0.5-5 mg yield from typical syngeneic tumor
# Per-sample sequencing depth: ≥500 reads/sgRNA at endpoint (lower than in vitro 300+)
# Multiple animals per condition (n=5-10) to account for clonal variation

# mageck count for in vivo
mageck count \
    --list-seq library.csv \
    --sample-label Plasmid,Animal1,Animal2,Animal3,Animal4,Animal5 \
    --fastq Plasmid.fq.gz A1.fq.gz A2.fq.gz A3.fq.gz A4.fq.gz A5.fq.gz \
    --norm-method median \
    --output-prefix in_vivo_screen
```

## Hit Calling for In Vivo

**Goal:** Identify per-gene fitness effects despite high inter-animal variability.

**Approach:** Each animal is a "replicate" with high variance due to clonal dynamics. Use MAGeCK MLE with animal-as-batch covariate, or run MAGeCK RRA per animal and meta-analyze.

```bash
# Option A: MAGeCK MLE with batch covariate
cat > in_vivo_design.txt <<EOF
Samples         baseline    tumor   animal_2  animal_3  animal_4  animal_5
Plasmid         1           0       0         0         0         0
Animal1         1           1       0         0         0         0
Animal2         1           1       1         0         0         0
Animal3         1           1       0         1         0         0
Animal4         1           1       0         0         1         0
Animal5         1           1       0         0         0         1
EOF

mageck mle \
    --count-table in_vivo_screen.count.txt \
    --design-matrix in_vivo_design.txt \
    --output-prefix in_vivo_mle
```

**Per-animal RRA + meta-analysis:**

```python
import pandas as pd

# Run mageck test on each animal vs plasmid
# Combine with Stouffer's Z method
from scipy.stats import norm

def meta_analyze_animals(per_animal_results):
    '''per_animal_results: list of MAGeCK gene_summary.txt per animal.'''
    merged = pd.concat([df.assign(animal=i) for i, df in enumerate(per_animal_results)])
    grouped = merged.groupby('id')
    meta = grouped.apply(lambda g: pd.Series({
        'mean_neg_score': g['neg|score'].mean(),
        'stouffer_z': norm.ppf(g['neg|p-value']).sum() / (len(g) ** 0.5),
        'animals_significant': (g['neg|fdr'] < 0.05).sum(),
        'n_animals': len(g)
    }))
    return meta.sort_values('stouffer_z')
```

## Failure Modes

### Clonal dominance from low complexity

**Trigger:** Implanted cells lack sufficient library complexity; a few clones dominate the tumor.
**Mechanism:** Inter-animal stochasticity in cell engraftment creates founder effects.
**Symptom:** Per-animal hit lists vary dramatically; no genes appear across all animals.
**Fix:** Use focused library to maintain coverage; increase animals per condition (n=10+); use CRISPR-StAR to delay bottleneck.

### Tumor DNA extraction yields no sgRNA reads

**Trigger:** Wrong library or library not amplified well from tumor DNA.
**Mechanism:** PCR primers don't match the sgRNA flanking; or insufficient DNA template.
**Symptom:** Low mapping rate (<10%); few sgRNAs detected per tumor.
**Fix:** Verify library plasmid sequence; design primers specific to lentiviral cassette; use 10-100 ng input DNA + 25 PCR cycles.

### In vivo PR-AUC against CEGv2 is poor

**Trigger:** Tumor biology differs from in vitro CEGv2 calibration; not all essentials are essential in animal context.
**Mechanism:** Cells in vivo have different growth conditions (nutrients, hypoxia, immune pressure) than in vitro; CEGv2 calibration assumes in vitro context.
**Symptom:** CEGv2 PR-AUC <0.5 in vivo despite high in vitro PR-AUC.
**Fix:** Use cell-type-and-context-specific essentialome (e.g., a corresponding in vitro screen of the same cell type) as a baseline; in vivo essentialome is biology-dependent.

### Pre-screen Cas9 selection failure

**Trigger:** Cas9-positive cells were not selected before implantation; library has Cas9-negative escapers.
**Mechanism:** Cas9-negative cells carry sgRNA but no editing; persist in tumor without biological perturbation.
**Symptom:** Specific essentiality signals weak; PR-AUC low.
**Fix:** Always select Cas9-positive cells (FACS or selection) before infection; verify by Cas9 IHC or flow.

### Inter-animal variability dominates hit calling

**Trigger:** Limited animals per condition (n=3-5); each has high variance.
**Mechanism:** Per-animal clonal dynamics produce different sgRNA distributions; no consistent signal across few animals.
**Symptom:** MAGeCK p-values inflated; FDR uncalibrated.
**Fix:** Increase animals per condition to 10+; use meta-analysis across animals (Stouffer); validate top hits in arrayed format with n=10 mice each.

### Tumor heterogeneity destroys screen signal

**Trigger:** Spontaneously arising mutations in some tumor regions create non-clonal heterogeneity.
**Mechanism:** Tumor heterogeneity is genuine biology; not all cells in tumor are descendants of original engrafted cells.
**Symptom:** Per-region sequencing shows different sgRNA distributions within same tumor.
**Fix:** Sample multiple tumor regions; or use whole-tumor genomic DNA pooling (averages out heterogeneity).

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Cells per animal | 1-5M typical for syngeneic; 5-10M for xenograft | Tumor model dependent |
| sgRNAs per gene in library (focused) | 4-6 | Standard convention |
| Library size for in vivo focused | 3,000-15,000 sgRNAs | Maintainable coverage |
| Coverage at endpoint | ≥50x, ideally 100-200x | Lower than in vitro 500x |
| Animals per condition | 10+ for hit-calling; 5 minimum | Inter-animal variability |
| Animals per condition for arrayed validation | 10 | Tighter signal needed |
| In vivo CEGv2 PR-AUC | >0.4 (context-dependent) | Lower than in vitro 0.7 |
| Late tumor sgRNA-per-gene | ~3.93 typical | Empirical from literature |
| Days to harvest (tumor) | 12-21 days post-implant | Time for selection to manifest |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No hits | Library complexity collapsed | Use focused library or CRISPR-StAR |
| Per-animal hit lists differ | Clonal dominance | Use focused library; increase animals |
| Low CEGv2 PR-AUC | Context-specific essentialome | Use in vivo-specific reference set |
| Low mapping rate | Wrong sequencing primers | Verify library lentiviral architecture |
| Coverage at endpoint <50x | Implantation bottleneck | Increase cells implanted; focused library |

## References

- Manguso RT et al. 2017. *Nature* 547:413. In vivo CRISPR screen for immune evasion; canonical focused-library design.
- Chen S et al. 2015. *Cell* 160:1246. Original in vivo Cas9 screening methodology.
- Aoki Y et al. 2023. *Cancer Gene Therapy* (doi 10.1038/s41417-023-00664-5). Clonal dynamics limit selection detection in in-vivo CRISPR screens. (Earlier "Sci Adv 9:eadg2451" attribution could not be verified.)
- Uijttewaal ECH et al. 2025. *Nat Biotechnol* 43:1848 (online Dec 2024). CRISPR-StAR intrinsic-control screening for in vivo models.
- Pacini C et al. 2024. *Nat Commun* 15:1230. DepMap screening benchmark including in vivo.

## Related Skills

- crispr-screens/library-design - Focused library design for in vivo
- crispr-screens/mageck-analysis - MAGeCK MLE with animal-as-batch covariate
- crispr-screens/hit-calling - Per-animal meta-analysis strategies
- crispr-screens/screen-qc - In-vivo-specific QC thresholds
- crispr-screens/batch-correction - Animal cohort as batch in MLE
- crispr-screens/combinatorial-screens - In vivo combinatorial screens
- crispr-screens/copy-number-correction - Cancer-line in vivo screens
- pathway-analysis/go-enrichment - Functional analysis of in vivo hits
