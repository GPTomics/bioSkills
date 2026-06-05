---
name: bio-immunoinformatics-neoantigen-prediction
description: Identify tumor neoantigens from somatic variants with pVACtools (pVACseq/pVACfuse/pVACbind/pVACvector/pVACview) for personalized cancer vaccines and checkpoint biomarkers. Encodes the field's hard truth that binding prediction is the easy, near-solved part and single-digit-percent PPV lives downstream — so it centers clonality/CCF, HLA LOH (the silent invalidator), expression, proximal-variant phasing, agretopicity/foreignness quality, and the predicted->presented->immunogenic validation tiers. Use when nominating vaccine targets, ranking neoantigens, or building a tumor-to-candidate pipeline. Binding details in mhc-binding-prediction; ranking in immunogenicity-scoring.
tool_type: mixed
primary_tool: pVACtools
---

## Version Compatibility

Reference examples tested with: Ensembl VEP 111+, pVACtools 4.1+, MHCflurry 2.1+, VAtools 5+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: pVACtools is now at 7.x; positional CLI args are stable across 4.x-7.x but defaults and the supported-algorithm list change between releases — always run `pvacseq run --help` against the installed build. pVACseq requires the Wildtype and Frameshift VEP plugins; the Downstream plugin was replaced by Frameshift in pVACtools 2.0, so 4.x+ pipelines must NOT use Downstream. A local IEDB install (`--iedb-install-directory`) is strongly preferred over the rate-limited public API for patient data.

# Neoantigen Prediction

**"Find neoantigens from my tumor mutations"** -> Translate somatic variants into mutant peptides, predict patient-HLA presentation, and rank by tumor-specific quality for vaccine/biomarker use.
- CLI: `pvacseq run` on a VEP-annotated, expression/readcount-annotated somatic VCF + patient HLA (pVACtools)
- CLI: `pvacfuse` (fusions via AGFusion/Arriba), `pvacbind` (arbitrary peptides), `pvacview` (manual re-tiering)
- Python: VAtools annotation, LOHHLA/CCF integration, aggregate-report parsing

## The Single Most Important Modern Insight -- binding is the easy part; PPV lives downstream

The visible surface of the field — NetMHCpan, MHCflurry, the IC50 column everyone sorts on — is the binding step, and binding is the one step the field has genuinely cracked. The positive predictive value of a binding-only neoantigen pipeline is single-digit percent: of peptides confidently called strong binders, the large majority are never presented, and of those presented, the large majority never elicit a T-cell response (TESLA; Wells 2020). This is structural, not a bad IC50 cutoff — each step of the presentation-and-recognition cascade multiplies a low conditional probability. The corrective: spend the analysis on the filters and features that govern the predicted->presented->immunogenic attrition (clonality/CCF, HLA LOH, expression, agretopicity, foreignness, processing, validation tiers) and treat the choice of binding algorithm as a near-afterthought with sane defaults. TESLA's five features that actually separated immunogenic peptides from binders: HLA binding affinity, source-gene expression ("tumor abundance"), peptide-HLA binding stability, hydrophobicity, and the two recognition features — agretopicity and foreignness.

## The pVACtools Suite

| Sub-tool | Input that defines the peptide | Use case |
|----------|--------------------------------|----------|
| pVACseq | VEP-annotated somatic VCF (SNV + indel/frameshift) | The workhorse: point mutations and frameshifts |
| pVACfuse | AGFusion / Arriba fusion output | Fusion-junction novel-ORF neoantigens |
| pVACbind | a plain peptide FASTA | Score arbitrary peptides (MS hits, splice peptides); no WT/agretopicity |
| pVACvector | chosen epitopes | Order epitopes into a vaccine string, minimizing junctional neo-epitopes |
| pVACview | `*.all_epitopes.aggregated.tsv` | Human-in-the-loop review and re-tiering (the decision step) |

## Upstream Chain (every link can silently poison the output)

| Step | Tool(s) | Failure if skipped/wrong |
|------|---------|--------------------------|
| Somatic calling (T/N) | Mutect2, Strelka2 (consensus) | Germline leak -> false neoantigens; indels matter most (frameshifts) |
| VEP annotation | VEP + Wildtype + Frameshift plugins, `--fasta`, `--tsl`, `--symbol` | Most error-prone step; wrong plugins -> no WT peptide / no frameshift ORF |
| HLA typing (I and II) | OptiType (class I, WES), arcasHLA (RNA), HLA-HD (II) | Wrong allele = confident garbage; type at 4-digit; reconcile DNA vs RNA |
| Expression | kallisto/salmon TPM, `vcf-expression-annotator` | `--expn-val` passes everything if unannotated -> ships unexpressed "neoantigens" |
| Read counts | bam-readcount, `vcf-readcount-annotator` | VAF/coverage filters pass everything if unannotated |
| Phasing | merge somatic+germline, WhatsHap / GATK ReadBackedPhasing | Proximal in-cis variants -> peptides the patient never makes (neoepiscope, Wood 2020) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| SNV + indel neoantigens | pVACseq, all_class_i | The workhorse; frameshifts via Frameshift plugin |
| Gene fusions | pVACfuse (AGFusion/Arriba, with STAR-Fusion read support) | Junction novel ORFs; demand junction read support |
| Proximal germline/somatic variants nearby | pVACseq `--phased-proximal-variants-vcf` | Otherwise the peptide sequence is wrong |
| Need quality features (DAI, foreignness, dissimilarity) | NeoFox / antigen.garnish on pVAC candidates | pVAC tiers; NeoFox computes the ~16 published features |
| Reproducible end-to-end | nextNEOpi (HLA + VEP + pVACseq + NeoFox + LOHHLA) | Wires the whole chain including LOHHLA and purity |
| Final candidate selection | pVACview manual re-tiering | Tiers say WHY a candidate failed; human triage |

## Run VEP, Then pVACseq

**Goal:** Produce the VEP annotation pVACseq actually consumes, then call neoantigens.

**Approach:** Run VEP with the Wildtype + Frameshift plugins and a protein FASTA; annotate expression and read counts with VAtools; supply a phased proximal-variants VCF; then `pvacseq run` with the patient HLA and sane filters.

```bash
pvacseq install_vep_plugin $VEP_PLUGINS          # installs Wildtype + Frameshift
vep --input_file somatic.vcf --output_file somatic.vep.vcf --format vcf --vcf \
    --symbol --terms SO --tsl --hgvs --fasta GRCh38.fa --offline --cache --dir_cache $VEP_CACHE \
    --plugin Frameshift --plugin Wildtype --pick

vcf-expression-annotator somatic.vep.vcf kallisto.tsv custom transcript -s TUMOR \
    --id-column target_id --expression-column tpm -o somatic.vep.expn.vcf

pvacseq run somatic.vep.expn.vcf TUMOR \
    "HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*44:02,HLA-C*07:02,DRB1*01:01" \
    all_class_i pvac_out/ \
    -e1 8,9,10,11 --iedb-install-directory $IEDB \
    --phased-proximal-variants-vcf phased.vcf.gz \
    --normal-vaf 0.02 --tdna-vaf 0.25 --trna-vaf 0.25 --expn-val 1.0 -t 8
```
Key flags: `-e1/-e2` epitope lengths; `-b/--binding-threshold` (default 500 nM); `--percentile-threshold` (recommend 2); `-m/--top-score-metric` median|lowest; `--allele-specific-binding-thresholds` (preferred over flat 500 nM); `--net-chop-method`/`--netmhc-stab` (processing + stability features).

## Compute Agretopicity (DAI) Correctly

**Goal:** Quantify how much more foreign the mutant looks than its wild-type counterpart.

**Approach:** Agretopicity (the fitness-model amplitude; Łuksza 2017) is the WT/MT binding ratio. A high value means the mutant binds while the WT does not — the surface is new to the immune system, so reactive T cells were not deleted in the thymus. The original differential agretopicity index (DAI; Duan 2014) is the difference form; both forms share the traps below. Requires the matched WT peptide (the Wildtype plugin), so pVACbind cannot compute it.

```python
import pandas as pd

def add_agretopicity(df, wt='Median WT IC50 Score', mt='Median MT IC50 Score'):
    '''Agretopicity (amplitude) = IC50_WT / IC50_MT (ratio > 1 = mutant binds better -> favorable).
    Anchor-position mutations inflate DAI without changing the TCR-facing surface, so
    pair DAI with anchor evaluation rather than trusting it alone.'''
    out = df.copy()
    out['agretopicity'] = out[wt] / out[mt]
    out['dai_favorable'] = out['agretopicity'] > 1
    return out
```

## Drop Candidates on Lost HLA Alleles (LOHHLA)

**Goal:** Remove neoantigens predicted to be presented by an HLA allele the tumor has deleted.

**Approach:** HLA LOH is an immune-escape mechanism in ~40% of NSCLC (McGranahan 2017) and is invisible to binding/expression/clonality filters. Run LOHHLA (or a subclonal-sensitive equivalent like DASH) with the HLA type and tumor purity/ploidy, then filter the aggregate report. This step sits outside pVACtools and errors silently if skipped.

```python
def drop_lost_allele_candidates(df, lost_alleles, allele_col='HLA Allele'):
    '''lost_alleles: set of alleles called as LOH-lost by LOHHLA. A peptide assigned
    to a lost allele is not weakly presented - it is not presented at all.'''
    return df[~df[allele_col].isin(set(lost_alleles))].copy()
```

## Per-Method Failure Modes

### HLA LOH silent invalidation
**Trigger:** ranking candidates without running LOHHLA. **Mechanism:** tumor deletes the haplotype that would present its neoantigens; upstream signals all look fine. **Symptom:** beautiful candidates on an absent allele. **Fix:** mandatory separate LOHHLA step; drop lost-allele candidates.

### Subclonal mis-tiering from raw VAF
**Trigger:** using VAF as clonality without purity/CN correction. **Mechanism:** clonality needs cancer cell fraction (CCF = f(VAF, purity, local CN)). **Symptom:** clonal mutation in low-purity sample read as subclonal (and vice versa in amplified regions). **Fix:** estimate purity (ASCAT/Sequenza/PURPLE) and CCF (PyClone) before tiering.

### Unphased proximal variants
**Trigger:** running pVACseq with only the somatic VCF when nearby in-cis variants exist. **Mechanism:** the translated peptide depends on both variants on the haplotype. **Symptom:** predicted/synthesized peptides the tumor never makes. **Fix:** supply `--phased-proximal-variants-vcf` (merge somatic+germline, phase with WhatsHap/GATK).

### Silent filter pass-through
**Trigger:** expression/VAF/coverage filters set but the values never annotated into the VCF. **Mechanism:** the filter passes everything when the field is absent. **Symptom:** unexpressed/low-coverage candidates in the output. **Fix:** annotate with VAtools first; confirm the FORMAT/INFO fields exist.

### Frameshift/fusion over-trust and MS-gap
**Trigger:** treating frameshift/fusion presentation scores like canonical SNV scores. **Mechanism:** EL/MS training is dominated by canonical 8-11mers from point mutations. **Symptom:** narrow-looking CIs on a poorly-supported class; no MS evidence misread as absence. **Fix:** widen confidence on these high-value classes; build a personalized MS search DB before claiming MS absence.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Binding threshold 500 nM (default) | pVACseq `-b` default | Entry gate only; prefer `--allele-specific-binding-thresholds` / %Rank |
| `--normal-vaf` 0.02 | pVACseq default | Germline-leak guard (esp. tumor-only-ish setups) |
| `--tdna-vaf` / `--trna-vaf` 0.25 | pVACseq default | Min tumor DNA/RNA VAF to keep |
| `--expn-val` 1.0 TPM | pVACseq default | Unexpressed mutation is not a neoantigen |
| Coverage normal/tDNA/tRNA 5/10/10 | pVACseq defaults | Below this, VAF/clonality calls are noise |
| Clonal CCF ~1 (clonal >> subclonal) | McGranahan 2016 | Subclonal targets select for resistant majority |
| Agretopicity/DAI > 1 favorable | Łuksza 2017; TESLA | WT binds poorly -> surface not tolerized |
| Validate beyond tier 1 | Wells 2020; Ott/Sahin 2017 | Predicted != presented != immunogenic |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| pVACseq misses frameshifts / no WT peptide | Used Downstream plugin or omitted Wildtype | Install Wildtype + Frameshift (Downstream dropped in pVACtools 2.0) |
| Everything passes the expression filter | TPM never annotated | `vcf-expression-annotator` before run |
| Non-overlapping neoantigen lists across labs | Different HLA typers/resolution | Type at 4-digit; reconcile DNA vs RNA; WES preferred |
| Candidates on a deleted allele | LOHHLA skipped | Run LOHHLA; drop lost-allele candidates |
| Wrong mutant peptide sequence | Proximal variants unphased | `--phased-proximal-variants-vcf` |
| Subclonal target promoted | Ranked by VAF/IC50, no CCF | Estimate purity + CCF; respect the Subclonal tier |

## References

- Wells DK, van Buuren MM, Dang KK, et al. 2020. Key parameters of tumor epitope immunogenicity revealed through a consortium approach improve neoantigen prediction (TESLA). *Cell* 183(3):818-834.
- Hundal J, Kiwala S, McMichael J, et al. 2020. pVACtools: a computational toolkit to identify and visualize cancer neoantigens. *Cancer Immunology Research* 8(3):409-420.
- McGranahan N, Furness AJS, Rosenthal R, et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. *Science* 351(6280):1463-1469.
- McGranahan N, Rosenthal R, Hiley CT, et al. 2017. Allele-specific HLA loss and immune escape in lung cancer evolution (LOHHLA). *Cell* 171(6):1259-1271.
- Łuksza M, Riaz N, Makarov V, et al. 2017. A neoantigen fitness model predicts tumour response to checkpoint blockade immunotherapy. *Nature* 551:517-520.
- Balachandran VP, Łuksza M, Zhao JN, et al. 2017. Identification of unique neoantigen qualities in long-term survivors of pancreatic cancer. *Nature* 551:512-516.
- Richman LP, Vonderheide RH, Rech AJ. 2019. Neoantigen dissimilarity to the self-proteome predicts immunogenicity and response to immune checkpoint blockade. *Cell Systems* 9(4):375-382.
- Wood MA, Nguyen A, Struck AJ, et al. 2020. neoepiscope improves neoepitope prediction with multivariant phasing. *Bioinformatics* 36(3):713-720.
- Lang F, Riesgo-Ferreiro P, Löwer M, Sahin U, Schrörs B. 2021. NeoFox: annotating neoantigen candidates with neoantigen features. *Bioinformatics* 37(22):4246-4247.
- Ott PA, Hu Z, Keskin DB, et al. 2017. An immunogenic personal neoantigen vaccine for patients with melanoma. *Nature* 547:217-221.

## Related Skills

- immunoinformatics/mhc-binding-prediction - the binding step (the solved, low-leverage part); EL abundance bias bites here
- immunoinformatics/mhc-class-ii-prediction - class II neoantigens for CD4 help (compounded uncertainty)
- immunoinformatics/immunogenicity-scoring - quality ranking (DAI, foreignness, dissimilarity) of the candidate list
- clinical-databases/hla-typing - the genotype substrate; wrong calls poison everything
- clinical-databases/somatic-signatures - clonal neoantigen burden predicts ICI response (McGranahan 2016)
- variant-calling/variant-calling - upstream somatic SNV/indel calls
- workflows/neoantigen-pipeline - the end-to-end orchestration
