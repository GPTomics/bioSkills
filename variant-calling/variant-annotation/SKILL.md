---
name: bio-variant-annotation
description: Annotates VCF variants with functional consequences, population frequencies, and pathogenicity scores using bcftools annotate/csq, Ensembl VEP, SnpEff, and ANNOVAR. Use when deciding which annotation engine and version to pin, which transcript set to report on (RefSeq vs Ensembl vs MANE Select/Plus Clinical, and why VEP --pick is dangerous clinically), how to reconcile HGVS 3'-shifting with VCF left-alignment, which consequence plus NMD status governs PVS1 eligibility, which single calibrated predictor to use for PP3/BP4 (REVEL, AlphaMissense, CADD, SpliceAI deltas), or how to read gnomAD v2/v3/v4 grpmax filtering allele frequency instead of one global AF cutoff. Not for ACMG combining rules or final classification (see variant-calling/clinical-interpretation).
tool_type: mixed
primary_tool: VEP
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, VEP 110+, SnpEff 5.2+, ANNOVAR 2020Jun07+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Note: SnpEff and SnpSift use single-dash `-version`, not `--version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: gnomAD v2 is GRCh37; v3/v4 are GRCh38. MANE transcripts exist only on GRCh38. Confirm the build of every annotation source matches the VCF before annotating.

# Variant Annotation

**"Annotate my variants with functional and clinical information"** -> Map each variant onto a transcript model, classify its coding/splice consequence, and attach population frequency and pathogenicity evidence.
- CLI: `vep` (Ensembl), `snpEff`/`SnpSift`, `table_annovar.pl` (ANNOVAR), `bcftools csq`/`annotate`
- Python: `cyvcf2` to parse VEP CSQ / SnpEff ANN strings; `bcftools +split-vep` to flatten them

## The governing principle: annotation is not deterministic

A variant's consequence is not a property of the variant. It is a property of the tuple (variant, transcript model, engine, engine version, parameter set). Change any element and the reported consequence, the HGVS string, and downstream the ACMG PVS1 eligibility can all change. The single most damaging naive belief in clinical genomics is that a VCF line has one true annotation. It does not. On any discordant result the first question is never "what does the variant do" but "which transcript and which tool+version produced that call."

The decision is therefore NOT to find the "right" tool. It is to PIN every axis and record it on the report: genome build, transcript set (prefer MANE Select on GRCh38), engine + version, predictor + version, gnomAD version. Reproducibility comes from pinning, not from picking. Never compare a variant annotated on RefSeq to one annotated on Ensembl.

Three independent axes of non-determinism: the transcript SET (RefSeq vs Ensembl/GENCODE vs MANE), the transcript-SELECTION heuristic (canonical, worst-consequence, `--pick`, MANE Select), and the ENGINE itself (VEP/SnpEff/ANNOVAR encode different splice-region widths, up/downstream windows, HGVS-shifting rules, and consequence-severity orderings). Discordance concentrates in indels, splice-region, and multi-transcript genes; loss-of-function calls that drive PVS1 are among the least concordant across engines (McLaren 2016 *Genome Biol* 17:122; Cingolani 2012 *Fly* 6(2):80-92; Wang 2010 *Nucleic Acids Res* 38(16):e164).

## Normalize before annotation, and never hand-derive HGVS from POS

Normalization is mandatory: the same variant represented differently produces different annotations.

```bash
# -m-any splits multiallelic records to biallelic so each ALT gets its own annotation
bcftools norm -f reference.fa -m-any input.vcf.gz -Oz -o normalized.vcf.gz
```

The HGVS 3'-rule vs VCF left-align clash is a genuine, still-live trap. VCF normalization requires indels **left-aligned** (most 5' on the forward genomic strand; Tan 2015 *Bioinformatics* 31(13):2202-2204). HGVS mandates the **opposite**: the 3'-rule places an indel in a repeat at the most 3' position with respect to the **transcript**. For a plus-strand gene, transcript-3' is the rightmost genomic position, the opposite end from VCF left-alignment; for a minus-strand gene the two can coincide by accident of strand. Net effect: a correctly left-aligned VCF POS and a correct HGVS `c.` string for the same indel can point to different repeat units.

Rules:
- Left-align + normalize the VCF BEFORE annotation (idempotent representation for matching/merging).
- Rely on the engine's HGVS generator to apply the transcript 3'-shift (VEP does; verify SnpEff/ANNOVAR per version). Never hand-derive an HGVS string from POS.
- When matching a patient variant to a ClinVar/literature HGVS, normalize BOTH to the same representation first. `dup` vs `ins` describe the same event but do not string-match. See variant-calling/variant-normalization.
- Adjacent SNVs in one codon (an MNV split by the caller) can be individually benign but jointly change the amino acid; per-SNV annotation silently misannotates. Use phase-aware annotation.

## Choosing the transcript set

| Set | What it is | When to report on it |
|-----|-----------|----------------------|
| MANE Select | One transcript/gene, byte-identical `NM_`/`ENST` on GRCh38 (Morales 2022 *Nature* 604:310-315) | Default for clinical reporting -- gives a single cross-database-stable c./p. |
| MANE Plus Clinical | Extra isoforms for genes where MANE Select misses known pathogenic variants | Add alongside MANE Select where assigned |
| RefSeq (`NM_`/`NR_`) | NCBI-curated; most ClinVar/HGMD/literature `c.` use these | Legacy clinical pins; may carry sequence absent from the reference |
| Ensembl/GENCODE (`ENST`) | Comprehensive, genome-aligned, more transcripts/gene | Research; NOT 1:1 with RefSeq -- c. positions differ |
| Worst-consequence across all | Most severe over every overlapping transcript | Discovery only -- inflates severity, manufactures false PVS1 |

Decision: report on MANE Select (plus MANE Plus Clinical where assigned), not worst-consequence. Worst-consequence reports a canonical splice change in a minor non-expressed isoform as "splice" even when MANE Select is intronic, manufacturing false PVS1 candidates; restricting to MANE Select alone can miss a variant that only hits a MANE Plus Clinical isoform -- which is exactly why that tier exists. The Ensembl "canonical" transcript is frequently NOT the MANE Select one, so migrating a pipeline to MANE changes some reported c. coordinates (expected, not erroneous). MANE is GRCh38-only; GRCh37 pipelines must lift over or maintain their own per-gene transcript pins.

**Why `--pick` is dangerous for clinical use.** VEP by default reports all consequences for all overlapping transcripts. `--pick` collapses to one block per variant using an ordered heuristic whose defaults are canonical status, biotype, consequence rank, then transcript length, then finally accession order. The late tiebreakers are not clinically motivated: when transcripts tie, `--pick` can let transcript length or alphanumeric `ENST` order decide which consequence a real patient gets, can pick per-variant (so variant A and variant B in one gene land on different transcripts, destroying coordinate consistency), and can silently hide a PVS1-eligible consequence behind a benign one. Defensible configurations: pin MANE Select (with Plus Clinical), or pin the lab's validated per-gene list. If `--pick`-style collapse is used at all, constrain it so MANE leads and length/accession never decide:

```bash
# --pick_order forces MANE first; length/accession can no longer choose the reported transcript
vep -i norm.vcf --vcf --cache --offline --assembly GRCh38 \
    --mane_select --pick --pick_order mane_select,canonical,biotype,rank -o out.vcf
```

## Consequence, impact, and NMD (the PVS1 hinge)

Anchor the vocabulary to Sequence Ontology (SO) terms (`missense_variant`, `stop_gained`, `frameshift_variant`, `splice_donor_variant`, `splice_acceptor_variant`, `start_lost`, `stop_lost`, `inframe_deletion`, `synonymous_variant`, ...), which VEP emits natively. SnpEff/ANNOVAR map to mostly-equivalent terms but differ on splice-region width and finer intronic terms (VEP adds `splice_donor_5th_base_variant`, `splice_polypyrimidine_tract_variant`), so term-level matching across engines is impossible.

Impact buckets are NOT evidence. SnpEff `HIGH/MODERATE/LOW/MODIFIER` and ANNOVAR `exonic;splicing` groupings are triage conveniences. `HIGH` lumps `stop_gained`, `frameshift`, and canonical splice together, but whether any of these earns PVS1 depends on NMD, exon location, and the gene's LOF mechanism -- none of which the bucket knows. Treating "HIGH impact" as "PVS1 met" is a classic error.

**The NMD 50-nt / last-exon rule** governs PVS1 strength (Abou Tayoun 2018 *Hum Mutat* 39(11):1517-1524). A premature termination codon (PTC) more than ~50-55 nt upstream of the last exon-exon junction triggers nonsense-mediated decay -> no protein -> strong LOF. A PTC in the **last exon**, or within ~50 nt of the final junction (3' end of the penultimate exon), **escapes NMD** -- the truncated protein IS made, so full-strength PVS1 on the "NMD -> no protein" logic is unjustified (downgrade to PVS1_Strong/Moderate/Supporting depending on how much protein / which domains are lost). Single-exon genes have no junctions, so a nonsense variant escapes NMD by construction; PTCs near the start can reinitiate downstream. PVS1 also requires that LOF is the established disease mechanism -- a null allele in a gain-of-function/dominant-negative gene must NOT fire PVS1. The ACMG combining lives in variant-calling/clinical-interpretation; this skill supplies the consequence + NMD inputs it needs.

## Annotation engines: run commands

Fix the transcript set and tool+version in the pipeline and record them; that, not tool choice, is what makes annotation reproducible.

### Ensembl VEP

**Goal:** Annotate consequence, HGVS, impact, frequencies, and plugin predictions against Ensembl/MANE.

**Approach:** Run offline against the cache with explicit assembly and transcript selection; add predictor plugins rather than relying on the built-in SIFT/PolyPhen.

```bash
# --everything enables --hgvs --symbol --canonical --af --af_gnomade --af_gnomadg --sift b
# --polyphen b --pubmed etc. Prefer --mane_select over the enabled --canonical for reporting.
vep -i norm.vcf.gz -o out.vcf --vcf --cache --offline \
    --species homo_sapiens --assembly GRCh38 --everything --mane_select --fork 4
```

```bash
# Calibrated predictors as plugins (one calibrated tool per evidence type -- see below)
vep -i norm.vcf.gz -o out.vcf --vcf --cache --offline \
    --plugin dbNSFP,dbNSFP4.3a.gz,REVEL_score,CADD_phred \
    --plugin AlphaMissense,file=AlphaMissense_hg38.tsv.gz \
    --plugin SpliceAI,snv=spliceai_snv.vcf.gz,indel=spliceai_indel.vcf.gz
```

### SnpEff / SnpSift

**Goal:** Fast batch effect annotation plus database cross-referencing.

**Approach:** `snpEff ann` against a prebuilt genome database, then chain `SnpSift annotate`/`filter`.

```bash
# Human GRCh38 DB expands to 3-4 GB in memory; give the JVM >= 8 GB or it OOMs/thrashes
snpEff -Xmx8g ann GRCh38.105 norm.vcf > out.vcf
snpEff -Xmx8g ann GRCh38.105 norm.vcf | SnpSift annotate clinvar.vcf.gz > annotated.vcf
```

### ANNOVAR

**Goal:** Table-driven gene/frequency/pathogenicity annotation.

**Approach:** `table_annovar.pl` with paired `-protocol`/`-operation` lists (g=gene, f=filter, r=region).

```bash
table_annovar.pl norm.vcf humandb/ -buildver hg38 -out annotated -remove \
    -protocol refGene,gnomad30_genome,clinvar_20230416,dbnsfp42a \
    -operation g,f,f,f -nastring . -vcfinput
```

### bcftools csq / annotate

**Goal:** Lightweight consequence prediction (csq) and database field transfer (annotate) without a full engine.

**Approach:** `csq` maps variants to a GFF3 and emits a `BCSQ` field; `annotate -c` copies ID/INFO columns from a position-matched source.

```bash
bcftools csq -f reference.fa -g genes.gff3.gz norm.vcf.gz -Oz -o csq.vcf.gz   # adds BCSQ
bcftools annotate -a dbsnp.vcf.gz -c ID norm.vcf.gz -Oz -o rsid.vcf.gz        # copy rsIDs
```

See usage-guide.md for BED/TAB annotation, field removal, `--set-id`, chromosome renaming, and database download recipes.

## Pathogenicity predictors: one calibrated tool, not a stack

Two structural problems pervade this literature. (1) Circularity: most predictors train on ClinVar/HGMD labels, so benchmarking or ACMG-calibrating on those same databases is partly self-referential; a headline "AUC 0.9x" is optimistic on truly novel variants. (2) Ensembles ingest each other: REVEL is a random forest over 13 component scores including SIFT and PolyPhen, so "REVEL agrees with PolyPhen" is not independent corroboration -- PolyPhen is inside REVEL. Independence between evidence lines is the load-bearing assumption of the ACMG points system; stacking correlated predictors silently over-calls pathogenic.

Therefore: use exactly ONE calibrated predictor per evidence type (missense; splicing), applied at its calibrated strength.

| Predictor | Scope | Use for PP3/BP4 |
|-----------|-------|-----------------|
| REVEL (Ioannidis 2016 *AJHG* 99(4):877-885) | rare missense | Best-calibrated single missense tool; calibrated thresholds below |
| AlphaMissense (Cheng 2023 *Science* 381(6664):eadg7492) | missense, proteome-wide | Not trained on ClinVar labels (uses population frequency + structure); use the CURRENT ClinGen SVI calibrated thresholds, not the developer class cutoffs |
| CADD (Kircher 2014 *Nat Genet* 46:310-315) | all variant types | Genome-wide/non-coding ranking, NOT missense PP3 -- see caveat |
| SIFT / PolyPhen-2 (Ng 2003 *NAR* 31:3812; Adzhubei 2010 *Nat Methods* 7(4):248) | missense | Do not use as standalone evidence -- see below |
| SpliceAI (Jaganathan 2019 *Cell* 176(3):535-548) | splice-altering | The splicing predictor; delta-score interpretation below |

**SIFT/PolyPhen alone are near-worthless now.** In the ClinGen SVI calibration (Pejaver 2022 *Am J Hum Genet* 109(12):2163-2177) neither reached even Supporting strength for PP3; both call a large fraction of all missense "damaging" (low positive predictive value on rare variants); and both are components of REVEL, so quoting them alongside it double-counts. Legacy pipelines surfacing "SIFT: deleterious, PolyPhen: probably damaging" prominently are decorative, not evidentiary.

**CADD >= 20 is not "pathogenic."** In Pejaver 2022 raw CADD did not reach Supporting for PP3, and the developer-recommended CADD >= 20 mapped to Moderate evidence for **benign** -- an inversion of how CADD 20 is casually used. Reserve CADD for its intended non-coding/genome-wide ranking.

**Calibrated REVEL thresholds (Pejaver 2022).** PP3_Supporting >= 0.644 and BP4_Supporting <= 0.290 are the well-reproduced values. The Moderate/Strong REVEL cutoffs (commonly quoted as PP3_Moderate >= 0.773, PP3_Strong >= 0.932; BP4_Moderate <= 0.183, BP4_Strong <= 0.016) come from the supplementary tables and are not uniformly reproduced -- verify against the Pejaver 2022 supplement / current ClinGen SVI recommendation table before hard-coding, rather than treating them as fixed. PP3 and BP4 are mutually exclusive by construction; only tools reaching >= Strong in the calibration qualify.

**SpliceAI delta scores.** Per variant, SpliceAI emits four deltas (acceptor gain/loss, donor gain/loss), each 0-1; the max is the headline. Developer-recommended interpretation: 0.2 high recall, 0.5 recommended, 0.8 high precision. Caveats: (i) know whether the pipeline uses the masked or raw model; (ii) the default scoring window is +/-50 bp -- deep-intronic/pseudoexon effects need a widened window (up to +/-10 kb) or are missed; (iii) a delta is a prediction, and converting it to PS3/PP3 strength needs the ClinGen splicing calibration, not the raw cutoffs; (iv) it does not report the RESULT (exon skip vs intron retention), which is what determines PVS1. Pangolin (Zeng 2022 *Genome Biol* 23:103) is an emerging tissue-aware alternative -- check current ClinGen splicing guidance. AlphaMissense and the splicing calibrations are still-evolving; verify the current ClinGen SVI approved-tool list before standardizing on one.

## Population frequency: grpmax filtering AF, not a global cutoff

gnomAD version + build is itself a trap. v2.1.1 is 141,456 individuals (125,748 exomes + 15,708 genomes) on **GRCh37** (Karczewski 2020 *Nature* 581(7809):434-443); v3 is genomes-only on GRCh38; v4 aggregates ~730k exomes + ~76k genomes on **GRCh38** (release totals per the gnomAD v4 release notes; the v4 genome constraint map is Chen 2024 *Nature* 625(7993):92-100). Comparing AF across versions requires liftover of the VARIANT (not just the coordinate), which can mis-map indels/segdups. "Absent" can mean "not callable here," not "not present in humans" -- always check site coverage/callability and the PASS/`AS_FilterStatus` flags, not just raw AF.

A single global AF cutoff ("AF > 1% -> benign") is wrong in both directions. The maximum credible population AF for a truly pathogenic allele is per-disease -- it depends on prevalence, allelic and genetic heterogeneity, inheritance, and penetrance (Whiffin 2017 *Genet Med* 19(10):1151-1158). Use the **filtering allele frequency (FAF)**: the lower bound of the 95% CI of the **grpmax** AF (the highest AF among genetic-ancestry groups, formerly "popmax"; v4 fields `grpmax`, `AF_grpmax`, `fafmax_faf95_max`, and the joint exome+genome VCF tag `fafmax_faf95_max_joint`). Global AF dilutes a variant common in one ancestry across the whole cohort; grpmax exposes it. Apply BA1/BS1 when FAF exceeds the disease's maximum credible AF -- too-lenient a global line benignizes nothing for ultra-rare high-penetrance disease, and too-strict a global line wrongly benignizes founder alleles that reach several percent in one ancestry.

**"In gnomAD therefore benign" is a fallacy.** Documented exceptions: recessive disease (healthy carriers -> pathogenic alleles present at carrier frequency, e.g. CFTR p.Phe508del); late-onset/reduced-penetrance disease (gnomAD adults can be pre-symptomatic carriers of adult-onset cancer/cardiomyopathy/neurodegeneration alleles, e.g. BRCA/Lynch); somatic/clonal-hematopoiesis contamination (low-AF calls in DNMT3A/TET2 can be somatic, not germline). Ancestry sampling is uneven, so "absent" is much weaker evidence for an under-represented ancestry than for a well-sampled one. The full BA1/BS1/PM2 combining lives in variant-calling/clinical-interpretation.

## Python: parse an annotated VCF

**Goal:** Flatten VEP CSQ (or SnpEff ANN) transcript blocks into per-transcript dicts for filtering.

**Approach:** Read the CSQ format from the header, then split each record's CSQ on commas (transcripts) and pipes (fields).

```python
from cyvcf2 import VCF

vcf = VCF('vep_output.vcf')
csq_fields = None
for h in vcf.header_iter():
    if h['HeaderType'] == 'INFO' and h['ID'] == 'CSQ':
        csq_fields = h['Description'].split('Format: ')[1].rstrip('"').split('|')
        break

for variant in vcf:
    csq = variant.INFO.get('CSQ')
    if not csq:
        continue
    for block in csq.split(','):
        ann = dict(zip(csq_fields, block.split('|')))
        # MANE_SELECT is populated only for the MANE transcript; prefer it over worst-consequence
        if ann.get('MANE_SELECT') and ann.get('IMPACT') in ('HIGH', 'MODERATE'):
            print(variant.CHROM, variant.POS, ann['SYMBOL'], ann['Consequence'])
```

`bcftools +split-vep -f '%CHROM\t%POS\t%SYMBOL\t%Consequence\n' -s worst out.vcf.gz` does the same at the CLI; `-s worst` and `-p` control which block(s) surface.

## Complete annotation pipeline

**Goal:** Normalize, then annotate on MANE Select with calibrated predictors, then triage.

**Approach:** Left-align/split multiallelics, run VEP with MANE-led selection and predictor plugins, filter to HIGH/MODERATE for review (triage, not classification).

```bash
#!/bin/bash
set -euo pipefail
INPUT=$1; REFERENCE=$2; VEP_CACHE=$3; OUT=$4

bcftools norm -f "$REFERENCE" -m-any "$INPUT" -Oz -o "${OUT}_norm.vcf.gz"
bcftools index "${OUT}_norm.vcf.gz"

# --mane_select + constrained --pick_order so MANE leads and length/accession never decide
vep -i "${OUT}_norm.vcf.gz" -o "${OUT}_vep.vcf" \
    --vcf --cache --offline --dir_cache "$VEP_CACHE" --assembly GRCh38 \
    --everything --mane_select --pick --pick_order mane_select,canonical,biotype,rank --fork 4

bgzip "${OUT}_vep.vcf" && bcftools index "${OUT}_vep.vcf.gz"
bcftools view -i 'INFO/CSQ~"HIGH" || INFO/CSQ~"MODERATE"' \
    "${OUT}_vep.vcf.gz" -Oz -o "${OUT}_review.vcf.gz"
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Two labs report different c. for one indel | HGVS 3'-shift vs VCF left-align, strand-dependent | Normalize both to one representation before matching; never hand-derive HGVS from POS |
| "HIGH impact stop_gained" assumed PVS1 | Impact bucket ignores NMD, exon location, LOF mechanism | Check the NMD 50-nt/last-exon rule and gene mechanism before invoking PVS1 |
| Consequence changed after MANE migration | Ensembl canonical != MANE Select for many genes | Expected, not an error; communicate the coordinate change |
| Empty gnomAD annotations | Build mismatch (v2=GRCh37 vs v3/v4=GRCh38) or chr naming (chr1 vs 1) | Match build; `bcftools annotate --rename-chrs`; check site callability |
| VEP/SnpEff/ANNOVAR disagree on consequence | Different transcript set / splice width / severity ranking | Not a bug; pin one engine+version+transcript set and record it |
| `--pick` picked a non-clinical transcript | Length/accession tiebreaker fell through | Constrain `--pick_order mane_select,...` or pin a per-gene list |
| Predictors "all agree it's damaging" | Correlated tools (REVEL contains SIFT/PolyPhen) double-counted | Use ONE calibrated predictor at its calibrated strength |

## Related Skills

- variant-calling/variant-normalization - Left-align and split multiallelics; the mandatory step before annotation and HGVS
- variant-calling/clinical-interpretation - ACMG/AMP combining rules, PVS1/PP3/BP4 strengths, ClinVar star ratings, final classification
- variant-calling/filtering-best-practices - Filter by annotation and quality fields
- variant-calling/vcf-basics - Query annotated INFO/CSQ fields
- variant-calling/vcf-manipulation - Merge and manipulate annotated VCFs
- database-access/entrez-fetch - Download annotation databases (ClinVar, dbSNP)

## References

- McLaren W, et al. The Ensembl Variant Effect Predictor. *Genome Biology*. 2016;17:122. doi:10.1186/s13059-016-0974-4
- Cingolani P, et al. A program for annotating and predicting the effects of SNPs, SnpEff. *Fly (Austin)*. 2012;6(2):80-92. doi:10.4161/fly.19695
- Wang K, Li M, Hakonarson H. ANNOVAR: functional annotation of genetic variants. *Nucleic Acids Research*. 2010;38(16):e164. doi:10.1093/nar/gkq603
- Morales J, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE). *Nature*. 2022;604:310-315. doi:10.1038/s41586-022-04558-8
- Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants. *Bioinformatics*. 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112
- Abou Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Human Mutation*. 2018;39(11):1517-1524. doi:10.1002/humu.23626
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *American Journal of Human Genetics*. 2022;109(12):2163-2177. doi:10.1016/j.ajhg.2022.10.013
- Ioannidis NM, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *American Journal of Human Genetics*. 2016;99(4):877-885. doi:10.1016/j.ajhg.2016.08.016
- Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. 2023;381(6664):eadg7492. doi:10.1126/science.adg7492
- Kircher M, et al. A general framework for estimating the relative pathogenicity of human genetic variants (CADD). *Nature Genetics*. 2014;46:310-315. doi:10.1038/ng.2892
- Ng PC, Henikoff S. SIFT: predicting amino acid changes that affect protein function. *Nucleic Acids Research*. 2003;31(13):3812-3814. doi:10.1093/nar/gkg509
- Adzhubei IA, et al. A method and server for predicting damaging missense mutations (PolyPhen-2). *Nature Methods*. 2010;7(4):248-249. doi:10.1038/nmeth0410-248
- Jaganathan K, et al. Predicting splicing from primary sequence with deep learning (SpliceAI). *Cell*. 2019;176(3):535-548. doi:10.1016/j.cell.2018.12.015
- Zeng T, Li YI. Predicting RNA splicing from DNA sequence using Pangolin. *Genome Biology*. 2022;23:103. doi:10.1186/s13059-022-02664-4
- Karczewski KJ, et al. The mutational constraint spectrum quantified from variation in 141,456 humans (gnomAD v2.1.1). *Nature*. 2020;581(7809):434-443. doi:10.1038/s41586-020-2308-7
- Chen S, et al. A genomic mutational constraint map using variation in 76,156 human genomes (gnomAD v4). *Nature*. 2024;625(7993):92-100. doi:10.1038/s41586-023-06045-0
- Whiffin N, et al. Using high-resolution variant frequencies to empower clinical genome interpretation. *Genetics in Medicine*. 2017;19(10):1151-1158. doi:10.1038/gim.2017.26
