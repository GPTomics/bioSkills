---
name: bio-clip-seq-m6a-clip
description: Map N6-methyladenosine (m6A) RNA modifications at single-nucleotide resolution using miCLIP (Linder 2015), miCLIP2 + m6Aboost machine learning (Kortel 2021), GLORI (Liu 2023, antibody-free chemical conversion), DART-seq (Meyer 2019, APOBEC1-YTH fusion), m6Anet (nanopore direct RNA), or MeRIP-seq with calibration. Use when distinguishing antibody-based from antibody-free m6A detection methods, applying the DRACH motif constraint, reconciling cross-method disagreements (DART 44% in DRACH vs GLORI), or detecting m6Am at the cap.
tool_type: mixed
primary_tool: miCLIP2
---

## Version Compatibility

Reference examples tested with: miCLIP2 pipeline (Kortel 2021), m6Aboost 1.0+, GLORI-tools (Liu 2023), Bullseye 1.0+, m6Anet 2.1+, EpiNano 1.2+, MeRIPSeq tools (exomePeak2 1.16+), nanocompore 1.0+, samtools 1.19+, bedtools 2.31+, R 4.3+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt the example to match the actual API rather than retrying.

# m6A CLIP (N6-Methyladenosine Profiling)

**"Map m6A modifications at single-nucleotide resolution"** -> Profile m6A on RNA using one of three orthogonal approaches: antibody-based UV-CL (miCLIP/miCLIP2), antibody-free chemical conversion (GLORI), or enzyme-fusion editing (DART-seq with APOBEC1-YTH). Nanopore direct RNA (m6Anet, nanocompore, EpiNano) provides a fourth modality. The DRACH consensus motif (D=A/G/U, R=A/G, A=m6A, C=C, H=A/C/U) constrains plausible sites but is not exclusive - only a fraction of DRACH instances are methylated; some m6A sites occur outside DRACH. Cross-method discordance is real: DART-seq C->U mutations are 44% within DRACH motifs (Liu 2023), suggesting many DART sites are not consensus m6A. GLORI is the new (2023) gold standard for stoichiometric single-base m6A.

- CLI (miCLIP2 antibody-based): `iCount` or custom pipeline through truncation + C->T mutation analysis; then m6Aboost ML scoring
- CLI (GLORI antibody-free): `GLORI-tools` Python pipeline; output is per-A m6A fraction (stoichiometric)
- CLI (DART-seq editing): `Bullseye` or `SAILOR` pipeline; identify C->U editing sites; filter by DRACH; cross-check against APOBEC1-only control
- CLI (m6Anet nanopore): `m6anet inference` on nanopolish eventalign output; per-site probability of m6A
- CLI (MeRIP-seq peak calling): `exomePeak2` in R for peak-level m6A from IP+input MeRIP libraries

The m6A field is rapidly evolving (2022-2026); single-base methods (GLORI, m6Anet) have largely replaced antibody-based miCLIP for new studies, but miCLIP2 remains the most common because of its eCLIP-like processing pipeline. Cross-method discordance means high-confidence m6A reporting should require concordance across at least two orthogonal methods.

## Methods Taxonomy

| Method | Detection chemistry | Resolution | Antibody | Stoichiometry | Strength | Fails when |
|--------|---------------------|------------|----------|---------------|----------|------------|
| MeRIP-seq (Dominissini 2012, Meyer 2012) | Anti-m6A IP + RNA-seq | Peak (50-300 nt) | Yes | No | Original m6A method; widely used | Low resolution; cannot distinguish m6A from m6Am |
| miCLIP (Linder 2015) | Anti-m6A + UV-CL + RT mutation | Single-nucleotide (some) | Yes | No | Single-nt subset of m6A peaks | Low yield of single-nt; high false-positive rate |
| miCLIP2 (Kortel 2021) | Anti-m6A + UV-CL + improved library | Single-nucleotide | Yes | No | Higher complexity; ML-classified (m6Aboost) | Antibody specificity remains issue |
| GLORI (Liu 2023) | Glyoxal + nitrite chemical conversion of A to N1-methyl-2-amino-6-oxopurine | Single-nucleotide | No (chemical) | Yes (stoichiometric) | Stoichiometric m6A fraction per site | New; less validated; harsh conversion may damage rare RNAs |
| DART-seq (Meyer 2019) | APOBEC1-YTH fusion edits C adjacent to m6A | Single-nucleotide (offset) | No | No | Antibody-free; in vivo | Only 44% of edits in DRACH motifs; high false positive |
| m6A-CLIP (Ke 2015) | Anti-m6A + UV-CL | Peak | Yes | No | Original UV-CL approach | Predecessor to miCLIP |
| m6Anet (Hendra 2022) | Nanopore direct RNA + neural net | Single-nucleotide (DRACH constraint) | No | Probability | Direct RNA; preserves isoform context | Restricted to DRACH; needs high coverage per site |
| EpiNano (Liu 2019) | Nanopore + SVM on signal features | Single-nucleotide | No | No | Pioneer nanopore m6A | Lower accuracy than m6Anet on benchmark |
| nanocompore (Leger 2021) | Nanopore + statistical test wt vs Mettl3-KO | Single-nucleotide | No | No | Comparative; high specificity | Requires KO control sample |
| DENA (Zhong 2024) | Nanopore + transformer model | Single-nucleotide | No | No | Single-sample tool | Newer; less validation |
| FTO/ALKBH5-aware methods | Eraser perturbation | Site | No | Indirect | Validates m6A regulation | Indirect |
| MAZTER-seq (Garcia-Campos 2019) | RNase T1 cleavage at unmodified ACA | Site (within ACA) | No | No | Antibody-free | Restricted to ACA context (subset of DRACH) |
| REF-seq (Werner 2020) | Endonuclease-cleavage | Site | No | No | Antibody-free | Restricted context |
| m6ACali (Yang 2024) | Calibrates MeRIP / miCLIP | Site | NA | Yes (calibration) | Cross-method calibration | Postprocessing only |

Methodology evolves; verify the latest benchmark publications (e.g., Tegowski 2024 Mol Cell review). The field is moving toward GLORI as the new gold standard but miCLIP2 remains the most-cited method because of its eCLIP-pipeline compatibility.

## Critical Choice: Antibody-Based vs Antibody-Free

**Antibody-based (MeRIP-seq, miCLIP, miCLIP2, m6A-CLIP):** Anti-m6A antibody (Abcam/Synaptic Systems) immunoprecipitates m6A-bearing RNA. The antibody is the only limitation - false positives from non-specific binding to long structured RNAs (especially poly-A) and false negatives at sites with low m6A stoichiometry. Mettl3 knockout calibration is recommended.

**Antibody-free chemical (GLORI):** Glyoxal + nitrite converts unmodified A to a nucleotide that reads as G; m6A is protected and reads as A. Sites are detected as A->G discrepancies post-conversion. Stoichiometric (the fraction of reads showing A vs G at a position = m6A fraction). Most rigorous but chemistry is harsh - degrades very long RNAs.

**Antibody-free enzymatic (DART-seq, APOBEC1-YTH):** APOBEC1 cytidine deaminase fused to YTH-domain (m6A reader) edits C residues adjacent to m6A. Editing pattern (C->U) marks m6A nearby but not exactly. 44% of DART edits in DRACH; many edits are off-target.

**Antibody-free nanopore (m6Anet, nanocompore, EpiNano):** Direct RNA sequencing detects m6A via current signal perturbation. Preserves isoform context. m6Anet AUC 0.83 on HEK293T; outperforms EpiNano and Tombo on benchmark (Liu 2024).

| Goal | Method |
|------|--------|
| Stoichiometric m6A fraction per site | GLORI |
| eCLIP-compatible processing pipeline | miCLIP2 + m6Aboost |
| Isoform-resolved m6A | m6Anet (nanopore) |
| Cell-line comparison (KO available) | nanocompore vs Mettl3-KO |
| High-throughput screening | DART-seq (in vivo, no UV) |
| Initial discovery (low cost) | MeRIP-seq (with calibration) |
| Variants in m6A context | GLORI + variant-effect analysis |
| Combined m6A + 5'-cap m6Am | miCLIP2 (detects both with separate motifs) |

## DRACH Motif Constraint

The DRACH consensus (D=A/G/U, R=A/G, A=m6A, C=C, H=A/C/U) is the dominant motif at m6A sites - 70-90% of high-confidence sites fall in DRACH context. But:
- Some m6A sites occur outside DRACH (~10-20% in calibrated datasets)
- Many DRACH instances are NOT methylated (only a subset)
- Filtering for DRACH-only loses 10-20% of sites; not-filtering inflates false positives

**miCLIP2 + m6Aboost (Kortel 2021)** trained on Mettl3 knockout calibration data to score sites without DRACH filtering. The m6Aboost ML model is the recommended approach when DRACH-blind detection is needed.

**GLORI** does not filter by DRACH; the per-A m6A fraction is reported regardless of context. The non-DRACH GLORI sites (10-20%) include genuine m6A in non-canonical context.

## Cross-Method Discordance

| Comparison | Concordance | Source |
|------------|-------------|--------|
| miCLIP vs miCLIP2 | ~70% | Kortel 2021 |
| miCLIP2 vs GLORI | ~60% (miCLIP2 calls in GLORI) | Liu 2023 |
| GLORI vs MeRIP-seq peaks | ~50% sites in MeRIP peaks | Liu 2023 |
| DART-seq vs GLORI | 44% DART edits in DRACH; 7-28% near GLORI sites | Liu 2023 |
| m6Anet vs miCLIP2 | ~75% concordance at high-coverage sites | Hendra 2022 |
| Antibody-based methods | High discordance between antibody lots | Tegowski 2024 |

**Reconciliation strategy:** Use GLORI as the new gold standard (2023+); cross-reference with m6Anet for nanopore isoform context; treat miCLIP2 + m6Aboost as a complementary in vivo perspective; treat DART-seq as a hypothesis-generating method. Three orthogonal methods agreeing on a site = high confidence.

## miCLIP2 Workflow

miCLIP2 (Kortel 2021) is the eCLIP-pipeline-compatible m6A method. It uses anti-m6A antibody + UV-CL + improved library prep that increases complexity 10x over miCLIP.

**Goal:** Produce a high-confidence single-nucleotide m6A site BED from anti-m6A miCLIP2 reads with antibody-false-positive suppression via m6Aboost machine learning.

**Approach:** Run the eCLIP-style preprocessing + STAR + UMI-dedup pipeline, call single-nt CL sites with PureCLIP using SMInput control, then apply m6Aboost (trained on Mettl3-KO calibration data) to discriminate genuine m6A sites from antibody false positives without requiring strict DRACH motif filtering.

```bash
# Step 1: Preprocessing (eCLIP-style - see clip-seq/clip-preprocessing)
umi_tools extract --bc-pattern=NNNNNNNNNN \
    --stdin=R1.fq.gz --read2-in=R2.fq.gz \
    --stdout=R1.umi.fq.gz --read2-out=R2.umi.fq.gz

cutadapt -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 6 -m 18 \
    -o R1.trim.fq.gz -p R2.trim.fq.gz \
    R1.umi.fq.gz R2.umi.fq.gz

# Step 2: Alignment (eCLIP-style)
STAR --runMode alignReads --genomeDir STAR_index \
    --readFilesIn R1.trim.fq.gz R2.trim.fq.gz --readFilesCommand zcat \
    --alignEndsType EndToEnd --outFilterMultimapNmax 1 --outFilterMismatchNoverReadLmax 0.04 \
    --outSAMtype BAM SortedByCoordinate

umi_tools dedup --method=unique --paired -I aligned.bam -S dedup.bam

# Step 3: Single-nt CL site detection - PureCLIP or custom
pureclip -i dedup.bam -bai dedup.bam.bai -g genome.fa \
    -ibam sminput.bam -ibai sminput.bam.bai \
    -o miCLIP2_sites.bed -or miCLIP2_regions.bed -nt 8

# Step 4: m6Aboost ML scoring (Kortel 2021)
# Requires: site BED + features (sequence context, C->T rate, truncation position)
# Trained on Mettl3 KO calibration data
# Output: m6A probability score per site
python m6aboost.py \
    --sites miCLIP2_sites.bed \
    --bam dedup.bam \
    --genome genome.fa \
    --output m6Aboost_predictions.bed

# Step 5: Filter at m6Aboost score >= 0.5 (default; tune per study)
awk '$5 >= 0.5' m6Aboost_predictions.bed > m6a_high_confidence.bed
```

## GLORI Workflow (Antibody-Free Stoichiometric)

GLORI (Liu 2023) is the new (2023) gold-standard for stoichiometric m6A. Chemistry: glyoxal + nitrite converts unmodified A; m6A is protected.

```bash
# GLORI-tools pipeline (Liu lab, github). GLORI-tools is a multi-step Python pipeline
# (`run_GLORI.py` is the typical orchestrator); the conceptual flow below is illustrative --
# verify the exact CLI against the GLORI-tools repo before scripting.
# 1. Pre-conversion sequencing (control)
# 2. Post-conversion sequencing (treated)
# 3. GLORI-tools computes per-A m6A fraction

python run_GLORI.py \
    --input pre_conversion.bam \
    --treated post_conversion.bam \
    --reference genome.fa \
    --output glori_sites.tsv

# Output columns: chr, pos, strand, m6A_fraction, coverage, p_value
# m6A_fraction: 0.0 = unmodified; 1.0 = fully methylated
# Filter at coverage >= 20 and m6A_fraction >= 0.1
awk 'NR>1 && $5 >= 20 && $4 >= 0.1' glori_sites.tsv > glori_high_confidence.tsv
```

## DART-seq Workflow (Editing-Based)

DART-seq (Meyer 2019) expresses APOBEC1-YTH fusion in cells; the YTH domain binds m6A, APOBEC1 edits nearby Cs.

```bash
# Bullseye pipeline (Meyer lab github)
# Requires APOBEC1-only (no YTH) control to subtract off-target editing
Bullseye \
    --ip dart_sample.bam \
    --control apobec1_only.bam \
    --reference genome.fa \
    --output dart_sites.bed

# Filter for DRACH motif overlap (44% of DART sites are in DRACH)
# Sites outside DRACH may be off-target editing
bedtools intersect -wa -u -s -a dart_sites.bed -b drach_motifs.bed > dart_drach_sites.bed
```

## m6Anet Workflow (Nanopore)

m6Anet (Hendra 2022) is the leading nanopore direct-RNA m6A detector. Uses signal-level features in a multiple-instance learning framework.

```bash
# Step 1: nanopolish eventalign on raw nanopore signal
# (assumes basecalled FASTQ, aligned BAM, raw FAST5/POD5)
nanopolish eventalign \
    --reads basecalled.fastq \
    --bam aligned.bam \
    --genome transcriptome.fa \
    --scale-events --signal-index --samples \
    > eventalign.tsv

# Step 2: m6Anet feature extraction
m6anet dataprep \
    --eventalign eventalign.tsv \
    --out_dir m6anet_features \
    --n_processes 8

# Step 3: m6Anet inference
m6anet inference \
    --input_dir m6anet_features \
    --out_dir m6anet_out \
    --pretrained_model HEK293T_RNA002

# Output: per-site probability of m6A
# Filter at probability_modified >= 0.9 (high confidence)
awk -F'\t' 'NR>1 && $5 >= 0.9' m6anet_out/data.indiv_proba.csv > m6Anet_high.tsv
```

## Per-Method Failure Modes

### miCLIP / miCLIP2 -- Antibody specificity

**Trigger:** Antibody lot variation; off-target binding to structured non-methylated RNAs.

**Mechanism:** Anti-m6A antibody (Abcam, Synaptic Systems) has variable specificity. Long structured RNAs (especially poly-A regions, snRNAs) capture non-specifically. False-positive rate without Mettl3-KO calibration is 30-50%.

**Symptom:** miCLIP sites overlap with snRNAs and long ncRNAs at unexpected rates; m6Aboost predicts < 30% of sites are true m6A.

**Fix:** Always include Mettl3-KO calibration (m6Aboost was trained on this). Apply m6Aboost ML; do not just filter by DRACH. Or switch to antibody-free GLORI.

### GLORI -- RNA degradation

**Trigger:** GLORI on long RNAs (> 5 kb); high glyoxal+nitrite concentration.

**Mechanism:** Harsh chemistry damages long RNAs; coverage at long transcripts drops 50-80% post-conversion.

**Symptom:** Long transcripts (e.g., Titin) have poor coverage post-GLORI; m6A sites in coding regions of long mRNAs under-called.

**Fix:** Use shorter conversion times for long-RNA studies (4 h vs 24 h); accept reduced power on long transcripts; cross-reference with miCLIP2 for long-RNA m6A.

### DART-seq -- Off-target editing

**Trigger:** APOBEC1-YTH expressed in cells; no APOBEC1-only control.

**Mechanism:** APOBEC1 has intrinsic C->U editing activity independent of YTH-m6A binding. Without APOBEC1-only control, 30-50% of edits are off-target.

**Symptom:** Many DART edits in non-DRACH context (44% in DRACH per Liu 2023); GO term enrichment of edited genes is non-specific.

**Fix:** Always run APOBEC1-only control in parallel; subtract its edits. Filter for DRACH motif overlap when reporting. Cross-validate with miCLIP2 or GLORI.

### m6Anet -- Coverage requirement

**Trigger:** Nanopore direct RNA on a low-input sample; per-site coverage < 20 reads.

**Mechanism:** m6Anet's multiple-instance learning needs >= 20 reads per DRACH position for stable probability estimate.

**Symptom:** Many "not enough coverage" sites in m6Anet output; gene-level coverage uneven.

**Fix:** Increase nanopore flowcell yield; pool replicates; restrict analysis to high-expression transcripts (TPM >= 5).

### MeRIP-seq -- Peak-level resolution

**Trigger:** MeRIP-seq on antibody-based platforms; peak width 100-300 nt.

**Mechanism:** MeRIP fragments are 100-300 nt; the peak captures a region containing m6A but cannot pinpoint the exact A.

**Symptom:** Peak BED width > 100 nt; downstream single-nt analysis impossible.

**Fix:** Combine MeRIP-seq with single-nt method (GLORI, miCLIP2). Or use m6ACali (Yang 2024) for cross-method calibration.

### DRACH-only filter -- Misses non-canonical m6A

**Trigger:** Filtered miCLIP2 / DART sites to DRACH-only.

**Mechanism:** 10-20% of validated m6A sites are outside DRACH context.

**Symptom:** Lost some validated sites; published m6A list shorter than expected.

**Fix:** Use m6Aboost (DRACH-blind ML) or GLORI (DRACH-blind chemical). Report both DRACH-filtered and unfiltered sets.

### Cross-method discordance frustration

**Trigger:** Three methods produce three different m6A site lists; user wants ONE truth.

**Mechanism:** Methods have different chemistries, sensitivities, and biases. They are not interchangeable. Discordance is real biology + technical.

**Symptom:** Two papers on the same RNA report different m6A sites.

**Fix:** Triangulate. Report (a) high-confidence sites from any single rigorous method (GLORI preferred); (b) consensus sites across 2+ methods. Acknowledge method limitations.

## Decision Tree by Use Case

| Scenario | Method | Why |
|----------|--------|-----|
| New 2024+ study, gold-standard single-base | GLORI | Stoichiometric, antibody-free |
| eCLIP-pipeline-compatible processing | miCLIP2 + m6Aboost | Uses eCLIP infrastructure |
| Isoform-resolved m6A | m6Anet (nanopore) | Long reads preserve isoforms |
| Mettl3 KO calibration available | miCLIP2 + m6Aboost; OR nanocompore | KO is the m6A negative control |
| In vivo, no UV | DART-seq | No UV CL needed |
| Initial discovery (low cost) | MeRIP-seq + exomePeak2 + m6ACali | Cheapest |
| Long RNAs (> 5 kb) | miCLIP2 or m6Anet (not GLORI) | GLORI degrades long RNAs |
| Variant in m6A context | GLORI single-base + variant-effect | Stoichiometric reveals dosage |
| m6Am at 5' cap | miCLIP2 (distinguishes via context) | The 5'-cap-adjacent A |
| Bacterial m6A | Custom methods | Mammalian DRACH irrelevant |
| Time-course m6A dynamics | GLORI per time point | Stoichiometric quantitation |
| Cross-species m6A | Use method validated in that species | Generalization not assumed |

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| miCLIP2 calls site; GLORI does not | Antibody false positive; or m6A fraction low | Trust GLORI for stoichiometry; flag miCLIP2 site for re-validation |
| GLORI calls site; miCLIP2 does not | Antibody false negative (saturation); or non-DRACH | Trust GLORI; check DRACH context of miCLIP2 site |
| DART edits not in DRACH | Off-target APOBEC1 editing | Subtract APOBEC1-only control; filter for DRACH |
| m6Anet calls site; miCLIP2 does not | Nanopore signal-specific detection; complementary | Cross-validate with GLORI; nanopore is orthogonal |
| MeRIP peak but no single-base call within | Peak captures multiple low-stoichiometry sites OR antibody non-specific | Use single-base method for confirmation |
| Discordance between antibody lots | Specificity variation | Use ENCODE-validated antibody; document lot |
| Cross-species method comparison | Method validated only in HEK293 / mouse | Re-validate before applying |
| Time-course shows decrease, methods disagree on magnitude | Stoichiometric (GLORI) vs fraction-based (miCLIP) | GLORI is quantitative; miCLIP is binary call |

**Operational rule for high-confidence m6A reporting:** (a) Use GLORI for stoichiometric single-base sites where chemistry permits; (b) Use miCLIP2 + m6Aboost where eCLIP-pipeline compatibility is required; (c) Use m6Anet for isoform-resolved or long RNAs; (d) Require concordance across at least two orthogonal methods for any m6A site claimed in publication.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| miCLIP2 sites > 100k - more than realistic m6A count | No m6Aboost ML scoring | Apply m6Aboost; expect 10-50k high-confidence |
| GLORI coverage uneven across transcripts | Glyoxal harshness on long RNAs | Shorter conversion; or use other methods for long RNAs |
| DART edits everywhere | No APOBEC1-only subtraction | Add APOBEC1-only control |
| m6Anet returns "no sites" | Coverage < 20 per DRACH | Pool replicates; restrict to high-expression transcripts |
| 10-20% sites outside DRACH | Real biology + some false positives | Report DRACH and non-DRACH separately |
| MeRIP peaks > 200 nt wide | Method resolution | Use single-base method for single-nt sites |
| Different methods give different sites | Method-specific biases | Triangulate; cross-validate |
| Antibody lot variation in miCLIP | Specificity drift | Document lot; use Mettl3-KO calibration |
| m6Am detection failing in miCLIP2 | Failed at 5'-cap | Check 5'-cap adjacent context filter |
| Cross-method calibration confusing | m6ACali heuristic | Apply pre-publication; verify with m6Aboost |

## References

- Dominissini D et al 2012 Nature 485:201 (MeRIP-seq)
- Meyer KD et al 2012 Cell 149:1635 (MeRIP-seq concurrent)
- Linder B et al 2015 Nat Methods 12:767 (miCLIP)
- Ke S et al 2015 Genes Dev 29:2037 (m6A-CLIP)
- Kortel N et al 2021 Nucleic Acids Res 49:e92 (miCLIP2 + m6Aboost)
- Liu C et al 2023 Nat Biotechnol 41:355 (GLORI)
- Meyer KD 2019 Nat Methods 16:1275 (DART-seq)
- Hendra C et al 2022 Nat Methods 19:1590 (m6Anet)
- Liu Y et al 2019 Nat Commun 10:4079 (EpiNano)
- Leger A et al 2021 Nat Commun 12:7198 (nanocompore)
- Garcia-Campos MA et al 2019 Cell 178:731 (MAZTER-seq)
- Tegowski M et al 2024 Mol Cell 84:1135 (m6A methods review)
- Yang Y et al 2024 Briefings Bioinform 25:bbae001 (m6ACali, benchmark)

## Related Skills

- clip-seq/clip-preprocessing - miCLIP2 uses eCLIP-style preprocessing
- clip-seq/clip-alignment - miCLIP2 uses eCLIP-style alignment
- clip-seq/crosslink-site-detection - miCLIP2 single-nt CL detection
- clip-seq/clip-peak-calling - MeRIP-seq exomePeak2 peak calling
- clip-seq/stamp-antibody-free - STAMP / DART-seq antibody-free approach
- long-read-sequencing/nanopore-methylation - Native nanopore m6A
- long-read-sequencing/basecalling - dRNA-seq basecalling
- epitranscriptomics/m6a-peak-calling - MeRIP-specific peak calling
- epitranscriptomics/m6a-differential - Differential m6A
- epitranscriptomics/m6anet-analysis - Nanopore m6Anet workflow
