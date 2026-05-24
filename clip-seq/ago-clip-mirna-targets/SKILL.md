---
name: bio-clip-seq-ago-clip-mirna-targets
description: Identify direct miRNA-target interactions from AGO HITS-CLIP, AGO-CLEAR-CLIP (chimeric reads), HEAP (Halo-Ago2 mouse), chimeric eCLIP / miR-eCLIP (deep miRNA-target profiling), or CLASH using chimeric-read processing pipelines, seed-pairing analysis, and 3' auxiliary pairing rules. Use when distinguishing direct miRNA targets from indirect, integrating CLIP-derived target maps with TargetScan / miRDB / DIANA predictions, applying canonical 7mer-8mer seed matching with 3' UTR context, or recovering miRNA-mRNA chimeras at scale.
tool_type: mixed
primary_tool: chimeric-eCLIP
---

## Version Compatibility

Reference examples tested with: eCLIP pipeline (Yeo lab), chimeric eCLIP analysis scripts (Manakov 2022), HEAP pipeline (Li 2020), Hyb pipeline (Travis 2014), TargetScanHuman 8.0, miRDB 6.0, samtools 1.19+, bedtools 2.31+, pyHyb 0.4+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt the example to match the actual API rather than retrying.

# AGO-CLIP and miRNA Target Identification

**"Identify direct miRNA-target interactions experimentally"** -> Use Argonaute (AGO1-4) CLIP-seq variants to map miRNA-binding sites on mRNAs, then resolve which miRNA pairs with each site. Three approaches: (a) standard AGO-CLIP recovers AGO-bound sites but cannot say which miRNA; (b) chimeric methods (CLEAR-CLIP, chimeric eCLIP / miR-eCLIP) ligate the miRNA to its target during library prep, producing miRNA-mRNA chimeric reads that unambiguously assign miRNA-target pairs; (c) HEAP uses HaloTag-Ago2 for in vivo profiling. The chimeric methods are the gold standard for direct miRNA-target identification; standard AGO-CLIP must be combined with computational seed-matching (TargetScan, miRDB) to infer miRNA pairing. Resolution: chimeric reads pinpoint single miRNA-target pairs; AGO-only CLIP identifies "AGO-binding sites" of which a subset are miRNA targets.

- CLI (chimeric eCLIP / miR-eCLIP processing): custom pipeline starting from eCLIP-style preprocessing + chimeric-read identification + miRNA-mRNA junction extraction
- CLI (CLEAR-CLIP custom Moore 2015 pipeline): Hyb (Travis 2014) for chimera analysis
- CLI (Hyb pipeline): `hyb run_hyb peaks.bam mature_miRNA.fa human.tab.gz` to find miRNA-mRNA chimeras
- CLI (HEAP analysis): standard HITS-CLIP processing pipeline + Halo-Ago2 capture details
- Python (seed-pairing analysis on AGO CLIP peaks): scan peaks for canonical 7mer-m8, 7mer-1A, 8mer, 6mer seeds + 3' UTR position + miRNA expression filter

The Yeo lab miR-eCLIP / chimeric eCLIP (Manakov 2022) is the modern depth-improved version of chimeric AGO-CLIP, enriching for chimeras of specific miRNAs of interest (30-175x enrichment via PCR or on-bead probe capture). For comprehensive miRNA-target mapping, miR-eCLIP combined with eCLIP-seq-style normalization is the state-of-the-art.

## Methods Taxonomy

| Method | Year | miRNA-target pairing | Chimera enrichment | Strength | Fails when |
|--------|------|---------------------|--------------------|----------|------------|
| HITS-CLIP for AGO | 2009 (Chi) | Indirect (computational seed) | None | Original; widely cited | Cannot assign miRNA without computational prediction |
| PAR-CLIP for AGO | 2010 (Hafner) | Indirect | None | T->C signature at CL position | Restricted to 4SU-permissive cells |
| AGO-CLEAR-CLIP (Moore 2015) | 2015 | Direct (chimera) | None (incidental) | First direct miRNA-target chimera method | Chimeric reads only 1-5% of library; deep sequencing needed |
| CLASH (Helwak 2013) | 2013 | Direct (chimera) | None | First general chimera method; pan-Argonaute | Lower chimera rate than CLEAR-CLIP |
| HEAP (Li 2020) | 2020 | Indirect (with chimeric step) | None | HaloTag-Ago2 in vivo mouse strain | Mouse only; requires transgenic model |
| chimeric eCLIP / miR-eCLIP (Manakov 2022) | 2022 | Direct (chimera) | 30-175x enriched | Deepest miRNA-target chimera profiling | Specialized library prep |
| AGO HITS-CLIP + targeted chimeric (Bracken et al; verify exact venue/year) | -- | Direct | Yes | Per-miRNA enrichment via probe capture | Older; superseded by miR-eCLIP |
| AGO-IP-RNA-seq (Karginov 2007) | 2007 | Indirect | None | Earliest; predecessor of CLIP for AGO | No crosslinking; misses transient targets |

Methodology evolves; verify the Manakov 2022 / Bracken 2024 papers for current chimeric eCLIP best practice. As of 2024, miR-eCLIP is the canonical approach for deep miRNA-target profiling.

## Critical Choice: Chimeric vs Computational miRNA-Target Pairing

Two fundamentally different strategies:

**Chimeric methods (CLEAR-CLIP, chimeric eCLIP / miR-eCLIP, CLASH):** During library prep, a ligation step covalently joins the miRNA to its target mRNA, producing chimeric reads (miRNA at 5' + target mRNA at 3'). The miRNA-target pair is read directly from the sequence. Pro: direct evidence of binding interaction; no inference. Con: chimera rate is 1-5% of library by default (enriched to 30-175x with miR-eCLIP for specific miRNAs); deep sequencing or enrichment needed.

**Computational pairing (HITS-CLIP / PAR-CLIP + seed-matching):** Standard AGO CLIP identifies AGO-bound peaks; downstream computational scanning matches each peak against canonical miRNA seeds (7mer-m8, 7mer-1A, 8mer) from TargetScan, miRDB, or DIANA databases. Pro: any AGO CLIP data can be analyzed; no special library prep. Con: indirect; assigns miRNAs based on canonical seed rules, missing non-canonical interactions (3' compensatory, central pairing).

The Bartel lab CLEAR-CLIP analysis revealed substantial 3' auxiliary pairing beyond canonical seeds: ~50% of miRNA-target interactions have weak or non-canonical seed matching but strong 3' supplementary pairing. Chimeric methods recover these; computational seed-only inference misses them.

| Goal | Method |
|------|--------|
| Direct miRNA-target pair identification | Chimeric eCLIP / miR-eCLIP |
| Specific miRNA's targets (deep) | miR-eCLIP with probe-capture enrichment |
| All AGO-binding sites (any miRNA) | Standard AGO eCLIP / HITS-CLIP |
| In vivo mouse tissue | HEAP (Halo-Ago2 mouse) |
| Pan-Argonaute interactome | CLASH or chimeric eCLIP |
| Initial discovery / cost-conscious | AGO HITS-CLIP + TargetScan |
| Non-canonical / 3'-compensatory miRNA pairing | Chimeric methods (CLEAR-CLIP) |
| Comparison across species | TargetScan + AGO HITS-CLIP (computational) |
| Validate specific miRNA-target prediction | miR-eCLIP with that miRNA's probe |

## miRNA-Target Pairing Rules

Computational seed-matching against TargetScan / miRDB requires understanding the canonical miRNA-target pairing rules:

| Seed type | Pairing positions (miRNA nt) | Position 1 | Pro / Con |
|-----------|-------------------------------|------------|-----------|
| 8mer | 2-7 + position 8 + A at position 1 | A required | Strongest; most conserved targets |
| 7mer-m8 | 2-7 + position 8 (no A1 requirement) | Any | Strong; common |
| 7mer-A1 | 2-7 (no position 8) + A at position 1 | A required | Moderate; common |
| 6mer | 2-7 | Any | Weak; very common (many false positives) |
| 6mer-A1 | 2-6 + A at position 1 | A required | Weak |
| 3'-compensatory | Weak 6mer + strong 3' UTR pairing 12-17 | Any | Discovered via CLEAR-CLIP; misses in seed-only methods |
| Central pairing | Positions 4-15 with no seed | Any | Rare; cleavage rather than translational repression |

For TargetScan integration: download the TargetScanHuman 8.0 conserved-site predictions; filter for 7mer-8mer (drop 6mer if too noisy); cross-reference with the CLIP peak BED of the analysis.

## Chimeric eCLIP / miR-eCLIP Workflow

**Goal:** Recover miRNA-mRNA chimeras from AGO chimeric eCLIP / miR-eCLIP libraries and produce a per-miRNA target list suitable for direct biological interpretation.

**Approach:** Apply eCLIP-style preprocessing, then run Hyb in chimera (`type=mim`) mode with bowtie2 alignment (required for short 21-23 nt miRNA sequences), filter chimeras to human mRNA targets, intersect with miRNA-expression atlas (filter > 100 TPM in matched cell type), and validate top targets against TargetScan conserved 7mer-m8 / 8mer predictions.

```bash
# Step 1: eCLIP-style preprocessing (see clip-seq/clip-preprocessing)
umi_tools extract --bc-pattern=NNNNNNNNNN \
    --stdin=R1.fq.gz --read2-in=R2.fq.gz \
    --stdout=R1.umi.fq.gz --read2-out=R2.umi.fq.gz

cutadapt -a AGATCGGAAGAGCACACGTCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -q 6 -m 18 -o R1.trim.fq.gz -p R2.trim.fq.gz \
    R1.umi.fq.gz R2.umi.fq.gz

# Step 2: Chimera-specific alignment
# Chimeric reads have miRNA sequence (21-23 nt) at 5' followed by target mRNA
# Step 2a: Trim 5' for miRNA portion + align miRNA part
# Step 2b: Trim 3' for target mRNA portion + align target part
# Use custom chimeric-eCLIP pipeline OR Hyb (Travis 2014)

# Hyb pipeline approach (CLEAR-CLIP and chimeric methods)
hyb \
    in=R1.trim.fq.gz \
    db=miRNA_and_human_mRNA.fa \
    align=blastall \
    type=mim   # multimer (miRNA-target) chimera mode

# Output: .blast and .hyb files with miRNA-mRNA chimera coordinates

# Step 3: Filter chimeras by miRNA + target alignment quality
# Hyb reports each chimera as: miRNA_id  target_id  miRNA_alignment  target_alignment
# Filter for:
#   - miRNA portion 18-25 nt
#   - Target portion 18-50 nt
#   - miRNA-target seed match (7mer-m8 / 7mer-A1 / 8mer)
#   - Target alignment unique
awk '$5 == "human_mRNA"' chimeras.hyb > chimeras_human_mRNA.tsv

# Step 4: Aggregate chimeras into per-miRNA target list
# Each miRNA -> targets (with read counts as binding affinity proxy)
awk '{print $3, $4}' chimeras_human_mRNA.tsv | sort | uniq -c | sort -rn > mirna_target_counts.tsv
```

## CLEAR-CLIP (Moore 2015) Analysis

CLEAR-CLIP was the first method to recover ~130k miRNA-target chimeras from mouse brain (Moore 2015). The analytical insight: AGO-CLIP reads ligated together during library prep produce chimeric reads at low rates that contain unambiguous miRNA-target pairs.

```bash
# CLEAR-CLIP analysis uses the Hyb pipeline (Travis 2014)
# Pre-requisite: AGO HITS-CLIP / PAR-CLIP BAM

# Extract candidate chimeric reads (reads that don't fully align to human mRNA)
samtools view -h dedup.bam | awk '$6 ~ /S/' | wc -l   # soft-clipped reads candidate chimeras

# Run Hyb in chimera mode
hyb \
    in=R1.trim.fq.gz \
    db=mature_miRNA_plus_human_mRNA.fa \
    align=blastall \
    type=mim

# Filter for canonical seed match
python analyze_chimeras.py \
    --chimeras chimeras.hyb \
    --mirna_db mature_human_miRNA.fa \
    --seed_types 7mer-m8 7mer-A1 8mer \
    --output validated_chimeras.tsv
```

## miR-eCLIP Probe Enrichment (Manakov 2022)

To recover deep coverage of one or a few miRNAs' targets, miR-eCLIP uses probe-based or PCR-based enrichment to amplify chimeras containing specific miRNAs.

```bash
# After chimera identification, filter for specific miRNA
# Example: enrich for hsa-miR-21 targets
grep "hsa-miR-21" chimeras_human_mRNA.tsv > mir21_targets.tsv

# Count unique target sites per miRNA
awk '{print $3}' mir21_targets.tsv | sort -u | wc -l

# Cross-reference with TargetScan conserved predictions for validation
bedtools intersect -wa -wb \
    -a mir21_targets_3utr_coords.bed \
    -b targetscan_mir21_conserved_targets.bed > mir21_validated_targets.bed
```

## Per-Method Failure Modes

### Standard AGO-CLIP -- Cannot assign miRNA

**Trigger:** Standard AGO eCLIP / HITS-CLIP run; user wants per-miRNA target list.

**Mechanism:** Standard AGO-CLIP enriches for AGO-bound RNAs but does not retain miRNA identity. Computational seed-matching infers which miRNAs are likely bound but each peak gets matched to dozens of candidate miRNAs.

**Symptom:** Peak BED has 100k peaks; seed-matching assigns 10-50 candidate miRNAs per peak.

**Fix:** Switch to chimeric method for direct pairing, OR filter computational predictions by miRNA expression in the same cell type (only consider miRNAs > 100 TPM in matched small-RNA-seq).

### Chimeric methods -- Low chimera rate

**Trigger:** Standard chimeric eCLIP without probe enrichment; expecting deep per-miRNA targets.

**Mechanism:** Chimeras are 1-5% of total reads in standard chimeric eCLIP. For a 30M-read library, only 300k-1.5M chimeras; distributed across 200+ miRNAs gives only ~5000-15000 per miRNA.

**Symptom:** Per-miRNA target count is sparse; rare miRNAs have < 100 chimeras.

**Fix:** Use miR-eCLIP with probe enrichment for specific miRNAs of interest (30-175x boost). Or sequence ultra-deep (200M+ reads) for global chimera profiling.

### Hyb -- BLAST sensitivity vs miRNA length

**Trigger:** miRNA sequences (21-23 nt) too short for BLAST default sensitivity.

**Mechanism:** BLAST defaults need >= 100 nt for reliable alignment. miRNA 21-23 nt hits below threshold; many true chimeras lost.

**Symptom:** Hyb returns few chimeras; rerun with `align=bowtie2` gives more.

**Fix:** Use bowtie2 mode for miRNA alignment (`hyb align=bowtie2 type=mim`); short-read aligners are designed for short sequences.

### Computational seed matching -- High false positive

**Trigger:** TargetScan predictions used as ground truth without CLIP validation.

**Mechanism:** TargetScan reports all potential 7mer-m8 / 8mer matches in 3' UTRs; many sites are not functional miRNA targets (no AGO binding observed).

**Symptom:** TargetScan predicts thousands of targets per miRNA; only a fraction are validated by CLIP.

**Fix:** Use CLIP overlap as the validation: TargetScan prediction AND AGO-CLIP peak = high-confidence target. Sites in TargetScan but not in CLIP = unfunctional predictions.

### Non-canonical miRNA-target pairing missed

**Trigger:** Seed-matching only; 3' compensatory pairing missed.

**Mechanism:** ~50% of miRNA-target interactions have weak seeds but strong 3' UTR pairing (positions 12-17). Seed-only matching loses these.

**Symptom:** Chimeric methods find targets that TargetScan misses; these have weak seeds.

**Fix:** Accept chimeric method's targets even with weak seeds (the chimera IS the evidence). Or use TargetScan + RNAhybrid (full miRNA-target duplex prediction) for non-canonical sites.

### HEAP -- Mouse-only

**Trigger:** Want HEAP-style in vivo AGO profiling in human tissue.

**Mechanism:** HEAP uses a transgenic mouse with Halo-Ago2 allele; not available in human or other species.

**Symptom:** Cannot replicate HEAP results in human.

**Fix:** Use eCLIP / chimeric eCLIP on human samples; HEAP is specifically for mouse tissue studies.

### miRNA expression filter forgotten

**Trigger:** Computational miRNA-target assignment without filtering by miRNA expression.

**Mechanism:** Many miRNA databases include rare or developmental-specific miRNAs. If the miRNA is not expressed in the cell type, it cannot bind anything.

**Symptom:** Per-miRNA target lists include miRNAs at < 1 TPM expression - implausible binding.

**Fix:** Cross-reference with matched small-RNA-seq from the same cell type; filter for miRNAs > 100 TPM. ENCODE-validated cell types have published miRNA atlases.

## Decision Tree by Use Case

| Scenario | Method | Why |
|----------|--------|-----|
| Direct miRNA-target identification, modern | chimeric eCLIP / miR-eCLIP | Direct chimeras; deep enrichment available |
| Specific miRNA's deep target list | miR-eCLIP with probe for that miRNA | 30-175x enrichment |
| Discover novel miRNA-target interactions | CLEAR-CLIP or chimeric eCLIP | Direct chimera, no seed prior |
| In vivo mouse tissue | HEAP (Halo-Ago2 mouse) | Mouse only |
| Initial AGO-binding site discovery | AGO HITS-CLIP / eCLIP | Cost-effective; no chimera |
| Compare miRNA targets across species | TargetScan + AGO HITS-CLIP each species | Computational + experimental |
| 3' compensatory / non-canonical | CLEAR-CLIP / chimeric eCLIP | Direct chimera captures non-canonical |
| miRNA-perturbation effects | KD/KO + AGO-CLIP + differential | See clip-seq/differential-clip |
| Cross-tissue miRNA profiling | AGO eCLIP each tissue | Tissue-specific cell-type |
| Validate single miRNA prediction | miR-eCLIP with that miRNA's probe | Direct experimental confirmation |

## Reconciliation: AGO-CLIP vs TargetScan vs Chimeric

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Chimera method finds targets TargetScan does not | Non-canonical / 3'-compensatory pairing | Trust chimera; novel target |
| TargetScan predicts; AGO eCLIP peak present; no chimera | Functional target without chimera in library | Likely real target; chimera capture stochastic |
| TargetScan predicts; no AGO eCLIP peak | Computational false positive | Not a functional target |
| AGO eCLIP peak; no TargetScan match | Non-canonical or rare miRNA seed | Investigate; may be 3' compensatory |
| Per-miRNA target counts vary 100x across miRNAs | miRNA expression varies | Filter by matched small-RNA-seq |
| Hyb chimeras 1% of library | Standard rate | Enrich with miR-eCLIP if needed |
| Different chimera tools give different counts | Algorithm sensitivity differs | Hyb is the most-cited; use it for canonical |
| HEAP and eCLIP discordant | Mouse vs human; in vivo vs cell line | Both correct in their context |
| miR-eCLIP enriched chimera count not 30x baseline | Probe inefficient | Verify probe design; use multiple probes per miRNA |

**Operational rule:** For publication-grade miRNA-target list: (a) chimeric eCLIP / miR-eCLIP for direct pairing; (b) cross-reference with TargetScan conserved predictions; (c) filter by miRNA expression > 100 TPM in matched small-RNA-seq; (d) validate top targets with reporter assay (luciferase / GFP fusion with target 3' UTR).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Hyb returns few chimeras | BLAST too stringent for short miRNAs | Use bowtie2 mode (`hyb align=bowtie2`) |
| Per-miRNA target list sparse | Low chimera rate without enrichment | Use miR-eCLIP probe enrichment |
| TargetScan predicts thousands per miRNA | No CLIP filter | Require CLIP peak overlap for high-confidence |
| miRNA assignments dominate by unexpressed miRNAs | No expression filter | Filter by matched small-RNA-seq > 100 TPM |
| Non-canonical sites missed | Seed-only matching | Use chimeric methods; or RNAhybrid full duplex |
| HEAP results don't replicate in human | Mouse-specific transgenic | Use eCLIP / chimeric eCLIP in human |
| 6mer matches dominate target list | Most weak seeds | Restrict to 7mer-m8 / 8mer; report 6mer separately |
| miRNA target chimeras unstrand-resolved | Strand information lost | Check BED column 6 throughout pipeline |
| Cross-tissue comparison naive | Tissue-specific miRNA expression | Match tissue-specific miRNA atlases |
| miR-eCLIP enrichment fails | Probe non-specific or low-affinity | Design multiple probes per miRNA; validate enrichment |

## References

- Chi SW et al 2009 Nature 460:479 (AGO HITS-CLIP)
- Hafner M et al 2010 Cell 141:129 (PAR-CLIP for AGO)
- Helwak A et al 2013 Cell 153:654 (CLASH; chimera method)
- Travis AJ et al 2014 Methods 65:263 (Hyb pipeline)
- Moore MJ et al 2015 Nat Commun 6:8864 (CLEAR-CLIP, 130k chimeras mouse brain)
- Bracken CP et al -- chimeric AGO-CLIP, targeted (consult current literature for verified venue/year; earlier "2016 Nat Methods 13:739" attribution could not be confirmed).
- Li K et al 2020 Mol Cell 80:1100 (HEAP, Halo-Ago2 in vivo mouse)
- Manakov SA et al 2022 bioRxiv 2022.02.13.480296 (chimeric eCLIP / miR-eCLIP, 30-175x enrichment)
- Agarwal V et al 2015 eLife 4:e05005 (TargetScan 7.0)
- Lewis BP et al 2003 Cell 115:787 (original 7mer/8mer seed rules)
- Bartel DP 2018 Cell 173:20 (miRNA target principles)
- McGeary SE et al 2019 Science 366:eaav1741 (TargetScan 8.0 / quantitative target prediction).

## Related Skills

- clip-seq/clip-peak-calling - AGO CLIP peak calls
- clip-seq/binding-site-annotation - 3' UTR annotation
- clip-seq/clip-motif-analysis - Seed motif scan
- clip-seq/differential-clip - miRNA perturbation experiments
- clip-seq/m6a-clip - DART-seq uses similar APOBEC1 fusion
- small-rna-seq/target-prediction - TargetScan / miRDB / DIANA
- small-rna-seq/differential-mirna - miRNA expression
- small-rna-seq/mirdeep2-analysis - miRNA discovery
