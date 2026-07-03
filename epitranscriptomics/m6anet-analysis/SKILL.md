---
name: bio-epitranscriptomics-m6anet-analysis
description: Detects m6A modifications from Oxford Nanopore direct-RNA-seq (DRS) signal using m6Anet (multiple-instance-learning over DRACH 5-mer signal). Covers the upstream pipeline (Dorado/Guppy basecalling -> minimap2 map-ont -> nanopolish eventalign -> m6anet dataprep -> m6anet inference), per-site vs per-read probability including the mod_ratio stoichiometry column, the DRACH-only constraint, minimum-coverage thresholds (20-50 reads/site), multi-condition comparison via xPore/Nanocompore/ELIGOS, Dorado native modification calling (RNA004, 2024+), and the cDNA-vs-DRS distinction (cDNA Nanopore CANNOT detect modifications). Use when calling m6A from ONT DRS without immunoprecipitation, choosing m6Anet vs xPore vs Nanocompore vs ELIGOS vs Dorado native, interpreting probability_modified vs mod_ratio vs per-read probabilities, deciding between m6Anet (known DRACH sites) and Dorado/Remora (genome-wide screening), pinning RNA002 vs RNA004 chemistry and basecaller versions, or troubleshooting eventalign/dataprep failures.
tool_type: python
primary_tool: m6Anet
---

## Version Compatibility

Reference examples tested with: m6anet 2.1+ (PyPI; project capitalisation `m6Anet`), nanopolish 0.14+, minimap2 2.26+, samtools 1.19+, Dorado 0.5+, xpore 2.1+, nanocompore 1.0.4+, ELIGOS 2 (GitHub `novoalab/Eligos2`), CHEUI (GitHub `comprna/CHEUI`), pandas 2.2+, pyranges 0.0.129+, pysam 0.22+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show m6anet` then `m6anet --help`, `m6anet dataprep --help`, `m6anet inference --help`
- CLI: `nanopolish --version`, `minimap2 --version`, `dorado --version`

If `m6anet` is invoked as a hyphenated command (`m6anet-dataprep`), the v1.x CLI is installed; v2.x uses subcommand syntax (`m6anet dataprep`). Pin m6anet version explicitly in any reproducible analysis. Dorado modification models are versioned independently of the basecaller (`m6A_DRACH@v0.1`, `m6A_DRACH@v1`, ...); model versions are NOT directly comparable across releases. RNA002 vs RNA004 nanopore chemistry use different signal characteristics — models trained on RNA002 do NOT transfer to RNA004.

# m6Anet Direct-RNA m6A Detection

**"Detect m6A from my Nanopore direct RNA data without IP"** -> Run the full pipeline from POD5 / FAST5 signal to per-site m6A modification probabilities: Dorado basecalling -> minimap2 alignment to TRANSCRIPTOME (not genome) with `-ax map-ont -uf -k14` -> nanopolish eventalign with `--scale-events --signal-index` (the m6Anet-required pair) -> m6anet dataprep -> m6anet inference. Report per-site `probability_modified` (model posterior that any reads at the site are modified) AND `mod_ratio` (per-site stoichiometry, the fraction of reads called modified) with coverage. CRITICAL: m6Anet is DRACH-only — non-DRACH sites are invisible to the model; cDNA-Nanopore data CANNOT be used (PCR erases the signal); RNA002 and RNA004 chemistry require different model versions.

- CLI: `dorado basecaller rna004_130bps_sup@v5.0.0 --emit-fastq pod5/` -- modern RNA004 basecalling
- CLI: `minimap2 -ax map-ont -uf -k14 --secondary=no transcriptome.fa reads.fastq` -- transcriptome alignment
- CLI: `nanopolish eventalign --reads reads.fastq --bam aligned.bam --genome transcriptome.fa --scale-events --signal-index` -- signal-to-event (m6Anet-required flag pair)
- CLI: `m6anet dataprep --eventalign eventalign.txt --out_dir m6anet_data/` -- feature extraction
- CLI: `m6anet inference --input_dir m6anet_data/ --out_dir m6anet_results/` -- per-site probability + mod_ratio

## The Single Most Important Modern Insight -- m6Anet scores per-site DRACH probabilities, but reliability depends on coverage and the model is DRACH-only

Hendra 2022 *Nat Methods* 19:1590 recommends a minimum 20-50 reads per site for stable per-site `probability_modified` estimates; sites with <20 reads at probability_modified > 0.9 are likely false positives driven by per-read noise. The per-site `mod_ratio` column (the fraction of reads called modified at the site) is the stoichiometry-aware signal — informative even when per-site `probability_modified` is borderline. CRITICAL: m6Anet only scores DRACH 5-mers (D=A/G/U, R=A/G, A=methylated, C=C, H=A/C/U); non-DRACH sites are invisible to the model regardless of methylation status. For non-DRACH discovery, switch to Dorado / Remora modification basecallers (genome-wide; pin the modification-model version explicitly) or to xPore / Nanocompore (which compare two conditions without an a-priori 5-mer restriction). cDNA-Nanopore data cannot be used for m6Anet because PCR amplification erases the modification signal — only direct-RNA-sequencing (ONT DRS, kit SQK-RNA002 or SQK-RNA004) preserves the modification signal in the ionic current. Cross-method comparison: m6Anet ~0.51 recall on the Pratanwanich/Hendra synthetic-mix benchmark at >=10% modification / >=10x coverage; Dorado RNA004 recall is markedly higher (reported ~0.9 at the same threshold in a 2024-2025 RNA004 benchmark) but with higher per-site FDR (~40% at low-prevalence sites) — the two figures are NOT from a head-to-head test set, so treat as bounds rather than a comparable pair. The systematic 10-tool benchmark Zhong 2023 *Nat Commun* 14:1906 documents broader precision-vs-recall tradeoffs. The right modern pipeline is Dorado native for first-pass discovery -> m6Anet (or CHEUI) for filtering high-confidence subset -> GLORI for orthogonal stoichiometry validation at named loci.

## Algorithmic Taxonomy

| Tool | Mechanism | RNA002 | RNA004 | Strength | Fails when |
|------|-----------|--------|--------|----------|------------|
| m6Anet (Hendra 2022 *Nat Methods* 19:1590) | Multiple-instance-learning NN over DRACH 5-mer ionic-current features from nanopolish eventalign | YES | YES (recent) | Best-in-class for DRACH m6A on RNA002; generalises across cell lines | DRACH-only; requires nanopolish eventalign upstream |
| CHEUI (Acera Mateos 2024 *Nat Commun* 15:3899) | Deep CNN on ionic-current signals predicting m6A AND m5C in single molecules | YES | Limited | Simultaneous m6A + m5C calling; single-molecule co-occurrence | Newer; smaller user base; m5C less validated than m6A |
| Nanocompore (Leger 2021 *Nat Commun* 12:7198) | 2-component GMM comparing current + dwell-time between two samples | YES | Limited | Generic (any modification with signal change); replicate-aware | Requires modification-free control sample (WT vs KO or IVT) |
| xPore (Pratanwanich 2021 *Nat Biotechnol* 39:1394) | Bayesian multi-sample GMM; estimates fraction-modified per site per sample | YES | Limited | No matched WT/KO needed; multi-sample design; replicate support | Per-site coverage threshold dropout |
| ELIGOS / Eligos2 (Jenjaroenpun 2021 *NAR* 49:e7) | Compares error profile between native dRNA and unmodified controls (IVT) | YES | Limited | General-purpose; validated on yeast rRNA (95% recall); both error AND signal modes | Requires IVT or cDNA-seq control; modification-type-blind |
| EpiNano (Liu H et al. 2019 *Nat Commun* 10:4079) | Basecalling-error features as classifier features for m6A | Albacore 2.1.7 only | NO | Historical relevance; the original error-feature approach | Tied to obsolete Albacore basecaller; do NOT use for new analyses |
| Tombo (`ont-tombo`; Stoiber 2017 bioRxiv) | Signal-level model comparison; de novo against canonical RNA model OR sample-compare | YES | Limited | First general modification-detection tool; multi-modification | Less actively maintained since 2020; modern Dorado / m6Anet preferred |
| Dorado native modification calls (ONT; RNA004 chemistry, 2024+) | ONT basecaller with built-in modification calling (m6A, Ψ, m5C, inosine) | Limited (ONT focus shifted to RNA004; verify against current Dorado release notes) | YES | Modification calls in same output as base calls; per-read probabilities; first-pass discovery | Per-site precision lower than m6Anet at borderline sites; FDR ~40% at low-prevalence sites; model-version-pinning mandatory |
| m6ABasecaller / m6Aiso / mAFiA / m6ATM (2023-2024) | Newer ML approaches integrated with Dorado / basecaller | Mixed | Mixed | Active development; some integrated with Dorado | Less benchmarked; verify model lineage before reporting |
| DRUMMER | Pipeline integrating multiple modification tools | YES | Limited | Convenience wrapper | Inherits each tool's limitations |
| DiffErr | Error-rate differential between samples | YES | Limited | Lightweight differential | Modification-type-blind |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| RNA002 chemistry, single-sample m6A discovery, DRACH context | m6Anet for first-pass; CHEUI for cross-check at high-confidence | EpiNano tied to obsolete basecaller; Tombo less maintained |
| RNA004 chemistry, first-pass screening | Dorado native modification calling -> m6Anet (or CHEUI) for filtering high-confidence subset | Single Dorado pass has ~40% FDR at low-prevalence sites |
| Two-condition comparison (WT vs KO) | xPore (Bayesian; no matched IVT needed) OR Nanocompore (with control sample); avoid single-condition tool with post-hoc differential | Single-tool per-condition then ad-hoc comparison loses statistical efficiency |
| Multi-condition (time course, multiple KOs) | xPore (multi-sample Bayesian); fall back to per-condition m6Anet + manual differential | Per-condition + ad-hoc comparison is less efficient |
| Need m6A AND m5C simultaneously | CHEUI (only published tool covering both) | Running m6Anet (m6A only) + separate m5C tool loses single-molecule co-occurrence |
| Non-DRACH site discovery | Dorado native modification calling (genome-wide) OR xPore / Nanocompore (5-mer-agnostic) | m6Anet / CHEUI are DRACH-only |
| Verified cDNA-Nanopore data (no DRS) | None of these tools apply — PCR erases modification signal; modification detection is impossible | Running m6Anet on cDNA returns near-zero modification probabilities everywhere |
| Cross-laboratory reproducibility | Pin: basecaller version, Dorado modification model version, m6Anet model SHA, nanopolish version, reference transcriptome version (GENCODE / Ensembl release) | Unpinned versions break cross-batch comparability |
| Absolute stoichiometry needed | NOT direct-RNA alone -- supplement with GLORI on a subset of sites for calibration | Per-read modification rate is a stoichiometry estimate but has ~5-15% per-read error; GLORI is the gold standard |
| Wanting to validate MeRIP peaks orthogonally | m6Anet at the MeRIP peak positions in DRACH context; report concordance | Cross-method validation by independent technologies is the strongest evidence |

Methodology evolves; before any high-stakes m6Anet / direct-RNA analysis, web-search "m6anet v2 release notes", "Dorado modification model release notes", "RNA004 m6A benchmark 2024" for current best practice.

## Full m6Anet Pipeline (RNA002 / RNA004 chemistry)

**Goal:** Take POD5 / FAST5 signal data through basecalling, transcriptome alignment, nanopolish event alignment, and m6anet inference to produce per-site m6A modification probabilities at DRACH 5-mers.

**Approach:** Basecall with Dorado (modern) or Guppy (legacy); align to TRANSCRIPTOME (not genome) with minimap2 `-ax map-ont -uf -k14 --secondary=no`; sort and index BAM; run nanopolish eventalign with the m6Anet-required `--scale-events --signal-index` flag pair (plus `--summary` / `--threads` for housekeeping); m6anet dataprep extracts features; m6anet inference produces per-site `probability_modified` and `mod_ratio`.

```bash
# Step 1: Basecalling (RNA004 chemistry example).
dorado download --model rna004_130bps_sup@v5.0.0

dorado basecaller \
    rna004_130bps_sup@v5.0.0 \
    pod5/ \
    --emit-fastq \
    > reads.fastq

# Step 2: Transcriptome alignment (NOT genome).
minimap2 \
    -ax map-ont \
    -uf \
    -k14 \
    --secondary=no \
    -t 12 \
    refs/transcriptome.fa \
    reads.fastq | \
samtools sort -@ 8 -o aligned.bam -

samtools index aligned.bam

# Step 3: Nanopolish eventalign with m6Anet-required --scale-events --signal-index pair.
# nanopolish uses --genome for both genome AND transcriptome FASTA (confusing nomenclature; the flag accepts either).
nanopolish index \
    -d pod5/ \
    reads.fastq

nanopolish eventalign \
    --reads reads.fastq \
    --bam aligned.bam \
    --genome refs/transcriptome.fa \
    --scale-events \
    --signal-index \
    --threads 12 \
    --summary nanopolish_summary.tsv \
    > eventalign.txt

# Step 4: m6Anet dataprep + inference. Output includes mod_ratio (per-site stoichiometry) as well as probability_modified.
m6anet dataprep \
    --eventalign eventalign.txt \
    --out_dir m6anet_data/ \
    --n_processes 8

m6anet inference \
    --input_dir m6anet_data/ \
    --out_dir m6anet_results/ \
    --n_processes 4 \
    --num_iterations 1000

ls m6anet_results/
```

The m6Anet-required nanopolish flags are `--scale-events --signal-index`; m6Anet was trained on this feature set. The additional `--samples --print-read-names` flags are sometimes needed by downstream tools (yanocomp, f5c-pipeline interop) but are NOT required by m6Anet itself. Transcriptome alignment with `-uf` forces the forward-strand interpretation (DRS is directional); `-k14` is the recommended k-mer for ONT DRS. Use `-ax map-ont` against a TRANSCRIPTOME reference (the reference is already spliced); switch to `-ax splice -uf -k14` only when aligning DRS reads to a GENOME reference (the lh3 cookbook pattern for SIRV / spike-in genomes).

## Filtering and Interpreting m6Anet Results

**Goal:** Apply minimum-coverage and probability thresholds to `data.site_proba.csv` to retain high-confidence m6A site calls; report `mod_ratio` (per-site stoichiometry) alongside `probability_modified` (per-site model posterior).

**Approach:** Read the CSV; filter by `n_reads >= 20` (conservative) or `>= 50` (stringent) AND `probability_modified >= 0.9` (high-precision threshold per Hendra 2022); inspect `mod_ratio` as the per-site stoichiometry estimate. The full column set is `transcript_id, transcript_position, n_reads, probability_modified, kmer, mod_ratio`.

```python
import pandas as pd

PROBABILITY_THRESHOLD = 0.9
MIN_COVERAGE = 20

sites = pd.read_csv('m6anet_results/data.site_proba.csv')

print(sites.columns.tolist())
print(f'Total DRACH sites tested: {len(sites)}')

filtered = sites[(sites['n_reads'] >= MIN_COVERAGE) &
                 (sites['probability_modified'] >= PROBABILITY_THRESHOLD)]

print(f'High-confidence m6A sites (n_reads >= {MIN_COVERAGE}, prob >= {PROBABILITY_THRESHOLD}): {len(filtered)}')

print(f'mod_ratio summary among high-confidence: min={filtered["mod_ratio"].min():.2f} median={filtered["mod_ratio"].median():.2f} max={filtered["mod_ratio"].max():.2f}')

per_transcript = (filtered
    .groupby('transcript_id')
    .agg(n_high_conf_sites=('transcript_position', 'count'),
         mean_mod_ratio=('mod_ratio', 'mean'),
         total_coverage=('n_reads', 'sum'))
    .reset_index()
    .sort_values('n_high_conf_sites', ascending=False))

per_transcript.head(20).to_csv('m6anet_results/top_modified_transcripts.tsv', sep='\t', index=False)
filtered.to_csv('m6anet_results/high_confidence_sites.tsv', sep='\t', index=False)
```

`probability_modified` is the per-site model posterior (multiple-instance-learning aggregation over reads); `mod_ratio` is the per-site stoichiometry (fraction of reads called modified above the model's internal per-read threshold). For stoichiometry claims, `mod_ratio` is the right column. For genuine per-read output (probability per individual read), enable `--read_proba_threshold <T>` during inference (NOT `--per_read_proba_threshold`); per-model defaults are very small (`0.033379376` for HCT116_RNA002), so 0.5 is an interpretation threshold on the output column, not the CLI default.

## xPore Two-Condition Differential

**Goal:** Compare m6A modification rate at each DRACH site between two conditions using a Bayesian multi-sample GMM; no matched IVT control needed.

**Approach:** Run nanopolish eventalign + xpore dataprep separately for each condition; build a YAML config listing conditions and per-condition runs; xpore diffmod tests for differential modification per site.

```bash
xpore dataprep \
    --eventalign eventalign_ctrl.txt \
    --out_dir xpore_ctrl/

xpore dataprep \
    --eventalign eventalign_treat.txt \
    --out_dir xpore_treat/

cat > xpore_config.yaml << 'EOF'
data:
  ctrl:
    rep1: xpore_ctrl/dataprep
  treat:
    rep1: xpore_treat/dataprep
out: xpore_diff_output/
EOF

xpore diffmod \
    --config xpore_config.yaml \
    --n_processes 8

head xpore_diff_output/diffmod.table
```

xPore reports `diff_mod_rate` (per-site rate difference) AND p-value AND posterior probability. Filter for high-confidence differential by combining `diff_mod_rate >= 0.1` (10 percentage points) AND `pval < 0.05`.

## Nanocompore Comparative Modification Detection

**Goal:** Compare current intensity and dwell time per position between WT and KO (or treated vs control) to detect modification differences via Gaussian mixture comparison.

**Approach:** nanopolish eventalign per condition; Nanocompore CLI compares the two; outputs per-position differential GMM logit + KS + MWU statistics.

```bash
nanocompore eventalign_collapse \
    --input eventalign_ctrl.txt \
    --output_dir nanocompore_ctrl/

nanocompore eventalign_collapse \
    --input eventalign_treat.txt \
    --output_dir nanocompore_treat/

nanocompore sampcomp \
    --file_list1 nanocompore_ctrl/out_eventalign_collapse.tsv \
    --file_list2 nanocompore_treat/out_eventalign_collapse.tsv \
    --label1 ctrl \
    --label2 treat \
    --fasta refs/transcriptome.fa \
    --outpath nanocompore_diff_output/ \
    --nthreads 8

head nanocompore_diff_output/outSampComp_results.tsv
```

Nanocompore reports `GMM_logit_pvalue`, `KS_dwell_pvalue`, `KS_intensity_pvalue`, `MW_dwell_pvalue`, `MW_intensity_pvalue`. The GMM_logit test is the primary modification signal; use as a starting filter, then inspect dwell / intensity for direction. CRITICAL: Nanocompore is modification-type-AGNOSTIC — it detects signal differences which may be any of m6A, m5C, Ψ, or other modifications. Type assignment requires orthogonal information.

## Per-Method Failure Modes

### cDNA-Nanopore data fed to m6Anet

**Trigger:** Sequencing run used cDNA-Nanopore (PCR-amplified) rather than DRS (direct RNA); m6Anet pipeline attempted on the output.

**Mechanism:** Modification detection from nanopore signals (m6Anet, xPore, Nanocompore, ELIGOS, Dorado native) requires the RAW ionic-current signal from direct-RNA sequencing. cDNA-Nanopore amplifies the cDNA via PCR, which erases all modification signal because PCR uses canonical bases. Only DRS protocols (kit SQK-RNA002 or SQK-RNA004; no PCR; RNA-native sequencing with the DNA RT-adapter) preserve the modification signal.

**Symptom:** `nanopolish eventalign` runs successfully but produces uniform-looking event distributions; m6anet inference returns near-zero `probability_modified` everywhere; no DRACH enrichment in flagged sites.

**Fix:** Confirm DRS protocol BEFORE running any modification-detection pipeline. Check sequencing-core run report for kit (SQK-RNA002 / SQK-RNA004) vs cDNA kit (SQK-PCS, SQK-LSK). If cDNA, modification detection is impossible; defer to MeRIP or chemistry-based methods.

### m6Anet run with low-coverage sites at high probability

**Trigger:** Reporting m6Anet `probability_modified > 0.9` at sites with `n_reads < 10`.

**Mechanism:** m6Anet's per-site probability is computed via multiple-instance aggregation over reads. With low coverage, per-read predictions are weakly aggregated and per-site probability is noisy. Hendra 2022 recommends 20-50 reads per site minimum for stable calls.

**Symptom:** Many "high-confidence" m6A sites at sparsely-covered transcripts; per-site probability scatter at low coverage; sites do not validate orthogonally.

**Fix:** Apply minimum-coverage filter (`n_reads >= 20` conservative; `>= 50` stringent) BEFORE the probability filter. For high-stakes sites at low coverage, report both per-site probability AND per-read modification rate; orthogonally validate with GLORI / SAC-seq.

### Genome alignment instead of transcriptome alignment

**Trigger:** minimap2 invoked with `-ax splice -uf` (genome splice-aware) before m6Anet pipeline.

**Mechanism:** m6Anet expects reads aligned to a TRANSCRIPTOME FASTA (per-transcript coordinates); nanopolish eventalign signal alignment is per-transcript. Genome-aligned reads with splice junctions break the per-transcript signal interpretation.

**Symptom:** m6anet dataprep fails with chromosome / transcript ID errors; OR runs but produces few site calls; OR per-site probabilities are nonsensical.

**Fix:** Use `minimap2 -ax map-ont -uf -k14 --secondary=no` against a TRANSCRIPTOME FASTA for the m6Anet pipeline. For downstream genome-coordinate visualisation, convert site predictions back to genome coordinates via the transcript-to-genome mapping (use the GTF + a custom script or `samtools` lifting).

### Nanopolish eventalign without m6Anet-required flags

**Trigger:** `nanopolish eventalign` invoked without `--scale-events --signal-index`.

**Mechanism:** m6Anet was trained on the specific eventalign output format produced by `--scale-events --signal-index`. Without these, the dataprep step fails or produces feature vectors that don't match the model's expected input. The additional `--samples --print-read-names` flags are commonly added because other downstream tools (yanocomp, f5c interop) need them, but m6Anet itself does NOT require them.

**Symptom:** m6anet dataprep fails with parse errors; OR runs but produces empty feature files; OR inference returns no calls.

**Fix:** Always pass `--scale-events --signal-index` for m6Anet; add `--samples --print-read-names` only if downstream tools need them. Verify against m6anet quickstart for the installed version.

### f5c eventalign substituted for nanopolish without validation

**Trigger:** GPU-accelerated f5c eventalign used in place of nanopolish for speed; m6Anet model run on f5c output.

**Mechanism:** f5c is a re-implementation of nanopolish; eventalign output is numerically very close but NOT bit-identical, particularly in `event_level_mean` and `model_kmer` columns under certain edge cases. m6Anet was trained on nanopolish output; f5c is empirically usable but not officially supported and can shift per-site probabilities by 5-10% in extreme cases.

**Fix:** Use nanopolish unless GPU acceleration is required AND the user has validated f5c output equivalence on a control dataset. Document the eventalign source in any published analysis.

### Dorado modification model version not pinned

**Trigger:** SKILL.md / pipeline uses Dorado for m6A modification calling without specifying the model version.

**Mechanism:** Dorado modification models are versioned independently of the basecaller (`m6A_DRACH@v0.1`, `m6A_DRACH@v1.0`, etc.). Calls from different model versions are NOT directly comparable; model retraining shifts per-site probabilities.

**Fix:** Pin the exact modification model SHA / version in every CLI invocation; record in pipeline metadata. For multi-batch projects, re-call older batches with the new model when upgrading.

### Reference transcriptome version drift

**Trigger:** m6Anet results in batch 1 computed against GENCODE v44 transcripts; batch 2 against GENCODE v45; direct transcript-coordinate comparison.

**Mechanism:** Transcript IDs include version suffixes (`ENST00000123456.1` vs `.2`); some transcripts change between releases. Site coordinates relative to transcript start may shift.

**Fix:** Pin GENCODE / Ensembl transcript release for ALL m6Anet runs in a project. Lift over results when upgrading.

### Per-read modification rate confused with per-site probability

**Trigger:** Reporting m6Anet `probability_modified` as "the fraction of molecules modified at this site".

**Mechanism:** `probability_modified` is the model's posterior probability that the site has any modified reads (multiple-instance learning aggregation). It is NOT the per-read modification rate. Per-read modification rate is the fraction of reads at the site whose per-read posterior exceeds 0.5 (or another threshold).

**Symptom:** Stoichiometry claims based on per-site probability; numbers don't match orthogonal GLORI / SAC-seq stoichiometry estimates.

**Fix:** For stoichiometry, compute per-read modification rate from per-read output; cross-validate against GLORI at named sites.

### RNA002 model run on RNA004 chemistry data

**Trigger:** Pipeline uses an m6Anet model trained on RNA002 to analyse RNA004 chemistry runs.

**Mechanism:** RNA002 and RNA004 nanopore chemistry produce different signal characteristics — different motor protein, different translocation kinetics. Models trained on RNA002 do NOT transfer to RNA004 (and vice versa).

**Fix:** Pin RNA chemistry version per project; use the m6Anet model trained on the matching chemistry. m6Anet v2+ has explicit RNA004 model support; verify against `m6anet --help`.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| m6Anet high-confidence; GLORI does not call | Site in low-stoichiometry tail; GLORI conservatism | Trust GLORI for absolute; m6Anet probability + per-read rate informative for prevalence |
| m6Anet calls; MeRIP peak in same gene but different position | MeRIP fragment-level vs m6Anet single-base; methylation position offset within peak | Map m6Anet site against MeRIP peak coordinates with tolerance (~100 nt); concordance at gene level |
| Dorado calls; m6Anet does not at same DRACH | Dorado broader / less conservative; m6Anet model trained more strictly | Use m6Anet for high-precision; Dorado for high-recall first-pass |
| xPore differential strong; Nanocompore weak | xPore Bayesian multi-sample; Nanocompore pairwise GMM less powered | Trust xPore for multi-sample designs |
| Per-site probability high but per-read rate low | Few reads contributed strongly; rest were near-threshold | Suspect over-fit per-site at low coverage; verify with higher-coverage replicate |
| All ELIGOS / EpiNano sites have low signal | Older tools tied to obsolete basecallers; signal interpretation broken | Switch to modern tools (m6Anet, Dorado native) |
| CHEUI m6A and m5C overlap at same site | Co-occurring modifications OR cross-talk in model | CHEUI is the only published tool that handles co-occurrence; interpret as biological signal if validated |
| Tombo de novo and sample-compare disagree | Different baseline assumptions | Prefer sample-compare with KO / IVT control |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| m6Anet per-site `probability_modified` (high confidence) | >= 0.9 | Hendra 2022 *Nat Methods* 19:1590 default; high precision, lower recall |
| m6Anet per-site `probability_modified` (discovery) | >= 0.7 | Hendra 2022 lower threshold for less-stringent screens |
| m6Anet per-site minimum coverage | >= 20 reads (conservative); >= 50 (stringent) | Hendra 2022 recommendation for stable per-site estimates |
| Per-read modification probability interpretation threshold | >= 0.5 on output column (NOT CLI default) | m6Anet per-read posterior column convention; the CLI flag is `--read_proba_threshold` and per-model defaults are very small (~0.003-0.03) |
| xPore differential modification rate | >= 0.1 (10 percentage points) | xPore convention for meaningful biological effect |
| Nanocompore GMM_logit p-value | < 0.05 | Standard convention; primary statistic |
| ELIGOS error-rate ratio threshold | >= 2 (modified vs control) | Jenjaroenpun 2021 convention |
| Dorado native modification probability (per-read) | >= 0.5 (default) | ONT convention; threshold tuned per model version |
| minimap2 k-mer for ONT DRS | 14 | ONT recommendation for direct-RNA |
| minimap2 strand flag for DRS | `-uf` (forward strand only) | DRS is directional; reverse-strand alignment is incorrect |
| Per-tool benchmark recall (>=10% mod, >=10x cov, GLORI ground truth) | m6Anet ~0.51 (Pratanwanich/Hendra synthetic-mix); Dorado RNA004 ~0.9 (reported in 2024-2025 RNA004 benchmarks; NOT head-to-head) | Hendra 2022 *Nat Methods* 19:1590; broader 10-tool tradeoff context in Zhong 2023 *Nat Commun* 14:1906 |
| Cross-method concordance (m6Anet + GLORI + SAC-seq) | ~70-85% at high-confidence sites | Approximate from per-method validation tables in Liu C 2023 *Nat Biotechnol* 41:355 and Hu L 2022 *Nat Biotechnol* 40:1210; verify against current cross-method benchmark |
| Common-core m6A sites (HEK293T cross-method) | ~6,000-15,000 | Intersection of orthogonal methods |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `m6anet-dataprep` not found | v1 hyphenated CLI not in installed v2 | Use `m6anet dataprep` (subcommand syntax) |
| nanopolish eventalign errors with "no event annotation" | FAST5 / POD5 not indexed; OR `--scale-events --signal-index` missing | Run `nanopolish index -d pod5/ reads.fastq` first; include the m6Anet-required flag pair |
| m6anet inference returns empty CSV | dataprep output empty; OR all sites filtered by default coverage threshold | Check dataprep log; verify feature TSV in dataprep output dir |
| Reads marked as cDNA in run summary | Library prep used cDNA kit (SQK-PCS / SQK-LSK), not DRS | Modification detection impossible; switch to MeRIP / chemistry-based methods |
| `minimap2 -ax splice` instead of `map-ont` | Genome-splice flag used for transcriptome alignment | Switch to `-ax map-ont -uf -k14` |
| `nanopolish eventalign` extremely slow | Single-threaded by default | Pass `--threads N` |
| f5c output produces shifted m6Anet probabilities | f5c is numerically near-identical but not bit-identical to nanopolish | Use nanopolish unless GPU acceleration required |
| Per-site probability high at 5 reads | Low coverage; per-site noise | Filter by `n_reads >= 20` before interpretation |
| xPore diffmod takes hours per condition | Single-threaded default | Pass `--n_processes 8` |
| Dorado m6A model not found | Model not downloaded | `dorado download --model m6A_DRACH@v1` |
| Transcriptome IDs missing from m6anet output | Reference mismatch between minimap2 input and m6anet config | Use same transcriptome FASTA throughout pipeline |
| Modification probabilities differ between Dorado releases | Model versioned independently from basecaller | Pin model version explicitly |
| Per-read CSV not generated | `--read_proba_threshold` not set (NOT `--per_read_proba_threshold`) | Add `--read_proba_threshold <T>` to `m6anet inference`; T is model-specific (default very small, e.g., 0.033 for HCT116_RNA002) |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why m6Anet over Dorado native?" | m6Anet for high-precision DRACH calling; Dorado for high-recall first-pass; both reported when feasible |
| "How were non-DRACH sites handled?" | m6Anet is DRACH-only by design; non-DRACH discovery via Dorado / xPore / Nanocompore; cited limitation |
| "What's the per-site coverage threshold?" | `n_reads >= 20` conservative per Hendra 2022; `>= 50` for stringent calls |
| "Was cross-validation against orthogonal methods done?" | High-stakes sites cross-checked against GLORI / m6A-SAC-seq / MeRIP peaks; cross-method concordance reported |
| "Why RNA002 vs RNA004?" | Chemistry version pinned per project; m6Anet model selected to match (RNA002 model on RNA002 data, RNA004 on RNA004) |
| "Was the Dorado modification model version pinned?" | Yes — version recorded in pipeline metadata; rerun batches with new model when upgrading |
| "Per-site probability vs mod_ratio?" | `probability_modified` is the per-site model posterior; `mod_ratio` is the per-site stoichiometry (fraction of reads called modified) — both reported |
| "Is mod_ratio the same as absolute stoichiometry?" | Approximately; per-read calls have ~5-15% error; for absolute stoichiometry at named loci, cross-validate with GLORI |
| "Why not just use Tombo / EpiNano?" | Tombo less maintained since 2020; EpiNano tied to obsolete Albacore basecaller; m6Anet / Dorado / CHEUI are modern |
| "How were f5c vs nanopolish handled?" | Used nanopolish (m6Anet-trained source); did not substitute f5c |

## References

- Hendra C, Pratanwanich PN, Wan YK, Goh WSS, Thiery A, Göke J (2022) Detection of m6A from direct RNA sequencing using a multiple instance learning framework. *Nat Methods* 19(12):1590-1598. doi:10.1038/s41592-022-01666-1
- Pratanwanich PN, Yao F, Chen Y et al (2021) Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. *Nat Biotechnol* 39(11):1394-1402. doi:10.1038/s41587-021-00949-w
- Leger A, Amaral PP, Pandolfini L et al (2021) RNA modifications detection by comparative Nanopore direct RNA sequencing. *Nat Commun* 12(1):7198. doi:10.1038/s41467-021-27393-3
- Jenjaroenpun P, Wongsurawat T, Wadley TD et al (2021) Decoding the epitranscriptional landscape from native RNA sequences. *Nucleic Acids Res* 49(2):e7. doi:10.1093/nar/gkaa620
- Liu H, Begik O, Lucas MC et al (2019) Accurate detection of m6A RNA modifications in native RNA sequences. *Nat Commun* 10(1):4079. doi:10.1038/s41467-019-11713-9
- Acera Mateos P, Sethi AJ, Ravindran A et al (2024) Prediction of m6A and m5C at single-molecule resolution reveals a transcriptome-wide co-occurrence of RNA modifications. *Nat Commun* 15:3899. doi:10.1038/s41467-024-47953-7
- Stoiber MH, Quick J, Egan R et al (2017) De novo identification of DNA modifications enabled by genome-guided nanopore signal processing. *bioRxiv* 094672. doi:10.1101/094672
- Zhong ZD, Xie YY, Chen HX et al (2023) Systematic comparison of tools used for m6A mapping from nanopore direct RNA sequencing. *Nat Commun* 14:1906. doi:10.1038/s41467-023-37596-5
- Liu C, Sun H, Yi Y et al (2023) Absolute quantification of single-base m6A methylation in the mammalian transcriptome using GLORI. *Nat Biotechnol* 41(3):355-366. doi:10.1038/s41587-022-01487-9
- Hu L, Liu S, Peng Y et al (2022) m6A RNA modifications are measured at single-base resolution across the mammalian transcriptome. *Nat Biotechnol* 40(8):1210-1219. doi:10.1038/s41587-022-01243-z
- Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094-3100. doi:10.1093/bioinformatics/bty191
- Loman NJ, Quick J, Simpson JT (2015) A complete bacterial genome assembled de novo using only nanopore sequencing data. *Nat Methods* 12(8):733-735. doi:10.1038/nmeth.3444
- Dominissini D, Moshitch-Moshkovitz S, Schwartz S et al (2012) Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. *Nature* 485(7397):201-206. doi:10.1038/nature11112
- Linder B, Grozhik AV, Olarerin-George AO, Meydan C, Mason CE, Jaffrey SR (2015) Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. *Nat Methods* 12(8):767-772. doi:10.1038/nmeth.3453

## Related Skills

- merip-preprocessing - Genome-aligned MeRIP for cross-validation against direct-RNA calls
- m6a-peak-calling - MeRIP fragment-level peaks for orthogonal validation comparison
- m6a-differential - Per-site direct-RNA modification rate comparable to MeRIP differential at named loci
- modification-visualization - Metagene and browser-track rendering of m6Anet site calls
- long-read-sequencing/basecalling - Upstream Dorado / Guppy basecalling fundamentals
- long-read-sequencing/long-read-alignment - General minimap2 alignment patterns
- long-read-sequencing/long-read-qc - Direct RNA QC (yield, length distribution, basecall accuracy)
- long-read-sequencing/nanopore-methylation - DNA methylation calling sibling (different chemistry, related framework)
- read-alignment/star-alignment - Splice-aware alignment for cross-method comparison context
- variant-calling/vcf-basics - General per-site variant call framework (analogue)
- rna-quantification/featurecounts-counting - For cross-method site-vs-MeRIP-peak comparison
