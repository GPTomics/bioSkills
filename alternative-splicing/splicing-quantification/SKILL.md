---
name: bio-splicing-quantification
description: Quantifies alternative splicing as PSI (percent spliced in) from RNA-seq using rMATS-turbo (BAM-based event), SUPPA2 (TPM-based event), MAJIQ V3 (LSV-based Bayesian), leafcutter (annotation-free intron clusters), VAST-TOOLS (cross-species with microexon support), Shiba (junction-imbalance-corrected, 2025 SOTA at low coverage), or IRFinder-S (intron retention coverage-aware). Distinguishes the five canonical event classes (SE, A5SS, A3SS, MXE, RI), special classes (microexons, exitrons, AFE/ALE), intron retention subtypes (canonical RI vs detained introns), and applies effective-length normalization. Use when measuring splice-site usage or isoform inclusion ratios from short-read RNA-seq.
tool_type: mixed
primary_tool: rMATS-turbo
---

## Version Compatibility

Reference examples tested with: rMATS-turbo 4.3+, SUPPA2 2.4+, leafcutter 0.2.9+, MAJIQ 3.0+, IRFinder-S 2.0+, kallisto 0.50+, Salmon 1.10+, pandas 2.2+, STAR 2.7.11+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Splicing Quantification

Quantify alternative splicing events as PSI (percent spliced in) from RNA-seq. PSI = inclusion read evidence / (inclusion + skipping read evidence), normalized for differential mapping opportunity between isoforms. The choice of *quantification unit* (event, intron cluster, LSV, transcript) determines which biological questions can be answered and which failure modes apply.

## Algorithmic Taxonomy

| Family | Unit | Reference tools | Fails when |
|--------|------|-----------------|------------|
| Event-based | Pre-defined SE/A5SS/A3SS/MXE/RI events from annotation | rMATS-turbo, SUPPA2, VAST-TOOLS | Event isn't in annotation; complex multi-junction events split arbitrarily; AFE/ALE confounded with splicing |
| LSV-based | Local Splice Variations at single source/target nodes | MAJIQ V3 | Memory-constrained environments; cohorts smaller than ~3 reps; non-academic users (license) |
| Junction-cluster | Annotation-free intron clusters by shared splice sites | leafcutter, leafcutter2 | Undersampled clusters lose power; topology biologically uninterpretable for novel events |
| Splice-graph | Graph nodes (non-overlapping exonic regions) | Whippet, Shiba | Whippet maintenance status uncertain since ~2022; complex multi-exon graphs |
| Coverage-aware (IR) | Intron body coverage + flanking junctions | IRFinder-S, S-IRFindeR, iREAD | Confounded by overlapping exons, repeats, low mappability regions |
| Isoform-based | Transcript abundance via EM | Salmon/kallisto + tximport | Salmon EM uncertainty propagates; many similar isoforms (TTN, MAPT) become indistinguishable |

The agent's first decision is **which family** the question requires, not which tool. Switching family is the response to a tool failing within a family — switching tool within a family rarely fixes systematic blind spots.

## Event Taxonomy (Beyond Standard SE/A5SS/A3SS/MXE/RI)

| Class | Code | Biology | Detection caveat |
|-------|------|---------|-------------------|
| Skipped exon (cassette) | SE | Default cassette-exon AS | Most common; well-handled by all tools |
| Alternative 5' splice site | A5SS | Alternative donor; intron 5' end varies | Sign convention tool-specific (see below) |
| Alternative 3' splice site | A3SS | Alternative acceptor; intron 3' end varies | Sensitive to BPS cancer mutations (SF3B1) — cryptic 3'ss ~10-30nt upstream |
| Mutually exclusive exons | MXE | Two exons paired, one included | Tool implementations vary; verify which "form 1" is yours |
| Retained intron | RI | Whole intron retained in mature mRNA | Junction-only quant systematically underdetects; needs IRFinder-S |
| **Microexon** | (sub-SE) | 3-27 nt; neural-enriched, SRRM4-regulated | Missed by default aligner anchor lengths (>=20-30 nt); needs VAST-TOOLS, MicroExonator, or long-read |
| **Exitron** | n/a | Intronic region within an annotated CDS exon | Mis-classified as A5SS/A3SS by most tools; use ExitronFinder, ScanExitron |
| **Alternative first exon (AFE)** | AFE | Alternative TSS / promoter use | Promoter-driven, NOT spliceosomal; confirm with FANTOM CAGE before reporting as splicing |
| **Alternative last exon (ALE)** | ALE | Alternative cleavage/polyadenylation | APA-driven, NOT spliceosomal; confirm with 3'-end-seq |
| **Detained intron (DI)** | (RI subtype) | Nuclear-retained on mature mRNA, regulated by Clk kinases | Distinct from cytoplasmic NMD-targeted RI (Boutz 2015 *Genes Dev*); requires fractionation to confirm |
| **Recursive splicing** | (intron subtype) | Long introns >50kb spliced via internal ratchet points | Sibley 2015 *Nature*; only detectable with long-read or nascent RNA-seq |

Tool-agnostic taxonomy reference: Wang 2008 *Nature*; Vaquero-Garcia 2016 *eLife* (LSV); Tapial 2017 *Genome Res* (VastDB).

## Tool Selection Matrix

| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| rMATS-turbo | Standard SE/A5SS/A3SS/MXE/RI in annotated organism, n>=3 | BAM + GTF | Fast, well-calibrated at n>=3, novel SS support | Junction read imbalance; novel multi-junction events; underdetects RI and microexons |
| SUPPA2 | Quick PSI from existing TPM, pilot analysis | Salmon/kallisto TPM + GTF | No alignment; fastest | Annotation-bound; high FDR (15-30%) at n<=2 vs n<=2 |
| MAJIQ V3 | Complex events, heterogeneous cohorts (HET) | BAM + GFF3 | Bayesian posterior PSI, complete LSV semantics | High memory (~50+ GB on cohorts); academic license; complex LSV interpretation needs care |
| leafcutter | Novel junction discovery, sQTL, low-memory environments | BAM (regtools junctions) | Annotation-free; ~400 MB memory; SOTA for unannotated organisms | Sensitive to read depth; cluster topology arbitrary for complex multi-junction events |
| Shiba (2025) | Low-coverage / few-replicate designs; junction imbalance correction | BAM + GTF | Best calibration in own benchmark at n=2 vs n=2 | New (2025); limited community calibration |
| VAST-TOOLS | Cross-species comparative AS, microexons | FASTQ | VastDB orthology, ExOrthist co-tool | Limited to species in VastDB |
| Whippet | Laptop-scale exploratory | FASTQ | Fast splice-graph-based PSI | Reduced active development since ~2022; underperforms on complex topologies |
| IRFinder-S | Intron retention specifically | FASTQ | Coverage + junction integration; CNN-based artifact filtering | IR-only; not for cassette events |
| S-IRFindeR | Replicate-stable IR ratio | BAM | Stable IR ratio metric | IR-only; less integrated than IRFinder-S |

Methodology evolves; verify benchmarks (Olofsson 2023 *Brief Bioinform*; Kubota 2025 *NAR*; Tran 2025 *WIREs RNA*) and tool docs before committing. Default 2026 recommendation: run rMATS-turbo + leafcutter and reconcile; add MAJIQ V3 for complex events / heterogeneous cohorts; switch to Shiba for n=2 vs n=2.

## PSI Definition and Effective Length Normalization

For a cassette exon, naive PSI ignores that the inclusion isoform contains more positions where a junction read can map than the skipping isoform. rMATS reports `IncFormLen` (= 2*(read_length - anchor) for the two flanking junctions, plus exon body bases) and `SkipFormLen` (= read_length - anchor) and computes:

```
PSI = (IJC / IncFormLen) / (IJC / IncFormLen + SJC / SkipFormLen)
```

where IJC = inclusion junction counts, SJC = skipping junction counts. **Skipping this normalization biases PSI by ~10-30% depending on read length and exon size.** SUPPA2 derives PSI from transcript TPMs, which the upstream Salmon/kallisto already accounts for. leafcutter operates on intron usage proportions within a cluster (a different statistic).

For long-read data, every read carries full isoform identity — effective-length normalization becomes unnecessary because each read counts as one isoform.

## Sign Conventions for Alternative Splice Sites

| Tool | A5SS interpretation | A3SS interpretation | "Inclusion" direction |
|------|----------------------|----------------------|------------------------|
| rMATS | "long" form = donor downstream of alternative donor | "long" form = acceptor upstream of alternative acceptor | PSI > 0 = more long form |
| SUPPA2 | Same as rMATS (long = inclusion of additional exon body) | Same | PSI > 0 = more long form |
| VAST-TOOLS | Encoded in event ID (`_D1` vs `_D2`) | Encoded in event ID (`_A1` vs `_A2`) | Document the chosen reference |
| MAJIQ | Per-junction within LSV; explicit donor/acceptor naming in VOILA | Same | PSI per junction in the LSV |

Always record which alternative form ΔPSI > 0 corresponds to in publication-grade reporting. Confusion is the most common reviewer comment for AS papers.

## rMATS-turbo Workflow

**Goal:** Quantify SE/A5SS/A3SS/MXE/RI events from BAMs aligned with STAR 2-pass.

**Approach:** Group BAMs by condition, run rMATS with `--statoff` for quantification only, then parse JC.txt files for per-replicate PSI.

```bash
rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --libType fr-firststrand \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp \
    --novelSS \
    --statoff
```

Key flags: `--novelSS` discovers junctions absent from the GTF (recommended with STAR 2-pass output). `--variable-read-length` allows mixed read lengths in the cohort. `--libType fr-firststrand` matches Illumina TruSeq stranded; verify with RSeQC `infer_experiment.py`. `--statoff` is for quantification-only runs; omit for differential testing.

```python
import pandas as pd

se_jc = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')

inc_cols = [c for c in se_jc.columns if c.startswith('IncLevel')]
se_jc['mean_PSI'] = se_jc[inc_cols].mean(axis=1)

per_rep_inc = se_jc['IJC_SAMPLE_1'].str.split(',').apply(lambda x: list(map(int, x)))
per_rep_skip = se_jc['SJC_SAMPLE_1'].str.split(',').apply(lambda x: list(map(int, x)))
min_inc = per_rep_inc.apply(min)
min_skip = per_rep_skip.apply(min)

reliable = se_jc[(min_inc + min_skip) >= 20]
```

### JC vs JCEC files

`SE.MATS.JC.txt` uses **only junction-spanning reads**. `SE.MATS.JCEC.txt` adds **reads contained within the alternative exon body** as inclusion evidence.

- Prefer **JC** for clean cassette-exon analysis when reads spanning junctions are sufficient.
- Use **JCEC** when alternative exons are short (<50nt) and junction-spanning reads are scarce.
- Avoid **JCEC** when intron retention overlaps the alternative exon body — exon-body reads may come from retained introns, not inclusion isoform.

## SUPPA2 Workflow

**Goal:** Compute event PSI from transcript TPM without alignment; useful when Salmon/kallisto TPMs already exist.

**Approach:** Generate IOE event definitions from GTF, then aggregate TPMs of transcripts including/excluding each event.

```bash
suppa.py generateEvents -i annotation.gtf -o events -f ioe -e SE SS MX RI AF AL

for ev in SE A5 A3 MX RI; do
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e transcript_tpm.tsv -o psi_${ev}
done
```

SUPPA2 is annotation-bound: events absent from the GTF cannot be quantified. Whether an event is detected depends entirely on which transcripts the upstream Salmon/kallisto index contains. Use GENCODE comprehensive over basic when SUPPA2 detection sensitivity matters.

## MAJIQ V3 Workflow

**Goal:** Quantify LSVs with Bayesian posterior PSI distributions; ideal for complex multi-junction events that don't fit canonical event types.

**Approach:** Build a splice graph from BAMs + GFF3, compute per-junction coverage with bootstrap, then run `majiq psi` for posterior PSI per LSV.

```bash
majiq build annotation.gff3 -c settings.ini -j 8 -o build_output
majiq psi build_output/sample1.majiq build_output/sample2.majiq -j 4 -o psi_output -n condition_psi
voila view -p 5000 -j 8 build_output/splicegraph.zarr psi_output/condition_psi.psi.voila -o voila_output
```

MAJIQ V3 (Aicher, Slaff, Jewell, Barash *bioRxiv* 2024; public release 2025) replaced V2's SQLite splicegraph (`splicegraph.sql`) with **Zarr storage** (`splicegraph.zarr`); the `.sql` is deprecated. V3 is ~3.2x faster than V2 via xarray/zarr/Dask parallelization. LSV output includes posterior mean PSI plus the full posterior distribution; this enables threshold-based testing (e.g. P(|ΔPSI| > 0.2)) rather than frequentist p-values.

## leafcutter Junction Quantification

**Goal:** Detect junctions and intron clusters annotation-free for downstream cluster-level usage.

**Approach:** Extract junctions per BAM with regtools, write filenames into a list, then cluster introns sharing splice sites.

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done
ls *.junc > juncfiles.txt

python leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -o leafcutter \
    -m 50 \
    -l 500000
```

`-a 8` = 8nt anchor minimum (raise to 12 for stricter; lower to 6 for microexon-friendly). `-m 50` = minimum junction reads per cluster. `-l 500000` = max intron length (relevant for long brain-gene introns; raise for genes like DSCAM, ROBO2, ANK3).

## Per-Tool Failure Modes

### rMATS-turbo: Junction Read Imbalance

**Trigger:** A cassette exon's flanking exons have unequal read mapping opportunity (very short upstream exon, repeat-overlapping flanks, or low-mappability regions).

**Mechanism:** rMATS' binomial model treats inclusion vs skipping junctions as having equal mappability. When mappability differs between the two junction types, the PSI estimate is biased.

**Symptom:** "Significant" rMATS calls with no concordant change in leafcutter or MAJIQ at the same locus; ΔPSI direction inconsistent with sashimi-plot intuition.

**Fix:** Run Shiba (Kubota 2025 *NAR*) which corrects junction-imbalance, or filter rMATS hits requiring concordant detection in leafcutter.

### SUPPA2: Sparse Empirical Null at Low Replicate Count

**Trigger:** n=2 vs n=2 (or n=3 vs n=2) design with `--method empirical`.

**Mechanism:** SUPPA2's empirical null is constructed from between-replicate ΔPSI distributions binned by transcript expression. With few replicates, the binned null is sparse and conservative-looking but actually under-calibrated.

**Symptom:** Inflated FDR (15-30% in benchmarks); many "significant" hits don't replicate or validate.

**Fix:** Switch to leafcutter or Shiba for n<=3 designs; or use `--method classical` (Wilcoxon) for very low replicate count; reconcile against orthogonal tool.

### MAJIQ V3: Complex LSV Interpretation

**Trigger:** A gene with 4+ alternative splice sites at one node (e.g. one source, multiple acceptors).

**Mechanism:** A complete LSV at a single node lists all observed junctions; PSI is per-junction within the LSV, not "PSI of one event."

**Symptom:** Reporting "PSI of the gene" doesn't make sense; per-junction PSIs sum to 1 across the LSV but no single number represents the gene.

**Fix:** Use VOILA to visualize the LSV graph and identify which junction(s) shifted; for cassette-style reporting, derive equivalent PSI from sum of inclusion-junctions / total junctions in the LSV.

### leafcutter: Cluster Topology Arbitrariness

**Trigger:** A cluster has 4+ introns sharing splice sites with non-canonical topology (e.g. mixed cassette + alternative donor + IR).

**Mechanism:** leafcutter clusters introns by shared splice sites; complex topologies don't map onto SE/A5SS/A3SS taxonomy and cluster-level "ΔPSI" hides which intron drove the change.

**Symptom:** Significant cluster-level p-value but multiple introns showing different effect-size directions.

**Fix:** Inspect the cluster in leafviz; report per-intron effect sizes (`effect_sizes.txt`); for canonical event reporting, map to SE/A5SS/A3SS via flanking exon coordinates manually.

## Reconciliation: When rMATS and leafcutter Disagree

The two most common short-read tools answer slightly different questions: rMATS classifies on annotated event templates; leafcutter classifies on observed cluster usage. Disagreement is informative.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| rMATS sig, leafcutter not sig | rMATS junction read imbalance OR rMATS event hits an annotation that leafcutter clustered differently | Inspect locus in IGV; check Shiba |
| leafcutter sig, rMATS not sig | Novel junction not in rMATS annotation; rMATS `--novelSS` may have missed it | Check `--novelSS` was on; rerun if not |
| Both sig, opposite ΔPSI direction | Event class mismatch (e.g. rMATS calls SE positive but leafcutter sees A5SS shift in same cluster) | Manually map cluster topology to event class |
| Both sig, same direction | High-confidence call | Report; cross-validate with sashimi-plot |

**Operational rule:** for high-confidence reporting, require concordant detection in two tools from different algorithmic families (event-based + cluster-based, or LSV + isoform-based).

## Intron Retention: Canonical vs Detained vs Co-Transcriptional Unspliced

Three biologically distinct states all called "IR" by generic tools:

1. **Canonical RI (cytoplasmic, NMD-substrate often)**: mature polyadenylated mRNA carries the intron; usually PTC-bearing and NMD-targeted, sometimes encoding an alternative protein.
2. **Detained intron (DI)** (Boutz 2015 *Genes Dev*): nuclear-localized, mature transcripts retaining a specific intron; a regulated reservoir released into translation upon signaling.
3. **Co-transcriptional unspliced**: nascent pre-mRNA captured before splicing complete; not a regulated state.

**Library prep determines which state(s) you see:**
- Poly(A) selection: enriches (1), depletes (2)/(3)
- rRNA depletion (cytoplasmic): captures (1)
- rRNA depletion (whole cell or nuclear): captures all three

**To distinguish DI from canonical RI:** subcellular fractionation (nuclear vs cytoplasmic RNA-seq), or NMD inhibitor (cycloheximide, NMDi-14) treatment — canonical RI mRNA increases under NMD inhibition; DI does not.

```bash
IRFinder -m FullAuto -r REF/ -d ir_output sample.fastq.gz
```

IRFinder-S (Lorenzi 2021 *Genome Biol*) uses CNN-based filtering of true IR vs noise; current SOTA for IR analysis. iREAD and S-IRFindeR (Broseus & Ritchie 2020 *bioRxiv*) are alternatives.

## Microexon Detection

Microexons (3-27 nt, neural-enriched, SRRM4-regulated; Irimia 2014 *Cell*) are missed by default short-read aligners requiring 20-30 nt anchors. Options:

| Approach | Tool | Notes |
|----------|------|-------|
| Curated database lookup | VAST-TOOLS + VastDB | Cross-species, microexon-aware (Tapial 2017 *Genome Res*) |
| De novo discovery | MicroExonator (Parada 2021 *Genome Biol*) | Snakemake pipeline |
| Tune the upstream aligner | `STAR --alignSJoverhangMin 6 --alignSJDBoverhangMin 1 --outFilterMismatchNoverReadLmax 0.04` | rMATS itself cannot recover microexons that STAR didn't pass through; lower DB-junction overhang to 1 (trusts annotated microexon coords) and combine with strict mismatch filter. Typical AS pipelines use STAR 8/3 which is too strict for microexons |
| Long-read sequencing | PacBio Iso-Seq, ONT | Solves the problem entirely; reads span microexons fully |

For brain / neural tissue or autism-spectrum studies, **microexon analysis must be explicit** — default short-read pipelines underdetect them by ~70%.

## Quality Thresholds

| Metric | Threshold | Source / Rationale |
|--------|-----------|---------------------|
| Junction reads per replicate | >=10-20 (per-replicate minimum) | Empirical PSI variance becomes <0.05 above this; below, PSI becomes a coin flip |
| PSI dynamic range | mean PSI 0.05-0.95 | Outside is near-constitutive; rMATS, SUPPA2 default filters drop these |
| Missing values | <50% of samples | Higher missingness indicates low expression — re-test with subset |
| Read length | >=75nt PE preferred; >=100nt for microexons | 50nt SE biases toward shorter exons (Wang 2008 *Nature*) |
| Library | rRNA depletion for IR analysis; poly(A) acceptable for cassette | Sims 2014 *Genome Res*; poly(A) loses pre-mRNA |
| STAR 2-pass | Cohort-style preferred over per-sample basic | Veeneman 2016 *Bioinformatics*: >=94% novel junction recovery |
| MAJIQ minreads / minpos | --minreads 10 --minpos 3 | Default; lower for low-coverage |
| leafcutter -m | 50 reads per cluster | Higher for rare events; lower for sQTL discovery |
| Anchor length | >=8 nt for short-read | Below this, false-positive junctions dominate (CIGAR-N noise) |

## Decision Tree by Scenario

| Scenario | Recommended tool(s) | Why |
|----------|----------------------|-----|
| Standard cassette analysis, n>=3, GENCODE-annotated | rMATS-turbo + leafcutter (concordance) | Default workflow; complementary algorithmic families |
| Non-model organism, no GENCODE-grade annotation | leafcutter + de novo discovery | Annotation-free |
| Heterogeneous cohort, n>=10 vs n>=10 (clinical, GTEx-style) | MAJIQ V3 with HET module | HET designed for between-sample variability dominance |
| Low coverage / few replicates (n=2 vs n=2) | Shiba | Junction-imbalance correction; SOTA at low coverage in 2025 benchmarks |
| Cross-species comparative (vertebrate panel) | VAST-TOOLS + VastDB | Orthology-aware events; ExOrthist co-tool |
| TPM-only available (no BAMs) | SUPPA2 | Annotation-bound but fast |
| Microexon focus (neural / ASD) | VAST-TOOLS or MicroExonator | Default tools systematically miss microexons |
| Intron retention focus | IRFinder-S (rRNA-depleted library) | Coverage-aware; CNN artifact filter |
| Detained introns specifically | IRFinder-S + nuclear/cytoplasmic fractionation | Required to separate DI from cytoplasmic RI |
| Long reads available | rMATS-long, FLAIR, IsoQuant | Full-isoform resolution; see long-read-splicing |
| Single-cell (full-length plate) | MARVEL, BRIE2 | See single-cell-splicing |
| Single-cell (10X 3') | Likely don't attempt; consider Sierra for APA | 10X 3' chemistry insufficient for AS |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `error: GTF gene_id parsing` (rMATS) | rMATS expects GENCODE-style gene_id; some Ensembl GTFs use different attribute order | `gffread input.gff3 -T -o standardized.gtf` |
| `KeyError: 'IJC_SAMPLE_1'` (rMATS parsing) | Output column missing; sometimes occurs when --statoff combined with novel events on older versions | Update rMATS-turbo to >=4.3.x; re-run |
| `MAJIQ: too few reads at junction` | Default `--minreads 10 --minpos 3` filters out the locus | Lower thresholds for low-coverage data; document filtering |
| `leafcutter: dispersion estimation failed` | Cluster has all-zero counts in one group | Pre-filter clusters with `--min-samps-feature-prop 3` |
| `SUPPA2: empirical p computed on N=4 nulls` | Insufficient replicates for empirical mode | Switch `--method classical` (Wilcoxon) for very low replicate count |
| `regtools: invalid CIGAR` | Non-BAM-spec read in input | Filter with `samtools view -h -F 0x100 -F 0x800` (drop secondary/supplementary) |
| `STAR: too many SJs` (in pass 2) | Cohort SJ.out.tab too large | Filter to junctions seen in >=3 samples or with >=3 unique reads before merging |

## Output Interpretation

PSI ranges 0 to 1: 1 = always included, 0 = always skipped, 0.5 = balanced. Sign of `IncLevelDifference` matches `--b1 minus --b2` group order — always document which is which in publications.

**NMD direction matters:** an increase in PSI of a poison exon (PTC-introducing) **decreases** functional protein due to NMD. Always check whether the alternative form is PTC-bearing using ORF-aware annotation (IsoformSwitchAnalyzeR consequences, or manual stop-codon distance check vs last exon-exon junction).

**Disease signatures:**
- **SF3B1** mutations (MDS, CLL, uveal melanoma): cryptic 3'ss ~10-30 nt upstream of canonical (Darman 2015 *Cell Rep*). Look for clustered A3SS hits.
- **U2AF1** mutations (lung adeno, MDS): altered preferences at 3'ss -3 position; cassette-exon shifts.
- **TDP-43 loss** (ALS/FTD): de novo cryptic exons in UNC13A, STMN2, ATG4B (Brown 2022 *Nature*; Klim 2019 *Nat Neurosci*) — annotation-free tools required (leafcutter denovo).

## Common Pitfalls

- Treating AFE/ALE as splicing — these are typically promoter-driven (AFE) or APA-driven (ALE), not spliceosomal. Confirm with FANTOM CAGE or 3'-end-seq.
- Confusing detained introns with NMD-targeted RI — both call as "IR" but have opposite biological fates.
- Using poly(A) libraries for IR analysis — biases toward mature transcripts, depletes pre-mRNA.
- Single-end short reads — junction-spanning reads need >=8nt overhang on both sides; biases toward shorter exons.
- Quoting "PSI of the gene" from MAJIQ LSV output — only per-junction PSI within an LSV is meaningful.
- Skipping STAR 2-pass — loses ~14% of novel junctions; matters for any non-canonical organism or condition.
- Trusting rMATS calls without `--novelSS` when STAR 2-pass found new junctions — rMATS will only quantify pre-annotated events.

## Related Skills

- differential-splicing - Compare PSI between conditions; use the same upstream alignment but switch to with-stat tools
- splicing-qc - Run BEFORE quantification to verify library, depth, strandedness, alignment quality
- isoform-switching - DTU framework with NMD/ORF/domain consequences; complementary to event-level PSI
- sashimi-plots - Visualize specific events for QC and reporting; concordance check across tools
- splice-variant-prediction - SpliceAI/Pangolin for variant impact predictions to test against PSI changes
- long-read-splicing - Full-isoform PSI without anchor-length limits; preferred for microexons and complex isoforms
- read-alignment/star-alignment - STAR 2-pass cohort-style alignment is required upstream
- rna-quantification/alignment-free-quant - Salmon/kallisto TPM is required for SUPPA2

## References

- Wang et al 2008 *Nature* - AS event taxonomy
- Vaquero-Garcia et al 2016 *eLife* - MAJIQ LSV framework
- Trincado et al 2018 *Genome Biol* - SUPPA2
- Li et al 2018 *Nat Genet* - leafcutter
- Tapial et al 2017 *Genome Res* - VAST-TOOLS / VastDB
- Wang et al 2024 *Nat Protoc* - rMATS-turbo
- Aicher, Slaff, Jewell, Barash 2024 *bioRxiv* - MAJIQ V3
- Kubota et al 2025 *NAR* - Shiba
- Lorenzi et al 2021 *Genome Biol* - IRFinder-S
- Boutz et al 2015 *Genes Dev* - detained introns
- Irimia et al 2014 *Cell* - SRRM4 microexons
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
- Olofsson et al 2023 *Brief Bioinform* - benchmark
- Tran et al 2025 *WIREs RNA* - methodology review
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing
- Veeneman et al 2016 *Bioinformatics* - STAR 2-pass benchmark
