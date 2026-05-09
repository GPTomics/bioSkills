---
name: bio-long-read-splicing
description: Analyzes alternative splicing from PacBio Iso-Seq (HiFi, Kinnex/MAS-Iso-seq) and Oxford Nanopore (direct cDNA, direct RNA, R10.4.1+) long-read RNA-seq with full-isoform resolution. Tools include FLAIR (correct/collapse/quantify/diffSplice for PacBio + ONT), IsoQuant (de-novo or annotation-guided isoform discovery 2024 SOTA), Bambu (annotation-aware Bayesian discovery + quantification with Novel Discovery Rate), SQANTI3/SQANTI-LR (isoform classification: FSM/ISM/NIC/NNC + artifact flags), rMATS-long (event calling on long-read isoforms), and minimap2 (-ax splice:hq for HiFi; -ax splice -k14 for ONT cDNA; add -uf only for direct RNA or stranded cDNA preps). Solves microexon detection, recursive splicing, complex multi-exon isoforms, and DTU without transcript-quantification uncertainty. Use when short-read AS limitations (anchor length, complex isoforms, microexons, recursive splicing, transcript ambiguity) demand full-isoform resolution.
tool_type: mixed
primary_tool: FLAIR
---

## Version Compatibility

Reference examples tested with: FLAIR 2.0+, IsoQuant 3.5+, Bambu 3.4+, SQANTI3 5.2+, minimap2 2.26+, samtools 1.19+, rMATS-long 0.2+, IsoSeq3 4.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Long-Read Splicing Analysis

Full-length long-read sequencing solves problems that short-read AS cannot: anchor-length-limited microexon detection, complex multi-exon isoform deconvolution, recursive splicing in long introns, and transcript-quantification uncertainty in DTU. The 2024-2026 transition: long-read is becoming the splicing default for high-resolution analysis.

## When Long-Read Wins

| Question | Why long-read wins |
|----------|---------------------|
| Microexon detection (3-27 nt) | Reads span the microexon entirely; no aligner anchor problem |
| Long-intron recursive splicing | Can detect ratchet point usage (Sibley 2015 *Nature*) |
| Complex isoform deconvolution (TTN, MAPT, NEFM) | Single read per isoform avoids EM ambiguity |
| DTU without quantification uncertainty | Transcript identity is read-level, not inferred |
| Novel transcript discovery | No annotation dependence |
| Phasing splicing with SNVs | Allele-resolved isoforms |
| Single-cell full-length isoforms | MAS-Iso-seq + 10X 5' is the practical SOTA |
| Cryptic splicing in TDP-43 ALS | Full-length reads confirm cryptic exon inclusion in target transcripts |

## Platform Selection Matrix

| Platform | Throughput | Accuracy (modal) | Best for | Fails when |
|----------|------------|------------------|----------|------------|
| PacBio Revio HiFi (Iso-Seq) | ~25M reads / SMRT cell | Q30+ (CCS) | Bulk transcript discovery; gold standard | Cost prohibitive for very large cohorts |
| PacBio Kinnex / MAS-Iso-seq | ~16x Iso-Seq via concatemer | Q30+ | High-throughput single-cell long-read | Kinnex de-array (skera) is an extra step |
| ONT direct cDNA (R10.4.1, PCS-114) | Millions / flowcell | ~98% simplex, ~99% duplex | Cost-effective; throughput | Minor higher error than HiFi |
| ONT direct RNA (RNA004, 2024+) | ~30M reads | ~96-98% | Native modifications (m6A, pseudo-U); no RT bias | Lower throughput; higher input |
| ONT pre-R10 (R9.4.1) | Same as R10 | ~85-90% | Legacy data | Pre-R10 not recommended for splicing analysis (false novel junctions) |

**Read length:** PacBio HiFi cdna typically 1-10 kb; ONT cdna 0.5-50+ kb (long-tailed). Both span typical mammalian transcripts. Direct RNA on ONT preserves true 5'/3' termini and modifications.

## Decision Tree by Use Case

| Use case | Recommended tools |
|----------|--------------------|
| Bulk Iso-Seq transcript discovery in well-annotated organism | minimap2 -ax splice:hq → IsoQuant or Bambu → SQANTI3 |
| Bulk ONT cDNA in well-annotated organism | minimap2 -ax splice -uf -k14 → IsoQuant or FLAIR → SQANTI3 |
| End-to-end pipeline for differential analysis | FLAIR (correct → collapse → quantify → diffSplice) |
| Joint discovery + quantification with calibrated novel rate | Bambu in R |
| De novo discovery for non-model organism | IsoQuant with --genedb omitted |
| Event-level differential splicing on long reads | rMATS-long |
| DTU on long-read transcript counts | DRIMSeq → DEXSeq/satuRn → stageR (no Salmon Gibbs needed) |
| Hybrid short+long for cohort | StringTie2 hybrid + FLAIR / IsoQuant |
| Single-cell full-length isoforms | MAS-Iso-seq + 10X 5' → FLAMES or scNanoGPS |
| Cryptic exon validation in ALS | minimap2 → FLAIR collapse → manual inspection of UNC13A, STMN2 |
| ASO design with full-isoform context | minimap2 → IsoQuant → SQANTI3 → ASO design (see splice-variant-prediction) |

## Splice-Aware Alignment

```bash
# PacBio HiFi (Iso-Seq) -> minimap2 splice:hq preset
minimap2 -ax splice:hq -uf --secondary=no \
    -t 16 \
    reference.fa \
    isoseq.fastq.gz | \
    samtools sort -@ 8 -o isoseq_aligned.bam
samtools index isoseq_aligned.bam

# ONT direct cDNA (PCS-114, PCB-114): unstranded by default; omit -uf
minimap2 -ax splice -k14 \
    -t 16 \
    reference.fa \
    ont_cdna.fastq.gz | \
    samtools sort -@ 8 -o ont_cdna_aligned.bam
samtools index ont_cdna_aligned.bam

# ONT direct RNA (RNA004): truly stranded (RNA molecule preserves direction); -uf is correct
minimap2 -ax splice -uf -k14 \
    -t 16 \
    reference.fa \
    ont_rna.fastq.gz | \
    samtools sort -@ 8 -o ont_rna_aligned.bam
samtools index ont_rna_aligned.bam
```

`-uf` forces all reads to the forward transcript strand — correct for direct RNA (single-stranded) and stranded cDNA library preps; **omit for unstranded cDNA** (default ONT PCS/PCB kits) or you lose ~half the reads. `--secondary=no` discards secondary alignments. For genomes with poorly-annotated splice sites, supplement with `--junc-bed gencode_junctions.bed`. uLTRA (Sahlin & Mäkinen 2021 *Bioinformatics*) and deSALT (Liu 2019 *Genome Biol*) are alternatives with higher precision on small/cryptic exons.

**Critical:** `splice:hq` is the preset for HiFi (Q30+ reads); plain `splice` is for ONT regardless of cDNA vs direct RNA. Using `splice` on HiFi data underuses the high quality; using `splice:hq` on ONT misses true junctions due to error-tolerance mismatch.

## FLAIR Workflow (correct → collapse → quantify → diffSplice)

**Goal:** Identify, quantify, and test full-length isoforms from long-read RNA-seq across conditions.

**Approach:** Correct splice junctions against short-read or annotation evidence, collapse isoforms, quantify per-sample expression, run diffSplice for differential isoform usage.

```bash
flair correct \
    --query aligned.bed \
    --genome reference.fa \
    --gtf gencode.v45.annotation.gtf \
    --shortread short_read_junctions.bed \
    --output flair_corrected \
    --threads 16

flair collapse \
    --query flair_corrected_all_corrected.bed \
    --reads sample.fastq.gz \
    --genome reference.fa \
    --gtf gencode.v45.annotation.gtf \
    --output flair_collapsed \
    --threads 16

flair quantify \
    --reads_manifest reads_manifest.tsv \
    --isoforms flair_collapsed.isoforms.fa \
    --output flair_quantified \
    --threads 16

flair diffSplice \
    --isoforms flair_collapsed.isoforms.bed \
    --counts_matrix flair_quantified_counts.tsv \
    --conditions_table conditions.tsv \
    --output flair_diffsplice \
    --threads 16
```

FLAIR (Tang 2020 *Nat Commun*) handles ONT and PacBio with the same workflow. Output includes per-event PSI, FDR, and visual sashimi-like plots. The `--shortread` flag for `flair correct` is **strongly recommended** when short-read RNA-seq is available — it dramatically improves splice junction precision.

## IsoQuant for Discovery + Quantification

**Goal:** De novo or annotation-guided isoform discovery and quantification with high precision.

**Approach:** Run `isoquant.py` with reference + reads + data type; output is GTF + counts.

```bash
isoquant.py \
    --reference reference.fa \
    --genedb gencode.v45.annotation.gtf \
    --fastq sample1.fastq.gz sample2.fastq.gz \
    --data_type pacbio_ccs \
    --output isoquant_output \
    --threads 16 \
    --model_construction_strategy default_pacbio
```

`--data_type` accepts `pacbio_ccs` (HiFi), `nanopore` (ONT), or `assembly`. As of v3.0+, `--genedb` is optional for de novo discovery. IsoQuant (Prjibelski 2023 *Nat Biotech*) is current SOTA for novel transcript reconstruction; pairs well with SQANTI3 for downstream classification.

Memory requirement: >=64 GB for atlas-scale runs.

## Bambu for Annotation-Aware Discovery + Quantification

**Goal:** Joint discovery and quantification with statistical filtering of novel isoforms.

**Approach:** R Bioconductor package; takes BAM + reference annotation + genome; outputs ranged SE objects of known + novel transcripts.

```r
library(bambu)

bam_files <- c('sample1.bam', 'sample2.bam', 'sample3.bam')
genome <- 'reference.fa'
gtf <- 'gencode.v45.annotation.gtf'

bambuAnnotations <- prepareAnnotations(gtf)

se <- bambu(
    reads = bam_files,
    annotations = bambuAnnotations,
    genome = genome,
    NDR = 0.1,
    ncore = 8
)

writeBambuOutput(se, path = 'bambu_output/')

tx_counts <- as.data.frame(assays(se)$counts)
gene_counts <- transcriptToGeneExpression(se)
```

Bambu (Chen 2023 *Nat Commun*) uses **NDR** (Novel Discovery Rate) as a single, calibrated parameter replacing per-sample heuristics:

| NDR | Interpretation |
|-----|----------------|
| 0.05 | Stringent; few novel transcripts; highest precision |
| 0.1 | Balanced (default) |
| 0.2-0.3 | Permissive; more novel discoveries; recall over precision |

Excellent for combined discovery + quantification when statistical filtering matters.

## SQANTI3 Classification

**Goal:** Classify discovered isoforms relative to reference; flag artifacts (intra-priming, RT-switching).

**Approach:** Run `sqanti3_qc.py` on the isoform GTF; review classification (FSM/ISM/NIC/NNC/antisense/genic/intergenic/fusion) and quality flags.

```bash
sqanti3_qc.py \
    isoforms.gtf \
    gencode.v45.annotation.gtf \
    reference.fa \
    --output sqanti3_qc \
    --aligner_choice minimap2 \
    --cage_peak refTSS_v3.3_human_coordinate.hg38.bed \
    --polyA_motif_list mouse_and_human.polyA_motif.txt \
    --cpus 8

sqanti3_filter.py rules \
    sqanti3_qc_classification.txt \
    --isoforms isoforms.fa \
    --gtf isoforms.gtf \
    --output sqanti3_filtered
```

| SQANTI category | Meaning |
|-----------------|---------|
| FSM (Full Splice Match) | All junctions match reference |
| ISM (Incomplete Splice Match) | Subset of reference junctions |
| NIC (Novel In Catalog) | Novel combination of known junctions |
| NNC (Novel Not in Catalog) | Contains novel junction |
| Antisense | Overlaps gene on opposite strand |
| Genic | Within gene but no junction match |
| Intergenic | Between genes |
| Fusion | Spans multiple genes |

SQANTI-LR (Pardo-Palacios 2024 *Nat Methods*) is the long-read-specific branch with QC tailored to ONT/PacBio error patterns. **Filter intra-priming and RT-switching** flags before reporting.

## rMATS-long for Differential Isoform Analysis on Long-Read Data

**Goal:** Apply differential isoform analysis to long-read transcript abundance with classification and visualization.

**Approach:** rMATS-long is a multi-script Python pipeline distributed via bioconda; entry point is `rmats-long` followed by the script name. It supports two modes: **abundance-based** (using ESPRESSO-style abundance estimates) and **ASM-based** (Alternative Splicing Modules — sets of isoforms sharing exon-junction structure). Run preprocessing scripts in order before `rmats_long.py`.

```bash
conda install -c conda-forge -c bioconda rmats-long

# Preprocessing pipeline (ASM mode); per-script flag names verified vs Xinglab/rmats-long
rmats-long organize_gene_info_by_chr.py --gtf annotation.gtf --out-dir gene_info_by_chr/

# simplify_alignment_info processes one BAM at a time -> one TSV
for bam in *.bam; do
    rmats-long simplify_alignment_info.py --in-file "$bam" --out-tsv "alignment_info/${bam%.bam}.tsv"
done

# organize_alignment_info_by_gene_and_chr requires a samples-tsv (sample_id<TAB>tsv_path)
rmats-long organize_alignment_info_by_gene_and_chr.py \
    --gtf-dir gene_info_by_chr/ \
    --out-dir organized/ \
    --samples-tsv samples.tsv

rmats-long detect_splicing_events.py --align-dir organized/ --out-dir events/
rmats-long create_gtf_from_asm_definitions.py --event-dir events/ --out-gtf asm.gtf
rmats-long count_reads_for_asms.py --align-dir organized/ --event-dir events/ --out-dir asm_counts/

# Main differential analysis (ASM mode)
rmats-long rmats_long.py \
    --group-1 ctrl1.bam,ctrl2.bam,ctrl3.bam \
    --group-2 trt1.bam,trt2.bam,trt3.bam \
    --event-dir events/ \
    --asm-counts-dir asm_counts/ \
    --align-dir organized/ \
    --gtf-dir gene_info_by_chr/ \
    --out-dir rmats_long_output/ \
    --adj-pvalue 0.05 \
    --delta-proportion 0.05 \
    --average-reads-per-group 10

# Alternative: abundance-based mode (when you already have ESPRESSO-style estimates)
rmats-long rmats_long.py \
    --abundance abundance.esp \
    --updated-gtf updated.gtf \
    --group-1 sample1,sample2,sample3 \
    --group-2 sample4,sample5,sample6 \
    --out-dir rmats_long_output/ \
    --no-splice-graph-plot
```

Key flags: `--adj-pvalue` (default 0.05), `--delta-proportion` (default 0.05), `--average-reads-per-group` (default 10), `--no-splice-graph-plot` (skip expensive splice-graph rendering).

rMATS-long is a separate tool from short-read rMATS-turbo. The predecessor `lr2rmats` used long reads only to *augment* the short-read rMATS GTF. The ASM framework treats AS as a set-of-isoforms problem, more natural for long-read data than rMATS-turbo's pre-defined event categories.

## DTU on Long-Read Counts

**Goal:** Apply DRIMSeq + DEXSeq + stageR DTU pipeline to long-read transcript counts (no quantification uncertainty).

**Approach:** Use FLAIR or Bambu transcript counts as input; long-read counts are read-level identities, so no Salmon Gibbs samples needed.

```r
library(DRIMSeq); library(DEXSeq); library(stageR)

counts <- read.table('flair_quantified_counts.tsv', header=TRUE, sep='\t')

samples <- data.frame(
    sample_id = c('s1', 's2', 's3', 's4', 's5', 's6'),
    condition = c('ctrl', 'ctrl', 'ctrl', 'trt', 'trt', 'trt')
)

d <- dmDSdata(counts = counts, samples = samples)
d <- dmFilter(
    d,
    min_samps_feature_expr = 3, min_feature_expr = 5,
    min_samps_feature_prop = 3, min_feature_prop = 0.1,
    min_samps_gene_expr = 6, min_gene_expr = 10
)
```

Then proceed with the standard DEXSeq + stageR DTU pipeline (see `isoform-switching` skill). IsoformSwitchAnalyzeR v2 has explicit long-read input support.

## Single-Cell Long-Read for Splicing

**Goal:** Combine cell typing (10X 5' short read) with full-length isoform structure (Kinnex / MAS-Iso-seq).

**Approach:** Split 10X library; sequence half short-read for cell typing, half PacBio Kinnex for isoforms; recover cell barcodes from long reads via FLAMES or skera (Kinnex de-array).

```bash
# Demultiplex MAS-Iso-seq reads
skera split \
    raw_kinnex.bam \
    mas12_primers.fasta \
    demuxed.bam

# Then proceed with lima -> isoseq3 refine -> isoseq3 cluster pipeline
# For barcode rescue from FLAMES:
match_cell_barcode \
    --bam demuxed.bam \
    --barcodes 10x_barcodes.tsv \
    --output flames_demuxed.bam
```

Joglekar 2024 *Nat Neurosci* used this approach for the mouse cortex isoform atlas. See `single-cell-splicing` for tools that work on the demultiplexed data.

## Per-Tool Failure Modes

### minimap2: Wrong Preset

**Trigger:** Using `-ax splice` for PacBio HiFi (instead of `-ax splice:hq`) or `-ax splice:hq` for ONT.

**Mechanism:** Presets configure k-mer size, error tolerance, and indel scoring; mismatched preset is sub-optimal.

**Symptom:** Lower alignment rate; missed junctions on HiFi, false novel junctions on ONT.

**Fix:** `splice:hq` for HiFi; `splice -k14` for ONT cDNA (unstranded); add `-uf` only for ONT direct RNA or stranded cDNA preps.

### IsoQuant: Memory Pressure

**Trigger:** Atlas-scale cohort or low-RAM environment.

**Mechanism:** IsoQuant builds graph structures across all reads simultaneously.

**Symptom:** OOM kill; very slow runtime.

**Fix:** Increase RAM to >=64 GB; or batch by chromosome.

### Bambu: NDR Mistuning

**Trigger:** NDR=0.5+ or NDR=0.01.

**Mechanism:** NDR controls the precision-recall tradeoff for novel transcripts.

**Symptom:** Too many spurious novel transcripts (high NDR) or missing real novel transcripts (low NDR).

**Fix:** Default NDR=0.1 is balanced; adjust based on validation expectations.

### SQANTI3: RT-Switching Flags

**Trigger:** PacBio/ONT cDNA libraries with template switching artifacts.

**Mechanism:** RT-switching produces chimeric reads spanning two unrelated transcripts; SQANTI3 flags these.

**Symptom:** Many "fusion" transcripts in non-cancer samples; biologically implausible.

**Fix:** Filter out RT-switching flags via `sqanti3_filter.py`; investigate library prep if rate >5%.

### FLAIR: Short-Read Augmentation Missing

**Trigger:** Running `flair correct` without `--shortread`.

**Mechanism:** FLAIR uses short-read junctions to correct long-read junction calls; without them, long-read errors persist as junction calls.

**Symptom:** Many false novel junctions; junction precision low.

**Fix:** Always include `--shortread short_read_junctions.bed` when short-read RNA-seq is available; generate with regtools junctions.

### rMATS-long: GTF-Only Input

**Trigger:** Trying to give rMATS-long raw long-read BAMs.

**Mechanism:** rMATS-long expects per-sample isoform GTFs (from FLAIR/IsoQuant collapse), not raw alignments.

**Symptom:** Confusing parsing errors.

**Fix:** Run FLAIR/IsoQuant per sample first; pass the resulting GTFs.

## Reconciliation: When Long-Read Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FLAIR has more isoforms than IsoQuant | FLAIR collapse less stringent; or IsoQuant filtered more aggressively | Both tools have valid pipelines; report based on use case |
| Bambu calls fewer novel than IsoQuant | Bambu NDR=0.1 is more conservative | Adjust NDR or trust Bambu's calibration |
| SQANTI3 classifies as NNC, FLAIR thinks FSM | GENCODE version mismatch | Verify both tools use same annotation |
| Long-read isoform calls don't match short-read events | Short-read EM ambiguity; or long-read coverage gap | Trust long-read for unambiguous; trust short-read for high-coverage events |

## Quality Control for Long-Read Splicing

| Metric | PacBio HiFi | ONT cDNA R10.4.1 |
|--------|-------------|-------------------|
| Read accuracy (modal) | Q30+ (>=99.9%) | ~98% simplex / ~99% duplex |
| Splice junction concordance to short-read truth | ~98% | 95-98% |
| Median read length (transcripts) | 1-4 kb | 0.5-3 kb |
| Throughput per run | ~25M HiFi reads | Tens of millions |
| Library input | 100-500 ng total RNA | 100-500 ng |
| Read direction | TSO + dT primed | TSO or random hexamer |

Pre-R10 ONT (R9.4.1) had ~85-90% junction concordance and is no longer recommended for splicing.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `minimap2: too many anchors` | Repeat-rich genome region | Use `-N 50` to limit secondary alignments |
| `IsoQuant: ssw-py not found` | Missing dependency | `pip install ssw-py` |
| `Bambu: prepareAnnotations failed` | GTF malformed | Validate GTF with `gffread -E` |
| `SQANTI3: kallisto not found` | sqanti3 expects kallisto for short-read overlap | `conda install -c bioconda kallisto` |
| `FLAIR: flair correct slow` | Genome FASTA not indexed | `samtools faidx reference.fa` |
| `skera: too many mismatches in adapter` | MAS primer mismatch | Verify primer fasta matches kit version |

## Quality Thresholds

| Metric | Recommendation | Source |
|--------|----------------|--------|
| Full-length non-chimeric (FLNC) % | >=80% (PacBio Iso-Seq) | PacBio convention |
| FSM% | >=50% in well-annotated genome | Tardaguila 2018 *Genome Res* |
| NNC% | <=30% (>30% suggests artifacts unless biologically interesting) | SQANTI3 convention |
| Junction support | >=2 reads (or >=3 with strict filtering) | Conservative |
| Bambu NDR | 0.1 default; 0.05 stringent | Chen 2023 *Nat Commun* |
| SQANTI3 RT-switching flag | filter out unless validated | SQANTI3 convention |
| SQANTI3 intra-priming flag | filter out | SQANTI3 convention |
| ONT R-version | R10.4.1+ for splicing | Splice junction concordance >=95% only with R10+ |
| HiFi CCS passes | >=3 | PacBio convention for Q30+ |

## Common Pitfalls

- **Ignoring reference annotation completeness** — SQANTI3 NNC categorization differs by GENCODE version; report version with results.
- **Not running isoseq3 refine** — concatemers and polyA artifacts inflate isoform counts.
- **Confusing FLAIR's 'collapse' with 'cluster'** — collapse merges similar isoforms post-alignment; cluster (in isoseq3) merges raw reads pre-alignment.
- **Treating ONT R9.x splice calls as reliable** — pre-R10.4.1 error patterns generate false novel junctions.
- **Skipping CAGE / polyA validation in SQANTI3** — TSS / TTS hallucination is common in long-read isoforms.
- **DTU on too few replicates** — long-read is expensive; n=2 vs n=2 is common but underpowered.
- **PacBio HiFi alignment with `-ax splice` (not `splice:hq`)** — use the HQ preset for HiFi data; default `splice` is for ONT.
- **Skipping `--shortread` in FLAIR correct** — long-read junction precision is much higher with short-read augmentation.

## Related Skills

- splicing-quantification - Short-read PSI for cross-validation
- isoform-switching - DTU framework on long-read counts
- single-cell-splicing - MAS-Iso-seq + 10X integration
- long-read-sequencing/isoseq-analysis - PacBio Iso-Seq general pipeline (CCS, lima, refine, cluster)
- long-read-sequencing/long-read-alignment - minimap2 splice:hq details
- long-read-sequencing/long-read-qc - QC for long-read data
- splice-variant-prediction - Cross-reference variant predictions with full isoforms

## References

- Tang et al 2020 *Nat Commun* - FLAIR
- Prjibelski et al 2023 *Nat Biotech* - IsoQuant
- Chen et al 2023 *Nat Commun* - Bambu
- Tardaguila et al 2018 *Genome Res* - SQANTI (original)
- Pardo-Palacios et al 2024 *Nat Methods* - SQANTI3 / LRGASP benchmark
- Wyman et al 2020 *bioRxiv* - TALON (note: not formally peer-reviewed)
- Li 2018 / 2021 *Bioinformatics* - minimap2
- Sahlin & Makinen 2021 *Bioinformatics* - uLTRA
- Sibley et al 2015 *Nature* - recursive splicing
- Al'Khafaji et al 2024 *Nat Biotech* - MAS-Iso-seq / Kinnex
- Joglekar et al 2024 *Nat Neurosci* - scISOr-Seq2 mouse cortex
- Tian et al 2021 *Nat Methods* - FLAMES
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing
