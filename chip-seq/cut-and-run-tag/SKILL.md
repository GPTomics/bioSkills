---
name: bio-chipseq-cut-and-run-tag
description: Analyzes CUT&RUN (Skene Henikoff 2017) and CUT&Tag (Kaya-Okur 2019) chromatin profiling data. Handles SEACR vs MACS2 peak calling (with the btaf375 2025 benchmark guidance), pA-MNase vs pA-Tn5 vs pAG-Tn5 chimera differences, E. coli spike-in carryover normalization, IgG-only control logic (no input), characteristic fragment-size signatures (25-75 bp for CUT&Tag), and lower depth requirements (5M reads typical vs 25M for ChIP). Use when calling peaks from CUT&RUN/CUT&Tag, scaling by E. coli spike-in carryover, choosing SEACR norm mode, or comparing CUT&RUN/Tag results to traditional ChIP.
tool_type: mixed
primary_tool: SEACR
---

## Version Compatibility

Reference examples tested with: SEACR 1.4+, MACS2 2.2.9+, MACS3 3.0.4+, samtools 1.19+, bowtie2 2.5+, bedtools 2.31+, deepTools 3.5+, GoPeaks 1.0+, LanceOtron (pip).

# CUT&RUN / CUT&Tag

**"Analyze CUT&RUN or CUT&Tag chromatin profiling data"** -> Use the lower-background, lower-input alternatives to traditional ChIP. CUT&RUN tethers MNase to an antibody via Protein A; CUT&Tag tethers Tn5 via Protein A/G. Both bypass cross-linking, fragmentation, and IP washes — producing 10-100× lower background, allowing 100-1000× lower cell input, and shifting the peak-calling problem from "find signal in noise" to "find signal in near-zero background."

- Aligner: bowtie2 (CUT&RUN/Tag standard) or bwa-mem; chromap optional
- Peak calling (CUT&RUN/Tag): SEACR (Meers 2019), MACS2 with `-f BAMPE --keep-dup all`, or both for consensus
- Spike-in: E. coli carryover from bacterially-produced pA-MNase/Tn5 (automatic, variable)
- Control: IgG-only (no input control; native chromatin has no meaningful "input")

CUT&RUN/CUT&Tag has different QC thresholds, different peak calling defaults, different spike-in protocols, and different antibody requirements than traditional ChIP. Treating it as ChIP fails silently.

## Protocol Variant Taxonomy

| Variant | Chimera | Year | Use case | Failure mode |
|---------|---------|------|----------|--------------|
| **CUT&RUN** (Skene Henikoff) | pA-MNase | 2017 | Native chromatin profiling; broad antibody compatibility | Native (no fixation) — gentler; MNase digest needs careful Ca²⁺ control |
| **CUT&Tag** (Kaya-Okur Henikoff) | pA-Tn5 (rabbit only) | 2019 | Lower cell input (~5000); faster; library-ready output | Rabbit-only antibody; PCR cycles can over-amplify |
| **CUT&Tag-IT** (Active Motif) | pA-Tn5 commercial | 2020 | Standardized lots; reproducible | Cost; vendor-locked |
| **pAG-Tn5 CUT&Tag** | pAG-Tn5 | 2020 | Binds both rabbit AND mouse IgG | More versatile; identical performance otherwise |
| **AutoCut&Tag** | pAG-Tn5 plate-based | 2021 | High-throughput (96-well) | Throughput at the cost of per-sample optimization |
| **CUTAC** (CUT&Tag-then-ATAC) | pAG-Tn5 + protocol modification | 2020 | Chromatin accessibility variant of CUT&Tag | Less common; not standard CUT&Tag |
| **scCUT&Tag** | pAG-Tn5 in droplets | 2021 | Single-cell histone mark profiling | Very sparse (~1000-5000 reads/cell) |

## Algorithmic Taxonomy

| Tool | Model | Strength | Fails when |
|------|-------|----------|------------|
| **SEACR** (Meers 2019) | Empirical threshold on signal block totals; IgG-aware "stringent" mode | Designed for sparse CUT&RUN data; "stringent + norm + IgG" is the recommended default | Wrong mode (top-X% without IgG; "non" mode if no upstream spike-in normalization); broad mark with very flat signal landscape |
| **MACS2 `-f BAMPE --keep-dup all`** | Local Poisson | Familiar; integrates well with downstream tools (DiffBind) | Default `-q 0.05` may be too lenient for low-background CUT&Tag; consider `-q 0.01` |
| **GoPeaks** (Yashar 2022) | Sliding-window thresholding | Broad-mark-oriented; faster than SEACR on broad data | Newer; smaller user base |
| **LanceOtron** (Hentges 2022) | CNN trained on ENCODE peaks | Parameter-free; handles both narrow and broad | Less validated for CUT&RUN/Tag specifically; web-only or pip |
| **MACS2 + SEACR consensus** | Intersection | Highest confidence; per 2025 btaf375 benchmark, best for cross-paper reproducibility | Most conservative; may miss true peaks at marginal regions |

**2025 benchmark (Bioinformatics btaf375):**
- MACS2 better for sharp peaks (H3K4me3, TFs)
- SEACR better for broad signal (H3K27me3)
- MACS2 + SEACR consensus best for publication-grade peak sets
- "norm stringent" + IgG is the recommended SEACR default

## SEACR Workflow (Canonical CUT&RUN/Tag Caller)

**Goal:** Call CUT&RUN/CUT&Tag peaks from aligned BAMs using SEACR with IgG-aware threshold.

**Approach:** Align with Henikoff parameters, convert BAM to fragment bedGraph via bamtobed-bedpe, then invoke SEACR with `norm stringent` mode and IgG control.

```bash
# 1. Align with bowtie2 (Henikoff lab standard parameters)
bowtie2 --local --very-sensitive --no-mixed --no-discordant \
    --phred33 -I 10 -X 700 \
    -x hg38 -1 reads_R1.fq -2 reads_R2.fq \
    -S aln.sam

# 2. Convert SAM to BAM, sort, index
samtools view -bS aln.sam | samtools sort -o aln.bam
samtools index aln.bam

# 3. Generate bedGraph for SEACR (paired-end fragments)
samtools view -bS -F 0x04 aln.bam | bedtools bamtobed -bedpe -i - > aln.bedpe
awk '$1==$4 && $6-$2 < 1000 {print $0}' aln.bedpe > aln.clean.bedpe
cut -f 1,2,6 aln.clean.bedpe | sort -k1,1 -k2,2n -k3,3n > aln.fragments.bed
bedtools genomecov -bg -i aln.fragments.bed -g hg38.chrom.sizes > aln.bedgraph

# 4. Same for IgG control
# (... produce igg.bedgraph similarly ...)

# 5. SEACR with stringent + norm + IgG control (recommended default)
bash SEACR_1.4.sh aln.bedgraph igg.bedgraph norm stringent target_peaks.bed
# Output: target_peaks.bed (BED of peak regions)

# Alternative: no IgG control, use top 1% of peaks
# bash SEACR_1.4.sh aln.bedgraph 0.01 non stringent target_peaks.bed
```

**SEACR mode selection:**
- `norm` (recommended): scales target to IgG distribution
- `non`: use ONLY if upstream spike-in normalization was applied; otherwise use `norm`
- `stringent`: top-half of signal blocks (recommended default)
- `relaxed`: full distribution (use only for very sparse signal)
- IgG control: `bash SEACR_1.4.sh target.bg igg.bg norm stringent out.bed`
- Top-X% without IgG: `bash SEACR_1.4.sh target.bg 0.01 non stringent out.bed`

## E. coli Spike-In Carryover (Automatic Spike-In)

CUT&RUN/CUT&Tag spike-in is "free" because the pA-MNase or pA-Tn5 carries E. coli DNA from bacterial production. Carryover is variable across batches but stable within a batch.

```bash
# Align reads to combined hg38 + E. coli genome (or sequential)
bowtie2 -x hg38_ecoli_combined -1 R1.fq -2 R2.fq -S aln.sam

# Count E. coli reads per sample
ECOLI_READS=$(samtools view -c -F 4 aln.bam ecoli_chr1)

# Scaling factor: smallest E. coli read count / per-sample E. coli reads
# Apply BEFORE peak calling for cross-condition comparison
```

**Expected E. coli alignment fractions:**
- Target ChIP samples: 0.5-2% of total reads
- IgG control: 2-5% (higher because no target chromatin to dilute)
- Below 0.1% E. coli: spike-in carryover lost; cross-condition normalization unreliable
- Above 10%: target ChIP failed (mostly E. coli)

The carryover is variable between batches of bacterial production; single-experiment carryover spike-in is noisier than deliberate Drosophila spike-in (ChIP-Rx). For high-stakes cross-condition claims, add deliberate Drosophila spike-in despite the E. coli carryover. See chip-seq/spike-in-normalization.

## QC Differences from Traditional ChIP

| Metric | Traditional ChIP | CUT&RUN/CUT&Tag |
|--------|------------------|------------------|
| FRiP (TF) | > 0.05 | > 0.10 (often > 0.25) |
| FRiP (histone) | > 0.10 | > 0.25 |
| Library size requirement | 20-50M | 3-10M (often sufficient) |
| Input control | Required | IgG only (no input meaningful) |
| Fragment size (CUT&Tag) | Sub-nucleosomal or mono-nucleosomal | Sharp peak at 25-75 bp (Tn5 staggered cuts) |
| Fragment size (CUT&RUN) | Variable | Mono- + di-nucleosomal pattern |
| Duplicates | Remove (MarkDuplicates) | **Keep for CUT&Tag** (low PCR cycles, dups have biology) |
| Spike-in alignment | Deliberate (Drosophila); 0.5-5% | Automatic E. coli; 0.5-5% |
| Cell input | 1-10M | 5,000-100,000 |

**Critical: `--keep-dup all` in MACS for CUT&Tag.** PCR cycles are 12-15 (vs 5-8 for ChIP); duplicates at high-coverage TF binding sites contain biology. Standard ChIP convention `--keep-dup auto` will over-deduplicate CUT&Tag data.

## Fragment Size as Diagnostic (Critical for CUT&Tag)

```bash
samtools view -f 0x2 sample.bam | awk '{print $9}' | awk '$1>0' \
    | sort -n | uniq -c | awk '{print $2, $1}' > frag_sizes.tsv
```

Expected for CUT&Tag:
- Sharp peak at 25-75 bp (Tn5 staggered insertion = ~9 bp + protein-DNA-protein interaction)
- Secondary peak at ~150-200 bp (mono-nucleosomal CUT&Tag from H3K4me3 etc)
- < 25 bp: Tn5 self-tagmentation noise; high abundance indicates over-tagmentation
- Flat distribution above 200 bp: poor enzyme activity or over-amplification

Expected for CUT&RUN:
- Mono-nucleosomal (~150 bp) for histone marks
- Sub-nucleosomal (~50-100 bp) for TFs (rare; CUT&RUN better for histones)
- Di-nucleosomal (~300 bp) secondary peak common

## Per-Tool Failure Modes

### SEACR -- Wrong mode for context

**Trigger:** Using `non` mode without prior spike-in normalization; using `relaxed` mode on standard CUT&Tag.

**Mechanism:** `non` assumes target is already scaled to IgG (typical for ChIP-Rx-style spike-in); on raw counts, it inflates false positives. `relaxed` includes the full distribution; appropriate only for sparse signal.

**Fix:** Default to `norm stringent` with IgG control. Use `non` only when upstream spike-in scaling has been applied.

### CUT&Tag MACS2 -- Default `--keep-dup auto` removes biology

**Trigger:** Using MACS2 default dedup settings on CUT&Tag.

**Mechanism:** Low PCR cycles (12-15) in CUT&Tag mean PCR duplicates contain real biology at high-coverage sites; auto-dedup over-filters.

**Symptom:** Peak counts much lower than published for same antibody / cell line.

**Fix:** `macs2 callpeak --keep-dup all -f BAMPE` for CUT&Tag. For CUT&RUN, dedup behavior depends on PCR cycles — verify with library complexity (NRF).

### pA-Tn5 vs pAG-Tn5 -- Antibody species mismatch

**Trigger:** Using pA-Tn5 (Henikoff original) with mouse primary antibody.

**Mechanism:** Protein A binds rabbit IgG much better than mouse IgG; mouse antibodies give weak signal with pA-Tn5.

**Fix:** Use pAG-Tn5 (binds both); or switch to a rabbit primary antibody for the same target.

### Digitonin permeabilization -- Wrong concentration

**Trigger:** Default 0.05% digitonin on all cell lines.

**Mechanism:** Optimal digitonin varies by cell type (some need 0.02%, some 0.1%); over-permeabilization releases chromatin into supernatant; under-permeabilization prevents antibody access.

**Symptom:** Inconsistent signal across cell lines; high IgG signal (under-permeabilized) or low target signal (over-permeabilized).

**Fix:** Titrate digitonin per cell line using a known-positive H3K4me3 antibody as control.

### ConA bead vs sepharose -- Volume / sample mismatch

**Trigger:** Switching bead type between protocols without adjusting volume.

**Mechanism:** ConA magnetic beads (e.g., Bangs) and sepharose ConA have different binding capacities; protocols designed for one give wrong cell loading for the other.

**Fix:** Follow Henikoff lab protocol exactly for the chosen bead; or titrate cell number per bead volume.

### Adapter readthrough in short fragments

**Trigger:** 100-150 bp reads on 25-75 bp CUT&Tag fragments.

**Mechanism:** Reads longer than fragments read through both adapters; downstream alignment loses the fragment.

**Symptom:** Many reads with adapter sequence at 3' end; alignment rate drops.

**Fix:** Aggressive adapter trimming with cutadapt: `-e 0.1 -O 5 --minimum-length 25`. Use 50 bp paired-end sequencing for CUT&Tag instead of 150 bp.

### MACS2 fragment-size modeling failure on CUT&Tag

**Trigger:** Running MACS2 without `-f BAMPE` on CUT&Tag PE data.

**Mechanism:** MACS2 in `-f BAM` mode tries to model fragment size from cross-correlation; CUT&Tag fragments are 25-75 bp, not the 200 bp ChIP expects; modeling fails or produces wrong estimate.

**Fix:** Always use `-f BAMPE` for CUT&Tag; MACS uses actual fragment spans from mate pairs.

## Reconciliation: When CUT&RUN/Tag Disagrees with ChIP

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Peak count much lower in CUT&Tag vs ChIP | Lower background reveals signal vs noise; CUT&Tag often has FEWER but cleaner peaks | Both correct; CUT&Tag specificity > sensitivity |
| Peak count much higher in CUT&Tag vs ChIP | `--keep-dup all` retained PCR duplicates as peaks | Verify NRF; if low, consider deduplicating with caution |
| Same antibody, different signal | Native chromatin vs cross-linked accessibility differs | Native CUT&RUN may miss DSG-dependent cofactors (BRD4); add brief fixation |
| FRiP very high (>50%) | Likely real for CUT&Tag (low background); confirm with motif enrichment | Verify motif enrichment at peaks; if missing, suspect technical artifact |
| H3K4me3 CUT&Tag peak count differs from ChIP | Expected; CUT&Tag has higher specificity | Trust CUT&Tag for sharp marks |
| H3K27me3 CUT&RUN/Tag misses regions | Broad domains require deeper sequencing; CUT&Tag was designed for sharp marks | Use CUT&RUN or traditional ChIP for very broad marks |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| SEACR "input file not bedgraph" | Wrong file format | Use `bedtools genomecov -bg`; not `bedGraphToBigWig` output |
| MACS2 modeling fails on CUT&Tag | Default `-f BAM` on PE | `-f BAMPE` |
| Peak count for mouse antibody | Used pA-Tn5 not pAG-Tn5 | Switch to pAG-Tn5 OR use rabbit primary |
| Very high adapter content in FASTQ | Read length > fragment length | Trim aggressively; consider 50 bp PE for CUT&Tag |
| Sample-to-sample carryover variability >5x | E. coli carryover variable between bacterial production batches | Add deliberate Drosophila spike-in for cross-condition |
| IgG signal as strong as target | Failed antibody / over-permeabilization | Validate antibody on positive control; titrate digitonin |

## References

- Skene PJ & Henikoff S 2017 eLife 6:e21856 (CUT&RUN)
- Skene PJ Henikoff JG & Henikoff S 2018 Nat Protoc 13:1006 (CUT&RUN protocol)
- Kaya-Okur HS et al 2019 Nat Commun 10:1930 (CUT&Tag)
- Kaya-Okur HS et al 2020 Nat Protoc 15:3264 (CUT&Tag protocol)
- Meers MP et al 2019 Epigenetics Chromatin 12:42 (SEACR)
- 2025 Bioinformatics btaf375 (CUT&RUN peak caller benchmark: MACS2 + SEACR consensus)
- Yashar A et al 2022 Bioinformatics Adv 2:vbac085 (GoPeaks)
- Hentges LD et al 2022 Bioinformatics 38:4255 (LanceOtron)
- Bartosovic M Kabbe M & Castelo-Branco G 2021 Nat Biotechnol 39:825 (scCUT&Tag)
- Henikoff S et al 2021 STAR Protocols 2:100822 (AutoCut&Tag)
- 2024 Nat Methods (scNanoSeq-CUT&Tag long-read scCUT&Tag)

## Related Skills

- chip-seq/peak-calling - Traditional ChIP peak calling (MACS3 + IDR vs naive overlap)
- chip-seq/chipseq-qc - QC battery (different thresholds for CUT&RUN/Tag)
- chip-seq/spike-in-normalization - Deliberate Drosophila spike-in beyond E. coli carryover
- chip-seq/differential-binding - DiffBind / csaw differential on CUT&RUN/Tag
- chip-seq/peak-annotation - Annotate CUT&RUN/Tag peaks (same tools as ChIP)
- chip-seq/super-enhancers - SE calling on CUT&Tag H3K27ac
- alignment-files/sam-bam-basics - BAM preparation
- read-qc/adapter-trimming - Aggressive trimming for short fragments
