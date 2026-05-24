---
name: bio-chipseq-allele-specific-binding
description: Detects allele-specific transcription factor or histone modification binding from heterozygous-variant ChIP-seq using WASP (reference-bias filter; mandatory upstream), RASQUAL (joint QTL + bias-corrected testing), BaalChIP (Bayesian beta-binomial with copy-number-aware overdispersion), and AlleleSeq (personalized diploid genome). Handles imprinted-locus awareness, X-inactivation artifacts, cancer copy-number imbalance, and integration with downstream caQTL / bQTL mapping. Use when identifying variants with allelic effects on TF binding, fine-mapping causal regulatory variants, validating deep-learning variant predictions, or characterizing cis-acting regulatory effects.
tool_type: mixed
primary_tool: WASP
---

## Version Compatibility

Reference examples tested with: WASP 0.3.4+, RASQUAL 1.1+, BaalChIP 1.30+ (Bioconductor), AlleleSeq 2.0+, samtools 1.19+, bcftools 1.19+, GATK 4.5+, pysam 0.22+.

# Allele-Specific Binding (ASB)

**"Identify variants that affect transcription factor or histone modification binding in cis"** -> Compare ChIP-seq read counts at the reference and alternate alleles of heterozygous variants in a single sample. Differential read counts (ALT vs REF at hetSNPs in peaks) reveal allele-specific binding.

- CLI (mandatory bias filter): WASP `mapping pipeline` to remove reference-allele mapping bias
- CLI (joint association): RASQUAL with `--n-permutations` for cis-QTL + ASB
- R (Bayesian beta-binomial): BaalChIP with copy-number-aware overdispersion
- CLI (personalized genome): AlleleSeq with phased diploid genome
- Statistical test: beta-binomial likelihood ratio or chi-squared on count tables

ASB analysis has three universal pitfalls: reference-allele mapping bias (universal across short-read aligners), imprinted loci (constitutively allele-skewed by biology), and copy-number variation (changes effective allele dose). All three must be addressed or results are unreliable.

## Method Taxonomy

| Method | Year | Approach | Strength | Fails when |
|--------|------|----------|----------|------------|
| **WASP** (van de Geijn 2015) | 2015 | Map reads, swap alleles, re-map, drop discordant | Universal first step; aligner-agnostic; mandatory preprocessing | Drops 22-31% of reads; reduces power; not an analysis method itself |
| **RASQUAL** (Kumasaka 2016) | 2016 | Joint genotype-phenotype association with per-feature `phi` bias parameter | Improves QTL mapping; integrates bias correction; works for ChIP/ATAC/RNA-seq | Computationally intensive; assumes binomial bias structure |
| **BaalChIP** (de Santiago 2017) | 2017 | Bayesian beta-binomial; copy-number-aware overdispersion | Cancer genomes (copy-number imbalance); rigorous inference | Slower; assumes copy-number known |
| **AlleleSeq** (Rozowsky 2011) | 2011 | Personalized diploid genome alignment | Avoids reference bias completely; conceptually cleanest | Requires phased genotype + diploid genome construction; computational cost |
| **MBASED** (Mayba 2014) | 2014 | Meta-analysis-based ASE; gene-level | RNA-seq oriented; adapted for ChIP gene-body binning | Gene-level not peak-level; less precise for narrow TF peaks |
| **AllelicImbalance** (R package) | — | Bioconductor multi-method | Easy R workflow | Requires variants and BAM; less rigorous than BaalChIP |
| **deepSEA / chromBPNet variant effects** | 2015 / 2024 | Deep-learning predictions | Sequence-only; no chromatin sample needed | Predictive not measurement; see chip-deep-learning |

## Universal First Step: WASP Reference-Bias Filter

**Goal:** Remove reads that show reference-allele mapping bias before any ASB testing.

**Approach:** Align reads, identify those overlapping heterozygous SNPs, swap alleles and re-align; reads that don't map consistently to the same position with both alleles are discarded. The output is a bias-corrected BAM at the cost of 22-31% read loss.

Reference-allele mapping bias is systematic: reads with the reference allele align more readily because the reference is the alignment target. This inflates REF allele frequency by 1-5% genome-wide. WASP fixes this:

```bash
# WASP mapping pipeline
# 1. Initial alignment
bowtie2 -x hg38 -1 R1.fq -2 R2.fq -S step1.sam
samtools view -bS step1.sam | samtools sort -o step1.bam
samtools index step1.bam

# 2. Identify reads overlapping hetSNPs; swap alleles; re-map
python /path/to/WASP/mapping/find_intersecting_snps.py \
    --is_paired_end \
    --is_sorted \
    --output_dir wasp_out/ \
    --snp_tab snps_tab.h5 \
    --snp_index snps_index.h5 \
    --haplotype haplotypes.h5 \
    --samples sample_list.txt \
    step1.bam

# 3. Re-map swapped reads
bowtie2 -x hg38 -1 wasp_out/step1.remap.fq.gz -S step2.sam
# (process step2.sam similarly)

# 4. Filter reads that don't map back consistently
python /path/to/WASP/mapping/filter_remapped_reads.py \
    step1.to.remap.bam step2.bam step1.keep.bam

# 5. Final WASP-filtered BAM (use this for all downstream ASB analysis)
samtools sort -o step1.wasp.bam step1.keep.bam
samtools index step1.wasp.bam
```

**WASP always drops 22-31% of reads.** This is the cost of bias correction; downstream power is reduced but ASB calls are trustworthy.

**Alternative to WASP filter:** RASQUAL's `phi` parameter models bias within the test rather than filtering reads. More sophisticated but assumes binomial bias structure.

## Workflow: BaalChIP (Recommended for Cancer / Copy-Number-Imbalanced Samples)

```r
library(BaalChIP)
library(BSgenome.Hsapiens.UCSC.hg38)

# Sample metadata
samples <- data.frame(
    SampleID = c('HCC1395_FOXA1_rep1', 'HCC1395_FOXA1_rep2'),
    Tissue = 'TNBC',
    Target = 'FOXA1',
    BAM = c('rep1.wasp.bam', 'rep2.wasp.bam'),
    Peaks = c('rep1_peaks.bed', 'rep2_peaks.bed'),
    Group = 'HCC1395'
)

# hetSNP file: VCF or BED with chrom, pos, ref, alt, allele frequencies
hetSNPs <- 'het_snps.bed'

# CNV file for copy-number-aware overdispersion (critical for cancer)
cnvs <- 'cnvs.bed'

# Initialize BaalChIP object
res <- BaalChIP(samplesheet = samples, hets = hetSNPs)

# Run filters and Bayesian test
res <- alleleCounts(res, min_base_quality = 10, min_mapq = 15)
res <- QCfilter(res, RegionsToFilter = c('blacklist_v2.bed'))
res <- mergePerGroup(res)
res <- filter1allele(res)
res <- getASB(res, Iter = 5000, conf_level = 0.95)
# Verify parameter names against the installed BaalChIP version (`?getASB`); some releases
# use `nIter` instead of `Iter`.

# Results
asb_table <- BaalChIP.report(res)
head(asb_table)
```

BaalChIP outputs per-hetSNP: allelic ratio, posterior, Bayes factor, ASB call.

## Workflow: RASQUAL (Joint cis-QTL + ASB)

```bash
# Prepare input
# - BAM filtered by WASP
# - Genotype VCF (phased)
# - Peak BED

# Run RASQUAL
rasqual \
    --y peak_counts.txt \
    --x covariates.txt \
    --k offsets.txt \
    --n N_samples \
    --p N_peaks \
    --j 0 -i 0 \
    --vcf genotypes.vcf \
    --window 250000 \
    --t 8 \
    > rasqual_results.txt

# Output columns: chrom, peak_id, n_RSNPs, n_FSNPs, n_imputed, summarized_phi,
#                 summarized_overdispersion, summarized_pi, beta, log10_BF, ...
```

RASQUAL's `phi` parameter is the per-feature bias estimate; `pi` is the allelic ratio.

## Workflow: AlleleSeq (Personalized Diploid Genome)

```bash
# Build personalized diploid genome from phased VCF
vcf2diploid -id SAMPLE -chr hg38.fa -vcf SAMPLE.phased.vcf -outDir personalized/

# Align reads to both maternal and paternal copies
bowtie2-build personalized/maternal.fa maternal_index
bowtie2-build personalized/paternal.fa paternal_index
bowtie2 -x maternal_index -1 R1.fq -2 R2.fq -S maternal.sam
bowtie2 -x paternal_index -1 R1.fq -2 R2.fq -S paternal.sam

# AlleleSeq pipeline
AlleleSeq2.pl SAMPLE maternal.sam paternal.sam genotype.vcf

# Output: per-hetSNP allelic counts and binomial test
```

Personalized genome avoids reference bias by construction. Cost: per-sample diploid genome generation and indexing.

## Three Universal Pitfalls

### Pitfall 1: Imprinted Loci Are Constitutively Skewed

Imprinted loci (H19, IGF2, MEG3, MEG8, KCNQ1OT1, etc.) show extreme allele bias by biology, not from differential binding.

```bash
# Filter imprinted loci before ASB analysis
wget https://imprintingdiseases.org/data/imprinted_loci_hg38.bed
bedtools intersect -v -a hetSNPs.bed -b imprinted_loci_hg38.bed > hetSNPs.non_imprinted.bed
```

### Pitfall 2: X-Inactivation in Females

In female samples, X-linked genes show extreme allele skew because each cell silences one X chromosome. This appears as ASB at every X-linked hetSNP.

```bash
# Filter chrX in female samples
awk '$1 != "chrX"' hetSNPs.bed > hetSNPs.autosomal.bed
# Or analyze chrX separately with imprinting-aware methods
```

### Pitfall 3: Copy-Number Imbalance (Cancer Genomes)

In cancer cells, copy-number gain of one allele alters effective allele dose; raw allelic ratios mix dose and binding effects. BaalChIP's copy-number-aware overdispersion handles this; other methods require pre-filtering CN-altered regions.

```bash
# Use ASCAT / Sequenza / FACETS to call allele-specific CNVs
# Exclude CN-altered regions from ASB analysis OR use BaalChIP
```

## Per-Tool Failure Modes

### WASP -- Reference panel mismatch

**Trigger:** Using a WASP SNP file from a different population than the sample.

**Mechanism:** WASP swaps alleles at known hetSNPs; if the variant isn't in the SNP file, no swap happens; reads retain reference bias.

**Symptom:** Sample-specific hetSNPs (not in 1KG) still show reference bias after WASP.

**Fix:** Build WASP SNP file from the sample's own genotype VCF, not a population panel; OR use RASQUAL which handles novel hetSNPs.

### WASP -- Excessive read loss

**Trigger:** WASP filter removes >40% of reads.

**Mechanism:** Many reads span multiple hetSNPs; each must re-map consistently after every allele swap; combinatorial loss.

**Fix:** Accept the loss (genuine bias correction) OR switch to AlleleSeq (personalized genome avoids the swap-and-remap step) OR RASQUAL (no read filtering).

### RASQUAL -- Convergence failure

**Trigger:** Sparse data (few hetSNPs per peak); strong copy-number imbalance.

**Mechanism:** EM convergence requires enough hetSNPs per feature; sparse data underspecifies the model.

**Fix:** Increase `--imputation-r2-cutoff` to require well-imputed SNPs; combine replicates; or switch to BaalChIP for sparse-data robustness.

### BaalChIP -- CN file mismatch

**Trigger:** CN BED uses different naming convention (chrX vs X) than BAMs.

**Mechanism:** BaalChIP silently doesn't apply CN-aware overdispersion if CN positions don't match BAM chromosomes.

**Symptom:** ASB calls at CN-altered regions look bimodal (one allele appears 100% bound).

**Fix:** Verify chromosome naming matches across CN file, BAM, hetSNP VCF.

### AlleleSeq -- Insufficient phasing

**Trigger:** Using unphased VCF for diploid genome construction.

**Mechanism:** AlleleSeq requires phased genotypes; without phasing, maternal and paternal genomes are randomly assigned.

**Fix:** Use trio or read-based phasing (HapCUT2, WhatsHap) before AlleleSeq.

### Imprinted loci not filtered

**Trigger:** Reporting ASB at H19 or IGF2.

**Mechanism:** These loci are biologically allele-skewed; the "ASB" call is correct but uninformative.

**Fix:** Always filter imprinted loci before reporting / interpreting ASB.

### Female chrX ASB artifacts

**Trigger:** Reporting ASB at chrX in female samples without X-inactivation correction.

**Mechanism:** Random X-inactivation silences one X per cell; population of cells shows extreme allele bias at any X-linked variant.

**Fix:** Filter chrX in female samples OR use methods that model X-inactivation (rare in standard ASB pipelines).

### Reference allele bias not corrected

**Trigger:** Running BaalChIP / chi-squared test directly without WASP or RASQUAL bias handling.

**Mechanism:** 1-5% genome-wide REF allele over-representation produces false-positive REF-favoring ASB calls.

**Symptom:** ASB calls skewed toward REF allele.

**Fix:** Always apply WASP (or RASQUAL's phi parameter) before testing.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| WASP filter applied; still REF-biased | Sample-specific hetSNPs not in WASP SNP file | Use sample's own genotype VCF for WASP |
| BaalChIP and RASQUAL disagree at sparse hetSNPs | Different sparse-data behavior | BaalChIP Bayesian more conservative for sparse; check posterior |
| ASB call at imprinted locus | Biology, not differential binding | Filter imprinted loci |
| ASB at chrX in female | X-inactivation | Filter chrX |
| ASB call where copy-number altered | Cancer dose effect | Use BaalChIP with CN file OR exclude CN-altered regions |
| chromBPNet predicts strong variant effect; ASB doesn't | Sample has low coverage at variant; chromBPNet predicts in counterfactual | Increase depth; ASB requires actual chromatin sample |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| WASP `find_intersecting_snps.py` fails | h5 SNP table format wrong | Convert from VCF via WASP's `extract_vcf.py` |
| BaalChIP "no overlap with peaks" | hetSNP and peak chrom naming mismatch | Standardize chrom prefixes |
| RASQUAL OOM | Window too large; too many features | Reduce `--window`; chunk feature list |
| AlleleSeq "diploid genome too large" | Many SVs in genome | Use small-variant only VCF; exclude SV-rich regions |
| ASB calls cluster at REF allele | WASP not applied OR insufficient | Re-run WASP with sample-specific SNP file |
| Many ASB at chrX in female | X-inactivation | Filter chrX |
| All "ASB" calls are at imprinted loci | Imprinting not filtered | Apply imprinted-loci BED |

## References

- Rozowsky J et al 2011 Mol Syst Biol 7:522 (AlleleSeq)
- van de Geijn B et al 2015 Nat Methods 12:1061 (WASP)
- Kumasaka N et al 2016 Nat Genet 48:206 (RASQUAL)
- de Santiago I et al 2017 Genome Biol 18:39 (BaalChIP)
- Mayba O et al 2014 Genome Biol 15:405 (MBASED)
- Chen J et al 2016 Nat Commun 7:11101 (1000 Genomes ASB / ASE survey)

## Related Skills

- chip-seq/peak-calling - Peak calling upstream
- chip-seq/chipseq-qc - QC before ASB analysis
- chip-seq/chip-deep-learning - Validate DL variant predictions against ASB
- chip-seq/peak-annotation - Annotate ASB variants to genes / cCREs
- atac-seq/allele-specific-accessibility - Parallel ATAC ASB workflow
- causal-genomics/fine-mapping - ASB as fine-mapping orthogonal evidence
- variant-calling/variant-annotation - Annotate hetSNPs before ASB
- phasing-imputation/haplotype-phasing - Required for AlleleSeq
