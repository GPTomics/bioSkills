---
name: bio-atac-seq-allele-specific-accessibility
description: Detect allele-specific chromatin accessibility from ATAC-seq using WASP, GATK ASEReadCounter, or RASQUAL. Use when mapping cis-regulatory genetic variants from heterozygous SNPs, separating cis from trans regulation, building chromatin QTL (caQTL) maps, validating GWAS variant function with allelic imbalance, or detecting reference allele mapping bias before downstream analysis.
tool_type: mixed
primary_tool: WASP
---

## Version Compatibility

Reference examples tested with: WASP 0.3.4+, GATK 4.4+, RASQUAL 1.1+, samtools 1.19+, bcftools 1.19+, vcftools 0.1.16+, plink 2.00+, MatrixEQTL 2.3+, QuASAR 0.1+, bowtie2 2.5+, bwa-mem2 2.2.1+, scipy 1.11+ (false_discovery_control), pandas 2+, pybedtools 0.10+.

Verify before use:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Allele-Specific Accessibility

**"Does this heterozygous SNP affect chromatin accessibility on its allele?"** -> Count ATAC reads supporting reference vs alternative allele at heterozygous sites in the same individual; significant deviation from 50:50 indicates cis-regulatory effect. Requires careful handling of reference-allele mapping bias (WASP filtering) and within-individual binomial testing.

- CLI: `WASP` (Geijn 2015) for de-biased reference mapping
- CLI: `gatk ASEReadCounter` for allele-specific count tables
- CLI: `RASQUAL` (Kumasaka 2016) for joint cis-mapping with allelic counts
- R: `QuASAR` (Harvey 2015) for genotype + ASE inference simultaneously

ASE/ASB analysis is fundamentally different from cohort-level differential. Statistical framework is binomial within-individual; sample size is the count of heterozygous SNPs in accessible regions, not the number of individuals.

## Algorithmic Taxonomy

| Tool | Method | Input | Strength | Fails when |
|------|--------|-------|----------|------------|
| WASP (Geijn 2015) | Realign reads where alt allele swap could change mapping; filter mapping-bias affected sites | BAM + VCF + reference | Mandatory for any ASE/ASB analysis; controls reference-allele bias | Slow on deep coverage; requires re-alignment step |
| GATK ASEReadCounter | Count REF and ALT reads at known heterozygous sites | BAM + VCF | Mature; integrates with GATK ecosystem; standard counter | Doesn't fix mapping bias (needs WASP first); single-sample |
| RASQUAL (Kumasaka 2016) | Joint cis-eQTL/caQTL model: total counts + allelic imbalance | BAM + VCF + peak counts | Best statistical power for caQTL when sample size is moderate (N=20-100); models genotype uncertainty | Complex setup; per-feature regression (slow); requires LD computation |
| QuASAR (Harvey 2015) | Genotype + ASE inference from RNA-seq or ATAC alone | BAM (no VCF needed) | Useful when genotypes are limited; integrates phasing | Less accurate than WASP+GATK when genotypes are known |
| MatrixEQTL on per-feature counts | Linear model fit on accessibility per peak | Genotypes + peak counts | Standard cohort-level caQTL; well-supported | No allelic imbalance information; needs sample size N >= 50 |
| Allelic Imbalance from Bayesian models (MAJIQ-style) | Bayesian beta-binomial | BAM + VCF | Models overdispersion appropriately | Niche; less mature than WASP+GATK |

Methodology evolves; verify against current Geijn 2015, Kumasaka 2016, Buchkovich 2015 before locking pipelines. Modern caQTL studies typically combine WASP + RASQUAL at modest N or WASP + GATK + linear caQTL at large N.

## Reference Allele Mapping Bias (The Single Most Important Issue)

When aligning reads to the reference genome, reads carrying the reference allele align with 0 mismatches; reads carrying the alternative allele have 1 mismatch and may fail to align (especially with `bwa-mem` -k or `bowtie2 --very-sensitive` thresholds). This **inflates the apparent reference-allele frequency** at every SNP and confounds ASE.

**Trigger:** Always (mapping bias is universal at heterozygous SNPs).

**Mechanism:** Aligners are mismatch-penalized. Without correction, ALT reads systematically under-align. Effect size: 1-5% bias toward reference at typical mismatch penalties.

**Symptom:** Per-SNP REF allele fraction skews above 50% even at sites without true allelic imbalance.

**Fix:** WASP. Pseudo-alleles every read at heterozygous sites; re-aligns swapped reads; keeps only reads that map identically to both haplotypes. Mandatory for any ASE/ASB analysis.

**Goal:** Remove reference-allele mapping bias before counting allele-specific ATAC reads.

**Approach:** Identify reads overlapping heterozygous SNPs, re-align allele-swapped versions, keep only reads consistent across both haplotypes, then count REF/ALT with GATK ASEReadCounter.

```bash
# WASP read-correction pipeline (Geijn 2015)
WASP_DIR=/path/to/WASP
PEAKS=peaks.bed                                  # ATAC peaks for filtering
SAMPLE=sample1
OUT=$SAMPLE.wasp_filtered.bam

# 1. Find reads at SNP sites
python $WASP_DIR/mapping/find_intersecting_snps.py \
    --is_paired_end \
    --is_sorted \
    --output_dir wasp_out/ \
    --snp_tab snp_h5/snp_tab.h5 \
    --snp_index snp_h5/snp_index.h5 \
    --haplotype snp_h5/haps.h5 \
    --samples $SAMPLE \
    $SAMPLE.bam

# 2. Re-align flipped-allele reads
bowtie2 -x hg38_idx -1 wasp_out/$SAMPLE.remap.fq1.gz \
                   -2 wasp_out/$SAMPLE.remap.fq2.gz \
                   -S wasp_out/$SAMPLE.remap.sam
samtools view -bS wasp_out/$SAMPLE.remap.sam | samtools sort -o wasp_out/$SAMPLE.remap.bam
samtools index wasp_out/$SAMPLE.remap.bam

# 3. Keep only consistently-mapped reads
python $WASP_DIR/mapping/filter_remapped_reads.py \
    wasp_out/$SAMPLE.to.remap.bam \
    wasp_out/$SAMPLE.remap.bam \
    wasp_out/$SAMPLE.kept.bam

# 4. Merge kept reads with non-overlapping reads
samtools merge $OUT \
    wasp_out/$SAMPLE.kept.bam \
    wasp_out/$SAMPLE.keep.bam

# After WASP filtering, GATK ASEReadCounter is safe
gatk ASEReadCounter \
    -I $OUT \
    -V heterozygous_snps.vcf \
    -R hg38.fa \
    -O $SAMPLE.ase_counts.tsv
```

## Per-Tool Failure Modes

### GATK ASEReadCounter without WASP -- Reference bias

**Trigger:** Running ASEReadCounter directly on a standard ATAC BAM without WASP filtering.

**Mechanism:** Reference allele over-counts due to alignment bias.

**Symptom:** Per-SNP reference fraction systematically > 0.5; aggregate plots show ~0.51-0.55 instead of 0.5.

**Fix:** WASP filter first, ALWAYS. There are no exceptions.

### Sample size for ASE per SNP

**Trigger:** Single-individual ATAC; per-SNP heterozygous coverage typically 10-100 reads.

**Mechanism:** Per-SNP binomial test has limited power; 10 reads at p=0.5 has 95% CI from 0.18 to 0.82 -- effectively no power for moderate effects.

**Fix:** Aggregate across many SNPs in the same peak (within-peak ASE); aggregate across replicates at same SNP; combine with cis-caQTL across cohort.

### RASQUAL -- LD computation requirement

**Trigger:** RASQUAL exits or returns NA for a feature.

**Mechanism:** RASQUAL estimates genotype/allelic correlation internally from the tabix-streamed VCF (no external LD matrix is needed or accepted). Failures instead come from the `-l` (testing SNP) and `-m` (feature SNP) counts not matching the SNPs actually present in the cis-window, or a feature with zero fSNPs.

**Fix:** Compute `-l`/`-m` from the actual VCF window rather than a fixed guess, and skip features with no feature SNPs; there is no LD-precompute step.

### WASP -- Phased vs unphased genotypes

**Trigger:** Using unphased genotypes for ASE.

**Mechanism:** ASE requires knowing which allele is on which haplotype to assign reads. Unphased het sites ambiguously assign reads.

**Fix:** Phase genotypes with SHAPEIT5, BEAGLE 5.4, or whatshap (read-based) before running WASP/ASE counter.

### Cohort-level caQTL without allelic info

**Trigger:** MatrixEQTL on peak counts without ASE.

**Mechanism:** MatrixEQTL maps cohort-level associations; misses cis-mode that ASE captures within individual.

**Symptom:** Power to detect caQTL is low (typical N=50-100 cohort gives ~hundreds of caQTLs vs ASE-augmented can give thousands).

**Fix:** Use RASQUAL (joint total + ASE) when N <= 100; or combine MatrixEQTL with separate ASE per individual.

### Read-deep peak coverage required

**Trigger:** Per-peak coverage < 30 reads at SNP site.

**Mechanism:** Binomial test power at p=0.5, n=30 yields detectable shifts only at |delta_p| >= 0.2.

**Fix:** Pool replicates if available; or restrict to peaks with sufficient coverage; or aggregate to per-individual peak-level ASE rather than per-SNP.

## Decision Tree by Setting

| Setting | Recommended pipeline |
|---------|---------------------|
| Single individual, ATAC + genotypes | WASP + GATK ASEReadCounter -> per-SNP and per-peak ASE; within-peak aggregation |
| Cohort N >= 100, want caQTL | WASP + GATK + MatrixEQTL on peak counts; supplement with ASE for cis-effects |
| Cohort N = 20-100 | WASP + RASQUAL (joint total + ASE) for max power |
| Cohort with no genotypes | QuASAR (infers genotypes from data) |
| Validating GWAS variant function | Look up het samples in cohort; aggregate ASE at the variant; ASB ratio |
| Trios or quartets | Per-trio phasing then ASE per individual |
| iPSC line genotype validation | Single-individual ASE at known SNPs |

## Cohort caQTL Pipeline

**Goal:** Build cohort-level chromatin QTLs by combining WASP-corrected per-sample counts with cis-genotype association.

**Approach:** WASP-correct each BAM, build consensus peakset, count reads in peaks per sample, then test cis-genotype association via MatrixEQTL (cohort) or RASQUAL (joint total + allelic).

```bash
# 1. WASP-correct each individual's BAM
for sample in $(cat samples.txt); do
    bash wasp_pipeline.sh $sample.bam $sample.vcf
done

# 2. Build consensus peakset (atac-seq/consensus-peakset)
# 3. Count reads in peaks per sample (featureCounts)
featureCounts -F SAF -a consensus.saf -o counts.tsv -p --countReadPairs *.wasp.bam

# 4. Cohort-level caQTL via MatrixEQTL
# (see R script in examples/)

# 5. Per-individual ASE (RASQUAL alternative)
# Parallel: gatk ASEReadCounter per sample, then merge for QuASAR meta-analysis
```

## Within-Peak ASE Aggregation

**Goal:** Boost per-SNP ASE power by pooling allele counts across heterozygous SNPs in the same peak.

**Approach:** Map each het SNP to its containing peak, sum REF and ALT counts per peak, run pooled binomial test against 50:50, apply BH FDR, and threshold on effect size.

```python
import pandas as pd, numpy as np
from scipy import stats

# Per-SNP allele counts at heterozygous sites
ase = pd.read_csv('sample.ase_counts.tsv', sep='\t')
ase = ase.rename(columns={'refCount': 'REF', 'altCount': 'ALT'})
ase['totalCount'] = ase['REF'] + ase['ALT']

# Map each SNP to its containing peak
ase['peak'] = map_snps_to_peaks(ase, 'consensus_peaks.bed')

# Aggregate within peak: pooled binomial test
def peak_ase(group):
    ref = group['REF'].sum()
    total = group['totalCount'].sum()
    if total < 30: return pd.Series({'ref_frac': np.nan, 'p_value': np.nan})
    p = stats.binomtest(ref, total, p=0.5).pvalue
    return pd.Series({'ref_frac': ref / total, 'p_value': p, 'snp_count': len(group)})

peak_ase_df = ase.groupby('peak').apply(peak_ase)
peak_ase_df['adj_p'] = stats.false_discovery_control(peak_ase_df['p_value'].fillna(1.0))
sig_ase = peak_ase_df[(peak_ase_df['adj_p'] < 0.05) & (abs(peak_ase_df['ref_frac'] - 0.5) >= 0.2)]
```

|ref_frac - 0.5| >= 0.2 is a 30:70 effect; smaller imbalances are detectable but biologically minor. Per-peak SNP count >= 2 strengthens the call.

## RASQUAL Joint Modeling

RASQUAL uses a non-standard CLI: it reads the VCF from stdin via tabix and uses single-letter flags. The canonical invocation pattern is:

**Goal:** Combine total accessibility counts and allele-specific counts into one joint cis-caQTL test per peak.

**Approach:** Pre-build binary count and offset files via rasqualTools, then iterate per-feature, streaming the cis-window VCF through tabix into RASQUAL with feature coordinates and SNP counts.

```bash
# 1. Pre-compute genotype offsets and binary count files (rasqualTools R package)
# (see rasqualTools::saveRasqualMatrices; produces .bin files for -y, -k)

# 2. Per-feature (peak) RASQUAL call: pipe tabix VCF in via stdin
# Per-line meaning: feature name, chromosome, start, end, n_testing_SNPs (cis-window rSNP candidates), n_feature_SNPs (fSNPs in the peak)
while IFS=$'\t' read -r name chr start end n_testing n_feature; do
    tabix cohort.vcf.gz $chr:$((start-500000))-$((end+500000)) | \
        rasqual -y counts.bin \
                -k offsets.bin \
                -n $N_SAMPLES \
                -j $FEATURE_INDEX \
                -l $n_testing -m $n_feature \
                -s $start -e $end \
                -f $name \
                > $name.rasqual.txt
done < features.tsv
```

`-y` is the binary count file; `-k` is the binary size-factor / offset file (both produced by `rasqualTools::saveRasqualMatrices` from R); `-j` is the row index of the feature; `-l` is the number of testing SNPs (cis-window rSNP candidates) and `-m` the number of feature SNPs (fSNPs overlapping the peak). RASQUAL does NOT use `--features`, `--counts`, `--vcf` flags; those are common in newer caQTL wrappers but not in stock RASQUAL.

For modern usage, the `rasqualTools` R wrapper (Kumasaka GitHub) handles this orchestration. Verify against `rasqual --help` because the flag set is unusual.

RASQUAL output includes joint p-values, total-only and ASE-only sub-tests; the joint test typically gains 1.5-3x power vs MatrixEQTL alone.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| GATK ASE shows REF > ALT systematically | WASP not run | Re-run with WASP filter; bias fixed |
| RASQUAL joint p << ASE-only p | Cis effect dominated by allelic component | Confirms cis-regulatory mechanism |
| ASE detected at SNP not in peak | Possible coding splice or 3' UTR effect | Check annotation; may not be regulatory |
| Cohort caQTL doesn't replicate per-individual ASE | Trans effect or technical artifact | Consider trans-effect; or single-individual outliers |

**Operational rule for high-confidence reporting:** A cis-regulatory variant must show (a) WASP-filtered allelic imbalance with adjusted p < 0.05 and effect size >= 0.2, AND (b) cohort-level caQTL p < 1e-5 (or RASQUAL joint p < 1e-5), AND (c) accessibility peak overlap. Validation against MPRA or CRISPRi-FlowFISH increases confidence.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Reference allele systematically over-represented | WASP not run | Mandatory WASP filter |
| ASE per-SNP coverage too low | Sparse coverage; rare alleles | Aggregate across SNPs in peak; or filter SNPs with allele freq 0.1-0.9 |
| RASQUAL crashes/NA on small features | Feature has no fSNPs, or `-l`/`-m` counts mismatch the VCF window | Skip features with 0 feature-SNPs; compute `-l`/`-m` from the actual window |
| caQTL replication low | Cohort-effect vs cis-effect confusion | Use RASQUAL joint or replicate ASE separately |
| Phased vs unphased confusion | Different software expectations | Phase with SHAPEIT5/BEAGLE before any ASE |
| GATK ASEReadCounter fails on multiallelic sites | Multi-allelic complications | Pre-filter VCF to biallelic only with bcftools |
| WASP runs slow | Re-alignment step is dominant | Parallelize per-chromosome; or use samtools faidx + region-based parallelism |

## References

- van de Geijn B et al 2015 Nat Methods 12:1061 (WASP; reference allele bias correction)
- Castel SE et al 2015 Genome Biol 16:195 (GATK ASEReadCounter; ASE framework)
- Kumasaka N et al 2016 Nat Genet 48:206 (RASQUAL; joint total + ASE caQTL)
- Harvey CT et al 2015 Bioinformatics 31:1235 (QuASAR)
- Buchkovich ML et al 2015 BMC Med Genomics 8:43 (reference mapping-bias correction with limited/no genotype data; AA-ALIGNER)
- Browning SR & Browning BL 2007 Am J Hum Genet 81:1084 (BEAGLE phasing)
- Patterson M et al 2015 J Comput Biol 22:498 (whatshap; read-based phasing)

## Related Skills

- atac-seq/atac-peak-calling - Generate peaks for ASE within-peak aggregation
- atac-seq/consensus-peakset - Cohort consensus peakset
- atac-seq/differential-accessibility - Cohort-level (cis + trans) differential
- atac-seq/deep-learning-atac - Predicted variant effects vs observed allelic imbalance
- atac-seq/enhancer-gene-linking - Map ASE-supported variants to target genes
- variant-calling/vcf-basics - VCF input
- variant-calling/joint-calling - Cohort genotype inputs
- phasing-imputation/haplotype-phasing - Phasing before WASP
- causal-genomics/fine-mapping - Use caQTL for fine-mapping
- population-genetics/association-testing - GWAS context
