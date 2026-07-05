---
name: bio-variant-calling-joint-calling
description: Joint genotype a cohort of per-sample gVCFs with GATK (HaplotypeCaller -ERC GVCF -> GenomicsDBImport or CombineGVCFs -> GenotypeGVCFs) or GLnexus for DeepVariant gVCFs, producing a squared-off sample-by-site genotype matrix. Use when deciding between joint genotyping and merging single-sample callsets (never bcftools merge as absent==hom-ref), choosing GenomicsDBImport vs CombineGVCFs by cohort size and memory, solving the N+1 problem so a new sample does not force re-calling everyone, understanding cohort rescue of low-coverage het sites, handling the spanning-deletion star allele and GQ/PL recomputation at the joint step, scaling to biobank cohorts by interval sharding, or picking DeepVariant+GLnexus over the GATK path on throughput. Not for single-sample calling (see variant-calling/gatk-variant-calling) or VQSR/hard-filter mechanism (see variant-calling/filtering-best-practices).
tool_type: cli
primary_tool: GATK
---

## Version Compatibility

Reference examples tested with: GATK 4.5+, GLnexus 1.4+, bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Joint Calling

**"Joint genotype my cohort samples"** -> Combine per-sample gVCFs into a single cohort callset with consistent genotyping across all sites, enabling cohort filtering and population-level analysis.
- CLI (GATK): `gatk HaplotypeCaller -ERC GVCF` -> `gatk GenomicsDBImport` (or `CombineGVCFs`) -> `gatk GenotypeGVCFs`
- CLI (DeepVariant cohorts): `deepvariant --output_gvcf` per sample -> `glnexus_cli --config DeepVariantWGS`

## The Governing Principle: Genotype Jointly, Never Merge Callsets

Joint genotyping is not the same operation as merging single-sample VCFs, and confusing the two silently corrupts a cohort callset. Two facts drive every decision below:

- **Joint genotyping rescues low-coverage het sites.** Evidence is borrowed across samples: when one carrier has a confident variant, the cohort allele frequency raises the Bayesian prior at that exact site for every other sample, so a second sample with only 2-3 supporting reads (which per-sample would fall below threshold) is rescued into a confident genotype. This "borrowing information" mechanism is strongest at low coverage and is why a cohort callset is more sensitive than N independent callsets (DePristo 2011 *Nat Genet* 43:491; Poplin 2018 *bioRxiv* 201178).
- **Joint genotyping produces a squared-off matrix** - a genotype at every variant site for every sample. The reference-confidence `<NON_REF>` allele in each gVCF lets GenotypeGVCFs distinguish confident homozygous reference (`0/0`) from no-data/no-call (`./.`) at a site another sample carries.

The decision this forces: **never `bcftools merge` single-sample callsets as if an absent record means hom-ref.** A single-sample VCF omits sites where that sample looked reference, so a naive merge fills those cells with `./.` (missing), NOT `0/0` - a sample genuinely hom-ref and a sample never assessed become indistinguishable, and downstream allele frequencies and association tests are wrong. `./.` != `0/0` is the load-bearing distinction (see variant-calling/vcf-manipulation for merge semantics and variant-calling/vcf-basics for the genotype-field grammar). Genotype from gVCFs so every cell is filled from evidence, not assumption.

## Why Joint Calling Matters

Single-sample calling discards cross-sample evidence that is critical for accurate genotyping:

- **Statistical power from shared evidence** - A site with 2 alt reads in one sample is borderline and would typically be missed. If 50 other samples in the cohort also show 2 alt reads at that site, the evidence is overwhelming and the variant is clearly real. Joint calling aggregates this weak-per-sample signal into strong cohort-level evidence.
- **Genotype refinement via cohort priors** - Individual genotype likelihoods are combined with cohort allele frequencies as a Bayesian prior. A heterozygous call at a common variant site (AF=0.3) receives more support than the same call at a site with no other carriers. This prior dramatically improves accuracy for low-coverage samples.
- **Consistent site representation** - All samples are genotyped at the same sites, producing homozygous-reference calls where applicable. Without joint calling, a missing genotype is ambiguous: it could mean homozygous-reference or simply insufficient coverage. This "missing = reference" assumption is a common source of false negatives in downstream analysis.
- **Cohort filtering eligibility** - Variant quality score recalibration (VQSR) and its successor VETS operate on the whole-cohort variant distribution and need cohort-scale variant counts: a single deep WGS genome supplies enough, but exomes are variant-poor so the ~30-sample floor applies to exomes/panels, not WGS (VQSR's Gaussian mixture needs enough variants to fit), so filtering is inherently a cohort operation, not a per-sample one (see Step 4).

## The N+1 Problem and Why gVCFs Solve It

Naive joint calling re-visits every BAM whenever the cohort changes: adding one genome forces re-calling all N. The gVCF workflow decouples expensive per-sample **discovery** (local assembly + PairHMM likelihoods, captured once per sample in the gVCF) from cheap cohort-wide **genotyping**. Adding sample N+1 then requires only generating that one gVCF plus re-running the cheap consolidation and GenotypeGVCFs - the assembly work for the existing N is never repeated. The gVCF is the reusable intermediate; GenomicsDB workspaces can even be updated in place (`--genomicsdb-update-workspace-path`). GATK frames this as decoupling "the initial identification of potential variant sites from the genotyping step, which is the only part that really needs to be done jointly" (see variant-calling/gatk-variant-calling for per-sample gVCF generation).

## Cohort Size Decision Table

| Cohort Size | Approach | Notes |
|---|---|---|
| <100 | CombineGVCFs or GenomicsDB | Either works; CombineGVCFs is simpler to manage |
| 100-10,000 | GenomicsDB + GenotypeGVCFs | Standard GATK Best Practices; shard by chromosome |
| 10,000-100,000 | GATK Biggest Practices | Heavily sharded and parallelized across intervals |
| >100,000 | DeepVariant + GLnexus, or Hail VDS | GATK becomes unwieldy at this scale; purpose-built tools required |

## GenomicsDBImport vs CombineGVCFs

Both produce a combined object that GenotypeGVCFs consumes; they differ in how they store it and how they scale.

| | GenomicsDBImport | CombineGVCFs |
|---|---|---|
| Storage | GenomicsDB workspace on a TileDB array backend; transposes sample-centric gVCFs into a locus-centric sparse 2-D array (samples x positions) | Pure-Java hierarchical merge into a single combined gVCF |
| Scaling | Best when N is large; the locus-centric transpose is what makes per-locus genotyping fast at scale | Fails when N grows - memory-hungry and slow; recommended only as a small-cohort fallback |
| Portability | Workspace is not a plain gVCF; genotype via `gendb://` | Output is a plain gVCF, portable and inspectable |
| Incremental | Add new samples with `--genomicsdb-update-workspace-path` (the N+1 win in practice) | No incremental mode; re-run over all samples |
| Best when | >100 samples, sharded by interval, biobank scale | <100 samples, or a small family/trio where simplicity wins |

Memory landmine specific to GenomicsDBImport: the heavy lifting runs in native C/C++ (TileDB), so cap the JVM heap (`--java-options -Xmx`) at ~80-90% of RAM. Over-allocating the JVM starves the native layer and causes a native out-of-memory failure that looks unrelated to heap size.

## Workflow Overview

```
Sample BAMs
    │
    ├── HaplotypeCaller (per-sample, -ERC GVCF)
    │   └── sample1.g.vcf.gz, sample2.g.vcf.gz, ...
    │
    ├── CombineGVCFs or GenomicsDBImport
    │   └── Combine into cohort database
    │
    ├── GenotypeGVCFs
    │   └── Joint genotyping
    │
    └── VQSR or Hard Filtering
        └── Final VCF
```

## Step 1: Per-Sample gVCF Generation

```bash
# Generate gVCF for each sample
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample1.bam \
    -O sample1.g.vcf.gz \
    -ERC GVCF

# With intervals (faster)
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample1.bam \
    -O sample1.g.vcf.gz \
    -ERC GVCF \
    -L intervals.bed
```

### Batch Processing

```bash
# Process all samples
for bam in *.bam; do
    sample=$(basename $bam .bam)
    gatk HaplotypeCaller \
        -R reference.fa \
        -I $bam \
        -O ${sample}.g.vcf.gz \
        -ERC GVCF &
done
wait
```

## Step 2a: CombineGVCFs (Small Cohorts)

For <100 samples:

```bash
gatk CombineGVCFs \
    -R reference.fa \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    -V sample3.g.vcf.gz \
    -O cohort.g.vcf.gz
```

### From Sample Map

```bash
# Create sample map file
# sample1    /path/to/sample1.g.vcf.gz
# sample2    /path/to/sample2.g.vcf.gz

ls *.g.vcf.gz | while read f; do
    echo -e "$(basename $f .g.vcf.gz)\t$f"
done > sample_map.txt

# Combine with -V for each
gatk CombineGVCFs \
    -R reference.fa \
    $(cat sample_map.txt | cut -f2 | sed 's/^/-V /') \
    -O cohort.g.vcf.gz
```

## Step 2b: GenomicsDBImport (Large Cohorts)

For >100 samples, use GenomicsDB:

```bash
# Create sample map
ls *.g.vcf.gz | while read f; do
    echo -e "$(basename $f .g.vcf.gz)\t$f"
done > sample_map.txt

# Import to GenomicsDB (per chromosome for parallelism)
gatk GenomicsDBImport \
    --sample-name-map sample_map.txt \
    --genomicsdb-workspace-path genomicsdb_chr1 \
    -L chr1 \
    --reader-threads 4

# Or all chromosomes
for chr in {1..22} X Y; do
    gatk GenomicsDBImport \
        --sample-name-map sample_map.txt \
        --genomicsdb-workspace-path genomicsdb_chr${chr} \
        -L chr${chr} &
done
wait
```

### Update GenomicsDB with New Samples

```bash
gatk GenomicsDBImport \
    --genomicsdb-update-workspace-path genomicsdb_chr1 \
    --sample-name-map new_samples.txt \
    -L chr1
```

### GenomicsDB Critical Caveats

GenomicsDB is powerful but has sharp edges that can cause data loss or silent failures:

- **No sample replacement** - Existing samples cannot be updated or overwritten. Only new samples with different names can be added. To fix a sample, the entire workspace must be recreated.
- **Intervals locked at import time** - The genomic intervals specified during the initial import cannot be changed on incremental updates. Adding new regions requires reimporting from scratch.
- **Fragment accumulation** - Each incremental batch creates a new database fragment. After thousands of incremental additions, file handle exhaustion becomes likely. Run `--consolidate` periodically to merge fragments.
- **Corruption risk on failed adds** - A failed incremental import can leave the datastore in an inconsistent state. Always backup the workspace directory before running `--genomicsdb-update-workspace-path`.
- **Batch size for memory** - Set `--batch-size 50` to control memory consumption. The default is `0`, which loads ALL samples in a single batch (maximum memory); a finite batch size trades a little speed for a bounded heap, so set it explicitly for large cohorts. Larger batches load more gVCFs simultaneously and can exhaust heap space.

## Step 3: GenotypeGVCFs

### From Combined gVCF

```bash
gatk GenotypeGVCFs \
    -R reference.fa \
    -V cohort.g.vcf.gz \
    -O cohort.vcf.gz
```

### From GenomicsDB

```bash
gatk GenotypeGVCFs \
    -R reference.fa \
    -V gendb://genomicsdb_chr1 \
    -O chr1.vcf.gz

# All chromosomes
for chr in {1..22} X Y; do
    gatk GenotypeGVCFs \
        -R reference.fa \
        -V gendb://genomicsdb_chr${chr} \
        -O chr${chr}.vcf.gz &
done
wait

# Merge chromosomes
bcftools concat chr{1..22}.vcf.gz chrX.vcf.gz chrY.vcf.gz \
    -Oz -o cohort.vcf.gz
```

### What GenotypeGVCFs Recomputes (and Why It Is Not a Copy)

GenotypeGVCFs re-derives genotypes jointly from the stored per-sample PL vectors (including the `<NON_REF>` likelihood) under a Bayesian model; it does not simply copy per-sample genotypes into a wider file. Two consequences matter when reading the output:

- **GQ and PL are recomputed against the finalized cohort allele set.** Once the real ALT alleles are known cohort-wide, the `<NON_REF>` likelihood is redistributed onto them and PLs are recomputed; GQ is then the difference of the two smallest PLs. A sample's genotype/GQ in the joint VCF can therefore differ from what its single-sample gVCF implied - this is the rescue mechanism working, not a bug.
- **The allele-frequency prior comes from `--heterozygosity`** (expected theta, ~0.001 for humans; verify in-tool) and `--indel-heterozygosity`, folding the cohort allele count into each sample's posterior. `--stand-call-conf` (~30; verify in-tool) drops sites below that QUAL.

### Multiallelics and the Spanning-Deletion `*` Allele

Joint genotyping across a cohort surfaces two representation issues absent from single-sample calling:

- **`--max-alternate-alleles`** caps the ALT alleles genotyped per site (most-supported kept; confirm the default with `gatk GenotypeGVCFs --help` for the installed version). Genotyping cost scales roughly exponentially in ALT count, so GATK caps it and discourages raising it; `--max-genotype-count` similarly bounds genotype configurations. Highly multiallelic sites are also where GenotypeGVCFs can blow past very large RAM at scale.
- **The `*` spanning/overlapping-deletion allele** (VCF 4.3 reserved) appears at a variant position that falls *inside an upstream deletion carried by some samples*. It means "for a sample carrying the upstream deletion, these bases are deleted/absent" - not reference, not the local ALT. Such a sample genotypes as `*/A` or `*/*`, which correctly keeps deletion-carriers from being called spuriously hom-ref at the interior site. Downstream tools must special-case it: VEP/SnpEff have no ref/alt sequence to predict a consequence on, and `bcftools norm` decomposition often splits or filters it, so it is a frequent source of annotation surprises after joint genotyping.

### With Allele-Specific Annotations

For larger cohorts where multiallelic sites are common, allele-specific annotations allow VQSR to evaluate each allele independently rather than penalizing a good allele because a co-occurring allele is poor:

```bash
gatk GenotypeGVCFs \
    -R reference.fa \
    -V gendb://genomicsdb \
    -O cohort.vcf.gz \
    -G StandardAnnotation \
    -G AS_StandardAnnotation
```

When allele-specific annotations are present, use `-AS` mode in VariantRecalibrator and ApplyVQSR for allele-level filtering.

## Step 4: Filtering Is a Cohort Operation

Filtering the joint VCF is a whole-cohort step, not a per-sample one, and it must run *after* joint genotyping. VQSR (and its GATK successor VETS) fit a model to the cohort-wide distribution of site annotations - VQSR's Gaussian mixture needs enough variants and enough overlap with the truth resources to converge, which is why it is unreliable on a single exome or a small panel. This is the decision:

| Cohort | Filter | Why |
|---|---|---|
| Single deep WGS, or a joint cohort (~30+ exomes) | VQSR, or VETS (isolation forest, VQSR's successor) | Enough variants to fit a stable multivariate model -- one WGS genome supplies millions of sites, but exomes are variant-poor so the ~30-sample floor applies to exomes/panels, not WGS; pass `-AS` when allele-specific annotations were emitted so each allele at a multiallelic site is filtered independently |
| Single exome, gene panel, or too few variants | Hard filters | GMM will not converge on too few variants; use fixed per-annotation thresholds instead |

The full VariantRecalibrator/ApplyVQSR/VETS invocations, training resources, tranche levels, and hard-filter thresholds live in variant-calling/filtering-best-practices - this skill does not duplicate that mechanism. A minimal hard-filter fallback for a small cohort:

```bash
gatk VariantFiltration \
    -R reference.fa \
    -V cohort.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    -O cohort.filtered.vcf.gz
```

## Batch Effects in Joint Calling

Joint genotyping mitigates most batch effects because it re-evaluates genotype likelihoods across all samples simultaneously, recalibrating quality scores against the full cohort distribution. However, certain batch effects persist through joint calling because they affect the underlying read data, not the genotyping model:

- **Different library prep protocols** - PCR-free vs PCR-based libraries produce different duplicate and error profiles
- **Different capture kits (WES)** - Exome kits target different regions; sites outside the intersection have systematically missing data in some batches
- **Significantly different coverage distributions** - 10x WGS samples mixed with 30x samples will have systematically different genotype quality at heterozygous sites
- **Different reference genome versions** - Mixing GRCh37 and GRCh38 alignments is not valid; all samples must use the same reference
- **Mixing WGS and WES** - Fundamentally different coverage profiles; off-target WES regions behave like very-low-coverage WGS

Mitigation: process all samples through an identical upstream pipeline (same aligner, same duplicate marking, same BQSR resources). If batches are unavoidable, include batch as a covariate in downstream association or differential analyses.

## When to Re-genotype

| Scenario | Action | Rationale |
|---|---|---|
| Adding new samples | Re-genotype (GenomicsDB incremental add + GenotypeGVCFs on full database) | New samples change cohort allele frequencies, improving all genotype calls |
| Changing reference genome | Full reprocess from alignment | gVCF coordinates are reference-specific |
| Updating caller version | Optional but recommended for consistency | Different caller versions may produce different quality scores; mixing versions adds noise |
| Adding new genomic intervals | Reimport from scratch | GenomicsDB intervals are locked at initial import; incremental update cannot expand them |

## DeepVariant + GLnexus Alternative

GLnexus is a scalable gVCF-merging/joint-genotyping engine (originally rocksdb-backed) that grows a cohort **incrementally** as samples are added, avoiding the full-cohort reprocessing the GenomicsDBImport + GenotypeGVCFs path requires. Yun et al. tuned its quality thresholds for DeepVariant output; the optimized presets ship in GLnexus v1.2.2+ as `DeepVariantWGS` (whole-genome) and `DeepVariantWES` (whole-exome). The original method was validated at cohort scale on ~50,000 exomes (Lin 2018 *bioRxiv* 343970). Prefer this path when DeepVariant is the caller (see variant-calling/deepvariant).

```bash
# Step 1: Run DeepVariant per sample to produce gVCFs
run_deepvariant --model_type=WGS \
    --ref=reference.fa --reads=sample.bam \
    --output_vcf=sample.vcf.gz --output_gvcf=sample.g.vcf.gz

# Step 2: Joint call with GLnexus (pre-tuned configs encode DeepVariant-tuned GQ + multiallelic handling)
# GLnexus emits BCF on stdout; pipe through bcftools to bgzip a VCF
glnexus_cli --config DeepVariantWGS --bed intervals.bed \
    sample1.g.vcf.gz sample2.g.vcf.gz ... | bcftools view - | bgzip -c > cohort.vcf.gz
```

### DeepVariant+GLnexus vs GATK GenomicsDB (representative numbers)

From the GLnexus benchmark (Yun et al. 2020 *Bioinformatics* 36:5582, GIAB 40x WGS and the 2,504-sample 1000 Genomes cohort). These are one study's figures at specific versions/coverage - treat as representative, not universal:

| Metric | DeepVariant + GLnexus | GATK (VQSR) |
|---|---|---|
| SNP F1 error | 0.07% | 1.23% |
| Indel F1 error | 1.14% | 2.92% |
| Cohort Mendelian violation rate | 1.7% | 5.0% |
| Cohort merge time, chr22 (2,504 samples) | 0.84 h | 6.83 h (GenomicsDBImport + GenotypeGVCFs) |
| Cohort gVCF footprint | 2.20 TB | 15.16 TB |

The throughput gap (GLnexus merge ~8x faster, DeepVariant gVCFs ~7x smaller on disk) is the practical reason large DeepVariant cohorts use GLnexus rather than routing DeepVariant gVCFs through GenotypeGVCFs.

## Scaling to Biobank Cohorts (tens of thousands+)

Naive GenotypeGVCFs does not scale to tens of thousands of samples: I/O and per-site QUAL computation dominate, GenotypeGVCFs can exceed very large RAM at highly multiallelic sites, and single-interval GenomicsDB workspaces plus fragment proliferation and open-file-descriptor limits become the recurring failures. The scaling levers:

- **Shard by interval.** Run one GenomicsDBImport + GenotypeGVCFs per chromosome (or finer) in parallel, then `bcftools concat`. A `--sample-name-map` file (sample<TAB>path, one per line) is mandatory at this scale - passing thousands of `-V` arguments is unmanageable and slow.
- **ReblockGVCF then GnarlyGenotyper (GATK "Biggest Practices").** ReblockGVCF drops uncalled/low-GQ alleles and re-bands reference blocks, shrinking files and merge time; GnarlyGenotyper approximates QUAL from a precomputed `QUALapprox` INFO field without iterating over all genotypes, the dominant cost saver above ~tens of thousands of samples. Broad switches production to reblocking around ~2,000 samples for cost. gnomAD v2.1 aggregated its callset in Hail and filtered with a custom random-forest model rather than VQSR (Karczewski 2020 *Nature* 581:434); later releases (v3+) ingest gVCFs directly via the Hail sparse combiner.
- **The DRAGEN / GLnexus route.** At biobank scale many projects avoid the vanilla GATK path entirely: DeepVariant + GLnexus (throughput above), or Illumina DRAGEN's integrated map-align-call engine. Verify a project's exact production pipeline rather than assuming it is GATK joint calling.

## Complete Pipeline Script

**Goal:** Run the full joint calling workflow from BAMs to filtered cohort VCF.

**Approach:** Generate per-sample gVCFs, import into GenomicsDB, joint genotype, then index and compute statistics.

```bash
#!/bin/bash
set -euo pipefail

REFERENCE=$1
OUTPUT_DIR=$2
THREADS=16

mkdir -p $OUTPUT_DIR/{gvcfs,genomicsdb,vcfs}

echo "=== Step 1: Generate gVCFs ==="
for bam in data/*.bam; do
    sample=$(basename $bam .bam)
    gatk HaplotypeCaller \
        -R $REFERENCE \
        -I $bam \
        -O $OUTPUT_DIR/gvcfs/${sample}.g.vcf.gz \
        -ERC GVCF &

    # Limit parallelism
    while [ $(jobs -r | wc -l) -ge $THREADS ]; do sleep 1; done
done
wait

echo "=== Step 2: Create sample map ==="
ls $OUTPUT_DIR/gvcfs/*.g.vcf.gz | while read f; do
    echo -e "$(basename $f .g.vcf.gz)\t$(realpath $f)"
done > $OUTPUT_DIR/sample_map.txt

echo "=== Step 3: GenomicsDBImport ==="
gatk GenomicsDBImport \
    --sample-name-map $OUTPUT_DIR/sample_map.txt \
    --genomicsdb-workspace-path $OUTPUT_DIR/genomicsdb \
    -L intervals.bed \
    --reader-threads 4

echo "=== Step 4: Joint genotyping ==="
gatk GenotypeGVCFs \
    -R $REFERENCE \
    -V gendb://$OUTPUT_DIR/genomicsdb \
    -O $OUTPUT_DIR/vcfs/cohort.vcf.gz

echo "=== Step 5: Index ==="
bcftools index -t $OUTPUT_DIR/vcfs/cohort.vcf.gz

echo "=== Statistics ==="
bcftools stats $OUTPUT_DIR/vcfs/cohort.vcf.gz > $OUTPUT_DIR/vcfs/cohort_stats.txt

echo "=== Complete ==="
echo "Joint VCF: $OUTPUT_DIR/vcfs/cohort.vcf.gz"
```

## Tips

### Memory for Large Cohorts

```bash
# Increase Java heap for GenotypeGVCFs (default 4g is often insufficient for >500 samples)
gatk --java-options "-Xmx64g" GenotypeGVCFs ...

# For GenomicsDBImport, --batch-size controls how many gVCFs are loaded simultaneously
gatk GenomicsDBImport --batch-size 50 ...
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Merged cohort has `./.` where samples are truly hom-ref; allele frequencies look wrong | `bcftools merge` of single-sample VCFs treats absent records as missing, not `0/0` | Genotype from gVCFs (GenotypeGVCFs/GLnexus) so every cell is filled from evidence; never merge single-sample callsets for a cohort matrix |
| GenomicsDBImport dies with a native/out-of-memory error despite a large `-Xmx` | Over-allocated JVM heap starves the native TileDB layer | Cap `--java-options -Xmx` at ~80-90% of RAM; the heavy lifting is native C/C++ |
| Cannot update an existing sample in GenomicsDB | GenomicsDB has no sample replacement; only new sample names can be added | Recreate the workspace to fix a sample; use `--genomicsdb-update-workspace-path` only to ADD |
| `--genomicsdb-update-workspace-path` cannot expand to new regions | Intervals are locked at initial import | Reimport from scratch to add genomic intervals |
| `*` alleles / genotypes like `*/A` break VEP/SnpEff or vanish after `bcftools norm` | Spanning-deletion symbolic allele has no ref/alt sequence to annotate | Expected after joint genotyping; special-case or split `*` records before annotation |
| VariantRecalibrator fails to converge or errors | Too few variants for the Gaussian mixture (a single exome/panel, not a single WGS) | Fall back to hard filters for exomes/panels below ~30 samples (see filtering-best-practices) |
| Fewer ALT alleles than expected at a multiallelic site | `--max-alternate-alleles` dropped the least-supported alts | Raise cautiously (cost scales ~exponentially); confirm the default with `--help` |

## Related Skills

- variant-calling/gatk-variant-calling - Single-sample HaplotypeCaller and per-sample gVCF generation (the N+1 intermediate)
- variant-calling/deepvariant - DeepVariant caller feeding the GLnexus pathway
- variant-calling/filtering-best-practices - VQSR/VETS and hard-filter mechanism (not duplicated here)
- variant-calling/vcf-manipulation - Merge/subset semantics and why single-sample merge != joint genotyping
- variant-calling/vcf-basics - Genotype-field grammar (`./.` vs `0/0`, the `*` allele)
- population-genetics/plink-basics - Population analysis of joint calls
- workflows/fastq-to-variants - End-to-end germline pipeline

## References

- DePristo MA, Banks E, Poplin R, et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. *Nature Genetics* 43(5):491-498 (2011). DOI 10.1038/ng.806.
- Poplin R, Ruano-Rubio V, DePristo MA, et al. Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv* 201178 (2018). DOI 10.1101/201178. Preprint; GATK's recommended cite for the GVCF/reference-confidence + joint-genotyping methodology.
- Yun T, Li H, Chang P-C, Lin MF, Carroll A, McLean CY. Accurate, scalable cohort variant calls using DeepVariant and GLnexus. *Bioinformatics* 36(24):5582-5589 (2020). DOI 10.1093/bioinformatics/btaa1081.
- Lin MF, Rodeh O, Penn J, et al. GLnexus: joint variant calling for large cohort sequencing. *bioRxiv* 343970 (2018). DOI 10.1101/343970. Preprint (original GLnexus method).
- Karczewski KJ, Francioli LC, Tiao G, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature* 581(7809):434-443 (2020). DOI 10.1038/s41586-020-2308-7.
