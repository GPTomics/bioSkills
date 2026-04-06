---
name: bio-variant-normalization
description: Normalize indel representation, decompose MNPs, and split multiallelic variants using bcftools norm. Use when comparing variants from different callers, preparing VCF for database annotation, or merging VCFs from multiple sources.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, cyvcf2 0.30+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

The `--atomize` flag requires bcftools 1.17+. Earlier versions require `vt decompose_blocksub` as an alternative.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Variant Normalization

Left-align indels, decompose MNPs, and split multiallelic sites using bcftools norm.

## When Normalization is Mandatory

Not normalizing before certain operations leads to missed matches and false discordance. Normalization is required:

- **Before comparing variants from different callers.** Each caller may represent the same indel at different positions or encode MNPs differently. Without normalization, identical variants appear discordant.
- **Before database annotation.** dbSNP, ClinVar, and gnomAD store variants in canonical left-aligned, parsimonious representation. A right-aligned or non-parsimonious indel will fail to match its database entry.
- **Before merging VCF files from different sources.** `bcftools merge` matches on CHROM/POS/REF/ALT; different representations of the same variant produce duplicate entries instead of a single merged record.
- **Before any variant set operations.** Intersection (`bcftools isec`), complement, and union operations all rely on exact positional matching. Non-normalized variants silently fall through set comparisons.

Normalization is generally safe to skip only when a single caller produced all variants and no cross-file comparison or database lookup is needed.

## Why Normalize?

The same variant can be represented multiple ways:

```
chr1  100  ATCG  A      (right-aligned)
chr1  100  ATC   A      (left-aligned, parsimonious -- the canonical form)
chr1  101  TCG   T      (shifted position, different anchor base)
```

The VCF specification mandates left-aligned, parsimonious representation, but not all callers comply. Normalization enforces this canonical form.

## Recommended Normalization Pipeline

The order of operations matters. Performing these steps out of order can produce incorrect results (e.g., left-aligning a multiallelic record may normalize differently than splitting first, then left-aligning each biallelic record independently).

The correct order:

1. **Decompose MNPs** into atomic SNPs (`--atomize`)
2. **Split multiallelic** sites into biallelic records (`-m-`)
3. **Left-align and trim** against the reference (`-f reference.fa`)

Combined as a piped pipeline:

```bash
bcftools norm --atomize input.vcf.gz | \
    bcftools norm -m- | \
    bcftools norm -f reference.fa -Oz -o normalized.vcf.gz
bcftools index normalized.vcf.gz
```

For VCFs without MNPs (e.g., GATK HaplotypeCaller output, which does not emit MNPs), the atomize step can be skipped:

```bash
bcftools norm -m- input.vcf.gz | \
    bcftools norm -f reference.fa -Oz -o normalized.vcf.gz
```

A single-pass `bcftools norm -f ref.fa -m-any` is acceptable for basic use cases but does not control the decomposition order and skips MNP atomization.

## Left-Alignment

**"Normalize my VCF before comparing callers"** -> Left-align indel representations and split multiallelic sites for consistent variant comparison.

```bash
bcftools norm -f reference.fa input.vcf.gz -Oz -o normalized.vcf.gz
```

Requires reference FASTA to determine the leftmost position. The reference must be the same genome build used during variant calling; mismatches between builds silently produce wrong results even when REF alleles happen to match locally.

### Check for Normalization Issues

```bash
bcftools norm -f reference.fa -c s input.vcf.gz > /dev/null
```

Check modes (`-c`):
- `w` - Warn on mismatch (default)
- `e` - Error on mismatch
- `x` - Exclude mismatches
- `s` - Set correct REF from reference

## Multiallelic Splitting

### Split Multiallelic to Biallelic

```bash
bcftools norm -m-any input.vcf.gz -Oz -o split.vcf.gz
```

Before:
```
chr1  100  .  A  G,T  30  PASS  .  GT  1/2
```

After:
```
chr1  100  .  A  G  30  PASS  .  GT  1/0
chr1  100  .  A  T  30  PASS  .  GT  0/1
```

### Splitting Caveats

Splitting creates artificial missing information. A sample with genotype 1/2 (compound heterozygous for two different ALT alleles) becomes 0/1 in both split records. The information that both alleles were present at the same site in the same individual is lost. This has consequences for:

- **Phasing and compound heterozygosity detection.** Clinical pipelines that identify compound hets (two damaging variants on different alleles of the same gene) can misinterpret split records as independent heterozygous calls rather than co-occurring alleles at one site.
- **Allele depth (AD) interpretation.** AD values are retained per allele in each split record, but the genotype relationship between alleles at the same site is gone.
- **Population allele frequency estimation.** Splitting followed by naive frequency calculation can double-count samples at multiallelic sites.

Decision guidance:

| Downstream tool | Splitting required? | Rationale |
|----------------|-------------------|-----------|
| PLINK, PLINK2 | Yes | PLINK requires biallelic records |
| Most GWAS tools | Yes | Expect biallelic sites |
| Hail | No | Handles multiallelics natively; splitting loses information |
| bcftools csq | No | Supports multiallelic consequence calling |
| VEP | Either | Handles both; multiallelic may give richer output |
| ClinVar matching | Yes | ClinVar entries are biallelic |

When a downstream tool does not require splitting, prefer keeping multiallelic sites intact to preserve genotype relationships.

### Split Options

| Option | Description |
|--------|-------------|
| `-m-any` | Split all multiallelic sites |
| `-m-snps` | Split multiallelic SNPs only |
| `-m-indels` | Split multiallelic indels only |
| `-m-both` | Split SNPs and indels separately |
| `-m+any` | Join biallelic sites into multiallelic |
| `-m+snps` | Join biallelic SNPs |
| `-m+indels` | Join biallelic indels |
| `-m+both` | Join SNPs and indels separately |

### Join Biallelic to Multiallelic

```bash
bcftools norm -m+any input.vcf.gz -Oz -o merged.vcf.gz
```

Rejoining after analysis can restore compound heterozygosity context, but only if the split records were not independently filtered (removing one allele of a 1/2 site makes the remaining record misleading).

## Atomize Complex Variants (MNP Decomposition)

Multi-nucleotide polymorphisms (MNPs) are adjacent substitutions reported as a single record (e.g., ATG->GCA). Not all callers emit MNPs:

| Caller | Emits MNPs? | Notes |
|--------|------------|-------|
| FreeBayes | Yes | Reports MNPs and complex events natively |
| Octopus | Yes | Local haplotype-aware, emits block substitutions |
| GATK HaplotypeCaller | No | Decomposes variants during calling; may emit nearby SNPs in the same haplotype block |
| DeepVariant | Rarely | Primarily emits SNPs and indels |

Decomposing MNPs is necessary when comparing output from callers that represent them differently. Without atomization, an MNP from FreeBayes will not match the equivalent individual SNPs from GATK.

### Atomize MNPs to SNPs

```bash
bcftools norm --atomize input.vcf.gz -Oz -o atomized.vcf.gz
```

Before:
```
chr1  100  .  ATG  GCA  30  PASS
```

After:
```
chr1  100  .  A  G  30  PASS
chr1  101  .  T  C  30  PASS
chr1  102  .  G  A  30  PASS
```

**Caveat:** Decomposition loses local phasing information. The original MNP record guarantees that A->G, T->C, and G->A occur on the same haplotype. After atomization, this co-occurrence is no longer explicit. If downstream analysis requires haplotype-aware interpretation (e.g., amino acid change prediction where the codon change matters), atomization may be inappropriate -- use `bcftools csq` on the un-atomized VCF instead.

### Atomize with Old Record Tag

```bash
bcftools norm --atomize --old-rec-tag ORIGINAL input.vcf.gz -Oz -o atomized.vcf.gz
```

Preserves the original record as an INFO annotation, enabling traceability back to the pre-atomized variant.

## Fixing Reference Alleles

**Goal:** Correct or remove variants whose REF allele does not match the reference genome.

**Approach:** Use bcftools norm -c with mode s (set correct REF) or x (exclude mismatches).

### Fix Mismatches from Reference

```bash
bcftools norm -f reference.fa -c s input.vcf.gz -Oz -o fixed.vcf.gz
```

This sets REF alleles to match the reference genome. Use with caution: REF mismatches often indicate a genome build mismatch, and silently "fixing" REF may mask a liftover error rather than correcting a trivial typo.

### Exclude Mismatches

```bash
bcftools norm -f reference.fa -c x input.vcf.gz -Oz -o clean.vcf.gz
```

Removes variants where REF does not match the reference. Safer than `-c s` when the cause of mismatch is unknown.

## Remove Duplicates After Splitting

```bash
bcftools norm -d exact input.vcf.gz -Oz -o deduped.vcf.gz
```

Duplicate removal options (`-d`):
- `exact` - Remove exact duplicates (same CHROM, POS, REF, ALT)
- `snps` - Remove duplicate SNPs only
- `indels` - Remove duplicate indels only
- `both` - Remove duplicate SNPs and indels
- `all` - Remove all duplicates at the same position
- `none` - Keep duplicates (default)

## Common Workflows

### Full Normalization for Caller Comparison

**Goal:** Make VCFs from different callers directly comparable.

**Approach:** Apply the same three-step normalization pipeline to each VCF, then use set operations.

```bash
for vcf in gatk.vcf.gz freebayes.vcf.gz; do
    base=$(basename "$vcf" .vcf.gz)
    bcftools norm --atomize "$vcf" | \
        bcftools norm -m- | \
        bcftools norm -f reference.fa -Oz -o "${base}.norm.vcf.gz"
    bcftools index "${base}.norm.vcf.gz"
done

bcftools isec -p comparison gatk.norm.vcf.gz freebayes.norm.vcf.gz
```

The `isec` output directories: `0000.vcf` = private to first file, `0001.vcf` = private to second, `0002.vcf`/`0003.vcf` = shared variants from each file.

### Before Database Annotation

```bash
bcftools norm --atomize variants.vcf.gz | \
    bcftools norm -m- | \
    bcftools norm -f reference.fa -Oz -o for_annotation.vcf.gz
bcftools index for_annotation.vcf.gz
```

### Prepare for GWAS (PLINK)

**Goal:** Produce a biallelic, SNP-only, deduplicated VCF suitable for PLINK import.

**Approach:** Normalize, split, restrict to SNPs, and remove duplicates.

```bash
bcftools norm -f reference.fa -m- input.vcf.gz | \
    bcftools view -v snps | \
    bcftools norm -d exact -Oz -o gwas_ready.vcf.gz
bcftools index gwas_ready.vcf.gz
```

## cyvcf2 Normalization Check

**Goal:** Assess how many variants require normalization before running bcftools norm.

**Approach:** Iterate with cyvcf2 and count multiallelic sites and complex (MNP) variants.

```python
from cyvcf2 import VCF

def needs_normalization(variant):
    if len(variant.ALT) > 1:
        return True
    ref, alt = variant.REF, variant.ALT[0]
    if len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
        return True
    return False

total, needs_norm, multiallelic, mnps = 0, 0, 0, 0
for variant in VCF('input.vcf.gz'):
    total += 1
    if len(variant.ALT) > 1:
        multiallelic += 1
    ref, alt = variant.REF, variant.ALT[0]
    if len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
        mnps += 1
    if needs_normalization(variant):
        needs_norm += 1

print(f'Total variants: {total}')
print(f'Needing normalization: {needs_norm} ({needs_norm/total*100:.1f}%)')
print(f'  Multiallelic sites: {multiallelic}')
print(f'  MNPs: {mnps}')
```

Note: this check does not detect indels requiring left-alignment, since that requires reference context. The count is a lower bound.

## Quick Reference

| Task | Command |
|------|---------|
| Left-align indels | `bcftools norm -f ref.fa in.vcf.gz` |
| Split multiallelic | `bcftools norm -m-any in.vcf.gz` |
| Join to multiallelic | `bcftools norm -m+any in.vcf.gz` |
| Atomize MNPs | `bcftools norm --atomize in.vcf.gz` |
| Fix REF alleles | `bcftools norm -f ref.fa -c s in.vcf.gz` |
| Remove duplicates | `bcftools norm -d exact in.vcf.gz` |
| Full pipeline | `bcftools norm --atomize \| bcftools norm -m- \| bcftools norm -f ref.fa` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `REF does not match` | Wrong reference or genome build mismatch | Verify the reference FASTA matches the build used during calling |
| `not sorted` | Unsorted input | Run `bcftools sort` first |
| `duplicate records` | Same position twice after splitting | Use `-d exact` to remove |
| `--atomize` unrecognized | bcftools < 1.17 | Upgrade bcftools, or use `vt decompose_blocksub` as alternative |

## Related Skills

- variant-calling/variant-calling - Generate VCF files from alignments
- variant-calling/filtering-best-practices - Filter after normalization
- variant-calling/vcf-manipulation - Merge, intersect, and compare VCFs
- variant-calling/variant-annotation - Annotate normalized variants against databases
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller workflow (does not emit MNPs)
- variant-calling/clinical-interpretation - ClinVar lookup requires normalized representation
- alignment-files/sam-bam-basics - BAM format and reference genome handling
