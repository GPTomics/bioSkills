---
name: bio-variant-normalization
description: Left-align and trim indels to parsimonious canonical form, decompose MNPs (atomize), and split multiallelic variants with bcftools norm. Use when comparing variants across callers or cohorts, preparing a VCF for database annotation or ClinVar/dbSNP matching, merging VCFs, reconciling vt-vs-bcftools representation discordance, or resolving the VCF-left-align vs HGVS-3'-rule clash.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, cyvcf2 0.30+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

The `--atomize`/`--old-rec-tag` flags require bcftools 1.12+ (`--keep-sum` requires 1.11+). Earlier versions require `vt decompose_blocksub` as an alternative.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Variant Normalization

Left-align indels, decompose MNPs, and split multiallelic sites using bcftools norm.

**"Put my VCF in canonical form before comparing, annotating, or merging"** -> Enforce one representation per biological event so tuple-keyed operations recognize the same variant across sources.
- CLI: `bcftools norm -m-any -f ref.fa` (split multiallelic, then left-align + parsimony each biallelic)
- Alternatives: `vt normalize` + `vt decompose`; GATK `LeftAlignAndTrimVariants`

## The Normalized Form (governing principle)

A variant is **normalized** if and only if it is BOTH (Tan, Abecasis & Kang 2015):

1. **Parsimonious (minimal representation).** No base can be trimmed from either end of REF/ALT without creating a zero-length allele or changing the variant. Indels keep exactly ONE anchor base (VCF forbids empty alleles), so a deletion of `A` is written `REF=CA ALT=C`, not `REF=A ALT=`.
2. **Left-aligned.** POS cannot be shifted further left while keeping the allele lengths and the reference-relative meaning unchanged.

Both are required: parsimony alone leaves an indel positionally ambiguous inside a repeat; left-alignment alone can leave redundant shared bases. Together they define a single canonical `(CHROM, POS, REF, ALT)` key, and normalization is **idempotent** -- normalizing an already-normalized VCF changes nothing. That idempotence is what makes the key safe to join, merge, and match on.

The algorithm (as in `vt normalize` and `bcftools norm`):
- **Right-trim:** while all alleles end in the same base AND all have length >= 2, drop the last base of each.
- **Left-extend:** while any allele has length 0, prepend the reference base immediately to the left to every allele and decrement POS.
- **Left-trim:** while all alleles start with the same base AND all have length >= 2, drop the first base of each and increment POS.

The right-trim-then-left-roll loop is exactly what shifts an indel through a repeat to its leftmost equivalent position.

## When Normalization is Mandatory

Not normalizing before certain operations leads to missed matches and false discordance. Normalization is required:

- **Before comparing variants from different callers.** Each caller may represent the same indel at different positions or encode MNPs differently. Without normalization, identical variants appear discordant.
- **Before database annotation.** dbSNP, ClinVar, and gnomAD store variants in canonical left-aligned, parsimonious representation. A right-aligned or non-parsimonious indel will fail to match its database entry.
- **Before merging VCF files from different sources.** `bcftools merge` matches on CHROM/POS/REF/ALT; different representations of the same variant produce duplicate entries instead of a single merged record.
- **Before any variant set operations.** Intersection (`bcftools isec`), complement, and union operations all rely on exact positional matching. Non-normalized variants silently fall through set comparisons.

Normalization is generally safe to skip only when a single caller produced all variants and no cross-file comparison or database lookup is needed.

## Why Left-Alignment Matters in Repeats (the silent miss)

The same variant can be written multiple ways:

```
chr1  100  ATCG  A      (right-aligned)
chr1  100  ATC   A      (left-aligned, parsimonious -- the canonical form)
chr1  101  TCG   T      (shifted position, different anchor base)
```

The ambiguity is worst in homopolymers and tandem repeats. In reference `...AAAAAAA...`, a single-base `A` deletion is positionally ambiguous: deleting ANY one of the seven A's yields the identical alternate haplotype, so different callers/aligners emit it at different POS. Left-alignment defines the canonical (leftmost) position, giving the one biological event one representation.

The load-bearing consequence: a one-base-off (non-left-aligned) indel in a homopolymer is a valid VCF line that produces a **different** `(CHROM, POS, REF, ALT)` tuple. Annotation and matching keyed on that tuple then silently **miss** the dbSNP/ClinVar/gnomAD entry that sits at the canonical coordinate -- a clinically actionable variant is reported as absent or novel, with **no error thrown**. The VCF is structurally valid; the numbers are wrong. `bcftools merge` of a normalized and an un-normalized cohort likewise emits duplicate rows for the same event, inflating counts and splitting allele frequencies.

Decision: always left-align + parsimony (`bcftools norm -f ref.fa`) against the EXACT downstream reference before annotation, database matching, set operations, or merging -- never rely on callers to emit canonical form.

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

## Tool Discordance: bcftools vs vt vs GATK

All three implement the same Tan-2015 left-align + parsimony core and agree on simple biallelic indels. They diverge on decomposition:

| Behavior | `vt` | `bcftools norm` | GATK `LeftAlignAndTrimVariants` |
|----------|------|-----------------|--------------------------------|
| Left-align + parsimony (biallelic) | yes (`normalize`) | yes (`-f`) | yes |
| Split multiallelic | `vt decompose` (separate step) | `-m-any` | `--split-multi-allelics` |
| Decompose MNP -> SNPs | `vt decompose_blocksub` splits MNPs **by default** (`vt decompose` only splits multiallelics) | only with `--atomize` (bcftools >= 1.12); **NOT by default** | does not decompose MNPs into SNPs |
| Block substitution / complex | `vt decompose_blocksub` | `--atomize` | limited |

The discordance that burns people: `vt decompose_blocksub` splits MNPs/block substitutions into SNPs by default, `bcftools norm` does not (it needs `--atomize`). Running the same normalization with the two tools yields a **different variant count** (bcftools keeps an MNP as one record; vt emits separate SNPs). If cohort A is atomized with vt and cohort B is left un-atomized with bcftools, every MNP is a systematic representation mismatch, manufacturing **spurious cohort-private "variants"** in any cross-cohort comparison.

Decision: standardize on ONE normalization tool + exact flag set across every cohort intended for comparison, and record the command. `bcftools norm` is the de-facto production standard (htslib-maintained; integrates split, atomize, and left-align in one pass). Note GATK's left-alignment window (`--max-leading-bases`) means indels inside repeats longer than the window may not fully left-align -- verify the installed default with `gatk LeftAlignAndTrimVariants --help` before trusting STR-embedded indels.

## Left-Alignment

**"Normalize my VCF before comparing callers"** -> Left-align indel representations and split multiallelic sites for consistent variant comparison.

```bash
bcftools norm -f reference.fa input.vcf.gz -Oz -o normalized.vcf.gz
```

Left-alignment cannot roll an indel leftward without the reference bases to its left, so `-f/--fasta-ref` is non-negotiable. Two traps make the reference choice load-bearing:
- The FASTA must be the **exact** reference the VCF was called against. A different build patch, `chr1`-vs-`1` contig naming, or a masked-vs-unmasked sequence yields **wrong** left-aligned coordinates. On a REF mismatch `bcftools norm` **errors and exits non-zero by default** (`-c e`) rather than warning -- do not paper over it with `-c w`; a single off-by-one in a contig shifts every indel.
- It must be the SAME reference used **downstream** (the annotation-database build, the other cohort). Normalizing to GRCh38 and matching against a GRCh37 dbSNP is a guaranteed miss even when local left-alignment is perfect.

### Check for Normalization Issues

```bash
bcftools norm -f reference.fa -c w input.vcf.gz > /dev/null
```

Check modes (`-c`):
- `e` - Error and exit on mismatch (default)
- `w` - Warn on mismatch and continue (use this to enumerate all mismatches)
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

### Field Reapportionment and Spanning Deletions

Splitting a multiallelic is NOT information-lossless. Per-allele fields must be reapportioned according to their header `Number` code, and `bcftools norm -m-` uses those codes to subset correctly:
- `Number=A` (one value per ALT, e.g. `AF`, `AC`) -> take the k-th element for the k-th ALT.
- `Number=R` (one per allele including REF, e.g. `AD` -> ref depth first) -> off-by-one relative to `A`.
- `Number=G` (one per genotype, e.g. `PL`, `GL`) -> subsetting needs the genotype-index formula.
- `Number=.` -> variable count; parsers **cannot** auto-subset, so these fields are silently carried whole onto every split record.

The trap: a custom INFO/FORMAT field mis-declared `Number=.` when it is really `A` is NOT subset on split, so every biallelic record keeps the full multiallelic vector and downstream tools read the **wrong allele's** value. Joining back (`-m+`) cannot always reconstruct the original PLs exactly. Rule: split once, early, and stay biallelic; join only for final delivery if a consumer requires it.

Spanning-deletion `*` alleles (VCF: "allele missing due to overlapping deletion") are meaningful only relative to the overlapping deletion recorded on another line. Splitting can strand a `*` allele from the deletion it references; bcftools handles this, but naive third-party splitters corrupt it. Use `--keep-sum AD` when the summed allele depth must be preserved across split records.

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

**Caveat -- decomposition destroys phase needed for functional annotation.** The original MNP record guarantees that its substitutions occur on the SAME haplotype. Atomization discards that guarantee. The concrete failure: two adjacent SNVs falling in one codon, annotated independently after decomposition, can each look **synonymous** while the true MNV (the haplotype) is **missense or nonsense** (or the reverse). VEP/SnpEff give the wrong amino-acid consequence on decomposed alleles because they no longer see the codon change.

This is the unresolved decompose-vs-atomic tension: decompose for variant **matching** (dbSNP/ClinVar/gnomAD lookup, allele-frequency comparison), but compute **functional consequence** on the haplotype-resolved (undecomposed / phased) representation -- run `bcftools csq` on the un-atomized VCF, which is codon-aware. Keep the atomized copy for matching and the un-atomized copy for annotation; do not feed atomized alleles to a per-record consequence caller. See variant-calling/variant-annotation.

### Atomize with Old Record Tag

```bash
bcftools norm --atomize --old-rec-tag ORIGINAL input.vcf.gz -Oz -o atomized.vcf.gz
```

Preserves the original record as an INFO annotation, enabling traceability back to the pre-atomized variant.

## The Normalization <-> HGVS 3'-Rule Clash

VCF normalization and HGVS nomenclature shift indels in **opposite** directions, so a correctly normalized VCF and a correct HGVS string for the same indel can name different repeat units:
- **VCF left-alignment** shifts an ambiguous indel to the most **5'** position on the **forward genomic strand**.
- **HGVS mandates the 3'-rule:** the indel is described at the most **3'** position with respect to the **transcript** (coding reading direction).

For a **plus-strand** gene, transcript-3' equals genomic-rightward -- the OPPOSITE end from VCF left-alignment. For a **minus-strand** gene, transcript-3' points genomic-leftward and can coincide with left-alignment, but only by accident of strand. The result: the VCF POS and the HGVS `c.` position for one deletion legitimately disagree, and a tool that generates HGVS by naively translating the left-aligned POS without re-shifting 3' on the transcript emits a non-compliant string. This is a leading cause of "two labs reported different c. positions for the same deletion."

Decision: **never hand-derive HGVS from POS.** Left-align + normalize the VCF for matching/merging (idempotent, reproducible), and rely on the annotation engine's HGVS generator to apply the transcript 3'-shift (VEP does this; verify SnpEff/ANNOVAR per version). When matching a patient variant to a ClinVar or literature HGVS string, normalize BOTH to the same representation first -- two strings that look different can be the same variant. See variant-calling/variant-annotation and variant-calling/clinical-interpretation.

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
| `--atomize` unrecognized | bcftools < 1.12 | Upgrade bcftools, or use `vt decompose_blocksub` as alternative |
| Split records carry wrong per-allele AF/AD | custom field mis-declared `Number=.` | Fix the header `Number` to `A`/`R`/`G` so bcftools subsets it on split |
| HGVS `c.` position disagrees with VCF POS | left-align (5' genomic) vs HGVS 3'-rule (transcript) | Expected; let the annotation engine emit HGVS, never hand-derive from POS |

## Related Skills

- variant-calling/variant-calling - Generate VCF files from alignments
- variant-calling/filtering-best-practices - Filter after normalization
- variant-calling/vcf-manipulation - Merge, intersect, and compare VCFs
- variant-calling/variant-annotation - Annotate normalized variants against databases
- variant-calling/gatk-variant-calling - GATK HaplotypeCaller workflow (does not emit MNPs)
- variant-calling/clinical-interpretation - ClinVar lookup requires normalized representation
- alignment-files/sam-bam-basics - BAM format and reference genome handling

## References

- Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants. *Bioinformatics.* 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112 (formal normalization definition: parsimony + left-alignment, and the `vt normalize` algorithm)
