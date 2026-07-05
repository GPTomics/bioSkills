---
name: bio-vcf-basics
description: View, query, and interpret VCF/BCF variant files with bcftools and cyvcf2. Use when inspecting variants, extracting fields with query format strings, converting VCF/BCF, or correctly reading a field -- QUAL (site) vs GQ (genotype) vs PL/GL likelihoods, AD vs DP and allele balance, GT phasing/ploidy/PS and missing-vs-hom-ref, INFO/FORMAT Number A/R/G semantics, symbolic alleles (<DEL>, <NON_REF>, spanning *) and END, or telling a raw gVCF apart from a filtered callset.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, cyvcf2 0.30+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# VCF/BCF Basics

**"Show me and extract fields from this VCF"** -> Parse the VCF/BCF format, then view, subset, or pull specific columns into a flat table.
- CLI: `bcftools view` / `bcftools query -f`
- Python: `cyvcf2.VCF` (iterate records with attribute access)

## The governing principle

A VCF field is only meaningful once its LEVEL and its `Number` are known. QUAL is a site-level property; GQ, PL, AD, DP, GT are per-sample. QUAL and GQ answer different questions and are NOT interchangeable. A field's header `Number` (A/R/G/.) dictates how many values it carries and how it must be re-subset after a multiallelic split. And several encodings are load-bearing traps: `.` (missing) is never `0/0` (hom-ref); a bare `*` ALT is a spanning-deletion placeholder, not an allele; a gVCF `<NON_REF>` record is a reference-confidence intermediate, not a filtered call. Read the header, read the Number, read the level -- a structurally valid VCF read at the wrong level silently produces wrong numbers with no error.

## Format Overview

| Format | Description | Use Case |
|--------|-------------|----------|
| VCF | Text format, human-readable | Debugging, small files |
| VCF.gz | Compressed VCF (bgzip) | Standard distribution |
| BCF | Binary VCF | Fast processing, large files |

## VCF Format Structure

```
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
chr1    1000    rs123   A       G       30      PASS    DP=50   GT:DP   0/1:25
```

### Header Lines (##)
- `##fileformat` - VCF version
- `##INFO` / `##FORMAT` - INFO / FORMAT field definitions (ID, Number, Type)
- `##FILTER` - Filter definitions
- `##contig` - Reference contigs (required for region indexing and contig order)
- `##reference` - Reference genome

### The Header Contract

Every INFO/FORMAT tag used in the body MUST be declared in a `##INFO`/`##FORMAT` line giving its ID, Number, and Type; parsers (bcftools, cyvcf2, pysam) read these declarations to know how many values a field holds and how to type it. An out-of-sync header -- a tag used but not declared, or declared with the wrong Number/Type -- silently breaks parsing: a `Number=1` declaration over data that holds a vector, or a missing `##contig`, makes tools drop, mistype, or mis-subset values with NO error thrown. After any hand-edit or annotation that adds a field, update the header to match (`bcftools +fill-tags` and `bcftools annotate` manage this automatically).

### Data Columns

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| POS | 1-based position of the first base in REF (contrast BED's 0-based half-open) |
| ID | Variant identifier (e.g., rs number) or `.` if novel |
| REF | Reference allele (matches the reference exactly) |
| ALT | Alternate allele(s), comma-separated. `*` = allele missing due to an overlapping deletion at this site |
| QUAL | Phred-scaled quality of the ALT assertion, `-10*log10 P(no variant)`; higher = more confident a variant exists (site-level, NOT per-sample) |
| FILTER | PASS or semicolon-separated filter names. `.` means filters were not applied |
| INFO | Semicolon-separated key=value pairs (site-level annotations) |
| FORMAT | Colon-separated format keys defining per-sample field order |
| SAMPLE | Colon-separated values matching FORMAT order |

## Critical Field Interpretation

What each field actually measures -- and what it does not -- drives every filtering and interpretation decision. QUAL, GQ, and PL answer three DIFFERENT questions.

### QUAL vs GQ vs PL/GL: three different confidences

| Field | Level | Scale | Question answered |
|-------|-------|-------|-------------------|
| QUAL (col 6) | Site | Phred: `-10*log10 P(no variant)` | "Is there ANY variant at this site?" |
| GQ (FORMAT) | Genotype | Phred, capped at 99 | "Is THIS sample's assigned genotype correct?" |
| PL (FORMAT) | Genotype | Phred, rebased to min=0 | Relative likelihood of every possible genotype |
| GL (FORMAT) | Genotype | log10, `<=0`, raw | Same info as PL, unscaled (`PL = -10*GL`, rebased) |

QUAL is computed once across all samples and SCALES with total depth, so a high-coverage artifact can carry a large QUAL -- hence QD (QUAL normalized by depth) is preferred for filtering. GQ is per-sample and does not scale with cohort size. They are NOT interchangeable: QUAL can be high while an individual genotype is uncertain (low GQ), and a sample can have a confident genotype (high GQ) at a site with only moderate QUAL. Filter site-level junk on QUAL/QD; no-call untrustworthy genotypes on GQ.

### PL/GL, and how GQ is derived

PL holds phred-scaled genotype likelihoods, rebased so the CALLED (most likely) genotype is exactly 0 and every other value is its phred penalty relative to that call. For a biallelic diploid site PL is ordered `[PL(0/0), PL(0/1), PL(1/1)]` -- the index of the 0 IS the genotype the caller assigned. GL is the same information as raw log10 likelihoods (`<=0`, larger is better). GQ = the difference between the two SMALLEST PL values, i.e. phred confidence in the call versus the next-best genotype; GQ=0 means the top two genotypes are tied (uninformative), GQ is capped at 99 by convention.

For a site with n alleles, diploid genotype `j/k` (j<=k) sits at PL index `k*(k+1)/2 + j` (this is the `Number=G` ordering). Getting this index formula wrong is the classic bug when re-parsing PL after a multiallelic split -- the vector must be re-subset by the formula, never sliced positionally.

### AD vs DP, and allele balance

AD (FORMAT, `Number=R`) is per-allele read depth `[ref_depth, alt1_depth, ...]`, REF first. DP is total depth. `sum(AD)` is often LESS than DP -- expected, not an error:
- DP counts all reads spanning the position, including uninformative reads (low base quality, ambiguous alignment, filtered reads).
- AD counts only reads that confidently support a specific allele.
- INFO/DP (site-level, summed across samples) differs from FORMAT/DP (per-sample).

Allele balance for a het is DERIVED from AD (GATK does not emit it directly): `AB = alt_AD / (ref_AD + alt_AD)`. A true het sits near 0.5; hets far from 0.5 (e.g. `<0.2` or `>0.8`) suggest a mapping artifact, CNV, or contamination.

### INFO/FORMAT Number semantics (A / R / G / .)

Every `##INFO`/`##FORMAT` header declares a `Number` telling a parser how many values a field holds AND how to re-subset it when a multiallelic record is split:

| Number | One value per | Examples | On multiallelic split |
|--------|---------------|----------|-----------------------|
| `A` | ALT allele | AF, AC | take the k-th element for the k-th ALT |
| `R` | allele incl. REF | AD | REF value first, then per-ALT (off-by-one vs A) |
| `G` | genotype | PL, GL | re-subset via the `k*(k+1)/2+j` index formula |
| `.` | variable/unknown | -- | parser CANNOT auto-subset; carried whole onto every record |
| `0` | flag (presence only) | -- | -- |

Load-bearing for correctness: `bcftools norm -m-` uses these codes to reapportion fields on split. A field mis-declared `Number=.` when it is really `A` keeps its full multiallelic vector on every split record, so downstream tools read the WRONG allele's value with no error. See variant-calling/variant-normalization for split/join reapportionment.

### Key INFO Annotations for Filtering

| Annotation | Meaning | What It Detects |
|-----------|---------|-----------------|
| QD | QUAL / allele depth | Low values suggest variant quality not supported by reads |
| FS | Fisher strand bias (phred-scaled) | Variant reads predominantly on one strand (artifact) |
| SOR | Strand odds ratio | Same as FS but handles high-depth sites better |
| MQ | Root mean square mapping quality | Low values indicate reads map ambiguously (paralogous regions) |
| MQRankSum | MQ difference: ref vs alt reads | Very negative = alt reads map much worse than ref (suspicious) |
| ReadPosRankSum | Read position: ref vs alt reads | Very negative = variant only at read ends (misalignment artifact) |

## Genotype Encoding

| Genotype | Meaning |
|----------|---------|
| `0/0` | Homozygous reference (confidently called ref) |
| `0/1` | Heterozygous |
| `1/1` | Homozygous alternate |
| `1/2` | Heterozygous for two different ALT alleles (compound het at multiallelic site) |
| `./.` | Missing genotype (no confident call) |
| `0\|1` | Phased heterozygous (allele before `\|` is on haplotype 1) |

### Phased vs Unphased

- `/` separates **unphased** alleles -- the two chromosomal copies are known, but which came from which parent is not
- `|` separates **phased** alleles -- haplotype assignment is known (read-backed phasing, trio analysis, or long-read sequencing)
- Phasing matters for compound heterozygosity: two variants in a gene are pathogenic together only if on different haplotypes (in *trans*), not the same haplotype (in *cis*)

### Phase Sets (PS)

A `|` is only meaningful WITHIN a phase set. The FORMAT/PS tag (an integer, usually the POS of the block's first variant) groups variants phased relative to EACH OTHER; `0|1` in two different PS blocks are not guaranteed to lie on the same physical haplotype. Read-backed phasers (WhatsHap) and trio phasing emit PS. A `|` with no consistent PS across records carries no global phase -- a subtle trap when merging phased VCFs.

### Ploidy and Missing vs Hom-Ref

Ploidy is read from the NUMBER of alleles in GT: `0/1` diploid, `0` haploid (chrY, chrM, male chrX outside the PAR), `0/1/1` triploid. Per-region ploidy (PAR, chrX in males, mito) must match the sample karyotype.

`.` (missing) is NOT reference. `./.` = no-call (genotype could not be determined, usually low depth); `0/0` = confidently called homozygous reference. Treating `./.` as `0/0` inflates the reference-allele count and biases allele frequencies, missingness, and burden tests. This is load-bearing: never impute `./.` as reference. In a gVCF, the ABSENCE of a record also does not mean reference -- see the reference-confidence model below.

### Multiallelic Genotypes

At multiallelic sites (e.g. ALT = G,T), allele indices reference the comma-separated ALT list: 0=REF, 1=first ALT, 2=second ALT. `1/2` means one copy of each ALT. Splitting multiallelics into biallelic records with `bcftools norm -m-` converts `1/2` into two `0/1` records, losing compound-heterozygosity information -- see variant-calling/variant-normalization for caveats.

## Symbolic Alleles, END, and Spanning Deletions

Not every ALT spells out a sequence. Symbolic alleles are angle-bracketed placeholders for events whose sequence is not given inline:

| ALT | Meaning |
|-----|---------|
| `<DEL>` `<DUP>` `<INS>` `<INV>` `<CNV>` | Structural-variant classes (sequence not spelled out) |
| `<NON_REF>` | gVCF: "any allele not yet observed" (reference-confidence model) |
| `<*>` | Same role as `<NON_REF>` in some callers' gVCF/mpileup output |
| `*` (bare) | Spanning deletion: allele MISSING because an upstream deletion on ANOTHER line overlaps this position |

Two parsing traps:
- `INFO/END` gives the end coordinate of a symbolic/large event. A tool that infers a record's span from `len(REF)` is WRONG for symbolic alleles -- it must read END. gVCF reference blocks also use END to mark the last position of the band.
- The bare `*` ALT is interpretable only relative to the overlapping deletion on another record; it is not a real alternate allele here. Splitting/subsetting can strand a `*` from the deletion it references (see variant-calling/variant-normalization).

### gVCF and the `<NON_REF>` Reference-Confidence Model

A gVCF (GATK HaplotypeCaller `-ERC GVCF`) is fundamentally different from a filtered callset: it emits a record for EVERY position or block, not just variant sites. Non-variant stretches are compressed into END-delimited blocks (bands) grouped by GQ, so a gVCF is not one line per base.

- Every record carries a symbolic `<NON_REF>` ALT with PL/AD computed against "any unseen allele." This lets joint genotyping evaluate a site in THIS sample even when the variant was only discovered in ANOTHER cohort sample -- the `<NON_REF>` likelihood supplies the evidence.
- Its purpose is to distinguish, at every site, confident homozygous reference from no-data/no-call -- solving the missing-vs-reference problem when squaring off a cohort matrix.
- A gVCF is NOT ready for analysis; it is an intermediate. It must be joint-genotyped (`GenomicsDBImport`/`CombineGVCFs` -> `GenotypeGVCFs`) to yield a normal VCF. Do NOT filter, annotate, or count variants on a raw gVCF, and never build a multi-sample callset by `bcftools merge`-ing single-sample project VCFs when gVCF joint-genotyping is available -- merging fabricates hom-ref genotypes. See variant-calling/joint-calling.

## bcftools view

**Goal:** View, subset, and convert VCF/BCF files from the command line.

**Approach:** Use `bcftools view` with flags for header control, region selection, sample extraction, and format conversion.

```bash
bcftools view input.vcf.gz | head           # full records
bcftools view -h input.vcf.gz               # header only
bcftools view -H input.vcf.gz | head        # skip header
bcftools view input.vcf.gz chr1:1000000-2000000   # region (needs index)
bcftools view -s sample1,sample2 input.vcf.gz      # keep samples
bcftools view -s ^sample3 input.vcf.gz             # exclude samples
```

## bcftools query

**Goal:** Extract specific fields from a VCF in a custom tabular format.

**Approach:** Use `bcftools query -f` with format specifiers for CHROM, POS, INFO, and FORMAT fields. Square brackets `[...]` loop over samples.

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf.gz
bcftools query -f '%CHROM\t%POS\t%INFO/DP\t%INFO/AF\n' input.vcf.gz
bcftools query -f '%CHROM\t%POS[\t%GT]\n' input.vcf.gz              # per-sample GT
bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT]\n' -s sample1 input.vcf.gz
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf.gz     # column header
```

### Common Format Specifiers

| Specifier | Description |
|-----------|-------------|
| `%CHROM` `%POS` `%ID` | Position fields |
| `%REF` `%ALT` | Alleles |
| `%QUAL` `%FILTER` | Site quality / filter status |
| `%INFO/TAG` | INFO field value |
| `%TYPE` | Variant type (snp, indel, etc.) |
| `[%GT]` `[%DP]` `[%AD]` `[%GQ]` | Per-sample FORMAT fields (loop in `[...]`) |
| `[%SAMPLE]` | Sample name |
| `\n` `\t` | Newline / tab |

## Format Conversion and Indexing

**Goal:** Convert between VCF, compressed VCF, and BCF, and index for region queries.

**Approach:** Use `bcftools view` output flags (`-Ov/-Oz/-Ou/-Ob`), then bgzip + index.

```bash
bcftools view -Ob -o output.bcf input.vcf.gz   # VCF -> BCF
bcftools view -Ov -o output.vcf input.bcf      # BCF -> VCF
bgzip input.vcf                                 # -> input.vcf.gz (bgzip, NOT gzip)
bcftools index input.vcf.gz                     # -> .csi index
bcftools index -t input.vcf.gz                  # -> .tbi (tabix) index
```

### Output Format Options

| Flag | Format |
|------|--------|
| `-Ov` | Uncompressed VCF |
| `-Oz` | Compressed VCF (bgzip) |
| `-Ou` | Uncompressed BCF (fast piping) |
| `-Ob` | Compressed BCF |

BCF is the binary encoding of VCF: faster to parse and smaller for large callsets. Region queries (`chr1:1-1000`) require a bgzipped+indexed VCF or a BCF -- plain `.gz` (gzip) is not seekable and fails.

## cyvcf2 Python Alternative

**Goal:** Read, query, and write VCF files programmatically in Python.

**Approach:** Use cyvcf2's `VCF` reader to iterate variants with attribute access to fields, and `Writer` to emit filtered output.

**"Parse this VCF in Python"** -> Open with cyvcf2 and iterate variant records.

### Open, Iterate, and Access Fields
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
for variant in vcf:
    # ALT is a list; QUAL is site-level and may be None
    print(variant.CHROM, variant.POS, variant.REF, variant.ALT)
    print(variant.ID, variant.QUAL, variant.FILTER, variant.var_type)
    dp = variant.INFO.get('DP')   # INFO field, None if absent
    af = variant.INFO.get('AF')
    break
vcf.close()
```

### Access Genotypes and Per-Sample Fields
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
samples = vcf.samples
for variant in vcf:
    # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN(missing), 3=HOM_ALT
    gts = variant.gt_types
    depths = variant.format('DP')   # numpy array, one row per sample
    gqs = variant.format('GQ')      # per-sample genotype quality
    ad = variant.format('AD')       # per-allele depth, Number=R
    print(dict(zip(samples, gts)))
    break
vcf.close()
```

Note: cyvcf2 codes missing genotypes as `gt_types == 2` (UNKNOWN) -- treat that as no-call, never as HOM_REF.

### Fetch Region and Read the Header
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
print(vcf.samples, vcf.seqnames)      # sample names, contig names
for info in vcf.header_iter():
    if info['HeaderType'] == 'INFO':
        print(info['ID'], info['Description'])
for variant in vcf('chr1:1000000-2000000'):   # requires an index
    print(variant.CHROM, variant.POS)
```

### Write Filtered VCF
```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('output.vcf', vcf)   # inherit the input header
for variant in vcf:
    if variant.QUAL is not None and variant.QUAL > 30:   # QUAL is site-level
        writer.write_record(variant)
writer.close()
vcf.close()
```

## Quick Reference

| Task | bcftools | cyvcf2 |
|------|----------|--------|
| View VCF | `bcftools view file.vcf.gz` | `VCF('file.vcf.gz')` |
| View header | `bcftools view -h file.vcf.gz` | `vcf.header_iter()` |
| Get region | `bcftools view file.vcf.gz chr1:1-1000` | `vcf('chr1:1-1000')` |
| Query fields | `bcftools query -f '%CHROM\t%POS\n'` | Loop with properties |
| Count variants | `bcftools view -H file.vcf.gz \| wc -l` | `sum(1 for _ in vcf)` |
| VCF to BCF | `bcftools view -Ob -o out.bcf in.vcf.gz` | Use Writer |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `no BGZF EOF marker` | Not bgzipped (plain gzip) | Recompress with `bgzip`, not `gzip` |
| `index required` / region query fails | Missing index | Run `bcftools index` (`-t` for tabix) |
| `sample not found` | Wrong sample name | Check with `bcftools query -l` |
| INFO/FORMAT field missing or mistyped | Header out of sync with body | Fix `##INFO`/`##FORMAT` Number/Type; use `bcftools +fill-tags` |
| Every hom-alt or missing site vanishes on filter | Treated `.`/`./.` as failing or as ref | Missing != hom-ref; make missing pass, never impute `0/0` |
| Wrong allele's AF/AD after split | `Number=.` field not re-subset | Declare the true `Number` (A/R/G) so bcftools reapportions |

## Related Skills

- variant-calling/variant-calling - Generate VCF from alignments
- variant-calling/variant-normalization - Split multiallelics, left-align, Number-code reapportionment
- variant-calling/filtering-best-practices - Filter variants by site (QUAL/QD) and genotype (GQ/DP)
- variant-calling/joint-calling - gVCF reference-confidence model and joint genotyping
- variant-calling/vcf-manipulation - Merge, concat, intersect VCFs
- alignment-files/pileup-generation - Generate pileup for calling

## References

- Danecek P, Auton A, Abecasis G, et al. The variant call format and VCFtools. *Bioinformatics.* 2011;27(15):2156-2158. doi:10.1093/bioinformatics/btr330 (VCF format definition)
- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008 (bcftools view/query/norm reference)
- The Variant Call Format Specification (VCFv4.3/4.4). GA4GH / samtools hts-specs. https://samtools.github.io/hts-specs/ (symbolic alleles, END, `*` overlapping-deletion allele, Number=A/R/G, PL/GL/GQ, gVCF `<NON_REF>`)
