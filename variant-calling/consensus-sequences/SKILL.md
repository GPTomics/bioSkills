---
name: bio-consensus-sequences
description: Generate consensus FASTA sequences by applying VCF variants onto a reference with bcftools consensus, or build viral/amplicon consensus with iVar. Use when reconstructing a sample-specific reference or haplotype, deciding -H haplotype vs IUPAC vs all-ALT projection, masking no-coverage sites so a consensus does not manufacture false reference calls, or setting iVar min-depth/min-frequency policy for surveillance genomes.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, samtools 1.19+, bedtools 2.31+, iVar 1.4+, minimap2 2.26+, BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: the `-H` argument vocabulary (N/R/A/I/LR/LA/SR/SA/NpIu, where `I` = IUPAC code for all genotypes) has grown across bcftools 1.x releases; IUPAC output is available both as the `-H I` code and the standalone `-I`/`--iupac-codes` flag. Always confirm the accepted letters with `bcftools consensus` on the installed version before scripting a selection.

# Consensus Sequences

**"Generate a consensus sequence from my VCF"** -> Apply called variants onto a reference FASTA, producing a sample-specific sequence, with a deliberate choice of haplotype projection and no-coverage masking.
- CLI (from a VCF): `bcftools consensus -f reference.fa input.vcf.gz`
- CLI (viral/amplicon from a BAM): `samtools mpileup ... | ivar consensus -p out`
- Python: `cyvcf2` + `Bio.SeqIO` for SNP-only prototypes

## The governing principle

`bcftools consensus` walks the reference and substitutes ALT alleles at the positions present in the VCF. Everything else is copied from the reference verbatim -- which drives three traps that ruin more consensus analyses than any tool bug:

1. **A consensus silently emits REFERENCE wherever the VCF is silent -- including positions with zero coverage.** No data and confidently-reference look identical in the output. An unmasked consensus therefore manufactures false confidence at exactly the sites where the sample was never observed. Mask no-coverage sites (below) or the FASTA lies.
2. **`-H 1` on an UNPHASED VCF yields a chimeric pseudo-haplotype.** Haplotype selection is only meaningful when genotypes are phased; on unphased data it mixes alleles from different real chromosomes into a sequence that exists in no cell. Verify `|` phasing before selecting a haplotype.
3. **A single FASTA cannot faithfully represent a diploid genome.** Every projection (`-H 1`, `-I`, `-H A`) is lossy in a different way; for phase-sensitive work keep the VCF, not the consensus.

The input VCF must be **bgzipped and indexed** (`bgzip` + `bcftools index`/`tabix`); plain-gzip or unindexed input errors out. The REF bases in the VCF must match the FASTA exactly or bcftools warns and skips those records. Normalize first (see Normalization).

## Basic Usage

`bcftools consensus` reads variants from a bgzipped, indexed VCF and writes FASTA:

```bash
bcftools index input.vcf.gz                                    # .csi index (or tabix -p vcf)
bcftools consensus -f reference.fa input.vcf.gz > consensus.fa
bcftools consensus -f reference.fa -o consensus.fa input.vcf.gz  # -o instead of redirect
```

For a multi-sample VCF, always pass `-s` -- without it, the applied genotypes are undefined:

```bash
bcftools query -l input.vcf.gz                                 # list samples
bcftools consensus -f reference.fa -s sample1 input.vcf.gz > sample1.fa
```

Restrict to a region with `-r` (the FASTA header is then `>chr:from-to`):

```bash
bcftools consensus -f reference.fa -r chr1:1000000-1010000 -s sample1 input.vcf.gz > gene.fa
```

## Haplotype Selection and the Phasing Trap

`-H` chooses which allele to apply from `FORMAT/GT`. The codes are case-insensitive:

| Option | Applies | Use when |
|--------|---------|----------|
| `-H 1` / `-H 2` | Allele at GT index 1 or 2 | Emitting one true chromosome -- **only valid on PHASED genotypes** |
| `-H A` | ALT allele in every genotype | Maximum divergence from reference; a chimera of both chromosomes |
| `-H R` | REF allele at heterozygous sites | Conservative consensus; discards het ALT alleles |
| `-H I` (or the standalone `-I` / `--iupac-codes` flag) | IUPAC ambiguity code | Retain heterozygosity in one sequence (see caveat below) |
| `-H LA`/`LR`/`SA`/`SR` | Longer/shorter allele, tie broken by ALT/REF | Length-driven selection; confirm the letter set on the installed version |

**The chimeric-haplotype footgun.** `-H 1`/`-H 2` are only meaningful when genotypes are phased (`0|1`, pipe separator). With **unphased** genotypes (`0/1`, slash), the assignment of "which allele is haplotype 1" is arbitrary *per site*, so `-H 1` across many heterozygous sites produces a **switch-error mosaic that corresponds to no real chromosome** -- while looking like a clean haplotype FASTA. This is the single most dangerous consensus mistake. Verify phasing before any `-H 1`/`-H 2`:

```bash
bcftools query -f '%CHROM\t%POS[\t%GT]\n' input.vcf.gz | head   # phased: 0|1 ; unphased: 0/1
```

If genotypes are unphased, phase first (read-backed WhatsHap/HapCUT2, trio, statistical SHAPEIT/Eagle -- accurate for common variants, poor for rare/singletons -- or native long-read phasing). See phasing-imputation/haplotype-phasing and variant-calling/vcf-basics for GT interpretation.

## What a Consensus Cannot Represent

A single consensus FASTA is a lossy projection of a diploid genome; the right projection depends on the downstream use, and some tasks need the VCF instead:

| Strategy | Flag | Best for | Loses |
|----------|------|----------|-------|
| Two haplotype sequences | `-H 1` + `-H 2` (phased) | Allele-specific expression, compound-het, HLA, cis-regulatory haplotypes | Nothing (if correctly phased) |
| IUPAC ambiguity codes | `-I` | Retaining het signal in one sequence | Phase/linkage; **many tree/alignment tools read IUPAC as N** |
| All ALT alleles | `-H A` | Max divergence, quick draft | Reality -- exists in no cell |
| REF at het sites | `-H R` | Conservative single sequence | Every heterozygous ALT allele |

Two hard boundaries:

- **For phase-sensitive work, keep the VCF (or two phased haplotype FASTAs), not a single consensus.** Collapsing hets to IUPAC or picking one allele discards linkage that the analysis needs -- treating a consensus FASTA as "the sample's genome" for compound-het or allele-specific analysis is a category error.
- **`bcftools consensus` cannot apply symbolic SV alleles** (`<DEL>`, `<INS>`, `<DUP>`, `<INV>`): those carry no ALT sequence, only INFO fields, so consensus has nothing to substitute. Short-read SV VCFs (Manta/DELLY) are mostly symbolic and are NOT directly consensus-able. Folding SVs into a consensus needs sequence-resolved records (long-read/assembly callers emit these) or an assembly-based approach -- see variant-calling/structural-variant-calling.

For phylogenetics specifically, prefer one clean phased haplotype or a homozygous-ALT-only sequence over IUPAC, because ambiguity codes are silently dropped by many tree builders:

```bash
bcftools view -i 'GT="AA"' input.vcf.gz | bcftools consensus -f reference.fa > hom_alt.fa
```

## Masking No-Coverage Sites (the load-bearing footgun)

Because unobserved positions are emitted as reference (trap 1), a consensus must mask sites with insufficient data. `-m mask.bed` replaces the listed regions (default char N via `--mask-with N`). The mask must be built from **callable depth**, and the depth step hides a silent bug:

**`samtools depth` WITHOUT `-a` OMITS zero-coverage positions** from its output -- so those positions never enter the low-depth BED, never get masked, and stay as reference: the exact false-confidence failure the mask was meant to prevent. Always use `-a` (report all positions) so no-coverage sites are captured:

```bash
# Build a mask of every position below the callable-depth threshold. -a is mandatory:
# without it, zero-coverage positions are absent from the output and escape masking.
samtools depth -a aligned.bam | awk '$3 < 10 {print $1"\t"$2-1"\t"$2}' | bedtools merge > lowcov.bed

bcftools consensus -f reference.fa -m lowcov.bed input.vcf.gz > consensus.fa
```

The `< 10` threshold is a minimum-callable-depth policy (10x is a common floor for confident base calls); set it to the depth below which the calls are not trusted. `bedtools genomecov -bga -ibam aligned.bam` is an equivalent zero-coverage-aware alternative that also emits 0-depth intervals.

Do NOT rely on `-M`/`-a` for this: `-M N` outputs N only for missing `./.` genotypes already present in the VCF, and `-a N` replaces every position absent from the VCF (which N-outs the entire non-variant genome). Neither distinguishes no-coverage from confident-reference -- only a depth-derived mask does.

## Normalization Before Consensus

**Goal:** Apply indels at the correct reference position and sequence.

**Approach:** Left-align and split multiallelics with `bcftools norm` so each record matches the reference context; consensus applies records positionally and mis-represented indels corrupt the output.

```bash
bcftools norm -f reference.fa input.vcf.gz -Oz -o norm.vcf.gz
bcftools index norm.vcf.gz
bcftools consensus -f reference.fa norm.vcf.gz > consensus.fa
```

Un-normalized or overlapping indels produce wrong sequence, and `bcftools consensus` only warns to stderr while still emitting output -- so the corruption is silent unless the stderr is inspected. Even after norm, two records whose REF spans collide remain a hazard; grep the run for warnings and inspect the region. See variant-calling/variant-normalization.

```bash
bcftools consensus -f reference.fa norm.vcf.gz 2>&1 >consensus.fa | grep -i 'overlap\|warn'
```

## Viral / Amplicon Consensus with iVar

For amplicon surveillance (SARS-CoV-2 and similar), `ivar consensus` builds a per-sample consensus directly from a pileup. Its two key thresholds are **epidemiological policy decisions, not defaults to accept blindly** -- they propagate into lineage assignment and transmission-cluster inference:

```bash
# Trim PCR primers FIRST -- primer-derived bases are not sample sequence and, at
# primer-binding-site mutations, cause reference-biased miscalls if left in.
ivar trim -b primers.bed -p trimmed -i aligned.bam
samtools sort -o trimmed.sorted.bam trimmed.bam

# -aa keeps all positions (so no-coverage becomes N), -A keeps orphan mates, -d 0 lifts the depth cap.
samtools mpileup -aa -A -d 0 -B -Q 0 trimmed.sorted.bam | ivar consensus -p sample -q 20 -t 0.5 -m 10 -n N
```

| Flag | Default | Decision |
|------|---------|----------|
| `-m` min depth | 10 | Below this, iVar emits N. Too low -> single-read sequencing errors become "mutations" that corrupt outbreak phylogenies. Too high -> excessive Ns, an unusably fragmented genome. |
| `-t` min frequency to call a base | 0 (majority) | 0 calls the most common base. For a strict majority consensus use 0.5. Too low bakes minority/within-host variants and contamination into the "genome", inflating diversity and creating phantom transmission links. Raise (e.g. 0.03) only deliberately for intrahost variant work, not for a reference consensus. |
| `-q` min base quality | 20 | Bases below this are not counted toward depth/frequency. |
| `-n` no-coverage char | N | Character emitted where depth `< -m`. |

Always report `-m` and `-t` alongside a surveillance consensus -- the genome is only as trustworthy as those two numbers. Alternatives: `bcftools consensus` from a called VCF, or ViralConsensus (Moshiri 2023) which calls consensus directly from the alignment without an intermediate VCF, faster and lower-memory for large batches.

## Filtering Before Consensus

Apply only trusted calls; pipe filtered VCF straight into consensus:

```bash
bcftools view -f PASS input.vcf.gz -Oz -o pass.vcf.gz && bcftools index pass.vcf.gz
bcftools consensus -f reference.fa pass.vcf.gz > consensus.fa

bcftools view -v snps input.vcf.gz -Oz -o snps.vcf.gz && bcftools index snps.vcf.gz  # SNPs only
```

Filtered VCFs must be re-bgzipped and re-indexed before `bcftools consensus` reads them.

## Chain Files and Naming

`-c chain.txt` writes a liftover chain mapping reference coordinates to consensus coordinates -- needed when indels shift positions and annotations must be lifted. `-p PREFIX` prepends a string to output sequence names (`>sample1_chr1`).

```bash
bcftools consensus -f reference.fa -c chain.txt -p "sample1_" input.vcf.gz > consensus.fa
```

## cyvcf2 Consensus (SNP-only prototypes)

For a quick SNP-only substitution in Python (production work should use `bcftools consensus`, which handles indels, phasing, and masking):

```python
from cyvcf2 import VCF
from Bio import SeqIO

ref = {rec.id: list(str(rec.seq)) for rec in SeqIO.parse('reference.fa', 'fasta')}
for v in VCF('input.vcf.gz'):
    if v.is_snp and len(v.ALT) == 1:
        ref[v.CHROM][v.POS - 1] = v.ALT[0]   # POS is 1-based; list index is 0-based
with open('consensus.fa', 'w') as fh:
    for chrom, seq in ref.items():
        fh.write(f'>{chrom}\n{"".join(seq)}\n')
```

## Verify the Consensus

```bash
minimap2 -a reference.fa consensus.fa | samtools view -b -o aln.bam   # inspect where it diverges
bcftools view -H input.vcf.gz | wc -l                                 # variants available to apply
```

## Common Errors

| Error / Symptom | Cause | Fix |
|-----------------|-------|-----|
| `the VCF file is not indexed` | Plain-gzip or missing index | `bgzip` then `bcftools index` (or `tabix -p vcf`) |
| `sequence "chr1" not found` | Chromosome names differ between FASTA and VCF | `bcftools annotate --rename-chrs map.txt` |
| `REF does not match` | Different reference than the caller used | Use the exact FASTA used for calling; normalize |
| Clean haplotype looks wrong | `-H 1` on an unphased VCF -> chimera | Verify `|` phasing; phase before `-H` |
| Consensus reference-identical over gaps | No-coverage sites emitted as reference | Mask with `samtools depth -a` derived BED and `-m` |
| Garbled indels, stderr overlap warnings | Un-normalized/overlapping records | `bcftools norm -f ref.fa` first; inspect warnings |
| `<DEL>`/`<INS>` not applied | Symbolic SV alleles carry no ALT sequence | Use sequence-resolved SV records; see structural-variant-calling |

## Related Skills

- variant-calling/variant-calling - Generate the VCF consensus is built from
- variant-calling/vcf-basics - Interpret GT and phasing (`|` vs `/`) before `-H`
- variant-calling/variant-normalization - Left-align indels before consensus
- variant-calling/filtering-best-practices - Restrict to trusted calls first
- variant-calling/structural-variant-calling - Sequence-resolved SVs for SV-aware consensus
- phasing-imputation/haplotype-phasing - Produce phased genotypes for true haplotypes
- phylogenetics/modern-tree-inference - Build trees from a consensus alignment

## References

- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008. (bcftools consensus / norm / mpileup.)
- Grubaugh ND, Gangavarapu K, Quick J, Matteson NL, De Jesus JG, Main BJ, et al. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. *Genome Biology.* 2019;20(1):8. doi:10.1186/s13059-018-1618-7. (iVar consensus/trim; depth `-m` and frequency `-t` thresholds.)
- Moshiri N. ViralConsensus: a fast and memory-efficient tool for calling viral consensus genome sequences directly from read alignment data. *Bioinformatics.* 2023;39(5):btad317. doi:10.1093/bioinformatics/btad317.
