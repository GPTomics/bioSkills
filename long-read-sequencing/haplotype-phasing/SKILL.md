---
name: bio-long-read-sequencing-haplotype-phasing
description: Phases small variants, SVs, and methylation from Oxford Nanopore and PacBio long reads (read-backed/physical phasing) with WhatsHap, LongPhase, or HiPhase, and haplotags the BAM (HP/PS tags) for allele-resolved downstream analysis. Covers why phase blocks break at het-sparse gaps (read length x heterozygosity), why phasing the VCF is useless until the BAM is haplotagged, the GT-pipe/PS and read HP/PS tag spec, reporting block N50 with switch error, the diploid-assumption/CNV/haploid-region traps, trio phasing as the gold standard, and the boundary to statistical panel phasing. Use when phasing long-read variants, haplotagging reads for allele-specific methylation/expression or phased SVs, choosing WhatsHap vs LongPhase vs HiPhase, trio phasing, or assessing phasing quality.
tool_type: cli
primary_tool: whatshap
---

## Version Compatibility

Reference examples tested with: whatshap 2.3+, longphase 1.7+, samtools 1.19+, tabix/htslib 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Behavior to record:
- `whatshap phase --reference` enables realignment mode (rescues indel phasing on error-prone long reads); omitting it falls back to lower-quality genotype-only phasing.
- Phasing the VCF and haplotagging the BAM are SEPARATE steps; downstream read-level tools need the BAM HP tag.
- `--max-coverage 15` (WhatsHap) is a runtime downsampling cap, not a minimum-depth requirement.

If code throws an error, introspect the installed tool (`whatshap phase --help`, `longphase --help`) and adapt the example to the actual API rather than retrying.

# Read-Backed Haplotype Phasing

**"Phase my long-read variants"** -> Reconstruct haplotypes directly from reads that span heterozygous sites, then haplotag the BAM so downstream tools can see the phase.
- CLI: `whatshap phase -o phased.vcf.gz --reference ref.fa --indels variants.vcf.gz aln.bam` then `whatshap haplotag -o haplotagged.bam --reference ref.fa phased.vcf.gz aln.bam`

This is read-backed (physical, panel-free) phasing of a single sample. Statistical/reference-panel phasing for imputation lives in phasing-imputation/haplotype-phasing; building phased haplotype contigs lives in genome-assembly/hifi-assembly.

## The Single Most Important Modern Insight -- Phasing the VCF Is Useless Until the BAM Is Haplotagged, and Blocks Break Where No Read Spans Two Hets

Read-backed phasing is sample-intrinsic and panel-free (the haplotypes are exactly what this individual's reads physically witness), with two load-bearing consequences:

1. **`phase` writes the VCF; `haplotag` writes the BAM - they are different products.** `whatshap phase` / `longphase phase` set GT pipe (`0|1`) and a PS phase-set in the VCF; they do NOT touch the BAM. Every read-level downstream tool (modkit `--partition-tag HP` for allele-specific methylation, pb-CpG-tools `--hap-tag HP`, Severus for phased SVs, IGV color-by-HP, `whatshap split`) keys on the per-read `HP` tag that ONLY `haplotag` writes. A user who runs `phase` and stops has a phased VCF and an un-haplotagged BAM, and the downstream step silently produces only an ungrouped partition. Triage: `samtools view haplotagged.bam | grep -m1 'HP:i:'`.
2. **Phase blocks break wherever no single read spans two adjacent hets.** Block length is capped by read length x heterozygosity: a long homozygous run (or a coverage/mapping dropout) ends a block no matter how long the reads are. A genome is phased into MANY blocks, not one haplotype per chromosome, and between-block phase is arbitrary until a long-range method (trio, Hi-C) stitches them. So "the genome is phased" is meaningless without block N50 AND switch error together.

## Tool Decision Tree

Switch-error accuracy is comparable across read-based tools (~0.1-0.4% on long reads); choose on speed, SV/mod co-phasing, platform, and pedigree.

| Scenario | Tool | Why |
|----------|------|-----|
| Careful default; indel phasing | WhatsHap (`--reference --indels`) | realignment mode rescues indels; widest downstream familiarity |
| Parents/pedigree sequenced | WhatsHap `--ped` (PedMEC) | the gold standard - chromosome-scale, lowest switch error |
| Whole-genome ONT speed | LongPhase (`--ont`) | ~10x faster; 30x human in ~1 min |
| Co-phase SVs / methylation into long blocks | LongPhase (`--sv-file`/`--mod-file`) | a phased SV bridges het-sparse gaps; block N50 ~25 Mbp |
| PacBio HiFi, joint small+SV+STR | HiPhase | PacBio-native one-pass phasing |
| Multi-tech (Hi-C / 10x) | HapCUT2 | models Hi-C/linked-read error |
| Inside PEPPER-Margin-DeepVariant | margin | legacy embedded haplotagger |

Clair3 uses WhatsHap (or LongPhase) internally to phase its het SNPs and haplotag the BAM feeding its full-alignment model - this skill owns that phase->haplotag mechanism (see clair3-variants).

## The Tags (the central distinction)

| Layer | Tag | Meaning |
|-------|-----|---------|
| VCF (per variant) | `GT` with `|` vs `/` | `0|1` phased (order = which haplotype carries ALT); `0/1` unphased |
| VCF (per variant) | `FORMAT/PS` (Integer) | phase-set / block id; variants sharing a PS are phased relative to each other (conventionally the first variant's position) |
| BAM (per read) | `HP:i:1` / `HP:i:2` | the haplotype this read was assigned to (written by `haplotag`) |
| BAM (per read) | `PS:i:<int>` | the phase set the read's assignment belongs to (matches the VCF PS) |

Unassigned reads carry NO HP tag (not `HP:i:0`). Do not confuse the VCF `HP` FORMAT tag (GATK style) with the BAM `HP` read tag.

## Phasing Quality - Report Block N50 AND Switch Error

| Metric | Tool | Trap |
|--------|------|------|
| phase-block N50/NG50 | `whatshap stats` | contiguity, not correctness; gameable by over-joining blocks (which raises switch errors) |
| phased fraction | `whatshap stats` | a tool can phase fewer easy sites to look better |
| switch error rate | `whatshap compare` | the primary accuracy number |
| switch vs flip decomposition | `whatshap compare` | a long switch propagates (damaging); a flip/short switch self-corrects (one wrong variant) - quote the decomposition |
| Hamming distance | `whatshap compare` | hypersensitive to switch position (a switch near a block start flips half the block) |

Long blocks with a high switch rate are worse, not better, than honest short blocks. Benchmark against a trio-/strand-seq-phased GIAB truth.

## Core Commands

```bash
# WhatsHap: phase (VCF), then haplotag (BAM). --reference enables realignment for indels.
whatshap phase -o phased.vcf.gz --reference ref.fa --indels variants.vcf.gz aln.bam
tabix -p vcf phased.vcf.gz
whatshap haplotag -o haplotagged.bam --reference ref.fa \
    --output-haplotag-list htlist.tsv.gz phased.vcf.gz aln.bam
samtools index haplotagged.bam

# Quality
whatshap stats --gtf blocks.gtf phased.vcf.gz                       # block N50, count, fraction
whatshap compare --names truth,mine truth.vcf.gz phased.vcf.gz      # switch error, flip decomposition

# Trio (gold standard) - --ped takes a PED file, not mother/father/child args
whatshap phase -o trio.vcf.gz --reference ref.fa --ped family.ped joint.vcf.gz mother.bam father.bam child.bam

# LongPhase: faster whole-genome, co-phase SNP+indel+SV(+5mC) into long blocks
longphase phase -s snps.vcf --indels --sv-file svs.vcf -b aln.bam -r ref.fa -o phased -t 16 --ont
longphase haplotag -s phased.vcf --sv-file phased_SV.vcf -b aln.bam -r ref.fa -o haplotagged -t 16

# Downstream consumer example: allele-specific methylation
modkit pileup haplotagged.bam asm/ --ref ref.fa --cpg --combine-strands --partition-tag HP
```

## Per-Method Failure Modes

### Phased VCF but no HP tags downstream
**Trigger:** running `phase` and pointing a read-level tool at the original BAM. **Mechanism:** `phase` writes the VCF only; the BAM HP tag comes from `haplotag`. **Symptom:** modkit returns only an ungrouped partition; IGV shows one color; Severus reports no phased SVs - all with no error. **Fix:** run `haplotag`; verify `samtools view ... | grep HP:i:`.

### Short blocks blamed on the tool
**Trigger:** a homozygosity-rich or inbred sample phasing into many short blocks. **Mechanism:** no intervening hets to link across a long homozygous run - intrinsic, not tool failure. **Symptom:** low block N50 despite good reads. **Fix:** expect it; use ultra-long reads or co-phase SVs (LongPhase) to bridge sparse-het gaps; only trio/Hi-C makes it chromosome-scale.

### Indels phased poorly
**Trigger:** `whatshap phase` without `--reference`. **Mechanism:** without realignment, allele support for indels in error-prone reads is noisy. **Symptom:** low indel phasing / errors. **Fix:** always pass `--reference ref.fa` (and `--indels`) on long reads.

### Confident phasing of a haploid/CNV region
**Trigger:** phasing chrX/Y/MT in an XY sample, or inside a CNV/segdup. **Mechanism:** the two-haplotype model is false there (hemizygous, >2 or 1 haplotype, or collapsed paralogs). **Symptom:** spurious micro-blocks, HP counts far from 50/50. **Fix:** treat phasing there as unreliable; do not interpret it as biology.

### Quoting N50 alone
**Trigger:** comparing phasers on block N50. **Mechanism:** N50 is inflated by over-joining, which raises switch errors. **Symptom:** "longer blocks" that are actually worse. **Fix:** report block N50 AND switch error together; use the flip decomposition.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Total depth ~15-20x for confident phasing | phasing practice | per-haplotype depth is ~half; below ~10x blocks fragment |
| `--max-coverage 15` is a runtime cap | WhatsHap | wMEC is exponential in per-site coverage; >15x is redundant, not required |
| long-read switch error ~0.1-0.4% | benchmarks vs trio truth | the achievable accuracy band |
| LongPhase SNP+SV block N50 ~25 Mbp | Lin 2022 | co-phasing SVs bridges het-sparse gaps (vs ~10-15 Mbp SNP-only) |
| ASM wants ~20x total | methylation practice | each haplotype must clear the ~10x per-site floor |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `modkit --partition-tag HP` has only an ungrouped partition | BAM never haplotagged | run `whatshap haplotag` / `longphase haplotag` |
| `--trio` flag not recognized | the flag is `--ped` | pass a PED file: `--ped family.ped` |
| Poor indel phasing | `--reference` omitted | add `--reference ref.fa --indels` |
| 0 reads usable in phase | BAM @RG sample != VCF sample | `--ignore-read-groups` (or fix sample names) |
| `longphase --platform ont` errors | platform is a bare flag | use `--ont` or `--pb` |
| Spurious phasing on chrX/CNV | diploid assumption violated | treat as unreliable; exclude haploid/CNV regions |

## References

- Patterson M, Marschall T, Pisanti N, et al. 2015. WhatsHap: weighted haplotype assembly for future-generation sequencing reads. *J Comput Biol* 22(6):498-509.
- Martin M, Patterson M, Garg S, et al. 2016. WhatsHap: fast and accurate read-based phasing. *bioRxiv* 085050.
- Garg S, Martin M, Marschall T. 2016. Read-based phasing of related individuals (PedMEC). *Bioinformatics* 32(12):i234-i242.
- Lin JH, Chen LC, Yu SC, Huang YT. 2022. LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants. *Bioinformatics* 38(7):1816-1822.
- Holt JM, Saunders CT, Rowell WJ, et al. 2024. HiPhase: jointly phasing small, structural, and tandem repeat variants from HiFi sequencing. *Bioinformatics* 40(2):btae042.
- Edge P, Bafna V, Bansal V. 2017. HapCUT2: robust and accurate haplotype assembly for diverse sequencing technologies. *Genome Res* 27(5):801-812.

## Related Skills

- clair3-variants - Produces the het VCF; Clair3 phases+haplotags internally via this mechanism
- long-read-alignment - Produces the BAM (keep `-Y` so supplementaries are taggable)
- nanopore-methylation - Allele-specific methylation via `modkit --partition-tag HP`
- structural-variants - Severus consumes a haplotagged BAM for phased/somatic SVs
- basecalling - LongPhase can co-phase 5mC from a modBAM
- phasing-imputation/haplotype-phasing - Statistical/reference-panel phasing for imputation
- genome-assembly/hifi-assembly - Phased de novo haplotype contigs (trio/Hi-C)
- hi-c-analysis/contact-pairs - Hi-C long-range phasing (orthogonal)
