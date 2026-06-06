---
name: bio-long-read-sequencing-medaka-polishing
description: Polishes Oxford Nanopore draft assemblies to higher consensus accuracy with medaka, a basecaller-model-specific neural consensus net, produces haploid variant calls (VCF) for microbial, mitochondrial, or viral samples, and generates amplicon/viral consensus sequences. Covers the model-matching footgun that silently degrades output, why Racon-first is obsolete and medaka runs directly on Flye output as a single pass, why HiFi must never be fed to medaka, the v1->v2 subcommand renames, and the precise medaka_variant deprecation. Use when polishing an ONT-only assembly, generating an amplicon/viral consensus, calling a haploid ONT consensus, or deciding whether medaka, dorado polish, or Clair3 is the right tool.
tool_type: cli
primary_tool: medaka
---

## Version Compatibility

Reference examples tested with: medaka 2.2+, minimap2 2.28+, samtools 1.19+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Results depend on inputs that outlive the binary version - record them:
- The medaka MODEL must match the basecaller (pore + chemistry + speed + mode + version), e.g. `r1041_e82_400bps_sup_v5.2.0`. A mismatch silently degrades output. Prefer auto-detection from the BAM; verify with `medaka tools list_models`.
- Default models advance with each release (consensus `..._sup_v5.2.0`, variant `..._sup_variant_v5.0.0` at time of writing); confirm with `medaka tools list_models`.
- medaka v2 renamed subcommands (`consensus`->`inference`, `stitch`->`sequence`, `variant`->`vcf`) and moved the backend to PyTorch; v1 tutorials fail.

If code throws an error, introspect the installed tool (`medaka --help`, `medaka_consensus --help`) and adapt the example to the actual API rather than retrying.

# Medaka Polishing

**"Polish my Nanopore assembly"** -> Run one medaka consensus pass directly on the assembler output, with the model that matches the basecaller - because a mismatched model silently makes the consensus worse.
- CLI: `medaka_consensus -i reads.fq -d draft.fa -o out/ -t 8` (model auto-detected from the basecaller annotation)

medaka is an Oxford Nanopore tool. For PacBio (HiFi/CLR) it is the wrong tool entirely - route to genome-assembly/assembly-polishing.

## The Single Most Important Modern Insight -- A Mismatched Model Silently Degrades; HiFi Must Never Be Fed to Medaka; Prove It on Held-Out Data

medaka is a basecaller-model-specific neural consensus net trained on one exact stack (pore + motor enzyme + speed + basecaller mode + basecaller version). Three consequences:

1. **The model must match the basecaller, and a mismatch fails silently.** Fed reads from a different stack, medaka applies corrections calibrated for an error fingerprint that is not there and misses the real one - the consensus gets WORSE, but medaka exits 0, writes a FASTA, and prints no warning. This is the #1 ONT-polishing footgun, sprung by ordinary acts (re-basecalling with newer Dorado, copying a 2020 model name, polishing a public assembly with the default). Prefer auto-detection (`medaka tools resolve_model --auto_model consensus reads.bam`); treat a stale model name as a reason to re-basecall, not to proceed.
2. **HiFi (and CLR) must never be fed to medaka.** It has no PacBio models; an ONT error-model net "corrects" HiFi toward errors HiFi does not make, and HiFi is already QV40+. If the reads are PacBio, medaka is simply wrong -> genome-assembly/assembly-polishing.
3. **Success is only real on held-out data.** medaka maximizes agreement between the consensus and its input pileup, so grading it on those same reads is circular and always looks good. medaka's "N changes" is a risk signal, not a success signal. Measure with reference-free Merqury QV before vs after on held-out / different-platform k-mers (design deferred to genome-assembly/assembly-polishing).

## What medaka Is For (three modes, same model rule)

| Mode | Input | medaka's role |
|------|-------|---------------|
| Assembly polishing | Flye/Canu draft + ONT reads | raise per-base QV (homopolymer-indel cleanup is the dominant win) |
| Haploid variant calling | ONT reads + reference (microbial, mito, viral) | `medaka_variant` wrapper -> haploid VCF (apply with `bcftools consensus` for a FASTA) |
| Amplicon / viral consensus | tiling-amplicon ONT reads | the non-signal consensus arm of ARTIC fieldbioinformatics / EPI2ME wf-artic |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| ONT-only Flye/Canu assembly | `medaka_consensus`, ONE pass, auto-detected model | model-matched consensus; racon pre-step is obsolete |
| Native bacterial isolate (modified DNA) | `medaka_consensus --bacteria` | bacterial-methylation model fixes methylation-motif errors |
| ONT small-variant (diploid/germline) calling | -> clair3-variants | medaka diploid calling deprecated in v2 (Clair3 surpassed it) |
| Haploid microbial/mito/viral VCF | `medaka_variant` (the renamed haploid wrapper) | still supported in v2 |
| Read-level / human polishing | `dorado polish` | ONT's emerging successor; identical bacterial weights to medaka today |
| PacBio HiFi/CLR | -> genome-assembly/assembly-polishing | medaka has no PacBio models; never ONT-polish HiFi |
| Unsure which basecaller model produced the reads | re-basecall, then auto-detect | a guessed model silently degrades the consensus |

## medaka_consensus Mechanics

The wrapper runs three steps: align (`mini_align`, a thin veil over `minimap2 -x map-ont`), infer (`medaka inference`, the neural net over the pileup), and stitch (`medaka sequence`, regions -> consensus FASTA).

```bash
# Canonical modern usage - model auto-detected from the basecaller annotation in the reads
medaka_consensus -i reads.fastq -d draft.fa -o medaka_out/ -t 8
# medaka_out/consensus.fasta is the polished assembly

# Native bacterial isolate: use the methylation-aware bacterial model
medaka_consensus -i reads.fastq -d draft.fa -o medaka_out/ -t 8 --bacteria

# Resolve / list models (do this when auto-detection cannot pick)
medaka tools resolve_model --auto_model consensus reads.bam
medaka tools list_models
```

medaka runs directly on the assembler (Flye) output as a SINGLE pass - do NOT pre-run Racon (contemporary models are trained on raw assembler output; v2 removed the bundled racon wrapper) and do NOT run medaka twice (iteration was racon's role; a second pass risks flipping correct bases).

### Haploid variant calling (v2 names)

medaka_variant emits a VCF only (no consensus FASTA); apply it to the reference with `bcftools consensus` to get a haploid consensus sequence.

```bash
# Wrapper form (renamed from medaka_haploid_variant in v2) - haploid samples only
medaka_variant -i reads.fastq -r reference.fa -o variant_out/

# Manual form - note v2 subcommand names and the hdf -> ref -> out argument order
minimap2 -ax map-ont reference.fa reads.fq | samtools sort -o aln.bam && samtools index aln.bam
medaka inference aln.bam probs.hdf --model r1041_e82_400bps_sup_variant_v5.0.0
medaka vcf probs.hdf reference.fa variants.vcf

# Optional: turn the VCF into a haploid consensus FASTA
bgzip variants.vcf && tabix -p vcf variants.vcf.gz
bcftools consensus -f reference.fa variants.vcf.gz > consensus.fasta
```

## Per-Method Failure Modes

### Silent model mismatch
**Trigger:** running medaka with a model that does not match the basecaller chemistry/version. **Mechanism:** the net corrects toward the wrong error fingerprint. **Symptom:** lower held-out QV; medaka exits 0 with no warning. **Fix:** auto-detect from the BAM; if forced to pick, derive from the actual basecaller and confirm in `list_models`; treat a stale name as a reason to re-basecall.

### HiFi fed to medaka
**Trigger:** polishing a PacBio assembly with medaka. **Mechanism:** ONT-only error model, no PacBio support, on already-QV40+ data. **Symptom:** degraded/homogenized consensus. **Fix:** do not; route to genome-assembly/assembly-polishing.

### Racon-first off-distribution
**Trigger:** running Racon before medaka out of habit. **Mechanism:** contemporary models are trained on raw assembler output; racon-polished input is off the training distribution. **Symptom:** no gain or mild harm. **Fix:** run medaka directly on the Flye output; one pass.

### Missing plasmid poisons the chromosome
**Trigger:** an assembly missing a small replicon (~80% identical to a chromosomal region). **Mechanism:** the absent plasmid's reads misalign onto the chromosome, and medaka "corrects" toward that spurious evidence. **Symptom:** clustered changes that introduce real errors. **Fix:** make the assembly structurally complete first; inspect medaka's changes for clustering (clustered = mapping artifact, not scattered homopolymer fixes).

### Validating on the polishing reads
**Trigger:** judging the polish by medaka's change count or by re-mapping the same reads. **Mechanism:** medaka optimizes agreement with its input pileup. **Symptom:** "improvement" that is circular. **Fix:** reference-free Merqury QV before vs after on held-out / different-platform k-mers.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| 1 medaka pass | medaka README | a single trained-model pass; iteration was racon's role, extra passes flip correct bases |
| model must match basecaller version | medaka model design | mismatch silently degrades; the #1 ONT-polishing error |
| inference threads ~2 | medaka inference behavior | the net is GPU-bound and scales poorly past ~2 CPU threads |
| HiFi QV40+ already | EBP/HiFi baseline | nothing for an ONT consensus net to gain; only harm |
| measure with held-out Merqury QV | Rhie 2020 | the only honest, reference-free before/after instrument |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `medaka consensus` not found / wrong args | v1 subcommand renamed | use `medaka inference` (or the `medaka_consensus` wrapper) |
| `medaka stitch` / `medaka variant` fail | v1 names | `medaka sequence` / `medaka vcf` |
| Polished assembly worse than draft | model mismatch | auto-detect the model; re-basecall if the model is stale |
| medaka errors on PacBio reads | no PacBio models | route to genome-assembly/assembly-polishing |
| Clustered, suspicious changes | missing/mis-structured contig in the draft | complete the assembly first; filter to high-identity alignments |
| Looking for diploid SNP calling | deprecated in v2 | use clair3-variants |

## References

- medaka. Oxford Nanopore Technologies. https://github.com/nanoporetech/medaka (no journal paper; cite the repository).
- Zheng Z, Li S, Su J, Leung AW, Lam TW, Luo R. 2022. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling (Clair3). *Nat Comput Sci* 2:797-803.
- Vaser R, Sović I, Nagarajan N, Šikić M. 2017. Fast and accurate de novo genome assembly from long uncorrected reads (Racon). *Genome Res* 27:737-746.
- Wick RR, Judd LM, Holt KE. 2023. Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing. *PLoS Comput Biol* 19(3):e1010905.
- Rhie A, Walenz BP, Koren S, Phillippy AM. 2020. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biol* 21:245.
- Wick RR. 2024. Medaka v2: progress and potential pitfalls. https://rrwick.github.io/2024/10/17/medaka-v2.html (blog; source of the missing-plasmid footgun).

## Related Skills

- basecalling - The basecaller model+version medaka's model must match
- clair3-variants - ONT small-variant (diploid/germline) calling; medaka diploid is deprecated
- long-read-alignment - minimap2 map-ont, the alignment medaka's mini_align wraps
- genome-assembly/assembly-polishing - Polishing strategy authority (HiFi doctrine, hybrid tiers, Merqury QV design)
- genome-assembly/long-read-assembly - Produces the Flye draft medaka polishes
- genome-assembly/assembly-qc - Merqury QV / BUSCO before-vs-after measurement
