---
name: bio-variant-calling-structural-variant-calling
description: Call structural variants (>=50 bp deletions, insertions, inversions, duplications, translocations) from short- or long-read data by reconstructing four orthogonal signals (discordant pairs, split reads via the SA tag, read depth, local assembly). Covers Manta, DELLY, LUMPY/smoove, GRIDSS2, SvABA for short reads and Sniffles2, cuteSV, pbsv, dipcall/PAV for long reads, each mapped to the signals it fuses and the blind spots that follow. Use when choosing an SV caller from its signal set and failure modes, decoding the SVLEN-sign / symbolic-vs-BND / CIPOS VCF representation minefield, force-genotyping a cohort matrix instead of unioning discovery VCFs, merging populations with sequence-aware Truvari vs position-only SURVIVOR, parameterizing a Truvari benchmark, or deciding when short-read insertion recall forces a switch to long reads. Not for pure copy-number dosage (see copy-number/cnvkit-analysis).
tool_type: cli
primary_tool: manta
---

## Version Compatibility

Reference examples tested with: Manta 1.6+, DELLY 1.2+, GRIDSS 2.13+, smoove 0.2+, SvABA 1.1+, bcftools 1.19+, samtools 1.19+, SURVIVOR 1.0.7+, Truvari 4.0+, Sniffles2 2.2+, cuteSV 2.0+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: DELLY 1.x and Truvari 4.x renamed flags/defaults from earlier lines (Truvari's alt-sequence-similarity param was `--pctsim`, now `--pctseq`); always confirm against `--help` before trusting old muscle memory.

# Structural Variant Calling

**"Call structural variants from my WGS data"** -> Reconstruct large rearrangements (DEL, INS, INV, DUP, BND) that no single read reports directly, by fusing evidence from how alignments break, disagree, or deviate from expectation.
- CLI (short read): `configManta.py` (Manta), `delly call`, `smoove call` (LUMPY), `gridss`, `svaba run`
- CLI (long read): `sniffles`, `cuteSV`, `pbsv`, `dipcall.kit`/`pav` (assembly-based)

## The governing principle: SV discovery is signal reconstruction, not base-calling

An SV caller is not a pileup genotyper. SVs are never observed directly by a read; they are *inferred* from four orthogonal signals, and each caller fuses only a subset. Its blind spots follow directly from which signals it omits (Mahmoud 2019 *Genome Biol* 20:246 is the canonical signal-to-blind-spot map).

1. Discordant read pairs (RP) - a pair whose insert size or orientation departs from the library distribution. A DEL inflates the apparent insert; a tandem DUP gives everted (`+/+`/`-/-`) pairs; an INV flips one strand; an interchromosomal event splits mates across chromosomes. RP never sees the junction base, so it localizes a breakpoint only to fragment-size uncertainty (~+/-300-500 bp).
2. Split reads / soft-clips (SR) - one read aligning in two pieces to non-contiguous positions, recorded via the BAM SA (supplementary alignment) tag. The clip point IS the breakpoint at base resolution (~+/-0-10 bp). Requires a read that straddles the junction with unique anchor on both sides - impossible inside a repeat longer than the read or for an INS larger than (read length - anchor).
3. Read depth (RD) - copy-number shifts local coverage (het DEL ~0.5x, DUP ~1.5x, hom DEL ~0). The ONLY signal for large CNVs whose breakpoints fall in unmappable repeats; blind to balanced events (INV, balanced BND) that leave dosage unchanged. Resolution is bin-size-limited (100 bp-1 kb), the worst of any signal.
4. Local assembly (AS) - reconstruct the breakpoint-spanning contig de novo from clipped/discordant reads, then realign it. Recovers the exact junction sequence, microhomology, and inserted bases; resolves events no single read spans. Most powerful and most expensive; separates GRIDSS/Manta/SvABA from LUMPY.

| Caller | RP | SR | RD | Assembly | Primary blind spot |
|--------|----|----|----|----------|--------------------|
| Manta | Yes | Yes | filter/score only | Yes (targeted, breakend-local) | Large INS beyond breakend-local contig; balanced events in repeats |
| DELLY | Yes | Yes (SR realign for refinement) | Yes (`delly cnv` module) | No de novo | INS (only small, SR-anchored); precision below assembly callers |
| LUMPY / smoove | Yes | Yes | No | No | INS entirely (no representation for novel inserted sequence) |
| GRIDSS2 | Yes | Yes | via GRIPSS/PURPLE downstream | Yes (genome-wide positional de Bruijn graph) | Pure RD-only CNVs with no junction reads; long-repeat interiors |
| SvABA | Yes | Yes | No | Yes (genome-wide SGA local assembly) | Large INS still assembly-limited; multi-kb DEL/DUP rely on RP linkage |
| CNVnator | No | No | Yes (mean-shift on RD bins) | No | Everything balanced; no breakpoint resolution; no small SVs |

### The real tradeoff is sensitivity vs breakpoint resolution (not vs specificity)

- RP-heavy calling maximizes SENSITIVITY (a pair spanning a junction is easy, works at low depth) but yields IMPRECISE breakpoints (+/-hundreds bp, `CIPOS=-289,301`).
- SR / assembly calling maximizes RESOLUTION (base-precise, `CIPOS=0,0`) but LOSES sensitivity wherever a read cannot cleanly straddle the junction (repeats, low depth, large INS).

Below ~30x the SR signal thins (fewer reads straddle any junction) and callers silently drift into the low-resolution RP-only regime. A caller demanding SR confirmation reports beautiful breakpoints and misses the repeat-mediated events that matter clinically; a caller accepting RP-only calls finds more but hands back wide confidence intervals.

## Caller selection

| Caller | Signals fused | Best when | Fails / do NOT use when |
|--------|---------------|-----------|--------------------------|
| Manta | RP+SR+targeted AS | Default germline/somatic; fastest assembly-capable caller (sub-hour on 30x) | WES/panel WITHOUT `--exome` (high-depth filter silently drops true SVs); large INS beyond local assembly |
| DELLY | RP+SR (+RD CNV module) | Cohorts (site-list-then-regenotype); INV/BND | General large-INS calling; needs base-precise breakpoints in repeats |
| smoove (LUMPY) | RP+SR (probabilistic) | Simple, sane-default germline pipeline at low-moderate depth | ANY INS-critical pipeline (structurally zero INS recall) |
| GRIDSS2 | RP+SR+genome-wide AS | Highest precision; somatic; single-breakends (viral integration, centromeric) | Using raw VCF as a callset (it is a breakpoint GRAPH - run GRIPSS/PURPLE/LINX); need speed |
| SvABA | RP+SR+genome-wide AS | The 20-300 bp indel/SV "twilight zone"; templated-insertion detection in cancer | Large balanced events; when Manta speed suffices |
| CNVnator / CNVpytor | RD only | Dosage CNVs with breakpoints in unmappable repeats | Balanced SVs; any breakpoint-precision need (see copy-number/cnvkit-analysis) |
| Sniffles2 / cuteSV / pbsv | long-read alignment | INS, repeat-mediated SVs, phased SVs, population `.snf` merging | Only short reads available |

Consensus recipe: run DELLY + Manta + SvABA + GRIDSS2 (the latter GRIPSS-filtered first, since its raw VCF is a breakpoint graph, not a callset) and require >=2/4 agreement - singletons are enriched for false positives, so intersection trades a little recall for a large precision gain. Caveat: this caps INS recall (LUMPY-class blind spots drag the intersection down); for INS-heavy work prefer one assembly caller + a graph genotyper, or long reads. Methods evolve - verify current ensemble practice against tool docs before committing.

## SV detection limits by platform

| SV type | Short read | Long read | Key limitation |
|---------|-----------|-----------|----------------|
| Deletion | Good (RP+SR) | Excellent | Short reads miss DELs buried in repeats |
| Duplication | Moderate (RP+RD) | Good | Tandem vs dispersed unreliable with short reads |
| Inversion | Moderate (RP) | Good | Breakpoints in repeats cause false negatives |
| Insertion | Poor (~30-50% recall) | Excellent (~90%+) | Physics limit: no short read spans an INS > read length |
| Translocation | Moderate (discordant) | Good | High FP rate near centromeres/telomeres |
| Complex/nested | Poor | Good (assembly) | Overlapping SVs confound short-read signals |

Insertions are a physics limit, not a tuning failure. Placing and sizing an INS requires reads spanning both junctions plus the novel bases; when the INS exceeds read length (Illumina 150 bp), no read spans it and the inserted sequence is unrecoverable without assembly - and short-read assembly fails again when the insertion is a repeat (mobile element, VNTR). GIAB HG002 Tier 1 contains MORE INS than DEL (7,281 INS vs 5,464 DEL >=50 bp; Zook 2020 *Nat Biotechnol* 38:1347), yet INS are the hardest short-read class. Ebert 2021 (*Science* 372:eabf7117) found 68% of 107,590 assembly-discovered SVs were missed by short reads. Do NOT benchmark short-read INS against a long-read truth set and blame the caller - ~30-50% is the ceiling by construction.

## The VCF representation minefield

The VCF spec represents SVs three ways and callers disagree on all three - this is where careful people lose days.

- SVTYPE - the class (DEL/INS/DUP/INV/BND).
- END (INFO) - the other breakpoint on the same chromosome. POS is the last unaffected anchor base, so a DEL at POS=1000, END=2000 removes 1001-2000. Off-by-one here corrupts every length calc. (VCF 4.4 began deprecating INFO/END for symbolic alleles; bcftools/GATK still emit and consume the classic field.)
- SVLEN - SIGNED length: DEL is NEGATIVE, INS/DUP positive by historical convention (Manta, DELLY follow it; some tools emit absolute values). NEVER filter on raw SVLEN without `ABS()` or a naive `SVLEN >= 50` drops every deletion.
- CIPOS / CIEND - confidence intervals encoding the RP-vs-SR resolution story: an SR-resolved breakpoint is `CIPOS=0,0`, an RP-only one is `CIPOS=-289,301`. This is the ground truth for "how much do I trust this breakpoint" and exactly what a population merger must respect.
- IMPRECISE (flag) - set when the breakpoint is RP/depth-derived, not junction-resolved; its absence implies PRECISE. An all-IMPRECISE callset cannot be merged tightly.

Breakend (BND) notation - the part everyone gets wrong. A rearrangement whose two sides are not a simple same-chromosome span is written as PAIRED BND records linked by `MATEID`, with an ALT string whose bracket direction encodes the join orientation:
```
2   321681  bnd_V  T  ]13:123456]T   MATEID=bnd_U
13  123456  bnd_U  A  A[2:321681[    MATEID=bnd_V
```
One reciprocal translocation is 2 BND records; chromothripsis is a graph of dozens. The SAME biological inversion can appear as `<INV>` in Manta, as 2+ BND records in GRIDSS, and as a DEL+DUP artifact pair in a naive RD caller - a merger that does not understand this triple-counts or drops the event.

## Manta

Manta builds a genome-wide breakend association graph from RP+SR, then does targeted local assembly of each candidate breakend and realigns the contig for base-resolution breakpoints. `candidateSmallIndels.vcf.gz` is the recommended indel-candidate input for Strelka2.

```bash
configManta.py --bam sample.bam --referenceFasta reference.fa --runDir manta_run
manta_run/runWorkflow.py -j 8
# results/variants/: diploidSV.vcf.gz (germline), candidateSV.vcf.gz (unscored superset),
#   candidateSmallIndels.vcf.gz (Strelka2 input)
```

WES/panel and RNA need explicit modes - default depth filtering assumes WGS-uniform coverage and silently drops true SVs at high-depth targeted loci:
```bash
# --exome disables the high-depth filter; --rna handles split alignments across splice junctions
configManta.py --bam sample.bam --referenceFasta reference.fa --exome \
    --callRegions regions.bed.gz --runDir manta_exome
```

Tumor-normal somatic mode adds `somaticSV.vcf.gz`:
```bash
configManta.py --tumorBam tumor.bam --normalBam normal.bam \
    --referenceFasta reference.fa --runDir manta_somatic
manta_somatic/runWorkflow.py -j 8
```

## DELLY

BCF output by default; the scalable cohort pattern is site-list-then-regenotype (Section below), not one big multi-BAM call.

```bash
delly call -g reference.fa -o sv_calls.bcf sample.bam
bcftools view sv_calls.bcf > sv_calls.vcf

# Per-type calling (INS recovers only small, SR-anchored insertions)
delly call -t DEL -g ref.fa -o deletions.bcf sample.bam   # also DUP, INV, BND, INS

# Somatic: call tumor+normal, then classify with a sample sheet
delly call -g reference.fa -o svs.bcf tumor.bam normal.bam
printf 'tumor\ttumor\nnormal\tcontrol\n' > samples.tsv
delly filter -f somatic -o somatic_svs.bcf -s samples.tsv svs.bcf
```

## smoove (LUMPY, with genotyping + depth annotation)

Nobody runs raw `lumpyexpress` anymore. smoove wraps LUMPY + svtyper + duphold with sane defaults and a cohort workflow. duphold annotates each DEL/DUP with three depth fold-changes: over the whole event (`DHFC`), over GC/mappability-matched bins (`DHBFC`), and over the immediate flanks (`DHFFC`) - a cheap, high-value filter for RP-only calls lacking depth support. duphold's own recommended false-positive filters are `DHFFC < 0.7` for deletions (flanking fold-change, its most robust metric) and `DHBFC > 1.3` for duplications.

```bash
smoove call --name sample --fasta reference.fa --outdir smoove_out -p 8 sample.bam
# smoove_out/sample-smoove.genotyped.vcf.gz

# Filter RP-only false positives that lack depth corroboration (duphold's recommended fields)
bcftools view -i '(SVTYPE!="DEL" && SVTYPE!="DUP") || (SVTYPE="DEL" && FMT/DHFFC<0.7) || (SVTYPE="DUP" && FMT/DHBFC>1.3)' \
    smoove_out/sample-smoove.genotyped.vcf.gz > smoove.dhfc.vcf
```

## GRIDSS: the raw VCF is a breakpoint graph, not a callset

GRIDSS2 is the only short-read caller doing genome-wide break-end assembly (positional de Bruijn graph) before calling, and the only one reporting single break-ends (SGL) where just one side maps - viral integrations, centromeric/telomeric junctions. But its raw VCF is deliberately noisy low-level breakpoints; it is NOT usable directly. The intended chain is GRIDSS -> GRIPSS (filtering/linkage) -> PURPLE (purity/ploidy/copy-number) -> LINX (interprets breakpoints into chromothripsis, breakage-fusion-bridge, fusions).

```bash
# --assembly must be a writable path; reference needs .fai and .dict
gridss --reference reference.fa --output gridss_raw.vcf \
    --assembly gridss_assembly.bam --threads 8 sample.bam

# Somatic: GRIDSS2 tumor+normal, then GRIPSS filtering
gridss --reference reference.fa --output gridss_raw.vcf --assembly asm.bam \
    --labels normal,tumor --threads 8 normal.bam tumor.bam
# GRIPSS (launched as `java -jar gripss.jar` in practice) also requires PON inputs
# (-pon_sgl_file, -pon_sv_file, -known_hotspot_file); confirm the full flag set with its docs.
gripss -sample tumor -reference normal -ref_genome reference.fa \
    -ref_genome_version 38 -pon_sgl_file sgl.pon -pon_sv_file sv.pon \
    -vcf gridss_raw.vcf -output_dir gripss_out/
```

## Genotyping is not discovery (do NOT union discovery VCFs into a population matrix)

Discovery answers "is there an SV here, and what/where?"; genotyping answers "what is each sample's GT (0/0, 0/1, 1/1) here?". These are different statistical problems. A sample recorded 0/0 in a discovery-VCF union may simply not have had that event DISCOVERED in it - a false "missing", not a true reference genotype. Building an allele-frequency-quality cohort matrix requires FORCE-GENOTYPING every sample at every discovered site.

Correct cohort workflow: discover per-sample -> merge to a non-redundant SITE list -> re-genotype ALL samples at ALL sites -> merge genotypes.

| Genotyper | Model | Use when |
|-----------|-------|----------|
| svtyper (in SpeedSeq, Chiang 2015) | Bayesian RP+SR at known breakpoints | LUMPY/smoove sites; DEL/DUP/INV/BND (NOT insertions) |
| Paragraph (Chen 2019) | realign reads to a per-SV sequence graph (ref+alt paths) | Modern short-read choice; genotypes INS (alt path contains the inserted bases) |
| GraphTyper2 (Eggertsson 2019) | joint SNV+SV over a pangenome graph | N in the tens of thousands (genotyped 49,962 Icelanders) |

```bash
# DELLY's native squared-off pattern: per-sample call -> merge SITES -> regenotype at sites
delly call -g ref.fa -o s1.bcf s1.bam        # ... one per sample
delly merge -o sites.bcf s1.bcf s2.bcf s3.bcf
delly call -g ref.fa -v sites.bcf -o s1.geno.bcf s1.bam   # regenotype each at the union sites
bcftools merge -m id -Ob -o cohort.bcf s1.geno.bcf s2.geno.bcf s3.geno.bcf
```

## Population merging: the merge parameters ARE the result

"The same SV" across samples/callers is a fuzzy concept because breakpoints disagree by CIPOS margins. Position-only merging inflates allele frequency by up to 2.2x versus sequence-aware merging, because it collapses distinct nearby alleles into one (English 2022 *Genome Biol* 23:271). For any AF-dependent analysis (constraint, association, catalogs), always report the merger and its parameters.

| Merger | Matching | Use when |
|--------|----------|----------|
| SURVIVOR (Jeffares 2017) | position + type + strand agreement, min-callers; NOT sequence-aware | Fast caller-consensus on ONE sample; 1000 bp `max_dist` cheerfully merges two different 300 bp DELs 800 bp apart |
| Jasmine (Kirsche 2023) | (breakpoint, length) proximity graph via KD-tree + constrained Kruskal | Long-read cohort merging at population scale |
| Truvari collapse (English 2022) | SEQUENCE-aware (refdist + size + alt-sequence similarity) | Any AF work or population catalog where allelic diversity must be preserved |

```bash
# SURVIVOR: max_dist=1000 min_callers=2 type_agree=1 strand_agree=1 est_dist=0 min_size=50
ls manta.vcf delly.vcf gridss.vcf smoove.vcf > vcf_list.txt
SURVIVOR merge vcf_list.txt 1000 2 1 1 0 50 merged.vcf   # single-sample caller consensus only
```

## Benchmarking: an SV F1 is meaningless without its Truvari parameters

Every SV F1/recall/precision figure is a function of the matching parameters; papers routinely report incomparable numbers. A Truvari `bench` true positive must match a truth variant under ALL of:

| Flag | Default | What loosening it does |
|------|---------|------------------------|
| `--refdist` | 500 | larger -> RP-only IMPRECISE calls start matching |
| `--pctsize` | 0.70 | lower -> size-sloppy calls pass (reciprocal size similarity) |
| `--pctseq` (was `--pctsim`) | 0.70 | `0` disables alt-sequence checking entirely - the quiet trick that inflates short-read INS scores |
| `--sizemin` / `--sizemax` | 50 / - | choosing the window can hide a caller weak at one size class |

```bash
# Report EVERY parameter; --pctseq 0 would make this uninterpretable
truvari bench -b truth.vcf.gz -c calls.vcf.gz -o bench_out/ --passonly \
    --refdist 500 --pctsize 0.70 --pctseq 0.70 --sizemin 50
```

Three escalating bars are routinely conflated: (1) event detection (loose, RP-only callers pass), (2) breakpoint accuracy (tight `--refdist`, only SR/assembly pass; matters at exon/splice boundaries), (3) genotype accuracy (add genotype-aware comparison; a caller can detect an event yet call het-as-hom, fatal for Mendelian work). Also stratify by region: run GIAB-CMRG (Wagner 2022 *Nat Biotechnol* 40:672), not just Tier 1 - CMRG covers the medically relevant repetitive genes Tier 1 EXCLUDES, and false duplications in GRCh37/38 (e.g. *CBS*, *KCNE1*) cause reference-specific misses that masking raised from 8% to 100% recall.

## Filter and annotate

```bash
bcftools view -i 'ABS(SVLEN) >= 50' svs.vcf > svs.min50.vcf   # ABS() is mandatory (DEL sign); see examples/svlen_sign_filter.py
bcftools view -i 'SVTYPE="DEL"' svs.vcf > deletions.vcf        # or INS/INV/DUP/BND
bcftools view -f PASS svs.vcf > svs.pass.vcf

AnnotSV -SVinputFile svs.vcf -genomeBuild GRCh38 -outputFile annotated_svs
# gene overlap, DGV/gnomAD-SV population AF, ClinVar pathogenicity
```

## Long reads: near-direct observation

A long read (HiFi ~15-25 kb, ONT tens of kb to Mb) physically spans the SV and both flanks in one molecule, converting inference-over-fragments into near-direct observation: INS become trivial (the read carries the inserted bases), repeat-mediated SVs resolve, and heterozygous variants on one read are natively phased. Switch to long reads when INS/complex/repeat SVs or phased haplotyping matter more than per-sample cost.

| Caller | Approach | Best for |
|--------|----------|----------|
| Sniffles2 (Smolka 2024) | signature + per-sample `.snf` population merge | ONT/HiFi general; mosaic/low-VAF SVs; linear-in-N cohorts |
| cuteSV (Jiang 2020) | signature clustering + refinement | Highest recall on noisy ONT |
| pbsv | official PacBio, tandem-repeat-aware | HiFi (pair with pbmm2 alignments) |
| SVIM (Heller 2019) | reports origin AND destination of duplications | tandem-vs-interspersed DUP discrimination |
| dipcall (Li 2018) / PAV (Ebert 2021) | assembly-vs-reference from phased haplotype assemblies | highest-quality callsets and truth sets |
| Severus (Keskus 2025) | phased somatic breakpoint graph | complex somatic rearrangements from long reads |

```bash
# minimap2/pbmm2 -> sort -> caller. Sniffles2 population design: per-sample .snf, then combine
sniffles --input sample.bam --reference reference.fa --snf sample.snf
sniffles --input s1.snf s2.snf s3.snf --vcf population.vcf   # combine scales linearly in N
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Default Manta on WES/panel loses true SVs | high-depth filter assumes WGS-uniform coverage | add `--exome` |
| Cohort "0/0" genotypes wrong; AF too low | took a UNION of discovery VCFs | force-genotype at merged sites (Paragraph/GraphTyper2/`delly -v`) |
| AF up to 2.2x too high vs another catalog | position-only merge (SURVIVOR-1000) | use sequence-aware Truvari for AF work |
| `ABS(SVLEN)>=50` filter drops all deletions | filtered raw SVLEN (DEL is negative) | always wrap in `ABS()` |
| Zero INS recall | LUMPY/smoove has no INS representation | use an assembly caller or long reads |
| Short-read INS recall "only ~40%" | benchmarked vs a long-read truth set | physics limit, not caller fault; ~30-50% is the ceiling |
| GRIDSS raw VCF looks like noise | it is a breakpoint GRAPH, not a callset | run GRIPSS -> PURPLE -> LINX |
| `-H 1` haplotype consensus is chimeric | genotypes were unphased | phase first (see variant-calling/consensus-sequences) |
| Truvari F1 not reproducible | reported without parameters | state refdist/pctsize/pctseq/sizemin and whether GT was compared |

## Related Skills

- copy-number/cnvkit-analysis - read-depth CNV detection for dosage changes with breakpoints in unmappable repeats (complements junction-based SV callers)
- long-read-sequencing/structural-variants - full long-read SV pipelines (Sniffles2, cuteSV, pbsv, assembly-based)
- variant-calling/consensus-sequences - applying variants to a reference; phasing before `-H` haplotype extraction; why symbolic SVs are not consensus-able
- variant-calling/vcf-manipulation - view, query, and reshape SV VCFs
- variant-calling/filtering-best-practices - general VCF filtering principles applicable to SV callsets
- variant-calling/variant-annotation - functional annotation of SVs (gene overlap, population AF, pathogenicity)
- alignment-files/alignment-filtering - BAM preparation and quality filtering before SV calling

## References

- Mahmoud M, Gobet N, Cruz-Davalos DI, Mounier N, Dessimoz C, Sedlazeck FJ. Structural variant calling: the long and the short of it. 2019 *Genome Biology* 20:246.
- Chen X, Schulz-Trieglaff O, Shaw R, Barnes B, Schlesinger F, Kallberg M, Cox AJ, Kruglyak S, Saunders CT. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. 2016 *Bioinformatics* 32:1220-1222.
- Rausch T, Zichner T, Schlattl A, Stutz AM, Benes V, Korbel JO. DELLY: structural variant discovery by integrated paired-end and split-read analysis. 2012 *Bioinformatics* 28:i333-i339.
- Layer RM, Chiang C, Quinlan AR, Hall IM. LUMPY: a probabilistic framework for structural variant discovery. 2014 *Genome Biology* 15:R84.
- Cameron DL, Schroder J, Penington JS, Do H, Molania R, Dobrovic A, Speed TP, Papenfuss AT. GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly. 2017 *Genome Research* 27:2050-2060.
- Cameron DL, Baber J, Shale C, Valle-Inclan JE, Besselink N, van Hoeck A, et al. GRIDSS2: comprehensive characterisation of somatic structural variation using single breakend variants and structural variant phasing. 2021 *Genome Biology* 22:202.
- Wala JA, Bandopadhayay P, Greenwald NF, O'Rourke R, Sharpe T, Stewart C, et al. SvABA: genome-wide detection of structural variants and indels by local assembly. 2018 *Genome Research* 28:581-591.
- Abyzov A, Urban AE, Snyder M, Gerstein M. CNVnator: an approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing. 2011 *Genome Research* 21:974-984.
- Chiang C, Layer RM, Faust GG, Lindberg MR, Rose DB, Garrison EP, Marth GT, Quinlan AR, Hall IM. SpeedSeq: ultra-fast personal genome analysis and interpretation. 2015 *Nature Methods* 12:966-968. (introduces the svtyper genotyper)
- Chen S, Krusche P, Dolzhenko E, Sherman RM, Petrovski R, Schlesinger F, et al. Paragraph: a graph-based structural variant genotyper for short-read sequence data. 2019 *Genome Biology* 20:291.
- Eggertsson HP, Kristmundsdottir S, Beyter D, Jonsson H, Skuladottir A, Hardarson MT, et al. GraphTyper2 enables population-scale genotyping of structural variation using pangenome graphs. 2019 *Nature Communications* 10:5402.
- Jeffares DC, Jolly C, Hoti M, Speed D, Shaw L, Rallis C, Balloux F, Dessimoz C, Bahler J, Sedlazeck FJ. Transient structural variations have strong effects on quantitative traits and reproductive isolation in fission yeast. 2017 *Nature Communications* 8:14061. (introduces SURVIVOR)
- Kirsche M, Prabhu G, Sherman R, Ni B, Battle A, Aganezov S, Schatz MC. Jasmine and Iris: population-scale structural variant comparison and analysis. 2023 *Nature Methods* 20:408-417.
- English AC, Menon VK, Gibbs RA, Metcalf GA, Sedlazeck FJ. Truvari: refined structural variant comparison preserves allelic diversity. 2022 *Genome Biology* 23:271.
- Zook JM, Hansen NF, Olson ND, Chapman L, Mullikin JC, Xiao C, et al. A robust benchmark for detection of germline large deletions and insertions. 2020 *Nature Biotechnology* 38:1347-1355. (GIAB HG002 SV Tier 1)
- Wagner J, Olson ND, Harris L, McDaniel J, Cheng H, Fungtammasan A, et al. Curated variation benchmarks for challenging medically relevant autosomal genes. 2022 *Nature Biotechnology* 40:672-680. (GIAB-CMRG)
- Li H, Bloom JM, Farjoun Y, Fleharty M, Gauthier L, Neale B, MacArthur D. A synthetic-diploid benchmark for accurate variant-calling evaluation. 2018 *Nature Methods* 15:595-597. (dipcall/syndip)
- Ebert P, Audano PA, Zhu Q, Rodriguez-Martin B, Porubsky D, Bonder MJ, et al. Haplotype-resolved diverse human genomes and integrated analysis of structural variation. 2021 *Science* 372:eabf7117. (PAV; 68% of SVs missed by short reads)
- Sedlazeck FJ, Rescheneder P, Smolka M, Fang H, Nattestad M, von Haeseler A, Schatz MC. Accurate detection of complex structural variations using single-molecule sequencing. 2018 *Nature Methods* 15:461-468. (Sniffles v1 + NGMLR)
- Smolka M, Paulin LF, Grochowski CM, Horner DW, Mahmoud M, Behera S, et al. Detection of mosaic and population-level structural variants with Sniffles2. 2024 *Nature Biotechnology*. doi:10.1038/s41587-023-02024-y.
- Jiang T, Liu Y, Jiang Y, Li J, Gao Y, Cui Z, et al. Long-read-based human genomic structural variation detection with cuteSV. 2020 *Genome Biology* 21:189.
- Heller D, Vingron M. SVIM: structural variant identification using mapped long reads. 2019 *Bioinformatics* 35:2907-2915.
- Keskus AG, et al. Severus detects somatic structural variation and complex rearrangements in cancer genomes using long-read sequencing. 2025 *Nature Biotechnology*. doi:10.1038/s41587-025-02618-8.
- pbsv - PacBio structural variant caller (no dedicated publication): github.com/PacificBiosciences/pbsv
- smoove (no dedicated publication): github.com/brentp/smoove
- Pedersen BS, Quinlan AR. Duphold: scalable, depth-based annotation and curation of high-confidence structural variant calls. 2019 *GigaScience* 8(4):giz040. (source of the DHFFC<0.7 / DHBFC>1.3 depth filters)
