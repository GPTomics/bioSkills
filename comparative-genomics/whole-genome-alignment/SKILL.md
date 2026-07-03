---
name: bio-comparative-genomics-whole-genome-alignment
description: Build whole-genome alignments using Progressive Cactus (Armstrong 2020 reference-free clade-level WGA), Minigraph-Cactus (Hickey 2024 pangenome-aware), LASTZ chain/net (UCSC pipeline), MUMmer4 (Marçais 2018 pairwise), minimap2 -x asm5/10/20 (Li 2018 fast pairwise), AnchorWave (Song 2022 WGD-aware), and Mauve / progressiveMauve (bacterial). Operates the HAL toolkit (Hickey 2013) for downstream extraction including halSynteny, halLiftover, halBranchMutations, and hal2maf. Use when constructing multi-species alignments for comparative-annotation projection (TOGA), synteny detection, conservation analyses (phyloP / PhastCons), or pangenome graph construction; selecting between reference-free (Cactus) and reference-anchored (LASTZ chains/nets) approaches; tuning sensitivity for closely vs distantly related genomes; or producing HAL files for genome-wide downstream tools.
tool_type: cli
primary_tool: Cactus
---

## Version Compatibility

Reference examples tested with: Progressive Cactus 2.9.1+ (ComparativeGenomicsToolkit/cactus; Armstrong 2020 Nature 587:246), Minigraph-Cactus (Hickey 2024 Nat Biotech 42:663; bundled with Cactus 2.5+), HAL toolkit 2.3+ (Hickey 2013 Bioinformatics 29:1341), LASTZ 1.04.22+, UCSC kentUtils for chain/net (Kent 2003 PNAS 100:11484), MUMmer 4.0.0+, minimap2 2.28+, AnchorWave 1.2.5+, progressiveMauve 2.4.0+, sibeliaz 1.2.5+, winnowmap 2.03+ (Jain 2022 Nat Methods 19:705). Toil workflow runner 6.0+ for Cactus on HPC/cloud.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `cactus --help`, `cactus-pangenome --help`, `halStats --help`, `lastz --version`, `nucmer --version`, `minimap2 --version`
- Python: `pip show toil`, `toil --version`

If code throws `Toil workflow restart failure`, `HAL file corrupted`, `WDL workflow missing`, the Cactus pipeline is Toil-based and requires careful checkpointing; failed runs must be restarted with `--restart`. HAL file versions differ across hal-toolkit releases; pin the version that produced the file.

# Whole Genome Alignment

**"Align these multiple genomes at the base-pair level"** -> Choose between **reference-free** progressive alignment (Cactus / Minigraph-Cactus: produces HAL, no privileged reference) and **reference-anchored** pairwise alignment (LASTZ chains/nets, MUMmer, minimap2: one genome is the reference, queries align to it). The fundamental tradeoff is **scale vs structure**: pairwise pipelines scale linearly per pair but lose multi-way relationships; progressive Cactus scales linearly with a tree but quadratically without and produces ancestrally-coherent alignments. For comparative genomics at vertebrate / mammal scale, Cactus is now the standard substrate (Zoonomia, Christmas 2023 Science 380:eabn3943; Bird10000 Genomes); for pangenome graph construction, Minigraph-Cactus (Hickey 2024) is the production pipeline.

- CLI: `cactus jobStore seqFile.txt output.hal --binariesMode local` -- reference-free progressive WGA
- CLI: `cactus-pangenome --reference ref name --vcf` -- pangenome graph from genomes
- CLI: `lastz target.fa[multiple] query.fa` then UCSC chain-net pipeline -- pairwise to a reference
- CLI: `minimap2 -ax asm5 ref.fa query.fa | samtools sort` -- fast pairwise for closely related
- CLI: `nucmer --maxmatch ref.fa query.fa` then `dnadiff` -- MUMmer4 pairwise
- CLI: `anchorwave proali --ploidy 4` -- WGD-aware sequence-level synteny alignment

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| Progressive Cactus (Armstrong 2020 Nature 587:246) | LASTZ pairwise -> CAF graph -> BAR -> ancestral assembly per node | HAL multi-genome alignment | Reference-free; scales to 1000s of genomes with phylogenetic tree; ancestral assembly at every internal node | Quadratic without tree; ALSO scaling problems if branch_scale and chaining parameters mis-tuned (Armstrong 2020 Supp Methods) |
| Minigraph-Cactus (Hickey 2024 Nat Biotech 42:663) | minigraph SV graph base + Cactus base-level resolution | Pangenome graph (VCF, GFA, BAM) | Combines structural variation (minigraph) with base-level alignment (Cactus); HPRC standard | Requires reference; designed for intra-species pangenomes (90 human haplotypes) |
| LASTZ + UCSC chain/net (Kent 2003 PNAS 100:11484) | LASTZ local alignments -> chains (collinear) -> nets (one-to-one selection) | Reference-anchored alignment (chain, net, MAF) | Standard for UCSC genome browser tracks; well-validated; pairwise interpretable | Reference-bias; quadratic scaling for N genomes; per-pair manual chain/net pipeline |
| MUMmer4 nucmer (Marçais 2018 PLoS Comp Biol 14:e1005944) | Maximal exact matches -> clustering -> alignment extension | Pairwise alignment (delta format) | Fast for closely related (>= 70% identity); robust | Slow at low identity; pairwise only |
| minimap2 (Li 2018 Bioinformatics 34:3094) | Minimizer-based seeding + base-level alignment | SAM / PAF pairwise | Very fast; preset `-x asm5/10/20` for divergence | Lower accuracy at < 70% identity vs LASTZ |
| AnchorWave (Song 2022 PNAS 119:e2113075119) | CDS/exon anchors -> wavefront alignment | Sequence-level synteny + SV | WGD-aware (ploidy parameter); plant-friendly | CDS-anchored; intergenic resolution limited |
| Progressive Mauve (Darling 2010 PLoS One 5:e11147) | Local collinear blocks (LCBs) -> progressive alignment | Multi-genome bacterial alignment | Bacterial-friendly; handles inversions natively | Slow for > 30 genomes; less common today |
| SibeliaZ (Minkin & Medvedev 2020 Nat Commun 11:6327) | de Bruijn graph + colinear blocks | Bacterial / closely related multi-genome alignment | Scales to 100+ genomes | Highly diverged genomes lose recall |
| HAL toolkit (Hickey 2013 Bioinformatics 29:1341) | Hierarchical Alignment Format | Index, slice, project, extract from HAL | Standard substrate for downstream Cactus tools (TOGA, halSynteny, halLiftover) | HAL file version pinning required |
| Winnowmap2 (Jain 2022 Nat Methods 19:705) | minimap2-derived with masked repetitive minimizers | SAM / PAF | Optimized for highly repetitive regions (centromeres, telomeres) | Niche; for repeat-rich pairwise alignment |
| BLASTZ (legacy) | LASTZ predecessor | Local alignments | Historical; superseded by LASTZ | Use LASTZ instead |
| Multiz (Blanchette 2004 Genome Res 14:708) | Reference-anchored chain-merging | Multi-species MAF | UCSC's older multi-species pipeline | Superseded by Cactus |
| Lagan / mLagan (Brudno 2003) | Local-global hybrid alignment | Pairwise / multi-species | Historical; not actively developed | Use Cactus or LASTZ chain/net |

Methodology evolves; Cactus / Minigraph-Cactus / HAL are now the dominant production pipelines for vertebrate / mammal / plant WGA. LASTZ chains/nets remain the standard for adding a single new species to a reference-anchored ecosystem (UCSC, Ensembl).

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Vertebrate clade WGA, 10-500 genomes | Progressive Cactus | Reference-free; phylogenetic-tree guided; HAL output |
| Mammalian Zoonomia-scale WGA, 500-2000 genomes | Progressive Cactus with `--branchScale` tuning + HPC | Cactus scales to 1000s with proper Toil config |
| Plant clade WGA, with WGD | AnchorWave proali OR Cactus | AnchorWave WGD-aware; Cactus more general |
| Bacterial / archaeal genome alignment | progressiveMauve OR SibeliaZ | Designed for compact, rearrangement-rich genomes |
| Single new species added to reference (e.g. UCSC, Ensembl track) | LASTZ + chain/net | Standard reference-anchored pipeline; interpretable |
| Closely related strains / haplotypes (~95% identity) | minimap2 `-x asm5` or MUMmer nucmer | Fast pairwise; sufficient accuracy |
| Cross-species ~70-90% identity | LASTZ with `--strategy` tuned | minimap2 loses accuracy below ~70% |
| Cross-species ~50% identity (e.g. human vs zebrafish) | LASTZ with HoxD55 matrix | Specialized parameters for distant genomes |
| Pangenome graph for variant calling | Minigraph-Cactus or PGGB | See [[pangenome-analysis]] for graph-based variant calling |
| Repeat-rich genome (large mammal, plant) | Cactus with masked input | Pre-mask repeats with RepeatMasker/RepeatModeler2 |
| Centromeric / telomeric alignment | Winnowmap2 | Optimized for repetitive minimizers |
| Comparative annotation projection | Cactus -> TOGA + CESAR (see [[comparative-annotation-projection]]) | Cactus HAL is TOGA input |
| Conservation analysis (phyloP, PhastCons) | Cactus -> hal2maf to MAF for phyloP | MAF format from HAL is the substrate |
| Pairwise alignment for synteny | minimap2 `-x asm5` -> SyRI | See [[synteny-analysis]] |
| Long-read assembly pairwise validation | nucmer + dnadiff | Standard for assembly QC |
| Reference-to-pangenome lift-over | halLiftover or vg paths | HAL coordinate-system tools |

## Per-Tool Failure Modes

### Cactus running quadratically without species tree

**Trigger:** Running `cactus` without `--guide-tree` or with a poorly-resolved guide tree.

**Mechanism:** Cactus is designed as a progressive aligner: it aligns siblings, then their ancestor, recursively up the species tree. Without a guide tree, it tries all-vs-all alignment, which is O(N^2) and explodes for > 10 genomes.

**Symptom:** Cactus job time scales as N^2; cluster runs hit wall time before completing.

**Fix:** Always provide a guide tree in `seqFile.txt`:
```
(((human:0.05, chimp:0.05):0.05, macaque:0.1):0.2, mouse:0.5);
human    /path/to/human.fa
chimp    /path/to/chimp.fa
macaque  /path/to/macaque.fa
mouse    /path/to/mouse.fa
```
Branch lengths in substitutions per site; need not be precise but must reflect relative divergence. For poorly-resolved trees, use STAR (concatenation) tree from a few hundred concatenated single-copy orthologs (cf. [[ortholog-inference]]).

### Branch-scale parameter mis-tuned

**Trigger:** Cactus alignment producing very few or very many alignment columns relative to expected.

**Mechanism:** Cactus uses `--branchScale` to scale internal-branch lengths for LASTZ chaining parameters. Default scaling assumes vertebrate-like substitution rates; plants, bacteria, and viruses need adjustment. Wrong scaling produces under-aligned (too low) or over-merged (too high) blocks.

**Symptom:** HAL file shows alignment column count vastly different from expected; halStats reports unusual coverage; comparative annotation downstream fails.

**Fix:** Adjust `--branchScale` per clade: vertebrates default 1.0; plants 0.5-0.7; bacteria 0.3-0.5; viruses 0.1-0.3. Validate against known orthologs (BUSCO single-copies should align across all genomes).

### Toil checkpoint restart failures

**Trigger:** Cactus job interrupted (cluster timeout, OOM, node failure) and restart.

**Mechanism:** Cactus uses Toil's job store for checkpointing. Restart requires identical `jobStore` path, identical config, and pointing to the same machine type (some cloud-stored jobs require AWS / Google Cloud credentials).

**Symptom:** Cactus restart errors with "job store inconsistent" or "missing intermediate file."

**Fix:** Use `cactus --restart jobStore seqFile.txt output.hal`. For HPC, ensure the jobStore path is on shared storage accessible to all nodes. For cloud, use S3 or Google Cloud Storage with proper credentials. If restart fails, restart from scratch with a new jobStore.

### Reference bias in LASTZ chain/net

**Trigger:** Reporting "human-chimp synteny" from LASTZ chains/nets with human as reference.

**Mechanism:** Chain/net is reference-anchored; the alignment is biased toward the reference's genomic context. Repeat-content differences and assembly-quality differences between reference and query produce alignment artifacts.

**Symptom:** Different chains/nets from different reference choices show inconsistent SV calls; e.g. mouse-as-reference and human-as-reference produce different INV counts.

**Fix:** For multi-genome questions, use Cactus (reference-free). For pairwise questions, run chain/net both directions (A-on-B and B-on-A) and report consistent calls. Document reference choice in methods. SyRI on minimap2 alignments has similar issues; cross-validate.

### minimap2 sensitivity loss at < 70% identity

**Trigger:** Running `minimap2 -x asm10` or `-x asm20` between distantly related species.

**Mechanism:** Minimap2's minimizer-based seeding loses sensitivity at lower identity; many regions never anchor and produce gaps in alignment.

**Symptom:** Alignment coverage low (< 50% of query mapped); SyRI / SV-calling on this alignment reports few or no SVs; high-identity regions align but mid-identity regions don't.

**Fix:** Switch to LASTZ for < 70% identity; LASTZ uses 6-bp seeds with scoring matrices designed for distant comparison (HoxD55 for human-fish). For minimap2, try `-x asm20` with manual sensitivity tuning (`-k 19 -w 19`) but expect inferior performance vs LASTZ. minimap2 documentation recommends presets up to 20% divergence.

### Cactus producing wrong topology in ancestor assembly

**Trigger:** Cactus is asked to assemble an ancestral sequence at an internal node, but the input tree has wrong topology.

**Mechanism:** Cactus uses the guide tree to determine which descendants align to which ancestor. Wrong topology produces incorrect ancestor sequence at internal nodes.

**Symptom:** Ancestor sequences inconsistent with biological expectation; downstream phylogenetic analyses on Cactus ancestor sequences disagree with concatenation phylogenies.

**Fix:** Validate species tree before Cactus alignment using `phylogenetics/modern-tree-inference` on concatenated single-copy orthologs. If topology is uncertain, use a more conservative resolution (e.g. polytomy resolved via STRIDE).

### Repeats inflating LASTZ chain count

**Trigger:** Running LASTZ on unmasked or weakly masked genomes.

**Mechanism:** Unmasked TEs produce millions of paralogous LASTZ alignments; chains form from these, creating false-positive synteny / alignment.

**Symptom:** Chain count >> 100,000 per pair; many short chains in TE-rich regions; chain/net pipeline produces excessive chains where biology predicts few.

**Fix:** Softmask both genomes (RepeatModeler2 + RepeatMasker) before LASTZ. LASTZ honors softmasked sequences in its seeding step (lowercase = ignored at seeding, full alignment after extension). Hardmasking (N's) is too aggressive and breaks LASTZ scoring.

### HAL version mismatch breaking downstream tools

**Trigger:** Using HAL files generated with one version of Cactus / HAL toolkit with downstream tools using a different version.

**Mechanism:** HAL file format is versioned; minor version mismatch usually OK, major version mismatch causes silent or loud failures.

**Symptom:** halStats / TOGA / halLiftover errors about "unrecognized HAL version."

**Fix:** Use the HAL toolkit version that wrote the file (check via `halStats --filter halVersion` or git log of the run). For shared HAL files, pin Cactus version in documentation.

### Insufficient masking before Cactus

**Trigger:** Cactus on unmasked vertebrate genomes.

**Mechanism:** Cactus's CAF/BAR algorithms identify candidate alignment columns; without masking, TE-derived blocks dominate. Memory consumption explodes and alignment results are mostly TE pseudo-alignments.

**Symptom:** Cactus runs out of memory; HAL file is huge (> 100 GB for what should be 10 GB); halStats shows enormous column count.

**Fix:** Mask all genomes with RepeatMasker (or species-specific RepeatModeler2 library) before Cactus. Soft-masking is acceptable; Cactus handles lowercase appropriately. Document masking strategy.

### AnchorWave WGD parameter mis-set

**Trigger:** AnchorWave proali with wrong `--ploidy` value for the query genome.

**Mechanism:** AnchorWave's WGD-aware alignment matches each query CDS to ploidy-many reference CDS regions; wrong ploidy causes either undermatching (low ploidy) or overmatching (high ploidy).

**Symptom:** Reported alignment count inconsistent with biology (e.g. hexaploid wheat reports diploid-level matches).

**Fix:** Set `--ploidy` to query genome's ploidy (1 for haploid reference, 2 for diploid query when reference is diploid, 4 for tetraploid query vs diploid reference). Check with `python -m jcvi.compara.synteny depth` on output anchors.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| Cactus runtime scaling | Linear with tree-balanced clade; quadratic without tree | Armstrong 2020 Supp Methods |
| Cactus memory per genome (vertebrate, masked) | 8-32 GB | Armstrong 2020 |
| Cactus minimum genome assembly N50 | >= 1 Mb for reliable WGA; chromosome-level preferred | Standard convention |
| Branch scale for vertebrates | 1.0 default | Cactus docs |
| Branch scale for plants | 0.5-0.7 | Tuning convention |
| Branch scale for bacteria | 0.3-0.5 | Tuning convention |
| LASTZ identity range | 60-99% practical; 50% with HoxD55 matrix | LASTZ docs |
| minimap2 -x asm5 | <= 5% divergence | minimap2 docs |
| minimap2 -x asm10 | <= 10% divergence | minimap2 docs |
| minimap2 -x asm20 | <= 20% divergence | minimap2 docs |
| MUMmer nucmer recommended identity | >= 70% | Marçais 2018 |
| Soft-masking required for WGA | >= 90% of TE families masked | Standard QC |
| HAL file size, vertebrate 100-genome | ~100-500 GB per HAL | Approximate; varies with masking |
| Cactus per-job CPU requirement | 4-16 cores per leaf | Cactus docs |
| Minigraph-Cactus haplotype count | tested at 90 human haplotypes | Hickey 2024 |
| LASTZ chain-net minimum chain length | 5000 bp default; tunable | UCSC convention |
| Toil retry limit | --retryCount 3 typical | Toil docs |
| Cactus dustmasking sensitivity | default sensitive; adjust per clade | Cactus configs |
| AnchorWave anchor minimum identity | 80% default; lower for distant species | AnchorWave docs |
| progressiveMauve maximum genomes (practical) | 20-30 | Darling 2010 |
| SibeliaZ minimum block size | 50 bp default; tunable | Minkin 2020 |
| HAL toolkit halLiftover precision | base-level for collinear; gap-aware for SV | Hickey 2013 |

## Progressive Cactus Standard Workflow

**Goal:** Build a multi-species reference-free WGA from N genomes.

**Approach:** Provide guide tree + per-species genome FASTA in `seqFile.txt` -> run `cactus` with Toil job store -> produce HAL -> extract MAF / synteny / annotation as needed.

```bash
# 1. Pre-mask genomes (softmask)
for fa in genomes/*.fa; do
    species=$(basename $fa .fa)
    BuildDatabase -name ${species}_DB $fa
    RepeatModeler -database ${species}_DB -threads 16
    RepeatMasker -lib ${species}_DB-families.fa -xsmall -pa 16 $fa
done

# 2. Prepare seqFile (tree + paths)
cat > seqFile.txt << 'EOF'
(((human:0.05, chimp:0.05):0.05, macaque:0.1):0.2, mouse:0.5);
human    genomes/human.fa.masked
chimp    genomes/chimp.fa.masked
macaque  genomes/macaque.fa.masked
mouse    genomes/mouse.fa.masked
EOF

# 3. Run Cactus
cactus jobStore_path seqFile.txt output.hal \
    --binariesMode local \
    --workDir /tmp/cactus_work \
    --maxCores 64 \
    --logFile cactus.log

# 4. Inspect HAL
halStats output.hal
# Reports: genome count, chromosome / contig count per genome, alignment column count

# 5. Extract MAF for downstream phyloP / conservation
hal2maf output.hal reference_genome ucsc.maf --refGenome human --chunkSize 1000000
```

For HPC deployment, use Toil's `--batchSystem slurm` or `--batchSystem kubernetes` for cluster-scale runs.

## LASTZ + UCSC Chain/Net Standard Pipeline

**Goal:** Pairwise WGA of a query genome to a reference, producing UCSC-format chain and net files.

**Approach:** LASTZ pairwise -> axtChain -> chainNet -> netSyntenic.

```bash
# 1. LASTZ pairwise alignment per chromosome
for chr in $(cut -f1 reference.fa.fai); do
    lastz reference.fa[multiple,unmask] query.fa[multiple,unmask] \
        --hspthresh=3000 --ydrop=9400 --gappedthresh=3000 \
        --ambiguous=iupac --inner=2000 --score=$LASTZ_PARAMS \
        --format=axt \
        > pairwise/${chr}.axt
done

# 2. Build chains. `axtChain` reads PSL when given `-psl`; for AXT input from lastz
# (`--format=axt`) drop the `-psl` flag, or emit PSL from LASTZ with `--format=psl` upstream.
for axt in pairwise/*.axt; do
    chr=$(basename $axt .axt)
    axtChain -linearGap=loose $axt reference.2bit query.2bit chains/${chr}.chain
done

# 3. Merge chains
chainMergeSort chains/*.chain > merged.chain
chainPreNet merged.chain reference.size query.size pre.chain

# 4. Build net (one-to-one selection)
chainNet pre.chain reference.size query.size \
    target.net query.net
netSyntenic target.net target.syntenic.net

# 5. Filter and convert
netToAxt target.syntenic.net pre.chain reference.2bit query.2bit final.axt
axtToMaf final.axt reference.size query.size final.maf
```

## Minigraph-Cactus for Pangenome Construction

**Goal:** Build a pangenome graph from related genomes (intra-species or sister-species).

**Approach:** Provide reference + alternative genomes -> minigraph for SV graph -> Cactus for base-level alignment.

```bash
# Prepare seqFile (similar to Cactus but designates reference)
cat > pangenome_seqs.txt << 'EOF'
reference   human.fa
hap1        haplotype1.fa
hap2        haplotype2.fa
hap3        haplotype3.fa
EOF

# Run minigraph-cactus
cactus-pangenome jobStore_path pangenome_seqs.txt \
    --outDir output_pangenome \
    --outName pangenome \
    --reference reference \
    --vcf \
    --gfa \
    --gbz \
    --indexCores 32 \
    --mapCores 32

# Outputs:
#   pangenome.gbz       Genome Browser Z compressed graph
#   pangenome.gfa.gz    Graph Fragment Assembly format
#   pangenome.vcf.gz    VCF for short variants from graph
#   pangenome.full.hal  Full Cactus HAL for downstream tools
```

For HPRC-scale (90 haplotypes), use `--mapCores 64 --indexCores 64` on a high-memory node. See [[pangenome-analysis]] for downstream graph analysis.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Cactus alignment column count >> LASTZ pairwise | Cactus includes ancestral alignment; pairwise is direct only | Cactus is expected to have more columns when integrating across siblings |
| minimap2 aligns 60%; LASTZ aligns 80% | Sensitivity differs at low identity | Use LASTZ for distant comparisons |
| LASTZ chain/net A-on-B differs from B-on-A | Reference bias | Use Cactus reference-free for unbiased multi-genome questions |
| Cactus HAL vs Multiz MAF disagree at deep node | Cactus reference-free; Multiz reference-anchored | Cactus is the modern standard; Multiz mostly superseded |
| AnchorWave finds collinearity in WGD region, LASTZ doesn't | AnchorWave WGD-aware; LASTZ unaware | AnchorWave correct for WGD lineages |
| Minigraph-Cactus VCF and Cactus pairwise differ | Pangenome graph integrates SV; pairwise reports relative to reference | Minigraph-Cactus is for variant-aware pangenome work |
| nucmer (MUMmer) and minimap2 disagree at 80% identity | nucmer more sensitive at higher identity; minimap2 faster | Either is fine; cross-validate |
| Cactus runs but HAL file empty for some genomes | Genome misnamed in seqFile or masking issue | Verify seqFile genome names match FASTA headers; verify masking |
| progressiveMauve / SibeliaZ vs Cactus on bacteria | Bacterial-specific aligners better-tuned for compact genomes | progressiveMauve / SibeliaZ preferred for bacterial work; Cactus is overkill |

**Operational rule for publication:** Multi-species WGA for comparative genomics uses Progressive Cactus with documented guide tree, masking, and Toil configuration. Pairwise reference alignments use LASTZ chain/net for UCSC tracks or minimap2 for fast pangenomics. Document choice of alignment tool, reference (if applicable), parameter values, and any non-default settings.

## Cohort Gotchas

- **Highly repetitive genomes (large mammals, plants):** memory consumption explodes; ensure >= 256 GB RAM per node
- **Centromeric / telomeric regions:** unalignable in most assemblies; Winnowmap2 for explicit telomere-to-telomere work; otherwise mask
- **Sex chromosomes:** different content; Cactus aligns where possible; downstream analyses must respect sex-chromosome specifics
- **Recently diverged inbred strains (e.g. mouse strains):** > 99% identity; minimap2 -x asm5 is fine; Cactus overkill
- **Pangenome with structural variants:** Minigraph-Cactus preferred over straight Cactus
- **Multi-region polyploids:** AnchorWave proali with ploidy specification; subgenome-aware

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Reference bias?" | Cactus is reference-free; or for pairwise pipelines, both directions reported |
| "Guide tree?" | Provided in seqFile; built from concatenated single-copy orthologs |
| "Masking?" | RepeatModeler2 species-specific + RepeatMasker softmask before WGA |
| "Tool choice?" | Cactus for multi-species (reference-free, scalable); LASTZ chain/net for single-pair (UCSC standard) |
| "Branch scale?" | Set to 1.0 for vertebrates; 0.5 for plants; documented |
| "Assembly quality?" | N50 reported; BUSCO completeness reported; chromosome-level preferred |
| "WGD?" | AnchorWave proali with explicit ploidy; or Cactus with subgenome assignment |
| "Reproducibility?" | Cactus version + Toil version + branch scale documented; jobStore path retained for restart |
| "HAL downstream?" | TOGA for annotation projection; hal2maf for phyloP / PhastCons; halSynteny for SV |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Cactus fails with "guide tree mismatch" | seqFile genome names don't match tree leaves | Verify exact name match in seqFile lines vs newick |
| Cactus memory exhaustion | Unmasked repeats | Pre-mask with RepeatMasker softmask |
| Toil restart fails | jobStore corruption | Restart from scratch with new jobStore |
| LASTZ produces 0 alignments | Sensitivity too strict | Lower `--hspthresh`; check seed-and-extend params |
| MAF from hal2maf huge | Default chunkSize too small | Increase `--chunkSize 10000000` for large genome |
| minimap2 SAM file empty | Wrong preset for divergence | Try `-x asm20` or LASTZ for distant comparison |
| HAL file unreadable in another tool | HAL version mismatch | Pin HAL toolkit version |
| Cactus produces many tiny alignment blocks | Branch scale too low | Increase `--branchScale` to 1.0 |
| nucmer maxmatch slow | Many MUMs (unmasked repeats) | Use `--mum` (unique) or mask repeats |
| AnchorWave reports few anchors | CDS GFF parsing issue | Verify GFF3 with `gene` and `CDS` features |
| TOGA fails with "no chain" | Cactus HAL not converted to chain | Use halSynteny + UCSC chain conversion |
| Cactus dead-locks on shared filesystem | Toil multi-node sync issue | Use AWS S3 / GCS jobStore for cloud |

## Tool Installation Notes

```bash
# Progressive Cactus (Docker / Singularity preferred for reproducibility)
docker pull quay.io/comparative-genomics-toolkit/cactus:latest
# Or via virtualenv
python3 -m venv cactus_env && source cactus_env/bin/activate && pip install cactus

# HAL toolkit (usually bundled with Cactus)
git clone https://github.com/ComparativeGenomicsToolkit/hal && cd hal && make

# LASTZ + UCSC kentUtils
conda install -c bioconda lastz ucsc-axt-chain ucsc-chain-merge-sort ucsc-chain-pre-net ucsc-chain-net ucsc-net-syntenic ucsc-net-to-axt

# MUMmer4
conda install -c bioconda mummer

# minimap2 + winnowmap2
conda install -c bioconda minimap2 winnowmap

# AnchorWave
conda install -c bioconda anchorwave

# progressiveMauve
conda install -c bioconda mauve

# SibeliaZ
conda install -c bioconda sibeliaz

# RepeatMasker / RepeatModeler2 (for masking)
conda install -c bioconda repeatmasker repeatmodeler

# Toil for HPC / cloud
pip install toil[all]
```

For Cactus on HPC, set up Toil with the appropriate batch system (`--batchSystem slurm`); for cloud, use Toil's S3/GCS support and a job store on shared cloud storage.

## References

- Armstrong J et al 2020 Nature 587:246 (Progressive Cactus)
- Hickey G et al 2024 Nat Biotech 42:663 (Minigraph-Cactus)
- Hickey G et al 2013 Bioinformatics 29:1341 (HAL toolkit)
- Kent WJ et al 2003 PNAS 100:11484 (UCSC chain/net)
- Schwartz S et al 2003 Genome Res 13:103 (LASTZ)
- Marçais G et al 2018 PLoS Comp Biol 14:e1005944 (MUMmer4)
- Li H 2018 Bioinformatics 34:3094 (minimap2)
- Li H 2021 Bioinformatics 37:4572 (minimap2 long-read updates)
- Song B et al 2022 PNAS 119:e2113075119 (AnchorWave)
- Darling AE et al 2010 PLoS One 5:e11147 (progressiveMauve)
- Minkin I & Medvedev P 2020 Nat Commun 11:6327 (SibeliaZ)
- Jain C et al 2022 Nat Methods 19:705 (Winnowmap2)
- Blanchette M et al 2004 Genome Res 14:708 (Multiz)
- Brudno M et al 2003 GR 13:721 (LAGAN)
- Liao W-W et al 2023 Nature 617:312 (HPRC draft pangenome)
- Sirén J et al 2024 Nat Methods 21:2017 (vg personalized pangenome references). The vg-Giraffe pangenome mapper is Sirén J et al 2021 Science 374:abg8871.
- Garrison E et al 2024 Nat Methods (PGGB)
- Paten B et al 2011 Genome Res 21:1512 (Cactus algorithm)
- ComparativeGenomicsToolkit Cactus + HAL toolkit documentation (https://github.com/ComparativeGenomicsToolkit)

## Related Skills

- comparative-genomics/synteny-analysis - Synteny detection from WGA (Cactus -> halSynteny)
- comparative-genomics/comparative-annotation-projection - TOGA + CESAR uses Cactus HAL
- comparative-genomics/pangenome-analysis - Minigraph-Cactus / PGGB for pangenome graph construction
- comparative-genomics/whole-genome-duplication - Ks-dating uses WGA-derived gene pair alignments
- alignment/multiple-alignment - MAFFT / MUSCLE for protein MSA from WGA-derived orthologs
- alignment/pairwise-alignment - LASTZ / minimap2 / MUMmer for sequence-level comparison
- alignment/structural-alignment - WGA-extracted protein orthologs for downstream alignment
- genome-assembly/assembly-qc - BUSCO / Compleasm checks before WGA
- variant-calling/structural-variant-calling - WGA-based SV detection from genome assemblies
- causal-genomics/heritability-partitioning - LDSC partitioning using phyloP scores from Cactus HAL
