---
name: bio-comparative-genomics-comparative-annotation-projection
description: Project gene annotations across genomes using TOGA (Kirilenko 2023 whole-genome-alignment chain-based projection with intactness classification), CESAR 2.0 (Sharma & Hiller 2017 codon-aware exon projection), LiftOff (Shumate & Salzberg 2020 reference-based annotation transfer), Liftover (UCSC), GeMoMa (Keilwagen 2019 evidence-based projection), and Comparative Annotation Toolkit (CAT). Use when transferring annotations from a well-annotated reference to query genome(s), classifying gene-loss vs gene-intact across many genomes at scale, building Zoonomia-style comparative annotations across hundreds of mammals or birds (Kirilenko 2023), detecting pseudogenization, projecting alternative isoforms, or selecting between WGA-anchored (TOGA) vs ortholog-based (LiftOff) annotation transfer strategies.
tool_type: cli
primary_tool: TOGA
---

## Version Compatibility

Reference examples tested with: TOGA 1.1.7+ (hillerlab/TOGA; Kirilenko 2023 Science 380:eabn3107), CESAR 2.0 (Sharma & Hiller 2017 NAR 45:8369), LiftOff 1.6.3+ (Shumate & Salzberg 2021 Bioinformatics 37(12):1639-1643), Comparative Annotation Toolkit (CAT) 2.4+, GeMoMa 1.9+ (Keilwagen 2019 Methods Mol Biol 1962:161), UCSC liftOver 2024+, Cactus 2.9.1+ (for HAL input), HAL toolkit 2.3+, NextFlow 24+ for TOGA pipeline, BUSCO 5.7+ / Compleasm 0.2.7+ for QC, Snakemake 8.0+ for CAT, R 4.4+. The current TOGA expects HAL from Cactus 2.5+; older HAL formats may fail.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `toga.py --help`; `cesar --help`; `liftoff --version`
- Python: `pip show liftoff`
- Java: `gemoma --help` (Java 11+)

If code throws `TOGA chain file missing`, `CESAR fragment not found`, `LiftOff annotation not parsed`, the toolchain expects specific input formats: TOGA needs HAL or chain/net files from Cactus / LASTZ; CESAR needs exon-level GFF; LiftOff needs reference GFF and aligned FASTA. Pre-process with the appropriate format conversion.

# Comparative Annotation Projection

**"Annotate this new genome using my well-annotated reference"** -> Annotation projection from a reference is the modern alternative to de novo gene prediction; it produces high-quality, comparable annotations across genomes by leveraging evolutionary conservation. The 2023-era standard is **TOGA + CESAR 2.0** (Kirilenko 2023 Science 380:eabn3107), which uses whole-genome alignment chains + ML classification + codon-aware exon projection to scale to hundreds of genomes (Zoonomia: 488 mammals; Bird10000 Genomes: 501 birds). For ortholog-based projection (no WGA needed), **LiftOff** (Shumate & Salzberg 2021 Bioinformatics 37(12):1639) is the standard. The critical decision is **WGA-anchored (TOGA) vs ortholog-anchored (LiftOff)**: TOGA explicitly classifies gene intactness vs loss using the alignment chains, LiftOff relies on reciprocal-best-hit equivalents.

- CLI: `toga.py --chain chain.bb --bed ref.bed --tDB target.2bit --qDB query.2bit --pn project_name` -- WGA-based projection
- CLI: `cesar -i exons.fa -d 4 -o output.aln` -- codon-aware exon alignment (used internally by TOGA)
- CLI: `liftoff -g ref.gff query.fa ref.fa -o query.gff` -- ortholog-based projection
- CLI: `gemoma` -- evidence-based comparative annotation (Java)
- CLI: `cat` (Comparative Annotation Toolkit) -- multi-species annotation projection

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| TOGA (Kirilenko 2023 Science 380:eabn3107) | Cactus HAL or LASTZ chains -> ML projection + intactness classifier | Per-gene I/PI/UL/L/M/PM codes; orthology classification; coding annotation via CESAR 2.0 | Modern paradigm; explicit gene-loss detection at scale; Zoonomia / Bird10000 standard | Requires Cactus WGA; not for prokaryotes |
| CESAR 2.0 (Sharma & Hiller 2017 NAR 45:8369) | HMM-based codon-aware exon projection | Aligned exons + frame preservation | Most accurate exon projection from WGA; preserves frame across indels | Used internally by TOGA; standalone use more rare |
| LiftOff (Shumate & Salzberg 2021 Bioinformatics 37(12):1639-1643) | Read-mapping-style ortholog detection + GFF transfer | Lifted GFF | Fast; no WGA required; standard for query-vs-reference pairs | Tandem duplicates ambiguous; not for gene loss detection |
| UCSC liftOver | Coordinate-based lift using chain files | Coordinate-lifted regions | Standard for coordinate transfers; not for gene annotations | Doesn't handle gene structure changes |
| Comparative Annotation Toolkit (CAT) | Snakemake workflow integrating Augustus + LiftOff + TransMap | Per-species comparative annotation | Integrates de novo + projection | Snakemake setup complex |
| GeMoMa (Keilwagen 2019 Methods Mol Biol 1962:161) | Reference protein homology + evidence integration | Comparative gene annotation | Combines multiple reference species evidence | Slower; less popular than TOGA / LiftOff |
| AUGUSTUS (Stanke 2008) | De novo prediction; not strictly projection | Per-genome annotation | Augments projection with de novo | Standalone de novo; lower comparative accuracy |
| BRAKER3 (Brůna 2024) | Augustus + GeneMark-ETP + RNA-Seq + protein | Comparative-aware de novo | Modern de novo with evidence | Not strictly projection |
| Funannotate (palmer lab) | Multi-evidence annotation including LiftOff | Funannotate annotations | Integrates evidence | Setup complex |
| Comparative Annotation Pipeline (CAP) | Earlier WGA-based annotation | Per-species per-gene | Historical; replaced by TOGA | Use TOGA |
| TransMap (Diekhans 2007) | UCSC genome browser annotation lifter | Per-locus lift | Tool for UCSC tracks | Tool-specific |
| Maker (Cantarel 2008) | Evidence-based de novo + projection | Per-genome annotation | Combines evidence | Maker is for novel genomes; LiftOff for transfer |

Methodology evolves; the Kirilenko 2023 TOGA paradigm (WGA-anchored + intactness classification) is the gold standard for vertebrate-scale comparative annotation. For pairwise transfers, LiftOff is the modern standard. Verify the current TOGA documentation (hillerlab/TOGA) before locking on a single approach.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| Annotate hundreds of mammal / bird genomes | TOGA with Cactus HAL | Scales to Zoonomia / Bird10000 |
| Annotate single new genome from reference | LiftOff | Fast; no WGA required |
| Detect gene loss across mammals | TOGA intactness classification | Explicit I/PI/UL/L/M/PM codes |
| Project alternative isoforms | TOGA (preserves multiple transcripts) | Standard |
| Project annotations to assembly with high N50 + chromosome-level | TOGA | Requires good assembly |
| Project to fragmented draft assembly | LiftOff (more tolerant) | LiftOff works on draft assemblies |
| Multi-species annotation pipeline | CAT (Snakemake) | Integrated workflow |
| Annotate plant genome from Arabidopsis | LiftOff with plant-specific options | Standard for plant work |
| Pseudogenization detection at scale | TOGA + intactness analysis | Designed for this |
| Reference-free gene prediction | BRAKER3 or AUGUSTUS | De novo; not projection |
| Comparative annotation of multiple references | GeMoMa | Multi-reference evidence integration |
| UCSC genome browser coordinate transfer | liftOver tool | Coordinate-specific |
| Annotation transfer to closely related strain (>95% ANI) | LiftOff | High accuracy at close divergence |
| Annotation transfer to deep divergence (mammal to fish) | TOGA + manual review | Requires WGA; expect lower coverage |
| Project annotations with WGD-aware handling | AnchorWave + TOGA-like or custom workflow | WGD-aware tools |
| Annotate non-coding RNAs | Specialized tools (Rfam, ncRNA-specific) | RNA detection different problem |
| Annotate immune / repetitive genes (MHC, OR) | PGR-TK MAP graph or manual | Repetitive regions; use [[pangenome-analysis]] |
| Annotate transposable elements | RepeatMasker / RepeatModeler | TE annotation different problem |
| Validate projected annotations | RNA-Seq alignment to projected | RNA-Seq evidence is gold standard |

## Per-Tool Failure Modes

### TOGA chain file missing or incompatible

**Trigger:** Running TOGA on Cactus HAL without proper chain file extraction.

**Mechanism:** TOGA requires UCSC-style chain files derived from Cactus HAL or LASTZ chains/nets pipeline. Cactus HAL doesn't directly produce chain files; conversion via halSynteny + chainNet + axtChain is required.

**Symptom:** TOGA fails with "chain file not found" or "no syntenic blocks for query."

**Fix:** Use `halSynteny` (HAL toolkit) to extract syntenic blocks; convert to chain format via `axtChain` and `chainNet`. The TOGA Nextflow wrapper handles this automatically; for manual runs, see UCSC kentUtils chain documentation.

### CESAR exon-fragment misalignment in highly divergent species

**Trigger:** Projecting from mouse to fish (~400 Myr divergence); many exons fail CESAR projection.

**Mechanism:** CESAR's HMM model is calibrated for vertebrate divergence (< 100 Myr typical). At deep divergence, exon boundaries shift; CESAR may misalign or fail to project.

**Symptom:** Many genes in mouse have TOGA "M" (missing) or "PI" (partial-intact) classification in fish; coverage of expected genes is low.

**Fix:** TOGA documentation recommends < 300 Myr divergence for reliable projection. For deeper divergence, manual review of failed exons; consider GeMoMa with multiple reference species. Some genes won't project because they're truly absent (orphan genes); others fail due to alignment limitations.

### LiftOff tandem duplicate ambiguity

**Trigger:** LiftOff on genomes with extensive tandem duplications (e.g., NLR clusters in plants, olfactory receptors in mammals).

**Mechanism:** LiftOff uses ortholog detection similar to OrthoFinder; tandem duplicates create many similar sequences, making reciprocal-best-hit identification ambiguous.

**Symptom:** LiftOff reports many "multimapped" genes; tandem clusters have one-to-many or many-to-one orthology calls.

**Fix:** Pre-collapse tandem duplicates manually; or use LiftOff with `-mismatch 5` and `-flank 0.5` for more relaxed mapping; or use TOGA which has tandem-aware classification.

### TOGA intactness classification false negatives

**Trigger:** TOGA classifies a gene as "Lost" when it is actually intact.

**Mechanism:** TOGA uses ML classifier on chain features + frame preservation; assembly gaps, short alignment fractions, or CESAR projection failures can cause false-loss calls.

**Symptom:** TOGA "Lost" gene actually present in independent validation (RNA-Seq, manual inspection); known biology contradicts loss.

**Fix:** Manual review of TOGA "Lost" calls in the loss_summ_data.tsv; cross-validate with RNA-Seq mapping; use ID + biology to verify. TOGA's classifier is calibrated for mammals/birds; non-canonical genome architectures may produce false losses.

### Reference choice bias

**Trigger:** Projecting from one reference (e.g., mouse) produces different annotation than from another (e.g., human).

**Mechanism:** Each reference's annotation has its own biases (gene structures, splice variants, missing genes). Projection inherits these biases; different references produce somewhat different annotations.

**Symptom:** Mouse-reference TOGA annotation has gene X missing in query; human-reference TOGA has it present; or different exon structures.

**Fix:** Project from multiple references; consensus annotation. CAT integrates multi-reference projection; manual review of inconsistencies. Document reference choice impact.

### Pseudogenization vs gene loss distinction

**Trigger:** Reporting a "lost" gene that retains coding sequence (frame may be conserved but expression lost).

**Mechanism:** TOGA detects loss of coding capacity (intactness), but doesn't directly identify pseudogenization (loss of expression). A pseudogene with intact reading frame may be classified "Intact" by TOGA.

**Symptom:** TOGA "Intact" annotation but the gene is pseudogene per RNA-Seq + Ribo-Seq evidence.

**Fix:** Combine TOGA intactness with expression data (RNA-Seq from species of interest); apply PseudoPipe (Zhang 2006 GR 16:1041) or RetroFinder (Baertsch 2008 GR 18:1675) for systematic pseudogene detection.

### Splice variant inconsistency across projections

**Trigger:** Projecting genes with extensive alternative splicing.

**Mechanism:** Reference may have alternative splice variants; projection of alternative isoforms requires per-isoform alignment which may fail for some variants.

**Symptom:** Query species has fewer projected isoforms than reference; canonical isoform present but alternatives missing.

**Fix:** TOGA projects each transcript independently; manual review of dropped isoforms. RNA-Seq from species of interest for novel splice variants.

### Annotation pipeline reference quality affecting projection

**Trigger:** Projecting from an outdated or buggy reference annotation.

**Mechanism:** Errors in reference annotation propagate to all projections. A wrong exon boundary in mouse is propagated to all mammals.

**Symptom:** Same exon boundary error appears across many projected annotations.

**Fix:** Verify reference annotation quality via BUSCO + manual gene model review; use updated Ensembl / NCBI releases (current 2024-Q4).

### Polyploid query genome handling

**Trigger:** Projecting annotations onto a polyploid query without subgenome consideration.

**Mechanism:** Projection sees multiple homeologous regions; ortholog detection ambiguous between subgenomes.

**Symptom:** Polyploid query annotation has ~2x the genes expected; many "redundant" projections from homeologs.

**Fix:** Assign subgenomes before projection (see [[whole-genome-duplication]]); project to each subgenome separately. AnchorWave proali handles ploidy; LiftOff doesn't natively.

### Chromosome-level vs scaffold-level reference

**Trigger:** Projecting from chromosome-level reference to scaffold-level query.

**Mechanism:** Scaffold-level query has gaps and ambiguous gene assignments; projection mostly succeeds but some genes split across scaffolds.

**Symptom:** Projected GFF has fragmented gene models; some genes have 2-3 entries across scaffolds.

**Fix:** Pre-scaffold query (Hi-C scaffolding if possible); or accept partial annotations; document fragmentation rate. TOGA reports "Partial Intact" for these cases.

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| TOGA intactness classes (loss_summ_data.tsv) | I (intact), PI (partial intact), UL (uncertain loss), L (lost), M (missing/assembly gap), PM (partial missing) | Kirilenko 2023 + TOGA repo |
| TOGA orthology relationships (orthology_classification.tsv) | one2one, one2many, many2one, many2many, PG (paralogous projection / no orthologous chain) | Kirilenko 2023 |
| TOGA "Intact" classification confidence | ML classifier posterior > 0.9 | Kirilenko 2023 supp |
| LiftOff coverage | >=80% of reference gene length aligned | Default |
| LiftOff identity | >=70% nucleotide identity (default) | Default |
| Maximum divergence for TOGA | ~300 Myr (vertebrate); validate per clade | Kirilenko 2023 |
| Maximum divergence for LiftOff | ~80% nucleotide identity | Empirical |
| Maximum divergence for CESAR | ~150 Myr (vertebrate) | Sharma & Hiller 2017 |
| Reference annotation BUSCO completeness | >= 95% | Standard QC |
| Assembly N50 for projection | >= 1 Mb; chromosome-level preferred | Standard |
| Annotation transfer success rate | 90-95% for closely related (< 50 Myr); 60-80% for moderately diverged | Empirical |
| Multi-reference consensus | >= 2 references agreeing | Manual standard |
| Pseudogene classification threshold | TOGA "Lost" + no RNA-Seq evidence | Operational |
| Tandem cluster window for LiftOff | 50 kb default | Default |
| GeMoMa minimum protein identity | 60% | Default |
| CAT pipeline runtime per genome | 1-5 hours on 16 cores | Empirical |
| TOGA per-genome runtime | 30 min - 5 hr on 16 cores | Empirical |
| Nextflow scaling | scales with cores | Standard |
| Reference annotation version | Ensembl / NCBI release 2024-Q4 minimum | Standard |
| Splice variant count | report; per-transcript projection | Variable |

## TOGA Standard Workflow

**Goal:** Project annotations from reference to query genome(s), classifying gene-loss / intactness.

**Approach:** Cactus WGA -> halSynteny + chainNet -> TOGA Nextflow pipeline.

```bash
# Prerequisites: Cactus HAL file from [[whole-genome-alignment]]
# OR LASTZ chain/net pipeline output

# 1. Extract syntenic blocks from HAL
halSynteny output.hal reference query --queryGenome query > query.synteny.psl

# 2. Convert PSL to UCSC chain format
axtChain -psl -linearGap=loose query.synteny.psl reference.2bit query.2bit chains/query.chain.gz
# Note: `-psl` is correct here only if the input is PSL. When the input comes from
# `lastz --format=axt`, drop the `-psl` flag (or emit PSL from LASTZ first).

# 3. Run TOGA. The canonical invocation is `python toga.py` from the TOGA checkout;
# the Nextflow-style command shown below mirrors the same arguments but may not be
# the standard entry point in your release -- verify against the hillerlab/TOGA README.
python toga.py \
    chains/query.chain.gz \
    reference_annotation.bed \
    reference.2bit \
    query.2bit \
    --pn project_name \
    --cpus 32

# Output:
#   project_name/loss_summ_data.tsv           Per-gene intactness call
#   project_name/orthology_classification.tsv One-to-one / one-to-many / many-to-many
#   project_name/query_annotation.bed         Lifted gene annotation
#   project_name/query_annotation.gff         GFF format
#   project_name/cesar_alignment/             CESAR exon alignments
```

```python
'''Parse TOGA loss summary to identify gene loss vs intact.'''
import pandas as pd


def load_toga_loss(loss_summary_path):
    '''loss_summary_data.tsv columns: TRANSCRIPT, STATUS, IS_INTACT, ...'''
    df = pd.read_csv(loss_summary_path, sep='\t')
    return df


def classify_genes(df):
    '''Standard TOGA classification: I/PI/UL/L/M/PM'''
    return df.groupby('STATUS')['TRANSCRIPT'].count()


def filter_high_confidence_intact(df):
    '''Filter to high-confidence intact (I) only. PI is partial-intact; L is Lost.'''
    return df[df['STATUS'] == 'I']
```

## LiftOff for Pairwise Annotation Transfer

**Goal:** Transfer reference annotation to query genome via ortholog mapping.

**Approach:** Minimap2-based ortholog detection -> per-gene transfer.

```bash
# Standard LiftOff (verify flags with `liftoff --help`)
liftoff -g reference.gff \
    query.fa reference.fa \
    -o query.gff \
    -u unmapped.txt \
    -copies \
    -overlap 0.5 \
    -mismatch 2 \
    -gap 5 \
    -threads 16

# Output:
#   query.gff             Lifted annotations
#   unmapped.txt          Genes failed to lift
```

For closely related species, default settings suffice. For divergent (75-90% identity), use `-mismatch 5 -gap 10`. For tandem-rich regions, `-copies` allows multiple projections.

## Comparative Annotation Toolkit (CAT)

```bash
# Setup
git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit
cd Comparative-Annotation-Toolkit && pip install .

# Edit config.yaml with reference + query genomes
# Then run Snakemake pipeline
snakemake --use-conda --cores 32 --configfile config.yaml
```

CAT integrates AUGUSTUS + LiftOff + TransMap; output is multi-species comparative annotation.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| TOGA "Lost" vs LiftOff "Mapped" | LiftOff more permissive; doesn't check frame preservation | Trust TOGA for loss claims; LiftOff for pairwise transfer |
| TOGA "Partial Intact" vs LiftOff "Mapped" | Frame disruption | Cross-validate; consider biology |
| Mouse-reference annotation vs human-reference | Reference choice bias | Multi-reference consensus |
| TOGA "Intact" but no RNA-Seq evidence | Possible pseudogene with intact ORF | Add expression evidence; reclassify |
| CESAR fails exon | Deep divergence or assembly gap | Manual review; consider GeMoMa alternative |
| GeMoMa vs LiftOff disagree | Evidence integration vs ortholog-based | GeMoMa for evidence-rich; LiftOff for fast pairwise |
| TOGA chain doesn't cover gene | Cactus alignment quality issue | Re-run Cactus with adjusted parameters; verify assembly |
| LiftOff produces duplicate transfers | Tandem cluster | Use `-copies` or restrict; manual review |
| CAT pipeline integrates de novo + projection | Combined evidence | Trust integrated; better than single tool |
| BUSCO of projected annotation low | Missing core genes | Likely tool failure; re-run with relaxed parameters |

**Operational rule for publication:** TOGA + Cactus HAL for clade-level annotation (Zoonomia-style); LiftOff for pairwise transfer to closely related (<100 Myr); BRAKER3 / Funannotate for de novo where projection fails; manual review of TOGA "Lost" / "PI" calls.

## Cohort Gotchas

- **Plant comparative annotation:** GENESPACE handles synteny-aware ortholog detection; project via plant-specific tools
- **Bacterial annotation:** different problem; use [[pangenome-analysis]] with Bakta consistent annotation
- **Single-cell expression data:** RNA-Seq evidence for projected genes essential
- **Repetitive genes (MHC, OR):** projection unreliable; use [[pangenome-analysis]] with PGR-TK
- **Recently diverged strains:** LiftOff with strict parameters; high accuracy
- **Polyploid query:** assign subgenomes first ([[whole-genome-duplication]])
- **Distantly related to reference (>300 Myr):** projection rate drops; consider de novo with comparative evidence
- **Fragmented draft genomes:** projection works on contigs but gene splits possible
- **Reference annotation quality:** verify BUSCO before propagating to all projections
- **Non-canonical genomes (B chromosomes, supernumerary):** typically excluded from projection

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Reference annotation quality?" | BUSCO completeness reported; updated to current Ensembl / NCBI release |
| "Maximum divergence for TOGA?" | < 300 Myr for vertebrate; documented and respected |
| "Tandem duplicate handling?" | LiftOff `-copies`; or pre-collapse for TOGA |
| "Gene loss detection?" | TOGA intactness classification (I/PI/UL/L/M/PM); cross-validated with RNA-Seq |
| "Pseudogenization?" | TOGA "Intact" classification verified against RNA-Seq + Ribo-Seq |
| "Cactus WGA quality?" | Pre-filtered repeats; BUSCO on reference and query; Toil reproducibility |
| "Multi-reference?" | Consensus annotation from 2-3 reference species; documented disagreements |
| "Polyploid?" | Subgenomes assigned via [[whole-genome-duplication]]; per-subgenome projection |
| "Annotation transfer rate?" | Per-species coverage reported (e.g. 85% projection success) |
| "Splice variants?" | Per-transcript projection; missing isoforms documented |
| "BUSCO of projected annotation?" | >= 90% complete; reported |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| TOGA "chain file not found" | Cactus HAL not converted to chain | Use halSynteny + axtChain |
| TOGA stops with "no syntenic blocks" | Genome unrelated to reference | Verify Cactus alignment is non-empty |
| LiftOff produces empty GFF | Reference / query genome mismatch | Verify both are FASTA; check identifier consistency |
| GeMoMa OOM | Java heap insufficient | Increase via `-Xmx32g` |
| CESAR fragment-not-found | Exon GFF malformed | Verify GFF format; convert if needed |
| TOGA Nextflow hangs | Cluster resource issue | Restart with `-resume`; check Nextflow config |
| BUSCO of projected annotation low | Missing core genes | Likely projection failure; manual review |
| Many "PI" classifications | Assembly fragmentation | Improve assembly; or accept partial intactness |
| CAT Snakemake step fails | Conda env issue | Re-create conda envs; check Snakemake logs |
| LiftOff "no overlap" | Reference / query coordinate systems differ | Verify same genome version |
| TOGA intactness disagrees with biology | Edge case; manual review needed | Inspect CESAR alignment; cross-validate with RNA-Seq |

## Tool Installation Notes

```bash
# TOGA
conda env create -f https://raw.githubusercontent.com/hillerlab/TOGA/master/toga_env.yml
# Requires Nextflow and nf-core

# CESAR 2.0 (bundled with TOGA)
git clone https://github.com/hillerlab/CESAR2.0

# LiftOff
pip install liftoff

# UCSC liftOver
conda install -c bioconda ucsc-liftover

# Comparative Annotation Toolkit (CAT)
git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit
cd Comparative-Annotation-Toolkit && pip install .

# GeMoMa
wget https://gemoma.de/jcag/gemoma.zip && unzip gemoma.zip

# Comparative tools
conda install -c bioconda cactus busco compleasm

# De novo annotation (alternative)
conda install -c bioconda braker3 funannotate maker
```

For TOGA pipeline at vertebrate-scale, use Nextflow with proper HPC config (Slurm / Kubernetes); allocate >= 32 cores per genome.

## References

- Kirilenko BM et al 2023 Science 380:eabn3107 (TOGA; Zoonomia + Bird10000 standard)
- Sharma V & Hiller M 2017 NAR 45:8369 (CESAR 2.0)
- Shumate A & Salzberg SL 2021 Bioinformatics 37(12):1639-1643 (LiftOff)
- Hickey G et al 2013 Bioinformatics 29:1341 (HAL toolkit)
- Keilwagen J et al 2019 Methods Mol Biol 1962:161 (GeMoMa)
- Brůna T et al 2024 NAR Genom Bioinform 6:lqae068 (BRAKER3)
- Cantarel BL et al 2008 Genome Res 18:188 (MAKER)
- Stanke M et al 2008 Bioinformatics 24:637 (AUGUSTUS)
- Diekhans M et al 2007 Nat Genet 39:583 (TransMap)
- Zhang Z et al 2006 Genome Res 16:1041 (PseudoPipe)
- Armstrong J et al 2020 Nature 587:246 (Progressive Cactus)
- Baertsch R et al 2008 Genome Res 18:1675 (RetroFinder)
- Liao W-W et al 2023 Nature 617:312 (HPRC draft pangenome; relevant context)
- Smith LP et al 2024 (Comparative Annotation Toolkit updates)
- Salzberg SL 2019 Bioinformatics 35:1844 (next-gen annotation pipelines)
- Comparative Genomics Toolkit (CGT) GitHub
- TimeTree (database) for divergence dates

## Related Skills

- comparative-genomics/whole-genome-alignment - Cactus WGA precedes TOGA
- comparative-genomics/synteny-analysis - Synteny detection from WGA
- comparative-genomics/ortholog-inference - TOGA orthology classification
- comparative-genomics/pangenome-analysis - PGR-TK for repetitive / clinical genes
- comparative-genomics/whole-genome-duplication - Subgenome assignment for polyploid query
- genome-annotation/eukaryotic-gene-prediction - BRAKER3 / Funannotate de novo alternative
- genome-annotation/functional-annotation - Function assignment downstream
- genome-annotation/annotation-transfer - Related skill on annotation transfer mechanisms
- genome-annotation/prokaryotic-annotation - Bakta for prokaryote annotation
- read-qc/rnaseq-qc - RNA-Seq evidence to validate projected annotations
- read-alignment/star-alignment - RNA-Seq alignment for validation
