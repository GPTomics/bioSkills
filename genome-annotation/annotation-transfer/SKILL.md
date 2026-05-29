---
name: bio-genome-annotation-annotation-transfer
description: Transfers gene annotations between genome assemblies via coordinate liftover (UCSC liftOver, CrossMap for same-species version updates) or feature/sequence projection (Liftoff for same/close species, miniprot for protein-level cross-species, TOGA/GeMoMa/CAT for distant clades). Covers the coordinate-vs-projection decision by divergence, why a successful lift is not biological confirmation, reference bias, the silent-dropping of unmapped features, build/PAR/MHC/inversion hazards, and transfer-vs-de-novo validation. Use when annotating a new assembly of a species with an existing reference, harmonizing coordinates across builds, or mapping annotations across related species.
tool_type: cli
primary_tool: Liftoff
---

## Version Compatibility

Reference examples tested with: Liftoff 1.6.3+, LiftoffTools 0.4+, miniprot 0.13+, CrossMap 0.7+, UCSC liftOver (current), BioPython 1.83+, gffutils 0.12+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

The **chain file must match the exact assembly pair** (build and patch); the **source and target build must be recorded** with every coordinate (a coordinate without a build is unusable). If code throws an error, introspect the installed tool and adapt rather than retrying.

# Annotation Transfer

**"Transfer annotations from a reference to my new assembly"** -> Map gene models from a well-annotated reference onto a target, by coordinate liftover (same-species, fast) or by re-aligning the actual gene sequence (cross-assembly/species, structure-aware), then validate against the target.
- CLI: `liftoff -g ref.gff3 -o out.gff3 -u unmapped.txt target.fa reference.fa` (note: target before reference), `liftOver in.bed map.chain out.bed unmapped` (intervals)

## The Single Most Important Modern Insight -- A Lift Is Geometry, Not Biology

Coordinate liftover and feature projection answer different questions, and choosing the wrong one is the dominant failure mode:

- **Coordinate liftover** (liftOver/CrossMap, on pre-computed chains) answers *"where does this interval sit in the other assembly's coordinate system?"* It moves numbers; it never re-examines the sequence.
- **Feature projection** (Liftoff/miniprot/TOGA) answers *"where is this gene, with its exon-intron structure intact, and is it still a functional gene?"* It re-aligns the biological sequence.

Three load-bearing consequences:

1. **A successful lift is not biological confirmation.** liftOver can place a gene at perfectly valid target coordinates that land in a pseudogenized, frameshifted, or collapsed-duplication region - the coordinate is right and the gene is dead, and the tool has no vocabulary to flag it. The output GFF inherits the reference's `gene`/`CDS` feature types verbatim. The "it mapped, ship it" culture is how a lifted GFF acquires the social status of a validated annotation while no one ever re-derived a model from sequence. Treat every lifted annotation as a hypothesis until target evidence (intact ORF, identity distribution, BUSCO, RNA-seq) has touched it.
2. **Transfer is reference-biased: it can only reproduce what the reference annotated.** Lineage-specific genes, target-specific expansions, novel isoforms, and orphan genes are *structurally invisible* - the new assembly's actual novelty (the reason it is interesting) is exactly what transfer cannot see. Always pair transfer with de novo + evidence; TOGA can call gene *loss* but is constitutionally incapable of calling gene *gain*.
3. **Read and classify the unmapped file - it is the most information-rich output.** liftOver writes failures to a side file, exits 0, and prints a clean shorter GFF with no visual tell; entire gene families can vanish silently. "Deleted in new" vs "Split in new" vs "Duplicated in new" are biologically distinct diagnoses (chain gap / rearrangement breakpoint / segmental duplication or gene family).

## Two Paradigms

| Paradigm | Tools | Operates on | Right for |
|----------|-------|-------------|-----------|
| A. Coordinate liftover | UCSC liftOver, CrossMap, segment_liftover, paftools | pre-computed chains; intervals (BED/GFF/VCF/BAM) | same-species version updates (hg19<->hg38, mm10<->mm39); variant/peak/CNV harmonization |
| B. Feature/sequence projection | Liftoff (nt), miniprot (protein), GeMoMa, TOGA, CAT, LiftOn | re-aligned gene sequence | cross-assembly/species; full gene models; polyploid/duplicated; no reliable chain |

## Decision Tree by Divergence

| Divergence | Recommended | Why |
|------------|-------------|-----|
| Same species, transfer intervals | liftOver / CrossMap | chain is dense; geometry suffices for variants/peaks |
| Same species, transfer gene models | Liftoff (`-chroms`) | structure-aware; per-interval liftOver fragments transcripts |
| Same genus (a few % divergence) | Liftoff + miniprot rescue for the divergent tail | nucleotide alignment robust; protein for the rest |
| Same family/order (tens-hundreds My) | TOGA or GeMoMa (multi-reference) | nucleotide saturates; orthology + gene-loss reasoning |
| Beyond family / lineage-specific content / heavy rearrangement | -> eukaryotic-gene-prediction (de novo) + transfer as evidence | reference too far; only de novo sees target-specific biology |
| Pan-genome / multi-haplotype | `vg annotate` onto the graph | avoids single-reference bias (tooling still maturing) |

Cross-species *coordinate* liftover is a methodological error (synteny fragments into thousands of short chains; most genes drop silently) - it is the wrong paradigm, not a tuning problem.

## Liftoff (Same / Close Species, Nucleotide)

```bash
liftoff -g reference.gff3 -o lifted.gff3 -u unmapped.txt -p 16 \
    -chroms chrom_map.txt -polish target.fasta reference.fasta
```

Positional args are **target first, then reference** (commonly swapped - a silent error). Liftoff extracts each gene's exon sequence, aligns with minimap2, and chooses the placement maximizing identity while preserving exon-intron structure. Key flags: `-a` (alignment coverage, default 0.5), `-s` (sequence identity, default 0.5), `-copies`/`-sc` (search for extra gene copies - a per-family decision, not a default), `-polish` (re-align to restore intact start/stop/splice, writes `*_polished.gff3`), `-exclude_partial`, `-chroms` (ordered chromosome mapping; reduces false cross-chromosome placements). LiftoffTools QCs the result (variants, synteny, copy-number changes). Same-species version updates should lift ≥99% - a 97% rate is a four-alarm signal of the wrong chain or coordinate-convention mismatch, not "pretty good."

## miniprot (Cross-Species, Protein)

```bash
miniprot -t 16 -d target.mpi target.fasta            # optional index
miniprot -Iut 16 --gff target.mpi proteins.faa > out.gff
```

Protein conserves far deeper than nucleotide (synonymous sites saturate), so miniprot works across species where Liftoff's nucleotide alignment fails. `-I` auto-sets max intron from genome length; `--gff` emits GFF3. **Frameshift and in-frame-stop tags in the output are the signal that the "gene" is pseudogenized in the target, not a clean ortholog** - inspect them; do not treat a miniprot hit as a functional gene by default. For a polished multi-reference, intron-aware annotation use GeMoMa (which reasons about intron-position conservation); for the DNA+protein hybrid use LiftOn.

## TOGA (Distant Species, Orthology + Gene Loss)

TOGA consumes a genome-alignment chain + reference BED12 and uses ML on chain features (including intronic/intergenic flanks - orthologs share flanking context, paralogs/retrocopies do not) to classify orthology (one2one ... one2zero) and **gene-loss/intactness** (intact / partially intact / lost / missing). It exists *precisely because* across deep time an inactivated gene still aligns - a coordinate lift reports the corpse as "present." Use TOGA for whole-clade ortholog projection; it does not discover target-specific novel genes (the reference-bias caveat of all of Paradigm B).

## Validating a Transfer with Python

**Goal:** Quantify transfer quality and, critically, check that lifted CDS are biologically intact, not just placed.

**Approach:** Compare gene counts for a transfer rate, then translate each lifted CDS from the target and check for a valid start, a single terminal stop, and correct length - coordinate success is not intactness.

```python
import gffutils
from Bio import SeqIO

def orf_integrity(lifted_gff, target_fasta):
    genome = SeqIO.to_dict(SeqIO.parse(target_fasta, 'fasta'))
    db = gffutils.create_db(lifted_gff, ':memory:', merge_strategy='merge')
    valid = total = 0
    for cds in db.features_of_type('CDS'):
        total += 1
        seq = genome[cds.seqid].seq[cds.start - 1:cds.end]
        if cds.strand == '-':
            seq = seq.reverse_complement()
        prot = seq.translate()
        if prot.startswith('M') and prot.endswith('*') and prot.count('*') == 1:
            valid += 1
    print(f'Intact ORFs: {valid}/{total} ({valid/total:.1%}) -- a clean lift can still land in a pseudogene')
    return valid, total
```

Also: read and **classify** the unmapped file (not just count it); run BUSCO on the *lifted protein set* and compare to the reference (a drop quantifies silently lost conserved genes); compare to a de novo annotation to expose reference bias.

## Hazards and the 1% That Matters

- **The 1% that does not lift trivially is enriched for the interesting parts.** ~99% of the human genome lifts hg19<->hg38 (liftOver ~99.99% over ~1.57M ClinVar variants), but failures concentrate in indels/large duplications (a pathogenic 8.1 kb *LDLR* duplication failed all three tools), PAR (double-maps), MHC, segmental duplications, centromeres/telomeres, and regions inverted between GRCh37/GRCh38 where coordinate conversion silently corrupts palindromic SNVs and imputation (the chr10 *MSMB*/rs10993994 prostate-cancer signal dropped p=2.86e-7 to 0.0011 across ~20,000 individuals; Sheng 2022). A lift touching any of these is suspect until validated.
- **Coordinate conventions and build identity bite.** BED is 0-based half-open; GFF/GTF/VCF are 1-based closed - a conversion at a format boundary shifts every start by 1 (passes smell tests, corrupts splice/start bases). VCF lifts must update the REF allele (CrossMap does; naive shifting does not). hg19 chrM (NC_001807) is *not* rCRS (NC_012920) - two files both labeled "hg19" can have incompatible mitochondrial coordinates.
- **Accumulation without revalidation.** Annotations get lifted v1->v2->v3 forever; a model wrong in 2009 stays byte-for-byte wrong in 2024 because lifting copies models without re-examining them. Re-run de novo + evidence at major assembly upgrades and reconcile.

## Per-Method Failure Modes

### Cross-species coordinate liftover
**Trigger:** liftOver/CrossMap to transfer genes between species. **Mechanism:** synteny fragments into short chains; most genes have no co-linear counterpart. **Symptom:** plausible-looking output that silently dropped most genes. **Fix:** miniprot/TOGA/GeMoMa (sequence/orthology), not chains.

### Not reading the unmapped file
**Trigger:** reporting a transfer complete from the success file alone. **Mechanism:** failures go to a side file; exit code 0. **Symptom:** a clean GFF missing entire gene families. **Fix:** read and classify unmapped (deletion/split/duplicated).

### Coordinate success treated as intactness
**Trigger:** trusting a lifted gene because it placed. **Mechanism:** the locus can be pseudogenized/frameshifted. **Symptom:** RNA-seq quantified against a gene with an internal stop at residue 40. **Fix:** `-polish` + ORF check; miniprot frameshift tags; TOGA intactness class.

### Swapped Liftoff positional args
**Trigger:** `liftoff ... reference.fa target.fa`. **Mechanism:** Liftoff is `target reference`. **Symptom:** nonsense mapping. **Fix:** target first, reference last.

### Loosening `-a`/`-s` to "rescue" features
**Trigger:** lowering coverage/identity to clear the unmapped pile. **Mechanism:** a 35%-coverage hit is usually a paralog/pseudogene/repeat match. **Symptom:** low-confidence placements laundered into the success file. **Fix:** treat default-threshold failures as signal; loosen only with a biological hypothesis and validate rescues individually.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Same-species mapping rate ≥99% | Liftoff/ClinVar studies | a 97% rate signals wrong chain / convention mismatch |
| Liftoff `-a`/`-s` default 0.5 | Liftoff | loosening manufactures false placements; failure is often signal |
| liftOver `-minMatch` default 0.95 (per-feature) | UCSC | a long feature with one chain gap fails silently |
| BUSCO on lifted set vs reference | completeness audit | a drop quantifies silently lost conserved genes |
| Divergence rule: species->liftOver/Liftoff; genus->+miniprot; family->TOGA/GeMoMa; beyond->de novo | lab convention | matches paradigm to where the chain/identity breaks |
| Every coordinate carries its build | reproducibility | a coordinate without a build is unusable |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Many unmapped features (same species) | wrong/patch-mismatched chain; contig naming (`chr1` vs `1`) | use the exact-pair chain; harmonize names |
| Mass gene loss, clean GFF | silent dropping | read/classify the unmapped file |
| Lifted genes with internal stops | landed in pseudogene/frameshift | `-polish`; ORF check; re-predict de novo in problem loci |
| Paralog/copy collapse or swap | no `-copies`, or mapped to the paralog | `-copies`/`-sc` per family; TOGA orthology graph |
| Most genes lost cross-species | coordinate liftover used across species | switch to miniprot/TOGA |
| mtDNA coordinates don't match | hg19 chrM != rCRS | record the exact MT record, not just "hg19" |

## References

- Shumate A, Salzberg SL. 2021. Liftoff: accurate mapping of gene annotations. *Bioinformatics* 37:1639-1643.
- Shumate A, Salzberg SL. 2024. LiftoffTools: a toolkit for comparing gene annotations mapped between genome assemblies. *F1000Research* 11:1230.
- Li H. 2023. Protein-to-genome alignment with miniprot. *Bioinformatics* 39:btad014.
- Zhao H, et al. 2014. CrossMap: a versatile tool for coordinate conversion between genome assemblies. *Bioinformatics* 30:1006-1007.
- Hinrichs AS, et al. 2006. The UCSC Genome Browser Database: update 2006. *Nucleic Acids Res* 34:D590-D598.
- Kent WJ, et al. 2003. Evolution's cauldron: duplication, deletion, and rearrangement in the mouse and human genomes (chains and nets). *PNAS* 100:11484-11489.
- Keilwagen J, et al. 2016. Using intron position conservation for homology-based gene prediction (GeMoMa). *Nucleic Acids Res* 44:e89.
- Otto TD, et al. 2011. RATT: Rapid Annotation Transfer Tool. *Nucleic Acids Res* 39:e57.
- Kirilenko BM, et al. 2023. Integrating gene annotation with orthology inference at scale (TOGA). *Science* 380:eabn3107.
- Fiddes IT, et al. 2018. Comparative Annotation Toolkit (CAT): simultaneous clade and personal genome annotation. *Genome Res* 28:1029-1038.
- Chao KH, et al. 2025. Combining DNA and protein alignments to improve genome annotation with LiftOn. *Genome Res* 35:311-325.
- Sheng X, Xia L, Cahoon JL, Conti DV, Haiman CA, Kachuri L, Chiang CWK. 2022. Inverted genomic regions between reference genome builds in humans impact imputation accuracy and decrease the power of association testing. *HGG Adv* 4:100159.

## Related Skills

- eukaryotic-gene-prediction - De novo annotation; the complement that captures target-specific genes
- annotation-qc - BUSCO on the lifted set and gene-structure integrity checks
- comparative-genomics/ortholog-inference - Orthology relationships across species
- comparative-genomics/synteny-analysis - Synteny context that validates or refutes a transfer
- genome-intervals/gtf-gff-handling - Parse and manipulate transferred annotations
