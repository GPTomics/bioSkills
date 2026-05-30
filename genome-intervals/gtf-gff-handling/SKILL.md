---
name: bio-genome-intervals-gtf-gff-handling
description: Parses, queries, converts, and extracts from GTF and GFF3 gene-model annotation files - walking the gene/transcript/exon/CDS hierarchy with gffutils (queryable SQLite DB), converting formats and extracting transcript/CDS/protein FASTA with gffread, slurping to dataframes with gtfparse/pyranges, and sanitizing malformed files with AGAT. Covers the 1-based-inclusive vs 0-based BED coordinate conversion (start-1 only), deriving implicit features (introns/UTRs/TSS), phase-not-frame, the stop-codon-in-or-out-of-CDS convention, and the chr1-vs-1 seqid and gene-ID-version mismatches that silently produce all-zero count matrices and dropped joins. Use when extracting features or sequences from an annotation, converting GTF<->GFF3 or GTF->BED, traversing the gene tree, or diagnosing a coordinate/provenance mismatch upstream of counting or DE.
tool_type: mixed
primary_tool: gffutils
---

## Version Compatibility

Reference examples tested with: gffutils 0.13+, gffread 0.12+, gtfparse 2.x, pyranges 0.1+ (or 1.0+ - see note), AGAT 1.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

Two version landmines specific to this skill: (1) **gtfparse changed its return type** - older releases returned a pandas DataFrame, gtfparse >=2.x returns a **polars** DataFrame by default; pass `result_type='pandas'` before chaining pandas idioms (`.copy()`, boolean masks). (2) **pyranges has a major-version API split** - pyranges 0.x and the 1.0 rewrite differ in method names and attribute access; check `import pyranges; pyranges.__version__` before pasting code. gffutils stores **1-based** coordinates while pyranges stores **0-based** - their `start` fields differ by one by design. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt rather than retrying.

# GTF/GFF Handling

**"Pull these features (or their sequences) out of my annotation, convert it, or find out why my counts are wrong."** -> Treat the file as a serialized gene-model tree: walk gene->transcript->exon/CDS, derive implicit features, and reconcile coordinate and namespace conventions before trusting any number.
- CLI: `gffread in.gtf -T -o out.gtf` (convert), `gffread -w tx.fa -g genome.fa in.gtf` (FASTA), `agat_convert_sp_gxf2gxf.pl` (sanitize)
- Python: `gffutils.create_db(...)` then `db.children(gene, featuretype='exon')` (tree query); `gtfparse.read_gtf(..., result_type='pandas')` / `pyranges.read_gtf(...)` (dataframe)

## The Single Most Important Modern Insight -- A GTF/GFF3 Is a Serialized Gene-Model Tree, Not a Table of Intervals

Almost every painful bug here comes from the file looking like a CSV while behaving like a tree, or from a coordinate/provenance mismatch the tools never warn about - and **every one of these failures is silent**: nothing throws, the wrong answer just propagates. Three load-bearing facts the tutorials skip:

1. **The coordinate conversion is asymmetric.** GTF/GFF3 are 1-based fully inclusive `[start, end]`; BED (and pyranges-internal) are 0-based half-open `[start-1, end)`. Convert to BED by **subtracting 1 from the start only - the end is unchanged** (the inclusive 1-based end and the exclusive 0-based end are the same integer). Doing `start-1` AND `end-1` shifts the feature one base left and is the classic over-correction: invisible in coverage/overlap, catastrophic in CDS translation (a one-base frameshift garbles the protein). gffutils keeps 1-based, pyranges stores 0-based, so their `start` fields differ by one *correctly* - never "fix" that discrepancy.

2. **The all-zero count matrix.** featureCounts/htseq-count match a read to a feature by **string equality on the chromosome name**, so `chr1` != `1` != `NC_000001.11` produces a perfectly well-formed matrix of **zeros with no error or warning** - the only signal is `~0%` assigned in the summary. Same bug one altitude up: gene-ID version suffixes (`ENSG00000223972.5` vs `ENSG00000223972`) silently drop rows on an annotation join. Audit every cross-file key (chromosomes between BAM/GTF/FASTA, gene IDs between GTF/count-matrix/annotation) by **set intersection, never by eye**, before any count or join.

3. **phase is not frame, and the stop codon is a 3-bp ghost.** Phase (column 8) is the strand-aware count of bases to trim from the segment's transcriptional 5' end to reach the next codon (0/1/2) - **recompute it on any CDS edit** (AGAT/gffread do; hand-editing coordinates without fixing phase frameshifts the translation). GTF (Ensembl/GENCODE) **excludes** the stop codon from CDS; GenBank/GFF3 often **include** it - so a CDS length off by exactly 3 nt (or a protein +/-1 stop) between two sources is a convention mismatch, not a bug.

## Tool Taxonomy

| Tool | Role | Mechanism | When |
|------|------|-----------|------|
| gffutils | Queryable gene-tree DB (Python) | builds a SQLite DB; `children`/`parents`/`region` traverse the hierarchy; keeps 1-based coords | walk gene->transcript->exon/CDS, derive introns, query by ID/coordinate |
| gffread | Converter + sequence extractor (CLI) | fast C++, genome-aware; one binary | GTF<->GFF3, extract transcript/CDS/protein FASTA, region filter |
| pyranges | Vectorized interval engine (Python) | PyRanges/pandas-like; stores 0-based half-open | overlap joins, set ops, dataframe-native interval work |
| gtfparse | GTF -> dataframe (Python) | one call explodes column 9 into attribute columns | quick column/filter work; NOT hierarchy-aware (flat table) |
| AGAT | GFF/GTF sanitizer (Perl CLI) | reconstructs the full tree; adds missing features, fixes IDs/phase, deflates attributes | a malformed/non-standard file - run FIRST, before parsing |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Walk the gene/transcript/exon hierarchy, derive introns | gffutils `create_db` + `children`/`parents` | the hierarchy is the point; flat parsers lose it |
| Convert GTF<->GFF3 or extract transcript/CDS/protein FASTA | gffread (`-T`, `-w`/`-x`/`-y -g`) | genome-aware, knows the stop-codon convention |
| Quick column/filter on a clean modern GTF | gtfparse (`result_type='pandas'`) | one-call dataframe; verify return type first |
| Overlap/set ops, large in-memory interval joins | pyranges | vectorized; route arithmetic -> interval-arithmetic |
| File malformed: no `##gff-version`, missing gene/exon lines, dup IDs, mixed conventions | AGAT `agat_convert_sp_gxf2gxf.pl` first | sanitize once vs writing a brittle parser around it |
| GTF -> BED for bedtools | `start-1`, end unchanged (-> bed-file-basics) | the off-by-one boundary is where it bites |
| Counts came out all-zero or DE join dropped rows | intersect seqid / gene-ID namespaces | string-equality match; no error is emitted |
| Counting reads per gene/feature | -> rna-quantification/featurecounts-counting | the seqid/strand landmines live there; set `-s` from chemistry |
| Judge whether the annotation itself is sound | -> genome-annotation/annotation-qc | this skill operates on the file, not its quality |

## Walk the Gene Tree and Derive Introns (gffutils)

**Goal:** Traverse gene -> transcript -> exon and reconstruct features (introns) that the file does not store explicitly.

**Approach:** Build a SQLite DB once (disabling gene/transcript inference when those lines already exist, for a ~100x speedup), then query children ordered by position and synthesize introns from the exon gaps.

```python
import gffutils

# disable_infer_* is GTF-only and applies when gene/transcript lines ALREADY exist (modern GENCODE/Ensembl) -> ~100x faster
db = gffutils.create_db('annotation.gtf', 'annotation.db', force=True,
                        disable_infer_genes=True, disable_infer_transcripts=True,
                        merge_strategy='create_unique')

gene = db['ENSG00000141510']                                    # gffutils returns 1-based coords (raw record)
for tx in db.children(gene, featuretype=['mRNA', 'transcript'], order_by='start'):
    exons = list(db.children(tx, featuretype='exon', order_by='start'))
    introns = list(db.interfeatures(exons, new_featuretype='intron'))   # introns are not stored - derived from exon gaps
    print(tx.id, len(exons), 'exons', len(introns), 'introns')
```

A modern GTF without `disable_infer_*` triggers the slow inference/merge machinery; an *older* minimal GTF lacking gene/transcript lines needs inference ON so gffutils reconstructs the envelopes. Match the flag to the file. introns, UTRs (`exon - CDS`), and TSS are derived, not stored - never infer biological absence from a missing feature line.

## Convert Formats and Extract Sequences (gffread)

gffread is genome-aware and respects the stop-codon convention, so it is the safe path for sequence extraction (naive coordinate math is not).

```bash
gffread annotation.gff3 -T -o annotation.gtf            # GFF3 -> GTF2 (default output is GFF3)
gffread -w transcripts.fa -g genome.fa annotation.gtf   # spliced exon (mature transcript) FASTA
gffread -x cds.fa         -g genome.fa annotation.gtf   # spliced CDS nucleotide FASTA
gffread -y proteins.fa    -g genome.fa annotation.gtf   # translated-CDS protein FASTA
gffread annotation.gtf -C -o coding.gtf                 # keep only coding transcripts
```

`-g` needs the genome FASTA (gffread auto-creates the `.fai`). `-w`/`-x`/`-y` splice the segments per transcript, so they handle multi-exon models correctly - do not concatenate exon FASTAs by hand.

## Convert GTF to BED with the Right Coordinate Shift

**Goal:** Emit a BED of a chosen feature type for bedtools, without the off-by-one frameshift.

**Approach:** Parse to a pandas frame, filter to the feature type, subtract 1 from the start *only*, leave the end untouched.

```python
import gtfparse

df = gtfparse.read_gtf('annotation.gtf', result_type='pandas')   # gtfparse >=2.x defaults to POLARS - force pandas
genes = df[df['feature'] == 'gene'].copy()
genes['start'] = genes['start'] - 1                              # 1-based inclusive -> 0-based half-open: START ONLY
bed = genes[['seqname', 'start', 'end', 'gene_id', 'score', 'strand']]
bed.to_csv('genes.bed', sep='\t', header=False, index=False)
```

For TSS/promoter derivation (strand-aware: `+` strand TSS = start, `-` strand TSS = end), route to proximity-operations - the promoter window is an imposed definition, not an annotated feature.

## Sanitize a Malformed File First (AGAT)

When a file lacks `##gff-version 3`, has non-Sequence-Ontology types, is missing `gene`/`exon`/UTR lines, has duplicate IDs, or mixes conventions, sanitize it once rather than coding around it:

```bash
agat_convert_sp_gxf2gxf.pl -g messy.gff3 -o clean.gff3   # adds missing ID/Parent + features, fixes dup IDs, recomputes phase, sorts
agat_convert_sp_gff2gtf.pl -g clean.gff3 -o clean.gtf    # GFF3 -> GTF (collapses level1->gene, level2->transcript)
```

AGAT *makes decisions* (which convention to standardize to, how to derive missing features) - usually a feature, but when a source's exact encoding must be preserved (e.g. auditing a submission), inspect what it changed rather than trusting blindly.

## Per-Method Failure Modes

### Over-correcting the coordinate conversion
**Trigger:** subtracting 1 from both start and end when converting to BED. **Mechanism:** only the start representation differs; the inclusive 1-based end equals the exclusive 0-based end. **Symptom:** every feature shifted one base left; invisible in coverage, frameshifts CDS translation. **Fix:** `start-1`, end unchanged.

### Comparing coordinates across gffutils and pyranges
**Trigger:** asserting equality on `start` fields from both libraries in one script. **Mechanism:** gffutils keeps 1-based, pyranges stores 0-based. **Symptom:** an off-by-one that looks like a bug; "fixing" it introduces a real error. **Fix:** confirm each library's convention; expect the difference.

### All-zero count matrix (seqid mismatch)
**Trigger:** BAM aligned to `chr1`, GTF annotated with `1`. **Mechanism:** counters match reads to features by chromosome-name string equality. **Symptom:** well-formed matrix of zeros, no error; `~0%` assigned in the summary. **Fix:** intersect the BAM `@SQ`/idxstats chromosome set with the GTF column-1 set programmatically; remap one namespace, re-confirm.

### Dropped rows on a gene-ID join
**Trigger:** count matrix keyed `ENSG...` joined to annotation keyed `ENSG....5`. **Mechanism:** exact string match on a versioned vs unversioned ID. **Symptom:** join returns a dataframe but rows vanish / annotation is NA. **Fix:** strip `.\d+$` on both sides for matching; keep the version in the stored annotation for provenance.

### CDS length off by exactly 3
**Trigger:** comparing CDS/protein length across two sources, or translating after a convention-flipping conversion. **Mechanism:** GTF excludes the stop codon from CDS; GenBank/GFF3 often include it. **Symptom:** length differs by 3 nt / 1 aa; protein does/does not end in `*`. **Fix:** suspect the convention before debugging code; extract CDS with gffread/AGAT, which know it.

### Editing CDS coordinates without recomputing phase
**Trigger:** trimming/merging/lifting CDS coordinates, leaving column 8 as-is. **Mechanism:** phase is a static integer; the chain of per-segment phases depends on cumulative coding length. **Symptom:** downstream translation (gffread `-y`, table2asn, EMBL) frameshifts or rejects. **Fix:** treat a CDS edit + phase recompute as one atomic operation; let AGAT/gffread recompute.

### gffutils pathologically slow on a modern GTF
**Trigger:** `create_db` on a GENCODE/Ensembl GTF without the infer flags. **Mechanism:** gffutils infers gene/transcript envelopes and runs the merge machinery. **Symptom:** create_db hangs for many minutes. **Fix:** `disable_infer_genes=True, disable_infer_transcripts=True` when those lines already exist (~100x faster).

## Quantitative Thresholds

| Convention / threshold | Source | Rationale |
|------------------------|--------|-----------|
| GTF/GFF3 1-based inclusive; convert to BED with start-1, end unchanged | UCSC/SO format specs | the inclusive 1-based end == the exclusive 0-based end; over-correcting both ends frameshifts CDS |
| CDS length differs by exactly 3 nt between sources | GTF vs GenBank/GFF3 stop-codon convention | GTF (Ensembl/GENCODE) excludes the stop from CDS; GenBank/GFF3 often include it |
| `disable_infer_*` -> ~100x create_db speedup | gffutils docs | inference/merge machinery is skipped when gene/transcript lines already exist |
| seqid intersection required before counting | featureCounts/htseq string-equality match | non-overlapping chromosome names -> all-zero matrix with no error |
| Strip `.\d+$` from gene IDs on both sides before a join | Ensembl/GENCODE/RefSeq versioned accessions | version suffix tracks model revision; mismatch drops rows silently |
| featureCounts default `-s 0` (unstranded) vs htseq-count `-s yes` (stranded) | tool defaults (Liao 2014; Anders 2015) | switching tools changes the counting model; set `-s` from library chemistry, not the default |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| All genes count zero | seqid mismatch (`chr1` vs `1` vs `NC_...`) | intersect BAM and GTF chromosome sets; remap one namespace |
| Counts low and flip when `-s` changes | wrong strandedness (featureCounts vs htseq defaults differ) | set `-s` from the library prep chemistry; verify assignment rate |
| Join drops rows / NA annotation | gene-ID version suffix (`ENSG....5` vs `ENSG...`) | strip `.\d+$` on both sides for matching |
| Biotype filter returns empty | attribute key differs by source: GENCODE uses `gene_type`, Ensembl/RefSeq use `gene_biotype` | check the actual key (it travels with the `chr1`-vs-`1` provenance split); query the present key |
| Translated protein is garbage | over-corrected coordinate (`start-1` AND `end-1`) | subtract 1 from start only |
| CDS/protein off by 3 nt / 1 aa | stop-codon-in-or-out-of-CDS convention | extract with gffread/AGAT; do not debug coordinate math |
| gtfparse pandas idioms raise AttributeError | gtfparse >=2.x returns polars | pass `result_type='pandas'` |
| pyranges AttributeError | 0.x vs 1.0 API mismatch | check `pyranges.__version__`; use matching method names |
| gffutils create_db hangs | infer machinery on a modern GTF | set `disable_infer_genes=True, disable_infer_transcripts=True` |
| `gffread -w/-x/-y` errors | missing or unindexed genome FASTA | pass `-g genome.fa` (gffread creates the `.fai`) |

## References

- Pertea G, Pertea M. 2020. GFF Utilities: GffRead and GffCompare. *F1000Research* 9:304.
- Stovner EB, Saetrom P. 2020. PyRanges: efficient comparison of genomic intervals in Python. *Bioinformatics* 36:918-919.
- Dale R. gffutils: GFF and GTF file manipulation and interconversion. Software, https://github.com/daler/gffutils (no journal publication).
- Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. doi:10.5281/zenodo.3552717.
- Rubinsteyn A, et al. gtfparse: parsing tools for GTF (gene transfer format) files. Software, https://github.com/openvax/gtfparse (no journal publication).
- Liao Y, Smyth GK, Shi W. 2014. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30:923-930.
- Anders S, Pyl PT, Huber W. 2015. HTSeq - a Python framework to work with high-throughput sequencing data. *Bioinformatics* 31:166-169.

## Related Skills

- bed-file-basics - BED format and the coordinate conversion this skill feeds into
- interval-arithmetic - Set operations on the features extracted here
- proximity-operations - Strand-aware TSS/promoter derivation from extracted features
- rna-quantification/featurecounts-counting - Consumes the GTF/GFF features; the seqid/strand landmines live there
- genome-annotation/functional-annotation - Downstream of feature/sequence extraction from the annotation
- genome-annotation/annotation-qc - Judges whether the annotation this skill parses is sound
- differential-expression/de-results - Map gene coordinates back to DE results
