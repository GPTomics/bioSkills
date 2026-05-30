# GTF/GFF Handling - Usage Guide

## Overview
GTF and GFF3 are the two serializations of a gene-model tree: gene -> transcript -> exon/CDS/UTR, with GTF using `key "value";` attributes and an implicit `gene_id`/`transcript_id` hierarchy and GFF3 using `key=value` with explicit `ID`/`Parent` links. This skill covers parsing, querying, converting, and extracting from these files - walking the hierarchy with gffutils, converting formats and extracting transcript/CDS/protein FASTA with gffread, slurping to dataframes with gtfparse/pyranges, and sanitizing malformed files with AGAT. Its central discipline is reconciling the *seams*: the 1-based-vs-BED coordinate conversion, phase-vs-frame, the stop-codon convention, and the chromosome-name and gene-ID namespace mismatches that silently corrupt counts and joins. Every failure mode in this domain is silent - nothing throws, the wrong answer just propagates.

## Prerequisites
```bash
# CLI tools
conda install -c bioconda gffread agat

# Python packages
pip install gffutils gtfparse pyranges
```
gtfparse >=2.x returns a polars DataFrame by default (pass `result_type='pandas'` for pandas idioms). pyranges has a 0.x-vs-1.0 API split - check `pyranges.__version__`. AGAT is a Perl toolkit installed via bioconda.

## Quick Start
Tell your AI agent what you want to do:
- "Extract all protein-coding gene coordinates from my GTF as BED"
- "Get every exon for TP53 and derive its introns"
- "Convert my GENCODE GTF to GFF3, or extract protein FASTA"
- "My count matrix is all zeros - figure out why"
- "Sanitize this malformed GFF3 before I parse it"

## Example Prompts

### Extracting Features and Sequences
> "Walk the gene tree in gencode.v44.annotation.gtf, list every transcript of TP53 with its exon count, and derive the introns."
> "Extract spliced transcript, CDS, and protein FASTA from my annotation against the genome FASTA, keeping only coding transcripts."
> "Pull all genes with gene_biotype lncRNA into a table."

### Coordinate and Format Conversion
> "Convert my gene features to BED for bedtools - make sure the 1-based-to-0-based shift is correct (start only, not both ends)."
> "Convert this GENCODE GTF to GFF3, then back to GTF, and tell me whether the stop codon stays in or out of the CDS."

### Diagnosing Silent Failures
> "featureCounts gave me a matrix of zeros - intersect the BAM chromosome names with the GTF seqids and tell me if it is a chr1-vs-1 mismatch."
> "My DESeq2 results lost half their genes when I joined gene symbols - check whether the gene IDs have version suffixes on one side."
> "A CDS I extracted is 3 nt shorter than the reference - is this the stop-codon convention?"

### Sanitizing
> "This GFF3 has no ##gff-version line and is missing exon features - run AGAT to reconstruct the full tree before I parse it."

## What the Agent Will Do
1. Identify the provider (Ensembl/GENCODE/RefSeq/UCSC), since it determines chromosome naming, the biotype attribute key (`gene_biotype` vs `gene_type`), and which features are present.
2. Choose the tool by task: gffutils for hierarchy traversal, gffread for conversion/FASTA, gtfparse/pyranges for dataframe work, AGAT to sanitize a malformed file first.
3. For tree queries, build a gffutils DB (disabling gene/transcript inference when those lines exist, for a ~100x speedup) and walk children/parents, deriving introns/UTRs as needed.
4. For coordinate conversion, subtract 1 from the start only when going to BED, leaving the end unchanged.
5. Before any count or join, intersect the relevant key namespaces (chromosomes, gene IDs) programmatically and report mismatches rather than proceeding silently.

## Tips
- Convert to BED with `start-1`, end unchanged - over-correcting both ends frameshifts CDS translation and is invisible in coverage.
- gffutils keeps 1-based coordinates, pyranges stores 0-based - their `start` fields differ by one by design; never "fix" that.
- Set gffutils `disable_infer_genes=True, disable_infer_transcripts=True` on modern GTFs that already have gene/transcript lines; leave inference on for older minimal GTFs that lack them.
- A CDS length off by exactly 3 nt between two sources is the stop-codon convention (GTF excludes it, GenBank/GFF3 often include it), not a code bug.
- Never trust hand-edited phase - recompute with AGAT or gffread; treat a CDS edit and a phase recompute as one operation.
- Reach for AGAT (`agat_convert_sp_gxf2gxf.pl`) the moment a file lacks `##gff-version`, has non-SO types, is missing gene/exon/UTR lines, or has duplicate IDs.
- Set featureCounts/htseq `-s` from the library chemistry, never the tool default - featureCounts defaults unstranded, htseq-count defaults stranded.

## Resources
- [GENCODE](https://www.gencodegenes.org/) - Human/mouse annotations and format FAQ
- [Ensembl FTP](https://ftp.ensembl.org/) - Multi-species annotations
- [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) - Sequence Ontology
- [gffread documentation](https://ccb.jhu.edu/software/stringtie/gff.shtml)
- [gffutils documentation](https://daler.github.io/gffutils/)

## Related Skills
- bed-file-basics - BED format and the coordinate conversion this skill feeds into
- interval-arithmetic - Set operations on the features extracted here
- proximity-operations - Strand-aware TSS/promoter derivation from extracted features
- rna-quantification/featurecounts-counting - Consumes the GTF/GFF features; the seqid/strand landmines live there
- genome-annotation/functional-annotation - Downstream of feature/sequence extraction from the annotation
- genome-annotation/annotation-qc - Judges whether the annotation this skill parses is sound
- differential-expression/de-results - Map gene coordinates back to DE results
