---
name: bio-genome-annotation-prokaryotic-annotation
description: Annotates bacterial and archaeal genomes (isolates, MAGs, plasmids) with Bakta (active versioned databases, NCBI-compliant output) or Prokka (legacy), producing GFF3/GenBank/EMBL/FASTA with INSDC locus tags. Covers Bakta-vs-Prokka-vs-PGAP-vs-DFAST choice, light-vs-full database tiers, translation-table selection (11/4/25), archaeal and leaderless-gene caveats, the small-ORF blind spot, pseudogene-vs-phase-variation, the pangenome re-annotation trap, and submission compliance. Use when annotating a newly assembled prokaryotic genome, choosing an annotation tool, re-annotating a collection for pangenomics, or preparing annotations for NCBI/DDBJ submission.
tool_type: cli
primary_tool: Bakta
---

## Version Compatibility

Reference examples tested with: Bakta 1.9+, Prokka 1.14.6 (legacy), tRNAscan-SE 2.0+, CheckM2 1.0+, gffutils 0.12+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

Annotation content tracks the **database version**, not just the binary: a Bakta full DB (schema v5+, ~30 GB zipped) and a Bakta light DB give different functional calls, and two runs months apart can differ purely from DB updates. Record the Bakta DB version in methods. Prokka's bundled databases are frozen (~2019-2021). If code throws an error, introspect the installed tool and adapt rather than retrying.

# Prokaryotic Genome Annotation

**"Annotate my bacterial genome"** -> Call protein-coding genes, tRNAs/rRNAs/ncRNAs, and other features, then assign function by database identity, and emit submission-ready files.
- CLI: `bakta --db db/ assembly.fa` (default), `prokka --outdir out assembly.fa` (legacy), NCBI PGAP (submission/RefSeq-grade)

## The Single Most Important Modern Insight -- Gene Calling Is Near-Solved; Function Is Not, and Wrong Labels Propagate

For a typical isolate, Prodigal/Pyrodigal recovers >95-99% of true coding genes. The unsolved problems are the **start codon** (translation initiation site), the **small/overlapping/recoded ORFs**, and above all the **function**. Three consequences a postdoc must internalize:

1. **A confidently wrong product name is worse than "hypothetical protein."** Names assigned by loose homology transfer across distant lineages and *self-amplify* through databases - there is no mechanism to retract a correction once it spreads, so annotation accuracy has gone *down*, not up, as sequencing scaled (Salzberg 2019 *Genome Biol* 20:92). Bakta's design counters this: gene-level identity only via exact (MD5/UniRef100) match, falling back to cluster level, then "hypothetical" rather than guessing. 25-50% hypothetical is healthy for a non-model organism; near 0% means over-confident transfer.

2. **Comparing gene counts across tools/versions/DBs is invalid.** Tool agreement is "highly dependent on the organism of study" and biased toward model organisms (Dimonaco 2022 *Bioinformatics* 38:1198) - there is no universal best tool. Differences between two annotations are mostly tool artifacts, not biology. For any collection (pangenomics), **re-annotate every assembly from FASTA with one pipeline + one pinned DB**; merging published annotations inflates the accessory genome ~10× (Tonkin-Hill 2020 *Genome Biol* 21:180).

3. **Annotation completeness ≠ assembly completeness.** A fragmented/contaminated assembly produces truncated partial CDS at every contig break, inflated counts, and missing rRNA operons. Run CheckM2 *before* trusting any annotation QC number.

## Tool Taxonomy

| Tool | Maintained | Database approach | Output | When |
|------|-----------|-------------------|--------|------|
| Bakta | Yes (active) | Curated, **versioned**, alignment-free (UniRef + AMRFinderPlus + expert systems) | GFF3/GBFF/EMBL/FASTA/TSV/JSON + plot | **Default for new work**; reproducible, MAG-aware |
| Prokka | Frozen ~2021 | BLAST hierarchy vs frozen UniProt/RefSeq + HMM | GFF3/GBK/FAA + tbl | Legacy only; pangenome pipelines (Roary) expect Prokka GFF |
| NCBI PGAP | Yes (NCBI) | RefSeq protein-family models + ProSplign | ASN.1/GenBank, submission-ready | **GenBank/RefSeq submission**; best pseudogene/frameshift/selenoprotein handling |
| DFAST | Yes (DDBJ) | DFAST reference DBs + swappable refs | INSDC files | DDBJ submission; fast, flexible reference swap |
| RAST / BV-BRC | Yes | SEED subsystems | GenBank/GFF + subsystems | Subsystem/metabolic framing; web; not for direct INSDC submission |

Prokka uses Aragorn (tRNA) + barrnap (rRNA); Bakta uses tRNAscan-SE 2.0 (tRNA, domain-specific) + Infernal/Rfam (rRNA/ncRNA) + PILER-CR (CRISPR arrays) + DeepSig (signal peptides). Bakta is deliberately conservative on naming, so it reports "hypothetical" where Prokka confidently (sometimes wrongly) names a gene - Bakta looking "less annotated" is it being honest.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Routine isolate, reproducible | Bakta, full DB | standardized, versioned, current functional calls |
| Submitting to GenBank/RefSeq | PGAP | RefSeq re-annotates with PGAP regardless; best pseudogene handling |
| Submitting to DDBJ | DFAST | INSDC-ready, DDBJ-aligned |
| MAG / metagenome bin | Bakta `--meta` + CheckM2 gate | anonymous-mode calling; QC before trust |
| Mycoplasma / Mollicutes | Bakta `--translation-table 4` | UGA=Trp; table 11 splits genes at every UGA |
| Archaeon | Bakta + verify tRNAscan-SE archaeal model; PGAP if N-termini matter | leaderless mRNAs degrade Prodigal TIS |
| Inherited Prokka pangenome pipeline | pin Prokka version, or re-call all in Bakta | tool consistency vs current biology |
| Subsystem/metabolic view | RAST / BV-BRC | SEED subsystem categories |
| Genome not yet assembled / poor QC | -> genome-assembly/assembly-qc | fix assembly before annotating |
| AMR for clinical report | -> run AMRFinderPlus/CARD-RGI directly | `--organism` context, point mutations, partials matter |

## Bakta (Default)

```bash
# Database (record the version): full ~30 GB for publishable annotation; light for triage/CI
bakta_db download --output db/ --type full

bakta --db db/ --output bakta_out --prefix ecoli_k12 \
    --genus Escherichia --species coli --strain K-12 \
    --locus-tag ECK12 --gram - --complete --threads 16 \
    assembly.fasta
```

Key flags: `--translation-table {11,4,25}` (default 11), `--gram {+,-,?}` (gates DeepSig signal-peptide calls; default `?`), `--complete` (all sequences are finished replicons; enables oriC detection - do NOT use on draft contigs), `--meta` (metagenome/MAG mode), `--compliant` (enforce INSDC structure), `--keep-contig-headers`, `--proteins <faa>` (trusted-protein transfer). Set `--genus`/`--species` from a GTDB-Tk classification, not a guess. Primary outputs: `.gff3`, `.gbff`, `.faa`, `.ffn`, `.fna`, `.tsv`, plus `.hypotheticals.tsv` and `.inference.tsv` (open the inference column to ask *why* a product was assigned).

## Prokka (Legacy)

```bash
prokka --outdir prokka_out --prefix my_genome --locustag MYORG \
    --genus Escherichia --species coli --cpus 8 --rfam assembly.fasta
```

Use only for tool-chain consistency with an existing Prokka-based pangenome workflow, and pin the version. Its bundled databases are frozen: a gene family characterized after ~2019 is "hypothetical" in Prokka but named by current Bakta/PGAP, so the *same gene* flips core/accessory purely on DB vintage.

## Coding Density and CDS Extraction with Python

**Goal:** Load Bakta/Prokka GFF3 into a queryable database and compute coding density, the first sanity number.

**Approach:** Build a gffutils database, sum CDS lengths, divide by genome length; flag values outside the expected band.

```python
import gffutils

CODING_DENSITY_LOW = 0.85   # <0.85 in a free-living bacterium suggests wrong table, fragmentation, or heavy pseudogenization
CODING_DENSITY_HIGH = 0.93  # >0.93 suggests ORF over-calling (spurious short hypotheticals)

def coding_density(gff_file, genome_length):
    db = gffutils.create_db(gff_file, ':memory:', merge_strategy='merge')
    coding_bp = sum(c.end - c.start + 1 for c in db.features_of_type('CDS'))
    density = coding_bp / genome_length
    if density < CODING_DENSITY_LOW or density > CODING_DENSITY_HIGH:
        print(f'WARNING: coding density {density:.1%} outside expected 85-93%')
    return density
```

## Hard Cases the Caller Gets Wrong

- **Wrong genetic code is silent and looks like fragmentation.** A Mycoplasma run under table 11 yields anomalously low coding density + many short "hypothetical" fragments because every internal UGA split a gene. Confirm the table from taxonomy (GTDB-Tk), never a guess. Table 4 (UGA=Trp, Mollicutes), table 25 (UGA=Gly, Gracilibacteria/SR1).
- **Pseudogene over-calling in reductive genomes is real biology.** *Mycobacterium leprae* has ~1,116 pseudogenes vs ~1,604 intact CDS (Cole 2001 *Nature* 409:1007) - high pseudogene fraction in a symbiont/host-restricted pathogen (*Rickettsia*, *Sodalis*) is a lifestyle signal, not a defect. PGAP (ProSplign aligns *through* frameshifts, emits `/pseudo`) is the best automated arbiter.
- **Phase variation is not a pseudogene.** A frameshift inside a homopolymer/SSR tract in a contingency locus (*Neisseria* ~65 candidate loci, *Haemophilus*, *Campylobacter* poly-G tracts) means the assembled cell was in the OFF phase - do NOT "correct" or polish it away. This confounds with ONT homopolymer indel error: a pseudogene spike is ambiguous between reductive biology, phase variation, and basecalling error; disambiguate with orthogonal (short-read/HiFi) data.
- **Programmed frameshifts make one gene look like two.** `prfB`/RF2 (a **+1** frameshift at a slippery `CTTT` + internal SD), `dnaX` (−1), and IS-element transposases encode one protein across two frames; naive callers emit two short ORFs. PGAP recognizes a curated set.
- **Small ORFs (<~50 aa) are systematically missed** - callers impose a ~30 aa minimum because short ORFs arise by chance and coding statistics are unreliable there. *E. coli* K-12 was missing dozens of 16-50 aa proteins. Treat the gene count as a lower bound on the small proteome; for sORF biology use a dedicated caller (smORFer, smORFinder) and, ideally, Ribo-seq (the experimental arbiter of translation).
- **Archaea and Actinobacteria use leaderless mRNAs** (no 5' UTR, no Shine-Dalgarno) - Prodigal's RBS-first scoring degrades TIS placement; GeneMarkS-2 (and PGAP) model leaderless transcription. Use tRNAscan-SE's archaeal model for archaeal intron-containing tRNAs.

## Per-Method Failure Modes

### Cross-tool gene-count comparison
**Trigger:** comparing counts from Bakta vs Prokka vs old RefSeq, or mixed-vintage records. **Mechanism:** tool/DB differences dominate biological differences. **Symptom:** "novel genes" or accessory-genome inflation. **Fix:** re-annotate every genome from FASTA with one pipeline + pinned DB.

### Wrong translation table
**Trigger:** table 11 on a Mollicute. **Mechanism:** UGA read as stop. **Symptom:** low coding density, doubled short gene count, high hypothetical fraction. **Fix:** `--translation-table 4`; confirm from GTDB-Tk.

### Annotating an unvetted assembly
**Trigger:** running Bakta before CheckM2. **Mechanism:** contig breaks truncate CDS; contamination mixes organisms. **Symptom:** partial CDS at ends, inflated/chimeric gene set, missing rRNA operons. **Fix:** CheckM2 first; contamination >5% -> decontaminate.

### Trusting an over-specific product name
**Trigger:** reading the product column as ground truth. **Mechanism:** loose-homology transfer below ~40% identity is unreliable, especially for promiscuous folds. **Symptom:** a "histidine-kinase expansion" that is one over-called domain. **Fix:** open `.inference.tsv`; prefer InterPro/Pfam architecture over free-text product when they conflict.

### Single-mode calling on a metagenome or tiny replicon
**Trigger:** Prodigal single mode on a MAG/plasmid/phage. **Mechanism:** self-training needs a full single-organism genome. **Symptom:** poor calls on short/mixed input. **Fix:** Bakta `--meta` (anonymous mode).

### Missing signal peptides on macOS
**Trigger:** Bakta on macOS. **Mechanism:** DeepSig was dropped from the default mac conda env (~v1.9.4). **Symptom:** silently absent signal-peptide calls. **Fix:** verify the run log; run on Linux if SPs matter.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Coding density ~88-90% (band 85-93%) | bacterial genome norm | <85% = wrong table/fragmentation/decay; >93% = over-calling |
| ~850-1,000 genes/Mb (~1 gene/kb) | bacterial gene density | far higher = over-call; far lower = under-call/wrong table |
| Hypothetical 25-50% | annotation norm | ~0% = over-confident transfer; >60-70% = novel lineage or wrong tool/table |
| tRNA count ≥ ~20 (often 40-60) | one isoacceptor set minimum | far fewer = fragmented assembly broke tRNA regions |
| rRNA operons ~1-15 (e.g. ~7 in *E. coli*) | copy-number norm | zero/fractional in a "complete" genome = short-read repeat collapse |
| Prodigal minimum ~30 aa / 90 nt | caller default | shorter ORFs need a dedicated sORF caller + Ribo-seq |
| CheckM2 contamination ≤5%, completeness ≥90% | MIMAG-aligned | above/below -> annotation QC numbers uninterpretable, fix assembly first |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Genes split, low coding density | wrong translation table | `--translation-table 4/25`; confirm from GTDB-Tk |
| Many hypothetical proteins | novel organism / light DB | full Bakta DB; add InterProScan/eggNOG; normal if 25-50% |
| Low gene / tRNA count, missing rRNA | fragmented assembly | CheckM2; prefer long-read/complete assembly |
| RefSeq record differs from the local GFF | RefSeq is PGAP, GenBank keeps the submitter's annotation | report WP_ accession + locus_tag; expect divergence |
| Pseudogene spike | ONT homopolymer indels vs real decay vs phase variation | check homopolymer context; corroborate with short-read/HiFi |
| Submission rejected on locus tags | invented prefix | register the prefix via BioSample (3-12 alnum, starts with a letter) |

## References

- Schwengers O, et al. 2021. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom* 7:000685.
- Seemann T. 2014. Prokka: rapid prokaryotic genome annotation. *Bioinformatics* 30:2068-2069.
- Hyatt D, et al. 2010. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11:119.
- Larralde M. 2022. Pyrodigal: Python bindings and interface to Prodigal. *J Open Source Softw* 7:4296.
- Lomsadze A, et al. 2018. Modeling leaderless transcription and atypical genes results in more accurate gene prediction in prokaryotes (GeneMarkS-2). *Genome Res* 28:1079-1089.
- Tatusova T, et al. 2016. NCBI prokaryotic genome annotation pipeline. *Nucleic Acids Res* 44:6614-6624.
- Li W, et al. 2021. RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. *Nucleic Acids Res* 49:D1020-D1028.
- Tanizawa Y, et al. 2018. DFAST: a flexible prokaryotic genome annotation pipeline for faster genome publication. *Bioinformatics* 34:1037-1039.
- Chan PP, et al. 2021. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Res* 49:9077-9096.
- Chklovski A, et al. 2023. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 20:1203-1212.
- Dimonaco NJ, et al. 2022. No one tool to rule them all: prokaryotic gene prediction tool annotations are highly dependent on the organism of study. *Bioinformatics* 38:1198-1207.
- Tonkin-Hill G, et al. 2020. Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biol* 21:180.
- Salzberg SL. 2019. Next-generation genome annotation: we still struggle to get it right. *Genome Biol* 20:92.
- Cole ST, et al. 2001. Massive gene decay in the leprosy bacillus. *Nature* 409:1007-1011.

## Related Skills

- functional-annotation - Add GO/KEGG/Pfam to hypothetical proteins with eggNOG-mapper/InterProScan
- ncrna-annotation - Detailed ncRNA identification with Infernal/Rfam beyond the built-in callers
- annotation-qc - CheckM2 completeness/contamination gate and gene-set sanity metrics
- genome-assembly/assembly-qc - Assess assembly quality before annotation
- genome-intervals/gtf-gff-handling - Parse and manipulate GFF3 output
- comparative-genomics/pangenome-analysis - Uniform re-annotation before pangenome clustering
