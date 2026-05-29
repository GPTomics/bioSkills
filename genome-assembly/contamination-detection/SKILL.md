---
name: bio-genome-assembly-contamination-detection
description: Detects and removes contamination in genome assemblies via two disjoint workflows - foreign-sequence screening of a single-organism (eukaryote/isolate) assembly with NCBI FCS-GX (GenBank-submission-mandatory), FCS-adaptor, and BlobToolKit blob plots; and MAG/bin quality assessment with CheckM2 plus GUNC (chimerism) plus GTDB-Tk taxonomy, judged against MIMAG. Covers why CheckM2 alone is blind to disjoint-marker chimeras, the FCS-GX RAM wall, organelle/NUMT triage, strain heterogeneity, and the HGT-vs-contamination (tardigrade) trap. Use when screening an assembly for foreign contamination before GenBank submission, assessing MAG completeness/contamination/chimerism, deciding which contigs to remove, or distinguishing real HGT from contaminant contigs.
tool_type: cli
primary_tool: CheckM2
---

## Version Compatibility

Reference examples tested with: FCS-GX 0.5+, FCS-adaptor 0.5+, BlobToolKit 4.3+, CheckM2 1.0+, GUNC 1.0+, GTDB-Tk 2.4+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

Database release drives results more than the binary version. Record: the **FCS-GX GX database** release (it grows each release; ~470 GiB and rising), the **CheckM2 DIAMOND DB** version, the **GUNC reference** (proGenomes 2.1 default vs GTDB), and the **GTDB-Tk reference package release** (e.g. R220) - GTDB-Tk fails loudly if the DB release does not match the binary. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Contamination Detection

**"Is my assembly contaminated?"** -> First decide which of two unrelated questions is being asked - "which contigs are not my organism?" (foreign-sequence screen) or "how complete/clean/chimeric is this bin?" (MAG quality) - then run the matching toolset; crossing them returns confidently wrong numbers.
- CLI (foreign): `run_fcsadaptor.sh` then `fcs.py screen genome --fasta asm.fa --gx-db <db> --tax-id <taxid>`, `blobtools create/add/view`
- CLI (MAG): `checkm2 predict -i bins/ -o out/` AND `gunc run -d bins/ -o out/`, `gtdbtk classify_wf`

## The Single Most Important Modern Insight -- "Contamination" Is Two Disjoint Problems With Disjoint Toolsets

There is no single "contamination" measurement. There are two questions that share a word and almost nothing else, and applying one toolset to the other problem runs to completion and prints a plausible, meaningless number.

1. **Foreign sequence in a single-organism assembly** (a eukaryotic nuclear genome or a cultured isolate with bacterial/human/vector contigs in it). The question is *which contigs are not my target organism?* Answer = a **set of contigs to remove, trim, or split out**. Tools: **NCBI FCS-GX** + **FCS-adaptor** (now GenBank-submission-mandatory), **BlobToolKit** (the GC x coverage x taxonomy blob plot). There is no "contamination %".

2. **MAG/bin quality** (a bin recovered from a metagenome). The question is *how complete is this bin and how much is from other organisms?* Answer = a **percentage pair** (completeness %, contamination %) plus a chimerism verdict, judged against **MIMAG**. Tools: **CheckM2** + **GUNC** together, **GTDB-Tk** for taxonomy.

**The cardinal category error:** running CheckM2/CheckM on a eukaryotic nuclear assembly (its bacterial/archaeal marker sets are meaningless for a eukaryote - eukaryote completeness is BUSCO; see assembly-qc), or running FCS-GX on a single MAG and reading the flagged-length fraction as "the MIMAG contamination %". Both wrong applications complete silently and emit a number. The only defense is knowing which question each tool answers.

**The second load-bearing fact - CheckM2 and GUNC are blind to different things, so report the pair, never the % alone.** CheckM2 estimates contamination from **marker-gene redundancy** (a single-copy marker seen twice = contamination signal). A **chimera of two organisms with disjoint marker complements** - organism A contributed markers 1-60, organism B markers 61-120 - shows *no* duplicated markers, so CheckM2 reports low contamination while the bin is biological nonsense. GUNC catches exactly this, because it scores whether taxonomic signal is *consistent across contigs* (clade separation score), not marker counts. Neither is a superset of the other.

## Tool Taxonomy

| Tool | Citation | Role | Problem |
|------|----------|------|---------|
| FCS-adaptor | Astashyn 2024 *Genome Biol* | adaptor/vector screen (VecScreen successor); tiny DB, trivial RAM | foreign (run FIRST) |
| FCS-GX | Astashyn 2024 *Genome Biol* | cross-taxon foreign-sequence screen vs a declared tax-id; GenBank-mandatory | foreign |
| BlobToolKit | Challis 2020 *G3* | GC x coverage x taxonomy blob plot for visual triage; the maintained successor to BlobTools | foreign |
| tiara | Karlicki 2022 *Bioinformatics* | deep-learning euk/prok/organelle/plastid/mito contig classifier (alignment-free) | foreign / organelle partitioning |
| CheckM2 | Chklovski 2023 *Nat Methods* | ML completeness + contamination (lineage-agnostic; handles novel/reduced lineages) | MAG (default) |
| CheckM (legacy) | Parks 2015 *Genome Res* | lineage-specific collocated markers on a placement tree; reports Strain heterogeneity | MAG (legacy; slow, ~40 GB RAM) |
| GUNC | Orakov 2021 *Genome Biol* | taxonomic-chimerism detection via clade separation score (CSS) | MAG (run WITH CheckM2) |
| GTDB-Tk | Chaumeil 2020 *Bioinformatics* | GTDB taxonomic classification of bacterial/archaeal genomes | MAG taxonomy |

BlobTools (Laetsch & Blaxter 2017 *F1000Res*) and CheckM are the legacy predecessors; BlobToolKit and CheckM2 are the maintained defaults. CAT/BAT (contig taxonomy) and MAGpurify (reference-free MAG cleaning) exist but are secondary; RefineM is deprecated by its author (use GUNC + manual curation, or MDMcleaner for SAGs/dark-matter lineages).

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Eukaryotic nuclear assembly headed for GenBank | FCS-adaptor -> FCS-GX (mandatory) -> BlobToolKit for triage | submission bounces without a clean FCS screen |
| Cultured bacterial/archaeal isolate | FCS-GX + FCS-adaptor; CheckM2 as a completeness sanity check | an isolate should be ~100% complete, low contam |
| MAG(s) binned from a metagenome | CheckM2 AND GUNC together; score vs MIMAG; GTDB-Tk for taxonomy | the % pair plus the orthogonal chimerism verdict |
| Unbinned metagenome | bin first (-> metagenome-assembly), THEN CheckM2 + GUNC per bin | CheckM2/GUNC consume bins, do not make them |
| Eukaryotic contigs mixed into a metagenome | tiara / Whokaryote / EukRep to partition euk from prok | CheckM2 markers are prokaryote-only |
| Organellar (mito/chloroplast) contigs in a nuclear assembly | tiara organelle classes + coverage spike -> separate, do NOT delete | organelles are real biology; submit as own record |
| Foreign-looking gene embedded in a host scaffold | investigate integration (-> comparative-genomics) before removing | could be real HGT/endosymbiont, not contamination |
| Read-level host removal before assembly | -> read-qc/contamination-screening, metagenomics/kraken-classification | this skill screens assembled sequence, not reads |
| ANI/species placement of the cleaned genome | -> comparative-genomics/genome-distance-and-species-delineation | taxonomy after decontamination |

## FCS-adaptor (Run First - Cheap and Unambiguous)

```bash
# Adaptor/vector screen; tiny DB, trivial RAM. --euk or --prok for the lineage.
run_fcsadaptor.sh --fasta-input assembly.fa.gz --output-dir ./adaptor_out --euk
```

Adaptor and vector hits are unambiguous - always trim or exclude. Running adaptor first removes the obvious junk and shrinks the input before the expensive GX pass.

## FCS-GX (The GenBank-Mandatory Foreign Screen)

```bash
# Screen against the GX database; --tax-id is the NCBI taxid of the SOURCE organism.
python3 ./fcs.py screen genome --fasta assembly.fa.gz --out-dir ./gx_out/ \
    --gx-db "$GXDB_LOC/gxdb" --tax-id 9606

# Apply ONLY the auto-clean actions (EXCLUDE/TRIM/FIX) to produce a cleaned FASTA.
zcat assembly.fa.gz | python3 ./fcs.py clean genome \
    --action-report ./gx_out/assembly.fa.9606.fcs_gx_report.txt \
    --output clean.fasta --contam-fasta-out contam.fasta
```

The action report (`<name>.<taxid>.fcs_gx_report.txt`, where `<name>` keeps the input name minus only `.gz` - e.g. `assembly.fa.gz` -> `assembly.fa.9606.fcs_gx_report.txt`) assigns each flagged region an action: **EXCLUDE** (drop the whole sequence), **TRIM** (remove a contaminated end), **FIX** (hard-mask an internal span), or **REVIEW** (flagged but NOT auto-cleaned - manual inspection required). A separate **INFO** action specifically marks sequence known to be integrated into host genomes (e.g. endosymbiont insertions). `clean genome` acts on EXCLUDE/TRIM/FIX only; it does **not** touch REVIEW or INFO. Those non-auto-cleaned tiers exist precisely because FCS-GX can flag legitimate HGT/endosymbiont sequence as contaminant - never bulk-delete every flagged contig without reading the report.

**The RAM wall (an operational, cloud-forcing fact).** The GX database is ~470 GiB on disk and the documented sweet spot is a **~512 GiB-RAM** host. Underprovisioned, the run does not fail - it crawls (minutes become days when the DB spills out of RAM). Best practice copies the DB into a tmpfs RAM disk; most labs rent a high-memory cloud VM for the screen. The first question before advising FCS-GX is "is there ~1/2 TB RAM or a cloud budget?".

## BlobToolKit (The Blob Plot - Visual Triage)

```bash
blobtools create --fasta assembly.fasta ./BlobDir
blobtools add --hits diamond.out --taxrule bestsumorder --taxdump /path/taxdump \
    --cov mapping.bam ./BlobDir
blobtools view --remote ./BlobDir   # interactive viewer; GC(x) vs coverage(y), sized by length, colored by taxonomy
```

Each contig is plotted by **GC fraction (x)** vs **coverage (y)**, sized by length, colored by best-hit taxonomy. The target organism forms one tight GC x coverage cloud; a free-living contaminant grows independently and forms its own cloud at a *different* GC *and* a *different* coverage. Inputs: a hit file (BLAST/DIAMOND vs nt/UniProt), a coverage BAM, and an NCBI taxdump. `blobtools filter` can extract or drop by taxon - but see the failure mode below before doing so.

## Merge CheckM2 and GUNC

**Goal:** Produce the joint CheckM2 x GUNC table that is the field standard for MAG QC, so a chimera invisible to CheckM2 is not reported as clean.

**Approach:** Load both reports, join on genome name, and apply MIMAG thresholds AND the GUNC pass flag together; never gate on contamination % alone.

```python
import pandas as pd

CONTAM_HQ = 5      # MIMAG high-quality: contamination < 5% (Bowers 2017); above this, gene-content inference is unreliable
COMPLETE_HQ = 90   # MIMAG high-quality: completeness > 90%
RRS_TRUST = 0.5    # GUNC pass is trustworthy only when reference_representation_score > 0.5; below it 'pass' means 'can't tell'

checkm = pd.read_csv('checkm2_out/quality_report.tsv', sep='\t')
gunc = pd.read_csv('gunc_out/GUNC.progenomes_2.1.maxCSS_level.tsv', sep='\t')
merged = checkm.merge(gunc, left_on='Name', right_on='genome', how='left')

merged['gunc_trustworthy'] = merged['reference_representation_score'] > RRS_TRUST
merged['high_quality'] = ((merged['Completeness'] > COMPLETE_HQ) &
                          (merged['Contamination'] < CONTAM_HQ) &
                          (merged['pass.GUNC'] == True) &
                          merged['gunc_trustworthy'])
merged.to_csv('combined_qc.tsv', sep='\t', index=False)
```

CheckM2 also runs as `checkm2 predict -i bins/ -o checkm2_out --threads 16 -x fa` (downloads its own DIAMOND DB), and GUNC as `gunc run -d bins/ -o gunc_out -t 16 -e .fa` against a downloaded reference (`gunc download_db`). GTDB-Tk taxonomy is `gtdbtk classify_wf --genome_dir bins/ --out_dir gtdbtk_out -x fa --cpus 16`.

## GUNC and the RRS Trap

`pass.GUNC = TRUE` at **CSS <= 0.45** is not a clean bill of health when `reference_representation_score` (RRS) is low. Low RRS means the genome is barely represented in GUNC's reference DB - there is not enough signal to *detect* chimerism, so "pass" degenerates into "can't tell, not clean". For genuinely novel lineages (microbial dark matter), a GUNC pass is weak evidence; trust it only when RRS > 0.5, and otherwise lean on manual contig inspection or MDMcleaner. Read the configurations: CheckM2-high-contam + GUNC-pass = duplicative contamination within one lineage; CheckM2-low-contam + GUNC-fail = the dangerous disjoint-marker chimera; both pass with high RRS = the only "clean" verdict.

## The Tardigrade Lesson (HGT vs Contamination)

Boothby 2015 (*PNAS* 112:15976) reported ~17% of tardigrade genes arrived by horizontal gene transfer; Koutsovoulos 2016 (*PNAS* 113:5053) showed almost all of it was undetected bacterial contamination in the assembly, dropping real HGT to ~1-2%. The original paper was **corrected, not retracted** - it stands as a cautionary monument. A foreign-looking gene has two explanations the *sequence alone cannot distinguish*: a contaminant contig, or a real HGT/endosymbiont gene integrated into the host genome. **The tell is physical integration, not taxonomy.** A real HGT gene sits *on a host contig* (host-gene flanks, host GC, host coverage, host introns); a contaminant sits *on its own contig* (foreign GC, its own coverage cloud, prokaryotic gene structure). Deleting everything that looks bacterial manufactures the opposite of Boothby's error - erasing real biology. Require integration evidence (host flanks, host-typical GC/coverage/introns, long reads spanning the host->foreign junction) before keeping or stripping.

## Organelles and NUMTs (A Third Category)

Mitochondrial and chloroplast sequence is legitimate biology that does not belong in the *nuclear* assembly record - it needs **separation, not deletion**. Deleting it loses the organelle genome (often the most-cited part of a non-model genome paper); leaving it inflates assembly size and creates fake duplicated content. Identify organellar contigs by tiara's plastid/mito classes plus a massive coverage spike and circular topology; assemble/submit them as their own record. The subtle trap: **NUMTs/NUPTs** (nuclear-integrated organellar fragments) are *real nuclear sequence that looks organellar* - do not strip them. Coverage discriminates: free organelle = huge coverage spike; NUMT = nuclear-level coverage.

## Per-Method Failure Modes

### CheckM2 on a eukaryotic nuclear assembly
**Trigger:** running CheckM/CheckM2 on a vertebrate/plant/insect genome. **Mechanism:** the marker sets are bacterial/archaeal single-copy genes; a eukaryote has none of the relevant context. **Symptom:** a plausible Completeness/Contamination pair that is pure noise. **Fix:** use BUSCO for eukaryote completeness (-> assembly-qc); CheckM2 is for prokaryotes/MAGs only.

### CheckM2 contamination read as the whole story (the chimera blind spot)
**Trigger:** "CheckM2 says 3%, the MAG is clean." **Mechanism:** CheckM2 counts marker redundancy and is blind to a chimera of two organisms with disjoint marker sets. **Symptom:** a 50/50 chimeric bin reported as low-contamination. **Fix:** run GUNC too; report the pair; trust a GUNC pass only when RRS > 0.5.

### FCS-GX % read as MIMAG contamination
**Trigger:** running FCS-GX on a single MAG and citing flagged-length as the contamination %. **Mechanism:** FCS-GX answers "foreign vs declared tax-id", not the intra-domain marker-redundancy MIMAG cares about. **Symptom:** a "contamination %" that means something different from CheckM2's. **Fix:** use CheckM2 for the MIMAG %; FCS-GX for foreign-sequence screening of single-organism assemblies.

### Auto-stripping every FCS-GX/blob foreign hit
**Trigger:** running `clean genome` and resubmitting, or `blobtools filter` on a whole taxon cloud, without reading the report. **Mechanism:** the screen flags legitimate HGT/endosymbiont/host genes with conserved-domain hits. **Symptom:** a real biological finding (HGT, symbiont) deleted from the assembly. **Fix:** treat EXCLUDE on near-complete-bacterial-genome contigs as real contamination; treat foreign flags on host-integrated, host-coverage, host-flanked sequence as REVIEW-and-verify.

### Trusting blob-plot taxonomy color over cluster geometry
**Trigger:** filtering out contigs by their BLAST-hit color. **Mechanism:** best-hit taxonomy is the noisiest axis - conserved domains (ribosomal proteins, HSP70) hit foreign taxa, and no-hit != contaminant. **Symptom:** real host genes deleted; large no-hit fractions panic-removed. **Fix:** decide on cluster geometry - off the main cloud on GC *and* coverage *and* consistent foreign taxonomy; coverage is the strongest single signal.

### Strain heterogeneity misread as contamination
**Trigger:** high CheckM contamination on a MAG with co-binned conspecific strains. **Mechanism:** near-identical single-copy markers from multiple strains register as duplicated. **Symptom:** inflated contamination % for what is one species' population. **Fix:** read CheckM's Strain heterogeneity column; high strain-heterogeneity = look closer (not ignore); GUNC + the strain flag disambiguate foreign-mixing from strain-mixing from duplication.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| MIMAG high-quality MAG: completeness > 90% AND contamination < 5% | Bowers 2017 *Nat Biotechnol* | above 5% contam, gene-content/metabolic inference is unreliable |
| MIMAG medium-quality: completeness >= 50% AND contamination < 10% | Bowers 2017 *Nat Biotechnol* | community-standard MQ floor |
| MIMAG HQ also requires 5S/16S/23S rRNA + >= 18 tRNAs | Bowers 2017 *Nat Biotechnol* | binners systematically lose rRNA; the most common reason a >90%/<5% MAG is only MQ |
| GUNC pass: CSS <= 0.45 | Orakov 2021 *Genome Biol* | benchmarked chimera/non-chimera cutoff; orthogonal to the % pair |
| GUNC pass trustworthy only when RRS > 0.5 | Orakov 2021 *Genome Biol* / GUNC docs | low RRS = too poorly represented to judge; pass means "can't tell" |
| Blob-plot contaminant call needs GC AND coverage AND taxonomy agreement | Laetsch & Blaxter 2017; field practice | any single axis alone is a false-positive generator |
| FCS-GX GX DB ~470 GiB, ~512 GiB RAM host | NCBI FCS wiki (grows per release) | underprovisioned = crawls, not fails; cloud-forcing |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| CheckM2 prints a Completeness/Contamination on a eukaryote | wrong tool for the organism | BUSCO for eukaryotes; CheckM2 is prokaryote/MAG only |
| "Clean" MAG (3% contam) is actually a chimera | CheckM2 blind to disjoint-marker chimeras | run GUNC; report the pair |
| GenBank submission bounced | FCS-GX/FCS-adaptor not run before submission | screen with FCS-adaptor then FCS-GX first |
| FCS-GX takes days | GX DB spilled out of RAM | provision ~512 GiB / tmpfs RAM disk / cloud VM |
| GTDB-Tk crashes on startup | DB release does not match the binary | install the reference package matching the GTDB-Tk version |
| Real HGT/symbiont gene deleted | auto-cleaned or taxon-filtered without review | check physical integration (host flanks/GC/coverage) before removal |
| Organelle genome lost | deleted as "contamination" | separate organellar contigs (coverage spike); submit as own record |

## References

- Astashyn A, et al. 2024. Rapid and sensitive detection of genome contamination at scale with FCS-GX. *Genome Biol* 25:60.
- Chklovski A, et al. 2023. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 20:1203-1212.
- Parks DH, et al. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Res* 25:1043-1055.
- Orakov A, et al. 2021. GUNC: detection of chimerism and contamination in prokaryotic genomes. *Genome Biol* 22:178.
- Bowers RM, et al. 2017. Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nat Biotechnol* 35:725-731.
- Chaumeil PA, et al. 2020. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics* 36:1925-1927.
- Challis R, et al. 2020. BlobToolKit - interactive quality assessment of genome assemblies. *G3 (Bethesda)* 10:1361-1374.
- Laetsch DR, Blaxter ML. 2017. BlobTools: interrogation of genome assemblies. *F1000Research* 6:1287.
- Karlicki M, Antonowicz S, Karnkowska A. 2022. Tiara: deep learning-based classification system for eukaryotic sequences. *Bioinformatics* 38:344-350.
- Boothby TC, et al. 2015. Evidence for extensive horizontal gene transfer from the draft genome of a tardigrade. *PNAS* 112:15976-15981.
- Koutsovoulos G, et al. 2016. No evidence for extensive horizontal gene transfer in the genome of the tardigrade Hypsibius dujardini. *PNAS* 113:5053-5058.

## Related Skills

- assembly-qc - BUSCO completeness for eukaryotes (not CheckM2) and the QC handoff
- metagenome-assembly - Binning produces the bins this skill scores with CheckM2 + GUNC
- hifi-assembly - False duplications inflate apparent content; distinct from contamination
- long-read-assembly - Produces the contigs screened here for foreign sequence
- metagenomics/kraken-classification - Read/contig taxonomic classification and pre-assembly host screening
- comparative-genomics/genome-distance-and-species-delineation - ANI/species placement after decontamination
- workflows/genome-assembly-pipeline - End-to-end assemble -> QC -> decontaminate
