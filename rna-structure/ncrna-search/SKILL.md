---
name: bio-rna-structure-ncrna-search
description: Searches for non-coding RNA homologs and classifies RNA families with Infernal covariance models against Rfam, scoring sequence AND secondary-structure conservation jointly. Use when deciding whether a covariance model is the right tool versus BLAST/nhmmer (structured ncRNA versus lncRNA or mature miRNA); choosing the Rfam gathering threshold over a flat E-value; resolving clan overlaps; building and calibrating a custom CM from a structure-annotated alignment; or preferring a family-specialized tool (tRNAscan-SE, barrnap) over a generic Rfam scan.
tool_type: cli
primary_tool: Infernal
---

## Version Compatibility

Reference examples tested with: Infernal 1.1.4+, BioPython 1.83+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# ncRNA Search

**"Search my sequences for known non-coding RNA families"** -> Score candidates against Rfam covariance models that capture both sequence and consensus secondary structure, or build a custom model for a novel family.
- CLI: `cmscan` for querying sequences against the Rfam CM database
- CLI: `cmsearch` for one CM against a sequence database
- CLI: `cmbuild` + `cmcalibrate` + `cmpress` for a custom covariance model

## The governing principle: a covariance model is only worth it for STRUCTURED RNA

A covariance model (CM) is a profile stochastic context-free grammar that models a family's consensus secondary structure AND primary sequence jointly. Its power source is covariation: a base pair is scored by a joint distribution over pair states, so a compensatory double mutation (C-G to G-C, or a G-U wobble) that PRESERVES the pair scores well even though both positions changed -- a BLAST or profile-HMM sees two mismatches and loses the signal. This is why Infernal detects remote homologs of structured RNAs (tRNA, rRNA, riboswitches, SRP, RNase P, ribozymes, snoRNA) past the sequence twilight zone.

The decisive corollary: a CM offers NO advantage when there is little conserved secondary structure to exploit. Many lncRNAs (R-scape finds no significant covariation in HOTAIR/Xist/SRA), mature miRNAs (~22 nt, no base-pairing in the mature strand), and primary-sequence-only motifs gain nothing from a CM -- it is slow and calibration-heavy, and a profile-HMM (nhmmer) or BLASTN is the correct, faster tool. A CM built from a structure-FREE alignment (no real `#=GC SS_cons` pairs) collapses to an HMM and buys nothing but runtime.

| Query / target property | Right tool | Why |
|---|---|---|
| Structured ncRNA, remote homology (tRNA, rRNA, riboswitch, ribozyme, SRP, snoRNA) | Infernal CM (Rfam) | covariation recovers pairs across sequence divergence |
| Structured RNA with a family-specialized tool | the specialist (table below) | tuned models + biology logic beat a generic scan |
| Close homolog, high identity, any RNA | BLASTN / nhmmer | sequence signal suffices, 100-1000x faster |
| lncRNA, mature miRNA, sequence-only motif | nhmmer / BLASTN | no conserved structure to exploit -> CM is overhead |
| Novel structured RNA, no Rfam family | build a custom CM, AFTER validating structure with R-scape | a CM is only as good as its SS_cons |

## Infernal toolchain

- `cmbuild model.cm aln.sto` -- build a CM from a structure-annotated Stockholm alignment; REQUIRES a `#=GC SS_cons` line.
- `cmcalibrate model.cm` -- fit E-value statistics; SLOW (minutes to hours) but MANDATORY before any E-value is meaningful. Without it, search still runs but only bit-score thresholding is valid.
- `cmpress model.cm` -- build the binary index `cmscan` requires (cmsearch does not need it).
- `cmsearch CM seqdb` -- one CM vs a sequence database (e.g. one family vs a genome).
- `cmscan CMdb seqs` -- query sequences vs a CM database (e.g. all of Rfam.cm); this is the genome/transcript ncRNA-annotation layout.
- `cmalign CM seqs` -- align/fold hits back to the CM consensus to recover each hit's secondary structure.

Rfam.cm ships PRE-CALIBRATED: run `cmpress` on it, but do NOT `cmcalibrate` it. Calibration is only for locally built custom models, and must be re-run after every rebuild.

## E-values depend on database size; gathering thresholds do not

A CM E-value scales linearly with the searched database size (Z): the SAME hit gets a different E-value depending on what was searched. `cmsearch` Z defaults to the target sequence-DB size counted on both strands (total residues x2); `cmscan` Z defaults to (query length x2 x number of models). So a cmscan E-value and a cmsearch E-value for the same locus are NOT comparable, and a flat `-E 1e-5` is exactly what curated thresholds replace. Fix Z explicitly with `-Z <Mb>` for reproducible cross-run comparison.

Bit scores are DB-size-INDEPENDENT, which is precisely why Rfam stores its curated cutoffs as bit scores:

| Flag | Threshold (bit-score cutoff stored in the CM) | When |
|---|---|---|
| `--cut_ga` | GA (gathering): the curated family-membership cutoff | the correct DEFAULT for Rfam annotation |
| `--cut_tc` | TC (trusted): lowest score of any known true positive | most conservative |
| `--cut_nc` | NC (noise): highest score of a known false positive | most permissive, exploratory |
| `-E` / `--incE` | flat E-value (reporting / inclusion) | custom CM, or a non-Rfam DB with no curated GA |
| `-T` / `--incT` | flat bit score | custom CM with no calibration, or DB-size-robust cut |

`--cut_ga` beats a flat E-value for Rfam because each family has a different signal-to-noise profile (a 70 nt tRNA vs a 2900 nt rRNA vs a short riboswitch); one flat cutoff over- or under-calls per family, and reintroduces the DB-size dependence GA was designed to avoid.

## The canonical Rfam annotation command

```bash
# One-time setup
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz && gunzip Rfam.cm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
cmpress Rfam.cm   # NOT cmcalibrate -- Rfam.cm is pre-calibrated

# Annotate a genome. -Z = 2 x genome size in Mb keeps E-values reproducible; --rfam is the
# large-DB strict filter; --nohmmonly forces CM mode so GA cutoffs stay valid for every model.
cmscan -Z 100 --cut_ga --rfam --nohmmonly --fmt 2 --clanin Rfam.clanin \
    --tblout genome.tblout Rfam.cm genome.fa > genome.cmscan

# Clan deoverlapping: --fmt 2 adds the 'olp' column; the documented Rfam filter drops hits
# marked '=' (dominated by a higher-scoring clanmate), keeping '^' (best of an overlap) and '*' (no overlap).
grep -v ' = ' genome.tblout > genome.deoverlapped.tblout
```

`--clanin` with `--fmt 2` plus the `grep -v ' = '` post-filter is the documented deoverlap path; `--oclan` is a valid in-tool alternative but not what the modern Rfam pipeline uses, so do not treat it as mandatory.

## Reading --fmt 2 output (the column shift that breaks fmt-1 parsers)

`--fmt 2` prepends an `idx` column and inserts a `clan name` column versus the default `--fmt 1`, so every downstream field index shifts -- a parser written for fmt 1 silently reads the wrong columns. Verified 0-based fmt-2 cmscan indices: idx 0, target/family name 1, accession 2, query/seq name 3, query accession 4, clan 5, mdl type 6, mdl_from 7, mdl_to 8, seq_from 9, seq_to 10, strand 11, trunc 12, pass 13, gc 14, bias 15, score 16, E-value 17, inc 18, olp 19. In cmscan the model name is column 1 and the sequence name is column 3; in cmsearch they are reversed -- a parser must branch on which tool produced the file.

Interpretive columns: `trunc` is `no`, `5'`, `3'`, or `5'&3'` for hits running off a contig end (often REAL incomplete genes worth keeping, not a discard flag; `--anytrunc` allows truncation at any internal position (E-values become less accurate), `--notrunc` disables it entirely); `bias` is the composition correction already subtracted (a large bias flags a low-complexity hit); mdl_from/mdl_to are consensus coordinates, so a small model span is a partial match even when `trunc` is `no`.

## Recovering the structure of a hit (the CM payoff over BLAST)

A significant CM hit yields not just a family label but an implied secondary structure -- the distinctive payoff a BLAST hit cannot give. Fold the hit sequences back to the model with `cmalign` to map the consensus structure onto them.

```bash
# Pull the family CM from Rfam, then align hits to recover their consensus structure
cmfetch Rfam.cm RF00005 > tRNA.cm
cmalign --outformat Pfam tRNA.cm hits.fa > hits.sto   # Stockholm with #=GC SS_cons per hit
```

The `#=GC SS_cons` line in the output is the structural hypothesis for each hit -- the prior to feed to probing (structure-probing) or to validate with covariation (covariation-analysis).

## A hit is a family assignment, not a functional call

A significant CM hit means the locus has sequence + structure consistent with the family; it does NOT prove the RNA is expressed, processed, or functional. Pseudogenes (tRNA-derived SINEs, rRNA pseudogenes) score well. tRNAscan-SE's high-confidence-set logic exists precisely to separate likely-functional tRNAs from numerous genomic pseudogenes. Frame a CM hit as a family assignment plus a structural hypothesis; confirm function with orthogonal evidence (expression, conservation, synteny, probing).

## Specialized tool beats a generic Rfam scan (when one exists)

| RNA class | Use this, not a generic Rfam scan | Why |
|---|---|---|
| tRNA | tRNAscan-SE 2.0 | tRNA-specific isotype/anticodon models, pseudogene filtering, high-confidence set |
| rRNA (5S/16S/18S/23S/28S) | barrnap or RNAmmer | per-kingdom HMM models tuned for rRNA; fast |
| C/D-box snoRNA | snoscan | models the guide-target duplex + box C/D a generic CM cannot |
| H/ACA + C/D snoRNA (ab initio) | snoReport 2.0 | RNA-fold + SVM on box/structure features |
| miRNA | miRBase lookup / miRDeep2 | CMs are poor on mature miRNA; precursors need read support |

Use the generic Rfam cmscan for broad "what ncRNA families are in here" sweeps and for classes without a dedicated tool.

## Building a custom CM well

```bash
# The Stockholm alignment MUST carry a #=GC SS_cons line in WUSS notation: <> or () = nested pairs,
# [] and {} = additional pseudoknot layers, . = unpaired (a structure-free alignment yields only an
# HMM-equivalent model). Validate the SS_cons covariation with R-scape FIRST (covariation-analysis).
cmbuild -n MYFAM myfam.cm alignment.sto
cmcalibrate --cpu 8 myfam.cm     # required for E-values; re-run after every rebuild
cmpress myfam.cm
cmsearch --cpu 8 -T 30 --tblout hits.tbl myfam.cm target.fa > hits.out
```

If `cmcalibrate` is skipped, thresholding MUST use bit score (`-T <bits>`), not E-value. The iterative search-align-rebuild loop (`cmsearch -A new_hits.sto`) grows a family but risks homology overextension and model drift -- gate each round on score and retained covariation, and recalibrate after every rebuild.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Parsed `score`/`evalue` are nonsense numbers | parsing `--fmt 2` output with fmt-1 indices | use fmt-2 indices (score 16, E-value 17), or run `--fmt 1` |
| `Error: failed to open ... .i1m` (cmscan) | Rfam.cm not pressed | `cmpress Rfam.cm` |
| E-values look meaningful on a custom CM but are not | `cmcalibrate` was skipped | calibrate, or threshold on bit score `-T` |
| Recalibrating Rfam.cm takes hours | Rfam.cm is already calibrated | press it, never calibrate it |
| Same hit, different E-value across runs | E-value scales with database size Z | use `--cut_ga`, or fix `-Z <Mb>` |
| Redundant overlapping hits from related families | clan overlap not resolved | `--fmt 2 --clanin` then `grep -v ' = '` |
| A CM search on a lncRNA finds nothing useful | no conserved structure to exploit | use nhmmer/BLASTN instead of a CM |
| tRNA search misses or over-calls genes | generic Rfam tRNA model lacks tRNA logic | use tRNAscan-SE 2.0 |

## Related Skills

- secondary-structure-prediction - Predict structure for novel ncRNA candidates with no Rfam hit
- covariation-analysis - Validate a custom CM's SS_cons with R-scape before building
- structure-probing - Experimental reactivities to corroborate a CM's consensus structure
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation pipelines
- alignment/msa-statistics - Evaluate alignment quality before CM building
- database-access/entrez-fetch - Fetch Rfam/RNAcentral records

## References

- Eddy SR, Durbin R. 1994. RNA sequence analysis using covariance models. Nucleic Acids Res 22(11):2079-2088. doi:10.1093/nar/22.11.2079
- Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29(22):2933-2935. doi:10.1093/bioinformatics/btt509
- Griffiths-Jones S, Bateman A, Marshall M, Khanna A, Eddy SR. 2003. Rfam: an RNA family database. Nucleic Acids Res 31(1):439-441. doi:10.1093/nar/gkg006
- Kalvari I, Nawrocki EP, Ontiveros-Palacios N, Argasinska J, Lamkiewicz K, Marz M, Griffiths-Jones S, Toffano-Nioche C, Gautheret D, Weinberg Z, Rivas E, Eddy SR, Finn RD, Bateman A, Petrov AI. 2021. Rfam 14: expanded coverage of metagenomic, viral and microRNA families. Nucleic Acids Res 49(D1):D192-D200. doi:10.1093/nar/gkaa1047
- Chan PP, Lin BY, Mak AJ, Lowe TM. 2021. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. Nucleic Acids Res 49(16):9077-9096. doi:10.1093/nar/gkab688
- Lagesen K, Hallin P, Rodland EA, Staerfeldt HH, Rognes T, Ussery DW. 2007. RNAmmer: consistent and rapid annotation of ribosomal RNA genes. Nucleic Acids Res 35(9):3100-3108. doi:10.1093/nar/gkm160
- Lowe TM, Eddy SR. 1999. A computational screen for methylation guide snoRNAs in yeast. Science 283(5405):1168-1171. doi:10.1126/science.283.5405.1168
- Rivas E, Clements J, Eddy SR. 2017. A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs. Nat Methods 14(1):45-48. doi:10.1038/nmeth.4066
