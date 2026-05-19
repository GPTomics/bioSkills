---
name: bio-clip-seq-crosslink-site-detection
description: Detect single-nucleotide crosslink (CL) sites in CLIP-seq data using truncation patterns (iCLIP/eCLIP CITS), crosslink-induced mutations (HITS-CLIP CIMS deletions, PAR-CLIP T-to-C), or HMM/kernel-density methods (PureCLIP, PARalyzer, CTK). Use when single-nucleotide resolution is required for motif registration (mCross), allele-specific binding (BEAPR), variant-effect prediction, or comparing crosslink chemistry across CLIP variants.
tool_type: cli
primary_tool: PureCLIP
---

## Version Compatibility

Reference examples tested with: PureCLIP 1.3.1+, CTK 1.1.4+, PARalyzer 1.5+, wavClusteR 2.34+, pyCRAC 1.5+, samtools 1.19+, bedtools 2.31+, pysam 0.22+, R 4.3+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws unexpected errors, introspect the installed tool and adapt the example to match the actual CLI rather than retrying.

# CLIP-seq Crosslink-Site Detection

**"Detect single-nucleotide crosslink sites in my CLIP data"** -> Identify the exact base where the protein-RNA UV adduct caused the reverse transcriptase to stop (truncation, in iCLIP/eCLIP), to read through with a mutation (deletion in HITS-CLIP, T->C in PAR-CLIP), or to leave a multi-base signature (PARalyzer kernel density for PAR-CLIP). Single-nucleotide resolution is the foundation of motif registration (mCross), allele-specific binding (BEAPR/ASPRIN), variant-effect prediction, and the most rigorous comparisons across CLIP variants. The detection chemistry differs by CLIP type: iCLIP/eCLIP read 5' end is one nucleotide downstream of the crosslink (CITS); PAR-CLIP reads contain T->C transitions at the crosslink (CIMS substitution); HITS-CLIP reads contain deletions at the crosslink (CIMS deletion).

- CLI (HMM, all CLIP variants): `pureclip -i dedup.bam -bai dedup.bam.bai -g genome.fa -ibam sminput.bam -ibai sminput.bam.bai -o crosslinks.bed -or regions.bed -nt 8 -dm 8`
- CLI (CTK CITS truncations, iCLIP/eCLIP): `parseAlignment.pl --map-qual 1 --min-len 18 --mutation-file mut.txt dedup.bam dedup.bed` then `tag2cluster.pl ... -cs5 5 -m 1` for truncation-cluster
- CLI (CTK CIMS deletions, HITS-CLIP): `getMutationType.pl dedup.bed mut.txt -type del` then `CIMS.pl dedup.bed mut.txt -big -c -p 0.01 cims.bed`
- CLI (CTK CIMS T->C, PAR-CLIP): `getMutationType.pl dedup.bed mut.txt -type sub -nuc t -mut c` then `CIMS.pl dedup.bed t2c.mut -p 0.001 t2c_cims.bed`
- CLI (PAR-CLIP kernel density): `PARalyzer params.ini` (parameters file defines read length, min reads per cluster, mutation rate threshold)
- CLI (PAR-CLIP wavClusteR R): `wavClusteR::filterClusters(cl, snps=NULL, filterFC=FALSE)` after wavelet clustering

Single-nt CL sites are the input to mCross motif registration (see clip-seq/clip-motif-analysis), to BEAPR/ASPRIN allele-specific binding analyses, and to RBPNet/DeepRiPe deep-learning models (see clip-seq/clip-deep-learning). They are NOT a replacement for broad peak calls; the two outputs are complementary. A common downstream error is passing a CLIPper peak BED to mCross instead of a PureCLIP single-nt BED - mCross requires the single-nt resolution to register motif position relative to the crosslink offset.

## Crosslink Chemistry by CLIP Variant

The detection method must match the underlying chemistry of how the reverse transcriptase encountered the protein-RNA adduct:

| CLIP variant | RT behavior at CL | Detection signature | Sensitivity (CL captured per read) | Sequence bias |
|--------------|-------------------|---------------------|-------------------------------------|---------------|
| iCLIP / iCLIP2 / iCLIP3 / eCLIP / irCLIP / FLASH | Truncates at adduct | 5' end of cDNA = CL site - 1 | ~80% of reads truncate, 20% read through with deletion | Strong U bias (~60-80% of CL events at U) |
| HITS-CLIP | Reads through with deletion (3-7% of crosslinked positions) | Single-base deletion in read | ~3-7% of reads contain deletion at CL site | U bias plus sequence-specific deletion rate |
| PAR-CLIP (4SU) | Reads through with T->C (20-50% of T positions in crosslinked reads) | T->C transition (or G->A on reverse strand) | Per-T rate 20-50%; per-read 1-5 conversions | Restricted to T positions; depends on 4SU incorporation rate |
| PAR-CLIP (6SG) | Reads through with G->A | G->A transition | Lower rate than 4SU T->C | Restricted to G positions |
| miCLIP / miCLIP2 (m6A) | Primarily truncation at the m6A site (RT stops at antibody-trapped m6A); secondary C->T transition observed in a subset of reads at m6A | 5' end (CL -1); C->T rate is a feature used by m6Aboost ML in addition to truncation, not a stand-alone primary signal | Mixed; m6Aboost ML integrates truncation + sequence context to score | DRACH motif context required |
| STAMP | C->T editing on target (not crosslink) | C->T editing in mRNA reads | NA (editing-based; not crosslink-based) | N/A |
| TRIBE | A->I editing on target | A->G in cDNA (I read as G) | NA | Adjacent to ADAR consensus |

## Algorithmic Taxonomy

| Tool | CLIP variant | Detection signal | Statistical model | Output resolution | Strength | Fails when |
|------|--------------|------------------|-------------------|-------------------|----------|------------|
| PureCLIP (Krakau 2017) | iCLIP/eCLIP/PAR-CLIP | Truncation + enrichment + CL motif | Non-homogeneous HMM | Single-nucleotide | The most comprehensive HMM model | F1 ~0.2 on broad-binding RBPs (only ~4 sites per CLIP on benchmark) |
| CTK CITS (Shah 2017) | iCLIP/eCLIP | Truncation only | Empirical FDR vs background | Single-nucleotide | Simple, well-validated for iCLIP | Less granular than PureCLIP; no HMM smoothing |
| CTK CIMS deletion (Shah 2017) | HITS-CLIP | Single-base deletions | Empirical FDR | Single-nucleotide | Standard for HITS-CLIP | Requires deletion-tolerant aligner (BWA -e) |
| CTK CIMS T->C (Shah 2017) | PAR-CLIP | T->C substitutions | Empirical FDR | Single-nucleotide | Alternative to PARalyzer | Less popular than PARalyzer |
| PARalyzer (Corcoran 2011) | PAR-CLIP | T->C transitions, kernel density | Kernel density estimation | Cluster (10-50 nt) | Field standard for PAR-CLIP | Cluster-level, not single-nt; needs careful parameter tuning |
| wavClusteR (Comoglio 2015) | PAR-CLIP | T->C with wavelet smoothing | Wavelet transform + clustering | Cluster | Robust to sequencing depth variability | Less single-nt; legacy R package |
| pyCRAC / kPLogo | CRAC (yeast) | Read truncation + deletion | Empirical | Single-nucleotide | Original CLIP analytics tool | Yeast-focused; perl/python legacy |
| Piranha bins | Any | Coverage in bins | ZTNB | Bin-level (50-200 nt) | Not crosslink-specific | Bin width too coarse for single-nt |
| HOMER tag2pos | Any | 5' end position | None | Read-end positions only | Quick truncation site dump | No statistical filtering |
| omniCLIP | Any | Coverage + variant pattern | Dirichlet-multinomial HMM | Region | Not focal | Too broad for single-nt |
| iCount (paths) | iCLIP | Crosslink site cluster | Empirical | Cluster | Used in nf-core/clipseq | Less single-nt than PureCLIP |

Methodology evolves; PureCLIP 2.0 and CTK 1.2 have minor flag changes; PARalyzer parameters need careful tuning per-RBP. Cross-validate single-nt sites with at least two tools when high-confidence reporting is required.

## Critical Choice: Truncation vs Mutation vs HMM

Three orthogonal approaches:

**Truncation-based (CITS, PureCLIP)** -- Find positions enriched in cDNA 5' ends. The RT enzyme stops at the adduct; the 5' end of the read maps to CL site - 1. Used for iCLIP/eCLIP. Pro: high yield (~80% of reads truncate). Con: U bias of crosslinking inflates U positions.

**Mutation-based (CIMS for deletions/substitutions; PARalyzer for T->C)** -- Find positions with crosslink-induced mutations. RT reads through the adduct with high mutation rate at the CL position. Used for HITS-CLIP (deletions, 3-7%) and PAR-CLIP (T->C, 20-50%). Pro: less U bias; PAR-CLIP signal is restricted to T positions. Con: lower yield (only mutated reads count).

**HMM-based (PureCLIP, omniCLIP)** -- Joint model of enrichment + CL signature + sequence context. State the most-likely posterior probability of CL state per nucleotide. Used across CLIP variants. Pro: integrates multiple signals; explicit input normalization. Con: F1 ~0.2 on bulk RBPs (very focal); single-nt resolution at the cost of breadth.

| Goal | Method | Tool |
|------|--------|------|
| Single-nt CL map for motif registration (mCross) | HMM | PureCLIP |
| Truncation-based ENCODE-compatible iCLIP | Truncation | CTK CITS |
| PAR-CLIP T->C signal | Mutation | PARalyzer (cluster) or CTK CIMS sub T->C (single-nt) |
| HITS-CLIP deletion signal | Mutation | CTK CIMS deletion |
| Cross-CLIP-variant comparable single-nt | HMM | PureCLIP (all variants) |
| Allele-specific CLIP (BEAPR) | Truncation + ASB | PureCLIP CL sites + WASP-filtered BAM |
| Variant effect (computational from CL) | HMM | PureCLIP + downstream DeepRiPe (see clip-deep-learning) |
| Yeast / CRAC | Truncation + deletion | pyCRAC |

## Per-Tool Failure Modes

### PureCLIP -- HMM convergence on sparse coverage

**Trigger:** RBP with low IP enrichment; SMInput depth uneven; PureCLIP run on full genome without -iv interval.

**Mechanism:** PureCLIP HMM convergence is sensitive to coverage depth. On regions with < 5 reads/100 bp the HMM defaults to background state; if half the genome is sparse, the state distributions become bimodal and the entire run takes > 24 h or runs out of memory.

**Symptom:** PureCLIP "Convergence not reached" warning; or output has 0 CL sites; or runtime > 24 h on standard server.

**Fix:** Restrict analysis to expressed transcripts: `-iv expressed_tx.bed`. Or filter the BAM to high-coverage regions first. Test on chr22 first to verify parameters before genome-wide run.

### PureCLIP -- Misses broad binding zones

**Trigger:** RBP with broad binding (PUM2 3' UTRs, SR proteins exonic enhancers); using PureCLIP exclusively.

**Mechanism:** PureCLIP HMM emits the CL state at very high stringency. Skipper 2023 benchmark reports F1 ~0.2 on bulk RBPs - only ~4 sites per CLIP on test data.

**Symptom:** PureCLIP CL site count is 100x lower than CLIPper / Skipper peak count.

**Fix:** Use PureCLIP for single-nt CL maps AND CLIPper/Skipper for broad-zone peak calls. They are complementary, not substitutes. PureCLIP `-or` regions output is broader than `-o` sites but still focal.

### CITS truncation -- Misapplied to PAR-CLIP

**Trigger:** CTK CITS run on PAR-CLIP data.

**Mechanism:** PAR-CLIP RT reads through the adduct (T->C mutation), not truncates. The 5' end of PAR-CLIP reads is at the fragment 5' end, not the CL site.

**Symptom:** CITS returns few sites; the sites that ARE returned do not align with the T->C mutation positions.

**Fix:** For PAR-CLIP use PARalyzer or CTK CIMS substitution mode (`-type sub -nuc t -mut c`). CITS is iCLIP/eCLIP only.

### CIMS deletion -- Misapplied to iCLIP

**Trigger:** CTK CIMS run with `-type del` on iCLIP data.

**Mechanism:** iCLIP RT truncates, doesn't delete. The deletion rate is < 1% in iCLIP (background level), so CIMS finds nothing significant.

**Symptom:** CIMS output BED has < 100 sites genome-wide.

**Fix:** Use CTK CITS for iCLIP (truncation-based) or PureCLIP. CIMS deletion is HITS-CLIP / older PAR-CLIP only.

### PARalyzer -- Parameter sensitivity

**Trigger:** PARalyzer default parameters on a new RBP; cluster count unexpected.

**Mechanism:** PARalyzer's kernel-density approach has 5+ parameters: minimum read depth, minimum mutation read count, mutation rate threshold, bandwidth, kernel type. Different combinations produce 5-100x different cluster counts.

**Symptom:** PARalyzer cluster count 5-100x off from literature reports for the same RBP.

**Fix:** Use the published PARalyzer parameter file for that RBP class (Hafner 2010 supplementary). Validate on a known-binding region with a reference dataset before publication.

### Aligner deletion-tolerance for HITS-CLIP

**Trigger:** HITS-CLIP aligned with STAR or bowtie2 in standard mode; CIMS finds few deletion sites.

**Mechanism:** Standard STAR/bowtie2 are aggressive about handling deletions; CIMS expects the deletions to be visible as CIGAR D operations in BAM. STAR/bowtie2 may soft-clip or fail to align deletion-bearing reads.

**Symptom:** Few CIGAR D operations in BAM (`samtools view dedup.bam | awk '$6 ~ /D/' | wc -l`); CIMS deletion site count low.

**Fix:** Re-align HITS-CLIP with BWA-aln (Yeo lab convention) or STAR `--scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1` to permit deletions. Or use Novoalign which is more deletion-tolerant.

### PAR-CLIP T->C mismatch ceiling

**Trigger:** PAR-CLIP aligned with STAR `--outFilterMismatchNoverReadLmax 0.04`; T->C reads discarded.

**Mechanism:** PAR-CLIP per-T conversion rate 20-50% means a 30 nt read with 8 Ts may have 4 T->C events = 13% mismatch rate, exceeding the 4% ceiling.

**Symptom:** 40-70% read loss at alignment for PAR-CLIP; downstream T->C site detection finds nothing.

**Fix:** Raise to 0.07 for PAR-CLIP only (see clip-seq/clip-alignment).

### Strand mis-assignment for truncation site

**Trigger:** Single-end CLIP; using read 5' end as CL site without strand consideration.

**Mechanism:** For plus-strand RNA, the cDNA 5' end is downstream of the CL position. For minus-strand RNA, it is upstream. Without strand-aware coordinate conversion, half of CL sites are mis-positioned by ~1 nt.

**Symptom:** mCross / motif analysis shows the motif positioned 1 nt off from expected.

**Fix:** Verify the BED output handles strand correctly. CTK CITS, PureCLIP, and most modern tools do this automatically; custom scripts must add `if strand == '-': cl_position = read_end_5p + 1`.

## Decision Tree by Use Case

| Scenario | Tool | Parameters |
|----------|------|-----------|
| Single-nt CL map for mCross motif registration | PureCLIP | `-i dedup.bam -ibam sminput.bam -dm 8` |
| iCLIP/eCLIP CITS truncation sites | CTK CITS | `tag2cluster.pl -cs5 5 -m 1` |
| HITS-CLIP deletion CIMS | CTK CIMS | `CIMS.pl -big -c -p 0.01` |
| PAR-CLIP T->C clusters | PARalyzer | Per-RBP params from Hafner 2010 |
| PAR-CLIP T->C single-nt | CTK CIMS substitution | `CIMS.pl -type sub -nuc t -mut c -p 0.001` |
| Yeast CRAC | pyCRAC | `pyCRAC.py` with HTP-tagged protein |
| Variant effect at CL sites | PureCLIP + DeepRiPe | See clip-seq/clip-deep-learning |
| Allele-specific binding | WASP-filtered BAM + PureCLIP | See clip-seq/clip-alignment for WASP |
| Cross-CLIP comparison | PureCLIP (consistent across variants) | Same parameters |
| m6A miCLIP2 sites | miCLIP2-specific (m6Aboost) | See clip-seq/m6a-clip |
| STAMP edit sites | Not crosslink-based | See clip-seq/stamp-antibody-free |

## Reconciliation: When Detection Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| PureCLIP << CTK CITS site count | PureCLIP HMM very focal; CITS empirical broader | Both correct; report which threshold |
| PARalyzer clusters << CIMS T->C sites for PAR-CLIP | PARalyzer is cluster-level; CIMS is single-nt | Aggregate CIMS to clusters for comparison |
| PureCLIP CL sites at unexpected positions | Strand mis-assignment OR aligner soft-clip | Check BED strand column; verify `--alignEndsType EndToEnd` |
| HITS-CLIP CIMS empty | Aligner not deletion-tolerant | Re-align with BWA-aln or STAR with adjusted scoring |
| PAR-CLIP CIMS << expected | Reads lost at mismatch ceiling | Raise STAR `--outFilterMismatchNoverReadLmax` to 0.07 |
| eCLIP CITS positions shifted by 1 nt | Strand handling off-by-one | Verify R2 5' end position handling |
| Multiple CITS sites within 5 nt | RT stops near each other; same CL event | Cluster within 10 nt window post-detection |
| Motif registers 5 nt off in mCross | CL positions wrong by RT-stop offset | Confirm tool reports "CL - 1" position consistently |

**Operational rule:** For motif registration and ASB, run PureCLIP with SMInput. For ENCODE-comparable truncation-based output, also run CTK CITS. Cross-validate by checking that mCross motif PWM is the same when fed PureCLIP sites vs CTK CITS sites. If they diverge by > 2 nt in motif position, suspect strand or off-by-one issue.

## Workflow: PureCLIP -> mCross -> ASB

```bash
# Step 1: PureCLIP single-nt sites
pureclip \
    -i sample.dedup.bam -bai sample.dedup.bam.bai \
    -g genome.fa \
    -ibam sminput.dedup.bam -ibai sminput.dedup.bam.bai \
    -o sample.crosslinks.bed \
    -or sample.regions.bed \
    -nt 8 -dm 8 -iv expressed_tx.bed

# Step 2: mCross motif registration
mCross -i sample.crosslinks.bed -g genome.fa -k 7 -n 5 -o mcross_out

# Step 3: Intersect with heterozygous SNPs for allele-specific binding
bedtools intersect -wa -wb -s -a sample.crosslinks.bed -b het_snps.vcf > cl_at_hets.bed

# Step 4: Test allele bias with BEAPR or ASPRIN
# (see allele-specific CLIP literature; not in this skill)
```

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| PureCLIP "Convergence not reached" | HMM didn't converge on sparse coverage | Restrict to expressed transcripts with `-iv expressed.bed` |
| PureCLIP > 24h runtime | Full genome on large library | Test on chr22 first; add `-iv` filter |
| CTK CIMS empty for HITS-CLIP | Aligner not deletion-tolerant | Re-align with BWA-aln |
| CITS truncation sites not single-nt | CIGAR S operations from soft-clip | Verify `--alignEndsType EndToEnd` upstream |
| PAR-CLIP CIMS empty | Reads lost at mismatch ceiling | Raise STAR mismatch tolerance |
| Off-by-one motif position | Strand handling differs | Verify tool documents whether output is CL or CL-1 |
| Single-nt motif lost in HOMER on PureCLIP sites | Site BED too narrow; need flanking | Extend `bedtools slop -b 15` before motif analysis |
| mCross "no crosslinks found" | Passed peak BED instead of CL BED | Use single-nt site output |
| PureCLIP output looks like uniform distribution | R2 5' end trimmed in preprocessing | Re-preprocess; -g and --trim_front2 banned for CLIP |
| Strand-specific motif inverted | BED column 6 wrong | Verify strand encoding; PureCLIP preserves correctly |

## References

- Konig J et al 2010 Nat Struct Mol Biol 17:909 (iCLIP truncation principle)
- Hafner M et al 2010 Cell 141:129 (PAR-CLIP T->C signature)
- Granneman S et al 2009 PNAS 106:9613 (CRAC, single-nt CL detection)
- Corcoran DL et al 2011 Genome Biol 12:R79 (PARalyzer)
- Comoglio F et al 2015 BMC Bioinformatics 16:32 (wavClusteR)
- Shah A et al 2017 Bioinformatics 33:566 (CTK / CIMS / CITS)
- Krakau S et al 2017 Genome Biol 18:240 (PureCLIP HMM)
- Sugimoto Y et al 2017 Nature 519:491 (iCLIP U-crosslink bias)
- Feng H et al 2019 Mol Cell 74:1189 (mCross requires CL sites)
- Wu B et al 2018 Nat Commun 9:5117 (BEAPR allele-specific binding)
- Van Nostrand EL et al 2020 Nature 583:711 (ENCODE eCLIP single-nt analyses)

## Related Skills

- clip-seq/clip-preprocessing - 5' base preservation critical for truncation detection
- clip-seq/clip-alignment - End-to-end alignment for crosslink-preserving BAM
- clip-seq/clip-peak-calling - Peak vs CL site is a complementary distinction
- clip-seq/clip-motif-analysis - mCross consumes CL sites
- clip-seq/clip-deep-learning - RBPNet trained on CL count distributions
- clip-seq/m6a-clip - miCLIP2 has its own CL detection
- clip-seq/stamp-antibody-free - STAMP uses C->U editing, not crosslink
- alignment-files/sam-bam-basics - CIGAR string semantics for D operations
