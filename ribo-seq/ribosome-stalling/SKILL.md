---
name: bio-ribo-seq-ribosome-stalling
description: Detect ribosome pausing and stalling at codon resolution from Ribo-seq, using local-relative occupancy metrics and A-site assignment. Use when studying elongation dynamics, codon dwell times, pause motifs, or ribosome collisions, and when judging whether a pause is real biology or a cycloheximide artifact.
tool_type: python
primary_tool: Plastid
---

## Version Compatibility

Reference examples tested with: plastid 0.6+, numpy 1.26+, scipy 1.12+, biopython 1.83+, twobitreader 3.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Ribosome Stalling Detection

**"Find ribosome pause sites in my data"** -> Detect codon positions where ribosomes dwell longer than the local average, attribute them to A-site decoding or nascent-chain effects, and judge whether the signal is real or a drug artifact.
- Python: `plastid` for A-site codon density, local-relative pause scoring, and motif context

## Read this first: cycloheximide destroys pause signal

A pause is only meaningful when footprint positions reflect in-vivo dwell times. Cycloheximide (CHX) pre-treatment of live cells violates this: arrest is not instantaneous, ribosomes run on after the drug, density redistributes downstream, codon-specific pausing is attenuated, and an artifactual start-codon peak appears. Hussmann 2015 showed CHX data report a WEAK NEGATIVE correlation between codon rate and tRNA abundance while flash-frozen data report a STRONG POSITIVE one -- the drug flips the conclusion. A pause analysis on CHX data largely measures the drug.

Decision rule before any dwell/pause analysis:
- Flash-frozen, no drug (or CHX only in lysis at high concentration) -> codon-resolution dwell-time analysis is valid.
- CHX pre-incubation of live cells -> restrict to qualitative/gene-level statements; do not report codon dwell times or tRNA correlations as biology.
- Harringtonine/lactimidomycin data are for initiation-site mapping, not elongation pausing (see initiation-site-mapping).

## A-site vs P-site: the offset choice changes the biology

The P-site offset (~12 nt from the 5' end for canonical 28-30 nt footprints) must be calibrated per read length, not hardcoded (see ribosome-periodicity). The relevant site depends on the mechanism: tRNA-availability/decoding pauses register at the A-SITE (A-site = P-site + 3), so codon-occupancy and tRNA work assign to the A-site. Nascent-chain effects (polyproline, charge) act at the P-site/exit tunnel and upstream. State which site is used; the peak position relative to A/P/E is itself diagnostic.

## Pause-metric selection

| Metric | Definition | Caveat |
|--------|-----------|--------|
| Per-transcript z-score | (density - gene mean)/gene SD | not the field standard; SD inflated by the peaks sought; arbitrary threshold |
| Pause score | local density / gene-mean density at that position | needs a per-gene coverage floor; the standard local-relative metric |
| Codon occupancy | mean over all instances of a codon of (position density / gene mean) | normalize each gene to its own mean FIRST, then pool; assign to A-site |
| RUST | binarize each position vs the gene mean, average the metafootprint | outlier-robust; resists a few high peaks dominating |
| Disome density | footprints from two stacked ribosomes (~58-62 nt) | the cleanest in-vivo strong-pause readout (Arpat 2020) |

The two normalization rules the naive z-score violates: never z-score across positions of differently-expressed genes (high-expression genes dominate) -- normalize each gene to its own mean first; and require a real per-gene coverage floor (a few hundred in-frame footprints), far above a `sum > 100` cutoff, or per-position metrics are noise.

## Calculate A-site codon density (plastid)

**Goal:** Get a per-codon occupancy vector for each CDS at the A-site.

**Approach:** Map footprints to the A-site offset, fetch the CDS count vector, and reduce each codon to its summed in-frame count.

```python
from plastid import BAMGenomeArray, GTF2_TranscriptAssembler, FivePrimeMapFactory
import numpy as np

def asite_codon_occupancy(bam_path, gtf_path, asite_offset=15):
    '''Per-codon A-site occupancy per CDS. A-site offset = P-site (~12) + 3.

    A single fixed offset is a simplification valid only when one read length
    dominates. For production, calibrate per length (ribosome-periodicity) and
    map with VariableFivePrimeMapFactory.from_file using A-site = P-site + 3.
    '''
    alignments = BAMGenomeArray(bam_path, mapping=FivePrimeMapFactory(offset=asite_offset))
    out = {}
    for tx in GTF2_TranscriptAssembler(gtf_path):
        if tx.cds_start is None:
            continue
        cds = tx.get_cds()
        counts = cds.get_counts(alignments)          # numpy vector over the CDS
        n_codons = len(counts) // 3
        # Sum the 3 positions of each codon into a SCALAR (one value per codon)
        per_codon = np.array([counts[i*3:i*3+3].sum() for i in range(n_codons)])
        out[tx.get_name()] = per_codon
    return out
```

`cds.get_counts(alignments)` is the count method on the SegmentChain; `BAMGenomeArray` has no `count_in_region`/`get_density`. Reducing each codon to a scalar (sum of its three positions) is essential -- storing the whole vector at each codon makes every downstream metric garbage.

## Score pauses with a local-relative metric

**Goal:** Flag codons where occupancy exceeds the gene's own average.

**Approach:** Divide each position by the gene mean (a pause score), require adequate coverage, and threshold.

```python
def pause_scores(per_codon_occupancy, min_total=500, score_threshold=5.0):
    '''Pause score = codon occupancy / gene-mean occupancy (local-relative).

    min_total: per-gene footprint floor; below this, scores are noise.
    score_threshold: fold-over-gene-mean to call a pause (tune per dataset).
    '''
    pauses = []
    for tx, occ in per_codon_occupancy.items():
        if occ.sum() < min_total:
            continue
        mean = occ.mean()
        if mean == 0:
            continue
        scores = occ / mean
        for pos in np.where(scores > score_threshold)[0]:
            pauses.append({'transcript': tx, 'codon': int(pos),
                           'pause_score': float(scores[pos])})
    return pauses
```

## Codon occupancy across genes

**Goal:** Estimate per-codon-type dwell, averaged across the transcriptome.

**Approach:** Normalize each gene to its own mean BEFORE pooling, then average per codon identity (the A-site codon).

Pool mean-of-ratios, not ratio-of-means: a raw average across genes is dominated by highly expressed genes. The per-codon occupancy is then a relative dwell estimate -- and only on no-drug data. The tRNA-availability correlation (codon occupancy vs tRNA adaptation index) is modest, sign- and protocol-dependent, and reflects charged-tRNA levels rather than gene copy number; report the effect size, not a presumed strong negative correlation.

## Known pause mechanisms

| Motif / feature | Mechanism |
|-----------------|-----------|
| Polyproline (PPP, PPG) | Rigid proline geometry stalls peptidyl transfer; rescued by eIF5A (eukaryotes) / EF-P (bacteria) |
| Poly-basic (Lys/Arg runs) | Basic nascent chain drags on the negatively-charged exit tunnel; poly-Lys also involves sliding on A-rich codons |
| Rare/low-tRNA codons | Slow A-site decoding; real but modest, and inflated in CHX data |
| Internal Shine-Dalgarno (bacteria) | Anti-SD base-pairing with 16S rRNA; real but contested (protocol-dependent) |

## Ribosome collisions and disome-seq (the modern readout)

When a ribosome stalls, the trailing ribosome collides into it, forming a disome whose ~58-62 nt footprint maps collision sites transcriptome-wide -- a cleaner in-vivo strong-pause readout than monosome relative density (Arpat 2020; ~10% of ribosomes can be in disomes). The collided-disome interface is the trigger for ribosome quality control: ZNF598 (mammals) / Hel2 (yeast) ubiquitinate small-subunit proteins, recruiting the splitting machinery and no-go decay. A monosome pause that coincides with a disome peak, replicates, and survives in no-drug data is strong evidence of a real, acted-upon stall.

## Extract pause-site sequence context

**Goal:** Find amino-acid motifs enriched at pause sites.

**Approach:** Build the per-transcript CDS sequences from a genome (plastid's `get_sequence` needs a genome, not a SegmentChain), then translate a window centered on the A-site codon of each pause.

```python
from Bio.Seq import Seq
import twobitreader

def cds_sequences_from_genome(gtf_path, twobit_path):
    '''Map transcript name -> spliced CDS nucleotide sequence.'''
    from plastid import GTF2_TranscriptAssembler
    genome = twobitreader.TwoBitFile(twobit_path)   # dict-like {chrom: seq}
    seqs = {}
    for tx in GTF2_TranscriptAssembler(gtf_path):
        if tx.cds_start is None:
            continue
        seqs[tx.get_name()] = tx.get_cds().get_sequence(genome)
    return seqs

def pause_motifs(pauses, cds_sequences, window_codons=5):
    '''Amino-acid context around each pause (centered on the A-site codon).'''
    motifs = []
    for p in pauses:
        seq = cds_sequences.get(p['transcript'])
        if not seq:
            continue
        c = p['codon']
        s, e = max(0, (c - window_codons) * 3), min(len(seq), (c + window_codons + 1) * 3)
        if (e - s) % 3 == 0:
            motifs.append(str(Seq(seq[s:e]).translate()))
    return motifs
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Strong start-codon "pause", odd tRNA correlation | CHX pre-treatment artifacts | Use flash-frozen no-drug data; restrict CHX data to gene-level claims |
| `AttributeError` on `count_in_region`/`get_density` | Not BAMGenomeArray methods | Use `cds.get_counts(alignments)` |
| Every codon occupancy identical | Whole count vector stored per codon | Store a scalar: sum the 3 positions of each codon |
| `TypeError` from `get_sequence` | Passed a SegmentChain, not a genome | Load a genome FASTA/2bit; call `cds.get_sequence(genome)` |
| Pauses dominated by one highly expressed gene | Global z-score / ratio-of-means | Normalize each gene to its own mean first; mean-of-ratios |
| Noisy, irreproducible pauses | Coverage floor too low (sum > 100) | Require a few hundred in-frame footprints per gene |
| tRNA correlation overstated | Assumed strong negative on CHX data | Report effect size; depends on charging and harvest |

## Related Skills

- ribosome-periodicity - Calibrate the A-site offset before scoring occupancy
- orf-detection - Locate the ORFs that pause sites fall within
- initiation-site-mapping - Distinguish initiation drugs from elongation pausing
- translation-efficiency - Gene-level translation context

## References

- Hussmann JA, Patchett S, Johnson A, Sawyer S, Press WH. 2015. Understanding biases in ribosome profiling experiments reveals signatures of translation dynamics in yeast. PLoS Genet 11(12):e1005732. doi:10.1371/journal.pgen.1005732
- Gerashchenko MV, Gladyshev VN. 2014. Translation inhibitors cause abnormalities in ribosome profiling experiments. Nucleic Acids Res 42(17):e134. doi:10.1093/nar/gku671
- O'Connor PBF, Andreev DE, Baranov PV. 2016. Comparative survey of the relative impact of mRNA features on local ribosome profiling read density. Nat Commun 7:12915. doi:10.1038/ncomms12915
- Arpat AB, Liechti A, De Matos M, Dreos R, Janich P, Gatfield D. 2020. Transcriptome-wide sites of collided ribosomes reveal principles of translational pausing. Genome Res 30(7):985-999. doi:10.1101/gr.257741.119
- Charneski CA, Hurst LD. 2013. Positively charged residues are the major determinants of ribosomal velocity. PLoS Biol 11(3):e1001508. doi:10.1371/journal.pbio.1001508
- Schuller AP, Wu CC, Dever TE, Buskirk AR, Green R. 2017. eIF5A functions globally in translation elongation and termination. Mol Cell 66(2):194-205. doi:10.1016/j.molcel.2017.03.003
- Li GW, Oh E, Weissman JS. 2012. The anti-Shine-Dalgarno sequence drives translational pausing and codon choice in bacteria. Nature 484(7395):538-541. doi:10.1038/nature10965
