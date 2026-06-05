#!/usr/bin/env python3
'''
cfDNA preprocessing with UMI/duplex consensus using fgbio.

Encodes the verified fgbio chain: extract UMIs -> align -> group (paired for duplex,
adjacency for simplex) -> call consensus permissively -> RE-align -> filter strictly.
Consensus reads are emitted unmapped by design, so the second alignment is mandatory.
'''
# Reference: bwa 0.7.17+, fgbio 2.1+, numpy 1.26+, pysam 0.22+, samtools 1.19+ | Verify API if version differs

import subprocess
import pysam
import numpy as np
from pathlib import Path


def _run(cmd):
    subprocess.run(cmd, check=True)


def preprocess_cfdna(input_bam, output_bam, reference, read_structure='6M11S+T',
                     duplex=False, threads=8):
    '''Run the cfDNA UMI/duplex consensus pipeline from an unmapped UMI-bearing BAM.

    read_structure tokens: M=UMI, S=skip/stem (omitting it bleeds UMI into template), T=template, +=all remaining.
    duplex=True uses paired grouping + CallDuplexConsensusReads + a true-duplex filter (2 1 1);
    duplex=False uses adjacency grouping + CallMolecularConsensusReads + a simplex filter (single value).
    '''
    work = Path(output_bam).parent
    stem = Path(output_bam).stem
    with_umis = work / f'{stem}_umis.bam'
    aligned = work / f'{stem}_aligned.bam'
    grouped = work / f'{stem}_grouped.bam'
    consensus_unmapped = work / f'{stem}_consensus_unmapped.bam'
    consensus_mapped = work / f'{stem}_consensus.bam'

    _run(['fgbio', 'ExtractUmisFromBam', '--input', str(input_bam), '--output', str(with_umis),
          '--read-structure', read_structure, read_structure, '--single-tag', 'RX'])

    # -Y soft-clips supplementaries so tag-bearing short-fragment sequence is not dropped.
    subprocess.run(f'bwa mem -t {threads} -Y {reference} {with_umis} | samtools sort -o {aligned} -',
                   shell=True, check=True)
    pysam.index(str(aligned))

    # paired grouping is mandatory for duplex; adjacency cannot reconstruct strand pairing.
    strategy = 'paired' if duplex else 'adjacency'
    _run(['fgbio', 'GroupReadsByUmi', '--input', str(aligned), '--output', str(grouped),
          '--strategy', strategy, '--edits', '1'])

    # Call permissively (--min-reads 1). For duplex this caller's --min-reads is a PRE-filter (fgbio #1009);
    # the real strictness is applied in FilterConsensusReads below.
    caller = 'CallDuplexConsensusReads' if duplex else 'CallMolecularConsensusReads'
    _run(['fgbio', caller, '--input', str(grouped), '--output', str(consensus_unmapped), '--min-reads', '1'])

    # RE-align: consensus reads are unmapped by design. ZipperBams transfers the per-base
    # consensus tags onto the new alignments; without it FilterConsensusReads has no tags to act on.
    realign = (f'samtools fastq {consensus_unmapped} | bwa mem -t {threads} -Y -p {reference} - '
               f'| fgbio ZipperBams --unmapped {consensus_unmapped} --ref {reference} '
               f'| samtools sort -o {consensus_mapped} -')
    subprocess.run(realign, shell=True, check=True)
    pysam.index(str(consensus_mapped))

    # min_reads "total strand1 strand2": "2 1 1" = true duplex (both strands seen). If values two and
    # three differ the more stringent must come first. Simplex uses a single value.
    min_reads = ['2', '1', '1'] if duplex else ['2']
    _run(['fgbio', 'FilterConsensusReads', '--input', str(consensus_mapped), '--output', str(output_bam),
          '--ref', str(reference), '--min-reads', *min_reads,
          '--max-read-error-rate', '0.025',   # fgbio default; per-read error ceiling
          '--max-base-error-rate', '0.1',     # fgbio default; per-consensus-base error ceiling
          '--min-base-quality', '40',         # mask consensus bases below Q40 to N
          '--reverse-per-base-tags'])         # per-base tags read in genomic orientation after alignment
    pysam.index(str(output_bam))
    return output_bam


def insert_size_qc(bam_path, max_size=600):
    '''Summarize cfDNA fragment lengths as a pre-analytical QC readout.

    Healthy: mode ~167 bp (mononucleosome) with a ~10.4 bp sawtooth below it.
    frac_over_250bp elevated signals leukocyte-lysis gDNA contamination (long fragments).
    A mode ~10 bp low is an adaptase (Swift/Accel-1S) chemistry signature, not contamination.
    '''
    bam = pysam.AlignmentFile(bam_path, 'rb')
    sizes = [r.template_length for r in bam.fetch()
             if r.is_proper_pair and not r.is_secondary and 0 < r.template_length <= max_size]
    bam.close()
    sizes = np.array(sizes)
    if len(sizes) == 0:
        return {'n': 0, 'mode_bp': 0, 'median_bp': 0.0, 'short_frac_90_150': 0.0, 'frac_over_250bp': 0.0}
    return {'n': len(sizes),
            'mode_bp': int(np.bincount(sizes).argmax()),
            'median_bp': float(np.median(sizes)),
            'short_frac_90_150': float(np.mean((sizes >= 90) & (sizes <= 150))),  # ctDNA-enriched window (Mouliere 2018)
            'frac_over_250bp': float(np.mean(sizes > 250))}


if __name__ == '__main__':
    print('cfDNA preprocessing pipeline')
    print('preprocess_cfdna(in_bam, out_bam, ref, duplex=False) - simplex UMI consensus')
    print('preprocess_cfdna(in_bam, out_bam, ref, duplex=True)  - duplex consensus (paired grouping, 2 1 1 filter)')
    print('insert_size_qc(bam) - fragment-length QC readout')
