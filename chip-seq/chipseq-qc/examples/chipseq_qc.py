import pysam
import pybedtools
import argparse

def calculate_frip(bam_file, peak_file):
    '''Calculate Fraction of Reads in Peaks'''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    total_reads = bam.count(read_callback=lambda r: not r.is_unmapped and not r.is_secondary)

    peaks = pybedtools.BedTool(peak_file)
    reads_in_peaks = 0
    for peak in peaks:
        reads_in_peaks += bam.count(peak.chrom, peak.start, peak.end)

    bam.close()
    return reads_in_peaks / total_reads if total_reads > 0 else 0

def calculate_nrf(bam_file):
    '''Calculate Non-Redundant Fraction'''
    bam = pysam.AlignmentFile(bam_file, 'rb')
    positions = set()
    total = 0
    for read in bam.fetch():
        if not read.is_unmapped and not read.is_secondary:
            total += 1
            pos = (read.reference_name, read.reference_start, read.is_reverse)
            positions.add(pos)
    bam.close()
    return len(positions) / total if total > 0 else 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ChIP-seq QC metrics')
    parser.add_argument('bam', help='Input BAM file')
    parser.add_argument('--peaks', help='Peak file for FRiP calculation')
    args = parser.parse_args()

    print(f'BAM: {args.bam}')
    if args.peaks:
        frip = calculate_frip(args.bam, args.peaks)
        print(f'FRiP: {frip:.4f}')

    nrf = calculate_nrf(args.bam)
    print(f'NRF: {nrf:.4f}')
