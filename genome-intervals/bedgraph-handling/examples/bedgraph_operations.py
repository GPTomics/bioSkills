import pyBigWig
import argparse
import sys

def bigwig_to_bedgraph(input_bw, output_bg, chrom=None, start=None, end=None):
    '''Convert bigWig to bedGraph format'''
    bw = pyBigWig.open(input_bw)

    with open(output_bg, 'w') as f:
        if chrom:
            intervals = bw.intervals(chrom, start, end)
            if intervals:
                for s, e, v in intervals:
                    f.write(f'{chrom}\t{s}\t{e}\t{v}\n')
        else:
            for chrom, size in bw.chroms().items():
                intervals = bw.intervals(chrom)
                if intervals:
                    for s, e, v in intervals:
                        f.write(f'{chrom}\t{s}\t{e}\t{v}\n')

    bw.close()

def bedgraph_stats(input_bg):
    '''Calculate statistics from bedGraph'''
    values = []
    total_bases = 0

    with open(input_bg) as f:
        for line in f:
            if line.startswith('track') or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, start, end, value = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            span = end - start
            values.append((value, span))
            total_bases += span

    weighted_sum = sum(v * s for v, s in values)
    mean_value = weighted_sum / total_bases if total_bases > 0 else 0
    max_value = max(v for v, s in values) if values else 0
    min_value = min(v for v, s in values) if values else 0

    print(f'Total bases: {total_bases}')
    print(f'Mean value: {mean_value:.4f}')
    print(f'Min value: {min_value:.4f}')
    print(f'Max value: {max_value:.4f}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='bedGraph operations')
    parser.add_argument('command', choices=['convert', 'stats'])
    parser.add_argument('input', help='Input file')
    parser.add_argument('--output', '-o', help='Output file')
    parser.add_argument('--chrom', help='Chromosome for region extraction')
    parser.add_argument('--start', type=int, help='Start position')
    parser.add_argument('--end', type=int, help='End position')

    args = parser.parse_args()

    if args.command == 'convert':
        if not args.output:
            print('Error: --output required for convert')
            sys.exit(1)
        bigwig_to_bedgraph(args.input, args.output, args.chrom, args.start, args.end)
    elif args.command == 'stats':
        bedgraph_stats(args.input)
