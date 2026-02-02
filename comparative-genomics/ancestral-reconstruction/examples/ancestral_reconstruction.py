'''Ancestral sequence reconstruction'''

import subprocess
import os


def write_asr_control(alignment, tree, outfile, output_dir='.'):
    '''Generate PAML control file for ancestral reconstruction

    RateAncestor = 1 enables reconstruction
    Output goes to RST file in working directory
    '''
    ctl = f'''      seqfile = {alignment}
     treefile = {tree}
      outfile = {outfile}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 2
        model = 2
    aaRatefile = lg.dat

        clock = 0
        Mgene = 0

    fix_alpha = 0
        alpha = 0.5
       Malpha = 0
        ncatG = 4

 RateAncestor = 1
    cleandata = 0
'''
    ctl_path = os.path.join(output_dir, 'asr.ctl')
    with open(ctl_path, 'w') as f:
        f.write(ctl)

    return ctl_path


def parse_rst_ancestors(rst_file):
    '''Parse ancestral sequences from PAML RST file

    Returns dict mapping node names to sequences
    Node numbering starts after extant sequences
    '''
    ancestors = {}
    current_node = None
    sequence_lines = []

    with open(rst_file) as f:
        in_sequence = False
        for line in f:
            if 'node #' in line.lower():
                if current_node and sequence_lines:
                    ancestors[current_node] = ''.join(sequence_lines).replace(' ', '')
                # Extract node number
                parts = line.split('#')
                if len(parts) > 1:
                    node_num = parts[1].split()[0]
                    current_node = f'Node{node_num}'
                    sequence_lines = []
                    in_sequence = True
            elif in_sequence and line.strip() and not line.startswith(' '):
                # Check if this is a sequence line
                if all(c in 'ACDEFGHIKLMNPQRSTVWY-' for c in line.strip().upper()):
                    sequence_lines.append(line.strip())
                elif sequence_lines:
                    in_sequence = False

    if current_node and sequence_lines:
        ancestors[current_node] = ''.join(sequence_lines).replace(' ', '')

    return ancestors


def extract_site_probabilities(rst_file):
    '''Extract posterior probabilities for ancestral states

    Confidence levels:
    - P > 0.95: High confidence (reliable for resurrection)
    - P 0.80-0.95: Moderate confidence
    - P < 0.80: Low confidence (consider alternatives)
    '''
    site_probs = []
    in_prob_section = False

    with open(rst_file) as f:
        for line in f:
            if 'Prob distribution' in line or 'posterior' in line.lower():
                in_prob_section = True
                continue

            if in_prob_section:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        site = int(parts[0])
                        state = parts[1]
                        prob = float(parts[2])
                        site_probs.append({
                            'site': site,
                            'state': state,
                            'probability': prob
                        })
                    except (ValueError, IndexError):
                        if site_probs:
                            in_prob_section = False

    return site_probs


def assess_reconstruction_quality(site_probs):
    '''Assess overall quality of ancestral reconstruction

    Returns quality metrics for experimental planning
    '''
    if not site_probs:
        return {'quality': 'unknown', 'message': 'No probability data found'}

    probs = [s['probability'] for s in site_probs]
    n_sites = len(probs)

    # Count sites by confidence
    high_conf = sum(1 for p in probs if p > 0.95)
    mod_conf = sum(1 for p in probs if 0.8 <= p <= 0.95)
    low_conf = sum(1 for p in probs if p < 0.8)

    mean_prob = sum(probs) / n_sites
    high_conf_frac = high_conf / n_sites

    # Overall quality assessment
    # For protein resurrection:
    # - >90% high confidence sites: Excellent, proceed with ML sequence
    # - 70-90% high confidence: Good, but test alternatives at uncertain sites
    # - <70% high confidence: Poor, may need better alignment/tree
    if high_conf_frac > 0.9:
        quality = 'excellent'
    elif high_conf_frac > 0.7:
        quality = 'good'
    elif high_conf_frac > 0.5:
        quality = 'moderate'
    else:
        quality = 'poor'

    return {
        'total_sites': n_sites,
        'high_confidence': high_conf,
        'moderate_confidence': mod_conf,
        'low_confidence': low_conf,
        'mean_probability': mean_prob,
        'high_conf_fraction': high_conf_frac,
        'quality': quality
    }


def identify_ambiguous_sites(site_probs, threshold=0.8):
    '''Find sites with ambiguous ancestral states

    These sites should be tested with alternative states
    in resurrection experiments
    '''
    return [s for s in site_probs if s['probability'] < threshold]


def summarize_asr_results(ancestors, quality_metrics, ambiguous_sites):
    '''Print summary of ancestral reconstruction'''
    print('Ancestral Sequence Reconstruction Results')
    print('=' * 50)

    print(f"\nReconstructed {len(ancestors)} ancestral nodes")
    for node, seq in list(ancestors.items())[:3]:
        print(f"  {node}: {seq[:50]}..." if len(seq) > 50 else f"  {node}: {seq}")

    print(f"\nReconstruction Quality: {quality_metrics['quality'].upper()}")
    print(f"  Total sites: {quality_metrics['total_sites']}")
    print(f"  High confidence (P>0.95): {quality_metrics['high_confidence']} "
          f"({quality_metrics['high_conf_fraction']*100:.1f}%)")
    print(f"  Low confidence (P<0.80): {quality_metrics['low_confidence']}")
    print(f"  Mean probability: {quality_metrics['mean_probability']:.3f}")

    if ambiguous_sites:
        print(f"\nAmbiguous sites requiring attention: {len(ambiguous_sites)}")
        for site in ambiguous_sites[:5]:
            print(f"  Site {site['site']}: {site['state']} (P={site['probability']:.3f})")
        if len(ambiguous_sites) > 5:
            print(f"  ... and {len(ambiguous_sites) - 5} more")


if __name__ == '__main__':
    print('Ancestral Sequence Reconstruction')
    print('=' * 50)

    # Example results (simulated)
    example_ancestors = {
        'Node6': 'MKFLILLFNILCLFPVLAADYKDDDDKGENLYFQG',
        'Node7': 'MKFLILLFNILCLFPVLAADYKDDDDKGENLYFQG',
        'Node8': 'MKFLVLLFNILCLFPVLAADYKDDDDKGDNLYFQG',
    }

    example_probs = [
        {'site': 1, 'state': 'M', 'probability': 0.999},
        {'site': 2, 'state': 'K', 'probability': 0.987},
        {'site': 3, 'state': 'F', 'probability': 0.923},
        {'site': 4, 'state': 'L', 'probability': 0.756},  # Ambiguous
        {'site': 5, 'state': 'I', 'probability': 0.812},
        {'site': 6, 'state': 'L', 'probability': 0.634},  # Ambiguous
    ]

    quality = assess_reconstruction_quality(example_probs)
    ambiguous = identify_ambiguous_sites(example_probs, threshold=0.8)

    summarize_asr_results(example_ancestors, quality, ambiguous)

    print('\n\nTo run on real data:')
    print('1. Prepare protein MSA (PHYLIP format)')
    print('2. Generate phylogenetic tree')
    print('3. ctl = write_asr_control(aln, tree, "asr.mlc")')
    print('4. subprocess.run(["codeml", ctl])')
    print('5. ancestors = parse_rst_ancestors("rst")')
    print('6. probs = extract_site_probabilities("rst")')
