'''Design HDR donor templates for CRISPR knock-ins'''

from Bio.Seq import Seq

# Common protein tags for knock-ins
TAGS = {
    'FLAG': 'GATTACAAGGATGACGATGACAAG',           # 24bp, 8aa
    '3xFLAG': 'GATTACAAGGATGACGATGACAAGGATTACAAGGATGACGATGACAAGGATTACAAGGATGACGATGACAAG',
    'HA': 'TACCCATACGATGTTCCAGATTACGCT',          # 27bp, 9aa
    'V5': 'GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG',  # 42bp, 14aa
    'MYC': 'GAACAAAAACTCATCTCAGAAGAGGATCTG',      # 30bp, 10aa
    '6xHIS': 'CATCACCATCACCATCAC',                # 18bp, 6aa
}

# Flexible linkers
LINKERS = {
    'GS': 'GGCAGC',                    # Short GS linker
    'GSGS': 'GGCAGCGGCAGC',            # Medium linker
    'GSGSG': 'GGCAGCGGCAGCGGC',        # Longer linker
}

# Optimal homology arm lengths
# ssODN: 30-60 nt each side (total ~120-200 nt)
# dsDNA: 200-800 bp each side for larger insertions
# Asymmetric: PAM-distal arm can be longer (improves HDR)
SSODN_ARM_LENGTH = (30, 60)
DSDNA_ARM_LENGTH = (200, 800)


def design_ssodn(target_seq, cut_site, insert_seq='', left_arm=50, right_arm=50):
    '''Design single-stranded oligodeoxynucleotide donor

    Args:
        target_seq: Genomic sequence around cut site
        cut_site: Position where Cas9 cuts (0-indexed)
        insert_seq: Sequence to insert
        left_arm: Left homology arm length (30-60nt optimal)
        right_arm: Right homology arm length (30-60nt optimal)

    Total length should be 100-200nt for efficient synthesis.
    Longer oligos have higher error rates and cost.
    '''
    left = target_seq[cut_site - left_arm:cut_site]
    right = target_seq[cut_site:cut_site + right_arm]

    ssodn = left + insert_seq + right
    ssodn_rc = str(Seq(ssodn).reverse_complement())

    return {
        'sense': ssodn,
        'antisense': ssodn_rc,
        'total_length': len(ssodn),
        'left_arm_length': left_arm,
        'right_arm_length': right_arm,
        'insert_length': len(insert_seq),
        'note': 'Test both sense and antisense - efficiency varies by locus'
    }


def design_asymmetric_ssodn(target_seq, cut_site, insert_seq, pam_side='right'):
    '''Design ssODN with asymmetric homology arms

    Asymmetric arms can improve HDR efficiency by ~2-3 fold.
    The PAM-distal arm should be longer (resected first during repair).

    PAM-proximal: 30-40nt
    PAM-distal: 60-90nt
    '''
    if pam_side == 'right':
        left_arm, right_arm = 70, 35  # PAM-distal longer
    else:
        left_arm, right_arm = 35, 70

    return design_ssodn(target_seq, cut_site, insert_seq, left_arm, right_arm)


def design_tag_insertion(target_seq, cut_site, tag_name, add_linker=True):
    '''Design ssODN for protein tagging

    Common use: Add epitope tag for immunoprecipitation or imaging.
    Position cut site at the codon boundary where tag will be inserted.
    '''
    tag = TAGS.get(tag_name.upper(), tag_name)

    if add_linker:
        insert = LINKERS['GSGS'] + tag  # Linker before tag
    else:
        insert = tag

    return design_ssodn(target_seq, cut_site, insert)


def design_point_mutation(target_seq, mutation_pos, new_base, arm_length=50):
    '''Design ssODN for a point mutation

    Centers the mutation in the ssODN for optimal incorporation.
    Consider also mutating the PAM to prevent re-cutting.
    '''
    mutant = list(target_seq)
    mutant[mutation_pos] = new_base
    mutant_seq = ''.join(mutant)

    left_start = mutation_pos - arm_length
    right_end = mutation_pos + arm_length + 1

    ssodn = mutant_seq[max(0, left_start):right_end]

    return {
        'sequence': ssodn,
        'length': len(ssodn),
        'mutation_position': min(arm_length, mutation_pos),
        'change': f'{target_seq[mutation_pos]}>{new_base}'
    }


def add_pam_mutation(donor_seq, pam_position):
    '''Add silent mutation to disrupt PAM and prevent re-cutting

    After HDR, the PAM (NGG) should be mutated to prevent Cas9
    from cutting the corrected allele.

    Strategy: Change NGG -> NGA (most common silent option)
    '''
    donor = list(donor_seq)

    # Check we have NGG at expected position
    if pam_position + 2 < len(donor):
        if donor[pam_position + 1] == 'G' and donor[pam_position + 2] == 'G':
            donor[pam_position + 2] = 'A'  # GG -> GA

    return ''.join(donor)


def design_dsdna_donor(target_seq, cut_site, insert_seq, arm_length=500):
    '''Design double-stranded DNA donor for larger insertions

    For insertions >50bp, dsDNA donors are more efficient.
    Can be delivered as PCR product or plasmid.

    Arm length: 200-800bp recommended for high efficiency
    '''
    left_arm = target_seq[cut_site - arm_length:cut_site]
    right_arm = target_seq[cut_site:cut_site + arm_length]

    donor = left_arm + insert_seq + right_arm

    return {
        'full_sequence': donor,
        'total_length': len(donor),
        'left_arm': left_arm,
        'right_arm': right_arm,
        'left_arm_length': arm_length,
        'right_arm_length': arm_length,
        'insert': insert_seq
    }


if __name__ == '__main__':
    # Example: Design ssODN for FLAG tag insertion
    # Simulated target sequence around cut site
    target = 'A' * 100 + 'ATGATCGATCGATCGATCGATCGAGG' + 'T' * 100
    cut_site = 110  # Position where Cas9 will cut

    print('HDR Template Design Examples')
    print('=' * 50)

    # 1. FLAG tag insertion
    print('\n1. FLAG tag insertion:')
    flag_donor = design_tag_insertion(target, cut_site, 'FLAG')
    print(f"   Sense ssODN length: {flag_donor['total_length']}nt")
    print(f"   Arms: {flag_donor['left_arm_length']}nt | {flag_donor['insert_length']}nt insert | {flag_donor['right_arm_length']}nt")

    # 2. Point mutation
    print('\n2. Point mutation (A->G):')
    mutation_donor = design_point_mutation(target, 115, 'G')
    print(f"   ssODN length: {mutation_donor['length']}nt")
    print(f"   Change: {mutation_donor['change']}")

    # 3. Asymmetric arms
    print('\n3. Asymmetric arm design:')
    asym_donor = design_asymmetric_ssodn(target, cut_site, TAGS['HA'])
    print(f"   Total length: {asym_donor['total_length']}nt")
    print(f"   Left arm: {asym_donor['left_arm_length']}nt (PAM-distal)")
    print(f"   Right arm: {asym_donor['right_arm_length']}nt (PAM-proximal)")

    # 4. dsDNA donor
    print('\n4. dsDNA donor for GFP insertion:')
    gfp_seq = 'ATG' + 'N' * 714 + 'TAA'  # ~720bp GFP CDS
    dsdna = design_dsdna_donor(target, cut_site, gfp_seq, arm_length=500)
    print(f"   Total length: {dsdna['total_length']}bp")
    print(f"   Homology arms: {dsdna['left_arm_length']}bp each")
