'''Map a two-chain interface and quantify it by buried surface area.'''
# Reference: biopython 1.83+ | Verify API if version differs
import sys
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

# Heavy-atom 4.5A is the midpoint of the 4-5A direct-vdW-contact regime and, by
# excluding hydrogens, is robust to whether H atoms were modeled. State the cutoff.
CONTACT_CUTOFF = 4.5
# 1.4A ~ radius of a water molecule; it MUST be identical across all three SASA
# terms below, otherwise BSA = SASA_A + SASA_B - SASA_complex is meaningless.
PROBE_RADIUS = 1.4

parser = PDBParser(QUIET=True)


def interface_residues(model, chain_a, chain_b, cutoff):
    atoms = [a for a in model.get_atoms() if a.element != 'H']
    ns = NeighborSearch(atoms)
    res_a, res_b = set(), set()
    for r1, r2 in ns.search_all(cutoff, level='R'):
        c1, c2 = r1.get_parent().id, r2.get_parent().id
        if c1 == chain_a and c2 == chain_b:
            res_a.add(r1); res_b.add(r2)
        elif c1 == chain_b and c2 == chain_a:
            res_b.add(r1); res_a.add(r2)
    return res_a, res_b


def chain_set_sasa(path, keep_chains):
    structure = parser.get_structure('s', path)
    model = structure[0]
    for chain in list(model):
        if chain.id not in keep_chains:
            model.detach_child(chain.id)
    ShrakeRupley(probe_radius=PROBE_RADIUS).compute(model, level='C')
    return sum(chain.sasa for chain in model)


def buried_surface_area(path, chain_a, chain_b):
    complex_sasa = chain_set_sasa(path, {chain_a, chain_b})
    a_sasa = chain_set_sasa(path, {chain_a})
    b_sasa = chain_set_sasa(path, {chain_b})
    return a_sasa + b_sasa - complex_sasa


def main(path, chain_a='A', chain_b='B'):
    model = parser.get_structure('complex', path)[0]
    res_a, res_b = interface_residues(model, chain_a, chain_b, CONTACT_CUTOFF)
    print(f'Interface at {CONTACT_CUTOFF}A heavy-atom contact:')
    print(f'  chain {chain_a}: {len(res_a)} residues')
    print(f'  chain {chain_b}: {len(res_b)} residues')

    bsa = buried_surface_area(path, chain_a, chain_b)
    print(f'Buried surface area: {bsa:.0f} A^2 total, {bsa / 2:.0f} A^2 per partner')
    # BSA magnitude is one probabilistic signal of a biological interface, not proof;
    # small biological and large crystal interfaces both exist (Levy 2010).
    verdict = 'suggestive of biological' if bsa / 2 > 800 else 'small - could be crystal packing'
    print(f'  (per-partner {bsa / 2:.0f} A^2: {verdict}; corroborate with PISA/conservation)')


if __name__ == '__main__':
    args = sys.argv[1:]
    if not args:
        print('usage: interface_bsa.py complex.pdb [chainA chainB]')
        sys.exit(1)
    main(*args)
