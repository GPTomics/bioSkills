'''Prepare a structure for docking/MD: fill missing atoms, drop heterogens, add hydrogens at a stated pH'''
# Reference: pdbfixer 1.9+, openmm 8.1+ | Verify API if version differs

import argparse
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def prepare(source, output, ph, keep_water):
    # PDBFixer(pdbid=...) fetches from RCSB; filename=... reads a local PDB/mmCIF.
    fixer = PDBFixer(pdbid=source) if len(source) == 4 and '.' not in source else PDBFixer(filename=source)

    # Finder order matters: residues and nonstandard first so addMissingAtoms sees both sets.
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=keep_water)  # drops buffer/cryoprotectant; keep catalytic groups deliberately
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()  # built residues/atoms are HYPOTHESES - reported below as provenance

    # pH is a CHOICE, not a safe default: it sets the protonation of titratable residues.
    # Standard states are only defensible for freely solvated residues; buried/pocket/metal-adjacent
    # residues need a pKa predictor (PROPKA/H++) and an H-bond optimizer (reduce), not this shortcut.
    fixer.addMissingHydrogens(pH=ph)

    with open(output, 'w') as out:
        PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)

    # missingResidues is a dict {(chain_index, residue_index): [resname, ...]}; flatten it for provenance
    built = [(ci, pos, name) for (ci, pos), names in fixer.missingResidues.items() for name in names]
    print(f'prepared {source} -> {output} at pH={ph}, keep_water={keep_water}')
    print(f'built {len(fixer.missingResidues)} segment(s), {len(built)} residue(s): {built if built else "none"}')
    print('provenance: record pH, protonation model, and every built atom/loop alongside this file')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='PDBFixer structure preparation')
    p.add_argument('source', help='PDB id (e.g. 1VII) or path to a local PDB/mmCIF')
    p.add_argument('-o', '--output', default='prepared.pdb')
    p.add_argument('--ph', type=float, default=7.0)  # 7.0 = a reasonable but explicit physiological choice
    p.add_argument('--keep-water', action='store_true')
    args = p.parse_args()
    prepare(args.source, args.output, args.ph, args.keep_water)
