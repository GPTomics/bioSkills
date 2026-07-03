'''Build a multi-species community with MICOM and predict community + per-member growth.

MICOM keeps each member compartmentalized and connects them through a shared extracellular medium,
avoiding the compartment-pooling artifact of a hand-merged single-pool "bag" model.
'''
# Reference: MICOM 0.33+, cobrapy 0.29+ | Verify API if version differs
# The __main__ guard is required: MICOM parallelizes and re-imports this module in workers.

from micom import Community
from micom.data import test_taxonomy


def main():
    # taxonomy: one row per taxon with id, file (per-taxon SBML), abundance. test_taxonomy() ships a
    # ready E. coli community; in practice abundances come from metagenomics (metagenomics/abundance-estimation).
    taxonomy = test_taxonomy().head(3)
    print(f'community members: {list(taxonomy["id"])}')

    community = Community(taxonomy, progress=False)   # compartmentalized multi-species model

    print('\n=== Cooperative tradeoff (community optimum, growth spread across members) ===')
    solution = community.cooperative_tradeoff(fraction=1.0)   # QP; needs HiGHS/CPLEX/Gurobi
    print(f'community growth rate: {solution.growth_rate:.4f} /h')
    members = solution.members.dropna(subset=['growth_rate'])
    for taxon, row in members.iterrows():
        print(f'  {taxon}: {row["growth_rate"]:.4f} /h')

    print('\nReminder: this community inherits every member model\'s errors; curate members first,')
    print('and confirm they share a namespace so metabolite exchange actually connects.')


if __name__ == '__main__':
    main()
