# Reference: pandas 2.2+, numpy 1.26+ | Verify API if version differs
# Randomize samples onto 8-position EPIC BeadChips so technical batch is orthogonal to
# phenotype BEFORE the bench. This is the single uncorrectable EWAS decision (the no-rescue
# theorem): a case-chip-A/control-chip-B layout cannot be fixed by any downstream correction.

import numpy as np
import pandas as pd

positions_per_chip = 8   # EPIC/EPICv2 BeadChips hold 8 samples (R01C01..R08C01); 450K holds 12 - set per array


def assign_balanced_layout(sheet, group_col='pheno', seed=0):
    'Interleave samples across chips so each chip holds a balanced mix of phenotype groups.'
    rng = np.random.default_rng(seed)
    sheet = sheet.copy()
    order = []
    for grp, idx in sheet.groupby(group_col).groups.items():
        shuffled = rng.permutation(list(idx))
        order.append(pd.Series(shuffled, name=grp))
    interleaved = [i for chunk in zip_longest(*order) for i in chunk if i is not None]
    slots = np.arange(len(interleaved))
    sheet.loc[interleaved, 'chip'] = slots // positions_per_chip + 1
    sheet.loc[interleaved, 'position'] = slots % positions_per_chip + 1
    return sheet


def zip_longest(*seqs):
    n = max(len(s) for s in seqs)
    for i in range(n):
        yield tuple(s.iloc[i] if i < len(s) else None for s in seqs)


if __name__ == '__main__':
    rng = np.random.default_rng(1)
    n = 96
    sheet = pd.DataFrame({'sample': [f's{i}' for i in range(n)],
                          'pheno': rng.choice(['case', 'control'], n),
                          'sex': rng.choice(['F', 'M'], n)})
    laid_out = assign_balanced_layout(sheet, group_col='pheno', seed=42).astype({'chip': int, 'position': int})

    balance = pd.crosstab(laid_out['chip'], laid_out['pheno'])
    single_group_chips = (balance == 0).any(axis=1).sum()
    print(balance)
    print(f'single-group chips: {single_group_chips} (want 0 so batch is orthogonal to phenotype)')
