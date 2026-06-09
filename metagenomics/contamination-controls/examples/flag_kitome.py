"""Prevalence-style contaminant flagging from extraction blanks (decontam's prevalence idea).

A contaminant is enriched in negative controls relative to true samples. This mirrors decontam's
prevalence signal without R: a taxon present in a higher fraction of blanks than of real samples, and
matching a canonical kitome genus, is a contaminant suspect. Use the real decontam package for the
statistical test; this is a transparent screen and a teaching example.
"""
# Reference: pandas 2.2+, numpy 1.26+ | Verify API if version differs
import sys
import numpy as np
import pandas as pd

KITOME_GENERA = ['Bradyrhizobium', 'Ralstonia', 'Burkholderia', 'Pseudomonas',
                 'Acinetobacter', 'Sphingomonas', 'Methylobacterium', 'Stenotrophomonas']


def prevalence(table, columns):
    """Fraction of the given samples in which each taxon is present (> 0)."""
    return (table[columns] > 0).mean(axis=1)


def flag_contaminants(table, blanks, samples, kitome=KITOME_GENERA):
    """table: taxa x samples abundances. Flag taxa more prevalent in blanks than in real samples."""
    prev_blank = prevalence(table, blanks)
    prev_sample = prevalence(table, samples)
    enriched_in_blanks = prev_blank > prev_sample
    is_kitome = table.index.to_series().apply(lambda t: any(g in t for g in kitome))
    flagged = table.index[enriched_in_blanks]
    return pd.DataFrame({
        'prev_blank': prev_blank, 'prev_sample': prev_sample,
        'enriched_in_blanks': enriched_in_blanks, 'kitome_genus': is_kitome,
    }).loc[flagged].sort_values('prev_blank', ascending=False)


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if path:
        table = pd.read_csv(path, sep='\t', index_col=0)
        blanks = [c for c in table.columns if 'blank' in c.lower()]
        samples = [c for c in table.columns if c not in blanks]
    else:
        rng = np.random.default_rng(0)
        taxa = ['Faecalibacterium_prausnitzii', 'Bacteroides_uniformis', 'Ralstonia_pickettii',
                'Bradyrhizobium_sp', 'Escherichia_coli']
        samples = [f'sample_{i}' for i in range(6)]
        blanks = ['blank_1', 'blank_2']
        data = pd.DataFrame(0.0, index=taxa, columns=samples + blanks)
        data.loc[['Faecalibacterium_prausnitzii', 'Bacteroides_uniformis', 'Escherichia_coli'], samples] = \
            rng.uniform(10, 500, size=(3, len(samples)))
        # Kitome: present in ALL blanks but only sporadically in real samples (low-biomass signature).
        data.loc[['Ralstonia_pickettii', 'Bradyrhizobium_sp'], blanks] = rng.uniform(5, 20, size=(2, len(blanks)))
        sporadic = rng.uniform(0, 15, size=(2, len(samples)))
        sporadic[sporadic < 8] = 0
        data.loc[['Ralstonia_pickettii', 'Bradyrhizobium_sp'], samples] = sporadic
        table = data

    flagged = flag_contaminants(table, blanks, samples)
    print('Contaminant suspects (enriched in blanks):')
    print(flagged.round(2).to_string())
    print(f'\nOf {len(flagged)} flagged, {int(flagged["kitome_genus"].sum())} match a canonical kitome genus.')
