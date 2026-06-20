#!/usr/bin/env python3
"""Formatted display table (great_tables) + a gene-symbol-safe Excel/CSV supplement."""
# Reference: great_tables 0.13+, pandas 2.2+, openpyxl 3.1+ | Verify API if version differs

import pandas as pd
from great_tables import GT
from openpyxl import Workbook

de = pd.DataFrame({
    'gene': ['SEPT9', 'MARCH1', 'TP53', 'BRCA1', 'DEC1'],
    'log2fc': [2.13, -1.04, 0.88, -2.50, 1.37],
    'padj': [1e-12, 3e-4, 0.04, 2e-8, 7e-3],
})


def display_table(df, path):
    """Formatted display table to PNG: scientific p-values, 2-decimal fold changes."""
    gt = (GT(df, rowname_col='gene')
          .fmt_number(columns='log2fc', decimals=2)
          .fmt_scientific(columns='padj')
          .tab_header(title='Differential Expression', subtitle='Top differentially expressed genes'))
    gt.save(path)


def write_gene_safe_xlsx(df, gene_col, path):
    """Write .xlsx forcing the gene column to Excel text so SEPT9 does not become a date."""
    wb = Workbook()
    ws = wb.active
    ws.append(list(df.columns))
    for record in df.itertuples(index=False):
        ws.append(list(record))
    gene_idx = list(df.columns).index(gene_col) + 1
    for column in ws.iter_cols(min_col=gene_idx, max_col=gene_idx, min_row=2):
        for cell in column:
            cell.number_format = '@'
    wb.save(path)


if __name__ == '__main__':
    display_table(de, 'de_table.png')
    # CSV is the robust machine-readable supplement (consumer imports the gene column as text)
    de.to_csv('de_supplement.csv', index=False)
    write_gene_safe_xlsx(de, 'gene', 'de_supplement.xlsx')
    print('Wrote de_table.png, de_supplement.csv, de_supplement.xlsx')
