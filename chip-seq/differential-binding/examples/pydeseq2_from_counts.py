# Reference: PyDESeq2 0.4+ | Verify API if version differs
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)
metadata = pd.DataFrame(
    {'condition': ['treated'] * 3 + ['control'] * 3},
    index=counts.columns
)

inference = DefaultInference(n_cpus=4)
# PyDESeq2 expects samples as rows, features as columns — transpose peaks-by-samples matrix
dds = DeseqDataSet(counts=counts.T, metadata=metadata, design='~ condition',
                   refit_cooks=True, inference=inference)
dds.deseq2()

ds = DeseqStats(dds, contrast=['condition', 'treated', 'control'], inference=inference)
ds.summary()
results = ds.results_df

results['peak_id'] = results.index
results['log2fc'] = results['log2FoldChange']
results['significant'] = results['padj'].lt(0.05).map({True: 'TRUE', False: 'FALSE'})
results[['peak_id', 'log2fc', 'pvalue', 'padj', 'significant']].to_csv(
    'differential.tsv', sep='\t', index=False
)
