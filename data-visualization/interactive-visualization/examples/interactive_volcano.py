import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np

np.random.seed(42)
n_genes = 1000

df = pd.DataFrame({
    'gene': [f'Gene{i}' for i in range(n_genes)],
    'log2FoldChange': np.random.normal(0, 1.5, n_genes),
    'pvalue': 10 ** (-np.random.uniform(0, 10, n_genes)),
    'baseMean': 10 ** np.random.uniform(1, 4, n_genes)
})
df['padj'] = np.minimum(df['pvalue'] * n_genes / np.arange(1, n_genes + 1), 1)
df['neg_log_pval'] = -np.log10(df['pvalue'])

df['significance'] = 'Not Significant'
df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] > 1), 'significance'] = 'Up'
df.loc[(df['padj'] < 0.05) & (df['log2FoldChange'] < -1), 'significance'] = 'Down'

color_map = {'Up': '#E64B35', 'Down': '#4DBBD5', 'Not Significant': '#999999'}

fig = px.scatter(df, x='log2FoldChange', y='neg_log_pval',
                 color='significance', hover_name='gene',
                 hover_data={'baseMean': ':.1f', 'padj': ':.2e',
                            'log2FoldChange': ':.2f', 'significance': False},
                 color_discrete_map=color_map,
                 title='Interactive Volcano Plot',
                 labels={'log2FoldChange': 'Log2 Fold Change',
                        'neg_log_pval': '-Log10 P-value'})

fig.add_hline(y=-np.log10(0.05), line_dash='dash', line_color='grey',
              annotation_text='p=0.05')
fig.add_vline(x=-1, line_dash='dash', line_color='grey')
fig.add_vline(x=1, line_dash='dash', line_color='grey')

fig.update_traces(marker=dict(size=6, opacity=0.7))
fig.update_layout(
    legend_title_text='Significance',
    hovermode='closest',
    template='plotly_white'
)

fig.write_html('volcano_interactive.html')
print(f"Saved: volcano_interactive.html")
print(f"Up-regulated: {sum(df['significance'] == 'Up')}")
print(f"Down-regulated: {sum(df['significance'] == 'Down')}")
