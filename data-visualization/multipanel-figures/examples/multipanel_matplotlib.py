'''Create multi-panel publication figures with matplotlib'''

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

# --- ALTERNATIVE: Use real data ---
# Load your analysis results:
# pca_df = pd.read_csv('pca_results.csv')
# de_results = pd.read_csv('de_results.csv')

np.random.seed(42)

# Create example data for demonstration panels
# Panel A: PCA
pca_x = np.concatenate([np.random.normal(-2, 0.8, 20), np.random.normal(2, 0.8, 20)])
pca_y = np.random.normal(0, 1, 40)
conditions = ['Control'] * 20 + ['Treatment'] * 20

# Panel B: Volcano
log2fc = np.random.normal(0, 1.5, 1000)
pvalue = 10 ** (-np.abs(log2fc) * np.random.uniform(0, 3, 1000))
sig = np.where((np.abs(log2fc) > 1) & (pvalue < 0.05),
               np.where(log2fc > 0, 'Up', 'Down'), 'NS')

# Panel C: Expression boxplot data
genes = ['Gene1', 'Gene2', 'Gene3', 'Gene4']
expr_ctrl = [np.random.normal(m, 0.5, 10) for m in [5, 8, 3, 6]]
expr_treat = [np.random.normal(m, 0.5, 10) for m in [7, 6, 5, 4]]

# Panel D: MA plot
base_mean = 10 ** np.random.uniform(0, 4, 1000)
ma_log2fc = np.random.normal(0, 1.5, 1000)

# Color schemes matching R examples
colors = {'Control': '#4DBBD5', 'Treatment': '#E64B35'}
sig_colors = {'Up': '#E64B35', 'Down': '#4DBBD5', 'NS': 'lightgray'}

# Figure 1: 2x2 grid with GridSpec
fig = plt.figure(figsize=(10, 8))
gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)

# Panel A: PCA
ax_pca = fig.add_subplot(gs[0, 0])
for cond in ['Control', 'Treatment']:
    mask = np.array(conditions) == cond
    ax_pca.scatter(pca_x[mask], pca_y[mask], c=colors[cond], label=cond, alpha=0.7, s=50)
ax_pca.set_xlabel('PC1 (45%)')
ax_pca.set_ylabel('PC2 (20%)')
ax_pca.set_title('PCA')
ax_pca.legend()

# Panel B: Volcano
ax_volcano = fig.add_subplot(gs[0, 1])
for s in ['NS', 'Down', 'Up']:  # Plot NS first so significant on top
    mask = sig == s
    ax_volcano.scatter(log2fc[mask], -np.log10(pvalue[mask]),
                       c=sig_colors[s], alpha=0.5, s=10, label=s)
ax_volcano.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1)
ax_volcano.axvline(-1, color='gray', linestyle='--', linewidth=1)
ax_volcano.axvline(1, color='gray', linestyle='--', linewidth=1)
ax_volcano.set_xlabel(r'$log_2$ Fold Change')
ax_volcano.set_ylabel(r'$-log_{10}$ p-value')
ax_volcano.set_title('Volcano Plot')

# Panel C: Boxplot
ax_box = fig.add_subplot(gs[1, 0])
positions = np.arange(len(genes)) * 2
width = 0.35
bp1 = ax_box.boxplot(expr_ctrl, positions=positions - width/2, widths=width,
                      patch_artist=True, boxprops=dict(facecolor=colors['Control']))
bp2 = ax_box.boxplot(expr_treat, positions=positions + width/2, widths=width,
                      patch_artist=True, boxprops=dict(facecolor=colors['Treatment']))
ax_box.set_xticks(positions)
ax_box.set_xticklabels(genes)
ax_box.set_ylabel('Expression (log2)')
ax_box.set_title('Top DE Genes')
ax_box.legend([bp1['boxes'][0], bp2['boxes'][0]], ['Control', 'Treatment'])

# Panel D: MA plot
ax_ma = fig.add_subplot(gs[1, 1])
ax_ma.scatter(np.log10(base_mean), ma_log2fc, alpha=0.3, s=10, c='gray')
ax_ma.axhline(0, color='red', linewidth=1)
ax_ma.set_xlabel('log10(Mean Expression)')
ax_ma.set_ylabel('log2 Fold Change')
ax_ma.set_title('MA Plot')

# Add panel labels
# Position relative to axes: (-0.1, 1.1) puts label above and left of plot
for ax, label in zip([ax_pca, ax_volcano, ax_box, ax_ma], 'ABCD'):
    ax.text(-0.1, 1.1, label, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top')

plt.suptitle('Figure 1: RNA-seq Analysis Summary', fontsize=14, y=1.02)
plt.savefig('multipanel_2x2.pdf', bbox_inches='tight')
plt.savefig('multipanel_2x2.png', dpi=300, bbox_inches='tight')
plt.close()

# Figure 2: Complex layout (large plot spanning rows)
fig = plt.figure(figsize=(12, 8))
gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

# Panel A: Large PCA spanning 2 rows, 2 columns
ax_pca = fig.add_subplot(gs[:, :2])
for cond in ['Control', 'Treatment']:
    mask = np.array(conditions) == cond
    ax_pca.scatter(pca_x[mask], pca_y[mask], c=colors[cond], label=cond, alpha=0.7, s=80)
ax_pca.set_xlabel('PC1 (45%)', fontsize=12)
ax_pca.set_ylabel('PC2 (20%)', fontsize=12)
ax_pca.set_title('PCA', fontsize=14)
ax_pca.legend(fontsize=10)
ax_pca.text(-0.05, 1.02, 'A', transform=ax_pca.transAxes,
            fontsize=16, fontweight='bold', va='bottom')

# Panel B: Volcano (top right)
ax_volcano = fig.add_subplot(gs[0, 2])
for s in ['NS', 'Down', 'Up']:
    mask = sig == s
    ax_volcano.scatter(log2fc[mask], -np.log10(pvalue[mask]),
                       c=sig_colors[s], alpha=0.5, s=8)
ax_volcano.set_title('Volcano', fontsize=12)
ax_volcano.text(-0.15, 1.05, 'B', transform=ax_volcano.transAxes,
                fontsize=14, fontweight='bold', va='bottom')

# Panel C: MA plot (bottom right)
ax_ma = fig.add_subplot(gs[1, 2])
ax_ma.scatter(np.log10(base_mean), ma_log2fc, alpha=0.3, s=8, c='gray')
ax_ma.axhline(0, color='red', linewidth=1)
ax_ma.set_title('MA Plot', fontsize=12)
ax_ma.text(-0.15, 1.05, 'C', transform=ax_ma.transAxes,
           fontsize=14, fontweight='bold', va='bottom')

plt.savefig('multipanel_complex.pdf', bbox_inches='tight')
plt.savefig('multipanel_complex.png', dpi=300, bbox_inches='tight')
plt.close()

print('Multi-panel figures saved: multipanel_2x2.pdf/png, multipanel_complex.pdf/png')
