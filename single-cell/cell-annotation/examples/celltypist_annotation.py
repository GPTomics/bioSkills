import scanpy as sc
import celltypist
import matplotlib.pyplot as plt

adata = sc.read_h5ad('adata_processed.h5ad')

celltypist.models.download_models(model='Immune_All_Low.pkl')
model = celltypist.models.Model.load(model='Immune_All_Low.pkl')

predictions = celltypist.annotate(adata, model=model, majority_voting=True)
adata = predictions.to_adata()

adata.obs['cell_type'] = adata.obs['majority_voting']
adata.obs['annotation_confidence'] = adata.obs['conf_score']

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color='cell_type', ax=axes[0], show=False, title='Cell Types')
sc.pl.umap(adata, color='annotation_confidence', ax=axes[1], show=False,
           title='Confidence Score', cmap='viridis')
plt.tight_layout()
plt.savefig('celltypist_annotation.png', dpi=150, bbox_inches='tight')
plt.close()

confidence_threshold = 0.5
adata.obs['high_confidence'] = adata.obs['annotation_confidence'] > confidence_threshold
low_conf_cells = adata.obs[~adata.obs['high_confidence']]
print(f'Low confidence cells: {len(low_conf_cells)} ({len(low_conf_cells)/len(adata)*100:.1f}%)')

cell_type_counts = adata.obs['cell_type'].value_counts()
cell_type_counts.to_csv('cell_type_counts.csv')

adata.write('adata_annotated.h5ad')
