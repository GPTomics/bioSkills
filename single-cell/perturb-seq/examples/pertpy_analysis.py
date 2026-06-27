'''Perturb-seq analysis with pertpy: mixture guide assignment, Mixscape escaper removal, E-distance.'''
# Reference: pertpy 0.9+, scanpy 1.10+ | Verify API if version differs
import scanpy as sc
import pertpy as pt

mdata = pt.dt.papalexi_2021()
adata = mdata.mod['rna']

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50)

# Guide assignment is a mixture problem, not a flat threshold: ambient contamination is biased to abundant guides
gdo = mdata.mod['gdo']
gdo.layers['counts'] = gdo.X.copy()
ga = pt.pp.GuideAssignment()
ga.assign_mixture_model(gdo, assigned_guides_key='assigned_guide')
print(gdo.obs['assigned_guide'].value_counts().head())

# Mixscape removes non-perturbed escapers before any DE (assignment != effective perturbation)
ms = pt.tl.Mixscape()
ms.perturbation_signature(adata, pert_key='perturbation', control='NT', n_neighbors=20)
ms.mixscape(adata, pert_key='gene_target', control='NT', layer='X_pert')   # pert_key renamed from labels
# An all-NP target confounds no-phenotype with no-editing: report the perturbed fraction, do not call the gene dead
print(adata.obs['mixscape_class_global'].value_counts())

# E-distance is the modern effect-size; pin the embedding and metric (sqeuclidean vs euclidean default changed)
dist = pt.tl.Distance(metric='edistance', obsm_key='X_pca')
pairwise = dist.pairwise(adata, groupby='gene_target')
print(pairwise.iloc[:5, :5])

# E-test: smallest p ~ 1/(n_perms+1), crushed by multiple testing across many perturbations
etest = pt.tl.DistanceTest('edistance', n_perms=1000)
results = etest(adata, groupby='gene_target', contrast='NT')
print(results.sort_values('pvalue').head())
