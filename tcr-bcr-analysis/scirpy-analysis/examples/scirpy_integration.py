'''Single-cell paired TCR/BCR analysis with scirpy on the awkward-array AIRR model.

Runs end-to-end against scirpy's bundled wu2020_3k dataset (a MuData with gex + airr
modalities), so no external files are needed. Demonstrates the modern data model:
AIRR lives in obsm['airr'] and is read via get.airr after pp.index_chains - NOT via
legacy per-chain obs columns like IR_VJ_1_junction_aa.
'''
# Reference: scirpy 0.24+ scanpy 1.10+ | Verify API if version differs

import warnings
import matplotlib
import scirpy as ir
import scanpy as sc

matplotlib.use('Agg')  # headless: build figure objects without opening a window
warnings.filterwarnings('ignore')

# wu2020_3k returns a MuData({'gex': ..., 'airr': ...}); real projects build this from
# ir.io.read_10x_vdj(...) (returns an AnnData) wrapped with the GEX AnnData in a MuData.
mdata = ir.datasets.wu2020_3k()
print('modalities:', list(mdata.mod.keys()), '| cells:', mdata.n_obs)

# index_chains is REQUIRED before any QC/clonotyping; it splits ragged per-cell chains
# into VJ/VDJ and flags multichain cells, writing obsm['chain_indices'].
ir.pp.index_chains(mdata)

# chain_qc writes airr:receptor_type / airr:receptor_subtype / airr:chain_pairing.
ir.tl.chain_qc(mdata)
print('\nchain pairing:')
print(mdata.obs['airr:chain_pairing'].value_counts())

# Always drop likely doublets (multichain, and TCR+BCR ambiguous cells).
# Orphan/extra cells are kept here; deleting all of them would preferentially remove
# small clones and inflate apparent clonal expansion, so match that choice to the question.
keep = mdata.obs['airr:receptor_type'].isin(['TCR', 'BCR']) & ~mdata.obs['airr:chain_pairing'].isin(['multichain'])
sub = mdata[keep].copy()
print(f'\nretained {sub.n_obs} cells after dropping doublets')

# TCR clonotypes: exact CDR3-nucleotide identity is correct for T cells (no somatic hypermutation).
# ir_dist metric/sequence MUST match the subsequent define_* call or the cache silently mismatches.
ir.pp.ir_dist(sub, metric='identity', sequence='nt', cutoff=0)
ir.tl.define_clonotypes(sub, receptor_arms='all', dual_ir='primary_only')
print('\nTCR clonotypes (identity):', sub.obs['airr:clone_id'].nunique())

# BCR-recommended settings (shown here on TCR data purely to exercise the call): identity
# clonotyping SHATTERS hypermutated B-cell lineages into fake singletons, so BCR needs
# distance clustering on nucleotides within same-V/same-J partitions. cutoff ~15 nt ~= 85% identity.
ir.pp.ir_dist(sub, metric='normalized_hamming', sequence='nt', cutoff=15)
ir.tl.define_clonotype_clusters(sub, sequence='nt', metric='normalized_hamming', receptor_arms='all', dual_ir='any', same_v_gene=True, same_j_gene=True)
print('distance clusters (normalized_hamming):', sub.obs['airr:cc_nt_normalized_hamming'].nunique())

# Clonal expansion: bins each cell singleton / 2 / >= 3 via breakpoints=(1, 2).
# target_col is resolved within the airr modality, so use the bare name 'clone_id'.
ir.tl.clonal_expansion(sub, target_col='clone_id')
print('\nclonal expansion:')
print(sub.obs['airr:clonal_expansion'].value_counts())

# Diversity and overlap depend entirely on the QC filter + clonotype definition above;
# report them alongside the number, and compare overlap only at equal sequencing depth.
div = ir.tl.alpha_diversity(sub, groupby='gex:sample', target_col='clone_id', metric='normalized_shannon_entropy', inplace=False)
print('\nnormalized Shannon per sample:')
print(div.head())

# Read AIRR fields with the get accessor - obsm['airr'] cannot be indexed like obs columns.
junction_vj = ir.get.airr(sub, 'junction_aa', 'VJ_1')
print('\nexample VJ junction_aa:', junction_vj.dropna().iloc[0])

# Overlay clonality on the transcriptome: cluster on GEX independently, then color the UMAP
# by a clonality column pushed into the GEX modality. Never cluster cells on receptor sequence.
gex = sub['gex']
sc.pp.normalize_total(gex, target_sum=1e4)
sc.pp.log1p(gex)
sc.pp.pca(gex, n_comps=30)
sc.pp.neighbors(gex)
sc.tl.umap(gex)
gex.obs['clonal_expansion'] = sub.obs['airr:clonal_expansion'].values
ax = sc.pl.umap(gex, color='clonal_expansion', show=False)  # figure built in memory, not saved
print('\nUMAP overlay of clonal expansion computed')

print('\nanalysis complete')
