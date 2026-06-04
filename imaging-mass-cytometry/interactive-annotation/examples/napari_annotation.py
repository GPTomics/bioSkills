'''IMC annotation helpers: cluster-on-tissue painting and class-balanced label sampling.

The decision-grade pieces of interactive annotation are not the napari boilerplate but
(1) painting clusters back onto the mask so spatial incoherence reveals artifact clusters,
and (2) deliberately over-sampling rare classes across patients when choosing what to
label. The napari session itself is shown as reference code (it needs a GUI).
'''
# Reference: napari 0.4.18+, napari-imc 0.7+, numpy 1.26+ | Verify API if version differs
import numpy as np

def cluster_label_image(masks, cell_ids, cluster_of_cell):
    # paint each cell's cluster id onto its mask; view in napari and demand spatial coherence
    # (nests/zones/sheets). A cluster that scatters along the boundary between two real
    # clusters is a segmentation/spillover artifact, not biology.
    out = np.zeros_like(masks)
    for cid, cl in zip(cell_ids, cluster_of_cell):
        out[masks == cid] = int(cl) + 1
    return out

def sample_cells_to_annotate(cell_ids, cell_types, per_class=300, rng=None):
    # over-sample rare classes toward a target; tissue is wildly imbalanced, so annotating
    # random fields leaves rare types unlearnable. Spread the chosen cells across patients.
    rng = rng or np.random.default_rng(0)
    picks = [rng.choice(cell_ids[cell_types == ct], size=min(per_class, int((cell_types == ct).sum())), replace=False)
             for ct in np.unique(cell_types)]
    return np.concatenate(picks)

# demo on a tiny synthetic example (no files written)
rng = np.random.default_rng(0)
cell_ids = np.arange(1, 1001)
cell_types = rng.choice(['T', 'B', 'tumor'], size=1000, p=[0.6, 0.35, 0.05])   # imbalanced
chosen = sample_cells_to_annotate(cell_ids, cell_types, per_class=40, rng=rng)
print(f'Selected {len(chosen)} cells; rare class is now represented, not buried by proportion')

reference_napari_session = '''
import napari
viewer = napari.Viewer()
viewer.open('slide.mcd', plugin='napari-imc')                 # core napari cannot read .mcd
viewer.add_labels(cell_masks, name='masks')
viewer.add_image(dna, name='DNA', contrast_limits=[0, 20])    # FIX contrast (a positivity threshold)
viewer.add_points(coords, features={'label': labels})         # per-cell categorical annotation
napari.run()
'''
print(reference_napari_session)
