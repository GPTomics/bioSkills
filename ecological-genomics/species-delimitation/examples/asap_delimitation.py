# Reference: biopython 1.83+, numpy 1.26+, scipy 1.12+ | Verify API if version differs
from Bio import AlignIO
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
from scipy.spatial.distance import squareform
import numpy as np
import matplotlib.pyplot as plt

# --- Step 1: Load aligned sequences ---
alignment = AlignIO.read('aligned_sequences.fasta', 'fasta')
n_seqs = len(alignment)
seq_length = alignment.get_alignment_length()
names = [record.id for record in alignment]

print(f'Sequences: {n_seqs}')
print(f'Alignment length: {seq_length}')

# --- Step 2: Compute pairwise p-distances ---
# p-distance: proportion of differing sites (uncorrected)
dist_matrix = np.zeros((n_seqs, n_seqs))
for i in range(n_seqs):
    for j in range(i + 1, n_seqs):
        seq_i = str(alignment[i].seq).upper()
        seq_j = str(alignment[j].seq).upper()

        valid = 0
        diff = 0
        for si, sj in zip(seq_i, seq_j):
            if si in 'ACGT' and sj in 'ACGT':
                valid += 1
                if si != sj:
                    diff += 1

        p_dist = diff / valid if valid > 0 else 0.0
        dist_matrix[i][j] = p_dist
        dist_matrix[j][i] = p_dist

# --- Step 3: Identify barcode gap ---
# Collect all pairwise distances
condensed = squareform(dist_matrix)

# Distance histogram to visualize barcode gap
# The gap between intra- and inter-specific distances is the barcode gap
# COI: typically 1-3% intraspecific, >3% interspecific for animals
plt.figure(figsize=(10, 5))
plt.hist(condensed, bins=100, color='steelblue', edgecolor='white', alpha=0.8)
plt.xlabel('Pairwise p-distance')
plt.ylabel('Frequency')
plt.title('Pairwise Distance Distribution (Barcode Gap)')
plt.axvline(x=0.03, color='red', linestyle='--', label='3% threshold')
plt.legend()
plt.tight_layout()
plt.savefig('barcode_gap_histogram.pdf')
plt.close()

# --- Step 4: Hierarchical clustering ---
Z = linkage(condensed, method='average')

# --- Step 5: Scan distance thresholds ---
# Test a range of thresholds to find the barcode gap region
# Stable plateaus in the number-of-groups curve indicate natural breaks
thresholds = np.arange(0.005, 0.15, 0.005)
n_groups_per_threshold = []

print('\nThreshold scan:')
for t in thresholds:
    clusters = fcluster(Z, t=t, criterion='distance')
    n_groups = len(set(clusters))
    n_groups_per_threshold.append(n_groups)
    print(f'  {t:.3f}: {n_groups} groups')

# Plot threshold vs number of groups
plt.figure(figsize=(10, 5))
plt.plot(thresholds, n_groups_per_threshold, 'o-', color='steelblue', markersize=4)
plt.xlabel('Distance threshold')
plt.ylabel('Number of groups')
plt.title('Species Count vs Distance Threshold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('threshold_scan.pdf')
plt.close()

# --- Step 6: Apply optimal threshold ---
# 0.03 (3%): standard COI barcode gap for animals
# Adjust based on barcode gap histogram and threshold scan
# For 12S: use 0.02 (2%) due to lower divergence rates
# For ITS fungi: use 0.03-0.05 depending on genus
optimal_threshold = 0.03
clusters = fcluster(Z, t=optimal_threshold, criterion='distance')
n_species = len(set(clusters))
print(f'\nOptimal threshold: {optimal_threshold}')
print(f'Putative species: {n_species}')

# --- Step 7: Assign species ---
species_assignments = {}
for name, cluster_id in zip(names, clusters):
    species_assignments.setdefault(cluster_id, []).append(name)

print('\nSpecies assignments:')
for sp_id in sorted(species_assignments.keys()):
    members = species_assignments[sp_id]
    print(f'  Species {sp_id} ({len(members)} individuals): {", ".join(members[:5])}{"..." if len(members) > 5 else ""}')

# --- Step 8: Intra- vs inter-specific distances ---
intra_dists = []
inter_dists = []
for i in range(n_seqs):
    for j in range(i + 1, n_seqs):
        if clusters[i] == clusters[j]:
            intra_dists.append(dist_matrix[i][j])
        else:
            inter_dists.append(dist_matrix[i][j])

if intra_dists:
    print(f'\nIntraspecific: mean={np.mean(intra_dists):.4f}, max={np.max(intra_dists):.4f}')
if inter_dists:
    print(f'Interspecific: mean={np.mean(inter_dists):.4f}, min={np.min(inter_dists):.4f}')

# Barcode gap present if max(intra) < min(inter)
if intra_dists and inter_dists:
    gap = np.min(inter_dists) - np.max(intra_dists)
    print(f'Barcode gap width: {gap:.4f} ({"present" if gap > 0 else "absent/overlap"})')

# --- Step 9: Dendrogram ---
plt.figure(figsize=(12, max(8, n_seqs * 0.2)))
dendrogram(Z, labels=names, orientation='right', leaf_font_size=6,
           color_threshold=optimal_threshold)
plt.axvline(x=optimal_threshold, color='red', linestyle='--',
            label=f'Threshold = {optimal_threshold}')
plt.xlabel('Distance')
plt.title('UPGMA Clustering for Species Delimitation')
plt.legend()
plt.tight_layout()
plt.savefig('species_dendrogram.pdf')
plt.close()

# --- Step 10: Export ---
with open('species_assignments.tsv', 'w') as f:
    f.write('individual\tspecies_id\n')
    for name, cluster_id in zip(names, clusters):
        f.write(f'{name}\t{cluster_id}\n')

print('\nResults written to species_assignments.tsv')
