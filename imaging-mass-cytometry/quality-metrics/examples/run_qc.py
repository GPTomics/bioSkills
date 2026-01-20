import numpy as np
import pandas as pd
from pathlib import Path
from scipy import ndimage

# === CONFIGURATION ===
output_dir = Path('qc_results/')
output_dir.mkdir(exist_ok=True)

# === 1. SIMULATE IMC DATA ===
print('Generating simulated IMC data...')

np.random.seed(42)
height, width = 256, 256
n_channels = 5
channel_names = ['CD45', 'CD3', 'CD68', 'panCK', 'DNA']

image_stack = np.zeros((n_channels, height, width), dtype=np.float32)

tissue_mask = np.zeros((height, width), dtype=bool)
tissue_mask[30:226, 30:226] = True
tissue_mask[100:120, 100:156] = False

for i in range(n_channels):
    background = np.random.normal(50, 20, (height, width))
    signal = np.random.normal(200, 50, (height, width))
    image_stack[i] = np.where(tissue_mask, signal, background)
    image_stack[i] = np.clip(image_stack[i], 0, 255)

image_stack[2][50:55, 50:55] = 255
image_stack[3, 100:110, :] *= 0.5

print(f'Image shape: {image_stack.shape}')

# === 2. SIGNAL-TO-NOISE RATIO ===
print('\n=== Signal-to-Noise Ratio ===')

def calculate_snr(image, mask=None):
    if mask is None:
        mask = image > np.percentile(image, 10)
    signal = np.mean(image[mask])
    noise = np.std(image[~mask])
    return signal / noise if noise > 0 else np.inf

snr_results = {}
for i, name in enumerate(channel_names):
    snr = calculate_snr(image_stack[i], tissue_mask)
    snr_results[name] = snr
    status = 'PASS' if snr > 3 else 'WARN' if snr > 1.5 else 'FAIL'
    print(f'{name}: SNR = {snr:.2f} [{status}]')

# === 3. CHANNEL CORRELATION ===
print('\n=== Channel Correlation ===')

flat_data = image_stack.reshape(n_channels, -1)
corr_matrix = np.corrcoef(flat_data)
corr_df = pd.DataFrame(corr_matrix, index=channel_names, columns=channel_names)

print(corr_df.round(2))

high_corr = []
for i, ch1 in enumerate(channel_names):
    for j, ch2 in enumerate(channel_names):
        if i < j and corr_matrix[i, j] > 0.7:
            high_corr.append((ch1, ch2, corr_matrix[i, j]))

if high_corr:
    print('\nHigh correlations detected:')
    for ch1, ch2, corr in high_corr:
        print(f'  {ch1} - {ch2}: {corr:.2f}')

# === 4. TISSUE INTEGRITY ===
print('\n=== Tissue Integrity ===')

dna_channel = image_stack[channel_names.index('DNA')]
dna_threshold = np.percentile(dna_channel, 50)
tissue_detected = dna_channel > dna_threshold

total_pixels = dna_channel.size
tissue_pixels = np.sum(tissue_detected)
coverage = tissue_pixels / total_pixels

labeled, n_fragments = ndimage.label(tissue_detected)
fragment_sizes = ndimage.sum(tissue_detected, labeled, range(1, n_fragments + 1))
largest_fragment = np.max(fragment_sizes) if len(fragment_sizes) > 0 else 0
fragmentation = 1 - (largest_fragment / tissue_pixels) if tissue_pixels > 0 else 1

print(f'Tissue coverage: {coverage:.1%}')
print(f'Number of fragments: {n_fragments}')
print(f'Fragmentation score: {fragmentation:.2f}')
print(f'Status: {"PASS" if coverage > 0.3 and fragmentation < 0.5 else "WARN"}')

# === 5. ACQUISITION ARTIFACTS ===
print('\n=== Acquisition Artifacts ===')

artifacts = []
for i, name in enumerate(channel_names):
    channel = image_stack[i]

    saturated_pct = np.sum(channel >= 254) / channel.size
    if saturated_pct > 0.01:
        artifacts.append({'channel': name, 'issue': 'saturation', 'severity': f'{saturated_pct:.2%}'})

    hot_threshold = np.percentile(channel, 99.9) * 1.5
    hot_pct = np.sum(channel > hot_threshold) / channel.size
    if hot_pct > 0.001:
        artifacts.append({'channel': name, 'issue': 'hot_pixels', 'severity': f'{hot_pct:.3%}'})

    row_means = np.mean(channel, axis=1)
    row_cv = np.std(row_means) / np.mean(row_means)
    if row_cv > 0.2:
        artifacts.append({'channel': name, 'issue': 'striping', 'severity': f'{row_cv:.2f}'})

if artifacts:
    print('Artifacts detected:')
    for a in artifacts:
        print(f"  {a['channel']}: {a['issue']} ({a['severity']})")
else:
    print('No major artifacts detected')

# === 6. DYNAMIC RANGE ===
print('\n=== Dynamic Range ===')

for i, name in enumerate(channel_names):
    channel = image_stack[i]
    low, high = np.percentile(channel, [1, 99])
    utilized = (high - low) / 255
    status = 'OK' if utilized > 0.3 else 'LOW'
    print(f'{name}: {utilized:.1%} range used [{status}]')

# === 7. SIMULATE SEGMENTATION AND CHECK ===
print('\n=== Segmentation QC ===')

n_cells = 100
seg_mask = np.zeros((height, width), dtype=np.int32)
for cell_id in range(1, n_cells + 1):
    cx, cy = np.random.randint(40, 216, 2)
    radius = np.random.randint(5, 12)
    y, x = np.ogrid[:height, :width]
    cell_mask = ((x - cx)**2 + (y - cy)**2) <= radius**2
    seg_mask[cell_mask] = cell_id

from skimage.measure import regionprops
props = regionprops(seg_mask)
areas = [p.area for p in props]
eccentricities = [p.eccentricity for p in props]

print(f'Cells detected: {len(props)}')
print(f'Mean area: {np.mean(areas):.1f} pixels')
print(f'Area CV: {np.std(areas)/np.mean(areas):.2f}')
print(f'Elongated cells: {np.sum(np.array(eccentricities) > 0.9)/len(props):.1%}')

# === 8. QC SUMMARY ===
print('\n=== QC SUMMARY ===')

qc_results = {
    'mean_snr': np.mean(list(snr_results.values())),
    'min_snr': min(snr_results.values()),
    'tissue_coverage': coverage,
    'n_artifacts': len(artifacts),
    'n_cells': len(props)
}

passed = (qc_results['min_snr'] > 1.5 and qc_results['tissue_coverage'] > 0.3 and qc_results['n_artifacts'] <= 1)

print(f"Mean SNR: {qc_results['mean_snr']:.2f}")
print(f"Min SNR: {qc_results['min_snr']:.2f}")
print(f"Tissue coverage: {qc_results['tissue_coverage']:.1%}")
print(f"Artifacts: {qc_results['n_artifacts']}")
print(f"Cells: {qc_results['n_cells']}")
print(f'\nOVERALL: {"PASS" if passed else "REVIEW REQUIRED"}')

# === 9. EXPORT ===
pd.DataFrame([qc_results]).to_csv(output_dir / 'qc_summary.csv', index=False)
corr_df.to_csv(output_dir / 'channel_correlations.csv')

if artifacts:
    pd.DataFrame(artifacts).to_csv(output_dir / 'artifacts.csv', index=False)

print(f'\nResults saved to {output_dir}/')
