library(flowCore)
library(ggplot2)
library(dplyr)

# === CONFIGURATION ===
output_dir <- 'results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. SIMULATE CYTOF DATA WITH BEADS ===
cat('Generating simulated CyTOF data with beads...\n')

set.seed(42)
n_cells <- 10000
n_beads <- 500

# Time course (acquisition order)
time <- 1:(n_cells + n_beads)

# Simulate drift (increasing over time)
drift_factor <- 1 + 0.0001 * time

# Cell data with drift
cell_cd45 <- rnorm(n_cells, 500, 50) * drift_factor[1:n_cells]
cell_cd3 <- rnorm(n_cells, 300, 40) * drift_factor[1:n_cells]

# Bead data (constant expected, but affected by drift)
bead_ce140_expected <- 600
bead_eu151_expected <- 550

bead_ce140 <- rnorm(n_beads, bead_ce140_expected, 30) * drift_factor[(n_cells+1):(n_cells+n_beads)]
bead_eu151 <- rnorm(n_beads, bead_eu151_expected, 25) * drift_factor[(n_cells+1):(n_cells+n_beads)]

# Cells have low bead channel signal
cell_ce140 <- rnorm(n_cells, 10, 5)
cell_eu151 <- rnorm(n_cells, 8, 4)

# Combine into data frames
cells <- data.frame(
    Time = 1:n_cells,
    CD45 = cell_cd45,
    CD3 = cell_cd3,
    Ce140 = cell_ce140,
    Eu151 = cell_eu151,
    is_bead = FALSE
)

beads <- data.frame(
    Time = (n_cells+1):(n_cells+n_beads),
    CD45 = rnorm(n_beads, 10, 5),  # Beads have low marker signal
    CD3 = rnorm(n_beads, 8, 4),
    Ce140 = bead_ce140,
    Eu151 = bead_eu151,
    is_bead = TRUE
)

# Interleave beads throughout acquisition
bead_positions <- sort(sample(1:(n_cells + n_beads), n_beads))
all_data <- rbind(cells, beads)
all_data <- all_data[order(c(setdiff(1:(n_cells+n_beads), bead_positions), bead_positions)), ]
all_data$Time <- 1:nrow(all_data)

cat('Total events:', nrow(all_data), '\n')
cat('Beads:', sum(all_data$is_bead), '\n')

# === 2. IDENTIFY BEADS ===
cat('\nIdentifying bead events...\n')

# Beads identified by high Ce140 and Eu151
bead_threshold <- 200
detected_beads <- all_data$Ce140 > bead_threshold & all_data$Eu151 > bead_threshold

cat('Detected beads:', sum(detected_beads), '\n')
cat('True positive rate:', round(sum(detected_beads & all_data$is_bead) / sum(all_data$is_bead) * 100, 1), '%\n')

# === 3. CALCULATE NORMALIZATION FACTORS ===
cat('\nCalculating normalization factors...\n')

bead_data <- all_data[detected_beads, ]

# Reference values (expected bead intensities)
reference <- c(Ce140 = bead_ce140_expected, Eu151 = bead_eu151_expected)

# Calculate time-dependent factors
n_bins <- 10
bead_data$time_bin <- cut(bead_data$Time, breaks = n_bins, labels = FALSE)

norm_factors <- bead_data %>%
    group_by(time_bin) %>%
    summarise(
        Ce140_median = median(Ce140),
        Eu151_median = median(Eu151),
        time_center = mean(Time),
        .groups = 'drop'
    ) %>%
    mutate(
        Ce140_factor = reference['Ce140'] / Ce140_median,
        Eu151_factor = reference['Eu151'] / Eu151_median,
        combined_factor = (Ce140_factor + Eu151_factor) / 2
    )

print(norm_factors)

# === 4. APPLY NORMALIZATION ===
cat('\nApplying normalization...\n')

# Interpolate factors for all time points
all_data$time_bin <- cut(all_data$Time, breaks = n_bins, labels = FALSE)
all_data <- merge(all_data, norm_factors[, c('time_bin', 'combined_factor')], by = 'time_bin', all.x = TRUE)

# Apply correction
all_data$CD45_normalized <- all_data$CD45 / all_data$combined_factor
all_data$CD3_normalized <- all_data$CD3 / all_data$combined_factor

# === 5. VISUALIZATION ===
cat('\nGenerating plots...\n')

# Bead drift over time
p1 <- ggplot(bead_data, aes(x = Time, y = Ce140)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = 'loess', color = 'red', se = FALSE) +
    geom_hline(yintercept = reference['Ce140'], linetype = 'dashed', color = 'blue') +
    theme_bw() +
    labs(title = 'Bead Drift (Ce140)', x = 'Time', y = 'Intensity')
ggsave(file.path(output_dir, 'bead_drift.png'), p1, width = 10, height = 4)

# Before/after CD45
cells_only <- all_data[!all_data$is_bead, ]

before_after <- data.frame(
    Time = rep(cells_only$Time, 2),
    CD45 = c(cells_only$CD45, cells_only$CD45_normalized),
    Status = rep(c('Before', 'After'), each = nrow(cells_only))
)

p2 <- ggplot(before_after, aes(x = Time, y = CD45, color = Status)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_smooth(method = 'loess', se = FALSE) +
    theme_bw() +
    labs(title = 'CD45 Before/After Normalization', x = 'Time', y = 'Intensity')
ggsave(file.path(output_dir, 'normalization_effect.png'), p2, width = 10, height = 5)

# Distribution comparison
p3 <- ggplot(before_after, aes(x = CD45, fill = Status)) +
    geom_histogram(bins = 50, alpha = 0.6, position = 'identity') +
    theme_bw() +
    labs(title = 'CD45 Distribution', x = 'Intensity', y = 'Count')
ggsave(file.path(output_dir, 'distribution_comparison.png'), p3, width = 8, height = 5)

# === 6. REMOVE BEADS AND EXPORT ===
cat('\nExporting normalized data...\n')

final_data <- all_data[!detected_beads, c('Time', 'CD45_normalized', 'CD3_normalized')]
colnames(final_data) <- c('Time', 'CD45', 'CD3')

write.csv(final_data, file.path(output_dir, 'normalized_cells.csv'), row.names = FALSE)

cat('Final cell count:', nrow(final_data), '\n')
cat('Results saved to', output_dir, '\n')
