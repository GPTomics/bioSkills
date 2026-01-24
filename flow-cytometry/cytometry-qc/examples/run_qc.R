library(flowCore)
library(ggplot2)
library(dplyr)

# === CONFIGURATION ===
output_dir <- 'qc_results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. SIMULATE FCS DATA ===
cat('Generating simulated flow cytometry data...\n')

set.seed(42)
n_events <- 20000

# Time course
time <- 1:n_events

# Simulate flow rate instability (clog at time 8000-9000)
flow_rate_factor <- rep(1, n_events)
flow_rate_factor[8000:9000] <- 0.3

# Simulate signal drift (gradual increase)
drift_factor <- 1 + 0.00005 * time

# Generate FSC/SSC
fsc_a <- rnorm(n_events, 100000, 15000) * drift_factor
ssc_a <- rnorm(n_events, 50000, 10000) * drift_factor

# Generate markers
cd45 <- rnorm(n_events, 500, 80) * drift_factor
cd3 <- rnorm(n_events, 300, 60) * drift_factor
cd4 <- rnorm(n_events, 200, 50) * drift_factor

# Add some dead cells (high viability dye)
viability <- rnorm(n_events, 100, 30)
dead_idx <- sample(1:n_events, n_events * 0.08)
viability[dead_idx] <- rnorm(length(dead_idx), 500, 100)

# Add margin events
margin_idx <- sample(1:n_events, n_events * 0.02)
fsc_a[margin_idx] <- 262143  # Max value

data <- data.frame(
    Time = time,
    FSC_A = pmax(fsc_a, 0),
    SSC_A = pmax(ssc_a, 0),
    CD45 = pmax(cd45, 0),
    CD3 = pmax(cd3, 0),
    CD4 = pmax(cd4, 0),
    Viability = pmax(viability, 0)
)

cat('Total events:', nrow(data), '\n')

# === 2. FLOW RATE CHECK ===
cat('\n=== Flow Rate Analysis ===\n')

n_bins <- 50
data$time_bin <- cut(data$Time, breaks = n_bins, labels = FALSE)
flow_rate <- data %>%
    group_by(time_bin) %>%
    summarise(events = n(), .groups = 'drop')

flow_cv <- sd(flow_rate$events) / mean(flow_rate$events) * 100
cat('Flow rate CV:', round(flow_cv, 1), '%\n')
# CV<20%: Standard flow rate stability threshold. Higher CV indicates clogs or bubbles.
cat('Status:', ifelse(flow_cv < 20, 'PASS', 'FAIL'), '\n')

# Plot flow rate
p1 <- ggplot(flow_rate, aes(x = time_bin, y = events)) +
    geom_line() +
    geom_hline(yintercept = mean(flow_rate$events), linetype = 'dashed', color = 'blue') +
    geom_hline(yintercept = mean(flow_rate$events) * 0.7, linetype = 'dashed', color = 'red') +
    theme_bw() +
    labs(title = 'Flow Rate Over Acquisition', x = 'Time Bin', y = 'Events per Bin')
ggsave(file.path(output_dir, 'flow_rate.png'), p1, width = 10, height = 4)

# === 3. SIGNAL DRIFT CHECK ===
cat('\n=== Signal Drift Analysis ===\n')

drift_check <- data %>%
    group_by(time_bin) %>%
    summarise(
        CD45_median = median(CD45),
        CD3_median = median(CD3),
        .groups = 'drop'
    )

# Linear trend
cd45_trend <- lm(CD45_median ~ time_bin, data = drift_check)
cd45_pct_change <- (predict(cd45_trend, newdata = data.frame(time_bin = n_bins)) -
                   predict(cd45_trend, newdata = data.frame(time_bin = 1))) /
                   predict(cd45_trend, newdata = data.frame(time_bin = 1)) * 100

cat('CD45 drift:', round(cd45_pct_change, 1), '%\n')
# drift<10%: Signal should not change >10% during acquisition. Indicates instrument instability.
cat('Status:', ifelse(abs(cd45_pct_change) < 10, 'PASS', 'WARNING'), '\n')

# Plot signal drift
p2 <- ggplot(drift_check, aes(x = time_bin, y = CD45_median)) +
    geom_point() +
    geom_smooth(method = 'lm', color = 'red', se = FALSE) +
    theme_bw() +
    labs(title = 'CD45 Signal Over Time', x = 'Time Bin', y = 'Median Intensity')
ggsave(file.path(output_dir, 'signal_drift.png'), p2, width = 10, height = 4)

# === 4. MARGIN EVENTS ===
cat('\n=== Margin Events ===\n')

max_val <- 262143
margin_mask <- data$FSC_A >= max_val * 0.99 | data$FSC_A <= 0

cat('Margin events:', sum(margin_mask), '(', round(mean(margin_mask) * 100, 2), '%)\n')
# <5% margin: Events at detector limits indicate saturation. >5% suggests voltage issues.
cat('Status:', ifelse(mean(margin_mask) < 0.05, 'PASS', 'WARNING'), '\n')

# === 5. DEAD CELL EXCLUSION ===
cat('\n=== Dead Cell Analysis ===\n')

viability_threshold <- quantile(data$Viability, 0.9)
dead_mask <- data$Viability > viability_threshold

cat('Dead cells:', sum(dead_mask), '(', round(mean(dead_mask) * 100, 1), '%)\n')
# <15% dead: Acceptable for most samples. >15% may indicate sample handling issues.
cat('Status:', ifelse(mean(dead_mask) < 0.15, 'PASS', 'WARNING'), '\n')

# Viability distribution
p3 <- ggplot(data, aes(x = Viability, fill = dead_mask)) +
    geom_histogram(bins = 100, alpha = 0.7) +
    geom_vline(xintercept = viability_threshold, linetype = 'dashed', color = 'red') +
    scale_fill_manual(values = c('gray50', 'red'), labels = c('Live', 'Dead')) +
    theme_bw() +
    labs(title = 'Viability Distribution', x = 'Viability Dye Intensity', y = 'Count', fill = 'Status')
ggsave(file.path(output_dir, 'viability.png'), p3, width = 8, height = 5)

# === 6. APPLY QC FILTERS ===
cat('\n=== Applying QC Filters ===\n')

# Combined filter
qc_pass <- !margin_mask & !dead_mask
cat('Events passing all QC:', sum(qc_pass), '(', round(mean(qc_pass) * 100, 1), '%)\n')

clean_data <- data[qc_pass, ]

# === 7. QC SUMMARY ===
cat('\n=== QC SUMMARY ===\n')
cat('Original events:', nrow(data), '\n')
cat('After margin removal:', sum(!margin_mask), '\n')
cat('After dead cell removal:', sum(!dead_mask), '\n')
cat('Final clean events:', nrow(clean_data), '\n')
cat('Total removed:', nrow(data) - nrow(clean_data), '(',
    round((1 - nrow(clean_data)/nrow(data)) * 100, 1), '%)\n')

# === 8. EXPORT ===
write.csv(clean_data, file.path(output_dir, 'qc_passed_data.csv'), row.names = FALSE)

# QC report
sink(file.path(output_dir, 'qc_report.txt'))
cat('=== FLOW CYTOMETRY QC REPORT ===\n\n')
cat('Date:', Sys.Date(), '\n')
cat('Total events:', nrow(data), '\n\n')
cat('--- Flow Rate ---\n')
cat('CV:', round(flow_cv, 1), '%\n')
cat('Status:', ifelse(flow_cv < 20, 'PASS', 'FAIL'), '\n\n')
cat('--- Signal Drift ---\n')
cat('CD45 change:', round(cd45_pct_change, 1), '%\n')
cat('Status:', ifelse(abs(cd45_pct_change) < 10, 'PASS', 'WARNING'), '\n\n')
cat('--- Event Filtering ---\n')
cat('Margin events removed:', sum(margin_mask), '\n')
cat('Dead cells removed:', sum(dead_mask), '\n')
cat('Final events:', nrow(clean_data), '\n')
sink()

cat('\nResults saved to', output_dir, '\n')
