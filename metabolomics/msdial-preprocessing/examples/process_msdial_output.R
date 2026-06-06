# Reference: R 4.3+ (base only) | Verify API if version differs
# Parses an MS-DIAL alignment-result export, splits metadata from per-sample
# intensities, and filters honestly on Fill% and MS/MS support. Generates a tiny
# synthetic MS-DIAL-style table (real exports carry 4 header rows above the column
# header) so the script runs with no external file. All output goes to tempdir().

make_synthetic_export <- function(path, n_features = 60, n_samples = 8) {
    set.seed(42)
    tags <- sample(c('Metabolite', 'Lipid', 'Suggested', 'Unknown'), n_features,
                   replace = TRUE, prob = c(0.25, 0.15, 0.2, 0.4))
    has_msms <- ifelse(tags %in% c('Metabolite', 'Lipid'),
                       sample(c('TRUE', 'FALSE'), n_features, replace = TRUE, prob = c(0.7, 0.3)),
                       'FALSE')
    meta <- data.frame(
        'Alignment ID' = seq_len(n_features),
        'Average Rt(min)' = round(runif(n_features, 0.5, 15), 3),
        'Average Mz' = round(runif(n_features, 80, 900), 4),
        'Metabolite name' = ifelse(tags == 'Unknown', 'Unknown', paste0(tags, '_', seq_len(n_features))),
        'Adduct type' = sample(c('[M+H]+', '[M+Na]+', '[M+NH4]+', '[M-H]-'), n_features, replace = TRUE),
        'Fill %' = sample(20:100, n_features, replace = TRUE),
        'MS/MS assigned' = has_msms,
        'Annotation tag (VS1.0)' = tags,
        check.names = FALSE
    )
    sample_names <- paste0('QC_', seq_len(2), c('', ''))
    sample_names <- c(paste0('Sample_', seq_len(n_samples)))
    intensities <- as.data.frame(matrix(round(rlnorm(n_features * n_samples, 11, 1.6)),
                                        nrow = n_features, dimnames = list(NULL, sample_names)),
                                 check.names = FALSE)
    body <- cbind(meta, intensities)
    # Real MS-DIAL exports prepend 4 metadata rows (class / file type / injection order / batch)
    # above the column header. Reproduce that offset so the skip=4 parse is exercised.
    header_block <- matrix('', nrow = 4, ncol = ncol(body))
    writeLines(apply(header_block, 1, paste, collapse = '\t'), path)
    suppressWarnings(write.table(body, path, sep = '\t', append = TRUE,
                                 row.names = FALSE, col.names = TRUE, quote = FALSE))
}

export_path <- file.path(tempdir(), 'AlignResult.txt')
make_synthetic_export(export_path)

msdial <- read.csv(export_path, sep = '\t', skip = 4, check.names = FALSE)
cat('Parsed', nrow(msdial), 'features x', ncol(msdial), 'columns\n')

meta_cols <- intersect(c('Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name',
                         'Adduct type', 'Fill %', 'MS/MS assigned', 'Annotation tag (VS1.0)'),
                       colnames(msdial))
last_meta <- max(match(meta_cols, colnames(msdial)))
sample_cols <- colnames(msdial)[(last_meta + 1):ncol(msdial)]

feature_info <- msdial[, meta_cols]
intensity <- as.matrix(msdial[, sample_cols])
rownames(intensity) <- msdial[['Alignment ID']]
cat('Split into', length(meta_cols), 'metadata cols and', length(sample_cols), 'sample cols\n')

# Fill% floor: below this a feature is mostly gap-filled noise-floor integrals, not
# measurements, so it would fabricate intensity for truly below-detection samples.
fill_floor <- 70
keep_fill <- feature_info[['Fill %']] >= fill_floor

# A named hit without MS/MS is at best a putative (MSI Level 3) ID; require MS/MS for identity.
has_msms <- feature_info[['MS/MS assigned']] == 'TRUE'
feature_info$msi_level <- ifelse(feature_info[['Annotation tag (VS1.0)']] %in% c('Metabolite', 'Lipid') & has_msms, 2L,
                          ifelse(feature_info[['Annotation tag (VS1.0)']] == 'Suggested', 3L, NA_integer_))

filtered_intensity <- intensity[keep_fill, ]
filtered_info <- feature_info[keep_fill, ]
cat('After Fill% >=', fill_floor, ':', nrow(filtered_intensity), '/', nrow(intensity), 'features\n')

cat('\nMSI confidence levels (survivors):\n')
print(table(filtered_info$msi_level, useNA = 'ifany'))

out_path <- file.path(tempdir(), 'msdial_filtered.tsv')
write.table(cbind(filtered_info, log2(filtered_intensity + 1)), out_path,
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('\nWrote filtered table to', out_path, '\n')

unlink(c(export_path, out_path))
cat('Cleaned up temp files.\n')
