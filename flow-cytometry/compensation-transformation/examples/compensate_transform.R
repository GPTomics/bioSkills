library(flowCore)

# Public flow cytometry data sources:
# - FlowRepository: https://flowrepository.org (FR-FCM-ZZPH for CyTOF)
# - ImmPort: https://www.immport.org (immunology datasets)
# - Bioconductor: flowWorkspaceData, flowStats have example FCS files
# - Example: fcs <- read.FCS(system.file('extdata', 'CytoTrol_CytoTrol_1.fcs', package='flowWorkspaceData'))

# Load data
fcs <- read.FCS('sample.fcs', transformation = FALSE)
cat('Loaded', nrow(fcs), 'events\n')

# Get compensation matrix from keywords
comp_kw <- keyword(fcs)$`$SPILLOVER`
if (!is.null(comp_kw)) {
    comp <- compensation(comp_kw)
    fcs_comp <- compensate(fcs, comp)
    cat('Applied compensation\n')
} else {
    cat('No compensation matrix in file\n')
    fcs_comp <- fcs
}

# Get marker channels (exclude scatter, time)
all_channels <- colnames(fcs)
marker_channels <- all_channels[!grepl('FSC|SSC|Time', all_channels)]
cat('Marker channels:', length(marker_channels), '\n')

# Apply logicle transformation
lgcl <- estimateLogicle(fcs_comp, marker_channels)
fcs_trans <- transform(fcs_comp, lgcl)

cat('Applied logicle transformation\n')

# Save
write.FCS(fcs_trans, 'sample_preprocessed.fcs')
cat('Saved preprocessed file\n')
