library(DiffBind)
library(rtracklayer)

samples <- data.frame(
    SampleID = c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'),
    Condition = c('control', 'control', 'treatment', 'treatment'),
    Replicate = c(1, 2, 1, 2),
    bamReads = c('ctrl1.bam', 'ctrl2.bam', 'treat1.bam', 'treat2.bam'),
    Peaks = c('ctrl1_peaks.narrowPeak', 'ctrl2_peaks.narrowPeak',
              'treat1_peaks.narrowPeak', 'treat2_peaks.narrowPeak'),
    PeakCaller = c('macs', 'macs', 'macs', 'macs')
)
write.csv(samples, 'samples.csv', row.names = FALSE)

dba_obj <- dba(sampleSheet = 'samples.csv')
dba_obj

dba_obj <- dba.count(dba_obj, summits = 250, minOverlap = 2)

dba_obj <- dba.normalize(dba_obj)

dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION)

dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

dba.show(dba_obj, bContrasts = TRUE)

dba.plotPCA(dba_obj, DBA_CONDITION, label = DBA_ID)
dba.plotMA(dba_obj)
dba.plotVolcano(dba_obj)
dba.plotHeatmap(dba_obj, contrast = 1, correlations = FALSE)

db_results <- dba.report(dba_obj, th = 0.05, fold = 1)
db_results

results_df <- as.data.frame(dba.report(dba_obj, th = 1))
write.csv(results_df, 'differential_binding_results.csv', row.names = FALSE)

db_sig <- dba.report(dba_obj, th = 0.05, fold = 2)
export(db_sig, 'significant_diff_peaks.bed', format = 'BED')
