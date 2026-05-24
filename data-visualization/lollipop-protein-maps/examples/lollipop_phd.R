# Reference: maftools 2.18+, trackViewer 1.38+ | Verify API if version differs

# PhD-level lollipop plot encoding the four correctness traps:
# (1) HGVSp_Short column verified, (2) explicit isoform, (3) counts printed,
# (4) hotspots labeled with canonical labels.

library(maftools)

# 1. INPUT -- MAF with HGVSp_Short column
maf <- read.maf(maf = 'cohort.maf', clinicalData = clinical)
stopifnot('HGVSp_Short' %in% colnames(maf@data))

# 2. VARIANT CLASS PALETTE -- CVD-safe
class_col <- c(Missense_Mutation = '#D55E00',
               Nonsense_Mutation = '#000000',
               Frame_Shift_Del   = '#0072B2',
               Frame_Shift_Ins   = '#56B4E9',
               Splice_Site       = '#CC79A7',
               In_Frame_Del      = '#009E73',
               In_Frame_Ins      = '#F0E442')

# 3. LOLLIPOP for TP53 with canonical hotspot labels
pdf('TP53_lollipop.pdf', width = 8, height = 4)
lollipopPlot(maf = maf,
             gene = 'TP53',
             AACol = 'HGVSp_Short',
             labelPos = c(175, 248, 273),                   # canonical hotspots
             labPosSize = 1.0,
             showMutationRate = TRUE,
             domainLabelSize = 1,
             printCount = TRUE,                             # CRITICAL: show counts
             colors = class_col,
             proteinID = 'P04637')                          # explicit canonical isoform
dev.off()

# 4. COHORT COMPARISON -- subtype A vs subtype B paired lollipop
maf_lumin <- subsetMaf(maf, clinQuery = 'Subtype == "Luminal"')
maf_basal <- subsetMaf(maf, clinQuery = 'Subtype == "Basal"')

pdf('TP53_lollipop_subtype.pdf', width = 10, height = 5)
lollipopPlot2(m1 = maf_lumin, m2 = maf_basal,
              gene = 'TP53',
              m1_name = 'Luminal',
              m2_name = 'Basal',
              AACol1 = 'HGVSp_Short', AACol2 = 'HGVSp_Short',
              colors = class_col)
dev.off()

# 5. trackViewer for custom domain coordinates (use when maftools Pfam cache lags)
library(trackViewer)
library(GenomicRanges)

# Build SNP GRanges -- one per mutation, size by count
mutation_summary <- maf@data[Hugo_Symbol == 'TP53',
                              .(count = .N, class = Variant_Classification[1]),
                              by = .(aa_pos = as.numeric(sub('p\\.[A-Z](\\d+).*', '\\1', HGVSp_Short)))]

snps <- GRanges('chr17',
                IRanges(mutation_summary$aa_pos, width = 1),
                color = class_col[mutation_summary$class],
                score = mutation_summary$count)
names(snps) <- mutation_summary$class

# Build domain feature GRanges from UniProt
# TP53 (P04637): Transactivation domain 1-42, DNA-binding 102-292, Tetramerization 323-356, Regulatory 363-393
features <- GRanges('chr17',
                    IRanges(c(1, 102, 323, 363),
                            width = c(41, 190, 33, 30),
                            names = c('TAD', 'DNA-binding', 'Tetramer', 'Reg')),
                    fill = c('#56B4E9', '#0072B2', '#009E73', '#CC79A7'),
                    height = 0.04)

pdf('TP53_trackviewer.pdf', width = 10, height = 4)
lolliplot(snps, features,
          ylab = 'Mutation count',
          xaxis = TRUE, yaxis = TRUE,
          legend = list(labels = names(class_col), col = class_col))
dev.off()

# 6. INTERACTIVE HTML supplement via g3viz
library(g3viz)
mutation_data <- hgvspChange2protein(maf, gene = 'TP53')
g3Lollipop(mutation_data,
           gene.symbol = 'TP53',
           protein.change.col = 'AA_Change',
           plot.options = g3Lollipop.theme(theme.name = 'nature'),
           output.filename = 'TP53_lollipop.html')

# 7. HOTSPOT VALIDATION reminder
# Novel hotspots require: (a) recurrence >2x background, (b) replication in TCGA Pan-Cancer + ICGC,
# (c) formal test via MutSig (Lawrence 2014) or statisticalhotspot (Chang 2016)
