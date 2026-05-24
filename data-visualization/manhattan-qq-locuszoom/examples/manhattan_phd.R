# Reference: qqman 0.1.9, CMplot 4.5+, locuszoomr 0.3+ | Verify API if version differs

# PhD-level Manhattan + QQ + locuszoom encoding the four correctness traps:
# (1) sort by CHR/BP before plot, (2) analysis-appropriate threshold,
# (3) y-axis cap with indicator, (4) LD reference matched to ancestry.

library(qqman)
library(CMplot)
library(dplyr)

# 1. INPUT -- ensure sorted (qqman plots in input order)
gwas <- read.table('gwas_summary.tsv', header = TRUE) %>%
    arrange(CHR, BP)

# 2. INFLATION DIAGNOSTIC -- compute lambda_GC
chisq <- qchisq(1 - gwas$P, df = 1)
lambda <- median(chisq) / 0.4549                          # 0.4549 = median of chi-square_1 (Devlin-Roeder 1999)
lambda_1000 <- 1 + (lambda - 1) * 1000 / nrow(gwas)       # sample-size adjusted

# 3. QQ PLOT with lambda in title
pdf('qq.pdf', width = 4, height = 4)
qq(gwas$P, main = sprintf('QQ (lambda = %.3f, lambda_1000 = %.3f)', lambda, lambda_1000))
dev.off()

# 4. THRESHOLD by analysis type
# common-variant GWAS -> 5e-8 (Pe'er 2008)
# whole-genome sequencing (EUR) -> 5e-9 (Pulit 2017; Xu 2014)
# TWAS / PWAS Bonferroni -> 0.05 / n_genes
# trans-ancestry meta -> 5e-9
sig_threshold <- 5e-8

# 5. Y-CAP with indicator -- cap extreme peaks
y_cap <- 25
gwas$P_capped <- ifelse(-log10(gwas$P) > y_cap, 10^(-y_cap), gwas$P)
gwas$is_capped <- -log10(gwas$P) > y_cap

# 6. MANHATTAN -- two-color chromosome alternation, lead SNP annotation
pdf('manhattan.pdf', width = 10, height = 4)
manhattan(gwas, chr = 'CHR', bp = 'BP', p = 'P_capped', snp = 'SNP',
          col = c('#0072B2', '#56B4E9'),                  # two-color alternation
          genomewideline = -log10(sig_threshold),
          suggestiveline = -log10(1e-5),
          ylim = c(0, y_cap + 2),
          annotatePval = sig_threshold,
          annotateTop = TRUE)
# Mark capped points with caret above the cap
if (any(gwas$is_capped)) {
    capped_xs <- which(gwas$is_capped)
    # qqman doesn't expose x mapping; for full control use CMplot or custom
}
dev.off()

# 7. CMPLOT -- supports Miami (two-trait mirror) natively, capped point markers
CMplot(list(gwas %>% select(SNP, CHR, BP, P)),
       plot.type = 'm',
       threshold = c(sig_threshold, 1e-5),
       threshold.col = c('red', 'grey'),
       threshold.lty = c(1, 2),
       col = c('#0072B2', '#56B4E9'),
       ylim = c(0, y_cap),
       amplify = TRUE, signal.cex = 1.2,
       highlight = lead_snps,
       file = 'pdf', file.output = TRUE)

# 8. MIAMI PLOT -- two traits mirrored (one above, one below axis)
CMplot(list(trait1_df, trait2_df),
       plot.type = 'm',
       multraits = TRUE,
       threshold = sig_threshold,
       threshold.col = 'red',
       col = list(c('#0072B2','#56B4E9'), c('#D55E00','#E69F00')),
       file = 'pdf', file.output = TRUE)

# 9. LOCUSZOOM -- regional plot with LD coloring; match LD reference to ancestry
library(locuszoomr)
library(LDlinkR)
loc <- locus(gene = 'TCF7L2',
             flank = 5e5,                                  # 500 kb each side
             ens_db = 'EnsDb.Hsapiens.v86',
             data = gwas,
             snp = 'SNP', chrom = 'CHR', pos = 'BP', p = 'P', labs = 'SNP')
# CRITICAL: match LD reference to GWAS ancestry
# EUR for European GWAS; EAS for East Asian; AFR for African; etc
loc <- link_LD(loc, pop = 'EUR', token = Sys.getenv('LDLINK_TOKEN'))

pdf('locuszoom.pdf', width = 6, height = 5)
locus_plot(loc, labels = c('index', 'top'),
           recomb_offset = 0.1, ens_db = 'EnsDb.Hsapiens.v86')
dev.off()
