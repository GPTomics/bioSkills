# Reference: ggseqlogo 0.2+ | Verify API if version differs

# PhD-level sequence logo encoding the four correctness traps:
# (1) bits not probability, (2) genome background composition, (3) explicit alphabet,
# (4) N annotated for small-sample awareness.

library(ggseqlogo)
library(ggplot2)

# 1. INPUT -- vector of aligned same-length DNA sequences
seqs <- readLines('aligned_motif.fa')
seqs <- seqs[!startsWith(seqs, '>')]
stopifnot(all(nchar(seqs) == nchar(seqs[1])))             # equal length required
n <- length(seqs)

# 2. BACKGROUND -- human genome composition; pass explicitly
human_bg <- c(A = 0.295, C = 0.205, G = 0.205, T = 0.295)

# 3. BITS LOGO with explicit background
p_bits <- ggseqlogo(seqs,
                    method = 'bits',                        # canonical encoding
                    bg_freq = human_bg) +
    labs(title = sprintf('Motif logo (bits; N = %d; human bg)', n)) +
    theme(plot.title = element_text(size = 10))

# 4. PROBABILITY logo for COMPARISON
p_prob <- ggseqlogo(seqs, method = 'probability') +
    labs(title = sprintf('Same motif (probability; N = %d)', n)) +
    theme(plot.title = element_text(size = 10))

# Side-by-side
library(patchwork)
ggsave('logo_comparison.pdf', p_bits / p_prob, width = 89, height = 80, units = 'mm',
       device = cairo_pdf)

# 5. MULTI-MOTIF STACK (TF-A vs TF-B)
multi <- list(`CTCF (N=200)` = ctcf_seqs,
              `REST (N=180)` = rest_seqs,
              `GATA1 (N=150)` = gata1_seqs)
p_stack <- ggseqlogo(multi, method = 'bits', bg_freq = human_bg, ncol = 1)
ggsave('multi_logo.pdf', p_stack, width = 89, height = 110, units = 'mm', device = cairo_pdf)

# 6. PROTEIN LOGO with custom functional-class palette
phospho_neighborhoods <- readLines('phospho_aligned.fa')
phospho_neighborhoods <- phospho_neighborhoods[!startsWith(phospho_neighborhoods, '>')]

protein_scheme <- make_col_scheme(
    chars = c('S','T','Y',           # phospho-acceptors
              'K','R','H',           # basic
              'D','E',               # acidic
              'A','V','L','I','M','F','W','C','G','P','N','Q'),  # hydrophobic/other
    cols  = c(rep('#D55E00', 3),
              rep('#0072B2', 3),
              rep('#CC79A7', 2),
              rep('#009E73', 12)))

p_prot <- ggseqlogo(phospho_neighborhoods,
                    method = 'bits',
                    seq_type = 'aa',
                    col_scheme = protein_scheme) +
    labs(title = sprintf('Phospho-substrate neighborhoods (N = %d)', length(phospho_neighborhoods)))
ggsave('protein_logo.pdf', p_prot, width = 130, height = 50, units = 'mm', device = cairo_pdf)

# 7. PWM INPUT instead of sequences
pwm_counts <- matrix(c(85, 5, 5, 5,    # position 1: A-conserved
                       5, 85, 5, 5,    # position 2: C-conserved
                       30, 20, 30, 20,  # position 3: A/G slight bias
                       5, 5, 5, 85),   # position 4: T-conserved
                     ncol = 4, byrow = FALSE,
                     dimnames = list(c('A', 'C', 'G', 'T'), NULL))
# ggseqlogo expects rows = letters, columns = positions (transpose if needed)
ggseqlogo(pwm_counts, method = 'bits', bg_freq = human_bg)
