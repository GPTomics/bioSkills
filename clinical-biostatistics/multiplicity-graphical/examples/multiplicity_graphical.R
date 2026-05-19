# Multiplicity control via graphical procedures and step-down methods.
#
# Reference: gMCP 0.8.16+ | Verify API if version differs
# Reference: graphicalMCP 0.2+ | Verify API if version differs
#
# Covers Bretz-Maurer 2009 graphical procedures, Holm/Hochberg/Hommel,
# serial and parallel gatekeeping, and the FDA 2022 Multiple Endpoints
# Guidance categories (co-primary, multiple primary, primary+key secondary).

# ----------------------------------------------------------------------
# 1. Standard Bretz-Maurer graphical procedure: primary + 2 key secondary
# ----------------------------------------------------------------------
library(gMCP)

# Hypothesis labels
hyps <- c('Primary', 'KeySec1', 'KeySec2')

# Initial alpha weights: all on Primary
weights <- c(1, 0, 0)

# Transition matrix: rows = source, cols = target
# When Primary rejects, alpha propagates 0.5 to KeySec1 and 0.5 to KeySec2
# When KeySec1 or KeySec2 rejects, all weight goes to the remaining hypothesis
transitions <- matrix(c(
    0,   0.5, 0.5,
    0,   0,   1,
    0,   1,   0
), nrow = 3, byrow = TRUE,
   dimnames = list(hyps, hyps))

graph_main <- matrix2graph(transitions, weights)
print(graph_main)

# Apply to trial p-values at alpha = 0.025 (one-sided regulatory)
p_vals <- c(Primary = 0.018, KeySec1 = 0.042, KeySec2 = 0.038)
result_main <- gMCP(graph_main, pvalues = p_vals, alpha = 0.025)
print(result_main)
# Primary rejects (0.018 < 0.025) -> alpha propagates 0.0125 to each secondary
# KeySec2 = 0.038 < 0.0125 * 2 = 0.025 (after propagation): rejects -> alpha to KeySec1
# KeySec1 = 0.042 > 0.025: does not reject


# ----------------------------------------------------------------------
# 2. Serial gatekeeping (pure hierarchical sequence)
# ----------------------------------------------------------------------
# H1 -> H2 -> H3 -> H4; full weight propagates only if predecessor rejects

hyps_serial <- c('H1', 'H2', 'H3', 'H4')
weights_serial <- c(1, 0, 0, 0)
trans_serial <- matrix(c(
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1,
    0, 0, 0, 0
), nrow = 4, byrow = TRUE,
   dimnames = list(hyps_serial, hyps_serial))

graph_serial <- matrix2graph(trans_serial, weights_serial)
p_serial <- c(H1 = 0.01, H2 = 0.02, H3 = 0.06, H4 = 0.001)
result_serial <- gMCP(graph_serial, pvalues = p_serial, alpha = 0.025)
print(result_serial)
# H1 rejects -> H2 tests at 0.025 -> H2 rejects -> H3 tests at 0.025 -> H3 does NOT reject -> H4 not tested


# ----------------------------------------------------------------------
# 3. Holm vs Hochberg vs Hommel
# ----------------------------------------------------------------------
# These are step-down/step-up procedures applied to a vector of p-values

p_vals_unstructured <- c(p1 = 0.018, p2 = 0.042, p3 = 0.038, p4 = 0.015, p5 = 0.029)
alpha <- 0.05

# Holm: any dependence; uniformly dominates Bonferroni
p_holm <- p.adjust(p_vals_unstructured, method = 'holm')

# Hochberg: requires PRDS; dominates Holm under PRDS
p_hochberg <- p.adjust(p_vals_unstructured, method = 'hochberg')

# Hommel: requires PRDS; dominates Hochberg by 1-3%
p_hommel <- p.adjust(p_vals_unstructured, method = 'hommel')

# Bonferroni (most conservative)
p_bonf <- p.adjust(p_vals_unstructured, method = 'bonferroni')

comparison <- data.frame(
    raw = p_vals_unstructured,
    Bonferroni = p_bonf,
    Holm = p_holm,
    Hochberg = p_hochberg,
    Hommel = p_hommel
)
print(round(comparison, 4))


# ----------------------------------------------------------------------
# 4. PRDS check (informal)
# ----------------------------------------------------------------------
# Hochberg/Hommel require PRDS (positive regression dependence on subsets).
# Under positive correlation among test statistics, PRDS typically holds.
# Under negative correlation (e.g., LDL vs HDL, complementary efficacy/safety),
# PRDS fails and Hochberg is anti-conservative.
#
# Practical check: estimate correlation among endpoint test statistics from data
# If all pairwise correlations >= 0, PRDS likely holds (use Hommel)
# If any correlation is meaningfully negative, fall back to Holm

# Example: assume p-values come from endpoints with positive correlation 0.4
# In confirmatory work, the correlation matrix should be estimated from the
# primary analysis residuals; sensitivity to PRDS should be documented in SAP


# ----------------------------------------------------------------------
# 5. Graphical procedure with Simes intersection tests (Bretz 2011)
# ----------------------------------------------------------------------
# When endpoints positively correlated, Simes-based intersection tests
# gain power over Bonferroni-based intersection tests

# Demonstration: same graph as Section 1, with Simes local tests
# In gMCP: specify the test type at the gMCP() call

# library(gMCP)
# result_simes <- gMCP(graph_main, pvalues = p_vals, alpha = 0.025, test = 'Simes')
# Compared to result_main above, more rejections expected when correlation > 0


# ----------------------------------------------------------------------
# 6. Subgroup alpha allocation (Dane 2019)
# ----------------------------------------------------------------------
# Allocate alpha across primary, key secondary, discovery subgroup

hyps_sg <- c('Primary', 'KeySec', 'SubgroupHR_OS')
weights_sg <- c(0.7, 0.2, 0.1)  # primary 70%, sec 20%, subgroup 10%
trans_sg <- matrix(c(
    0,    0.5,  0.5,   # primary rejects -> half to sec, half to subgroup
    0.5,  0,    0.5,
    0.5,  0.5,  0
), nrow = 3, byrow = TRUE,
   dimnames = list(hyps_sg, hyps_sg))

graph_sg <- matrix2graph(trans_sg, weights_sg)
print(graph_sg)


# ----------------------------------------------------------------------
# 7. Co-primary endpoints (NO alpha split needed)
# ----------------------------------------------------------------------
# For co-primary (both must reach significance), each test at full alpha
# The multiplicity issue is reversed: power loss because joint power < marginal
# Inflate n so that joint power (e.g., 0.81 = 0.9 * 0.9) meets target

# Example: power calculation
p_marginal_target <- 0.90  # each endpoint
joint_power <- p_marginal_target^2  # if independent
cat(sprintf('\nCo-primary marginal power 0.90 -> joint power %.3f under independence\n',
            joint_power))
# Under positive correlation, joint power > product; sponsor should estimate empirically


# ----------------------------------------------------------------------
# 8. Win-ratio composite (Pocock 2012) -- alternative to multiplicity
# ----------------------------------------------------------------------
# Win-ratio compares each pair of patients (one per arm) on a hierarchical
# composite of outcomes. Single test; no multiplicity.
# Pre-specify component hierarchy in SAP (e.g., CV death > HF hospitalisation > NYHA worsening)
# Implementation: library(WINS) or library(WWR) in R

# library(WINS)
# wins_result <- win.stat(data, ep_type = c('time-to-event', 'binary'), ...)
