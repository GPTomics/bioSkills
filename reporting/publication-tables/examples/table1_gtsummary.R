# Reference: gtsummary 2.0+, flextable 0.9+ | Verify API if version differs
# Descriptive Table 1 (SMD, not p-values) + an inferential regression table, exported to Word.

library(gtsummary)
library(flextable)

# Descriptive Table 1 by treatment arm. For a randomized trial, report SMD for balance
# instead of a p-value column: a baseline p-value tests a null already known to be true.
# missing = 'ifany' shows a row whenever a value is absent (never hide missingness).
tbl1 <- trial |>
    tbl_summary(by = trt, missing = 'ifany',
                statistic = list(all_continuous() ~ '{median} ({p25}, {p75})',
                                 all_categorical() ~ '{n} ({p}%)')) |>
    add_overall() |>
    add_difference(test = everything() ~ 'smd') |>
    bold_labels()

tbl1 |> as_flex_table() |> save_as_docx(path = 'table1.docx')

# Inferential table: logistic regression as odds ratios with 95% CIs. The effect estimate
# and CI are primary; exponentiate = TRUE reports OR rather than log-odds.
model <- glm(response ~ trt + age + grade, data = trial, family = binomial)
tbl_reg <- tbl_regression(model, exponentiate = TRUE) |> bold_p()

tbl_reg |> as_flex_table() |> save_as_docx(path = 'regression_table.docx')

# Side-by-side univariable + multivariable (a common reporting layout)
tbl_uv <- tbl_uvregression(trial, y = response, method = glm,
                           method.args = list(family = binomial), exponentiate = TRUE)
tbl_merge(list(tbl_uv, tbl_reg), tab_spanner = c('**Univariable**', '**Multivariable**')) |>
    as_flex_table() |> save_as_docx(path = 'uv_mv_table.docx')

cat('Tables written: table1.docx, regression_table.docx, uv_mv_table.docx\n')
