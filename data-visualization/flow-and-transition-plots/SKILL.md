---
name: bio-data-visualization-flow-and-transition-plots
description: Build Sankey, alluvial, river, and CONSORT-style flow diagrams to visualize cohort transitions, cell-state changes, or pipeline filtering using ggalluvial, networkD3, plotly, and consort. Use when showing how entities move between categories across timepoints (cell states, drug response classes, patient flow through a trial) or filtering pipelines (variants filtered through QC stages).
tool_type: mixed
primary_tool: ggalluvial
---

## Version Compatibility

Reference examples tested with: ggalluvial 0.12+, networkD3 0.4+, plotly 4.10+, consort 0.2+ (CONSORT diagrams), pySankey 0.0.1+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Flow and Transition Plots

**"Show how things flow between categories"** -> Render entities as ribbons whose width encodes count, flowing between ordered columns of categories. Sankey emphasizes total flow magnitude; alluvial emphasizes per-entity continuity (each row's path is traceable); CONSORT formalizes the trial-filtering convention. The decision space: which method (Sankey vs alluvial vs CONSORT), how to order categories within each column, and whether to highlight specific entity trajectories.

- R: `ggalluvial::geom_alluvium`, `networkD3::sankeyNetwork`, `consort::consort_plot`
- Python: `plotly.graph_objects.Sankey`, `pySankey`

## The Single Most Important Modern Insight -- Sankey vs Alluvial Are Different

**Sankey** plots show flow from sources to sinks; each ribbon represents an aggregate count. The horizontal direction is "flow." Use for energy flows, web-traffic funnels, cohort dropouts.

**Alluvial** plots track *individual entities* through multiple ordered category columns (axes). Each row of input data becomes a continuous ribbon; intersections at each axis show counts in each category. Use for cell-state transitions across timepoints, drug-response trajectories, longitudinal class changes.

A Sankey shows "100 cells became neuron, 50 became glia"; an alluvial shows "of the 100 that became neurons at t2, 80 came from the proliferating pool at t1." Different encoding, different scientific story.

## Decision Tree by Use Case

| Use case | Recommended | Tool |
|----------|-------------|------|
| Single timepoint, source-to-sink flow | Sankey | networkD3, plotly |
| Multi-timepoint entity trajectories | Alluvial | ggalluvial |
| Clinical trial patient flow | CONSORT (formal vertical box-and-arrow) | consort R package |
| Variant filtering pipeline | CONSORT-style flow | consort or manual diagrammeR |
| Cell-state transitions (scRNA timepoints) | Alluvial OR Sankey if 2 timepoints | ggalluvial |
| Drug response class changes | Alluvial | ggalluvial |
| Gene-set membership across conditions | UpSet (alternative) | data-visualization/upset-plots |

## ggalluvial -- Modern R Default for Alluvial

**Goal:** Visualize entity (e.g., cell, patient) trajectories across multiple ordered axes with ribbon-width = count.

**Approach:** Reshape to "lodes" (long) format with one row per entity-stratum, or "alluvia" (wide) format with one row per entity; use `geom_alluvium` for ribbons and `geom_stratum` for column boxes.

```r
library(ggalluvial)
library(ggplot2)

# Wide (alluvia) format: one row per entity
df_wide <- data.frame(
    entity_id = 1:1000,
    t1 = sample(c('A', 'B', 'C'), 1000, replace = TRUE),
    t2 = sample(c('A', 'B', 'C'), 1000, replace = TRUE),
    t3 = sample(c('A', 'B', 'C'), 1000, replace = TRUE))

ggplot(df_wide, aes(axis1 = t1, axis2 = t2, axis3 = t3)) +
    geom_alluvium(aes(fill = t1), alpha = 0.7, width = 1/6) +
    geom_stratum(width = 1/6, fill = 'grey90', color = 'black') +
    geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c('t1', 't2', 't3')) +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#009E73')) +
    labs(y = 'Entities', x = NULL) +
    theme_classic()
```

`aes(fill = t1)` colors each ribbon by its starting class — common pattern for "where did this end up cluster come from?" stories.

## networkD3 -- Interactive Sankey

```r
library(networkD3)

# Nodes and links
nodes <- data.frame(name = c('Source A', 'Source B', 'Sink X', 'Sink Y', 'Sink Z'))
links <- data.frame(source = c(0, 0, 1, 1),
                    target = c(2, 3, 3, 4),
                    value  = c(40, 30, 50, 20))

sankeyNetwork(Links = links, Nodes = nodes,
              Source = 'source', Target = 'target', Value = 'value',
              NodeID = 'name',
              colourScale = JS('d3.scaleOrdinal(d3.schemeCategory10);'),
              fontSize = 12, nodeWidth = 30, height = 400, width = 700)
```

networkD3 produces interactive HTML — drag nodes, hover for values. For static publication figure, screenshot or export via `webshot2`.

## plotly Sankey (Python)

```python
import plotly.graph_objects as go

fig = go.Figure(go.Sankey(
    node=dict(label=['Source A', 'Source B', 'Sink X', 'Sink Y', 'Sink Z'],
              color=['#0072B2', '#56B4E9', '#D55E00', '#E69F00', '#009E73']),
    link=dict(source=[0, 0, 1, 1],
              target=[2, 3, 3, 4],
              value=[40, 30, 50, 20],
              color=['rgba(0,114,178,0.4)'] * 4)))
fig.update_layout(title='Flow', font_size=12)
fig.write_html('sankey.html')
fig.write_image('sankey.pdf')           # requires Kaleido (NOT orca; orca is EOL)
```

## CONSORT Diagrams -- The Formal Trial-Flow Standard

CONSORT 2010 (Schulz 2010 *BMJ* 340:c332) is the canonical clinical-trial flow diagram. The `consort` R package implements the structure:

```r
library(consort)

# Trial enrollment flow
g <- add_box(txt = c('Assessed for eligibility (n=200)'))
g <- add_side_box(g, txt = c('Excluded (n=50)\n  - Not meeting criteria (n=30)\n  - Declined (n=15)\n  - Other (n=5)'))
g <- add_box(g, txt = c('Randomized (n=150)'))
g <- add_split(g, txt = c('Allocated to intervention (n=75)\n  - Received as allocated (n=70)\n  - Did not receive (n=5)',
                          'Allocated to control (n=75)\n  - Received as allocated (n=73)\n  - Did not receive (n=2)'))
g <- add_box(g, txt = c('Lost to follow-up (n=2)\nDiscontinued (n=3)',
                        'Lost to follow-up (n=1)\nDiscontinued (n=2)'))
g <- add_box(g, txt = c('Analysed (n=75)\nExcluded from analysis (n=0)',
                        'Analysed (n=75)\nExcluded from analysis (n=0)'))
plot(g)
```

CONSORT is a *required* element in randomized trial publication (CONSORT 2010 statement, item 13a).

## Per-Method Failure Modes

### Sankey used when alluvial is appropriate

**Trigger:** Multi-timepoint cell-state data plotted as Sankey instead of alluvial.

**Mechanism:** Sankey collapses to source-sink summary; loses entity-trajectory continuity.

**Symptom:** Reader sees "cluster A -> 50% to B, 50% to C" but cannot trace individual trajectories.

**Fix:** Use ggalluvial for multi-axis trajectories; Sankey for single-step source-to-sink.

### Category ordering within column not specified

**Trigger:** Default ggalluvial ordering by frequency.

**Mechanism:** Categories shuffle position across columns; ribbons cross excessively.

**Symptom:** Visual spaghetti; hard to follow.

**Fix:** Set explicit factor levels (`factor(t1, levels = c('A', 'B', 'C'))`) AND consider ggalluvial's `lode.guidance` to minimize crossings.

### Ribbon coloring by destination instead of origin

**Trigger:** `geom_alluvium(aes(fill = t3))` for a "where did these come from" story.

**Mechanism:** Color encodes the wrong axis; readers misinterpret.

**Symptom:** Story is "where did final cluster Z come from" but ribbons are colored by Z — every ribbon to Z is the same color.

**Fix:** `aes(fill = t1)` if origin matters; `aes(fill = t3)` if destination matters.

### CONSORT diagram missing required boxes

**Trigger:** Skipping "Lost to follow-up" or "Excluded from analysis" boxes.

**Mechanism:** CONSORT 2010 requires reporting at each stage.

**Symptom:** Submission flagged for non-compliance with CONSORT 2010 item 13a.

**Fix:** Use `consort` package which scaffolds the required structure; cross-check against CONSORT 2010 statement.

### plotly Sankey value sum mismatch

**Trigger:** Source-to-target sums don't balance.

**Mechanism:** plotly Sankey requires conservation: sum of in-flows = sum of out-flows at each non-terminal node.

**Symptom:** Layout renders but node sizes look wrong; ribbons stretch/compress incorrectly.

**Fix:** Verify upstream data: per-node `sum(value where target=node) == sum(value where source=node)` for internal nodes.

### Static export of plotly Sankey fails silently

**Trigger:** `fig.write_image('sankey.pdf')` without kaleido installed.

**Mechanism:** plotly defaults to Kaleido for static export since orca EOL; kaleido is optional dependency.

**Symptom:** No file written; no error in some plotly versions.

**Fix:** `pip install kaleido`; verify with `import kaleido`.

## Reconciliation: When Implementations Differ

| Pattern | Cause | Action |
|---------|-------|--------|
| ggalluvial and networkD3 show different orderings | Different default stratum/node ordering | Set explicit factor levels / node order |
| CONSORT box counts don't sum | Box-content arithmetic error | Audit each box; consort package enforces structure |
| Cells appear/disappear between alluvial axes | Missing data in some timepoints | Decide: drop entities with NA; OR add "Missing" category |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Max categories per column for legibility | 5-7 | Visualization practical |
| Max axes for alluvial | 4-5 | Above this ribbons too crossed |
| CONSORT requirement | Required for RCTs | Schulz 2010 CONSORT 2010 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Excessive ribbon crossing | Categories unordered | Explicit factor levels; lode.guidance |
| Trajectories not traceable | Sankey used instead of alluvial | Switch to ggalluvial |
| CONSORT non-compliant | Missing required boxes | Use consort package |
| Sankey node sizes wrong | Flow not conserved | Audit source-target sums |
| plotly static export blank | kaleido not installed | `pip install kaleido` |
| Color story unclear | Wrong axis for fill | Decide origin vs destination story |

## References

- Brunson J. 2020. ggalluvial: Layered grammar for alluvial plots. *J Open Source Softw* 5(49):2017.
- Sankey MH. 1898. The thermal efficiency of steam engines. *Proc Inst Civil Eng* 134:278-312. (origin)
- Schulz KF, Altman DG, Moher D; CONSORT Group. 2010. CONSORT 2010 Statement: updated guidelines for reporting parallel group randomised trials. *BMJ* 340:c332.
- Riehmann P, Hanfler M, Froehlich B. 2005. Interactive Sankey diagrams. *IEEE Symp Information Visualization*.

## Related Skills

- data-visualization/upset-plots - Alternative for set-intersection rather than flow
- clinical-biostatistics/trial-reporting - CONSORT diagrams in trial publication
- single-cell/trajectory-inference - Cell-state transition data for alluvial
- workflows/biomarker-pipeline - Pipeline filtering flows
