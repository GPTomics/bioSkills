# bioSkills Skill Writing Operational Reference

Complete cookbook for writing skills in this repo. Generated from CLAUDE.md, Anthropic best-practices doc, and recent decision-grade exemplars.

## A. FRONTMATTER

Required:
- `name`: <=64 chars, lowercase + hyphens, no reserved words ("anthropic", "claude")
- `description`: <=1024 chars, third-person, must include "Use when..." trigger
- `tool_type`: one of python | r | cli | mixed
- `primary_tool`: single value only (no commas)

Workflow skills also need: `workflow: true`, `depends_on: [...]`, `qc_checkpoints: {...}`

## B. SKILL.md SECTION ORDER

1. Frontmatter
2. `## Version Compatibility` (mandatory if code present)
3. Main narrative / algorithm taxonomy
4. Decision tree / per-tool failure modes
5. Common Errors table (for code-heavy skills)
6. References block (in-line citations + bibliography)
7. `## Related Skills` (LAST, plain text, qualified paths)

Length cap: 500 lines. Voice: third-person throughout (including code comments).

## C. DECISION-GRADE PATTERN

Recent postdoc-grade skills include:

**Algorithmic Taxonomy table:**
| Tool | Model | Input | Output | Strength | Fails when |
|------|-------|-------|--------|----------|------------|

**Decision Tree by Scenario:**
| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|

**Per-Tool Failure Modes** (one per failure):
```
### [Tool] -- [Pitfall]
**Trigger:** condition activating failure
**Mechanism:** why the tool produces wrong output
**Symptom:** observable signal
**Fix:** how to detect and correct
```

**Quantitative Thresholds:**
| Quantity | Threshold | Source / Rationale |
|----------|-----------|--------------------|

**Reconciliation paragraph** when methods disagree.

**Goal/Approach intent-first code** for multi-step pipelines:
```markdown
**Goal:** What this section achieves.
**Approach:** How the code achieves it.
```python
...code...
```
```

## D. VERSION COMPATIBILITY (verbatim template)

```markdown
## Version Compatibility

Reference examples tested with: <package> <version>+, <package> <version>+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.
```

NO version annotations in heading titles. Version Compatibility block is single source of truth.

## E. usage-guide.md SECTIONS

1. `## Overview` (2-3 sentences)
2. `## Prerequisites`
3. `## Quick Start` (bullets only: `- "..."`)
4. `## Example Prompts` (blockquotes only: `> "..."`)
5. `## What the Agent Will Do`
6. `## Tips`
7. `## Related Skills` (plain text, no bold)

## F. examples/ RULES

- At least one runnable script per skill
- Version header as first line:
  `# Reference: <pkg> <ver>+, <pkg> <ver>+ | Verify API if version differs`
- Magic numbers documented with rationale comments

## G. Category README

```markdown
# category-name
## Overview
brief description.
**Tool type:** ... | **Primary tools:** ...
## Skills
| Skill | Description |
## Example Prompts
- "..."
## Requirements
```bash
pip install ...
```
## Related Skills
- **other-category** - description    ← BOLD here (only place)
```

Bold ONLY in category READMEs. SKILL.md and usage-guide.md use plain text.

## H. 19-CHECK VALIDATION SUITE

```bash
# 1. Names <=64
grep -h "^name:" */*/SKILL.md | awk -F': ' '{print $2}' | awk 'length > 64'
# 2. Names lowercase
grep -h "^name:" */*/SKILL.md | awk -F': ' '{print $2}' | grep '[A-Z]'
# 3. "Use when" present
grep -L "Use when" */*/SKILL.md
# 4. Third-person description
grep -l -E "^description:.*\b(you|your|You|Your)\b" */*/SKILL.md
# 5. Under 500 lines
find . -name "SKILL.md" -exec wc -l {} \; | awk '$1 > 500'
# 6. Related Skills present
grep -L "## Related Skills" */*/SKILL.md
# 7. primary_tool single value
grep -r "^primary_tool:.*,.*" --include="SKILL.md" .
# 8. Related Skills no bold
grep -A20 "## Related Skills" */*/SKILL.md | grep "\*\*"
# 9-12. usage-guide sections
grep -L "## Overview" */*/usage-guide.md
grep -L "## Prerequisites" */*/usage-guide.md
grep -L "## Quick Start" */*/usage-guide.md
grep -L "## Example Prompts" */*/usage-guide.md
# 13. Quick Start uses bullets not blockquotes
grep -A10 "## Quick Start" */*/usage-guide.md | grep "^>"
# 14. No empty examples/
find . -type d -name "examples" -exec sh -c 'test -z "$(ls -A "$1")" && echo "$1"' _ {} \;
# 15. No broken skill names
grep -r "vcf-filtering\|clustering-annotation\|enrichment-analysis\|multi-panel-figures\|contact-matrices\|quality-assessment\|dexseq-analysis\|variant-to-annotation" --include="*.md" . | grep -v CLAUDE.md
# 16. Cross-category refs qualified
grep -A30 "## Related Skills" */*/SKILL.md | grep "^[^:]*:- " | grep -v "/" | grep -v "^--$"
# 17. usage-guide Related Skills no bold
grep -A10 "## Related Skills" */*/usage-guide.md | grep "\*\*"
# 18. Version Compatibility present
for f in $(grep -rl '```' */*/SKILL.md); do grep -L "## Version Compatibility" "$f"; done
# 19. No version annotations in headings
grep -n '^#' */*/SKILL.md | grep -v 'Version Compatibility' | grep -E '\([^)]*[0-9]+\.[0-9]+\+' | head -20
```

## I. KNOWN HALLUCINATION CATEGORIES TO WATCH FOR

- Wrong column names in R packages (`$FDR` vs `$adj.P.Val`)
- NULL vs FALSE in R (many functions check `!is.null()`)
- Deprecated parameter names (verify with help() / ?fn)
- Sklearn metric strings (precompute distances, use `metric='precomputed'`)
- CLI flag fabrication (always check `--help`)
- Wrong default methods (e.g. `multipletests` defaults to Holm-Sidak not Holm)
- Table orientation in crosstab (scipy reorders columns)

## J. CHEMOINFORMATICS-SPECIFIC LESSONS

- RDKit moved from MolStandardize (deprecated Q1 2024) to rdMolStandardize
- Open Babel 3.x: `from openbabel import pybel` (not `import pybel`)
- ETKDGv3 is default conformer embedding (replaces ETKDGv2)
- AutoDock Vina 1.2 Python API: `from vina import Vina`
- ChEMBL pipeline uses RDKit + MolVS; canSARchem uses ChemAxon for canonical tautomer
- ADMETlab 3.0 API: https://admetlab3.scbdd.com/api/predict (POST with `{"smiles": [...]}`)
- GNINA = SMINA fork = Vina fork (Vina < SMINA < GNINA); GNINA adds CNN scoring
- DiffDock-L recommended hybrid with GNINA rescoring (not standalone)
- PoseBusters 7.5% vdW overlap cutoff
- Boltz-2 affinity module approaches FEP accuracy 1000x faster
- REINVENT 4: 4 generators (de novo, scaffold decoration, linker, optimization) x 3 algorithms (TL, RL, CL)
- AiZynthFinder 4.0: MO-MCTS + Chemformer template-free
- OpenFE has ABFE protocol (relative + absolute binding free energies)
- ECFP4 has radius=2 (diameter 4); ECFP6 has radius=3 (diameter 6)
