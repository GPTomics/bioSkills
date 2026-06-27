---
name: bio-small-rna-seq-target-prediction
description: Predicts and prioritizes miRNA target genes with seed-based tools (miRanda, TargetScan, miRDB) and experimentally validated databases (miRTarBase, multiMiR). Use when deciding that a predicted target is a hypothesis not a finding; ranking by the right score (weighted context++, mirSVR, miRDB); raising confidence by intersecting predictions with inversely-correlated mRNA DE; weighing validated (CLIP/reporter) over predicted evidence; or avoiding the circular enrichment of unfiltered target lists.
tool_type: mixed
primary_tool: miRanda
---

## Version Compatibility

Reference examples tested with: miRanda 3.3a+, BioPython 1.83+, pandas 2.2+, gseapy 1.1+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# miRNA Target Prediction

**"Predict target genes for my miRNAs"** -> Generate candidate mRNA targets by seed complementarity and thermodynamics, then raise confidence with conservation, validated databases, and matched expression.
- CLI: `miranda miRNA.fa UTR.fa -sc 140 -en -20 -strict` for de novo prediction
- Python: query TargetScan/miRDB downloads and miRTarBase for validated interactions

## The governing principle: a predicted target is a hypothesis, not a finding

Seed-based prediction has a false-positive rate near 50% even for conserved sites (Pinzon 2017), because a 6mer seed match occurs by chance roughly once per 4 kb of sequence, so a multi-kb 3' UTR carries many spurious matches. It also has a RECALL problem in the opposite direction: AGO-CLIP and CLASH show that roughly 60% of real interactions are noncanonical, and the bulged, seedless, and 3'-compensatory sites among them are missed entirely by seed-only tools (3'-supplementary sites keep a canonical seed and are still found), so a clean seed list is incomplete, not just imprecise (Helwak 2013). Worse, a single miRNA represses most of its real targets only modestly - typically less than two-fold at the protein level (Baek 2008; Selbach 2008) - so miRNAs are rheostats, not switches, and large single-target claims should be distrusted. The decisive move is therefore not running more predictors (five seed-based tools agreeing is pseudo-replication, not independent evidence) but climbing an evidence ladder and, above all, intersecting predictions with INVERSELY-correlated differentially expressed mRNA or protein from the SAME samples. Prediction proposes; expression disposes.

| Evidence tier (low to high) | What it means |
|------------------------------|----------------|
| seed match alone | weakest; common by chance |
| + conservation (TargetScan PCT) | precision up, recall down (misses species-specific targets) |
| + multiple independent tools / ML (miRDB) | modestly higher precision |
| + AGO-CLIP footprint | the miRNA's complex bound there (but CLIP is cell-type/state-specific - a peak from another tissue is weak evidence) |
| + CLASH/CLEAR-CLIP chimera | direct miRNA-target duplex |
| + anti-correlated matched miRNA/mRNA(protein) DE | functional in YOUR system |
| + reporter / seed-mutation rescue (miRTarBase "strong") | causal, gold standard |

## Decision: which predictor, and how it scores

| Tool | Scoring philosophy | Use for | Caveat |
|------|--------------------|---------|--------|
| TargetScan (context++) | conservation + 14-feature regression of repression | conserved-site prioritization | v7 = context++; v8 = a different Kd/biochemical model |
| miRanda + mirSVR | thermodynamic alignment + expression-trained regression | non-conserved / non-canonical sites | permissive; tune thresholds |
| miRDB / MirTarget | ML (SVM) on CLIP + overexpression data | data-driven ranking (score 0-100) | score >= 80 is the conventional high-confidence cut |
| RNAhybrid | pure MFE hybridization, no seed constraint | exploratory, no-seed sites | most false positives without filters |
| miRTarBase / TarBase (ENCORI) | experimentally validated interactions | the gold tier; anchor claims here | "less strong" CLIP/NGS entries are not individually validated |
| multiMiR | unifies predicted + validated sources | one-call aggregation | inherits each source DB's errors |

## De novo prediction with miRanda

**Goal:** Predict miRNA-mRNA target sites by complementarity and duplex energy.

**Approach:** Align miRNA sequences against 3' UTRs with a minimum score and a maximum (negative) energy, requiring strict seed pairing.

```bash
miranda miRNA.fa UTRs.fa -sc 140 -en -20 -strict -out predictions.txt

# -sc 140: minimum alignment score (keep alignments with score >= 140; default 140)
# -en -20: maximum free energy in kcal/mol (keep energies <= -20; value is negative)
# -strict: require canonical seed pairing at positions 2-8 (no gaps/wobble in seed)
# Tunable: many use -sc 150 / -en -7 (looser) up to -sc 155 / -en -20 (stringent)
```

## Parse miRanda output

**Goal:** Extract interaction records into a DataFrame.

**Approach:** Read the lines miRanda prefixes with '>' (per-hit summary) and pull miRNA, target, score, and energy.

```python
import pandas as pd

def parse_miranda(output_file):
    rows = []
    with open(output_file) as f:
        for line in f:
            if line.startswith('>') and not line.startswith('>>'):
                p = line.strip().split('\t')
                if len(p) >= 5:
                    rows.append({'mirna': p[0].lstrip('>'), 'target': p[1],
                                 'score': float(p[2]), 'energy': float(p[3])})
    return pd.DataFrame(rows)
```

## TargetScan context-scores lookup

**Goal:** Retrieve conserved-site predictions and rank a miRNA's targets.

**Approach:** Read the downloadable per-site context-scores file and rank by the weighted context++ score (more negative = stronger predicted repression across the whole UTR).

```python
import pandas as pd

def query_targetscan(mirbase_id, ts_file='Predicted_Targets_Context_Scores.default_predictions.txt'):
    # Verified column names: the miRNA column is 'Mirbase ID' (NOT 'miRNA family'),
    # the gene column is 'Gene ID', and 'weighted context++ score' aggregates a UTR's sites.
    df = pd.read_csv(ts_file, sep='\t')
    hits = df[df['Mirbase ID'] == mirbase_id]
    return hits.sort_values('weighted context++ score')   # ascending: most negative first
```

## miRDB (machine-learning) lookup

**Goal:** Retrieve ML-based target predictions above the conventional confidence cut.

**Approach:** Read the miRDB prediction download and keep targets with score >= 80.

```python
def query_mirdb(mirna_id, mirdb_file='miRDB_v6.0_prediction_result.txt'):
    df = pd.read_csv(mirdb_file, sep='\t', header=None, names=['mirna', 'refseq', 'score'])
    hits = df[df['mirna'] == mirna_id]
    return hits[hits['score'] >= 80].sort_values('score', ascending=False)
```

## Validated targets and unified lookup

**Goal:** Anchor target claims in experimental evidence rather than prediction.

**Approach:** Query miRTarBase for validated interactions and weight by evidence type; use multiMiR (R) to unify predicted and validated sources in one call.

```python
def get_validated_targets(mirna, mirtarbase_file='miRTarBase_MTI.xlsx'):
    df = pd.read_excel(mirtarbase_file)
    hits = df[df['miRNA'] == mirna]
    # 'Support Type' separates strong (reporter/western/qPCR) from less-strong (CLIP/NGS)
    return hits[['Target Gene', 'Experiments', 'Support Type']]
```

## The confidence move: intersect with anti-correlated mRNA DE

**Goal:** Keep only targets that behave functionally in the actual experiment.

**Approach:** Intersect predicted (or CLIP-supported) targets of UP miRNAs with DOWN mRNAs from matched samples; note the blind spot that translation-only targets may not move at the mRNA level.

```python
def functional_targets(predicted_targets, mrna_de, mirna_direction):
    # mrna_de: DataFrame with index = gene, column 'log2FC' from matched mRNA-seq.
    # Anti-correlation: an UP miRNA should repress -> targets DOWN (and vice versa).
    # Blind spot: miRNAs also act translationally, so some real targets stay flat at
    # the mRNA level (need ribosome profiling / proteomics to see those).
    want_down = mirna_direction == 'up'
    moved = mrna_de[(mrna_de['log2FC'] < 0) == want_down].index
    return [g for g in predicted_targets if g in set(moved)]
```

## Seed match analysis and site types

**Goal:** Locate seed matches in a UTR and classify site strength.

**Approach:** The seed is miRNA positions 2-7; a site is the reverse complement of the seed in the 3' UTR. Canonical sites by decreasing efficacy: 8mer > 7mer-m8 > 7mer-A1 > 6mer.

```python
from Bio.Seq import Seq

def find_seed_matches(mirna_seq, utr_seq):
    # 7mer-m8 site = reverse complement of miRNA positions 2-8 found in the UTR
    seed = str(Seq(mirna_seq)[1:8])
    site = str(Seq(seed).reverse_complement())
    matches, start = [], 0
    while True:
        pos = utr_seq.find(site, start)
        if pos == -1:
            break
        matches.append(pos)
        start = pos + 1
    return matches
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `KeyError: 'miRNA family'` on TargetScan file | Wrong column name for the context-scores file | The miRNA column is `Mirbase ID`; rank by `weighted context++ score` |
| Hundreds of "targets", almost none real | Treating seed prediction as truth | Intersect with anti-correlated mRNA DE; anchor in miRTarBase strong evidence |
| Every miRNA "regulates cancer pathways" | Enrichment on an unfiltered predicted target list (circular) | Build the list from validated/CLIP or expression-filtered targets before enrichment |
| Five tools "agree" so a target is trusted | All five use the seed (pseudo-replication) | Require an orthogonal evidence tier (CLIP/validated/expression), not more seed tools |
| A strong single-target claim | miRNAs repress most targets < 2-fold | Treat large single-target effects skeptically; demand validation |
| ceRNA/sponge mechanism asserted | Stoichiometry usually too low to matter (Denzler) | Require absolute abundance (miRNA copies vs added sites) before accepting it |

## Related Skills

- differential-mirna - Source of the DE miRNAs to predict targets for
- pathway-analysis/go-enrichment - Enrich a target list (only after evidence-filtering)
- database-access/entrez-fetch - Fetch UTR/gene sequences and identifiers
- clip-seq/ago-clip-mirna-targets - AGO-CLIP / CLASH direct target evidence

## References

- Agarwal V, Bell GW, Nam JW, Bartel DP. 2015. Predicting effective microRNA target sites in mammalian mRNAs. *eLife* 4:e05005. doi:10.7554/eLife.05005
- Betel D, Koppal A, Agius P, Sander C, Leslie C. 2010. Comprehensive modeling of microRNA targets predicts functional non-conserved and non-canonical sites. *Genome Biol* 11:R90. doi:10.1186/gb-2010-11-8-r90
- Chen Y, Wang X. 2020. miRDB: an online database for prediction of functional microRNA targets. *Nucleic Acids Res* 48:D127-D131. doi:10.1093/nar/gkz757
- Huang HY, Lin YC, Cui S, et al. 2022. miRTarBase update 2022: an informative resource for experimentally validated miRNA-target interactions. *Nucleic Acids Res* 50:D222-D230. doi:10.1093/nar/gkab1079
- Ru Y, Kechris KJ, Tabakoff B, et al. 2014. The multiMiR R package and database: integration of microRNA-target interactions. *Nucleic Acids Res* 42:e133. doi:10.1093/nar/gku631
- Pinzón N, Li B, Martinez L, et al. 2017. microRNA target prediction programs predict many false positives. *Genome Res* 27:234-245. doi:10.1101/gr.205146.116
- Baek D, Villén J, Shin C, et al. 2008. The impact of microRNAs on protein output. *Nature* 455:64-71. doi:10.1038/nature07242
- Selbach M, Schwanhäusser B, Thierfelder N, et al. 2008. Widespread changes in protein synthesis induced by microRNAs. *Nature* 455:58-63. doi:10.1038/nature07228
- Helwak A, Kudla G, Dudnakova T, Tollervey D. 2013. Mapping the human miRNA interactome by CLASH reveals frequent noncanonical binding. *Cell* 153:654-665. doi:10.1016/j.cell.2013.03.043
- Denzler R, Agarwal V, Stefano J, Bartel DP, Stoffel M. 2014. Assessing the ceRNA hypothesis with quantitative measurements of miRNA and target abundance. *Mol Cell* 54:766-776. doi:10.1016/j.molcel.2014.03.045
