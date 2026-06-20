---
name: bio-read-qc-rnaseq-qc
description: Runs RNA-seq-specific post-alignment QC - strandedness inference, gene-body 5'-3' coverage, read distribution (exonic/intronic/intergenic), rRNA/globin/mitochondrial rate, transcript integrity (TIN), and saturation - with RSeQC, Qualimap, RNA-SeQC, and Picard. Use when validating RNA-seq libraries before quantification or differential expression, diagnosing degradation or gDNA contamination, or determining library strandedness. For raw-FASTQ QC use quality-reports; for UMI dedup use umi-processing.
tool_type: mixed
primary_tool: RSeQC
---

## Version Compatibility

Reference examples tested with: RSeQC 5.0+, Qualimap 2.3+, RNA-SeQC 2.4+, Picard 3.1+, salmon 1.10+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# RNA-seq QC -- post-alignment metrics that have no DNA analogue

Assess strandedness, integrity, feature distribution, and enrichment on the ALIGNED BAM, using RSeQC / Qualimap / RNA-SeQC / Picard against a gene model.

**"Run RNA-seq QC"** -> Infer strandedness, gene-body coverage, exonic/intronic/intergenic distribution, rRNA rate, and TIN from the BAM.
- CLI: `infer_experiment.py -i aligned.bam -r genes.bed12` (strandedness)
- CLI: `picard CollectRnaSeqMetrics` / `qualimap rnaseq` / `rnaseqc collapsed.gtf in.bam out/`

Scope: this skill OWNS transcriptome QC on the aligned BAM. Raw-FASTQ QC (adapters, base quality) -> read-qc/quality-reports. UMI dedup -> read-qc/umi-processing. Quantification -> rna-quantification/featurecounts-counting. OUT OF SCOPE: differential expression (differential-expression/deseq2-basics).

## The Single Most Important Modern Insight

1. **These are POST-ALIGNMENT QC: every metric needs an aligned BAM AND a gene model (BED12 / GTF / refFlat / collapsed-GTF), which is the line that separates them from FastQC.** FastQC answers "is the sequencer output clean?"; RNA-seq QC answers "did I sequence the transcriptome I think I sequenced, in the orientation I think, with the integrity I think?" The metrics below (strand, exonic rate, rRNA rate, 5'-3' bias) have no DNA analogue because DNA has no exons, no strand of transcription, and no rRNA fraction. The most common setup error is feeding the wrong gene-model format (RSeQC wants BED12; Qualimap a GTF; RNA-SeQC a COLLAPSED GTF; Picard a refFlat + ribosomal_intervals).

2. **Getting strandedness wrong SILENTLY HALVES OR ZEROS the counts -- no error is thrown.** dUTP (TruSeq Stranded mRNA, most rRNA-depletion kits) is fr-firststrand = REVERSE = featureCounts `-s 2` = htseq `reverse` = salmon `ISR` = STAR ReadsPerGene column 4. Run it as "forward" and reads land on the antisense gene: counts collapse toward zero and the antisense neighbor inflates (running stranded data as UNSTRANDED, by contrast, roughly doubles counts). The tell is a huge "assigned to no feature" fraction or counts ~2x below the unstranded run. ALWAYS infer strandedness empirically (`infer_experiment.py`, salmon `-l A`, or how_are_we_stranded_here) before quantifying -- never assume from the kit name.

3. **In standard bulk RNA-seq WITHOUT UMIs, do NOT mark or remove duplicates.** A highly expressed gene legitimately produces many fragments sharing identical coordinates; at the read level a PCR duplicate and a natural duplicate are INDISTINGUISHABLE. Coordinate dedup (Picard MarkDuplicates) preferentially deletes reads from the most abundant and shortest transcripts, introducing an expression- and length-dependent bias. This is the OPPOSITE of DNA-seq. Duplication rate is a DIAGNOSTIC ("low complexity / over-sequenced / low input"), never a remove step. The only correct way to remove RNA PCR duplicates is UMIs (read-qc/umi-processing); UMI-protocol RNA-seq (QuantSeq, 10x) inverts the rule.

Integrity bonus: RIN is an electrophoresis estimate measured BEFORE library prep; gene-body coverage and TIN are the post-hoc TRUTH measured from the aligned reads. Use DV200 (% fragments >200 nt), not RIN, for FFPE/archival. In a cohort with variable quality, regress medTIN out as a covariate rather than discarding samples.

## Tool Taxonomy

| Tool | Gene model | Role |
|------|-----------|------|
| RSeQC | BED12 | The script suite: infer_experiment, geneBody_coverage, read_distribution, tin, junction_saturation, read_duplication |
| Qualimap 2 | GTF | `qualimap rnaseq`: feature distribution + transcript 5'-3' profile + junctions in one HTML (bamqc is the generic, non-RNA mode) |
| RNA-SeQC 2 | COLLAPSED GTF | GTEx/TOPMed tool; scales to tens of thousands of samples; exonic/intronic/intergenic + rRNA rate + TPM |
| Picard CollectRnaSeqMetrics | refFlat + ribosomal_intervals | PCT_CODING/UTR/INTRONIC/INTERGENIC/RIBOSOMAL, MEDIAN_5PRIME_TO_3PRIME_BIAS (cannot compute rRNA without the intervals) |
| SortMeRNA | rRNA database | Filter/quantify rRNA reads directly |

QC-gate order: (1) FastQC on raw FASTQ -> (2) align (STAR/HISAT2) -> (3) post-alignment QC: strandedness FIRST (it gates correct quantification), then read distribution, gene-body + TIN, rRNA/globin/MT, duplication + saturation -> (4) aggregate with MultiQC and judge each sample against the cohort.

## Strandedness -- infer, then set every tool to match

```bash
infer_experiment.py -i aligned.bam -r genes.bed12     # samples reads, reports the two fractions
salmon quant -i index -l A -r sample.fq.gz -o quant/  # -l A auto-detects; see lib_format_counts.json
```

| Protocol | infer_experiment dominant fraction | salmon -l (PE/SE) | featureCounts -s | htseq | STAR ReadsPerGene col |
|----------|------------------------------------|-------------------|------------------|-------|-----------------------|
| Unstranded | both ~0.5 | IU / U | 0 | no | 2 |
| fr-secondstrand (forward) | "1++,1--,2+-,2-+" | ISF / SF | 1 | yes | 3 |
| fr-firststrand (reverse, dUTP -- common) | "1+-,1-+,2++,2--" | ISR / SR | 2 | reverse | 4 |

Single-end infer_experiment drops the read-number prefix: forward = "++,--", reverse = "+-,-+". A STAR sanity check: the ReadsPerGene column with the most counts and fewest N_noFeature is the correct strand (the wrong column makes N_noFeature blow up). Picard STRAND_SPECIFICITY is a notorious inversion: NONE / FIRST_READ_TRANSCRIPTION_STRAND (= forward/fr-secondstrand) / SECOND_READ_TRANSCRIPTION_STRAND (= dUTP/reverse/fr-firststrand, the common case).

## Gene-body coverage and integrity

```bash
geneBody_coverage.py -i aligned.bam -r genes.bed12 -o coverage   # 5'->3' uniformity curve
tin.py -i aligned.bam -r genes.bed12 > tin.txt                   # per-transcript integrity; medTIN = sample score
```

3' bias (coverage piling at the 3' end) = RNA degradation OR oligo-dT priming of degraded/FFPE RNA -- which is why poly-A protocols fail on FFPE and rRNA-depletion + random priming is preferred there. 5' bias is rarer (5'-capture protocols / artifacts). Flat = intact RNA. RIN/DV200/TIN: RIN (1-10, pre-prep, electrophoresis) predicts degradation; DV200 (% >200 nt) is the FFPE metric because fragmented RNA has no rRNA peaks for RIN; TIN is measured from the data and can be used as a DE covariate.

## Read distribution and enrichment

```bash
read_distribution.py -i aligned.bam -r genes.bed12 > distribution.txt
```

- High INTRONIC = pre-mRNA / nuclear RNA or gDNA contamination (in snRNA-seq it is SIGNAL, not a fail).
- High INTERGENIC = gDNA contamination or annotation gaps. gDNA drives intronic AND intergenic up together; an annotation gap drives only intergenic.
- rRNA rate = the readout of poly-A-selection / rRNA-depletion efficiency (high = wasted reads, failed depletion).
- Globin (HBA/HBB) crowds whole-blood PAXgene libraries -- deplete (GLOBINclear); globin% is the readout.
- Mitochondrial %: high = degradation (bulk) or dying cells / ambient contamination (single-cell; in snRNA-seq it should be LOW).

## Duplication and saturation -- diagnostic, not a remove step

```bash
# Duplication as a DIAGNOSTIC only -- do NOT remove duplicates in non-UMI bulk RNA-seq
read_duplication.py -i aligned.bam -o dup                         # sequence- and mapping-based curves
junction_saturation.py -i aligned.bam -r genes.bed12 -o junc_sat  # enough depth for splicing?
```

## Complete QC pipeline

**Goal:** Produce a per-sample RNA-seq QC summary covering strandedness, distribution, integrity, and Picard metrics.

**Approach:** Infer strandedness first, run the RSeQC suite, then Picard with STRAND_SPECIFICITY set to the inferred protocol, and append to one report (do NOT dedup).

```bash
#!/bin/bash
set -euo pipefail
SAMPLE=$1; BAM=$2; BED12=$3; REFFLAT=$4; RRNA_INTERVALS=$5
STRAND=${6:-SECOND_READ_TRANSCRIPTION_STRAND}   # SECOND = dUTP/reverse (common); FIRST = forward; NONE = unstranded

REPORT="${SAMPLE}_rnaseq_qc.txt"
echo "=== RNA-seq QC: $SAMPLE ===" > "$REPORT"

echo "--- Strandedness (set downstream tools to match) ---" >> "$REPORT"
infer_experiment.py -i "$BAM" -r "$BED12" >> "$REPORT"

echo "--- Read distribution ---" >> "$REPORT"
read_distribution.py -i "$BAM" -r "$BED12" >> "$REPORT"

geneBody_coverage.py -i "$BAM" -r "$BED12" -o "${SAMPLE}_genebody"
tin.py -i "$BAM" -r "$BED12"                       # writes <bam>.summary.txt (mean/median TIN) + <bam>.tin.xls
echo "--- TIN (medTIN = median column of the summary) ---" >> "$REPORT"
cat *.summary.txt >> "$REPORT" 2>/dev/null

echo "--- Picard RNA-seq metrics (STRAND=$STRAND) ---" >> "$REPORT"
picard CollectRnaSeqMetrics I="$BAM" O="${SAMPLE}_picard.txt" \
    REF_FLAT="$REFFLAT" STRAND_SPECIFICITY="$STRAND" RIBOSOMAL_INTERVALS="$RRNA_INTERVALS"

cat "$REPORT"
```

## The collapsed gene model

A standard GTF lists many overlapping isoforms per gene, so a read that is exonic in isoform A but intronic in B is ambiguous and overlapping isoforms double-count the same base. RNA-SeQC 2 REQUIRES a COLLAPSED model (one flattened transcript per gene, inter-gene overlaps excluded), built with GTEx `collapse_annotation.py`. Mismatched or un-collapsed models are a leading cause of "my exonic rate looks wrong". Picard `PCT_*` metrics are FRACTIONS (0-1), not percentages, despite the name.

## Quantitative Thresholds

| Metric | Anchor | Source / rationale |
|--------|--------|--------------------|
| Mapping rate | > 0.2 exclude below (GTEx); > 85% typical | GTEx v8 RNA-SeQC gate |
| Intergenic rate | < 0.3 | GTEx; above = gDNA / annotation |
| rRNA rate | < 0.3 (GTEx); <5% polyA, <10% depleted in practice | depletion efficiency |
| Uniquely mapped reads | >= 30M (ENCODE human) | ENCODE long-RNA standard |
| medTIN | > 70 good, 50-70 moderate, < 50 poor | RSeQC TIN |
| 5'-to-3' bias | near 1 flat; > 2 strong degradation | Picard MEDIAN_5PRIME_TO_3PRIME_BIAS |

Thresholds are protocol-specific: an intronic rate that fails a poly-A bulk sample is normal/required for snRNA-seq (nuclei are >50% intronic); a 3' bias that condemns fresh poly-A is expected for FFPE. Apply cohort-relative outlier logic on top.

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Counts ~halved / huge "no feature" fraction | Wrong strandedness | Infer first; set featureCounts/htseq/salmon/Picard to match |
| RNA-seq DE has odd length bias | Marked duplicates on non-UMI bulk RNA-seq | Do not dedup; report duplication as a diagnostic |
| Exonic rate looks wrong in RNA-SeQC | Un-collapsed multi-isoform GTF | Use a collapsed GTF (GTEx collapse_annotation.py) |
| Picard rRNA metric is 0/blank | No ribosomal_intervals supplied | Build the interval list from rRNA features + BAM dict |
| snRNA-seq "fails" high intronic rate | Bulk gate applied to nuclear RNA | Intronic reads are signal in snRNA; use an intron-inclusive reference |
| Picard percentages look 100x too small | PCT_* are fractions (0-1) | Multiply by 100 for display |

## References

Wang L, Wang S, Li W. 2012. RSeQC: quality control of RNA-seq experiments. Bioinformatics 28(16):2184-2185.
Okonechnikov K, Conesa A, Garcia-Alcalde F. 2016. Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics 32(2):292-294.
Graubert A, Aguet F, Ravi A, Ardlie KG, Getz G. 2021. RNA-SeQC 2: efficient RNA-seq quality control and quantification for large cohorts. Bioinformatics 37(18):3048-3050.
Schroeder A, Mueller O, Stocker S, et al. 2006. The RIN: an RNA integrity number for assigning integrity values to RNA measurements. BMC Molecular Biology 7:3.
Wang L, Nie J, Sicotte H, et al. 2016. Measure transcript integrity using RNA-seq data. BMC Bioinformatics 17:58.
Smith T, Heger A, Sudbery I. 2017. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research 27(3):491-499.

## Related Skills

read-qc/quality-reports - Raw-FASTQ QC before alignment
read-qc/umi-processing - Molecule-accurate dedup for UMI RNA-seq
read-qc/contamination-screening - rRNA and cross-species contamination
read-alignment/star-alignment - Aligner that emits ReadsPerGene strandedness columns
rna-quantification/featurecounts-counting - Strand-aware quantification after QC
differential-expression/deseq2-basics - Use medTIN as a covariate in the design
