---
name: bio-workflows-genome-annotation-pipeline
description: Orchestrates genome annotation from assembled contigs to functional annotation, forking prokaryotic (Bakta one-step, genetic-code table from GTDB-Tk) vs eukaryotic (RepeatMask -> BRAKER3 -> functional -> ncRNA), then eggNOG/InterProScan functional assignment and Infernal/tRNAscan ncRNA. Use when committing the pro-vs-eukaryotic path and the genetic-code table from taxonomy (never guessing), annotating ONLY a decontaminated QC-passed assembly (CheckM2 before prokaryotic annotation is non-negotiable), committing the evidence set (RNA-seq + protein drives BRAKER3 training), soft-masking with a curated repeat library before gene prediction, or pinning the tool + DB version for any pangenome comparison. Hands mechanism to the genome-annotation component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: Bakta
goal_approach_exempt: true
workflow: true
depends_on:
  - genome-annotation/prokaryotic-annotation
  - genome-annotation/eukaryotic-gene-prediction
  - genome-annotation/repeat-annotation
  - genome-annotation/functional-annotation
  - genome-annotation/ncrna-annotation
  - genome-annotation/annotation-qc
  - genome-assembly/assembly-qc
qc_checkpoints:
  - after_repeat_masking: "Repeat content within expected range for taxon"
  - after_gene_prediction: "Gene count plausible, BUSCO completeness >90%"
  - after_functional_annotation: ">60% of genes with functional assignment"
---

## Version Compatibility

Reference examples tested with: BRAKER3 3.0+, BUSCO 5.5+, Bakta 1.9+, Infernal 1.1+, InterProScan 5.66+, Prokka 1.14+, RepeatMasker 4.1+, RepeatModeler 2.0.4+ (-threads replaced -pa in 2.0.4), eggNOG-mapper 2.1+, pandas 2.2+, tRNAscan-SE 2.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Genome Annotation Pipeline

**"Annotate my genome assembly"** -> Orchestrate prokaryotic (Bakta) or eukaryotic (BRAKER3) gene prediction, repeat masking (RepeatMasker), functional annotation (eggNOG-mapper, InterProScan), and ncRNA annotation (Infernal).

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

A gene set is ~95% right and 100% confident; its trustworthiness is decided at four seams, not inside the gene-finder.

1. **Pro- vs eukaryotic is THE fork, committed from taxonomy up front, and it fixes the genetic-code table.** Prokaryote -> Bakta one-step (verify the genetic-code TABLE from GTDB-Tk classification, never guess — a Mycoplasma under table 11 splits every gene at internal UGA). Eukaryote -> multi-step RepeatMask -> BRAKER3 -> functional -> ncRNA. There is no general-purpose eukaryote annotator: alternative genetic codes, trans-splicing, and polycistronic transcription break standard pipelines.
2. **Annotate ONLY a decontaminated, QC-passed assembly.** CheckM2 before prokaryotic annotation is non-negotiable: contamination >5% mixes two organisms' genes into one chimeric set; a gene-finder trained on a contaminated/fragmented assembly produces confidently-wrong models genome-wide that are invisible in the GFF3. Annotation quality is bounded above by assembly quality.
3. **The evidence set is a committed input, not an afterthought.** Eukaryotic: RNA-seq BAM + protein (OrthoDB) evidence drives BRAKER3's high-confidence training-set mining (the real advance — learning from loci where transcripts AND homology agree). Committing RNA-seq (ideally Iso-Seq for isoforms+UTRs) is decided at project design; without it the annotation is one-isoform, CDS-only, UTR-less and silently poisons AS/3'-tag/APA analyses.
4. **Tool + DB version + date is a reproducibility commitment.** Bakta's DB is versioned (record it); Prokka's is frozen ~2019 (a post-2019 gene is "hypothetical" in Prokka, "named" in Bakta — accessory-vs-core flips on tool vintage alone). For any comparison, re-annotate everyone with ONE pipeline + ONE DB version from FASTA.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Pro- vs eukaryotic path + genetic-code table (from taxonomy) | The whole tool chain; a wrong code table splits genes at recoded stops |
| Decontaminated, QC-passed assembly (CheckM2 gate) | Chimeric gene set / genome-wide corrupt training if skipped; annotation quality is bounded by assembly quality |
| Evidence set (RNA-seq + protein) | BRAKER3 training quality; without RNA-seq the annotation is isoform-naive, UTR-less |
| Tool + DB version | Named-vs-hypothetical and accessory-vs-core flip on tool vintage; re-annotate all with one version for comparison |

## Pipeline Overview

```
Assembled contigs
    |
    v
[0. Assembly QC] ----------> QUAST, BUSCO (confirm assembly quality)
    |
    +----- Prokaryotic? -----> Path A: Bakta (one-step annotation)
    |                                |
    |                                v
    |                          Annotated genome (GFF3, GenBank, FASTA)
    |
    +----- Eukaryotic? ------> Path B: Multi-step pipeline
                                    |
                                    v
                              [1. Repeat Masking] ----> RepeatModeler + RepeatMasker
                                    |
                                    v
                              [2. Gene Prediction] ---> BRAKER3 (RNA-seq + protein evidence)
                                    |
                                    v
                              [3. Functional Annotation] -> eggNOG-mapper + InterProScan
                                    |
                                    v
                              [4. ncRNA Annotation] ---> Infernal + tRNAscan-SE
                                    |
                                    v
                              Annotated genome (GFF3, proteins, functional tables)
```

## Path A: Prokaryotic Annotation (Bakta)

Bakta provides comprehensive one-step annotation for bacteria and archaea. Preferred over Prokka for new projects.

### Database Setup

```bash
bakta_db download --output /path/to/bakta_db --type full
```

### Run Bakta

```bash
bakta \
    --db /path/to/bakta_db \
    --output bakta_out \
    --prefix my_genome \
    --locus-tag MYORG \
    --genus Escherichia --species "coli" \
    --strain K12 \
    --gram - \
    --translation-table 11 \
    --threads 8 \
    assembly.fasta
# Set --translation-table from the GTDB-Tk classification, never a guess: table 11 for most
# bacteria, but --translation-table 4 for Mycoplasma/Spiroplasma (UGA = Trp, not stop) --
# annotating a Mycoplasma under table 11 splits every gene at its internal UGA codons.
# Add --complete ONLY for finished replicons; omit it for draft contigs (the common input).
```

### Prokaryotic QC Checkpoint

```python
import subprocess
import json

def validate_prokaryotic_annotation(bakta_dir, prefix, expected_cds_range=(500, 8000)):
    '''
    QC gates for prokaryotic annotation.
    - CDS count in expected range for genome size
    - tRNA count >= 20 (typical minimum for free-living bacteria)
    - rRNA operons detected
    '''
    gff_file = f'{bakta_dir}/{prefix}.gff3'

    feature_counts = {'CDS': 0, 'tRNA': 0, 'rRNA': 0, 'tmRNA': 0, 'ncRNA': 0}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3 and fields[2] in feature_counts:
                feature_counts[fields[2]] += 1

    qc_pass = True
    if not (expected_cds_range[0] <= feature_counts['CDS'] <= expected_cds_range[1]):
        print(f'WARNING: CDS count {feature_counts["CDS"]} outside expected range {expected_cds_range}')
        qc_pass = False
    if feature_counts['tRNA'] < 20:
        print(f'WARNING: Only {feature_counts["tRNA"]} tRNAs detected (expect >= 20)')
        qc_pass = False

    print(f'Feature summary: {feature_counts}')
    return qc_pass, feature_counts
```

## Path B: Eukaryotic Annotation

### Step 1: Repeat Masking

```bash
# Build the RepeatModeler database FIRST, then the species-specific library
BuildDatabase -name mygenome assembly.fasta
RepeatModeler -database mygenome -threads 8 -LTRStruct

# CURATE the de novo library against a protein DB before masking, or real multi-copy gene
# families (NLR/R-genes, zinc-fingers) get masked and "discovered" as a gene-poor repertoire.
# Then soft-mask with the curated library (RepeatMasker uses the bundled Dfam DB in addition).
RepeatMasker \
    -lib mygenome-families.fa \
    -pa 8 \
    -xsmall \
    -gff \
    -dir repeat_out \
    assembly.fasta
```

#### Repeat Masking QC Checkpoint

```python
def check_repeat_content(repeatmasker_tbl, taxon='vertebrate'):
    '''
    Verify repeat content is within expected range for taxon.
    Typical ranges:
    - Vertebrate: 30-60%
    - Insect: 15-45%
    - Plant: 20-85%
    - Fungus: 3-20%
    '''
    expected_ranges = {
        'vertebrate': (30, 60), 'insect': (15, 45),
        'plant': (20, 85), 'fungus': (3, 20)
    }
    low, high = expected_ranges.get(taxon, (5, 80))

    with open(repeatmasker_tbl) as f:
        for line in f:
            if 'total interspersed' in line.lower():
                pct = float(line.strip().split()[-1].replace('%', ''))
                break

    qc_pass = low <= pct <= high
    if not qc_pass:
        print(f'WARNING: Repeat content {pct:.1f}% outside expected range ({low}-{high}%) for {taxon}')
    return qc_pass, pct
```

### Step 2: Gene Prediction with BRAKER3

```bash
# BRAKER3 combines GeneMark-ETP, AUGUSTUS, and TSEBRA
# Uses both RNA-seq and protein evidence for best results
braker.pl \
    --genome=repeat_out/assembly.fasta.masked \
    --bam=rnaseq_sorted.bam \
    --prot_seq=proteins.fa \
    --softmasking \
    --threads 8 \
    --species=my_species \
    --gff3 \
    --workingdir=braker_out

# If only RNA-seq evidence available
braker.pl \
    --genome=repeat_out/assembly.fasta.masked \
    --bam=rnaseq_sorted.bam \
    --softmasking \
    --threads 8 \
    --species=my_species \
    --gff3

# If only protein evidence available (use OrthoDB proteins)
braker.pl \
    --genome=repeat_out/assembly.fasta.masked \
    --prot_seq=orthodb_proteins.fa \
    --softmasking \
    --threads 8 \
    --species=my_species \
    --gff3
```

#### Gene Prediction QC Checkpoint

```bash
# BUSCO completeness on predicted proteins. Use the DEEPEST applicable clade dataset
# (e.g. insecta_odb10 / embryophyta_odb10), NOT the shallow eukaryota_odb10.
# The diagnostic that matters: compare this proteome BUSCO to a genome-mode BUSCO on the
# same assembly -- a large gap means the predictor missed present genes (see genome-annotation/annotation-qc).
busco \
    -i braker_out/braker.aa \
    -l <clade>_odb10 \
    -o busco_annotation \
    -m proteins \
    --cpu 8
```

```python
def check_gene_prediction(braker_gff, busco_summary, expected_genes_range=(15000, 35000)):
    '''
    QC gates after gene prediction.
    - Gene count within expected range for genome
    - BUSCO completeness > 90%
    - Mean exons per gene > 1 (spliced genes expected in eukaryotes)
    '''
    gene_count = 0
    exon_count = 0
    with open(braker_gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            feature = line.strip().split('\t')[2] if len(line.strip().split('\t')) >= 3 else ''
            if feature == 'gene':
                gene_count += 1
            elif feature == 'exon':
                exon_count += 1

    mean_exons = exon_count / gene_count if gene_count > 0 else 0

    with open(busco_summary) as f:
        for line in f:
            if line.strip().startswith('C:'):
                completeness = float(line.strip().split('C:')[1].split('%')[0])
                break

    issues = []
    if not (expected_genes_range[0] <= gene_count <= expected_genes_range[1]):
        issues.append(f'Gene count {gene_count} outside expected range {expected_genes_range}')
    if completeness < 90:
        issues.append(f'BUSCO completeness {completeness:.1f}% < 90%')
    if mean_exons < 2:
        issues.append(f'Mean exons/gene {mean_exons:.1f} is low for eukaryote')

    print(f'Genes: {gene_count}, Mean exons/gene: {mean_exons:.1f}, BUSCO: {completeness:.1f}%')
    return len(issues) == 0, issues
```

### Step 3: Functional Annotation

```bash
# eggNOG-mapper for comprehensive functional annotation
emapper.py \
    -i braker_out/braker.aa \
    --output eggnog_results \
    --cpu 8 \
    -m diamond \
    --tax_scope auto \
    --go_evidence non-electronic \
    --target_orthologs all \
    --seed_ortholog_evalue 1e-5 \
    --override

# InterProScan for domain annotation (complementary to eggNOG)
interproscan.sh \
    -i braker_out/braker.aa \
    -b interpro_results \
    -f tsv,gff3 \
    -goterms \
    -pa \
    -cpu 8
```

#### Functional Annotation QC Checkpoint

```python
import pandas as pd

def check_functional_annotation(eggnog_annotations, total_genes):
    '''
    QC gate: > 60% of genes should have functional assignment.
    Below 50% suggests database issues or highly divergent organism.
    '''
    cols = ['query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl',
            'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko']
    df = pd.read_csv(eggnog_annotations, sep='\t', comment='#', header=None)
    df.columns = (cols + [f'c{i}' for i in range(len(df.columns) - len(cols))])[:len(df.columns)]
    annotated = len(df[df['Description'] != '-'])
    pct_annotated = annotated / total_genes * 100

    has_go = len(df[df['GOs'] != '-'])
    has_kegg = len(df[df['KEGG_ko'] != '-'])

    print(f'Annotated: {annotated}/{total_genes} ({pct_annotated:.1f}%)')
    print(f'With GO terms: {has_go}, With KEGG: {has_kegg}')

    if pct_annotated < 60:
        print('WARNING: <60% annotated. Check database version or use broader taxonomy scope.')
    return pct_annotated >= 60
```

### Step 4: ncRNA Annotation

```bash
# tRNAscan-SE for tRNA genes
tRNAscan-SE \
    -E \
    --thread 8 \
    -o trna_results.txt \
    --gff trna.gff \
    assembly.fasta

# Infernal for Rfam-based ncRNA annotation. Rfam.cm ships pre-calibrated: cmpress it, never recalibrate.
# --cut_ga uses the curated per-family bit-score gathering thresholds (the correct Rfam default over a
# flat E-value); --rfam is the large-DB strict filter; --nohmmonly keeps GA valid for every model.
cmpress Rfam.cm
cmscan \
    --cpu 8 \
    --cut_ga --rfam --nohmmonly \
    --tblout rfam_results.tbl \
    --fmt 2 \
    --clanin Rfam.clanin \
    Rfam.cm \
    assembly.fasta
# Clan deoverlapping: drop hits marked '=' (dominated by a higher-scoring clanmate)
grep -v ' = ' rfam_results.tbl > rfam_results.deoverlapped.tbl
```

## Merging Annotations

```python
def merge_annotations(braker_gff, trna_gff, rfam_tbl, eggnog_tsv, output_gff):
    '''Merge gene predictions, ncRNAs, and functional annotations into final GFF3.'''
    import subprocess

    # Use AGAT for GFF merging and validation
    subprocess.run([
        'agat_sp_merge_annotations.pl',
        '--gff', braker_gff,
        '--gff', trna_gff,
        '-o', output_gff
    ], check=True)

    # Validate final GFF3
    subprocess.run([
        'agat_sp_statistics.pl',
        '--gff', output_gff,
        '-o', output_gff.replace('.gff3', '_stats.txt')
    ], check=True)

    print(f'Merged annotations written to {output_gff}')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Chimeric gene set; training corrupted genome-wide | Annotated a contaminated assembly | CheckM2/FCS-GX gate BEFORE annotation; decontaminate first |
| Genes split at recoded stops; low coding density, high hypothetical | Wrong genetic-code table | Set the table from GTDB-Tk taxonomy, not a guess |
| Real NLR/immune gene families deleted; suspiciously gene-poor | Over-masking with an uncurated repeat library | Filter the RepeatModeler library against a protein DB; confirm conserved families survive; soft-mask `-xsmall` |
| Accessory genome inflated ~10x in a pangenome | Frozen-DB annotation drift (Panaroo) | Re-annotate all assemblies with ONE pipeline + DB version from FASTA; existing GenBank annotations are unusable for pangenomics |
| BUSCO looks great but models are wrong | BUSCO-only quality claim (certifies ~1000 easy genes) | Proteome-mode BUSCO on delivered proteins; add mono-exonic fraction + length distribution + mRNA:gene ratio |
| Low gene count | Repeat masking too aggressive | Soft-mask (`-xsmall`) not hard-mask; curate the TE library first |
| Isoform-naive annotation (mRNA:gene = 1.00) | No RNA-seq evidence | Add RNA-seq (Iso-Seq for isoforms/UTRs) to BRAKER3 |

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

GENOME="assembly.fasta"
RNASEQ_BAM="rnaseq_sorted.bam"
PROTEINS="orthodb_proteins.fa"
BAKTA_DB="/path/to/bakta_db"
THREADS=8

# Determine organism type
ORGANISM_TYPE="${1:-eukaryotic}"  # prokaryotic or eukaryotic

if [ "$ORGANISM_TYPE" == "prokaryotic" ]; then
    echo "Running prokaryotic annotation with Bakta"
    echo "Step 0: Contamination/completeness gate (CheckM2 before annotation is non-negotiable)"
    checkm2 predict --input $GENOME --output-directory checkm2_out --threads $THREADS
    # Inspect checkm2_out/quality_report.tsv: proceed only if Completeness is high and Contamination < 5%
    # --translation-table comes from the GTDB-Tk classification, never a guess (rule 1). Bakta silently
    # assumes table 11; Mycoplasma/Spiroplasma need 4 (UGA = Trp, not stop) or every gene is truncated.
    bakta --db $BAKTA_DB --output bakta_out --prefix genome --translation-table 11 \
          --locus-tag MYORG --threads $THREADS $GENOME
    echo "Done. Results in bakta_out/"

else
    echo "Running eukaryotic annotation pipeline"

    echo "Step 0: Assembly QC (contiguity + completeness before committing to annotation)"
    quast.py $GENOME -o quast_out --threads $THREADS
    busco -i $GENOME -l "${LINEAGE:?set the DEEPEST applicable clade dataset, e.g. insecta_odb10 / embryophyta_odb10; eukaryota_odb10 is too shallow and inflates completeness}" -o busco_asm -m genome --cpu $THREADS

    echo "Step 1: Repeat masking"
    BuildDatabase -name mygenome $GENOME
    RepeatModeler -database mygenome -threads $THREADS -LTRStruct
    RepeatMasker -lib mygenome-families.fa -pa $THREADS -xsmall -gff -dir repeat_out $GENOME

    echo "Step 2: Gene prediction with BRAKER3"
    braker.pl --genome=repeat_out/$(basename $GENOME).masked \
              --bam=$RNASEQ_BAM --prot_seq=$PROTEINS \
              --softmasking --threads $THREADS --gff3 --workingdir=braker_out

    echo "Step 3: BUSCO QC (use the deepest applicable clade dataset, not eukaryota_odb10)"
    busco -i braker_out/braker.aa -l "${LINEAGE:?set the DEEPEST applicable clade dataset, e.g. insecta_odb10 / embryophyta_odb10; eukaryota_odb10 is too shallow and inflates completeness}" -o busco_check -m proteins --cpu $THREADS

    echo "Step 4: Functional annotation"
    emapper.py -i braker_out/braker.aa --output eggnog_out --cpu $THREADS -m diamond

    echo "Step 5: ncRNA annotation"
    tRNAscan-SE -E --thread $THREADS -o trna_out.txt --gff trna.gff $GENOME
    cmscan --cpu $THREADS --cut_ga --rfam --nohmmonly --tblout rfam.tbl --fmt 2 --clanin Rfam.clanin Rfam.cm $GENOME

    echo "Done. Check braker_out/, eggnog_out*, trna.gff, rfam.tbl"
fi
```

## Related Skills

- genome-annotation/prokaryotic-annotation - Bakta and Prokka details
- genome-annotation/eukaryotic-gene-prediction - BRAKER3 and AUGUSTUS options
- genome-annotation/repeat-annotation - Soft-masking before gene prediction
- genome-annotation/functional-annotation - eggNOG-mapper and InterProScan
- genome-annotation/ncrna-annotation - Infernal/Rfam and tRNAscan-SE detail
- rna-structure/ncrna-search - Covariance-model search, gathering thresholds, and clan resolution
- genome-annotation/annotation-qc - BUSCO genome-vs-proteome, OMArk, CheckM2 gates
- genome-assembly/assembly-qc - Pre-annotation assembly quality checks
- genome-intervals/gtf-gff-handling - GFF3/GTF hierarchy traversal, AGAT sanitizing/validation, coordinate conversion, and seqid-consistency checks on the merged annotation
- workflows/genome-assembly-pipeline - Upstream: hands off the decontaminated, QC-passed FASTA (with its QV/BUSCO)

## References

- Salzberg SL (2019) Next-generation genome annotation: we still struggle to get it right. *Genome Biology* 20:92. DOI 10.1186/s13059-019-1715-2. (error propagation.)
- Gabriel L, Bruna T, Hoff KJ, et al (2024) BRAKER3: fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA. *Genome Research* 34:769-777. DOI 10.1101/gr.278090.123. (high-confidence training-set mining.)
- Tonkin-Hill G, MacAlasdair N, Ruis C, et al (2020) Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biology* 21:180. DOI 10.1186/s13059-020-02090-4. (annotation-drift accessory inflation.)
- Schwengers O, Jelonek L, Dieckmann MA, et al (2021) Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics* 7:000685. DOI 10.1099/mgen.0.000685.
