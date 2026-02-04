# New Categories Implementation Plan

## Categories to Add
1. **alternative-splicing** - Isoform and splicing analysis
2. **chemoinformatics** - Small molecule and drug discovery
3. **liquid-biopsy** - cfDNA and circulating tumor DNA analysis

## New Workflows to Add
4. **workflows/splicing-pipeline** - FASTQ to differential splicing
5. **workflows/drug-discovery-pipeline** - Target to virtual screening
6. **workflows/liquid-biopsy-pipeline** - cfDNA to tumor monitoring

---

# Category 1: alternative-splicing

## Overview
Alternative splicing analysis for RNA-seq data. Covers splicing quantification, differential splicing between conditions, isoform switching analysis, and visualization.

**Tool type:** mixed | **Primary tools:** rMATS, SUPPA2, leafcutter, IsoformSwitchAnalyzeR

## Skills (6 total)

### 1. splicing-quantification
```yaml
name: bio-splicing-quantification
description: Quantifies alternative splicing events (PSI/percent spliced in) from RNA-seq alignments using SUPPA2. Calculates inclusion levels for skipped exons, alternative 5'/3' splice sites, mutually exclusive exons, and retained introns. Use when measuring splice site usage or isoform ratios from BAM files.
tool_type: python
primary_tool: SUPPA2
```

**Key functions:**
- SUPPA2 generateEvents (extract events from GTF)
- SUPPA2 psiPerEvent (calculate PSI from TPM)
- rMATS event extraction
- PSI calculation formulas with documentation

**Example script:** `quantify_splicing.py`
- Input: BAM/TPM files, GTF annotation
- Output: PSI matrices per event type

### 2. differential-splicing
```yaml
name: bio-differential-splicing
description: Detects differential alternative splicing between conditions using rMATS, leafcutter, or SUPPA2 diffSplice. Identifies exon skipping, intron retention, and other splicing changes with statistical significance. Use when comparing splicing patterns between treatment groups, tissues, or disease states.
tool_type: mixed
primary_tool: rMATS
```

**Key functions:**
- rMATS (turbo) differential analysis
- leafcutter cluster + differential
- SUPPA2 diffSplice
- Multiple testing correction (FDR)
- Effect size (delta PSI) thresholds

**Example scripts:** `diff_splicing_rmats.sh`, `diff_splicing_leafcutter.R`

### 3. isoform-switching
```yaml
name: bio-isoform-switching
description: Analyzes isoform switching events and their functional consequences using IsoformSwitchAnalyzeR. Predicts protein domain changes, NMD sensitivity, and open reading frame alterations. Use when investigating how splicing changes affect protein function between conditions.
tool_type: r
primary_tool: IsoformSwitchAnalyzeR
```

**Key functions:**
- IsoformSwitchAnalyzeR workflow
- Domain prediction (Pfam)
- NMD sensitivity
- ORF analysis
- Consequence summary plots

**Example script:** `isoform_switch_analysis.R`

### 4. sashimi-plots
```yaml
name: bio-sashimi-plots
description: Creates sashimi plots showing RNA-seq read coverage and splice junction usage using ggsashimi or rmats2sashimiplot. Visualizes differential splicing events with junction read counts. Use when visualizing specific splicing events or validating differential splicing results.
tool_type: python
primary_tool: ggsashimi
```

**Key functions:**
- ggsashimi (Python wrapper)
- rmats2sashimiplot
- IGV sashimi
- Grouping by condition
- Custom color schemes

**Example script:** `plot_sashimi.py`

### 5. splicing-qc
```yaml
name: bio-splicing-qc
description: Assesses splicing analysis quality including junction saturation, splice site strength, and event coverage metrics using RSeQC and custom metrics. Use when evaluating RNA-seq data quality for splicing analysis or troubleshooting low event detection.
tool_type: python
primary_tool: RSeQC
```

**Key functions:**
- junction_saturation.py (RSeQC)
- Splice site strength (MaxEntScan)
- Event coverage filtering
- Sample correlation of PSI

**Example script:** `splicing_qc.py`

### 6. single-cell-splicing
```yaml
name: bio-single-cell-splicing
description: Analyzes alternative splicing at single-cell resolution using BRIE2 or Leafcutter2. Identifies cell-type-specific splicing patterns and splicing heterogeneity. Use when analyzing isoform usage in scRNA-seq data or finding splicing differences between cell populations.
tool_type: python
primary_tool: BRIE2
```

**Key functions:**
- BRIE2 for scRNA-seq
- Leafcutter2 pseudobulk
- Cell-type PSI aggregation
- Splicing velocity (preliminary)

**Example script:** `sc_splicing_brie2.py`

---

# Category 2: chemoinformatics

## Overview
Computational chemistry and drug discovery. Covers molecular representations, property prediction, similarity searching, virtual screening, and ADMET analysis.

**Tool type:** python | **Primary tools:** RDKit, Open Babel, DeepChem

## Skills (7 total)

### 1. molecular-io
```yaml
name: bio-molecular-io
description: Reads, writes, and converts molecular file formats (SMILES, SDF, MOL2, PDB) using RDKit and Open Babel. Handles chemical structure parsing, canonicalization, and format conversion. Use when loading chemical libraries, converting between formats, or standardizing molecular representations.
tool_type: python
primary_tool: RDKit
```

**Key functions:**
- Chem.MolFromSmiles, MolToSmiles
- Chem.SDMolSupplier, SDWriter
- Open Babel conversion
- SMILES canonicalization
- Mol sanitization and error handling

**Example script:** `molecule_io.py`

### 2. molecular-descriptors
```yaml
name: bio-molecular-descriptors
description: Calculates molecular descriptors and fingerprints (ECFP, MACCS, physicochemical properties) using RDKit. Computes Lipinski properties, topological indices, and 2D/3D descriptors. Use when featurizing molecules for machine learning or filtering by drug-likeness criteria.
tool_type: python
primary_tool: RDKit
```

**Key functions:**
- Descriptors.MolecularDescriptorCalculator
- AllChem.GetMorganFingerprintAsBitVect (ECFP)
- MACCSkeys.GenMACCSKeys
- Lipinski.* (HBD, HBA, LogP, MW)
- QED score

**Example script:** `calculate_descriptors.py`

### 3. similarity-searching
```yaml
name: bio-similarity-searching
description: Performs molecular similarity searches using Tanimoto coefficient on fingerprints via RDKit. Finds structurally similar compounds in chemical libraries using ECFP or MACCS keys. Use when finding analogs of a query compound or clustering molecules by structural similarity.
tool_type: python
primary_tool: RDKit
```

**Key functions:**
- DataStructs.TanimotoSimilarity
- BulkTanimotoSimilarity
- Similarity matrix computation
- Clustering (Butina)
- MCS (maximum common substructure)

**Example script:** `similarity_search.py`

### 4. substructure-search
```yaml
name: bio-substructure-search
description: Searches molecular libraries for substructure matches using SMARTS patterns with RDKit. Filters compounds by pharmacophore features or functional groups. Use when finding compounds containing specific chemical moieties or filtering libraries by structural features.
tool_type: python
primary_tool: RDKit
```

**Key functions:**
- Chem.MolFromSmarts
- mol.HasSubstructMatch
- GetSubstructMatches
- Common SMARTS patterns (rings, functional groups)
- Highlighting matches

**Example script:** `substructure_search.py`

### 5. admet-prediction
```yaml
name: bio-admet-prediction
description: Predicts ADMET (absorption, distribution, metabolism, excretion, toxicity) properties using DeepChem or ADMETlab models. Estimates bioavailability, CYP inhibition, hERG liability, and Ames mutagenicity. Use when filtering compounds for drug-likeness or prioritizing leads by predicted safety.
tool_type: python
primary_tool: DeepChem
```

**Key functions:**
- DeepChem ADMET models
- SwissADME API (if available)
- ADMETlab2.0 predictions
- Tox21 models
- PAINS filter

**Example script:** `predict_admet.py`

### 6. virtual-screening
```yaml
name: bio-virtual-screening
description: Performs structure-based virtual screening using AutoDock Vina or molecular docking workflows. Docks compound libraries against protein targets and ranks by predicted binding affinity. Use when screening chemical libraries against a protein structure to find potential binders.
tool_type: python
primary_tool: AutoDock Vina
```

**Key functions:**
- Vina docking setup
- Receptor preparation (PDBQT)
- Ligand preparation
- Batch docking
- Results parsing and ranking

**Example script:** `virtual_screen.py`

### 7. reaction-enumeration
```yaml
name: bio-reaction-enumeration
description: Enumerates chemical libraries through reaction SMARTS transformations using RDKit. Generates virtual compound libraries from building blocks using defined reactions. Use when creating combinatorial libraries or enumerating products from synthetic routes.
tool_type: python
primary_tool: RDKit
```

**Key functions:**
- AllChem.ReactionFromSmarts
- RunReactants
- Library enumeration
- Reaction validation
- Product filtering

**Example script:** `enumerate_reactions.py`

---

# Category 3: liquid-biopsy

## Overview
Cell-free DNA (cfDNA) and circulating tumor DNA (ctDNA) analysis for non-invasive cancer detection and monitoring. Covers fragment analysis, tumor fraction estimation, mutation detection, and methylation-based detection.

**Tool type:** mixed | **Primary tools:** ichorCNA, cfDNApipe, DELFI

## Skills (6 total)

### 1. cfdna-preprocessing
```yaml
name: bio-cfdna-preprocessing
description: Preprocesses cell-free DNA sequencing data including adapter trimming, alignment, and duplicate removal optimized for short cfDNA fragments. Applies cfDNA-specific QC thresholds. Use when processing plasma cfDNA sequencing data before downstream analysis.
tool_type: python
primary_tool: cfDNApipe
```

**Key functions:**
- cfDNApipe preprocessing
- Fragment length filtering (100-200bp typical)
- Duplicate marking (UMI-aware if applicable)
- cfDNA-specific mapping (short fragments)
- Coverage uniformity

**Example script:** `preprocess_cfdna.py`

### 2. fragment-analysis
```yaml
name: bio-fragment-analysis
description: Analyzes cfDNA fragment size distributions and fragmentomics features using DELFI or Griffin. Extracts nucleosome positioning and fragmentation patterns for cancer detection. Use when leveraging fragment length patterns for tumor detection or tissue-of-origin analysis.
tool_type: python
primary_tool: DELFI
```

**Key functions:**
- Fragment length distribution
- DELFI fragmentomics features
- Griffin coverage patterns
- Nucleosome footprinting
- Ratio of short to long fragments

**Example script:** `fragment_analysis.py`

### 3. tumor-fraction-estimation
```yaml
name: bio-tumor-fraction-estimation
description: Estimates circulating tumor DNA fraction from shallow whole-genome sequencing using ichorCNA. Detects copy number alterations and calculates ctDNA percentage in plasma. Use when quantifying tumor burden from liquid biopsy samples or monitoring treatment response.
tool_type: r
primary_tool: ichorCNA
```

**Key functions:**
- ichorCNA workflow
- Read depth binning
- GC correction
- HMM segmentation
- Tumor fraction from CNA

**Example script:** `estimate_tumor_fraction.R`

### 4. ctdna-mutation-detection
```yaml
name: bio-ctdna-mutation-detection
description: Detects somatic mutations in circulating tumor DNA using specialized variant callers optimized for low allele fractions. Applies error suppression and UMI-based deduplication. Use when identifying tumor mutations from plasma DNA at variant allele frequencies below 1%.
tool_type: python
primary_tool: VarDict
```

**Key functions:**
- VarDict low-AF calling
- UMI consensus with fgbio/Picard
- Error suppression strategies
- Panel-of-normals filtering
- VAF tracking over time

**Example script:** `detect_ctdna_mutations.py`

### 5. methylation-based-detection
```yaml
name: bio-methylation-based-detection
description: Analyzes cfDNA methylation patterns for cancer detection and tissue-of-origin using cfMeDIP-seq or bisulfite sequencing. Identifies cancer-specific methylation signatures. Use when using methylation biomarkers for early cancer detection or minimal residual disease monitoring.
tool_type: python
primary_tool: methylDackel
```

**Key functions:**
- cfMeDIP-seq processing
- Methylation deconvolution
- Tissue-of-origin prediction
- Cancer-specific DMR detection
- MCED (multi-cancer early detection) panels

**Example script:** `cfdna_methylation.py`

### 6. longitudinal-monitoring
```yaml
name: bio-longitudinal-monitoring
description: Tracks ctDNA dynamics over time for treatment response monitoring using serial liquid biopsy samples. Analyzes tumor fraction trends and mutation clearance kinetics. Use when monitoring patients during therapy or detecting molecular relapse before clinical progression.
tool_type: python
primary_tool: pandas
```

**Key functions:**
- Serial sample tracking
- Tumor fraction trend analysis
- Mutation persistence/clearance
- Molecular response criteria
- Time-to-progression modeling

**Example script:** `longitudinal_monitoring.py`

---

# New Workflows

## workflows/splicing-pipeline
```yaml
name: bio-splicing-pipeline
description: End-to-end alternative splicing analysis from FASTQ to differential splicing results. Aligns with STAR, quantifies events with rMATS/SUPPA2, and visualizes with sashimi plots. Use when performing comprehensive splicing analysis from raw RNA-seq data.
tool_type: mixed
primary_tool: rMATS
```

**Steps:**
1. Read QC and trimming (fastp)
2. STAR alignment (2-pass)
3. Splicing quantification (rMATS or SUPPA2)
4. Differential splicing analysis
5. Isoform switching (IsoformSwitchAnalyzeR)
6. Visualization (sashimi plots)
7. Functional annotation

## workflows/drug-discovery-pipeline
```yaml
name: bio-drug-discovery-pipeline
description: Virtual screening pipeline from target structure to ranked compound hits. Prepares receptor, filters compound library by drug-likeness, docks candidates, and predicts ADMET. Use when performing computational drug discovery against a protein target.
tool_type: python
primary_tool: RDKit
```

**Steps:**
1. Target preparation (PDB cleanup)
2. Binding site identification
3. Library filtering (Lipinski, PAINS)
4. Docking (AutoDock Vina)
5. Rescoring and filtering
6. ADMET prediction
7. Hit prioritization

## workflows/liquid-biopsy-pipeline
```yaml
name: bio-liquid-biopsy-pipeline
description: Cell-free DNA analysis pipeline from plasma sequencing to tumor monitoring. Processes cfDNA reads, estimates tumor fraction, detects mutations, and tracks longitudinal dynamics. Use when analyzing liquid biopsy samples for cancer detection or monitoring.
tool_type: mixed
primary_tool: ichorCNA
```

**Steps:**
1. cfDNA preprocessing
2. Fragment analysis (DELFI)
3. Tumor fraction estimation (ichorCNA)
4. Mutation detection (if panel/WES)
5. Methylation analysis (if applicable)
6. Longitudinal tracking
7. Clinical reporting

---

# Deprecation Checks Required

## Alternative Splicing
- [ ] rMATS-turbo vs legacy rMATS (use turbo)
- [ ] SUPPA2 current API
- [ ] leafcutter vs leafcutter2
- [ ] IsoformSwitchAnalyzeR version compatibility

## Chemoinformatics
- [ ] RDKit 2023+ API changes
- [ ] DeepChem current models
- [ ] Open Babel 3.x compatibility
- [ ] AutoDock Vina vs Vina 1.2

## Liquid Biopsy
- [ ] ichorCNA current version
- [ ] cfDNApipe maintenance status
- [ ] DELFI availability
- [ ] fgbio UMI tools current API

---

# Cross-Category References

## alternative-splicing references:
- read-alignment/star-alignment - STAR 2-pass for junctions
- differential-expression/deseq2-basics - Compare with gene-level
- rna-quantification/salmon-quantification - Transcript TPM input
- data-visualization/ggplot2-fundamentals - Plot customization

## chemoinformatics references:
- structural-biology/pdb-parsing - Protein structure handling
- structural-biology/modern-structure-prediction - Target structures
- machine-learning/biomarker-discovery - ML for compound activity

## liquid-biopsy references:
- variant-calling/variant-calling - Somatic mutation calling
- copy-number/cnvkit-analysis - CNV detection principles
- methylation-analysis/bismark-alignment - Methylation processing
- clinical-databases/clinvar-queries - Variant annotation

---

# Validation Checklist

Before implementation, verify:
- [ ] All names ≤64 chars, lowercase, hyphens only
- [ ] All descriptions include "Use when..."
- [ ] All descriptions are third-person
- [ ] primary_tool is single value
- [ ] SKILL.md will be under 500 lines
- [ ] usage-guide.md has all required sections
- [ ] examples/ has at least one script per skill
- [ ] All cross-references use qualified paths
- [ ] No deprecated functions used
