---
name: bio-proteomics-peptide-identification
description: Peptide-spectrum matching and protein identification from MS/MS data. Use when identifying peptides from tandem mass spectra. Covers database searching, spectral library matching, and FDR estimation using target-decoy approaches.
tool_type: mixed
primary_tool: pyOpenMS
---

# Peptide Identification

## Database Search with pyOpenMS

```python
from pyopenms import MSExperiment, MzMLFile, FASTAFile, ProteaseDigestion
from pyopenms import ModificationsDB, AASequence

# Load FASTA database
fasta_entries = []
FASTAFile().load('uniprot_human.fasta', fasta_entries)

# In-silico digestion
digestion = ProteaseDigestion()
digestion.setEnzyme('Trypsin')
digestion.setMissedCleavages(2)

peptides = []
for entry in fasta_entries:
    seq = AASequence.fromString(entry.sequence)
    result = []
    digestion.digest(seq, result)
    peptides.extend([(entry.identifier, str(p)) for p in result])
```

## Working with Search Results (idXML)

```python
from pyopenms import IdXMLFile, ProteinIdentification, PeptideIdentification

protein_ids = []
peptide_ids = []
IdXMLFile().load('search_results.idXML', protein_ids, peptide_ids)

for pep_id in peptide_ids:
    rt = pep_id.getRT()
    mz = pep_id.getMZ()
    for hit in pep_id.getHits():
        sequence = hit.getSequence()
        score = hit.getScore()
        charge = hit.getCharge()
```

## FDR Estimation (Target-Decoy)

```python
def calculate_fdr(scores, is_decoy, score_threshold):
    above_threshold = scores >= score_threshold
    n_target = ((~is_decoy) & above_threshold).sum()
    n_decoy = (is_decoy & above_threshold).sum()
    fdr = n_decoy / n_target if n_target > 0 else 1.0
    return fdr

def find_score_at_fdr(scores, is_decoy, target_fdr=0.01):
    sorted_scores = np.sort(scores)[::-1]
    for threshold in sorted_scores:
        fdr = calculate_fdr(scores, is_decoy, threshold)
        if fdr <= target_fdr:
            return threshold
    return sorted_scores[-1]
```

## R: Search Result Processing

```r
library(MSnbase)

# Read mzIdentML results
psms <- readMzIdData('results.mzid')

# Filter to 1% FDR
psms_filtered <- psms[psms$qvalue <= 0.01, ]

# Unique peptides per protein
peptide_counts <- table(psms_filtered$accession)
```

## Spectral Library Search

```python
from pyopenms import SpectraSTSearchAlgorithm, MSExperiment

# Load spectral library
library = MSExperiment()
MzMLFile().load('spectral_library.mzML', library)

# Match query spectra against library
# Returns similarity scores and library matches
```

## Related Skills

- data-import - Load raw MS data before identification
- protein-inference - Group peptides to proteins
- ptm-analysis - Identify modified peptides
