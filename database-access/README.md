# database-access

## Overview

Programmatic access to NCBI, EMBL-EBI (Ensembl, ENA, UniProt, BioStudies), and interaction-database resources. Covers Entrez E-utilities, the modern NCBI Datasets v2 CLI, Ensembl REST and BioMart, UniProt's post-2022 REST API, SRA/ENA/cloud-mirror sequencing data download, GEO expression data with SuperSeries handling, BLAST (remote and local with BLAST+ v5 databases), modern remote-homology methods (PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek), pre-computed ortholog retrieval, and protein-interaction networks across STRING, BioGRID, IntAct, SIGNOR, OmniPath.

**Tool type:** mixed | **Primary tools:** Bio.Entrez, NCBI Datasets CLI, BLAST+, sra-tools, Foldseek, MMseqs2, DIAMOND, HMMER, pybiomart, requests (UniProt REST, Ensembl REST, STRING, BioGRID, OmniPath, SIGNOR)

## Skills

| Skill | Description |
|-------|-------------|
| entrez-search | Search NCBI with ESearch/EInfo/EGQuery; history server, retmax caps, query translator, weekly index lag |
| entrez-fetch | Retrieve records with EFetch/ESummary; rettype matrix, GI deprecation, XML schema brittleness |
| entrez-link | Cross-database ELink; linkname semantics, asymmetric round-trip, neighbor_history for batches |
| batch-downloads | Bulk EFetch via history server; WebEnv TTL, EPost 200-ID limit, retry/resume; Datasets defection |
| ncbi-datasets-cli | NCBI Datasets v2 CLI (2023+); the modern bulk path for genome/gene data |
| sra-data | SRA + ENA + STRIDES cloud mirrors; fasterq-dump scratch trap, 10x technical reads |
| geo-data | GEO + BioStudies/ArrayExpress; SuperSeries trap, processed-vs-raw decision, GEOparse/GEOquery |
| blast-searches | Remote BLAST; Karlin-Altschul stats, max_target_seqs misuse (Shah 2019), CBS modes |
| local-blast | Local BLAST+ with v5 databases; -task taxonomy, soft masking, thread saturation |
| remote-homology | PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek (van Kempen 2024) for distant homology |
| ortholog-inference | Pre-computed orthologs from OrthoDB, Compara, OMA, eggNOG, PANTHER, KEGG, HomoloGene |
| uniprot-access | UniProt REST (post-2022 endpoint); JSON schema navigation, ID mapping async, stream vs search |
| ensembl-rest | Ensembl REST API; archive endpoints for reproducibility, VEP, Compara homology |
| biomart-queries | Ensembl BioMart for bulk ID mapping and coordinate / ortholog wide-tables (>5K rows) |
| interaction-databases | STRING, BioGRID, IntAct, SIGNOR, Reactome, HuRI, HuMAP, OmniPath; signed vs unsigned, license matrix |

## Example Prompts

- "Find PubMed articles about CRISPR published in 2024 with field-qualified query"
- "Download the human reference genome via NCBI Datasets CLI with auto MD5 verification"
- "Resolve GSE123456 to SRA runs via pysradb, then download FASTQ from the ENA mirror with md5 check"
- "Detect SuperSeries before processing a GEO accession; recommend per-SubSeries handling"
- "BLAST a protein with hitlist_size=500 to avoid the max_target_seqs trap; sort by bit-score not E-value"
- "Build a BLAST+ v5 database with parse_seqids; use dc-megablast for cross-species DNA"
- "Find structural homologs of a protein via Foldseek against AlphaFoldDB; use ProstT5 for sequence-only access"
- "Run MMseqs2 with -s 7.5 (HMMER-equivalent sensitivity) for distant protein homology"
- "Get high-confidence Compara orthologs of human BRCA1 in mouse, with confidence score"
- "Fetch UniProt P04637 as JSON; navigate the new (post-2022) nested schema for gene/sequence/PDB xrefs"
- "Map 5,000 Ensembl Gene IDs to HGNC + RefSeq + UniProt in one BioMart query (not REST loop)"
- "Pull STRING v12 interactions at confidence 700+ with per-channel breakdown (experiments vs textmining)"
- "Build a signed signaling network with SIGNOR (the only signed/mechanism-typed major DB)"
- "Aggregate STRING + OmniPath + BioGRID with per-edge provenance for cross-resource consensus"

## Requirements

```bash
# NCBI Entrez + UniProt + Ensembl REST (Python clients)
pip install biopython requests pandas pybiomart pysradb GEOparse networkx

# NCBI Datasets v2 CLI
conda install -c conda-forge ncbi-datasets-cli

# SRA toolkit
conda install -c bioconda sra-tools

# BLAST+
conda install -c bioconda blast

# Modern homology tools
conda install -c bioconda hmmer mmseqs2 diamond hhsuite foldseek

# Orthology
conda install -c bioconda orthofinder busco

# R alternatives (Bioconductor)
# BiocManager::install(c('GEOquery', 'biomaRt', 'STRINGdb', 'OmnipathR'))
```

## Related Skills

- **sequence-io** - Read/write sequence files after downloading
- **sequence-manipulation** - Work with downloaded sequences (slicing, translation, motif search)
- **alignment-files** - Process SAM/BAM after read alignment
- **variant-calling** - Variant calling and annotation from downloaded data
- **comparative-genomics** - De novo orthology computation (vs the pre-computed access here)
- **gene-regulatory-networks** - Network analysis from interaction data
- **data-visualization** - Visualize interaction networks (Cytoscape/networkx)
- **pathway-analysis** - GO/KEGG/Reactome enrichment using IDs mapped here
- **structural-biology** - PDB / AlphaFoldDB structures cross-referenced from UniProt
