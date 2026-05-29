# Contamination Detection - Usage Guide

## Overview

Contamination detection covers two unrelated problems that share a word. The first is **foreign sequence in a single-organism assembly** - bacterial, human, or vector contigs in a eukaryotic nuclear genome or a cultured isolate - answered by NCBI FCS-adaptor, FCS-GX (now mandatory for GenBank submission), and the BlobToolKit blob plot, and resolved as a set of contigs to remove, trim, or split. The second is **MAG/bin quality** - how complete and how chimeric a metagenome-assembled genome is - answered as a completeness/contamination percentage pair plus a chimerism verdict from CheckM2 and GUNC run together, with GTDB-Tk taxonomy, judged against the MIMAG standard. Crossing the two toolsets (CheckM2 on a eukaryote, FCS-GX percentage as MIMAG contamination) produces confidently wrong numbers. The deepest judgment is that contamination and real horizontal gene transfer are the same signal read two ways - the tell is physical integration, not taxonomy.

## Prerequisites

```bash
# FCS-adaptor + FCS-GX (Docker/Singularity wrappers; large GX DB)
curl -LO https://github.com/ncbi/fcs/raw/main/dist/fcs.py
curl -LO https://github.com/ncbi/fcs/raw/main/dist/run_fcsadaptor.sh
# Download the GX database (~470 GiB; needs ~512 GiB RAM or a cloud VM)
python3 fcs.py db get --mft "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest" --dir "$GXDB_LOC/gxdb"

# CheckM2
pip install checkm2
checkm2 database --download

# GUNC
conda install -c bioconda gunc
gunc download_db ./gunc_db

# GTDB-Tk (reference package release must match the binary)
conda install -c bioconda gtdbtk
download-db.sh

# BlobToolKit + tiara
pip install blobtoolkit tiara
```

## Quick Start

Tell your AI agent what you want to do:
- "Screen my insect assembly for foreign contamination before GenBank submission"
- "Run CheckM2 and GUNC on my MAGs and tell me which pass MIMAG high-quality"
- "Make a blob plot of my assembly and tell me which contigs are contaminants"
- "This bin is 95% complete and 3% contaminated - is it actually a chimera?"
- "Is this bacterial-looking contig real HGT or contamination?"

## Example Prompts

### Foreign-sequence screening (single-organism assembly)
> "I have a draft fungal nuclear assembly going to GenBank. Run FCS-adaptor then FCS-GX with the correct source tax-id, and walk me through the action report - which contigs are EXCLUDE/TRIM/FIX versus REVIEW, and which REVIEW items I should inspect rather than auto-clean."
> "Build a BlobToolKit BlobDir for my assembly from a DIAMOND hit file and a coverage BAM, render the blob plot, and tell me whether the off-cloud contigs are contaminants based on GC, coverage, and taxonomy together - not color alone."
> "FCS-GX flagged a bacterial-looking gene island in my aphid genome as EXCLUDE. It sits on a host scaffold with host GC and coverage. Should I keep it as possible Wolbachia-derived HGT?"

### MAG/bin quality
> "Run CheckM2 and GUNC on my 200 metagenome bins, merge the reports, and classify each into MIMAG high/medium/low quality - flag any bin that passes CheckM2 contamination but fails GUNC as a likely chimera."
> "One of my MAGs is GUNC-pass but the reference_representation_score is 0.3. Can I trust the pass, or is the lineage too novel to judge?"
> "CheckM reports 18% contamination but high strain heterogeneity for this bin. Is that contamination or co-binned strains?"

### Organelle and category-error triage
> "My nuclear assembly has a few contigs with 50x the average coverage that hit mitochondrial genes. Are these the organelle genome, contamination, or NUMTs, and what should I do with each?"
> "I ran CheckM2 on my plant nuclear genome and got 12% completeness. Why is that meaningless, and what should I run instead?"

## What the Agent Will Do

1. Determine which problem applies - foreign sequence in a single-organism assembly, or MAG/bin quality - before choosing tools.
2. For foreign-sequence screening: run FCS-adaptor first (cheap), then FCS-GX against the GX database with the source tax-id, and parse the action report into EXCLUDE/TRIM/FIX/REVIEW tiers.
3. Build a BlobToolKit blob plot (GC x coverage x taxonomy) for visual triage when contigs need investigation.
4. For MAGs: run CheckM2 and GUNC together, merge the reports, add GTDB-Tk taxonomy, and apply MIMAG thresholds plus the GUNC pass flag (gated on RRS).
5. Triage organellar contigs into a third category (separate, not delete) and distinguish free organelle from NUMTs by coverage.
6. Flag the HGT-vs-contamination decision by physical integration rather than taxonomy, and never report a contamination % without its GUNC verdict.

## Tips

- Run FCS-adaptor before FCS-GX: adaptor/vector hits are unambiguous and the adaptor DB is tiny.
- FCS-GX needs the ~470 GiB GX database resident in RAM (~512 GiB host); underprovisioned it crawls rather than fails. Plan a high-memory or cloud machine.
- `fcs.py clean genome` applies EXCLUDE/TRIM/FIX only - it does NOT act on REVIEW items, which are exactly where real HGT/symbiont biology hides.
- Always report CheckM2 and GUNC together. A 3%-contamination MAG can still be a disjoint-marker chimera that only GUNC catches.
- A GUNC pass at low reference_representation_score (< 0.5) means "can't tell", not "clean".
- The blob-plot decision is cluster geometry (off the main GC + coverage cloud), never taxonomy color alone; no-hit is not a contaminant.
- Organelles are a third category: separate and submit as their own record; do not delete, and do not strip genuine NUMTs (coverage discriminates).
- Never call CheckM2/CheckM on a eukaryotic nuclear assembly - use BUSCO. Never read an FCS-GX flagged-length fraction as the MIMAG contamination %.

## Related Skills

- assembly-qc - BUSCO completeness for eukaryotes (not CheckM2) and the QC handoff
- metagenome-assembly - Binning produces the bins this skill scores with CheckM2 + GUNC
- hifi-assembly - False duplications inflate apparent content; distinct from contamination
- long-read-assembly - Produces the contigs screened here for foreign sequence
- metagenomics/kraken-classification - Read/contig taxonomic classification and pre-assembly host screening
- comparative-genomics/genome-distance-and-species-delineation - ANI/species placement after decontamination
- workflows/genome-assembly-pipeline - End-to-end assemble -> QC -> decontaminate
