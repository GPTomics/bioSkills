# ncRNA Search - Usage Guide

## Overview
Search for non-coding RNA homologs and classify RNA families with Infernal covariance models (CMs) against Rfam. A CM scores both sequence and consensus secondary-structure conservation jointly, giving higher sensitivity than BLAST or profile HMMs for STRUCTURED RNAs whose sequence has diverged but whose structure is conserved. For RNAs without conserved structure (lncRNAs, mature miRNAs), a CM offers no advantage and a faster sequence method is correct.

## Prerequisites
```bash
# Infernal
conda install -c bioconda infernal

# Rfam database (pre-calibrated: press it, never recalibrate it)
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz && cmpress Rfam.cm
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin

# Family-specialized tools (often the better choice for a known class)
conda install -c bioconda trnascan-se barrnap

# Python dependencies
pip install biopython pandas
```

## Quick Start
Tell your AI agent what you want to do:
- "Scan my transcripts against Rfam to classify ncRNA families"
- "Find rRNAs, tRNAs, and snoRNAs in my genome"
- "Build a covariance model for my novel RNA family"
- "Should I use a covariance model or BLAST for this RNA?"
- "Find all tRNAs in my genome assembly"

## Example Prompts

### Rfam Classification
> "I have non-coding transcripts. Scan them against Rfam and classify by family, resolving clan overlaps."

> "Search this single RNA against Rfam and tell me what it is and whether the structure matches."

### Genome-wide ncRNA Search
> "Annotate ncRNA families in my bacterial genome using the documented Rfam pipeline."

> "Find all tRNAs in my genome assembly with the right tool."

### Tool Choice
> "My RNA is a lncRNA with no obvious conserved structure. Is a covariance model worth it?"

### Custom Covariance Models
> "I have a Stockholm alignment of a novel RNA family with a consensus structure. Validate the covariation, then build and search a CM."

### Parsing
> "Parse the Infernal --fmt 2 output, deoverlap clans, and summarize family assignments."

## What the Agent Will Do
1. Decide whether a covariance model is the right tool (structured ncRNA) or a sequence method/specialist is better
2. Run cmscan (sequences vs Rfam) or cmsearch (one CM vs database) with gathering thresholds and clan resolution
3. Deoverlap clan hits and parse the --fmt 2 tabular output with the correct column layout
4. For custom families, validate the consensus structure (R-scape), then build, calibrate, and search
5. Summarize family assignments and flag truncated or pseudogene-like hits

## Tips
- **CM vs sequence search** - Covariance models win only for structured ncRNAs with diverged sequence. For lncRNAs, mature miRNAs, and sequence-only motifs use nhmmer/BLASTN; a CM is slow overhead with no gain.
- **Gathering thresholds** - Use `--cut_ga` for Rfam; the curated per-family bit-score cutoffs are database-size-independent. A flat `-E 1e-5` reintroduces the DB-size dependence GA was built to avoid.
- **E-values depend on database size** - The same hit gets a different E-value depending on what was searched (Z); fix it with `-Z <Mb>` for reproducibility, and prefer bit-score thresholds for comparison.
- **Clan overlap** - Use `--fmt 2 --clanin` then drop hits marked `=` (`grep -v ' = '`); this is the documented Rfam deoverlap. `--oclan` is a valid alternative but not the documented path.
- **Calibration** - Run `cmcalibrate` on any custom CM before trusting E-values (or threshold on bit score). Rfam.cm is already calibrated; only `cmpress` it.
- **Specialized tools** - Prefer tRNAscan-SE 2.0 for tRNA, barrnap/RNAmmer for rRNA, snoscan/snoReport for snoRNA. They add biology a generic scan lacks.
- **A hit is a hypothesis** - A CM hit is a family assignment plus a structural hypothesis, not proof of transcription or function; pseudogenes score well.
- **SS_cons is what makes it a CM** - A custom model needs a real `#=GC SS_cons`; a structure-free alignment yields only an HMM, and the structure should pass R-scape covariation before you build.

## Related Skills
- secondary-structure-prediction - Predict structure for novel ncRNA candidates with no Rfam hit
- covariation-analysis - Validate a custom CM's consensus structure with R-scape before building
- structure-probing - Experimental reactivities to corroborate a CM's consensus structure
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation pipelines
- alignment/msa-statistics - Evaluate alignment quality before CM building
- database-access/entrez-fetch - Fetch Rfam/RNAcentral records
