---
name: hpc-slurm-nextflow
category: custom-zam-skills
primary_tool: slurm
description: Use when creating, reviewing, or troubleshooting SLURM, Nextflow, nf-core, Apptainer/Singularity, micromamba, and module-system workflows on university HPC clusters.
---

# HPC SLURM + Nextflow Skill

## Use when

- Creating `sbatch` scripts for bioinformatics tools.
- Running nf-core or custom Nextflow pipelines on SLURM.
- Debugging `.command.sh`, `.command.err`, `.command.log`, and work directories.
- Handling Apptainer/Singularity containers, Wave/Seqera containers, offline HPC environments, and module systems.
- Preparing commands for RAPTOR/Lanta-style project directories using `/project`, `/scratch`, and SLURM accounts.

## Core response pattern

For any HPC command/script request, provide:

1. **Assumptions**: cluster account, partition, threads, memory, input paths.
2. **Copy-paste script**: complete `sbatch` or command block.
3. **Safety checks**: `set -euo pipefail`, path checks, log directories.
4. **Expected output**: files created and where to inspect logs.
5. **Troubleshooting**: common failure modes and specific commands to inspect.

## SLURM script template

```bash
#!/usr/bin/env bash
#SBATCH -J JOB_NAME
#SBATCH -A ACCOUNT_ID
#SBATCH -p compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err

set -euo pipefail

mkdir -p logs results work

# Optional: module system
# module purge
# module load samtools/1.20

# Optional: micromamba
# source ~/.bashrc
# micromamba activate ENV_NAME

THREADS=${SLURM_CPUS_PER_TASK:-16}

INPUT="/path/to/input"
OUTDIR="results"

if [[ ! -e "$INPUT" ]]; then
    echo "ERROR: Input not found: $INPUT" >&2
    exit 1
fi

# Run command here
# tool --threads "$THREADS" --input "$INPUT" --outdir "$OUTDIR"
```

## Nextflow on SLURM

Prefer:

```bash
nextflow run PIPELINE \
  -profile singularity,slurm \
  -work-dir /scratch/$USER/nxf_work/PROJECT \
  --outdir /project/PROJECT/results \
  -resume
```

For local config snippets:

```groovy
process {
  executor = 'slurm'
  queue = 'compute'
  cpus = 8
  memory = '32 GB'
  time = '24h'
}

singularity {
  enabled = true
  autoMounts = true
  cacheDir = '/project/containers/singularity_cache'
}
```

## Debugging Nextflow

Always inspect the failing process directory:

```bash
cd /path/to/work/ab/cdef...
ls -lah
cat .command.sh
cat .command.err
cat .command.log
bash .command.run
```

Use these commands:

```bash
nextflow log
nextflow log <run_name> -f process,tag,status,exit,submit,duration,realtime,cpus,memory,workdir
nextflow run ... -resume
```

## Common issues

### `/etc/bashrc: BASHRCSOURCED: unbound variable`

Cause: `set -u` with a site bashrc expecting unset variables.

Fix options:

```bash
set +u
source /etc/bashrc
set -u
```

or avoid sourcing global bashrc inside strict-mode scripts.

### Internet only on head node

- Pull containers/databases on the head node first.
- Set Apptainer/Singularity cache to a shared project directory.
- Avoid database download inside compute jobs.

### Apptainer bind issues

Use explicit binds:

```bash
apptainer exec \
  -B /project:/project \
  -B /scratch:/scratch \
  container.sif command
```

## Version Compatibility

Reference versions:

- SLURM 22+
- Nextflow 23.10+
- Apptainer 1.2+
- SingularityCE 3.11+
- micromamba 1.5+

Verify syntax and cluster-specific options if versions or local policies differ.
