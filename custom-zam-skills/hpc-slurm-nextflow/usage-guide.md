# HPC Slurm + Nextflow Usage Guide

This skill supports two separate HPC profiles. Select the cluster first, then load the matching reference.

## Cluster Selection

- LANTA: use `references/lanta.md`.
- RAPTOR / MedCMU-HPC: use `references/raptor.md`.
- nf-core/rnaseq shared flow: use `references/nfcore-rnaseq.md` after selecting the cluster.

If the user has not named the cluster, ask whether they mean LANTA or RAPTOR before giving commands with site-specific paths, partitions, account flags, modules, or container locations.

## Do Not Mix Profiles

- LANTA-only conventions: `lanta.nstda.or.th`, `transfer.lanta.nstda.or.th`, `/project/<PROJECT_ID-SHORTNAME>`, `/scratch/<PROJECT_ID-SHORTNAME>`, `myquota`, `myqueue`, `#SBATCH -A <LANTA_ACCOUNT>`, `#SBATCH -p compute`.
- RAPTOR-only conventions: `ssh <username>@<host>`, standard `squeue`/`sacct` flow, MedCMU module inventory, `/common/sif/`, `module load apptainer`, `module load nextflow`.

## Examples

- `examples/nextflow.lanta.config`: LANTA Slurm + Apptainer config.
- `examples/nextflow.raptor.config`: RAPTOR Slurm + Apptainer config.
- `examples/run_nfcore_rnaseq_lanta.sh`: copy-paste LANTA launcher.
- `examples/run_nfcore_rnaseq_raptor.sh`: copy-paste RAPTOR launcher.
- `examples/rnaseq_samplesheet.csv`: nf-core/rnaseq samplesheet skeleton.

## Default Safety Rules

- Keep account IDs, usernames, hostnames, and project names masked until the user provides them.
- Keep large work directories and containers outside `$HOME`.
- Use Slurm for compute work; do not run heavy analysis directly on login nodes.
- Include `-resume`, report, trace, timeline, and DAG outputs for production Nextflow runs.
- Verify module names and pipeline versions on the target cluster before production execution.
