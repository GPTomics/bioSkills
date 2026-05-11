# LANTA HPC Profile

Source of truth: ThaiSC LANTA user guide, https://thaisc.atlassian.net/wiki/spaces/LANTA/overview

Use this profile only when the user explicitly targets LANTA. Keep LANTA commands separate from RAPTOR/MedCMU-HPC commands.

## Identity and Access

- Login/submission host: `lanta.nstda.or.th`.
- Data transfer host: `transfer.lanta.nstda.or.th`.
- LANTA uses two-factor authentication for SSH login.
- Use the transfer node for large transfers, downloads, container pulls, and environment preparation that may take more than about 5 minutes.
- Do not run heavy bioinformatics computation on the frontend/login node.

```bash
ssh <LANTA_USER>@lanta.nstda.or.th
scp -r <LOCAL_DIR> <LANTA_USER>@transfer.lanta.nstda.or.th:/project/<PROJECT_ID-SHORTNAME>/
rsync -rvz <LOCAL_DIR>/ <LANTA_USER>@transfer.lanta.nstda.or.th:/project/<PROJECT_ID-SHORTNAME>/<DEST>/
```

Prefer `rsync -rvz` for LANTA project directories. Avoid `rsync -a` unless the user confirms permissions and group ownership are appropriate for the project directory.

## Storage Model

Use the documented LANTA storage tiers:

- `$HOME`: personal storage, small files only, documented quota is 100 GB.
- `/project/<PROJECT_ID-SHORTNAME>`: shared project space for data, results, and runnable project files; documented default project quota is 5 TB.
- `/scratch/<PROJECT_ID-SHORTNAME>`: shared temporary high-performance space for calculation intermediates; documented as temporary and subject to deletion after 30 days since last access.

Recommended layout:

```text
/project/<PROJECT_ID-SHORTNAME>/
  data/
  refs/
  results/
  logs/
  containers/

/scratch/<PROJECT_ID-SHORTNAME>/<USER>/
  nextflow-work/
  tmp/
  apptainer-cache/
```

Never place large Nextflow `work/` directories, FASTQ/BAM files, or container caches in `$HOME`.

Check storage before large work:

```bash
myquota
```

## Slurm Defaults

LANTA examples use the `compute` partition and an account ID such as `lt123456`; keep the public skill masked with `<LANTA_ACCOUNT>`.

Minimal Slurm header pattern:

```bash
#!/usr/bin/env bash
#SBATCH -J <JOB_NAME>
#SBATCH -p compute
#SBATCH -A <LANTA_ACCOUNT>
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o logs/%x.%j.out
#SBATCH -e logs/%x.%j.err
```

Useful commands:

```bash
sbatch submit.sh
myqueue
cat slurm-<JOB_ID>.out
```

Use Slurm and Nextflow resource directives instead of running compute-heavy commands directly on the login host.

## Module System

LANTA uses Lmod environment modules.

```bash
module avail      # or: ml av
module load <module>
module list       # or: ml
module swap <old> <new>
module unload <module>
module restore
module purge
```

Verify exact module names on LANTA before production runs:

```bash
module avail nextflow
module avail apptainer
module avail singularity
module avail java
```

## Containers and Offline Compute Nodes

Treat LANTA compute/GPU/memory nodes as environments without direct internet access unless the user confirms otherwise. Pull containers and download reference assets on `transfer.lanta.nstda.or.th`, then reuse the shared cache/path from compute jobs.

Recommended cache variables:

```bash
export APPTAINER_CACHEDIR=/scratch/<PROJECT_ID-SHORTNAME>/<USER>/apptainer-cache
export APPTAINER_TMPDIR=/scratch/<PROJECT_ID-SHORTNAME>/<USER>/apptainer-tmp
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
```

For explicit binds:

```bash
apptainer exec \
  -B /project:/project \
  -B /scratch:/scratch \
  <container.sif> <command>
```

## Nextflow Pattern

Use a site config such as `examples/nextflow.lanta.config` and submit from `lanta.nstda.or.th`:

```bash
nextflow run nf-core/rnaseq \
  -r <NFCORE_RNASEQ_VERSION> \
  -c nextflow.lanta.config \
  --input rnaseq_samplesheet.csv \
  --outdir /project/<PROJECT_ID-SHORTNAME>/results/rnaseq \
  --fasta /project/<PROJECT_ID-SHORTNAME>/refs/genome.fa \
  --gtf /project/<PROJECT_ID-SHORTNAME>/refs/genes.gtf \
  -work-dir /scratch/<PROJECT_ID-SHORTNAME>/<USER>/nextflow-work/rnaseq \
  -resume \
  -with-report /project/<PROJECT_ID-SHORTNAME>/results/rnaseq/report.html \
  -with-trace /project/<PROJECT_ID-SHORTNAME>/results/rnaseq/trace.txt \
  -with-timeline /project/<PROJECT_ID-SHORTNAME>/results/rnaseq/timeline.html \
  -with-dag /project/<PROJECT_ID-SHORTNAME>/results/rnaseq/flow.svg
```

## LANTA Troubleshooting

1. Confirm the command uses LANTA paths, not RAPTOR `/common/sif` or MedCMU examples.
2. Run `myquota` before rerunning failed large jobs.
3. Check queue state with `myqueue`; use Slurm tools if available.
4. Inspect `.nextflow.log` in the launch directory.
5. Inspect the failed work directory under `/scratch/<PROJECT_ID-SHORTNAME>/<USER>/nextflow-work/`.
6. Check `.command.err`, `.command.out`, `.command.sh`, `.command.run`, and `.exitcode`.
7. If container or reference downloads failed on compute nodes, pre-stage them via `transfer.lanta.nstda.or.th` and rerun with `-resume`.
