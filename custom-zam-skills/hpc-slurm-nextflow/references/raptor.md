# RAPTOR / MedCMU-HPC Profile

Source of truth: MedCMU-HPC documentation, https://medcmu-hpc.netlify.app/docs/intro/

Use this profile only when the user explicitly targets RAPTOR, MedCMU, or MedCMU-HPC. Keep RAPTOR commands separate from LANTA commands.

## Identity and Access

- Login pattern: `ssh <username>@<host>`.
- The user must provide the actual RAPTOR/MedCMU host if it is not already known.
- Windows users can use MobaXTerm; Linux/macOS users can use a normal terminal.
- SSH key login is recommended by the MedCMU documentation.

```bash
ssh <RAPTOR_USER>@<RAPTOR_HOST>
scp <LOCAL_FILE> <RAPTOR_USER>@<RAPTOR_HOST>:<REMOTE_DIR>/
rsync -rvz <LOCAL_DIR>/ <RAPTOR_USER>@<RAPTOR_HOST>:<REMOTE_DIR>/
```

## Policy and Job Conduct

- Login nodes are for job management, file checks, and lightweight setup, not computation.
- Submit compute work through Slurm using `sbatch` or `srun`.
- Keep running jobs manageable; the MedCMU guidance says ideally under 4 running jobs per user.
- Use GPU/specialized resources only when the workflow actually needs them.
- Project leaders manage project directories; keep only research data on the HPC.

## Slurm Defaults

MedCMU-HPC Slurm examples use `short` in sample commands. Keep the partition configurable unless the user gives a project-specific default.

Minimal Slurm header pattern:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --partition=<PARTITION>
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
```

Useful commands:

```bash
sinfo
sbatch script.sbatch
squeue
sacct
scancel <JOB_ID>
srun -p short -t 01:00:00 --mem=4G -c 1 --pty bash
```

Common resource flags:

```text
--partition / -p
--time / -t
--mem
--mem-per-cpu
--cpus-per-task / -c
--nodes / -N
--gpus / -G
--job-name / -J
--output / -o
--error / -e
```

## Module System

RAPTOR/MedCMU-HPC uses environment modules.

```bash
module avail        # or: ml av
module load <name>  # or: ml load <name>
module list         # or: ml
module swap <old> <new>
module unload <name>
module restore
module purge
ml spider <name>
```

Documented useful modules/environments include:

- `apptainer/1.3.5`
- `nextflow/24.10.6`, `nextflow/25.04.6`, `nextflow/25.10.3`, `nextflow/24.10.0`
- `nf-core` as a Micromamba environment
- `fastp/1.0.0`, `fastp/0.23.4`
- `fastqc/0.12.1`
- `multiqc/1.25.1`
- `star/2.7.11b` and `STAR/2.7.11b-GCC-13.2.0`
- `samtools/1.21`
- `subread/2.1.1`
- `java-jdk/21.0.3+9`

Verify current module names on the cluster before production runs:

```bash
module avail nextflow
module avail apptainer
module avail fastp
module avail fastqc
module avail multiqc
module avail star
module avail subread
```

If using Micromamba environments:

```bash
module load mamba
micromamba env list
micromamba activate <environment_name>
```

## Apptainer and SIF Images

MedCMU-HPC documents a central SIF image collection under `/common/sif/`.

```bash
module load apptainer
ls /common/sif/*/*.sif
apptainer pull image_name.sif docker://repository_name:tag
apptainer shell image_name.sif
apptainer exec image_name.sif <command>
apptainer exec --nv image_name.sif <gpu_command>
```

Use `--nv` only for GPU workflows. Do not assume `/common/sif/` exists on LANTA.

## Nextflow Pattern

RAPTOR examples can load `apptainer` and `nextflow` modules before running a workflow.

```bash
module purge
module load apptainer
module load nextflow

nextflow run nf-core/rnaseq \
  -r <NFCORE_RNASEQ_VERSION> \
  -c nextflow.raptor.config \
  --input rnaseq_samplesheet.csv \
  --outdir <PROJECT_DIR>/results/rnaseq \
  --fasta <REFERENCE_DIR>/genome.fa \
  --gtf <REFERENCE_DIR>/genes.gtf \
  -work-dir <SCRATCH_DIR>/nextflow-work/rnaseq \
  -resume
```

## RAPTOR Troubleshooting

1. Confirm the command uses RAPTOR/MedCMU paths and module names, not LANTA `myqueue`, `myquota`, `transfer.lanta.nstda.or.th`, or `#SBATCH -A <LANTA_ACCOUNT>`.
2. Confirm compute work is submitted through Slurm, not run directly on a login node.
3. Use `squeue`, `sacct`, and logs under `logs/` to inspect scheduler state.
4. Inspect `.nextflow.log` in the launch directory.
5. Inspect failed work directories and `.command.*` files.
6. Check `module list`, `which nextflow`, `nextflow -version`, and `apptainer --version`.
7. If using a central SIF path, confirm it exists with `ls /common/sif/*/*.sif` or use a project-local SIF/cache path.
8. Rerun with `-resume` after fixing inputs, resources, modules, or containers.
