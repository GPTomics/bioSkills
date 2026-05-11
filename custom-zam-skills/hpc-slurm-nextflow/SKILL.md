---
name: hpc-slurm-nextflow
category: custom-zam-skills
tool_type: cli
primary_tool: Nextflow
description: Use when creating, reviewing, or troubleshooting Slurm, Nextflow, nf-core, Apptainer/Singularity, module-system, and RNA-seq workflows on LANTA or RAPTOR/MedCMU-HPC clusters. Use this for cluster-specific commands, storage paths, job submission, container cache setup, nf-core/rnaseq launches, and debugging failed HPC or Nextflow jobs.
---

# HPC Slurm + Nextflow

## Cluster Routing

Identify the target HPC before giving site-specific commands.

- If the user says LANTA, use `references/lanta.md` plus shared workflow guidance in `references/nfcore-rnaseq.md`.
- If the user says RAPTOR, MedCMU, or MedCMU-HPC, use `references/raptor.md` plus shared workflow guidance in `references/nfcore-rnaseq.md`.
- If the user does not identify the cluster, ask whether the target is LANTA or RAPTOR before using paths, partitions, accounts, module names, or container cache locations.
- Never mix site-specific conventions. Do not use RAPTOR `/common/sif` paths in LANTA commands. Do not use LANTA `myquota`, `myqueue`, `transfer.lanta.nstda.or.th`, or `#SBATCH -A <LANTA_ACCOUNT>` in RAPTOR commands.

## Default Response Pattern

For command, config, or troubleshooting requests, provide:

1. Assumptions: cluster, project path, scratch path, account/partition placeholder, module/container strategy, input type, and output directory.
2. Runnable block: a complete `nextflow run`, `sbatch`, or diagnostic command sequence.
3. Safety checks: avoid compute on login nodes, keep large intermediates out of `$HOME`, create log/work/cache directories, and keep credentials/account IDs masked.
4. Expected outputs: final result directories, reports, Nextflow trace/report/timeline/DAG files, and where logs appear.
5. Troubleshooting: commands for `.nextflow.log`, failed work directories, `.command.*` files, Slurm accounting, quota, and queue state.

## Reusable Resources

Load only the files needed for the cluster and task:

- `references/lanta.md`: LANTA login, storage, transfer node, Slurm account/partition style, Lmod, no-compute-internet behavior, quota, scratch policy, and troubleshooting.
- `references/raptor.md`: RAPTOR/MedCMU-HPC SSH, Slurm, modules, Apptainer, `/common/sif`, installed bioinformatics modules, job policy, and troubleshooting.
- `references/nfcore-rnaseq.md`: shared nf-core/rnaseq samplesheet, parameters, command pattern, output expectations, and cluster-specific launch notes.

Use these examples when the user asks for concrete files or copy-paste commands:

- `examples/nextflow.lanta.config`
- `examples/nextflow.raptor.config`
- `examples/run_nfcore_rnaseq_lanta.sh`
- `examples/run_nfcore_rnaseq_raptor.sh`
- `examples/rnaseq_samplesheet.csv`

## Shared Nextflow Standards

Prefer a Nextflow Slurm profile over wrapping all pipeline work inside one large `sbatch` allocation. Use process labels and per-process resources, then tune from `trace.txt`.

Include production run metadata when practical:

```bash
-resume \
-with-report report.html \
-with-trace trace.txt \
-with-timeline timeline.html \
-with-dag flow.svg
```

For nf-core/rnaseq, prefer explicit references and a pinned pipeline release when the user knows the reference build. Do not invent genome paths or account IDs. Keep placeholders such as `<PROJECT>`, `<LANTA_ACCOUNT>`, `<RAPTOR_HOST>`, and `<REFERENCE_DIR>` until the user supplies concrete values.

## Troubleshooting Order

1. Confirm the cluster profile matches the command path and scheduler style.
2. Check quota and storage location before debugging tool behavior.
3. Inspect `.nextflow.log` for orchestration and config errors.
4. Inspect the failed work directory: `.command.sh`, `.command.run`, `.command.err`, `.command.out`, and `.exitcode`.
5. Check Slurm state with site-appropriate commands: LANTA `myqueue`/Slurm tools, RAPTOR `squeue`, `sacct`, `scancel`.
6. Check container cache, bind paths, and offline compute-node assumptions.
7. Use `-resume` after fixing inputs, config, resources, or cache paths.

## Guardrails

- Use masked practical defaults in public examples. Do not commit real private project names, account IDs, sample paths, or credentials.
- Do not recommend running heavy computation on login/frontend nodes.
- Do not store Nextflow `work/`, large FASTQ/BAM files, or container caches in `$HOME`.
- Verify current module names and pipeline versions on the cluster before presenting final production commands.
- Separate observation from certainty when a site policy is not documented in the loaded reference.
