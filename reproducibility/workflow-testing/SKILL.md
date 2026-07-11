---
name: bio-reproducibility-workflow-testing
description: Ships reference and synthetic test data plus expected outputs with a bioinformatics workflow so it can be regression-tested and a silent breakage caught automatically, the community's weakest reproducibility practice. Use when choosing a test tier (stub/dry-run wiring, tiny smoke data, full-size CI); generating shareable synthetic data with the flavour of real data without disclosing sensitive samples (dwgsim/ART/InSilicoSeq/Polyester simulators, or samtools/seqtk downsampling); deciding an assertion strategy for outputs that are not bytewise reproducible (normalize-then-compare VCF/BAM, snapshot testing, numeric tolerance, structural checks) instead of a naive md5 that fails on header timestamps and nondeterministic ordering; selecting a framework (nf-test, pytest-workflow, snakemake --generate-unit-tests, cwltest, miniwdl); and wiring CI (GitHub Actions matrix, scheduled runs) that alerts when a container/tool version change silently alters the expected output.
tool_type: mixed
primary_tool: nf-test
---

## Version Compatibility

Reference examples tested with: nf-test 0.9+, pytest-workflow 2.1+, Snakemake 8+, samtools 1.19+, dwgsim 0.1.13+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pytest --workflow-... --help`; `pip show pytest-workflow`

If code throws an error, introspect the installed tool and adapt the example to match the
actual API rather than retrying.

Note: nf-test config moved from `nextflow.config`-style to `nf-test.config` and the assertion DSL
(`assertAll`, `snapshot(...).match()`) evolves between releases; confirm against `nf-test --help` and
the current docs. `snakemake --generate-unit-tests` writes to `.tests/unit/`.

# Workflow Testing

**"Ship a sample input and its expected output so anyone can check the workflow still works"** -> Provide reference data (real, downsampled, or synthetic) and a machine-checkable expectation, so a user verifies the pipeline in minutes and CI flags the moment a tool change silently alters a result, instead of the failure being discovered by accident mid-experiment.

## The governing principle: testing is the weakest link, and "no error" is not "correct"

Providing a sample input and the corresponding expected output is not standard practice in bioinformatics, and it should be (Cohen-Boulakia 2017 *Future Gener Comput Syst* 75:284-298). A pipeline that exits 0 has proved only that no process crashed, not that the numbers are right. Two failure modes this catches that nothing else does:

- A dependency or container version bumps and the *output changes* while every process still exits 0. Without a regression test, this is discovered by accident, months later, and the reflex is to never reuse anything again.
- A refactor rewires a channel/rule and quietly drops samples. The run succeeds; the expectation (record count, per-sample outputs) fails.

Reference data and reference results should be part of the workflow artifact, versioned alongside the code, exactly as tests are in professional software engineering.

## Decision: the three test tiers

Run them in ascending cost; the cheap tiers gate the expensive ones in CI.

| Tier | Question it answers | Data | Runs the tools? | Mechanism |
|------|---------------------|------|-----------------|-----------|
| Stub / dry-run (wiring) | Is the DAG/channel graph valid; does every step connect? | none | no | `nextflow run -stub-run` / `-preview`; `snakemake -n` (dry run); `cwltool --validate` |
| Smoke (tiny data) | Do the real tools run end-to-end and produce the expected small output? | tiny fixture or simulated reads | yes, on toy input | nf-core `-profile test`; `nf-test`; `pytest-workflow`; snakemake `.tests` |
| Integration (full/near-full) | Does a realistic run still match the reference result? | downsampled or full reference dataset | yes | scheduled CI (nightly/weekly), larger runners |

Most day-to-day protection comes from the smoke tier: it is fast enough for every pull request and runs the actual tools. nf-core ships tiny inputs in the `nf-core/test-datasets` repo and wires them through `-profile test`; mirror that pattern for a bespoke pipeline.

## Decision: where the test data comes from

The keynote's key enabler is that you can now generate data "with the flavour of real data without disclosing anything sensitive." Prefer, in order: a tiny hand-curated fixture, a downsample of shareable real data, or a simulator seeded for determinism.

| Source | Tool | Use when | Watch |
|--------|------|----------|-------|
| Downsample real data | `samtools view -s SEED.FRAC`, `seqtk sample -s SEED` | shareable public data exists (e.g. a 1000G BAM) | the seed is part of the fixture; record it |
| Simulate short DNA reads | dwgsim, wgsim (BWA), ART | need reads with known truth (variants, coverage) and no sample disclosure | fix the RNG seed or the reads are not reproducible |
| Simulate metagenomic reads | InSilicoSeq (`iss generate`) | community profiles with known composition | model choice (HiSeq/NovaSeq) changes error profile |
| Simulate RNA-seq | Polyester (R) | known differential-expression truth | seed via `set.seed()` before `simulate_experiment` |
| Simulate long reads | pbsim, badread | ONT/PacBio error modelling | error model must match the chemistry claimed |
| Hand-curated fixture | manual | a 3-read, 2-variant edge case | keep it minimal and documented |

Whatever the source, the fixture is an input to reproducibility itself: pin the seed, keep it tiny (CI-sized), and version it with the tests.

## Decision: how to assert (the naive-md5 trap)

Genomic outputs are usually NOT bytewise reproducible: VCF/BAM headers carry timestamps, command lines, and program versions; record order can be nondeterministic across threads; floats vary across architectures. A plain `md5sum` of the whole file gives false failures. Choose the assertion to the output type.

| Output | Wrong assertion | Right assertion |
|--------|-----------------|-----------------|
| VCF | md5 of the file (header timestamp/`##bcftools_command` differ every run) | `bcftools view --no-version | grep -v '^##'` then compare records; or sort + normalize then md5 |
| BAM | md5 of the binary (compression/order differ) | compare `samtools flagstat` + sorted read records; or `samtools view` piped to a stable projection |
| Tables (TSV/CSV) | md5 when row order is nondeterministic | sort rows/columns first, then md5; or column-wise compare |
| Floating-point metrics | exact equality | tolerance compare (`abs(a-b) < eps`); nf-test `assert ... closeTo` |
| Large/opaque output | nothing | structural checks: non-empty, expected row count, schema/column names, key invariants |

Snapshot testing formalizes this: capture a known-good, normalized snapshot once, then diff every future run against it. When a container digest bumps and the biology changes, the snapshot diff is the alert. nf-test writes snapshots to a `.snap` file; treat a snapshot change in a PR as a claim that behavior changed on purpose, and review it.

## Decision: the test framework

| Framework | Engine | Test spec | Levels | Snapshot | Best when |
|-----------|--------|-----------|--------|----------|-----------|
| nf-test | Nextflow | Groovy `.nf.test` | process / workflow / pipeline | yes (`.snap`) | Nextflow/nf-core (the de facto standard) |
| pytest-workflow | any (Snakemake/CWL/shell) | YAML | pipeline (files/exit/stdout) | no (md5/contains) | engine-agnostic, declarative file checks |
| snakemake unit tests | Snakemake | generated Python | per-rule | no | `snakemake --generate-unit-tests` bootstraps from a real run |
| cwltest | CWL | YAML job + expected | tool / workflow | no | CWL conformance and regression |
| miniwdl / wdl-tests | WDL | Python / testkit | task / workflow | no | WDL, Terra ecosystem |

## nf-test: a process test with normalized comparison

```groovy
// tests/modules/mark_duplicates.nf.test -- runs the real tool on tiny simulated input,
// then asserts on NORMALIZED output so header timestamps never cause a false failure.
nextflow_process {
    name "Test MARK_DUPLICATES"
    process "MARK_DUPLICATES"

    test("marks duplicates on a tiny simulated BAM") {
        when {
            process {
                """
                input[0] = [ [ id:'test' ], file("\${projectDir}/tests/data/tiny.sorted.bam") ]
                """
            }
        }
        then {
            assertAll(
                { assert process.success },
                // stable projection: flagstat is order- and timestamp-independent
                { assert snapshot(
                    path(process.out.bam[0][1]).linesGzip.findAll { !it.startsWith('@PG') }
                  ).match("bam_body") },
                { assert path(process.out.metrics[0]).text.contains("PERCENT_DUPLICATION") }
            )
        }
    }
}
```

## pytest-workflow: engine-agnostic YAML

```yaml
# tests/test_pipeline.yml -- declarative; content checks avoid brittle whole-file md5.
- name: rnaseq smoke test on simulated reads
  command: nextflow run main.nf -profile test,docker --outdir out
  files:
    - path: out/salmon/test/quant.sf
      contains: ["NumReads"]           # schema/invariant, not a fragile md5
      must_not_contain: ["NaN"]
    - path: out/multiqc/multiqc_report.html
  exit_code: 0
```

## CI: fail fast, and alert on silent output drift

```yaml
# .github/workflows/ci.yml -- cheap tiers on every PR; full run on a schedule.
on: [pull_request]
jobs:
  wiring:                                   # tier 1: seconds, no tools
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: nextflow run main.nf -stub-run -profile test
  smoke:                                    # tier 2: minutes, real tools on tiny data
    needs: wiring
    runs-on: ubuntu-latest
    strategy:
      matrix: { profile: [docker, singularity] }   # a container-runtime change is a real risk
    steps:
      - uses: actions/checkout@v4
      - run: nf-test test --profile ${{ matrix.profile }}    # snapshot diff = the silent-drift alert
# tier 3 (full reference dataset) goes in a separate `schedule:` workflow, not per-PR.
```

The snapshot diff is the mechanism that "alerts automatically when a tool change silently breaks the expected output": pin the container by digest, and when someone bumps it, the smoke test's snapshot changes and the PR goes red instead of the surprise landing in someone's experiment.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| md5 test fails on every run, output looks identical | header timestamp / `##command` / program version in the file | strip headers or `--no-version`; compare records |
| Test passes locally, fails in CI (or vice versa) | nondeterministic record order across thread counts | sort before comparing; fix threads in the test profile |
| Simulated reads differ every run | RNG seed not fixed | set the simulator seed (`-z`, `set.seed()`, `-s`) and version the fixture |
| Snapshot changes and nobody knows why | a container/tool version drifted | pin by digest; review the snapshot diff as an intentional-change claim |
| Floating-point metric test is flaky | exact equality across architectures | tolerance compare (`closeTo`, `abs(a-b) < eps`) |
| "Green" pipeline still produced wrong biology | tested exit code only, not output content | add content/structural assertions, not just `exit_code: 0` |
| CI too slow to run per PR | full dataset in the smoke tier | move full runs to a scheduled job; keep PR data tiny |
| Refactor silently dropped samples, tests passed | no per-sample or record-count assertion | assert expected sample set / row count |

## Related Skills

- reproducibility/by-design - Why testing is pillar 3 and how it fits the reproducibility doctrine
- reproducibility/workflow-provenance - Provenance records which tool versions produced a tested output
- workflow-management/nf-core-pipelines - `-profile test` and `nf-core/test-datasets` conventions to mirror
- workflow-management/nextflow-pipelines - `-stub-run` wiring tests and pinning containers by digest
- read-qc/quality-reports - QC outputs (FastQC/MultiQC) that make good structural test targets

## References

- Cohen-Boulakia S, Belhajjame K, Collin O, et al. 2017. Scientific workflows for computational reproducibility in the life sciences: status, challenges and opportunities. *Future Gener Comput Syst* 75:284-298.
- Beaulieu-Jones BK, Greene CS. 2017. Reproducibility of computational workflows is automated using continuous analysis. *Nat Biotechnol* 35(4):342-346.
- Grüning B, Chilton J, Köster J, et al. 2018. Practical computational reproducibility in the life sciences. *Cell Syst* 6(6):631-635.
- Sandve GK, Nekrutenko A, Taylor J, Hovig E. 2013. Ten simple rules for reproducible computational research. *PLoS Comput Biol* 9(10):e1003285.
- Huang W, Li L, Myers JR, Marth GT. 2012. ART: a next-generation sequencing read simulator. *Bioinformatics* 28(4):593-594.
- Frazee AC, Jaffe AE, Langmead B, Leek JT. 2015. Polyester: simulating RNA-seq reads from differential expression experiments. *Bioinformatics* 31(17):2778-2784.
