#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: RNA-seq quantification workflow
doc: |
  Trims reads with fastp and quantifies with Salmon.
  Typical runtime: ~30 min per sample at 8 threads.

requirements:
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  sample_ids:
    type: string[]
    doc: Sample identifiers matching FASTQ file prefixes
  fastq_1_files:
    type: File[]
    doc: Array of R1 FASTQ files (gzipped)
  fastq_2_files:
    type: File[]
    doc: Array of R2 FASTQ files (gzipped)
  salmon_index:
    type: Directory
    doc: Pre-built Salmon index directory
  threads:
    type: int
    default: 8
    doc: Threads per sample (8-16 typical for Salmon)

outputs:
  trimmed_reads_1:
    type: File[]
    outputSource: fastp/trimmed_1
  trimmed_reads_2:
    type: File[]
    outputSource: fastp/trimmed_2
  fastp_reports:
    type: File[]
    outputSource: fastp/json_report
  quant_dirs:
    type: Directory[]
    outputSource: salmon_quant/quant_dir

steps:
  fastp:
    run: fastp.cwl
    scatter: [reads_1, reads_2, sample_id]
    scatterMethod: dotproduct
    in:
      reads_1: fastq_1_files
      reads_2: fastq_2_files
      sample_id: sample_ids
      threads: threads
    out: [trimmed_1, trimmed_2, json_report]

  salmon_quant:
    run: salmon_quant.cwl
    scatter: [reads_1, reads_2, sample_id]
    scatterMethod: dotproduct
    in:
      reads_1: fastp/trimmed_1
      reads_2: fastp/trimmed_2
      sample_id: sample_ids
      index: salmon_index
      threads: threads
    out: [quant_dir]

---

# fastp.cwl - inline tool definition
cwlVersion: v1.2
class: CommandLineTool
id: fastp
baseCommand: fastp

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/fastp:0.23.4--hadf994f_2
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 4000

inputs:
  reads_1:
    type: File
    inputBinding:
      prefix: -i
  reads_2:
    type: File
    inputBinding:
      prefix: -I
  sample_id:
    type: string
  threads:
    type: int
    default: 4
    inputBinding:
      prefix: --thread

arguments:
  - prefix: -o
    valueFrom: $(inputs.sample_id)_trimmed_R1.fq.gz
  - prefix: -O
    valueFrom: $(inputs.sample_id)_trimmed_R2.fq.gz
  - prefix: --json
    valueFrom: $(inputs.sample_id)_fastp.json

outputs:
  trimmed_1:
    type: File
    outputBinding:
      glob: "*_trimmed_R1.fq.gz"
  trimmed_2:
    type: File
    outputBinding:
      glob: "*_trimmed_R2.fq.gz"
  json_report:
    type: File
    outputBinding:
      glob: "*_fastp.json"

---

# salmon_quant.cwl - inline tool definition
cwlVersion: v1.2
class: CommandLineTool
id: salmon_quant
baseCommand: [salmon, quant]

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0
  ResourceRequirement:
    coresMin: $(inputs.threads)
    # Salmon typically needs 8-16GB depending on index size
    ramMin: 16000

inputs:
  index:
    type: Directory
    inputBinding:
      prefix: -i
  reads_1:
    type: File
    inputBinding:
      prefix: "-1"
  reads_2:
    type: File
    inputBinding:
      prefix: "-2"
  sample_id:
    type: string
  threads:
    type: int
    default: 8
    inputBinding:
      prefix: --threads

arguments:
  - prefix: -l
    valueFrom: A
  - prefix: -o
    valueFrom: $(inputs.sample_id)_salmon
  - --validateMappings

outputs:
  quant_dir:
    type: Directory
    outputBinding:
      glob: "*_salmon"
