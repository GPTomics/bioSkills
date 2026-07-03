#!/usr/bin/env cwl-runner
# Validate this document (no compute) with: cwltool --validate rnaseq.cwl
# Self-contained via $graph: one Workflow (id: main) plus its CommandLineTools in a single file.
# Run the workflow entry point with: cwltool rnaseq.cwl#main job.yml
cwlVersion: v1.2

$graph:

  - id: main
    class: Workflow
    label: RNA-seq trim / align / index / archive workflow
    doc: |
      Per sample: fastp trims -> STAR aligns to a prebuilt genomeDir -> samtools indexes the
      coordinate-sorted BAM -> samtools converts to a reference-compressed CRAM. Demonstrates
      secondaryFiles on both an indexed reference File input and a BAM output. Digest-pinned
      images keep the run reproducible; :latest would silently drift and break reproducibility.

    requirements:
      ScatterFeatureRequirement: {}

    inputs:
      sample_ids:
        type: string[]
        doc: Sample identifiers matching the FASTQ arrays position-for-position
      fastq_1_files:
        type: File[]
        doc: R1 FASTQ files (gzipped)
      fastq_2_files:
        type: File[]
        doc: R2 FASTQ files (gzipped)
      genome_dir:
        type: Directory
        doc: Prebuilt STAR genomeDir (a Directory, not an indexed File)
      reference:
        type: File
        doc: Genome FASTA; its .fai and .dict companions travel with it as a type property
        secondaryFiles:
          # samtools/CRAM needs the .fai; Picard/GATK tools need the ^.dict. The caret ^ STRIPS one
          # extension: genome.fasta -> genome.dict (NOT genome.fasta.dict). required:false = optional.
          - {pattern: .fai, required: true}
          - {pattern: ^.dict, required: false}
      threads:
        type: int
        default: 8
        doc: Threads per sample; 8-16 typical for STAR/fastp on a human genome

    outputs:
      fastp_reports:
        type: File[]
        outputSource: fastp/json_report
      indexed_bams:
        type: File[]
        outputSource: samtools_index/indexed_bam   # each carries its .bai as a secondaryFile
      cram_files:
        type: File[]
        outputSource: samtools_cram/cram

    steps:
      fastp:
        run: '#fastp'
        scatter: [reads_1, reads_2, sample_id]
        scatterMethod: dotproduct                  # equal-length arrays zipped 1:1 (R1[i] with R2[i])
        in:
          reads_1: fastq_1_files
          reads_2: fastq_2_files
          sample_id: sample_ids
          threads: threads
        out: [trimmed_1, trimmed_2, json_report]

      star_align:
        run: '#star_align'
        scatter: [reads_1, reads_2, sample_id]
        scatterMethod: dotproduct
        in:
          reads_1: fastp/trimmed_1
          reads_2: fastp/trimmed_2
          sample_id: sample_ids
          genome_dir: genome_dir
          threads: threads
        out: [sorted_bam]

      samtools_index:
        run: '#samtools_index'
        scatter: bam
        in:
          bam: star_align/sorted_bam
        out: [indexed_bam]

      samtools_cram:
        run: '#samtools_cram'
        scatter: [bam, sample_id]
        scatterMethod: dotproduct
        in:
          bam: samtools_index/indexed_bam          # arrives with its .bai staged (declared on the tool input)
          reference: reference
          sample_id: sample_ids
          threads: threads
        out: [cram]

  - id: fastp
    class: CommandLineTool
    baseCommand: fastp
    requirements:
      DockerRequirement:
        # Pinned tag for reproducibility; digest-pin (@sha256:) in production. :latest would drift.
        dockerPull: quay.io/biocontainers/fastp:0.23.4--hadf994f_2
      ResourceRequirement:
        coresMin: $(inputs.threads)
        ramMin: 4000                               # MB; fastp is light, 4 GB covers paired FASTQ streaming
    inputs:
      reads_1: {type: File, inputBinding: {prefix: -i}}
      reads_2: {type: File, inputBinding: {prefix: -I}}
      sample_id: string
      threads: {type: int, default: 4, inputBinding: {prefix: --thread}}
    arguments:
      - {prefix: -o, valueFrom: $(inputs.sample_id)_trimmed_R1.fq.gz}
      - {prefix: -O, valueFrom: $(inputs.sample_id)_trimmed_R2.fq.gz}
      - {prefix: --json, valueFrom: $(inputs.sample_id)_fastp.json}
    outputs:
      trimmed_1: {type: File, outputBinding: {glob: '*_trimmed_R1.fq.gz'}}
      trimmed_2: {type: File, outputBinding: {glob: '*_trimmed_R2.fq.gz'}}
      json_report: {type: File, outputBinding: {glob: '*_fastp.json'}}

  - id: star_align
    class: CommandLineTool
    baseCommand: [STAR, --runMode, alignReads]
    requirements:
      DockerRequirement:
        dockerPull: quay.io/biocontainers/star:2.7.11b--h43eeafb_1   # pinned; digest-pin in production
      ResourceRequirement:
        coresMin: $(inputs.threads)
        ramMin: 32000                              # MB; STAR holds the human genome index in RAM (~30 GB)
    inputs:
      reads_1: {type: File, inputBinding: {prefix: --readFilesIn, position: 1}}
      reads_2: {type: File, inputBinding: {position: 2}}
      genome_dir: {type: Directory, inputBinding: {prefix: --genomeDir}}
      sample_id: string
      threads: {type: int, default: 8, inputBinding: {prefix: --runThreadN}}
    arguments:
      - {prefix: --readFilesCommand, valueFrom: zcat}
      - --outSAMtype                               # STAR takes two tokens: --outSAMtype BAM SortedByCoordinate
      - BAM
      - SortedByCoordinate
      - {prefix: --outFileNamePrefix, valueFrom: $(inputs.sample_id)_}
    outputs:
      sorted_bam:
        type: File
        outputBinding: {glob: '*_Aligned.sortedByCoord.out.bam'}

  - id: samtools_index
    class: CommandLineTool
    # Runs `samtools index BAM` and emits the SAME bam carrying its freshly built .bai as a secondaryFile,
    # so downstream random-access tools receive the index without hand-wiring it.
    baseCommand: [samtools, index]
    requirements:
      DockerRequirement:
        dockerPull: quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1   # pinned; digest-pin in production
      InitialWorkDirRequirement:
        listing: [$(inputs.bam)]                   # stage the bam into the workdir so index writes .bai beside it
    inputs:
      bam: {type: File, inputBinding: {position: 1}}
    outputs:
      indexed_bam:
        type: File
        secondaryFiles: [.bai]                     # OUTPUT secondaryFile: the .bai is collected with the bam
        outputBinding: {glob: $(inputs.bam.basename)}

  - id: samtools_cram
    class: CommandLineTool
    baseCommand: [samtools, view, -C]
    requirements:
      DockerRequirement:
        dockerPull: quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1   # pinned; digest-pin in production
      ResourceRequirement:
        coresMin: $(inputs.threads)
        ramMin: 4000
    inputs:
      bam:
        type: File
        secondaryFiles: [.bai]                     # INPUT secondaryFile: the .bai must arrive co-staged
        inputBinding: {position: 1}
      reference:
        type: File
        secondaryFiles:                            # CRAM compression reads the reference via its .fai; .dict travels along
          - {pattern: .fai, required: true}
          - {pattern: ^.dict, required: false}
        inputBinding: {prefix: -T}
      sample_id: string
      threads: {type: int, default: 8, inputBinding: {prefix: -@}}
    arguments:
      - {prefix: -o, valueFrom: $(inputs.sample_id).cram}
    outputs:
      cram:
        type: File
        outputBinding: {glob: '*.cram'}
