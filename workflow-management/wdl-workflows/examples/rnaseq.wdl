version 1.0

# RNA-seq quantification workflow
# Trims reads with fastp and quantifies with Salmon
# Typical runtime: ~30 min per sample at 8 threads
#
# Reproducibility: docker images below are pinned to immutable BioContainers build tags
# (name--buildhash), not :latest. For the strongest guarantee pin by digest
# (org/image@sha256:...): a digest fixes reproducibility AND stabilizes Cromwell call caching,
# which hashes the resolved image identity into the cache key. Lint: `miniwdl check rnaseq.wdl`.

workflow rnaseq_pipeline {
    input {
        Array[String] sample_ids
        Array[File] fastq_1_files
        Array[File] fastq_2_files
        File salmon_index
        Int threads = 8
    }

    # Process each sample in parallel
    scatter (idx in range(length(sample_ids))) {
        call fastp {
            input:
                sample_id = sample_ids[idx],
                reads_1 = fastq_1_files[idx],
                reads_2 = fastq_2_files[idx],
                threads = threads
        }

        call salmon_quant {
            input:
                sample_id = sample_ids[idx],
                reads_1 = fastp.trimmed_1,
                reads_2 = fastp.trimmed_2,
                index = salmon_index,
                threads = threads
        }
    }

    output {
        Array[File] trimmed_r1 = fastp.trimmed_1
        Array[File] trimmed_r2 = fastp.trimmed_2
        Array[File] fastp_reports = fastp.json_report
        Array[File] quant_files = salmon_quant.quant_sf
    }

    meta {
        author: "bioSkills"
        description: "RNA-seq quantification with fastp and Salmon"
    }
}

task fastp {
    input {
        String sample_id
        File reads_1
        File reads_2
        Int threads = 4
    }

    # Dynamic disk: (paired inputs + trimmed outputs) approx 3x, +10 GiB headroom.
    # ceil() rounds UP so disk never under-sizes and kills the task mid-run.
    Int disk_gb = ceil((size(reads_1, "GiB") + size(reads_2, "GiB")) * 3) + 10

    command <<<
        fastp \
            -i ~{reads_1} \
            -I ~{reads_2} \
            -o ~{sample_id}_trimmed_R1.fq.gz \
            -O ~{sample_id}_trimmed_R2.fq.gz \
            --json ~{sample_id}_fastp.json \
            --html ~{sample_id}_fastp.html \
            --thread ~{threads}
    >>>

    output {
        File trimmed_1 = "~{sample_id}_trimmed_R1.fq.gz"
        File trimmed_2 = "~{sample_id}_trimmed_R2.fq.gz"
        File json_report = "~{sample_id}_fastp.json"
        File html_report = "~{sample_id}_fastp.html"
    }

    runtime {
        docker: "quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
        cpu: threads
        memory: "4 GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: 3
    }
}

task salmon_quant {
    input {
        String sample_id
        File reads_1
        File reads_2
        File index
        Int threads = 8
    }

    # Salmon needs ~8-16GB RAM depending on index; disk for index + reads + outputs
    Int disk_gb = ceil(size(index, "GiB") + size(reads_1, "GiB") * 2) + 20

    command <<<
        salmon quant \
            -i ~{index} \
            -l A \
            -1 ~{reads_1} \
            -2 ~{reads_2} \
            -o ~{sample_id}_salmon \
            --threads ~{threads}
    >>>

    output {
        File quant_sf = "~{sample_id}_salmon/quant.sf"
        File quant_genes = "~{sample_id}_salmon/quant.genes.sf"
        File cmd_info = "~{sample_id}_salmon/cmd_info.json"
    }

    runtime {
        docker: "quay.io/biocontainers/salmon:1.10.0--h7e5ed60_0"
        cpu: threads
        # Salmon memory scales with index size; 16GB safe for most transcriptomes
        memory: "16 GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: 3
    }
}
