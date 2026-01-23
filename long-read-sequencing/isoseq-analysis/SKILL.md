---
name: bio-long-read-sequencing-isoseq-analysis
description: Analyze PacBio Iso-Seq data for full-length isoform discovery and quantification. Use when characterizing transcript diversity or identifying novel splice variants.
tool_type: cli
primary_tool: IsoSeq3
---

# Iso-Seq Analysis

## IsoSeq3 Pipeline

```bash
# CCS reads (circular consensus)
ccs input.subreads.bam ccs.bam --min-rq 0.9

# Remove primers
lima ccs.bam primers.fasta demux.bam --isoseq

# Refine (remove polyA, concatemers)
isoseq3 refine demux.primer_5p--primer_3p.bam primers.fasta refined.bam

# Cluster into isoforms
isoseq3 cluster refined.bam clustered.bam --verbose
```

## SQANTI3 Quality Control

```bash
# Classify isoforms against reference
sqanti3_qc.py \
    clustered.hq_transcripts.fasta \
    reference.gtf \
    reference.fa \
    -o sqanti_output

# Categories: FSM (full splice match), ISM (incomplete),
# NIC (novel in catalog), NNC (novel not in catalog)
```

## Collapse and Quantify

```bash
# Collapse redundant isoforms
isoseq3 collapse clustered.bam reference.fa collapsed.gff

# Quantify with minimap2
minimap2 -ax splice:hq -uf reference.fa reads.fq > aligned.sam
```

## TAMA for Annotation

```bash
# Merge with reference annotation
tama_merge.py \
    -f file_list.txt \
    -p merged \
    -a 50 -z 50
```

## Related Skills

- **long-read-sequencing/basic-analysis** - ONT/PacBio basics
- **rna-quantification** - Expression analysis
- **genome-annotation** - GTF/GFF handling
