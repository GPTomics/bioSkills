# Amplicon Processing with DADA2

## Overview

DADA2 infers exact amplicon sequence variants (ASVs) from amplicon sequencing data, providing single-nucleotide resolution without clustering.

## ASVs vs OTUs

| Feature | ASVs (DADA2) | OTUs (97%) |
|---------|--------------|------------|
| Resolution | Single nucleotide | 3% similarity |
| Reproducibility | Exact sequences | Cluster-dependent |
| Comparability | Can merge studies | Requires re-clustering |
| Sensitivity | Higher | Lower |

## Workflow Steps

1. **Inspect quality profiles** - Determine trim points
2. **Filter and trim** - Remove low-quality bases
3. **Learn error rates** - Model sequencing errors
4. **Denoise** - Infer true sequences
5. **Merge pairs** - Combine forward/reverse
6. **Remove chimeras** - Eliminate PCR artifacts

## Key Parameters

### filterAndTrim

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| truncLen | Truncate reads at position | Based on quality plots |
| maxEE | Maximum expected errors | 2 |
| maxN | Maximum Ns allowed | 0 |
| truncQ | Truncate at quality score | 2 |

### Choosing truncLen

- Examine quality profiles
- Trim where median quality drops below Q30
- Ensure sufficient overlap for merging (20bp minimum)

## Common Issues

### Low merge rate
- Check amplicon length vs read length
- Increase minimum overlap
- Reads may not overlap (too long amplicon)

### High chimera rate
- Normal: 10-25% of ASVs, but few reads
- Check PCR cycles (lower is better)

### Few reads passing filter
- Loosen maxEE (try 2, 5)
- Adjust truncLen
- Check raw data quality

## Output Files

```r
# Save for downstream analysis
saveRDS(seqtab_nochim, 'seqtab_nochim.rds')

# Export as FASTA
seqs <- getSequences(seqtab_nochim)
headers <- paste0('ASV', seq_along(seqs))
writeXStringSet(DNAStringSet(seqs), 'asv_seqs.fasta')
```
