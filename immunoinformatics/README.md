# immunoinformatics

## Overview

Computational immunology: MHC class I and class II binding/presentation, tumor neoantigen identification and immunogenicity ranking, B/T-cell epitope prediction, and TCR-antigen specificity for vaccine design and cancer immunotherapy. Skills are decision-grade and encode where each prediction is trustworthy and where it silently fails.

**Tool type:** mixed | **Primary tools:** MHCflurry, NetMHCpan, NetMHCIIpan, pVACtools, NeoFox, BepiPred, DiscoTope, tcrdist3

## Skills

| Skill | Description |
|-------|-------------|
| mhc-binding-prediction | Class I peptide-MHC binding/presentation (MHCflurry, NetMHCpan-4.1); BA vs EL, %Rank vs nM, abundance bias |
| mhc-class-ii-prediction | Class II (DR/DQ/DP) CD4 epitope prediction (NetMHCIIpan-4.3, MixMHC2pred); register, DQ pairing, why it is less reliable |
| neoantigen-prediction | Tumor neoantigens with pVACtools; clonality/CCF, HLA LOH, phasing, the downstream-attrition reality |
| immunogenicity-scoring | Rank candidates by likely T-cell response (NeoFox, PRIME2.0, fitness model); rank-don't-threshold |
| epitope-prediction | B-cell (BepiPred-3.0/DiscoTope-3.0) and T-cell epitopes; the maturity asymmetry, conformational dominance |
| tcr-epitope-binding | TCR specificity by clustering (tcrdist3/GLIPH2) + lookup; why de-novo prediction for unseen epitopes fails |

## Example Prompts

- "Predict class I presentation for these peptides with my patient's HLA genotype"
- "Predict CD4 epitopes against this DRB1/DQ genotype and tell me how much to trust DQ calls"
- "Find neoantigens from my somatic VCF and drop candidates on HLA-LOH-lost alleles"
- "Rank these neoantigens by immunogenicity within the patient, keeping features visible"
- "Identify conformational B-cell epitopes on this antigen from its AlphaFold model"
- "Cluster my TCR repertoire and annotate clusters with known specificities"

## Requirements

```bash
pip install mhcflurry pvactools vatools neofox bepipred3 tcrdist3
mhcflurry-downloads fetch
# NetMHCpan-4.1, NetMHCIIpan-4.3, MixMHCpred, MixMHC2pred, DiscoTope-3.0, PRIME, BigMHC
# are standalone academic/lab downloads (DTU Health Tech / Gfeller lab / Karchin lab) or
# wrapped by the IEDB REST API. Ensembl VEP + an HLA typer (OptiType/arcasHLA/HLA-HD) are
# needed for the neoantigen pipeline.
```

## Related Skills

- **clinical-databases** - HLA typing, ClinVar/gnomAD variant annotation, somatic signatures
- **variant-calling** - somatic SNV/indel calls feeding neoantigen prediction
- **tcr-bcr-analysis** - TCR/BCR repertoire sequencing upstream of specificity analysis
- **structural-biology** - AlphaFold structures enabling conformational B-cell epitope prediction
- **workflows** - the neoantigen-pipeline and tcr-pipeline end-to-end orchestrations
