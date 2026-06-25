---
name: bio-rna-structure-structure-probing
description: Processes experimental RNA structure probing data (SHAPE-MaP, DMS-MaPseq) into per-nucleotide reactivity profiles with ShapeMapper2, then uses them as soft restraints on thermodynamic folding. Covers reagent and readout choice (SHAPE vs DMS, mutational-profiling vs RT-stop), the three control samples, per-transcript normalization, the Deigan vs Zarringhalam pseudo-energy models, in-cell versus in-vitro interpretation, and multi-conformation deconvolution. Use when converting probing reads to reactivities; deciding SHAPE versus DMS parameters; judging whether low reactivity means base-paired or protein-bound; or detecting whether an RNA populates more than one structure.
tool_type: cli
primary_tool: ShapeMapper2
---

## Version Compatibility

Reference examples tested with: ShapeMapper2 2.1.5+, ViennaRNA 2.6+, SEISMIC-RNA 0.20+, matplotlib 3.8+, pandas 2.2+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structure Probing

**"Process my SHAPE-MaP experiment to get RNA reactivity profiles"** -> Convert per-nucleotide chemical-modification signal into reactivities, normalize them, then feed them as SOFT restraints into thermodynamic folding.
- CLI: `shapemapper` (ShapeMapper2) for end-to-end SHAPE-MaP / DMS-MaP processing
- CLI: `RNAfold --shape` (ViennaRNA) for reactivity-restrained folding
- CLI: `seismic` (SEISMIC-RNA) for DMS-MaPseq and multi-conformation clustering

## The governing principle: reactivity is flexibility, not pairing, and it is a constraint, not a structure

SHAPE reagents acylate the ribose 2'-OH at a rate set by local nucleotide FLEXIBILITY / conformational dynamics; DMS methylates the Watson-Crick face of A (N1) and C (N3). High reactivity means flexible/accessible, low means constrained. The routine over-interpretation is "low reactivity = base-paired": a nucleotide can be unreactive because it is base-paired OR because it is tertiary-contacted, protein-bound, ligand-occluded, or stacked. Reactivity probes CONSTRAINT, and it cannot by itself distinguish pairing from other protection. Before interpreting any low-reactivity region, rule out the TECHNICAL causes first: a long flat run can be low read depth / no-data or global undermodification, not structure -- check effective depth and confirm low-depth positions are carried as -999, not 0.

Two consequences govern everything below:
- Reactivities are a RESTRAINT on folding, not a structure. They enter as pseudo-free-energy terms that bias the thermodynamic fold and raise accuracy substantially, but they do not yield a structure on their own and cannot disambiguate pairing from other protection.
- A per-position profile is a POPULATION AVERAGE. If the RNA samples more than one structure (riboswitches, dynamic mRNAs), the averaged profile can match no real structure. Detecting and deconvolving multiple conformations requires per-read data (mutational profiling) and clustering (SEISMIC-RNA / DREEM), not a single profile.

## Reagent and readout choices

SHAPE acts on the backbone 2'-OH, so it reports all four bases; DMS reads only A/C by default (G/U carry no Watson-Crick-face signal and must be masked to no-data before folding). Newer DMS-MaPseq protocols recover G(N1)/U(N3) signal for a four-base readout (via mutation-signature filtering plus optimized buffer, not the buffer alone) -- only treat DMS as four-base if the protocol and analysis explicitly enable it; otherwise mask G/U.

| Goal | Reagent | Reads | Note |
|------|---------|-------|------|
| In-vitro, all-base, fast | 1M7 (SHAPE) | A/C/G/U flexibility | default; m/b = 1.8/-0.6 |
| In-cell SHAPE (membrane-permeable) | NAI, NAI-N3, 5NIA, 2A3 | A/C/G/U | NAI-N3 -> icSHAPE click enrichment; 2A3 among the strongest in vivo |
| In-cell, cheap, A/C-resolved | DMS | A(N1)/C(N3) | works in vivo; mask G/U to no-data |
| Tertiary-contact flagging | 1M6/NMIA vs 1M7 (differential SHAPE) | A/C/G/U | report the DIFFERENCE, not absolute |
| Fill G/U coverage | CMCT (G/U), kethoxal (G) | G,U / G | low throughput, rarely MaP-coupled |

Mutational profiling (MaP) vs RT-stop is a fundamental analysis fork, not a detail: in MaP the reverse transcriptase reads THROUGH the adduct (Mn2+/TGIRT/Marathon-RT) and misincorporates, encoding each modification as a point mutation; in RT-stop the adduct truncates the cDNA and the read 5'-end is counted. They need different scoring, and a pipeline tuned for one is wrong for the other.

| Axis | MaP (mutation) | RT-stop (truncation) |
|------|----------------|----------------------|
| Adduct encoded as | misincorporation/deletion | RT drop-off (read 5'-end) |
| Reads/molecule informative | many | one |
| Single-molecule / correlated analysis | yes (per-read mutation strings) | no |
| Tools | ShapeMapper2, SEISMIC-RNA, rf-count -sm 3/4 | rf-count -sm 1/2, icSHAPE, StructureFold |
| Methods | SHAPE-MaP, DMS-MaPseq | Mod-seq, Structure-seq, DMS-seq, icSHAPE |

icSHAPE / Structure-seq / DMS-seq are RT-stop methods; do not push them through a MaP mutation-rate pipeline (ShapeMapper2/SEISMIC).

## ShapeMapper2: the three samples and verified flags

The three sample roles are distinct: MODIFIED is the signal; UNTREATED subtracts background (SNPs, RT errors, intrinsic damage); DENATURED normalizes sequence-dependent reactivity bias (divide modified by denatured when no good no-data reference exists). The untreated control is mandatory; the denatured control improves normalization.

```bash
# ShapeMapper2 is Linux-only; on macOS use Docker/Singularity (see usage-guide).
shapemapper \
    --target target_rna.fa \
    --name my_rna \
    --modified --R1 mod_R1.fastq.gz --R2 mod_R2.fastq.gz \
    --untreated --R1 unmod_R1.fastq.gz --R2 unmod_R2.fastq.gz \
    --out results/ \
    --nproc 8 \
    --min-depth 5000
```

| Option | Effect (verified default) |
|--------|---------------------------|
| `--target` / `--name` | reference FASTA / output basename |
| `--modified` / `--untreated` / `--denatured` | the three samples, each followed by `--R1/--R2` |
| `--amplicon` | primer-trimmed amplicon mode |
| `--min-depth` | minimum effective depth to report a nt (default 5000) |
| `--min-qual-to-count` | minimum basecall quality in a mutation (default 30, not 20) |
| `--max-bg` | max untreated mutation frequency (default 0.05) |
| `--star-aligner` | use STAR instead of Bowtie2 (recommended for long targets) |
| `--nproc` / `--overwrite` | threads / overwrite output |

ShapeMapper2 writes a `results/` tree. The reactivity table is `<name>_<RNA>_profile.txt`; the folder-ready files are SEPARATE: `<name>_<RNA>.shape` (2 columns: position, normalized reactivity; excluded = -999) and `<name>_<RNA>.map` (4 columns: position, normalized reactivity, stderr, base). There is no combined `_map.shape` file.

Key `profile.txt` columns: `Nucleotide` is the POSITION integer (1-based); `Sequence` is the base character; `Reactivity_profile` is raw; `Norm_profile` is after normalization. Fold with `Norm_profile` (or the `.shape`/`.map` file), never the raw `Reactivity_profile` -- the pseudo-energy parameters assume normalized input.

## Normalization: per-transcript, and why raw reactivities are not comparable

Raw reactivity (background-subtracted modified mutation rate) sits on an arbitrary, experiment-specific scale set by reagent dose, RT efficiency, and depth. The standard 2-8% / box-plot normalization excludes outliers (top ~2% as a whisker cap), then scales by the mean of the next most-reactive ~8-10% of nucleotides, so normalized values mostly fall ~0-2 with ~1.0 = average reactivity. This scale factor is PER TRANSCRIPT: raw reactivities from different transcripts or experiments are NOT comparable, so never pool or compare raw reactivities across them. To compare two conditions (e.g. +/- ligand), use delta-SHAPE at matched positions with the per-nt standard errors, not raw subtraction. Low-depth nucleotides must become no-data (-999), not zero.

## Reactivity-restrained folding and the pseudo-energy model

The Deigan model adds a soft pseudo-energy to every nucleotide in a stacked pair: deltaG = m * ln(1 + reactivity) + b. It is a restraint, not a hard constraint -- a nucleotide can still pair against the data if the global fold demands it.

| Situation | Model / ViennaRNA flag | Parameters |
|-----------|------------------------|------------|
| Standard SHAPE (1M7/NAI) | Deigan `--shapeMethod="Dm1.8b-0.6"` | m=1.8, b=-0.6 (Hajdin 2013) |
| Noisy data / probabilistic target | Zarringhalam `--shapeMethod="Z"` | target pairing probability |
| Penalize unpaired only | Washietl `--shapeMethod="W"` | perturbation vector |
| DMS-MaPseq | Deigan-style, A/C only, G/U set to -999 | no DMS-specific standard; commonly reuse 1.8/-0.6, or tune |

The m=1.8, b=-0.6 pair is the Hajdin et al. 2013 standard and the ViennaRNA "Deigan" DEFAULT -- it is NOT Deigan et al. 2009's own values (m=2.6, b=-0.8); cite it correctly. For DMS, apply the restraint ONLY to A/C and set G/U to -999, or the model invents constraints at bases that carry no signal. The folded energy a tool reports after a SHAPE/DMS restraint INCLUDES the pseudo-energy bonus, so it is not comparable to an unrestrained MFE -- judge the result by structure agreement, not by a more-negative energy.

```bash
# Fold directly from the ShapeMapper2 .shape file (already normalized)
RNAfold --shape=results/my_rna_my_rna.shape --shapeMethod="Dm1.8b-0.6" --noPS < target_rna.fa
```

## In-cell vs in-vitro: the occupancy trap

In-vitro (refolded, deproteinized) RNA reports pure thermodynamics; in-cell reports the RNA as it exists, with bound proteins, ligands, and chaperone-remodeled states all altering reactivity. An in-cell PROTECTED nucleotide may be protein-bound or ligand-occluded, not base-paired -- the single biggest in-cell misinterpretation. Cells also actively unfold mRNA: genome-wide in-vivo DMS showed mRNAs are MORE unfolded in vivo than in vitro (Rouskin 2014). The in-cell-minus-in-vitro difference is itself the signal for protein/ligand footprints (Spitale 2015).

| Want | Condition | Reagent | Caveat |
|------|-----------|---------|--------|
| De-novo thermodynamic structure | in-vitro refolded | 1M7 SHAPE / DMS | MFE-like, no proteins |
| Functional in-cell state | in-cell | NAI/2A3/5NIA, DMS | protected != paired (occupancy) |
| Protein/ligand footprint | in-cell vs in-vitro delta | matched reagent | needs both, matched depth |

## Multiple conformations: cluster, do not average

If a profile looks inconsistent with any single structure, the RNA may populate more than one. DREEM (Tomezsko 2020) and its maintained successor SEISMIC-RNA cluster MaP reads by co-occurring mutations (expectation-maximization) to deconvolve coexisting conformers; RING-MaP (Homan 2014) and PAIR-MaP (Mustoe 2019) use correlated mutations between positions to detect through-space communication and direct base pairs. These need per-read mutation data, which only MaP provides.

```bash
# DMS-MaPseq processing and multi-conformation clustering with SEISMIC-RNA
# Subcommand names vary by version (released: align/relate/mask/cluster; recent dev renames
# relate->idmut, mask->filter). Run `seismic --help` to confirm before scripting.
seismic align target.fa reads_R1.fq.gz reads_R2.fq.gz --out seismic_out
seismic relate seismic_out target.fa --out seismic_out
seismic mask seismic_out --out seismic_out
seismic cluster seismic_out --max-clusters 3 --out seismic_out
```

For RNA Framework (rf-count -> rf-norm), the reference is `-f` and the BAM/SAM files are positional; choose the MaP-vs-RT-stop scoring in rf-norm with `-sm` (1 Ding RT-stop, 2 Rouskin RT-stop, 3 Siegfried MaP, 4 Zubradt MaP) and the normalization with `-nm` (1 = 2-8% default, 3 = box-plot); restrict reactive bases for DMS with `-rb AC`.

```bash
rf-count -f reference.fa modified.bam untreated.bam -o rf_out/
rf-norm -i rf_out/index.rci -t rf_out/modified.rc -u rf_out/untreated.rc -sm 3 -nm 1 -rb AC
```

## Quality thresholds

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| Effective depth | >= 5000 | reliable per-nt mutation-rate estimation for MaP |
| Untreated mutation rate | < 0.5% | overall expectation; higher suggests SNP, RT-prone motif, or damage (individual nt above `--max-bg`=5% are auto-excluded) |
| Modified mutation rate | ~1-10% | too low = undermodified; too high = degraded |
| No-data marker | -999 | low-depth/high-background nt; carry through folding, do not treat as 0 |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `KeyError: 'Reactivity_profile'` or garbled sequence | reading base from `Nucleotide` (it is the position integer) | read the base from `Sequence`; fold with `Norm_profile` |
| `FileNotFoundError: my_rna_map.shape` | no combined file is written | use the separate `<name>_<RNA>.shape` (2-col) and `.map` (4-col) |
| Folding barely changes with SHAPE data | folding with raw `Reactivity_profile`, or vector mis-indexed | use `Norm_profile`; vector is 1-indexed (prepend -999), -999 = no data |
| DMS constraints look noisy at G/U | G/U carry no Watson-Crick DMS signal | mask G/U to -999 before folding and before normalization |
| `rf-count -t target.fa -r mod.bam -rc unt.bam` errors | wrong flags | reference is `-f`; BAMs are positional; there is no `-r`/`-rc` |
| Two conditions disagree but raw reactivities were compared | raw values are per-transcript, non-comparable | compare normalized profiles (delta-SHAPE) with standard errors |
| Profile fits no single structure | RNA populates multiple conformations | cluster MaP reads with SEISMIC-RNA / DREEM |
| In-cell protected region called "paired" | protection may be protein/ligand occupancy | compare in-cell vs in-vitro; do not equate protection with pairing |
| A long unreactive stretch read as a stable hairpin | could be low depth/no-data or undermodification, not pairing | check effective depth (>=5000) and that low-depth nt are -999 before interpreting |

## Related Skills

- secondary-structure-prediction - The folding engine the reactivities restrain
- ncrna-search - Identify the RNA family and a CM consensus structure to probe against
- covariation-analysis - Independent (evolutionary) evidence for the pairs probing suggests
- epitranscriptomics/m6a-peak-calling - RNA modifications that confound DMS/SHAPE reactivity
- clip-seq/binding-site-annotation - In-cell protection as an RBP footprint
- read-qc/quality-reports - QC of the underlying sequencing reads

## References

- Merino EJ, Wilkinson KA, Coughlan JL, Weeks KM. 2005. RNA structure analysis at single nucleotide resolution by selective 2'-hydroxyl acylation and primer extension (SHAPE). J Am Chem Soc 127(12):4223-4231. doi:10.1021/ja043822v
- Mortimer SA, Weeks KM. 2007. A fast-acting reagent for accurate analysis of RNA secondary and tertiary structure by SHAPE chemistry. J Am Chem Soc 129(14):4144-4145. doi:10.1021/ja0704028
- Deigan KE, Li TW, Mathews DH, Weeks KM. 2009. Accurate SHAPE-directed RNA structure determination. Proc Natl Acad Sci USA 106(1):97-102. doi:10.1073/pnas.0806929106
- Zarringhalam K, Meyer MM, Dotu I, Chuang JH, Clote P. 2012. Integrating chemical footprinting data into RNA secondary structure prediction. PLoS ONE 7(10):e45160. doi:10.1371/journal.pone.0045160
- Hajdin CE, Bellaousov S, Huggins W, Leonard CW, Mathews DH, Weeks KM. 2013. Accurate SHAPE-directed RNA secondary structure modeling, including pseudoknots. Proc Natl Acad Sci USA 110(14):5498-5503. doi:10.1073/pnas.1219988110
- Cordero P, Kladwang W, VanLang CC, Das R. 2012. Quantitative dimethyl sulfate mapping for automated RNA secondary structure inference. Biochemistry 51(36):7037-7039. doi:10.1021/bi3008802
- Rouskin S, Zubradt M, Washietl S, Kellis M, Weissman JS. 2014. Genome-wide probing of RNA structure reveals active unfolding of mRNA structures in vivo. Nature 505(7485):701-705. doi:10.1038/nature12894
- Homan PJ, Favorov OV, Lavender CA, Kursun O, Ge X, Busan S, Dokholyan NV, Weeks KM. 2014. Single-molecule correlated chemical probing of RNA. Proc Natl Acad Sci USA 111(38):13858-13863. doi:10.1073/pnas.1407306111
- Siegfried NA, Busan S, Rice GM, Nelson JAE, Weeks KM. 2014. RNA motif discovery by SHAPE and mutational profiling (SHAPE-MaP). Nat Methods 11(9):959-965. doi:10.1038/nmeth.3029
- Smola MJ, Rice GM, Busan S, Siegfried NA, Weeks KM. 2015. Selective 2'-hydroxyl acylation analyzed by primer extension and mutational profiling (SHAPE-MaP) for direct, versatile and accurate RNA structure analysis. Nat Protoc 10(11):1643-1669. doi:10.1038/nprot.2015.103
- Spitale RC, Flynn RA, Zhang QC, Crisalli P, Lee B, Jung JW, Kuchelmeister HY, Batista PJ, Torre EA, Kool ET, Chang HY. 2015. Structural imprints in vivo decode RNA regulatory mechanisms. Nature 519(7544):486-490. doi:10.1038/nature14263
- Zubradt M, Gupta P, Persad S, Lambowitz AM, Weissman JS, Rouskin S. 2017. DMS-MaPseq for genome-wide or targeted RNA structure probing in vivo. Nat Methods 14(1):75-82. doi:10.1038/nmeth.4057
- Busan S, Weeks KM. 2018. Accurate detection of chemical modifications in RNA by mutational profiling (MaP) with ShapeMapper 2. RNA 24(2):143-148. doi:10.1261/rna.061945.117
- Incarnato D, Morandi E, Simon LM, Oliviero S. 2018. RNA Framework: an all-in-one toolkit for the analysis of RNA structures and post-transcriptional modifications. Nucleic Acids Res 46(16):e97. doi:10.1093/nar/gky486
- Mustoe AM, Lama NN, Irving PS, Olson SW, Weeks KM. 2019. RNA base-pairing complexity in living cells visualized by correlated chemical probing (PAIR-MaP). Proc Natl Acad Sci USA 116(49):24574-24582. doi:10.1073/pnas.1905491116
- Tomezsko PJ, Corbin VDA, Gupta P, Swaminathan H, Glasgow M, Persad S, Edwards MD, Rouskin S. 2020. Determination of RNA structural diversity and its role in HIV-1 RNA splicing (DREEM). Nature 582(7812):438-442. doi:10.1038/s41586-020-2253-5
