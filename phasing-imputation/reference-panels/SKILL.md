---
name: bio-phasing-imputation-reference-panels
description: Selects and prepares the reference panel that phasing/imputation copies haplotypes from (1000 Genomes, HRC, TOPMed, HGDP+1kGP/gnomAD, CAAPA), matching panel ancestry to the target, reconciling genome build and chromosome naming, and running the strand/allele harmonization gate. Covers why ancestry-match beats panel size (imputation can only copy haplotypes the panel contains), why palindromic A/T and C/G SNPs flip strand without erroring, why liftover is a strand-flip generator in between-build inverted regions, that HRC is SNP-only and TOPMed is never downloadable (governance can override accuracy), and panel formats (msav, bref3, imp5). Use when choosing a panel for a target ancestry, preparing or converting a panel, aligning study data, or deciding between downloadable and server-only panels. Phasing is haplotype-phasing; imputation is genotype-imputation; PCA for ancestry is population-genetics/population-structure; HLA panels are clinical-databases/hla-typing.
tool_type: cli
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, PLINK 1.9+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

A reference panel is DATA, not a tool: it has a dated release and a fixed genome build (GRCh37 or GRCh38). Record the exact panel name, version, and build with every result; "1000 Genomes" without a version is unreproducible because Phase 3 (GRCh37, low-coverage) and the high-coverage NYGC 3202 release (GRCh38, 30x) are different call sets. HRC and TOPMed are server-only (not downloadable); panel-build tools (Minimac4 `--compress-reference`, `bref3.jar`, `imp5Converter`) and the Will Rayner harmonization check are separate downloads.

# Reference Panels -- Choosing and Preparing the Prior Imputation Copies From

**"Which reference panel should I use, and how do I prepare it?"** -> Match the panel's ancestry to the target population, reconcile build and strand, then convert to the engine's format - because imputation can only copy haplotypes the panel contains, so the panel IS the prior, and a mismatched ancestry or a flipped strand corrupts the result without any error.
- CLI: `bcftools norm -m -any -f ref.fa` then the strand/allele harmonization check against the panel sites, then `minimac4 --compress-reference` / `bref3.jar` / `imp5Converter` to build the engine format

Scope: panel selection (ancestry-match, build, access), strand/allele harmonization, build/liftover, and format conversion. The phasing engine that consumes the panel -> haplotype-phasing. Imputation -> genotype-imputation. PCA to establish the target ancestry -> population-genetics/population-structure. Classical HLA-allele imputation needs a dedicated HLA panel -> clinical-databases/hla-typing. VCF normalization mechanics -> variant-calling/variant-normalization.

## The Single Most Important Modern Insight -- Ancestry Match Beats Panel Size, Because Imputation Copies Haplotypes and a Mismatched Panel Has None Worth Copying

The reflex "TOPMed has 97k samples, HRC has 32k, so use TOPMed" is right for a European cohort and wrong for an ancestry-mismatched one, and the reason is mechanical, not statistical: imputation copies haplotype segments from panel samples that resemble the target, so if no panel sample carries the target population's haplotypes there is nothing to copy, and adding ten thousand more European haplotypes does nothing for an East African sample (Marchini & Howie 2010 *Nat Rev Genet* 11:499). The binding resource is not panel size but how many panel samples share the target's ancestry. Three facts follow:

1. **Ancestry composition often dominates size.** HRC (32k, European-heavy) imputes African-ancestry rare variants poorly while TOPMed lifts them dramatically - not because TOPMed is 3x bigger but because it contains African-American and Hispanic/Latino haplotypes. Admixed samples need a panel with both ancestral components and admixed individuals, because admixed haplotypes are mosaics a single-ancestry panel cannot reconstruct at the switch points. The honest answer is "large AND matched"; where both are not available, which one wins depends on whether the target is common or rare variants (size matters more for rare).
2. **The self-graded INFO/R2 cannot see ancestry mismatch.** When the panel lacks the target's haplotypes the model still finds some template and reports a confident-looking quality about a copy that is systematically wrong (full metric theory -> imputation-qc). High R2 in an under-represented ancestry is not reassurance; validate against masked truth stratified by frequency.
3. **A panel encodes whose haplotypes are trusted to fill the gaps.** HRC is European-heavy because its component cohorts were; TOPMed is diverse because it was designed to be. No panel represents everyone, and the imputation inherits exactly the panel population's representation, including its gaps. This is the central equity problem of imputation, and low-coverage WGS (unbiased ascertainment -> genotype-imputation) is the main route around it.

## The Major Panels

Numbers are routinely misquoted; these are the verified figures. State the exact version and build in any method.

| Panel | Samples | Build | Indels | Diversity | Access |
|-------|---------|-------|--------|-----------|--------|
| 1000G Phase 3 (Auton 2015) | 2,504 | GRCh37 (GRCh38 lift) | yes | 26 pops, broad but shallow | public download |
| 1000G high-cov NYGC (Byrska-Bishop 2022) | 3,202 (incl. 602 trios) | GRCh38 | yes | same 26 pops, 30x | public download |
| HRC r1.1 (McCarthy 2016) | 32,470 (64,940 haps) | GRCh37 only | NO - SNP-only, MAF floor ~5e-4 | European-heavy | server-only (Michigan) |
| TOPMed r2 (Taliun 2021) | 97,256 (~308M sites) | GRCh38 only | yes | very diverse (large AA, Hispanic) | server-only, never downloadable |
| HGDP+1kGP / gnomAD (Koenig 2024) | 4,094 (76-80 pops) | GRCh38 | yes | maximally diverse per-sample | public download |
| CAAPA (Mathias 2016) | 883 African-ancestry | GRCh37 | (SNP) | African / African-American | server-supported |

The "1000G" trap: Phase 3 (Auton 2015, GRCh37, low-coverage) and the high-coverage NYGC 3202 release (Byrska-Bishop 2022, GRCh38, 30x, with 602 trios) are different call sets; the NYGC release is strictly better for rare variants. State which.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Study is GRCh37 | HRC r1.1, 1000G Phase 3, or CAAPA (all GRCh37) | match the build; avoid liftover |
| Study is GRCh38 | TOPMed, 1000G NYGC, or HGDP+1kGP (all GRCh38) | match the build; avoid liftover |
| Build mismatch unavoidable | lift over ONCE, strand-aware, then re-run the harmonization check | liftover flips strand in inverted regions (see failure modes) |
| European cohort, common-variant GWAS | HRC (server) or 1000G | large and European-rich |
| African / admixed / Hispanic / multi-ancestry | TOPMed (server) | diversity wins; contains the matching haplotypes |
| Need a downloadable, diverse, local panel | HGDP+1kGP (gnomAD) | global, jointly-called, and not server-gated |
| Data cannot leave the institution / country | downloadable panels only (1000G, HGDP+1kGP) | governance overrides accuracy (see failure modes) |
| Need indels imputed | 1000G, TOPMed, or HGDP+1kGP | HRC is SNP-only |
| Classical HLA alleles | -> clinical-databases/hla-typing (dedicated HLA panel) | standard SNP panels cannot impute HLA alleles |
| Establish the target ancestry first | -> population-genetics/population-structure | PCA, not a panel operation |

The governing principle: ancestry match beats panel size, and governance (can the data be uploaded to a US server?) often narrows the field before accuracy does.

## The Strand / Allele Harmonization Gate

This is the step that silently corrupts results when skipped. The job: align every study variant's alleles to the panel REF/ALT, fix strand, and drop the variants that cannot be safely resolved. Normalize first (`bcftools norm -m -any -f ref.fa`), because the same indel represented two ways will not match.

The field-standard gate is Will Rayner's check (`HRC-1000G-check-bim.pl` and its bgen/VCF variants): it compares a QC'd PLINK `.bim` plus an allele-frequency file against the panel's sites list and EMITS a `Run-plink.sh` that updates positions, ref/alt, and strand, removes unresolvable SNPs, and splits by chromosome. The check diagnoses; the script fixes - both are required, and a surprising number of pipelines run the check and never execute the script.

The allele-frequency concordance plot (study AF vs panel AF) is the visual gate, not decoration: a tight diagonal is good; points on the `y = 1 - x` anti-diagonal are strand flips; a general smear is sample mislabeling or the wrong panel ancestry. Compare against the ancestry-MATCHED sub-panel's frequencies - an African cohort vs a European panel's AF smears even with perfect strand. Read this plot before uploading, every time.

## Panel Formats and the Genetic-Map Pairing

| Engine | Native format | Build command |
|--------|---------------|---------------|
| Minimac4 | `.msav` (current) or legacy `.m3vcf` | `minimac4 --compress-reference ref.vcf.gz > ref.msav` (legacy: `Minimac3 --processReference`) |
| Beagle 5.x | `.bref3` or plain VCF | `java -jar bref3.jar ref.vcf.gz > ref.bref3` (needs fully phased, non-missing, `|`-separated, per chromosome) |
| IMPUTE5 | `.imp5` or VCF/BCF | `imp5Converter --h ref.vcf.gz --r chr20 --o ref.chr20.imp5` |

A panel is half the input; phasing/imputation also needs a genetic (recombination) map matched to the panel build. A panel VCF comes site-only (the legend, used for the harmonization check) or full (the haplotypes, used to impute) - obtain both. The map is build-specific: an hg19 map with a GRCh38 panel silently mis-places recombination rates.

## chrX and MHC

- **chrX**: split into PAR1 / nonPAR / PAR2 using build-correct coordinates. PAR and female nonPAR are diploid; male nonPAR is HAPLOID and must be coded haploid (a het call there is an error). Mixed ploidy in one file crashes most tools; the Michigan/TOPMed servers split, impute, and re-merge automatically.
- **MHC (chr6 ~28-34 Mb)**: extreme LD and polymorphism. Standard panels impute SNPs there but unreliably and cannot impute classical HLA alleles at all - that needs a dedicated HLA panel and tool -> clinical-databases/hla-typing.

## Per-Method Failure Modes

### Palindromic (A/T, C/G) SNP strand flip
**Trigger:** keeping strand-ambiguous SNPs without a frequency-based strand check. **Mechanism:** A/T and C/G alleles are their own reverse complements, so opposite-strand study and panel still "match" on alleles; the variant passes every join and imputes cleanly while allele-swapped. **Symptom:** flipped effect direction at that locus and everything imputed in LD with it; no error. **Fix:** resolve strand by allele frequency; drop palindromic SNPs with MAF > 0.4 (cannot be disambiguated near 0.5); treat "I kept all palindromic SNPs" as proof strand was never checked.

### Liftover across builds
**Trigger:** running `liftOver` to reach a panel in the other build and imputing without re-checking. **Mechanism:** ~2-5 Mb of the genome is inverted between GRCh37 and GRCh38 (BBIS regions); lifting a variant there changes its strand, and the allele-based check cannot see it on a palindrome (Sheng & Chiang 2023 *HGG Adv* 4:100159). **Symptom:** silent allele-swaps in inverted regions; the TOPMed server's own internal conversion had this bug. **Fix:** prefer a panel native to the study build; if forced, lift once with a strand-aware method, then re-run the harmonization check against the new build.

### Chromosome-naming mismatch
**Trigger:** `1` vs `chr1` between study and panel. **Mechanism:** GRCh38/TOPMed use `chr` prefixes, GRCh37 panels do not, so the join matches nothing. **Symptom:** "0 variants matched" and a wasted day; does not corrupt, just fails. **Fix:** `bcftools annotate --rename-chrs` before any check.

### Expecting a SNP-only or rarity-floored panel to carry a variant
**Trigger:** imputing an indel against HRC, or a variant rarer than the panel floor. **Mechanism:** HRC is SNP-only; its MAC>=5 cutoff means nothing below MAF ~5e-4 is in the panel; un-present variants cannot be imputed at any quality. **Symptom:** the variant is absent or near-zero R2, misread as "imputed poorly." **Fix:** use a panel that contains the variant class (1000G/TOPMed/gnomAD for indels; a WGS panel for rarer variants).

### Server-only panel blocked by governance
**Trigger:** planning to use TOPMed/HRC for data that cannot be uploaded. **Mechanism:** TOPMed is never downloadable and both are server-only; consent/data-residency/IRB rules may forbid uploading participant genotypes to a US server. **Symptom:** the best panel is legally unusable. **Fix:** use a downloadable panel (1000G, HGDP+1kGP) and impute locally.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Drop palindromic (A/T, C/G) SNPs with MAF > 0.4 | Rayner check default | strand unresolvable from alleles and frequency too near 0.5 to disambiguate |
| Allele-frequency concordance flag > 0.2 (stringent 0.1) | Rayner check default | a large study-vs-panel AF gap signals a strand, build, or ancestry problem |
| HRC MAF floor ~5e-4 (MAC>=5 / 32,470) | McCarthy 2016 *Nat Genet* 48:1279 | nothing rarer is in the panel and cannot be imputed |
| Male nonPAR chrX coded haploid | biological ploidy | a het call in male nonPAR is an error and mis-models every male |
| Genetic map must match the panel build | Li & Stephens 2003 *Genetics* 165:2213 | an hg19 map on a GRCh38 panel mis-places recombination silently |
| Match AF comparison to the ancestry-matched sub-panel | Marchini & Howie 2010 *Nat Rev Genet* 11:499 | true frequencies differ by ancestry; a mismatch smears the plot even with perfect strand |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "0 variants matched" the panel | chr naming (`1` vs `chr1`) | `bcftools annotate --rename-chrs` |
| Flipped effect direction at some loci | unresolved palindromic strand | run the harmonization check; drop A/T,C/G MAF>0.4; execute `Run-plink.sh` |
| Indels missing after imputation | HRC is SNP-only | use 1000G/TOPMed/gnomAD |
| AF concordance plot smears off-diagonal | wrong panel ancestry or sample mislabel | compare to the matched sub-panel; check sample labels |
| Cannot download HRC/TOPMed (403) | server-only / access-controlled | use the imputation server, or a downloadable panel |
| Engine errors building bref3 | unphased or missing genotypes in the panel VCF | bref3 needs fully phased, non-missing, `|`-separated input |
| Imputation degraded in one region after liftover | BBIS inverted region strand flip | use a native-build panel; re-check after any liftover |
| Tempted to impute ancestry subgroups of one cohort separately | re-creates the differential-imputation confound (batch-differential quality) | impute all samples together against one large diverse panel (TOPMed or HGDP+1kGP), not per-stratum -> imputation-qc |

## References

- Auton A, Brooks LD, Durbin RM, et al. (1000 Genomes Project Consortium). 2015. A global reference for human genetic variation. *Nature* 526:68-74.
- Byrska-Bishop M, Evani US, Zhao X, et al. 2022. High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. *Cell* 185:3426-3440.
- McCarthy S, Das S, Kretzschmar W, et al. 2016. A reference panel of 64,976 haplotypes for genotype imputation. *Nat Genet* 48:1279-1283.
- Taliun D, Harris DN, Kessler MD, et al. 2021. Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program. *Nature* 590:290-299.
- Koenig Z, Yohannes MT, Nkambule LL, et al. 2024. A harmonized public resource of deeply sequenced diverse human genomes. *Genome Res* 34:796-809.
- Mathias RA, Taub MA, Gignoux CR, et al. 2016. A continuum of admixture in the Western Hemisphere revealed by the African Diaspora genome. *Nat Commun* 7:12522.
- Marchini J, Howie B. 2010. Genotype imputation for genome-wide association studies. *Nat Rev Genet* 11:499-511.
- Das S, Forer L, Schonherr S, et al. 2016. Next-generation genotype imputation service and methods. *Nat Genet* 48:1284-1287.
- Sheng X, Xia L, Cahoon JL, et al. 2023. Inverted genomic regions between reference genome builds in humans impact imputation accuracy and decrease the power of association testing. *HGG Adv* 4:100159.
- Luo Y, Kanai M, Choi W, et al. 2021. A high-resolution HLA reference panel capturing global population diversity enables multi-ancestry fine-mapping in HIV host response. *Nat Genet* 53:1504-1516.

## Related Skills

- haplotype-phasing - The phasing engine that consumes the panel; the genetic-map pairing
- genotype-imputation - Impute untyped variants once the panel is prepared
- imputation-qc - INFO/R2 quality, which cannot detect ancestry mismatch
- variant-calling/variant-normalization - Split multiallelics and left-align before harmonization
- population-genetics/population-structure - PCA to establish target ancestry for panel choice
- clinical-databases/hla-typing - Classical HLA-allele imputation with a dedicated panel
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate
