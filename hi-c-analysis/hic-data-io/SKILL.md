---
name: bio-hi-c-analysis-hic-data-io
description: Loads, converts, and manipulates Hi-C contact matrices in cooler format (.cool/.mcool/.scool) and Juicer .hic, using cooler (Python + CLI), hic2cool, and hictk. Covers the single-resolution mcool URI (file.mcool::/resolutions/<bp>), the load-bearing divisive-vs-multiplicative weight-naming rule (KR/VC/VC_SQRT auto-divisive vs cooler's multiplicative weight), what survives .hic<->.cool conversion (FRAG matrices and norm vectors do not), raw-vs-balanced coarsening, the .pairs upper-triangle/chromsize-order contract, and chrom-naming/bin-table provenance. Use when loading a cooler, converting .hic to .mcool, selecting a resolution, building a cooler from pairs or a matrix, coarsening/zoomifying, importing Juicer norm vectors, or debugging all-NaN balanced matrices and chr1-vs-1 empty fetches.
tool_type: mixed
primary_tool: cooler
---

## Version Compatibility

Reference examples tested with: cooler 0.10+, hic2cool 1.0+, hictk 1.0+, bioframe 0.7+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Two version boundaries change BEHAVIOUR, not just signatures: hic2cool >= 0.5.0 stores Juicer norms un-inverted (divisive); < 0.5.0 inverted them to multiplicative, so two coolers made from one .hic by different hic2cool versions disagree numerically. cooler standardized creation/balance signatures around 0.8-0.10 (`balance_cooler` keyword-only after `clr`; `store=False` default). Record the converter version in provenance.

# Hi-C Data I/O

**"Load my Hi-C contact matrix, convert it, and pull out a region."** -> Open the cooler at one resolution, fetch raw or balanced pixels, and convert across .hic/.cool/.mcool while knowing what does not survive the round trip.
- Python: `cooler.Cooler('file.mcool::/resolutions/10000').matrix(balance=True).fetch('chr1')`
- CLI: `hic2cool convert in.hic out.mcool -r 0` / `hictk convert in.hic out.mcool` / `cooler cload pairs ...`

## The Single Most Important Modern Insight -- cooler Stores Observed Counts; .hic Bakes In Norms, So Conversion Is Never Lossless Both Ways

A `.cool` is a thin, open, HDF5-native *store*: three tables (`chroms`, `bins`, `pixels`) holding raw observed integer counts in an upper-triangle COO layout, plus optional cached `weight` bias columns. Balancing, expected, O/E, eigenvectors -- everything else is computed downstream on demand. Juicer's `.hic` is the opposite philosophy: a sealed binary *deliverable* with all resolutions, precomputed normalization vectors, and expected vectors welded in. Every footgun in this skill descends from that split:

- **Conversion translates between two philosophies, not two encodings.** `.hic -> .cool` loses FRAG (restriction-fragment) matrices (cooler has no FRAG concept) and Juicer's precomputed expected; a norm whose vector is missing arrives as all-NaN. `.cool -> .hic` loses asymmetric matrices, arbitrary extra pixel value columns, and non-Hi-C labeled arrays.
- **The `weight` column NAME is load-bearing.** cooler's own ICE `weight` is applied MULTIPLICATIVELY by `matrix(balance=True)`. Juicer KR/VC norms are DIVISIVE. `Cooler.matrix(divisive_weights=None)` (the default) decides by column name: weights named **KR, VC, or VC_SQRT are auto-treated as divisive**; everything else (including `weight`) is multiplicative. Import a KR vector under the name `weight` and it is applied the wrong way -- garbage, no error.
- **`.mcool` is a container of resolutions, not a matrix.** Every downstream tool wants a single-resolution URI `file.mcool::/resolutions/<bp>`, never the bare `.mcool`. Each resolution is balanced from scratch; the 100kb `weight` is NOT derivable from the 10kb `weight`.

## Format and Tool Taxonomy

| Format / Tool | Role | Mechanism | When |
|---------------|------|-----------|------|
| `.cool` | single-resolution store | HDF5 chroms/bins/pixels, raw counts + optional weight | one resolution; the analysis unit |
| `.mcool` | multi-resolution container | `/resolutions/<bp>/` each a full cooler | HiGlass tilesets; pick a resolution via URI |
| `.scool` | single-cell container | shared bins, `/cells/<name>/pixels` | scHi-C; avoids duplicating the bin table per cell |
| `.hic` (Juicer) | sealed deliverable | binary, all resolutions + baked norm/expected | Juicer/Juicebox ecosystem; BP or FRAG bins |
| `.pairs` (4DN) | upstream contact list | bgzip + pairix index; upper-triangle, flipped | input to `cooler cload`; provenance of the matrix |
| cooler (CLI+Py) | the open standard store/API | pandas/scipy selectors; `cload`/`balance`/`zoomify`/`dump` | the default; cooltools integration |
| hic2cool | .hic -> .cool/.mcool | 4DN-canonical norm handler (>=0.5.0 un-inverted) | importing Juicer norms; BP only |
| hictk | fast cross-format convert/dump | C++; reads .hic v6-9 + cooler, writes .hic v9 + cooler | large files; faster than hic2cool/straw; no FRAG, no asymmetric |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Bare `.mcool`, tool errors / wrong scale | append `::/resolutions/<bp>` URI | the .mcool is a container; tools need one resolution |
| `.hic` -> `.mcool`, want Juicer KR/SCALE preserved | `hic2cool convert -r 0` | 4DN-canonical norm handling; keeps divisive names |
| `.hic` -> cooler, large file, speed matters | `hictk convert` | C++, order-of-magnitude faster; but no FRAG |
| `.hic` was FRAG-binned | re-bin from pairs in BP | FRAG does not survive any converter -> contact-pairs |
| Build a cooler from `.pairs` | `cooler cload pairs -c1 -p1 -c2 -p2 sizes.txt:bp` | needs flipped/deduped pairs -> contact-pairs |
| Need lower resolution | `cooler zoomify --balance` (sum RAW, re-ICE) | cannot sum balanced pixels; re-balance per resolution |
| Imported KR/VC norms | keep their original names | the KR/VC/VC_SQRT auto-divisive rule fires only on those names |
| `matrix(balance=True)` raises / all NaN | balance first; NaN rows = masked bins | unbalanced file has no weight; masked bins are expected NaN |
| `fetch('1')` returns empty on a `chr1` file | harmonize chrom naming first | chr1-vs-1 silently zeros every join, no error |
| Balance/normalize for analysis | -> matrix-operations | ICE/KR mechanics, O/E, expected live there |
| Two coolers, compare pixels | confirm bin tables byte-identical | different contigs/order shift every bin_id |
| scHi-C many cells | `cooler.create_scool` (`.scool`) -> single-cell | shared bin table; per-cell pixels |

## Load a Cooler and Select a Resolution

```python
import cooler

cooler.fileops.list_coolers('matrix.mcool')                      # e.g. ['/resolutions/1000', '/resolutions/10000', ...]
clr = cooler.Cooler('matrix.mcool::/resolutions/10000')          # single-resolution URI, never the bare .mcool
clr.binsize, clr.chromnames, clr.info['sum']                     # info is unvalidated metadata, not a guarantee
'weight' in clr.bins().columns                                   # is this resolution balanced?
```

`clr.matrix().fetch(region)` takes a UCSC region string or a bare chrom name; one arg -> symmetric square, two args -> rectangular submatrix (incl. trans). `clr.bins().fetch('chr1')`, `clr.pixels().fetch(region)` slice the tables.

## Fetch Balanced vs Raw Pixels

**Goal:** Pull a chromosome submatrix as a dense array, correctly balanced.

**Approach:** Confirm the file carries a `weight` column, then `matrix(balance=True)`; on a balanced file, all-NaN rows are masked low-coverage bins (correct), not a bug. For Juicer-imported KR/VC/VC_SQRT weights, name them correctly and the divisive auto-rule fires; force it with `divisive_weights=` only if a custom column is misnamed.

```python
balanced = clr.matrix(balance=True).fetch('chr1')                # multiplicative cooler weight
raw = clr.matrix(balance=False).fetch('chr1')                    # observed counts
kr = clr.matrix(balance='KR').fetch('chr1')                      # KR/VC/VC_SQRT auto-treated as divisive by name
sparse = clr.matrix(balance=True, sparse=True).fetch('chr1')     # scipy COO for large chromosomes
```

## Convert .hic to cooler (and Back)

```bash
hic2cool convert in.hic out.mcool -r 0                           # -r 0 = all resolutions -> .mcool; norm vectors un-inverted (>=0.5.0)
hic2cool convert in.hic out.cool -r 10000                        # single resolution -> .cool
hic2cool extract-norms in.hic out.mcool                          # add Juicer norm vectors to an existing matching cooler
hictk convert in.hic out.mcool                                   # fast C++ path; --resolutions to subset (single -> .cool)
hictk convert in.mcool out.hic                                   # cooler -> .hic v9 ONLY; drops asymmetric/extra columns
```

Neither converter reads/writes FRAG-binned matrices; a FRAG `.hic` yields only its BP resolutions. Record the hic2cool version: the 0.5.0 inversion boundary changes weight values.

## Build a Cooler from a Matrix (Vectorized)

**Goal:** Turn an in-memory numpy contact matrix into a cooler without a hand-rolled O(n^2) loop.

**Approach:** Binnify the chromsizes, take only the upper triangle (cooler stores `symmetric_upper`), pull the nonzero coordinates with a single vectorized `np.triu` + `np.nonzero`, and assemble the pixel DataFrame in one shot.

```python
import cooler
import bioframe
import numpy as np
import pandas as pd

chromsizes = bioframe.fetch_chromsizes('hg38')                   # pin the assembly + chrom naming up front
bins = cooler.binnify(chromsizes, 10000)                         # 10kb bins; bin table identity defines pixel comparability

upper = np.triu(matrix)                                          # cooler stores the upper triangle only
i, j = np.nonzero(upper)                                         # vectorized; never loop over all bin pairs
pixels = pd.DataFrame({'bin1_id': i, 'bin2_id': j, 'count': upper[i, j]})
cooler.create_cooler('new.cool', bins, pixels, assembly='hg38', symmetric_upper=True)
```

For pairs, prefer the CLI: `cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 chromsizes.txt:10000 in.pairs out.cool` (the pairs must already be flipped/deduped -> contact-pairs).

## Coarsen Correctly (Sum Raw, Re-balance Per Resolution)

**Goal:** Produce a lower-resolution or multi-resolution file whose weights are valid.

**Approach:** Coarsen the RAW counts then re-run ICE at each new resolution; `zoomify` does exactly this. Never sum balanced pixels -- weights are resolution-specific and summed balanced values are silently wrong.

```python
cooler.zoomify_cooler('hires.cool', 'out.mcool', resolutions=[10000, 50000, 100000, 500000], chunksize=10_000_000)
cooler.coarsen_cooler('hires.cool', 'lowres.cool', factor=5, chunksize=10_000_000)   # raw sum; re-balance after
```

```bash
cooler zoomify -r 10000,50000,100000,500000 --balance hires.cool   # raw-coarsen then ICE afresh per level
```

## Export, Merge, and Inspect

```python
cooler.merge_coolers('merged.cool', ['rep1.cool', 'rep2.cool'], mergebuf=20_000_000)   # bin tables MUST match
np.save('chr1.npy', clr.matrix(balance=True).fetch('chr1'))
```

```bash
cooler dump -t pixels --join --balanced in.cool > pixels.tsv      # --balanced requires a balanced file
cooler dump -t bins in.cool > bins.tsv
cooler info in.cool ; cooler ls -l in.mcool
```

## Reference, Blacklist, and Chrom-Naming Provenance

The bin table IS the identity of a cooler. Two coolers are pixel-comparable only if their bin tables are byte-identical: same chroms, same order, same contigs present, same binsize. Dropping `chrM`/scaffolds in one pipeline shifts every `bin_id` and makes pixel comparison nonsense. `chr1` (UCSC) vs `1` (Ensembl) silently zeros every cross-tool join and `fetch('1')` on a `chr`-named file returns nothing with no error -- harmonize naming (and the assembly) across the cooler, the genome FASTA, any phasing/annotation track, and any blacklist BED (genome-intervals/bed-file-basics). `info['genome-assembly']` is unvalidated metadata, not a checksum. The `.pairs` `#chromsize` header ORDER defines the upper-triangle convention and the cooler's `assembly`/chromsizes must match the order used to flip the pairs upstream (contact-pairs), or bin assignment and the triangle disagree.

## Multi-Way and Single-Cell Contacts (Decompose or Defer)

Standard Hi-C is pairwise. Multi-way assays (Pore-C, SPRITE, GAM) record higher-order co-occurrence; the common path is to DECOMPOSE concatemers into pairwise contacts and store them in a normal cooler, deferring true multi-way analysis to assay-specific tools rather than forcing it into the COO model. For single-cell Hi-C, `.scool` shares one bin table across `/cells/<name>/pixels` (`cooler.create_scool`); per-cell sparsity and imputation are a distinct world -> single-cell/scatac-analysis for the single-cell chromatin context.

## Per-Method Failure Modes

### Bare .mcool passed where a single resolution is required
**Trigger:** `cooler.Cooler('f.mcool')` or a cooltools call on the bare `.mcool`. **Mechanism:** an `.mcool` is a group of resolutions, not one matrix. **Symptom:** KeyError, wrong/aggregate resolution, or a tool error. **Fix:** use `f.mcool::/resolutions/<bp>`; list with `cooler.fileops.list_coolers`.

### KR/VC norm imported under the name `weight`
**Trigger:** renaming a divisive Juicer norm to `weight`. **Mechanism:** `divisive_weights=None` treats only KR/VC/VC_SQRT as divisive; `weight` is multiplicative. **Symptom:** balanced values are wrong, no error. **Fix:** keep the original KR/VC/VC_SQRT name, or pass `divisive_weights=True` explicitly.

### hic2cool version mismatch across coolers
**Trigger:** two coolers from one `.hic` made by hic2cool <0.5.0 and >=0.5.0. **Mechanism:** pre-0.5.0 inverted norms to multiplicative; >=0.5.0 keeps them divisive. **Symptom:** the same norm gives different balanced values. **Fix:** regenerate both with one version; record the version in provenance.

### Summed balanced pixels to coarsen
**Trigger:** building a coarse balanced matrix by summing finer balanced values. **Mechanism:** balancing weights are resolution-specific. **Symptom:** plausible-looking but wrong coarse values. **Fix:** sum RAW then re-ICE per resolution (`cooler zoomify --balance`).

### all-NaN balanced matrix
**Trigger:** `matrix(balance=True)` raises or returns all NaN. **Mechanism:** an unbalanced file has no `weight`; or those bins were masked during balancing. **Symptom:** error (unbalanced) or NaN rows/cols (masked). **Fix:** balance first (`cooler balance` / `balance_cooler(..., store=True)`); accept masked-bin NaNs as correct.

### FRAG `.hic` "lost resolution" after conversion
**Trigger:** converting a FRAG-binned `.hic`. **Mechanism:** cooler/hictk/hic2cool have no FRAG concept. **Symptom:** only BP resolutions appear; FRAG matrix missing. **Fix:** re-bin in BP from the pairs.

### chrom-naming mismatch
**Trigger:** cooler is `chr1`, a track/blacklist is `1` (or vice versa). **Mechanism:** chrom names never match. **Symptom:** empty `fetch`, zero overlap, no error. **Fix:** harmonize naming across cooler, FASTA, tracks, blacklist.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| hic2cool >= 0.5.0 (un-inverted norms) | hic2cool README | the 0.5.0 boundary flips divisive-vs-multiplicative storage; pin it |
| Cooler weight named KR/VC/VC_SQRT -> divisive | cooler `divisive_weights` rule | only these names auto-trigger divisive; all else multiplicative |
| `ignore_diags=2` (balance default) | `cooler.balance_cooler` default | drop the main + first diagonal (self/near-diagonal artifacts) before ICE |
| `mad_max=5` (balance default) | `cooler.balance_cooler` default | mask bins >5 MAD from the median coverage marginal |
| `min_nnz=10` (balance default) | `cooler.balance_cooler` default | mask sparse bins with <10 nonzero pixels |
| mcool resolution ladder = integer multiples of base | HiGlass tiling | non-integer-multiple levels break tile aggregation; 4DN uses a nice-number series |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `clr.matrix(balance=True)` raises / all NaN | unbalanced file, or masked bins | balance first; masked-bin NaN is expected |
| Empty / wrong-resolution result on `.mcool` | bare `.mcool` passed | use `f.mcool::/resolutions/<bp>` |
| Balanced values look wrong, no error | KR/VC norm renamed to `weight` (treated multiplicative) | keep KR/VC/VC_SQRT name or set `divisive_weights=True` |
| Two coolers from one `.hic` disagree | hic2cool 0.5.0 inversion boundary | regenerate with one version; record it |
| `fetch('1')` returns nothing | chr1-vs-1 naming mismatch | harmonize chrom naming across all inputs |
| FRAG resolutions missing after convert | no FRAG concept in cooler | re-bin from pairs in BP |
| Coarse matrix values wrong | summed balanced pixels | sum raw then re-ICE (`zoomify --balance`) |
| `AttributeError` on a cooler function | pre-0.8/0.10 API change | `help(cooler.<fn>)`; `balance_cooler` is keyword-only after `clr` |

## References

- Abdennur N, Mirny LA. 2020. Cooler: scalable storage for Hi-C data and other genomically labeled arrays. *Bioinformatics* 36(1):311-316.
- Open2C, Abdennur N, Abraham S, Fudenberg G, Flyamer IM, Galitsyna AA, et al. 2024. Cooltools: enabling high-resolution Hi-C analysis in Python. *PLoS Comput Biol* 20(5):e1012067.
- Rossini R, Paulsen J. 2024. hictk: blazing fast toolkit to work with .hic and .cool files. *Bioinformatics* 40(7):btae408.
- Durand NC, Shamim MS, Machol I, Rao SSP, Huntley MH, Lander ES, Aiden EL. 2016. Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. *Cell Syst* 3(1):95-98.
- Imakaev M, Fudenberg G, McCord RP, Naumova N, Goloborodko A, Lajoie BR, Dekker J, Mirny LA. 2012. Iterative correction of Hi-C data reveals hallmarks of chromosome organization. *Nat Methods* 9(10):999-1003.
- Knight PA, Ruiz D. 2013. A fast algorithm for matrix balancing. *IMA J Numer Anal* 33(3):1029-1047.
- Kerpedjiev P, Abdennur N, Lekschas F, et al. 2018. HiGlass: web-based visual exploration and analysis of genome interaction maps. *Genome Biol* 19:125.

## Related Skills

- contact-pairs - Produces the flipped/deduped .pairs that cooler cload bins into a matrix
- matrix-operations - Balancing (ICE/KR), expected, and O/E that operate on the loaded cooler
- hic-visualization - Render the matrices loaded here
- compartment-analysis - Consumes the balanced cooler at compartment resolution
- read-alignment/bwa-alignment - Produces the BAM that pairtools converts to .pairs
- genome-intervals/bed-file-basics - Chrom-naming and BED handling for blacklists/annotation joins
- single-cell/scatac-analysis - Single-cell chromatin context for .scool / scHi-C
