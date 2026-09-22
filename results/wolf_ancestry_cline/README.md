# Red wolf, Great Lakes wolf and eastern wolf on the gray-wolf ↔ coyote continuum

**Question.** North America's three "problem" wolves — the Great Lakes wolf, the eastern
(Algonquin) wolf and the red wolf — are each argued to be either distinct lineages or
products of gray-wolf × coyote hybridization. How much of each genome traces to the
gray-wolf lineage versus the coyote lineage, and do they form a west-to-south cline?

![Wolf ancestry cline](wolf_ancestry_cline.png)

## Findings

| Group | Genomes | Gray-wolf ancestry (f4-ratio, GQ ≥ 20) | Same, all hard calls | D(X, western wolf; coyote, jackal) |
|---|---|---|---|---|
| Western gray wolf (control, leave-one-out) | 4 | 0.91–1.05 per genome | 0.98–1.02 | — |
| **Great Lakes wolf** | Wolf18, Wolf40 (Isle Royale) | **74% ± 5%** | 74% ± 4% | 0.081, Z = 3.4 |
| **Eastern wolf, Algonquin only** | 2 | **56% ± 4%** | 65% ± 3% | 0.094, Z = 4.1 |
| Eastern wolf, all three as labelled | + QuebecWolf | 45% ± 4% | 49% ± 2% | 0.129, Z = 6.0 |
| **Red wolf** | Wolf25, Wolf26 | **36% ± 4%** | 36% ± 3% | 0.197, Z = 8.3 |
| Coyote (control, leave-one-out) | 3 | −0.05–0.09 per genome | −0.03–0.04 | — |

(± is one chromosome-jackknife standard error; 38 autosome blocks.)

1. **Red wolves stand apart; Great Lakes and Algonquin wolves overlap.** Point estimates
   run as a cline — ~74% gray wolf in Great Lakes wolves, ~56–65% in Algonquin eastern
   wolves, ~36% in red wolves (64% coyote-lineage) — but only the red-wolf step is decisive
   in direct tests that need no reference model:
   D(red, Algonquin; coyote, jackal) = 0.113 (Z = 4.8) and
   D(red, Great Lakes; coyote, jackal) = 0.134 (Z = 4.8).
   Great Lakes versus Algonquin is **not** resolved: D(Algonquin, Great Lakes; coyote,
   jackal) = 0.023 (Z = 0.9) with two genomes each. (Pooling in the coyote-like Quebec
   genome would have made that step look significant, Z = 2.6 — one reason it is excluded.)
2. **All three share significantly more alleles with coyotes than western gray wolves do**
   (every individual D > 0, Z = 2.2–8.2), so none is simply a gray-wolf population.
3. **Four independent views agree on the order.** The reference-based f4-ratio (panel a/b),
   unsupervised sNMF with no reference groups (panel d), PCA (panel c: Great Lakes →
   Algonquin → red wolf → coyote along PC1) and the D-statistics (panel e) all give the
   same point order, with red wolves clearly the most coyote-like of the three.
4. **The result survives the genotype-quality filter** (panel f). At GQ ≥ 20 the <10×
   genomes lose 50–80% of calls, mostly homozygous-reference; re-running on all hard
   calls moves the red and Great Lakes estimates by <1 point. The Algonquin wolves are the
   most filter-sensitive (56% → 65%), so their value is best read as a range.
5. **The "QuebecWolf" genome is a coyote.** It scores ~0% gray-wolf ancestry, has the
   largest coyote D (Z = 6.1), and clusters with coyotes in PCA and sNMF. Its NCBI record
   (SAMN05736724) already calls it an "apparent gray wolf / coyote hybrid"; in this panel it
   is indistinguishable from coyotes (plausibly an eastern coyote, or a sample mix-up). It is kept in the labelled eastern-wolf group
   for transparency, but the Algonquin-only estimate is the one to quote.
6. **Isle Royale shows its bottleneck.** Among high-coverage genomes (panel g), the Isle
   Royale wolf (Wolf40) has the lowest heterozygosity of any North American wolf
   (0.74 het/kb vs 0.93–1.06 for Yellowstone and Great Lakes 01), consistent with that
   island population's documented inbreeding.

### How this fits the literature

The order and rough magnitudes agree with vonHoldt et al. 2016 (doi:10.1126/sciadv.1501714),
who reported red wolves as roughly three-quarters coyote ancestry, Algonquin wolves as
intermediate, and Great Lakes wolves as predominantly gray wolf with a coyote component.
Our Great Lakes value (74% wolf) sits at the more coyote-rich end of that picture; see the
caveat on what "coyote lineage" means below — any eastern-wolf ancestry in Great Lakes
wolves is counted on the coyote side here.

## What this can and cannot say

- **"Coyote lineage" means the non-gray-wolf side of the tree.** The f4-ratio measures
  drift shared with coyotes versus the gray-wolf stem. It **cannot distinguish recent
  coyote admixture from an old, distinct North American lineage that is sister to
  coyotes** (the *Canis lycaon* / *C. rufus* "distinct species" hypothesis). Both predict
  exactly this cline. Separating them needs divergence-time modeling (e.g. SMC++, momi2,
  or ancestry-tract lengths), which these summary statistics do not provide.
- **Tiny samples.** Two or three genomes per target group, and red wolves from the
  captive breeding program. These are genome-wide estimates for *these individuals*, not
  population-wide parameters.
- **Low coverage.** Five of seven target genomes are 5–7×. The f4-ratio is robust to this
  (panel f), heterozygosity is not, so panel g shows only ≥ 15× genomes.
- **Reference bias.** Reads were mapped to the dog assembly (CanFam3.1), which is closer to
  wolves than to coyotes.
- **Panel-relative heterozygosity.** Heterozygous sites per kb of fetched window, scaled by
  the site density of the joint callset; comparable between samples here, not an absolute
  genome-wide π.
- This is exploratory comparative genomics, not a taxonomic or conservation-management
  determination.

## Data and methods

- **Source:** the NHGRI Dog Genome Project 722-genome joint callset
  (`722g.990.SNP.INDEL.chrAll.vcf.gz`, BioProject PRJNA448733, Plassais et al. 2019,
  doi:10.1038/s41467-019-09373-w), read by indexed byte range only — **2.0 GB transferred
  of a 309 GB file**. Sample identities come from the source's `722g_Metadata_May2018.xlsx`
  and NCBI BioSample records (see [`wolf_cline_samples.csv`](../../configs/examples/wolf_cline_samples.csv)).
  Target genomes: vonHoldt et al. 2016 (PRJNA255370; red wolves Wolf25/26, Great Lakes 01),
  PRJNA299506 (Isle Royale), PRJNA331222 (Algonquin and Quebec).
- **Panel:** 500 contiguous windows (~4 MB compressed each) evenly spaced across the 38
  autosomes. Every PASS biallelic SNP in a window is kept, so the panel is **not
  ascertained** on dog SNP arrays — **75,762 SNPs** variable in this cohort.
- **f4-ratio:** α = f4(Eurasian wolves, golden jackal; X, coyotes) /
  f4(Eurasian wolves, golden jackal; western gray wolves, coyotes) (Patterson et al. 2012).
  Eurasian wolves: Altai, Chukotka, Bryansk, Croatia, Portugal, Spain. Western gray
  wolves: Yellowstone ×3, Alaska. Controls are scored leave-one-out. Validated on msprime
  simulations with known admixture (`tests/unit/test_f_statistics.py`).
- **Uncertainty:** weighted delete-one block jackknife with each autosome as a block.
- **PCA / sNMF:** North American genomes, all hard calls, ≥ 90% call rate, MAF ≥ 5%,
  thinned to one SNP per 2 kb (7,215 SNPs); sNMF K = 2, best of 5 restarts. sNMF tends to
  pull intermediate genomes toward a pole (red wolves 15–18% wolf there vs 36% by f4-ratio),
  which is why the f4-ratio is the headline estimate.

## Files

| File | Contents |
|---|---|
| `wolf_ancestry_cline.png` / `.svg` | The figure above |
| `f4_ratio_wolf_ancestry.csv` | Gray-wolf ancestry per group and genome, both call sets |
| `d_statistics.csv` | D(X, western wolf; coyote, jackal) per group and genome |
| `pca_snmf.csv` | PC1/PC2 and sNMF K = 2 proportions |
| `sample_qc_heterozygosity.csv` | Coverage, missingness under both filters, heterozygosity |
| `summary.json` | Everything above plus the between-target D contrasts (Algonquin-only eastern wolves) |

## Reproduce

```bash
python scripts/fetch_wolf_cline_panel.py --windows 500 --window-bytes 4000000 --threads 6
python scripts/wolf_ancestry_cline.py
```

The fetch writes `data/wolf_cline/panel.npz` (gitignored); the analysis runs in ~15 s.
