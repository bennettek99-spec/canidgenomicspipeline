# Red wolf, Great Lakes wolf and eastern wolf on the gray-wolf ↔ coyote continuum

**Question.** North America's three "problem" wolves — the Great Lakes wolf, the eastern
(Algonquin) wolf and the red wolf — are each argued to be either distinct lineages or
products of gray-wolf × coyote hybridization. How much of each genome traces to the
gray-wolf lineage versus the coyote lineage, and do they form a west-to-south cline?

![Wolf ancestry cline](wolf_ancestry_cline.png)

## Findings

Allele frequencies are estimated from genotype likelihoods (see *Genotype calling*
below); the hard-call columns are shown for comparison.

| Group | Genomes | Gray-wolf ancestry (f4-ratio, genotype likelihoods) | Hard calls, GQ ≥ 20 | All hard calls | D(X, western wolf; coyote, jackal) |
|---|---|---|---|---|---|
| Western gray wolf (control, leave-one-out) | 4 | 0.98–1.02 per genome | 0.91–1.05 | 0.98–1.02 | — |
| **Great Lakes wolf** | Wolf18, Wolf40 (Isle Royale) | **75% ± 4%** | 74% ± 5% | 74% ± 4% | 0.103, Z = 5.1 |
| **Eastern wolf, Algonquin only** | 2 | **64% ± 3%** | 56% ± 4% | 65% ± 3% | 0.145, Z = 7.4 |
| Eastern wolf, all three as labelled | + QuebecWolf | 49% ± 2% | 45% ± 4% | 49% ± 2% | 0.191, Z = 10.6 |
| **Red wolf** | Wolf25, Wolf26 | **36% ± 3%** | 36% ± 4% | 36% ± 3% | 0.230, Z = 13.4 |
| Coyote (control, leave-one-out) | 3 | −0.02–0.05 per genome | −0.05–0.09 | −0.03–0.04 | — |

(± is one chromosome-jackknife standard error; 38 autosome blocks.)

1. **A three-step cline, with the red-wolf step the largest.** Gray-wolf ancestry falls
   from ~75% in Great Lakes wolves to ~64% in Algonquin eastern wolves to ~36% in red
   wolves (64% coyote-lineage). Direct tests that need no reference model:
   D(red, Algonquin; coyote, jackal) = 0.101 (Z = 5.3),
   D(red, Great Lakes; coyote, jackal) = 0.145 (Z = 6.9),
   D(Algonquin, Great Lakes; coyote, jackal) = 0.053 (Z = 3.0).
   The Algonquin vs Great Lakes step is the weakest, and it is only **suggestive**. It is
   Z = 2.5 on all hard calls and not significant (Z = 0.9) on GQ ≥ 20 calls. Tested one
   genome pair at a time, all four pairs point the same way but reach only Z = 0.8–2.8.
   One of the two Algonquin genomes also looks compromised (see *Follow-up tests*). The
   cleaner Algonquin genome (AlgonquinWolf13467) gives ~61% gray wolf.
2. **All three share significantly more alleles with coyotes than western gray wolves do**
   (every individual D > 0, Z = 3.5–12.5), so none is simply a gray-wolf population.
3. **Four independent views agree on the order.** The reference-based f4-ratio (panel a/b),
   unsupervised sNMF with no reference groups (panel d), PCA (panel c: Great Lakes →
   Algonquin → red wolf → coyote along PC1) and the D-statistics (panel e) all give the
   same order.
4. **The "QuebecWolf" genome is a coyote.** It scores ~3% gray-wolf ancestry, has the
   largest coyote D (Z = 13.6), and clusters with coyotes in PCA and sNMF. Its NCBI record
   (SAMN05736724) already calls it an "apparent gray wolf / coyote hybrid"; in this panel
   it is indistinguishable from coyotes (plausibly an eastern coyote, or a sample mix-up).
   It is kept in the labelled eastern-wolf group for transparency, but the Algonquin-only
   estimate is the one to quote.
5. **Isle Royale shows its bottleneck.** Among high-coverage genomes (panel g), the Isle
   Royale wolf (Wolf40) has the lowest heterozygosity of any North American wolf
   (0.74 het/kb vs 0.93–1.06 for Yellowstone and Great Lakes 01), consistent with that
   island population's documented inbreeding.

## Robustness tests

![Robustness tests](robustness/wolf_cline_robustness.png)

Produced by [`scripts/wolf_cline_robustness.py`](../../scripts/wolf_cline_robustness.py);
tables in [`robustness/`](robustness/). All tests use genotype-likelihood frequencies.

### 1. Is dog ancestry inflating the "gray wolf" estimates? Mostly no.

Dogs sit inside the gray-wolf clade, so the f4-ratio counts any dog ancestry as gray
wolf. qpAdm (Haak et al. 2015) fits each target as western gray wolf + dog (five North
American breeds) + coyote. The outgroups are golden jackal, Eurasian wolves, Asian
wolves, village dogs and dhole.

| Target | Gray wolf | Dog | Coyote | Fit p (3-way) | Fit p (wolf + coyote only) |
|---|---|---|---|---|---|
| Yellowstone wolf (control) | 101% ± 5% | −2% ± 3% | 1% ± 3% | 0.59 | 0.73 |
| Great Lakes wolf | 75% ± 5% | 2% ± 2% | 24% ± 4% | 0.36 | 0.28 |
| Eastern wolf (Algonquin) | 59% ± 3% | **4% ± 2%** (Z = 2.2) | 37% ± 3% | 0.63 | 0.21 |
| Red wolf | 35% ± 4% | 1% ± 2% | 65% ± 3% | 0.89 | 0.95 |

- **Red and Great Lakes wolves:** no detectable dog ancestry. Their gray-wolf estimates
  stand as they are.
- **Algonquin wolves (and the Quebec coyote):** the ~4–6% "dog" weight (individual
  Z = 2.4–3.0) is **not dog-specific**. A model using a Mexican-wolf-like source instead
  of dogs fits just as well, and with both offered the weights become unidentifiable
  (*Follow-up tests*). There is no evidence of dog ancestry here.
- **Mexican wolves (control):** the 10% ± 4% "dog" weight is an artifact. The direct test
  shows no excess dog sharing, which agrees with Fitak et al. 2018 (next section).
- Wolf18 (Great Lakes 01) is the one genome no wolf + dog + coyote model fits
  (p = 0.005–0.03). *Follow-up tests* explains why.

### 2. Where did the coyote-side ancestry come from? A shared eastern component.

Recent local hybridization predicts that each target's coyote component matches the
coyotes living near it. A long-separate lineage predicts it matches no living coyote
better than another. qpAdm was run with each coyote genome as the source and the
other two as outgroups (panel b):

| Target | California coyote as source | Alabama coyote | Midwest coyote |
|---|---|---|---|
| Yellowstone wolf (control) | p = 0.54 | 0.42 | 0.40 |
| Mexican wolves (control) | 0.11 | 0.18 | 0.07 |
| Great Lakes wolf | 0.004 | 0.007 | **0.73** |
| Eastern wolf (Algonquin) | <1e-4 | 1e-4 | **0.51** |
| Red wolf | <1e-4 | <1e-4 | **0.81** |

- **Only the Midwest coyote fits.** For every target and every individual genome
  (including QuebecWolf), the Midwest coyote is the only coyote that works as a source.
  The controls fit with any coyote, so this is not an artifact of that genome's coverage.
- **The Midwest coyote is itself admixed.** Modelled as gray wolf + another coyote, it
  needs 13–18% extra ancestry and still fails (p ≤ 1e-4). The California and Alabama
  coyotes fit each other cleanly (p = 0.22). So the Midwest coyote carries a component
  that is not western-gray-wolf-like, and it shares that component with all three
  eastern wolves.
- **Interpretation.** These wolves don't derive their coyote ancestry *from* Midwest
  coyotes. Rather, there is a **shared eastern North American component**: present in
  Great Lakes wolves, Algonquin wolves, red wolves and at least one Midwest coyote,
  absent from California and Alabama coyotes and from western wolves.
- **Not what simple local hybridization predicts for red wolves.** The Alabama coyote,
  geographically closest to the red wolf's historic range, fits red wolves worst.
  Caveat: coyotes only reached the Southeast in the 20th century, so a modern Alabama
  coyote is an imperfect proxy for the coyotes red wolves met.
- **The simple D-tests don't show this on their own.** The D(coyote i, coyote j; X,
  jackal) tests (panel c) are all non-significant (|Z| ≤ 1.8 for the target groups), and
  if anything they lean slightly toward the California coyote. The qpAdm signal comes
  from combining many f4 statistics, including those against the wolf outgroups. It
  rests on a single 8× Midwest coyote genome, so it needs replication with more coyotes
  before it carries much weight.

### 3. Does the answer depend on which reference genomes were chosen? No.

The f4-ratio was recomputed for 96 combinations per target (panel d):

- **Eurasian reference:** all six Eurasian wolves / Europe only / Asia only.
- **North American wolf reference:** Yellowstone + Alaska / Yellowstone / Alaska / Mexican.
- **Coyote reference:** all three / each coyote alone.
- **Outgroup:** golden jackal / dhole.

Results:

- **Medians:** Great Lakes 75% (range 69–88%), Algonquin 64% (58–74%), red wolf 36%
  (30–44%), QuebecWolf 3% (−6–9%).
- **The ordering holds in all 96 of 96 models.**
- **The one systematic shift:** using Mexican wolves as the North American wolf reference
  raises every target by 4–8 points.

### 4. Follow-up tests: Mexican wolves, Wolf18, and one suspect genome

All values are in [`robustness/followup_tests.csv`](robustness/followup_tests.csv).

**Mexican wolves have no dog ancestry, and they reveal what the "dog" source was absorbing.**

- **No dog signal in the direct test.** Relative to Yellowstone wolves, Mexican wolves
  share no extra alleles with breed dogs (D = +0.004, Z = 0.1) or village dogs (Z = 0.1).
- **What actually differs.** Yellowstone wolves share *more* with Eurasian and Asian wolves
  than Mexican wolves do (Z = 1.9 and 2.1). That fits Mexican wolves as a distinct, older
  North American lineage. In qpAdm, the "dog" source was simply a convenient stand-in for
  "wolf ancestry unlike Yellowstone's".
- **Consistent with the literature.** Fitak et al. 2018 genotyped 87 Mexican wolves at
  >172,000 SNPs and found no biologically significant dog ancestry
  ([doi:10.1093/jhered/esy009](https://doi.org/10.1093/jhered/esy009)). Taron et al.
  2021 found historical and modern Mexican wolves form a discrete unit, distinct from
  other North American wolves ([doi:10.1111/mec.16037](https://doi.org/10.1111/mec.16037)).

**Wolf18 (Great Lakes 01) carries Mexican-wolf-like gray-wolf ancestry.**

- **Leaning toward Mexican wolves.** In a test where coyote ancestry cancels out,
  D(Yellowstone, Wolf18; Mexican wolf, Eurasian wolf) = −0.071 (Z = −2.4): Wolf18 is
  closer to Mexican wolves than Yellowstone wolves are. The Great Lakes (Z = −2.1) and
  Algonquin (Z = −1.9) groups lean the same way; red wolves do not (Z = +0.2).
- **Only a Mexican-like source fixes the model.** Wolf18 is rejected as western wolf +
  coyote (p = 0.005). Adding dog does not help (p = 0.007). Adding Mexican wolves does
  (p = 0.20, ~18% Mexican-like, though that weight is poorly determined at ±35%).
- **Best guess:** Wolf18 is ~75–80% gray wolf and ~20–25% coyote-lineage, as the f4-ratio
  says. Part of its gray-wolf ancestry comes from an older North American wolf lineage
  (today best represented by Mexican wolves) rather than from Yellowstone/Alaska-type
  wolves.
- **Caution.** This rests on one genome and Z ≈ 2 signals; it is a hypothesis to test
  with more Great Lakes genomes, not a result.

**AlgonquinWolf13470 shares an extraordinary amount with the Isle Royale genome.**

- **The anomaly.** D(AlgonquinWolf13467, AlgonquinWolf13470; Wolf40, jackal) = −0.47
  (Z = −29.5). One Algonquin genome shares far more with the Isle Royale wolf than its
  packmate does; qpAdm models it as ~60% "other Algonquin wolf" + ~40% "Isle Royale wolf".
- **Not a sequencing-batch effect.** Other genomes from the same BioProject (the Alaska
  wolf and QuebecWolf) show nothing (|Z| ≤ 1.5).
- **Not a duplicated sample.** Genotype concordance with Wolf40 is 46% at variable sites,
  versus ~90%+ expected for a duplicate.
- **Likely explanation.** The remaining possibilities are contamination of this 7× library
  with DNA from an Isle-Royale-related wolf, a labelling error, or improbably close kinship
  across ~700 km. Treat AlgonquinWolf13470 as suspect.
- **Effect on the results.** It pulls the Algonquin estimate toward Isle Royale's value
  (69% vs 61% for AlgonquinWolf13467), so ~61% is the more trustworthy Algonquin figure.
  It makes Algonquin and Great Lakes wolves look *more* alike, so it does not create the
  cline. It weakens the evidence that they differ.

### Genotype calling

The GQ ≥ 20 filter used in the first version removed 50–80% of calls from the 5–7×
genomes, mostly homozygous-reference ones. Frequencies are now fitted by EM directly from
the genotype likelihoods (Kim et al. 2011), which uses every read without forcing a
hard call. Evidence that this is the better choice:

- The likelihood estimates agree with all-hard-call estimates to within ~1 point
  (panel f of the main figure).
- They also place the controls closest to their true values.
- The GQ ≥ 20 set was the outlier, and it was biased mainly for the Algonquin wolves.

`tests/unit/test_f_statistics.py` checks the EM on simulated 2× reads, where hard calls
are biased and the EM is not.

### How this fits the literature

The order and rough magnitudes agree with vonHoldt et al. 2016 (doi:10.1126/sciadv.1501714).
They reported red wolves as roughly three-quarters coyote ancestry, Algonquin wolves as
intermediate, and Great Lakes wolves as predominantly gray wolf with a coyote component.
Our Great Lakes value (75% wolf) sits at the more coyote-rich end of that picture. Any
eastern-wolf ancestry in Great Lakes wolves is counted on the coyote side here (see the
first caveat below).

## What this can and cannot say

- **"Coyote lineage" means the non-gray-wolf side of the tree.** The f4-ratio measures
  drift shared with coyotes versus the gray-wolf stem. It **cannot distinguish recent
  coyote admixture from an old, distinct North American lineage that is sister to
  coyotes** (the *Canis lycaon* / *C. rufus* "distinct species" hypothesis). The shared
  eastern component in robustness test 2 fits either story: an old eastern lineage whose
  ancestry persists across these canids, or a history of hybridization that has spread
  one gene pool through them. Separating the two needs ancestry-tract lengths or
  divergence-time modeling (SMC++, momi2), not these summary statistics.
- **Tiny samples.** Two or three genomes per target group, a single genome per coyote
  locality, and red wolves from the captive breeding program. These are genome-wide
  estimates for *these individuals*, not population-wide parameters.
- **Low coverage.** Five of seven target genomes are 5–7×. Genotype-likelihood
  frequencies handle this for the f-statistics. Heterozygosity does not, so panel g
  shows only ≥ 15× genomes.
- **Reference bias.** Reads were mapped to the dog assembly (CanFam3.1), which is closer to
  wolves than to coyotes.
- **Panel-relative heterozygosity.** Heterozygous sites per kb of fetched window, scaled by
  the site density of the joint callset. It is comparable between samples here, but it is
  not an absolute genome-wide π.
- This is exploratory comparative genomics, not a taxonomic or conservation-management
  determination.

## Data and methods

- **Source:** the NHGRI Dog Genome Project 722-genome joint callset
  (`722g.990.SNP.INDEL.chrAll.vcf.gz`, BioProject PRJNA448733, Plassais et al. 2019,
  doi:10.1038/s41467-019-09373-w), read by indexed byte range only: **2.0 GB transferred
  of a 309 GB file**. Sample identities come from the source's `722g_Metadata_May2018.xlsx`
  and NCBI BioSample records (see [`wolf_cline_samples.csv`](../../configs/examples/wolf_cline_samples.csv)).
  Target genomes: vonHoldt et al. 2016 (PRJNA255370; red wolves Wolf25/26, Great Lakes 01),
  PRJNA299506 (Isle Royale), PRJNA331222 (Algonquin and Quebec).
- **Panel:** 500 contiguous windows (~4 MB compressed each) evenly spaced across the 38
  autosomes, for 38 genomes: the targets plus coyotes, western/Mexican/Eurasian/Asian gray
  wolves, dogs, village dogs, golden jackal and dhole. Every PASS biallelic SNP in a window
  is kept, so the panel is **not ascertained** on dog SNP arrays: **87,538 SNPs** variable
  in this cohort.
- **f4-ratio:** α = f4(Eurasian wolves, golden jackal; X, coyotes) /
  f4(Eurasian wolves, golden jackal; western gray wolves, coyotes) (Patterson et al. 2012).
  Eurasian wolves: Altai, Chukotka, Bryansk, Croatia, Portugal, Spain. Western gray
  wolves: Yellowstone ×3, Alaska. Controls are scored leave-one-out.
- **qpAdm:** least-squares mixture weights with a block-jackknife χ² test of the model's
  remaining constraints (`canidae.analysis.f_statistics.qpadm`). This is a light
  reimplementation, not ADMIXTOOLS.
- **Validation:** msprime simulations with known two- and three-way admixture confirm that
  the f4-ratio and qpAdm recover the true weights, that correct qpAdm models fit, and that
  a model missing a source is rejected (`tests/unit/test_f_statistics.py`).
- **Uncertainty:** weighted delete-one block jackknife with each autosome as a block.
- **PCA / sNMF:** North American genomes, all hard calls, ≥ 90% call rate, MAF ≥ 5%,
  thinned to one SNP per 2 kb (7,215 SNPs); sNMF K = 2, best of 5 restarts. sNMF tends to
  pull intermediate genomes toward a pole (red wolves 15–18% wolf there vs 36% by
  f4-ratio), which is why the f4-ratio is the headline estimate.

## Files

| File | Contents |
|---|---|
| `wolf_ancestry_cline.png` / `.svg` | Main figure |
| `f4_ratio_wolf_ancestry.csv` | Gray-wolf ancestry per group and genome, all three genotype treatments |
| `d_statistics.csv` | D(X, western wolf; coyote, jackal) per group and genome |
| `pca_snmf.csv` | PC1/PC2 and sNMF K = 2 proportions |
| `sample_qc_heterozygosity.csv` | Coverage, missingness, heterozygosity |
| `summary.json` | Everything above plus the between-target D contrasts |
| `robustness/wolf_cline_robustness.png` / `.svg` | Robustness figure |
| `robustness/qpadm_dog_ancestry.csv` | Three-way and two-way qpAdm fits |
| `robustness/qpadm_coyote_source.csv` | qpAdm with each coyote genome as the source |
| `robustness/qpadm_coyote_as_target.csv` | Each coyote modelled as gray wolf + another coyote |
| `robustness/d_coyote_pairs.csv` | D(coyote i, coyote j; X, jackal) |
| `robustness/reference_rotation.csv` | All 96 f4-ratio models per target |
| `robustness/summary.json` | Summary of the robustness tests |

## Reproduce

```bash
python scripts/fetch_wolf_cline_panel.py --windows 500 --window-bytes 4000000 --threads 6
python scripts/wolf_ancestry_cline.py
python scripts/wolf_cline_robustness.py
```

The fetch writes `data/wolf_cline/panel.npz` (gitignored); each analysis runs in under a
minute.
