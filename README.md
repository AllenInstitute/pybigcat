# pybigcat

Python port of the R package **[scrattch.bigcat](https://github.com/AllenInstitute/scrattch.bigcat)** —
iterative, embedding-based single-cell transcriptomic clustering following the HiCAT scheme
(approximate-KNN → Jaccard/SNN graph → Leiden/Louvain → DE-based merge, applied recursively).

`pybigcat` is the **bigcat** counterpart to
[`transcriptomic_clustering`](https://github.com/AllenInstitute/transcriptomic_clustering) (which mirrors
**scrattch.hicat**): it shares the core HiCAT algorithm but aligns the graph construction, DE scoring, and
merge behavior to `scrattch.bigcat`, and adds the big-data (on-disk / parquet) machinery that has no hicat
equivalent.

## Installation

```bash
pip install -e .        # or add the repo root to PYTHONPATH
```

Dependencies are in `requirements.txt`.

## Quick start (first-step clustering on a precomputed embedding)

```python
import scanpy as sc, numpy as np, pandas as pd
import pybigcat as tc
from pybigcat.iterative_clustering import onestep_clust, OnestepKwargs

adata = sc.read(h5ad_path)                       # normalized (ln CPM) expression
adata.obsm['latent'] = np.asarray(pd.read_csv(latent_csv, index_col=0))   # e.g. scVI

thr = {'q1_thresh':0.5,'q2_thresh':None,'cluster_size_thresh':10,'qdiff_thresh':0.7,
       'padj_thresh':0.01,'lfc_thresh':0.6931472,'score_thresh':300,'low_thresh':0.6931472,
       'min_genes':5,'merge_mode':'aligned'}
kw = OnestepKwargs(
    latent_kwargs={'latent_component':'latent'},
    cluster_louvain_kwargs={'k':15,'annoy_trees':50,'louvain_method':'vtraag',
                            'weighting_method':'jaccard_snn','annoy_seed':1,'resolution':1.0},
    merge_clusters_kwargs={'thresholds':thr,'k':4,'de_method':'ebayes'})
clusters, markers = onestep_clust(adata, kw, random_seed=1)
```

Key bigcat-aligned options: `weighting_method='jaccard_snn'` (SNN graph), `annoy_seed` (seed-invariant
KNN), `merge_mode='aligned'` (R-faithful merge) or `'fast'` (faster, slightly coarser).

---

## What changed vs stock `transcriptomic_clustering`, and their measured effects

Quantified impact of each change on the **first clustering step** (194,221 cells / 32-dim scVI, top-level
`onestep`, no recursion). ARI = Adjusted Rand Index; "self-consist" = same-pipeline seed7-vs-seed99 ARI;
"cross" = R↔Python ARI.

| Change | Test / where | Measured effect |
|---|---|---|
| **SNN graph** (`weighting_method='jaccard_snn'`) vs stock KNN-edge | graph-edge overlap vs R | edge Jaccard **0.31 → 0.90**; weight Pearson **0.23 → 0.988** |
| **DE-score cap = 20/gene** (`calc_de_score`) | DE-score engine test; three-way | with cap, R `fast_limma` ≡ Py `ebayes` (score r=**1.000**, 100% merge-decision agreement). Without cap (stock) → **over-split to 78–84 clusters** vs R's 45–49 |
| **Merge #1 (per-round `best`+extras<th/2) + #3 (`num<min_genes`)** | merge-on-shared vs R | ARI **0.61 → 0.68**; count 50 → 48 (R=42) |
| **Merge #4 (Euclidean candidate metric)** | merge-on-shared vs R | ARI **0.68 → 0.95**; count 48 → **41** (matches R 42) — the dominant merge fix |
| **`merge_mode` `fast` vs `aligned`** | merge-on-shared vs R | fast: 0.906 / 52 cl / ~9 min · aligned: **0.954** / 41 cl / ~30–50 min |
| **Fixed `annoy_seed`** | first-step diff-seed | Python self-consist **0.51 → 0.77** |
| **`annoy_trees` 17 → 50** | first-step diff-seed | Python self-consist **0.51 → 0.57** |
| **Leiden vs Louvain** (both sides) | first-step | self-consist Py 0.77→**0.84**, R 0.91→**0.93**; cross 0.47→**0.52** |
| **`max_cl_size` disabled** | no-cap test | negligible: R self 0.93→0.92, Py 0.84→0.80, cross 0.52→0.52 |
| **ln-converted thresholds** (`lfc`×ln2, `low`×ln2, `padj` 0.01, `min.cells` 10) | rationale (log2-vs-ln); not isolated | no standalone ARI test — corrects log-base scale + over-splitting bias |
| **All changes combined** (`pybigcat` vs stock, same params) | three-way (below) | **cross R↔Py 0.23 → 0.885**; self-consist 0.869 → 0.910 (≈ R's 0.915); counts 78–84 → 47 (R 45–49) |
| **Edge + merge alignment only** (Leiden held) | first-step | **cross 0.52 → 0.88** |

Full methodology and the R-side changes are documented in the alignment study (`changes.md`,
`firststep_comparison.md`, `three_way_comparison.md`).

---

## Three-way validation: R `scrattch.bigcat` vs stock `transcriptomic_clustering` vs `pybigcat`

First-step (top-level `onestep`, no recursion) on the same 194,221 cells / 32-dim scVI embedding, at two
seeds (7, 99):

1. **R `scrattch.bigcat`** — Leiden, no cap, with the `sample_cells` seed fix (reference).
2. **Stock `transcriptomic_clustering`** — pristine upstream (KNN-edge graph, seeded Annoy, no DE cap, stock merge).
3. **`pybigcat`** — all bigcat-aligned changes (SNN graph, fixed Annoy, DE cap, aligned merge).

**The two Python packages use identical tuning parameters** (`k=15`, 50 Annoy trees, Leiden, `resolution=1.0`,
the DE/merge thresholds, merge `k=4`, `de_method='ebayes'`); the differences are therefore **code only** —
which is exactly what this comparison isolates.

### Confusion heatmaps

**Self-consistency — seed 7 vs seed 99 (each pipeline against itself)**
![Three-way self-consistency](figures/3way_selfconsist.png)

**R `scrattch.bigcat` (rows) vs stock `transcriptomic_clustering` (cols)**
![R vs stock](figures/3way_R_vs_stock.png)

**R `scrattch.bigcat` (rows) vs `pybigcat` (cols)**
![R vs pybigcat](figures/3way_R_vs_pybigcat.png)

Heatmaps show the Jaccard overlap between each row-cluster and each column-cluster; columns are ordered to
place best matches near the diagonal. A crisp diagonal = strong agreement.

### Results

| Pipeline | clusters (s7 / s99) | self-consistency (7 vs 99) | vs R (cross, avg) |
|---|---|---|---|
| **R `scrattch.bigcat`** | 45 / 49 | 0.915 | — |
| **stock `transcriptomic_clustering`** | **78 / 84** | 0.869 | **≈ 0.23** |
| **`pybigcat`** | 47 / 47 | 0.910 | **≈ 0.885** |

Same-seed ARI is 1.00 for all three (each is deterministic at a fixed seed).

### Findings

1. **`pybigcat` agrees with R ~4× better than stock does** — R↔Python cross ARI jumps from **≈0.23** to
   **≈0.885**. Matching the graph construction (SNN) and the DE-merge is what closes the gap.
2. **Stock over-splits** to 78–84 clusters (vs R's 45–49); `pybigcat` produces 47, inside R's seed range.
   The over-splitting comes from the KNN-edge graph + missing DE-score cap + stock merge.
3. **`pybigcat` matches R's reproducibility** — self-consistency 0.910 ≈ R's 0.915 (stock 0.869).
4. **The cross now approaches the within-pipeline seed-noise ceiling** — the residual R↔`pybigcat`
   difference is on the order of ordinary seed noise, not a systematic algorithmic gap. A cell-perfect
   match is not attainable (different Annoy forests + different Leiden implementations between R and Python).

### What drives self-consistency (Python)

Python's cross-seed stability is driven mainly by **`annoy_trees=50` + Leiden**, *not* by fixing the
Annoy seed. More trees make the approximate-KNN graph nearly seed-invariant, and Leiden is far more
seed-stable than Louvain; together they take Python different-seed ARI from **0.51 → ~0.91** (matching
R's 0.915). Fixing `annoy_seed` adds only ~0.04 *under Leiden* — it was the dominant lever only under
Louvain (+0.20). Keep `annoy_seed` fixed anyway for a *guaranteed* identical graph across seeds
(reproducibility at ~zero cost), but it is not the main driver.

> **Note on `annoy_trees`:** if not set explicitly it defaults to `max(1, int(log2(n_cells)))` — e.g.
> only **~17** trees for 194k cells (`clustering.py:357–358`), *not* a fixed value. That auto-default is
> the low-tree, less-stable regime, so **set `annoy_trees=50` (or higher) explicitly** for the
> self-consistency above; leaving it unset silently gives the weaker default.

Self-consistency is a **separate axis** from cross-agreement with R: the latter (0.23 → 0.885) required
the SNN graph + DE-merge alignment, not seed handling. Details in `self_consistency.md`.

## Level of support

Research code, updated on no fixed schedule. Issues and pull requests welcome.
