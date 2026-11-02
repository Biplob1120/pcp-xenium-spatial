"""Step 4: spatial neighbourhoods / niches around fungal lesions.
- kNN neighbourhood composition per cell -> k-means "niches" (pooled across samples)
- squidpy neighbourhood enrichment + co-occurrence per sample
- lesion definition WITHOUT fungal probes (the panel has no Cryptococcus genes):
  granuloma cores = dense macrophage-rich regions (macrophage + giant/fused-macrophage program),
  then distance of every cell to the nearest lesion core, binned into rings.
  TODO: validate lesion masks against H&E/morphology images (MODE=all download) or the paper's figures.
Output: results/xenium_niches.h5ad, tables + figures."""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN, MiniBatchKMeans

sys.path.insert(0, str(Path(__file__).parent))
from common import get_config  # noqa: E402

cfg = get_config()
out = Path(cfg["out_dir"])
nc = cfg["niche"]
x = sc.read_h5ad(out / "xenium_annotated.h5ad")
types = x.obs["cell_type"].cat.categories
x.obs["is_mac"] = x.obs["lineage"].isin(["Macrophage", "Alveolar_macrophage", "Monocyte"])

comp = np.zeros((x.n_obs, len(types)), dtype=np.float32)
dist_lesion = np.full(x.n_obs, np.nan)
lesion_id = np.full(x.n_obs, -1)
codes = x.obs["cell_type"].cat.codes.to_numpy()
for s in x.obs["sample"].cat.categories:
    idx = np.where(x.obs["sample"].to_numpy() == s)[0]
    xy = x.obsm["spatial"][idx]
    tree = cKDTree(xy)
    _, nn = tree.query(xy, k=nc["n_neighbors"] + 1)
    nn = nn[:, 1:]
    oh = np.eye(len(types), dtype=np.float32)[codes[idx]]
    comp[idx] = oh[nn].mean(1)

    # lesion cores: dense, macrophage-rich neighbourhoods (thresholds in config.yaml -> niche:)
    # NOTE: alveolar macrophages in airspaces give high macrophage fraction at LOW density,
    # hence the density criterion. Tune on uninfected samples (they should give ~no lesions).
    r = nc["lesion_bandwidth_um"]
    nb = tree.query_ball_point(xy, r)
    mac = x.obs["is_mac"].to_numpy()[idx]
    frac_mac = np.array([mac[n].mean() for n in nb])
    dens = np.array([len(n) for n in nb])
    core = (frac_mac >= nc["lesion_min_mac_frac"]) & (dens >= np.quantile(dens, nc["lesion_density_quantile"]))
    if nc.get("lesion_require_giant_program", False) and "prog_Giant_fused_mac" in x.obs:
        giant = x.obs["prog_Giant_fused_mac"].to_numpy()[idx]
        core &= giant >= np.quantile(giant, 0.5)
    if core.sum() >= nc["lesion_min_cells"]:
        lab = DBSCAN(eps=r / 2, min_samples=nc["lesion_min_cells"] // 2).fit_predict(xy[core])
        good = lab >= 0
        if good.any():
            core_idx = np.where(core)[0][good]
            lesion_id[idx[core_idx]] = lab[good]
            d, _ = cKDTree(xy[core_idx]).query(xy)
            dist_lesion[idx] = d
    print(s, "lesion-core cells:", int((lesion_id[idx] >= 0).sum()))

x.obsm["nhood_comp"] = comp
x.obs["lesion_id"] = lesion_id
x.obs["dist_to_lesion_um"] = dist_lesion
bins = [-0.1, 0, 25, 50, 100, 200, 400, np.inf]
x.obs["lesion_ring"] = pd.cut(x.obs["dist_to_lesion_um"], bins,
                              labels=["core", "0-25", "25-50", "50-100", "100-200", "200-400", ">400"])

km = MiniBatchKMeans(n_clusters=nc["n_niches"], random_state=cfg["seed"], n_init=10).fit(comp)
x.obs["niche"] = pd.Categorical(km.labels_.astype(str))
pd.DataFrame(km.cluster_centers_, columns=types).to_csv(out / "niche_composition.tsv", sep="\t")

# per-sample squidpy statistics (graph per sample; coordinates in microns)
for s in x.obs["sample"].cat.categories:
    xs = x[x.obs["sample"] == s].copy()
    sq.gr.spatial_neighbors(xs, coord_type="generic", radius=nc["radius_um"], delaunay=False)
    sq.gr.nhood_enrichment(xs, cluster_key="cell_type", seed=cfg["seed"], show_progress_bar=False)
    pd.DataFrame(xs.uns["cell_type_nhood_enrichment"]["zscore"],
                 index=xs.obs["cell_type"].cat.categories,
                 columns=xs.obs["cell_type"].cat.categories).to_csv(out / f"nhood_enrichment_{s}.tsv", sep="\t")

# program scores vs distance to lesion (core of the Pneumocystis-relevant question)
prog = [c for c in x.obs if c.startswith("prog_")]
ring = x.obs.groupby(["sample", "lesion_ring"], observed=True)[prog].mean()
ring.to_csv(out / "program_by_lesion_ring.tsv", sep="\t")
ring_ct = x.obs.groupby(["sample", "lesion_ring", "cell_type"], observed=True).size().unstack(fill_value=0)
ring_ct.div(ring_ct.sum(1), axis=0).to_csv(out / "celltype_frac_by_lesion_ring.tsv", sep="\t")

# figures: one spatial map per sample coloured by niche / ring
for key in ["niche", "lesion_ring", "cell_type"]:
    S = list(x.obs["sample"].cat.categories)
    fig, axes = plt.subplots(3, 4, figsize=(20, 15))
    for ax, s in zip(axes.flat, S):
        m = (x.obs["sample"] == s).to_numpy()
        cats = x.obs[key].astype("category")
        ax.scatter(*x.obsm["spatial"][m].T, c=cats.cat.codes[m], s=0.2, cmap="tab20", rasterized=True)
        ax.set_title(s); ax.set_aspect("equal"); ax.axis("off")
    for ax in axes.flat[len(S):]:
        ax.axis("off")
    plt.tight_layout(); plt.savefig(out / f"figures/04_spatial_{key}.png", dpi=150); plt.close()

x.write_h5ad(out / "xenium_niches.h5ad", compression="gzip")
print("done")
