"""Step 1: load all 11 Xenium samples (GSE328931), QC, normalise, cluster.
Output: results/xenium_qc.h5ad, QC tables and figures. No biological claims here."""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from common import LINEAGE_MARKERS, get_config, keep_present, load_xenium_sample  # noqa: E402

cfg = get_config()
out = Path(cfg["out_dir"])
meta = pd.read_csv(cfg["xenium_samples"], sep="\t")

adatas = []
for r in meta.itertuples():
    a = load_xenium_sample(cfg["xenium_dir"], r.gsm, r.sample)
    a.obs["group"] = r.group
    a.obs["day"] = r.day_post_infection
    adatas.append(a)
    print(f"{r.sample}: {a.n_obs} cells x {a.n_vars} genes")
adata = ad.concat(adatas, join="inner", merge="same")
adata.obsm["spatial"] = np.vstack([a.obsm["spatial"] for a in adatas])
for c in ["sample", "group"]:
    adata.obs[c] = adata.obs[c].astype("category")

# ---- QC metrics ----
sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
neg = adata.obs[["control_probe_counts", "control_codeword_counts", "unassigned_codeword_counts"]].sum(1)
adata.obs["neg_frac"] = neg / adata.obs["total_counts"].clip(lower=1)
q = cfg["qc"]
keep = (
    (adata.obs["total_counts"] >= q["min_counts"])
    & (adata.obs["n_genes_by_counts"] >= q["min_genes"])
    & (adata.obs["neg_frac"] <= q["max_neg_control_frac"])
    & adata.obs["cell_area"].between(q["min_cell_area"], q["max_cell_area"])
)
qc_tab = adata.obs.groupby("sample", observed=True).agg(
    n_cells=("total_counts", "size"),
    median_counts=("total_counts", "median"),
    median_genes=("n_genes_by_counts", "median"),
    median_area=("cell_area", "median"),
    median_neg_frac=("neg_frac", "median"),
)
qc_tab["n_pass"] = keep.groupby(adata.obs["sample"], observed=True).sum()
qc_tab.to_csv(out / "qc_per_sample.tsv", sep="\t")
print(qc_tab)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, k in zip(axes, ["total_counts", "n_genes_by_counts", "cell_area"]):
    adata.obs.boxplot(column=k, by="sample", ax=ax, rot=90, showfliers=False)
    ax.set_yscale("log")
plt.suptitle("")
plt.tight_layout()
plt.savefig(out / "figures/01_qc_boxplots.png", dpi=150)

adata = adata[keep.values].copy()

# ---- normalise / embed / cluster ----
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=None)  # Xenium: scale to median counts
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=50, random_state=cfg["seed"])
# TODO: if batch effects across slides are strong, try harmony (scanpy.external.pp.harmony_integrate)
sc.pp.neighbors(adata, n_neighbors=15, random_state=cfg["seed"])
sc.tl.umap(adata, random_state=cfg["seed"])
for res in (0.3, 0.6, 1.0):
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2,
                 random_state=cfg["seed"])

markers = {k: keep_present(v, adata.var_names) for k, v in LINEAGE_MARKERS.items()}
dp = sc.pl.dotplot(adata, markers, groupby="leiden_0.6", standard_scale="var", return_fig=True)
dp.savefig(out / "figures/01_lineage_dotplot.png", dpi=150)
sc.pl.umap(adata, color=["sample", "group", "leiden_0.6"], ncols=3, show=False)
plt.savefig(out / "figures/01_umap.png", dpi=150, bbox_inches="tight")

adata.write_h5ad(out / "xenium_qc.h5ad", compression="gzip")
print("saved", out / "xenium_qc.h5ad", adata.shape)
