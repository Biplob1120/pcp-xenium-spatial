"""Step 2: build a Pneumocystis (PCP) mouse-lung scRNA-seq reference.
Inputs: GSE157627 (CD45+ lung cells, weeks 0-5, 6 samples) and GSE225246 (lung, control/PCP/
dexamethasone/dex+PCP, 12 samples), both Cell Ranger filtered matrices (mm10).
Output: results/pcp_reference.h5ad with `ref_label` (coarse + PCP-relevant states).
Labels come from marker scoring + YOUR manual review (TODO) - check them against the papers:
  Respir Res 2024 doi:10.1186/s12931-024-02709-1 (Trem2hi PAM) and
  Front Immunol 2023 doi:10.3389/fimmu.2023.1179094 (Mmp12+ macrophages)."""
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
from common import LINEAGE_MARKERS, PROGRAMS, get_config, keep_present, load_scrna_sample  # noqa: E402

cfg = get_config()
out = Path(cfg["out_dir"])
meta = pd.read_csv(cfg["scrna_samples"], sep="\t")

ads = []
for r in meta.itertuples():
    a = load_scrna_sample(cfg["scrna_dir"], r.dataset, r.sample)
    a.obs["condition"], a.obs["week"] = r.condition, r.week
    ads.append(a)
    print(r.dataset, r.sample, a.shape)
ref = ad.concat(ads, join="inner")
ref.var["mt"] = ref.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(ref, qc_vars=["mt"], percent_top=None, inplace=True)
# thresholds follow GSE157627 processing notes (>=200 genes, <25% mito); tune after plotting
ref = ref[(ref.obs.n_genes_by_counts >= 200) & (ref.obs.pct_counts_mt < 25)].copy()
# TODO: doublet removal, e.g. sc.pp.scrublet(ref, batch_key="sample")

ref.layers["counts"] = ref.X.copy()
sc.pp.normalize_total(ref, target_sum=1e4)
sc.pp.log1p(ref)
sc.pp.highly_variable_genes(ref, n_top_genes=3000, batch_key="sample", flavor="seurat")
sc.pp.pca(ref, n_comps=50, mask_var="highly_variable", random_state=cfg["seed"])
try:  # batch correction across samples (harmonypy called directly; output orientation differs by version)
    import harmonypy

    ho = harmonypy.run_harmony(ref.obsm["X_pca"], ref.obs, ["dataset", "sample"], random_state=cfg["seed"])
    Z = np.asarray(ho.Z_corr)
    ref.obsm["X_pca_harmony"] = Z if Z.shape[0] == ref.n_obs else Z.T
    rep = "X_pca_harmony"
except Exception as e:  # harmonypy missing/failing -> fall back, but say so
    print("harmony failed, using uncorrected PCA:", e)
    rep = "X_pca"
sc.pp.neighbors(ref, use_rep=rep, random_state=cfg["seed"])
sc.tl.umap(ref, random_state=cfg["seed"])
sc.tl.leiden(ref, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2,
             random_state=cfg["seed"])

# ---- marker-score based coarse labels ----
for name, genes in {**LINEAGE_MARKERS, **PROGRAMS}.items():
    g = keep_present(genes, ref.var_names)
    if g:
        sc.tl.score_genes(ref, g, score_name=f"score_{name}", random_state=cfg["seed"])
lin_cols = [f"score_{k}" for k in LINEAGE_MARKERS if f"score_{k}" in ref.obs]
cl_lin = ref.obs.groupby("leiden", observed=True)[lin_cols].mean().idxmax(axis=1).str.replace("score_", "")
ref.obs["coarse"] = ref.obs["leiden"].map(cl_lin).astype(str)

# PCP-relevant sub-states inside macrophages / CD4 T cells (cluster-level, thresholded on z-score)
lab = ref.obs["coarse"].copy()
mac = ref.obs["coarse"].isin(["Macrophage", "Alveolar_macrophage", "Monocyte"])
cl = ref.obs.groupby("leiden", observed=True)[["score_PAM_Trem2", "score_Mmp12_mac", "score_Th17_IL17_IL22"]].mean()
z = (cl - cl.mean()) / cl.std()
pam_cl = z.index[z["score_PAM_Trem2"] > 1]
mmp_cl = z.index[z["score_Mmp12_mac"] > 1]
th17_cl = z.index[z["score_Th17_IL17_IL22"] > 1]
lab[mac & ref.obs["leiden"].isin(pam_cl)] = "Mac_Trem2hi_PAM"
lab[mac & ref.obs["leiden"].isin(mmp_cl)] = "Mac_Mmp12"
lab[ref.obs["coarse"].eq("T_CD4") & ref.obs["leiden"].isin(th17_cl)] = "T_CD4_Th17"
ref.obs["ref_label"] = lab.astype("category")
# TODO (important): open results/figures/02_* and manually rename clusters where marker scoring is wrong.
#   e.g. ref.obs["ref_label"] = ref.obs["leiden"].map({...})

ref.obs.groupby(["ref_label", "condition"], observed=True).size().unstack(fill_value=0).to_csv(
    out / "reference_label_counts.tsv", sep="\t")
sc.pl.umap(ref, color=["dataset", "condition", "leiden", "ref_label"], ncols=2, show=False)
plt.savefig(out / "figures/02_reference_umap.png", dpi=150, bbox_inches="tight")
dp = sc.pl.dotplot(ref, {k: keep_present(v, ref.var_names) for k, v in PROGRAMS.items()},
                   groupby="ref_label", standard_scale="var", return_fig=True)
dp.savefig(out / "figures/02_reference_programs.png", dpi=150)
ref.write_h5ad(out / "pcp_reference.h5ad", compression="gzip")
print("saved reference", ref.shape)
