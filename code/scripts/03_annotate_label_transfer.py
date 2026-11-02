"""Step 3: annotate Xenium cells.
(a) broad lineage for every cell by marker scoring on Xenium clusters (panel-aware markers);
(b) immune sub-states transferred from the PCP scRNA reference with CellTypist, trained only on
    genes shared by the reference and the Xenium panel (~5k genes);
(c) optional scANVI transfer (TODO, GPU) as a second method for the robustness check.
Output: results/xenium_annotated.h5ad"""
import sys
from pathlib import Path

import celltypist
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from common import LINEAGE_MARKERS, PROGRAMS, get_config, keep_present  # noqa: E402

cfg = get_config()
out = Path(cfg["out_dir"])
x = sc.read_h5ad(out / "xenium_qc.h5ad")
ref = sc.read_h5ad(out / "pcp_reference.h5ad")

# ---- (a) broad lineage on Xenium ----
for name, genes in LINEAGE_MARKERS.items():
    g = keep_present(genes, x.var_names)
    if g:
        sc.tl.score_genes(x, g, score_name=f"score_{name}", random_state=cfg["seed"])
cols = [c for c in x.obs if c.startswith("score_")]
cl_lin = x.obs.groupby("leiden_1.0", observed=True)[cols].mean().idxmax(axis=1).str.replace("score_", "")
x.obs["lineage"] = x.obs["leiden_1.0"].map(cl_lin).astype("category")
# TODO: review figures/01_lineage_dotplot.png and correct cluster->lineage mapping by hand.

# ---- (b) CellTypist on shared genes ----
shared = sorted(set(ref.var_names) & set(x.var_names))
print("shared genes reference/panel:", len(shared))
r = ref[:, shared].copy()
r.X = r.layers["counts"].copy()
sc.pp.normalize_total(r, target_sum=1e4)
sc.pp.log1p(r)
model_path = out / "celltypist_pcp_reference.pkl"
if model_path.exists():  # reuse the model trained in an earlier run (training takes ~1 h)
    model = celltypist.models.Model.load(str(model_path))
else:
    model = celltypist.train(r, labels="ref_label", feature_selection=True, top_genes=300,
                             n_jobs=-1, use_SGD=True, mini_batch=True, balance_cell_type=True)
    model.write(str(model_path))

q = x[:, shared].copy()
q.X = q.layers["counts"].copy()
sc.pp.normalize_total(q, target_sum=1e4)
sc.pp.log1p(q)
pred = celltypist.annotate(q, model=model, majority_voting=False)
prob = pred.probability_matrix
x.obs["ct_label"] = pred.predicted_labels["predicted_labels"].values
x.obs["ct_conf"] = prob.max(axis=1).values

# Two-tier annotation: the Xenium lineage (marker based, robust at low counts) decides the broad
# cell class; the reference only chooses a sub-state that is compatible with that lineage.
# Cells whose best compatible sub-state is weak keep their lineage name (no "lowconf" bin).
ALLOWED = {
    "Macrophage": ["Mac_Trem2hi_PAM", "Mac_Mmp12", "Alveolar_macrophage", "Monocyte"],
    "Alveolar_macrophage": ["Alveolar_macrophage", "Mac_Trem2hi_PAM", "Mac_Mmp12"],
    "Monocyte": ["Monocyte", "Mac_Trem2hi_PAM", "Mac_Mmp12"],
    "T_CD4": ["T_CD4", "T_CD4_Th17"],
    "Treg": ["T_CD4", "T_CD4_Th17"],
}
MIN_SUB = cfg.get("annot", {}).get("min_substate_prob", 0.3)
lin = x.obs["lineage"].astype(str).values
cell_type = lin.copy()
sub_conf = np.full(x.n_obs, np.nan)
for L, labs in ALLOWED.items():
    labs = [l for l in labs if l in prob.columns]
    m = lin == L
    if not labs or not m.any():
        continue
    P = prob.loc[m, labs].values
    best = P.argmax(axis=1)
    bestp = P.max(axis=1)
    lab = np.array(labs)[best]
    cell_type[m] = np.where(bestp >= MIN_SUB, lab, L)
    sub_conf[m] = bestp
# Marker-checked relabelling (see figures/03_celltype_marker_dotplot.png from the previous run):
#  - the "Plasma" lineage cluster expresses Cd19/Ms4a1/Cd79a >> Jchain/Mzb1 -> B cells
#  - the reference "T_CD4_Th17" state lands on Mki67/Top2a-high T cells without Rorc/Il17a ->
#    call it proliferating T; Th17 is assessed separately with the Il17a/Il17f/Il22/Rorc program score
#  - Neutrophil lineage has almost no Csf3r/Cxcr2 signal and <400 cells/sample -> merge into Myeloid_other
RELABEL = {"Plasma": "B", "T_CD4_Th17": "T_proliferating", "Neutrophil": "Myeloid_other"}
cell_type = np.array([RELABEL.get(c, c) for c in cell_type])
x.obs["cell_type"] = pd.Categorical(cell_type)
x.obs["substate_conf"] = sub_conf
print("lineage counts:\n", x.obs["lineage"].value_counts())

# ---- program scores (for niches / time course) ----
for name, genes in PROGRAMS.items():
    g = keep_present(genes, x.var_names)
    if g:
        sc.tl.score_genes(x, g, score_name=f"prog_{name}", random_state=cfg["seed"])

# ---- (c) optional scANVI (TODO, run on a GPU partition) ----
# import scvi
# scvi.model.SCVI.setup_anndata(r_counts, layer="counts", batch_key="sample", labels_key="ref_label")
# ... train SCVI -> SCANVI.from_scvi_model(..., unlabeled_category="Unknown") -> load_query_data(x[:, shared])
# store as x.obs["scanvi_label"] and compare with cell_type in 05 (agreement = robustness metric).

x.obs.groupby(["sample", "cell_type"], observed=True).size().unstack(fill_value=0).to_csv(
    out / "celltype_counts_per_sample.tsv", sep="\t")
sc.pl.umap(x, color=["lineage", "cell_type", "substate_conf"], ncols=3, show=False)
plt.savefig(out / "figures/03_annotation_umap.png", dpi=150, bbox_inches="tight")
mk = {k: keep_present(v, x.var_names) for k, v in LINEAGE_MARKERS.items()}
mk = {k: v for k, v in mk.items() if v}
sc.pl.dotplot(x, mk, groupby="cell_type", standard_scale="var", show=False)
plt.savefig(out / "figures/03_celltype_marker_dotplot.png", dpi=150, bbox_inches="tight")
x.write_h5ad(out / "xenium_annotated.h5ad", compression="gzip")
print(x.obs["cell_type"].value_counts())
