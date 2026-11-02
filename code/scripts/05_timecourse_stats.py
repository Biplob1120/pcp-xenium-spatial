"""Step 5: time course + robustness.
Biological replicate = mouse/section (n = 2-3 per group). Statistics are therefore descriptive /
exploratory: per-sample summaries, Spearman trend over granuloma day (30/60/90; n=7), and
uninfected vs infected contrasts reported as effect sizes. Never test cells as replicates.

Robustness checks (mirroring the template project):
 R1 segmentation spillover: in T cells, fraction of counts from macrophage/epithelial-only genes;
    compare whole-cell vs nucleus-only counts (needs transcripts.parquet, MODE=transcripts download).
 R2 label-transfer confidence: results restricted to ct_conf >= 0.7 should agree with all cells.
 R3 leave-one-sample-out: niche/ring trends recomputed without each sample.
 R4 panel coverage: list key genes absent from the panel (docs/panel_coverage.tsv).
Output: results/timecourse_*.tsv, results/robustness_*.tsv, figures."""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from common import get_config, gunzip_once  # noqa: E402

cfg = get_config()
out = Path(cfg["out_dir"])
x = sc.read_h5ad(out / "xenium_niches.h5ad")
obs = x.obs.copy()
obs["day"] = pd.to_numeric(obs["day"], errors="coerce")

# ---- time course: per-sample composition and programs ----
frac = obs.groupby(["sample", "cell_type"], observed=True).size().unstack(fill_value=0)
frac = frac.div(frac.sum(1), axis=0)
frac.index = frac.index.astype(str)
smeta = obs.groupby("sample", observed=True)[["group", "day"]].first()
smeta.index = smeta.index.astype(str)
frac.join(smeta).to_csv(out / "timecourse_celltype_fraction.tsv", sep="\t")
prog = [c for c in obs if c.startswith("prog_")]
near = obs[obs["lesion_ring"].isin(["core", "0-25", "25-50"])]
if near.empty:
    print("WARNING: no lesion cells detected in any sample - tune niche: thresholds in config.yaml")
p_near = near.groupby("sample", observed=True)[prog].mean().reindex(smeta.index).join(smeta)
p_near.to_csv(out / "timecourse_programs_near_lesion.tsv", sep="\t")

gran = smeta[smeta["group"] == "granuloma"].index
rows = []
for col in list(frac.columns) + prog:
    v = (frac[col] if col in frac else p_near[col]).reindex(gran)
    ok = v.notna()
    rho, p = spearmanr(smeta.loc[gran, "day"][ok], v[ok]) if ok.sum() >= 4 else (np.nan, np.nan)
    rows.append((col, rho, p, len(gran)))
trend = pd.DataFrame(rows, columns=["feature", "spearman_rho_vs_day", "p_raw", "n_samples"])
trend["note"] = "exploratory; n=7 sections, 3 timepoints"
trend.sort_values("p_raw").to_csv(out / "timecourse_trend_granuloma.tsv", sep="\t", index=False)

fig, axes = plt.subplots(1, len(prog), figsize=(3 * len(prog), 3), sharex=True)
for ax, c in zip(np.atleast_1d(axes), prog):  # noqa
    ax.scatter(p_near.loc[gran, "day"], p_near.loc[gran, c])
    ax.set_title(c.replace("prog_", ""), fontsize=8); ax.set_xlabel("day")
plt.tight_layout(); plt.savefig(out / "figures/05_programs_vs_day.png", dpi=150)

# ---- R2: confidence-filter agreement ----
hi = obs[obs["ct_conf"] >= 0.7]
f_hi = hi.groupby(["sample", "cell_type"], observed=True).size().unstack(fill_value=0)
f_hi = f_hi.div(f_hi.sum(1), axis=0).reindex(columns=frac.columns, fill_value=0)
f_hi.index = f_hi.index.astype(str)
r2 = pd.Series({s: spearmanr(frac.loc[s], f_hi.loc[s])[0] for s in f_hi.index}, name="rho_all_vs_highconf")
r2.to_csv(out / "robustness_R2_confidence.tsv", sep="\t")

# ---- R3: leave-one-sample-out for the day trend ----
loo = []
for drop in gran:
    keep = gran.drop(drop)
    for c in prog:
        v = p_near.loc[keep, c]
        ok = v.notna()
        loo.append((drop, c, spearmanr(smeta.loc[keep, "day"][ok], v[ok])[0] if ok.sum() >= 4 else np.nan))
pd.DataFrame(loo, columns=["dropped", "feature", "rho"]).to_csv(out / "robustness_R3_loo.tsv", sep="\t", index=False)

# ---- R1: spillover in T cells (whole cell vs nucleus-only), needs transcripts.parquet ----
foreign = [g for g in ["Cd68", "Adgre1", "Fcgr1", "Mertk", "Marco", "Ager", "Lamp3", "Epcam"] if g in x.var_names]
t_mask = obs["lineage"].isin(["T_CD4", "T_CD8", "Treg"]).to_numpy()
cnt = x.layers["counts"]
fg = np.asarray(cnt[t_mask][:, [x.var_names.get_loc(g) for g in foreign]].sum(1)).ravel()
tot = np.asarray(cnt[t_mask].sum(1)).ravel()
r1 = obs.loc[t_mask, ["sample"]].assign(foreign_frac_wholecell=fg / np.maximum(tot, 1))
meta = pd.read_csv(cfg["xenium_samples"], sep="\t")
for r in meta.itertuples():
    tp = Path(cfg["xenium_dir"]) / f"{r.gsm}_{r.sample}_transcripts.parquet.gz"
    if not tp.exists():
        print("R1 nucleus-only skipped (no transcripts file):", r.sample)
        continue
    t = pd.read_parquet(gunzip_once(tp), columns=["cell_id", "feature_name", "overlaps_nucleus", "qv"])
    t = t[(t["qv"] >= 20) & (t["overlaps_nucleus"] == 1)]
    t["cell_id"] = r.sample + ":" + t["cell_id"].astype(str)
    t["feature_name"] = t["feature_name"].astype(str)
    sub = t[t["cell_id"].isin(r1.index)]
    nuc_tot = sub.groupby("cell_id").size()
    nuc_fg = sub[sub["feature_name"].isin(foreign)].groupby("cell_id").size()
    r1.loc[nuc_tot.index, "foreign_frac_nucleus"] = (nuc_fg.reindex(nuc_tot.index, fill_value=0) / nuc_tot).values
r1.groupby("sample", observed=True).mean(numeric_only=True).to_csv(out / "robustness_R1_spillover.tsv", sep="\t")
# TODO (optional, like the template): re-segment with Proseg or Baysor from transcripts.parquet and
#   re-run 01-04 on the new cells; report how much each headline result changes.

# ---- R4: panel coverage reminder ----
cov = pd.read_csv(Path(__file__).parents[1] / "docs/panel_coverage.tsv", sep="\t")
cov[cov["in_GSE328931_panel"] == "NO"].to_csv(out / "robustness_R4_missing_genes.tsv", sep="\t", index=False)
print("done; see results/*.tsv - interpret as exploratory (small n)")
