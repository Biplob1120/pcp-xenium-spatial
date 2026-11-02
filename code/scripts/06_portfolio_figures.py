"""Clean, portfolio-ready figures from results/xenium_niches.h5ad.
Reads only obs + obsm (not the expression matrix), so it is fast and light on memory.
Output: results/figures/portfolio/*.png"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import get_config  # noqa: E402

try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem

cfg = get_config()
out = Path(cfg["out_dir"])
fig_dir = out / "figures" / "portfolio"
fig_dir.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(cfg["seed"])
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 200, "savefig.bbox": "tight"})

# ---------- load obs + coordinates only ----------
with h5py.File(out / "xenium_niches.h5ad", "r") as f:
    obs = read_elem(f["obs"])
    xy = read_elem(f["obsm/spatial"])
    umap = read_elem(f["obsm/X_umap"]) if "X_umap" in f["obsm"] else None
obs = obs.copy()
obs["sample"] = obs["sample"].astype(str)
obs["condition"] = obs["sample"].str.split("_").str[0]          # UN, WT, G30, G60, G90
COND_ORDER = [c for c in ["UN", "WT", "G30", "G60", "G90"] if c in set(obs["condition"])]
COND_NAME = {"UN": "Uninfected", "WT": "WT infected", "G30": "Granuloma d30",
             "G60": "Granuloma d60", "G90": "Granuloma d90"}
SAMPLE_ORDER = sorted(obs["sample"].unique(), key=lambda s: (COND_ORDER.index(s.split("_")[0]), s))
print("cells:", len(obs), "samples:", SAMPLE_ORDER)

# ---------- one fixed colour per cell type (same in every figure) ----------
CT_ORDER = ["AT1", "AT2", "Ciliated", "Endothelial", "Fibroblast", "SMC_pericyte",
            "Alveolar_macrophage", "Mac_Trem2hi_PAM", "Macrophage", "Monocyte", "Myeloid_other",
            "T_CD4", "T_proliferating", "B"]
present = list(obs["cell_type"].astype(str).unique())
CT_ORDER = [c for c in CT_ORDER if c in present] + sorted(c for c in present if c not in CT_ORDER)
PALETTE_MAP = {  # distinct hues: epithelium blues, stroma browns/greys, myeloid greens/reds/oranges, lymphoid purples
    "AT1": "#a6cee3", "AT2": "#1f78b4", "Ciliated": "#17becf", "Endothelial": "#fdbf6f",
    "Fibroblast": "#cab2d6", "SMC_pericyte": "#8c6d31", "Alveolar_macrophage": "#b2df8a",
    "Mac_Trem2hi_PAM": "#33a02c", "Macrophage": "#e31a1c", "Monocyte": "#ff7f00",
    "Myeloid_other": "#000000", "T_CD4": "#6a3d9a", "T_proliferating": "#f781bf", "B": "#b15928"}
EXTRA = ["#7f7f7f", "#bcbd22", "#9edae5", "#c49c94"]
CT_COL = {c: PALETTE_MAP.get(c, EXTRA[i % len(EXTRA)]) for i, c in enumerate(CT_ORDER)}
ct = obs["cell_type"].astype(str).to_numpy()


def ct_legend(ax_or_fig, **kw):
    h = [plt.Line2D([], [], marker="o", ls="", color=CT_COL[c], markersize=8, label=c.replace("_", " "))
         for c in CT_ORDER]
    return ax_or_fig.legend(handles=h, frameon=False, **kw)


# ---------- 1. UMAP on a 150k subsample ----------
if umap is not None:
    n = min(150_000, len(obs))
    idx = rng.choice(len(obs), n, replace=False)
    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5), gridspec_kw={"wspace": 0.55})
    ax = axes[0]
    ax.scatter(*umap[idx].T, c=[CT_COL[c] for c in ct[idx]], s=0.4, lw=0, rasterized=True)
    ax.set_title("Cell types"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ct_legend(ax, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=10)
    ax = axes[1]
    cond_col = dict(zip(COND_ORDER, ["#7f7f7f", "#d62728", "#fdae6b", "#e6550d", "#7f2704"]))
    cc = obs["condition"].to_numpy()[idx]
    ax.scatter(*umap[idx].T, c=[cond_col[c] for c in cc], s=0.4, lw=0, rasterized=True)
    ax.set_title("Condition"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1")
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=cond_col[c], markersize=8, label=COND_NAME[c])
                       for c in COND_ORDER], frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.suptitle(f"Xenium mouse lung, {len(obs):,} cells ({n:,} shown)", y=1.01)
    plt.savefig(fig_dir / "fig1_umap_celltype_condition.png"); plt.close()

# ---------- 2. spatial cell-type maps: one representative section per condition ----------
reps = []
for c in COND_ORDER:
    cand = [s for s in SAMPLE_ORDER if s.startswith(c + "_")]
    reps.append(max(cand, key=lambda s: (obs["sample"] == s).sum()))
fig, axes = plt.subplots(1, len(reps), figsize=(4.6 * len(reps), 5.2))
for ax, s in zip(np.atleast_1d(axes), reps):
    m = (obs["sample"] == s).to_numpy()
    ax.scatter(xy[m, 0], xy[m, 1], c=[CT_COL[c] for c in ct[m]], s=0.3, lw=0, rasterized=True)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    ax.set_title(f"{COND_NAME[s.split('_')[0]]}\n({s})")
    x0, y1 = xy[m, 0].min(), xy[m, 1].max()
    ax.plot([x0, x0 + 1000], [y1 + 150, y1 + 150], color="k", lw=2)
    ax.text(x0 + 500, y1 + 450, "1 mm", ha="center", va="top", fontsize=9)
ct_legend(fig, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=10)
plt.tight_layout()
plt.savefig(fig_dir / "fig2_spatial_celltypes_by_condition.png"); plt.close()

# ---------- 3. zoomed immune view of the largest lesion in a granuloma and WT section ----------
IMMUNE = ["Alveolar_macrophage", "Mac_Trem2hi_PAM", "Macrophage", "Monocyte", "T_CD4", "T_proliferating", "B"]
zoom_samples = [s for s in reps if s.split("_")[0] in ("WT", "G30", "G60", "G90")]
if "lesion_id" in obs and zoom_samples:
    fig, axes = plt.subplots(1, len(zoom_samples), figsize=(4.8 * len(zoom_samples), 5))
    for ax, s in zip(np.atleast_1d(axes), zoom_samples):
        m = (obs["sample"] == s).to_numpy()
        lid = obs["lesion_id"].to_numpy()[m]
        pts = xy[m]
        if (lid >= 0).any():
            big = pd.Series(lid[lid >= 0]).value_counts().index[0]
            cen = pts[lid == big].mean(0)
        else:
            cen = np.median(pts, 0)
        half = 700
        w = (np.abs(pts[:, 0] - cen[0]) < half) & (np.abs(pts[:, 1] - cen[1]) < half)
        sub_ct = ct[m][w]
        imm = np.isin(sub_ct, IMMUNE)
        ax.scatter(*pts[w][~imm].T, c="#e0e0e0", s=1.5, lw=0, rasterized=True)
        ax.scatter(*pts[w][imm].T, c=[CT_COL[c] for c in sub_ct[imm]], s=2.5, lw=0, rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
        ax.set_title(f"{COND_NAME[s.split('_')[0]]} ({s})\nimmune cells, 1.4 mm window")
    h = [plt.Line2D([], [], marker="o", ls="", color=CT_COL[c], markersize=8, label=c.replace("_", " "))
         for c in IMMUNE if c in CT_COL]
    h.append(plt.Line2D([], [], marker="o", ls="", color="#e0e0e0", markersize=8, label="non-immune"))
    fig.legend(handles=h, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(fig_dir / "fig3_lesion_zoom_immune.png"); plt.close()

# ---------- 4. composition: stacked bars per sample ----------
frac = (obs.groupby(["sample", "cell_type"], observed=True).size().unstack(fill_value=0))
frac.columns = frac.columns.astype(str)
frac = frac.reindex(index=SAMPLE_ORDER, columns=CT_ORDER, fill_value=0)
frac = frac.div(frac.sum(axis=1), axis=0) * 100
fig, ax = plt.subplots(figsize=(11, 5))
bottom = np.zeros(len(frac))
for c in CT_ORDER:
    ax.bar(frac.index, frac[c], bottom=bottom, color=CT_COL[c], width=0.85, label=c.replace("_", " "))
    bottom += frac[c].to_numpy()
ax.set_ylabel("% of cells"); ax.set_ylim(0, 100); ax.set_title("Cell-type composition per section")
plt.xticks(rotation=45, ha="right")
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=9)
plt.savefig(fig_dir / "fig4_composition_stacked.png"); plt.close()
frac.round(3).to_csv(out / "portfolio_celltype_percent.tsv", sep="\t")

# ---------- 5. key immune populations by condition (dots = sections) ----------
key = [c for c in ["Monocyte", "Macrophage", "Mac_Trem2hi_PAM", "Alveolar_macrophage", "T_CD4", "B"] if c in frac]
cond_of = pd.Series([s.split("_")[0] for s in frac.index], index=frac.index)
fig, axes = plt.subplots(1, len(key), figsize=(3.1 * len(key), 3.8), sharex=True)
for ax, c in zip(np.atleast_1d(axes), key):
    for i, cond in enumerate(COND_ORDER):
        v = frac.loc[cond_of == cond, c]
        ax.bar(i, v.mean(), color=CT_COL[c], alpha=0.35, width=0.7)
        ax.scatter(np.full(len(v), i) + rng.uniform(-0.12, 0.12, len(v)), v, color=CT_COL[c], edgecolor="k", s=28, zorder=3)
    ax.set_xticks(range(len(COND_ORDER))); ax.set_xticklabels(COND_ORDER, rotation=45)
    ax.set_title(c.replace("_", " "), fontsize=11)
axes[0].set_ylabel("% of cells")
plt.tight_layout()
plt.savefig(fig_dir / "fig5_immune_fraction_by_condition.png"); plt.close()

# ---------- 6. gene programs by condition (whole section mean, dots = sections) ----------
prog = [c for c in obs.columns if c.startswith("prog_")]
if prog:
    pm = obs.groupby("sample", observed=True)[prog].mean().reindex(SAMPLE_ORDER)
    ncol = min(5, len(prog)); nrow = int(np.ceil(len(prog) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.2 * nrow), squeeze=False)
    for ax, c in zip(axes.flat, prog):
        for i, cond in enumerate(COND_ORDER):
            v = pm.loc[[s for s in SAMPLE_ORDER if s.startswith(cond + "_")], c]
            ax.scatter(np.full(len(v), i) + rng.uniform(-0.1, 0.1, len(v)), v, s=30, color="#3182bd", edgecolor="k")
            ax.hlines(v.mean(), i - 0.25, i + 0.25, color="k")
        ax.set_xticks(range(len(COND_ORDER))); ax.set_xticklabels(COND_ORDER, rotation=45)
        ax.set_title(c.replace("prog_", "").replace("_", " "), fontsize=10)
    for ax in list(axes.flat)[len(prog):]:
        ax.axis("off")
    fig.supylabel("mean program score per section")
    plt.tight_layout()
    plt.savefig(fig_dir / "fig6_programs_by_condition.png"); plt.close()

# ---------- 7. cell-type mix vs distance from lesion, per condition ----------
if "lesion_ring" in obs:
    RING_ORDER = ["core", "0-25", "25-50", "50-100", "100-200", "200-400", ">400"]
    inf_conds = [c for c in COND_ORDER if c != "UN"]
    fig, axes = plt.subplots(1, len(inf_conds), figsize=(4.2 * len(inf_conds), 4.6), sharey=True)
    for ax, cond in zip(np.atleast_1d(axes), inf_conds):
        o = obs[(obs["condition"] == cond) & obs["lesion_ring"].notna()]
        t = o.groupby([o["lesion_ring"].astype(str), o["cell_type"].astype(str)]).size().unstack(fill_value=0)
        t = t.reindex(index=[r for r in RING_ORDER if r in t.index], columns=CT_ORDER, fill_value=0)
        n_cells = t.sum(axis=1)
        t = t.div(n_cells.replace(0, np.nan), axis=0).fillna(0) * 100
        bottom = np.zeros(len(t))
        for c in CT_ORDER:
            ax.bar(t.index, t[c], bottom=bottom, color=CT_COL[c], width=0.85)
            bottom += t[c].to_numpy()
        ax.set_title(COND_NAME[cond]); ax.set_ylim(0, 100)
        ax.set_xticks(range(len(t))); ax.set_xticklabels(t.index, rotation=45, ha="right")
        ax.set_xlabel("distance from lesion (µm)")
    np.atleast_1d(axes)[0].set_ylabel("% of cells")
    ct_legend(fig, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=10)
    fig.suptitle("Cell types from lesion core outward", y=1.02)
    plt.tight_layout()
    plt.savefig(fig_dir / "fig7_celltypes_by_lesion_distance.png"); plt.close()

# ---------- 8. niche composition heatmap ----------
f_niche = out / "niche_composition.tsv"
if f_niche.exists():
    nc = pd.read_csv(f_niche, sep="\t", index_col=0)
    nc = nc.reindex(columns=[c for c in CT_ORDER if c in nc.columns])
    counts = obs["niche"].astype(str).value_counts()
    nc.index = [f"Niche {i} (n={counts.get(str(i), 0):,})" for i in nc.index]
    fig, ax = plt.subplots(figsize=(10, 0.45 * len(nc) + 2))
    im = ax.imshow(nc.to_numpy(), cmap="viridis", aspect="auto")
    ax.set_xticks(range(nc.shape[1])); ax.set_xticklabels([c.replace("_", " ") for c in nc.columns], rotation=45, ha="right")
    ax.set_yticks(range(nc.shape[0])); ax.set_yticklabels(nc.index)
    plt.colorbar(im, ax=ax, label="mean neighbour fraction")
    ax.set_title("Spatial niches: which cell types live together")
    plt.savefig(fig_dir / "fig8_niche_composition.png"); plt.close()

print("DONE. Figures in", fig_dir)
for p in sorted(fig_dir.glob("*.png")):
    print(" ", p.name)
