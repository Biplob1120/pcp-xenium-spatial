"""Shared helpers: config, Xenium loading from GEO per-sample files, scRNA loading.

GEO GSE328931 provides per-sample XOA outputs *without* experiment.xenium / gene_panel.json,
so spatialdata_io.xenium() cannot read them directly. `load_xenium_sample` therefore reads
cell_feature_matrix.h5 + cells.parquet(.gz) itself. If you ever have a full XOA `outs/`
folder (e.g. from the authors or the 10x fallback dataset), `load_xenium_outs` uses
spatialdata_io.xenium().
"""
from __future__ import annotations

import argparse
import gzip
import io
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

XENIUM_PIXEL_SIZE_UM = 0.2125  # morphology image pixel size for Xenium (check per run if using images)


def get_config() -> dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    Path(cfg["out_dir"]).mkdir(parents=True, exist_ok=True)
    (Path(cfg["out_dir"]) / "figures").mkdir(exist_ok=True)
    return cfg


def read_parquet_maybe_gz(path: Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as fh:
            return pd.read_parquet(io.BytesIO(fh.read()))
    return pd.read_parquet(path)


def gunzip_once(path: Path) -> Path:
    """Decompress a large .gz (e.g. transcripts.parquet.gz) to disk once; return the plain path."""
    path = Path(path)
    out = path.with_suffix("")
    if not out.exists():
        with gzip.open(path, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst, length=64 << 20)
    return out


def load_xenium_sample(xenium_dir: str | Path, gsm: str, sample: str):
    """Return AnnData (cells x genes, raw counts) with obsm['spatial'] in microns."""
    import scanpy as sc

    d = Path(xenium_dir)
    prefix = f"{gsm}_{sample}_"
    adata = sc.read_10x_h5(d / f"{prefix}cell_feature_matrix.h5", gex_only=True)
    adata.var_names_make_unique()
    cells_path = d / f"{prefix}cells.parquet.gz"
    if not cells_path.exists():
        cells_path = d / f"{prefix}cells.parquet"
    cells = read_parquet_maybe_gz(cells_path)
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id")
    missing = adata.obs_names.difference(cells.index)
    if len(missing):
        raise ValueError(f"{sample}: {len(missing)} matrix barcodes absent from cells.parquet")
    cells = cells.loc[adata.obs_names]
    adata.obs = cells.drop(columns=[c for c in ["x_centroid", "y_centroid"] if c in cells]).copy()
    adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    adata.obs["sample"] = sample
    adata.obs["gsm"] = gsm
    adata.obs_names = [f"{sample}:{c}" for c in adata.obs_names]
    return adata


def load_xenium_outs(outs_dir: str | Path):
    """Full XOA outs folder -> SpatialData via spatialdata_io (needs experiment.xenium)."""
    from spatialdata_io import xenium

    return xenium(str(outs_dir), cells_as_circles=True)


def find_mtx_dirs(root: Path) -> list[Path]:
    return sorted({p.parent for p in Path(root).rglob("matrix.mtx*")})


def load_scrna_sample(scrna_dir: str | Path, dataset: str, sample: str):
    import scanpy as sc

    root = Path(scrna_dir) / dataset / sample
    dirs = find_mtx_dirs(root)
    if len(dirs) != 1:
        raise FileNotFoundError(f"expected one matrix.mtx under {root}, found {dirs}")
    ad = sc.read_10x_mtx(dirs[0], var_names="gene_symbols", cache=False)
    ad.var_names_make_unique()
    ad.obs_names = [f"{dataset}:{sample}:{b}" for b in ad.obs_names]
    ad.obs["dataset"], ad.obs["sample"] = dataset, sample
    return ad


# Broad lineage markers verified present in the GSE328931 panel (5K mouse + 46 custom genes).
# Absent from this panel (do not use): Sftpc, Sftpb, Scgb1a1, Krt5, Cldn5, Acta2, C1qa/b,
# Siglecf, Ly6c2, Lyz2, S100a8/9, Nkg7, Trdc, Apoe, Spp1, Clec7a, Il6, Lcn2, Retnla, Ear2.
LINEAGE_MARKERS = {
    "AT1": ["Ager", "Hopx"],
    "AT2": ["Lamp3", "Abca3", "Etv5"],
    "Ciliated": ["Foxj1"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr"],
    "Fibroblast": ["Col1a1", "Col1a2", "Pdgfra"],
    "SMC_pericyte": ["Myh11", "Pdgfrb"],
    "Macrophage": ["Cd68", "Adgre1", "Fcgr1", "Mertk"],
    "Alveolar_macrophage": ["Marco", "Pparg", "Chil3", "Itgax"],
    "Monocyte": ["Ccr2", "Cd14", "Itgam"],
    "Neutrophil": ["Csf3r", "Cxcr2", "Ly6g"],
    "DC": ["Flt3", "Xcr1", "Clec9a", "Ccr7"],
    "T_CD4": ["Cd3e", "Cd3d", "Cd4"],
    "T_CD8": ["Cd3e", "Cd8a", "Cd8b1"],
    "Treg": ["Foxp3"],
    "NK": ["Ncr1", "Klrb1c", "Gzma"],
    "B": ["Cd19", "Ms4a1", "Cd79a"],
    "Plasma": ["Jchain", "Mzb1"],
    "Proliferating": ["Mki67", "Top2a"],
}

# Programs of interest (all verified present in the panel).
PROGRAMS = {
    "PAM_Trem2": ["Trem2", "Gpnmb", "Lpl", "Cd9", "Ctsk"],
    "Mmp12_mac": ["Mmp12", "Mmp9"],
    "Giant_fused_mac": ["Dcstamp", "Ocstamp", "Acp5", "Ctsk"],
    "Th17_IL17_IL22": ["Rorc", "Il17a", "Il17f", "Il22", "Il23r"],
    "Th1_IFNg": ["Ifng", "Cxcl9", "Cxcl10", "Gbp2", "Irf1"],
    "Type2": ["Gata3", "Il5", "Il13", "Il1rl1", "Arg1"],
    "M1_like": ["Nos2", "Il1b", "Tnf", "Cxcl2"],
    "TLS_iBALT": ["Cxcl13", "Ccl19", "Ccr7", "Ms4a1"],
    "CLR_fungal_sensing": ["Clec4n", "Clec4e", "Cd209a"],
}


def keep_present(genes, var_names):
    s = set(var_names)
    return [g for g in genes if g in s]
