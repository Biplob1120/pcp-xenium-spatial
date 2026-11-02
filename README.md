# Spatial niches of antifungal lung immunity: mapping Pneumocystis-associated macrophage and Th17 programs onto a fungal granuloma Xenium atlas

*Portfolio project (public data only). Modeled on the analysis style of
https://hossainms.github.io/xenium-human-lung-cancer-ffpe/ (QC, then annotation, then spatial niches, then robustness checks).*
Author: Md Biplob Hosen, Kolls Lab, Tulane University (BMS PhD program)

## Scientific question
Do the host programs that single-cell studies of *Pneumocystis* pneumonia (PCP) identified, namely
**Trem2-hi "PCP-associated macrophages" (PAM)**, **Mmp12+ protective macrophages**, and **Th17 / IL-17 / IL-22
CD4 T cells**, form organized spatial niches around fungal lesions in the lung? If so, how do these niches
change as a fungal lesion matures from early infection to a 90-day granuloma?

No public spatial transcriptomics dataset of Pneumocystis-infected lung exists (searched GEO, 10x, Zenodo and
PubMed/bioRxiv, Oct 2026). This project therefore **transfers PCP single-cell programs onto the closest public fungal
lung spatial atlas**, a mouse *Cryptococcus neoformans* granuloma time course profiled with Xenium Prime 5K. The
project also produces a ready-to-use pipeline for future Kolls lab Pneumocystis Xenium data.

## Rationale
- **Pneumocystis and the Kolls lab.** CD4 T cells, Th17 cells and the IL-17/IL-22 axis are central to host defense against
  Pneumocystis and other fungi. Single-cell data place Trem2-hi interstitial macrophages and effector CD4 T cells at the
  center of PCP (GSE157627), and show that glucocorticoids suppress protective Mmp12+ macrophages (GSE225246).
  Single-cell data cannot show *where* these cells sit relative to the fungus. Spatial data can.
- **Translational (CCTS) framing.** PCP is a leading opportunistic infection in people with HIV (a large burden in
  the US South) and in patients on immunosuppression or steroids. Knowing which immune niches surround and contain fungal lesions in the lung could
  point to tissue-level biomarkers and host-directed therapy targets. It also complements the lab's cell-free RNA
  diagnostics work: genes expressed in lesion niches are candidate sources of circulating host RNA.

## Data (verified on GEO, 31 Oct 2026)
| Dataset | What | Samples | Files (per sample) |
|---|---|---|---|
| **GSE328931** (Del Poeta lab, Stony Brook; Nat Commun 2026, PMID 42209491, PMC13392052) | Xenium Prime, **Mouse 5K Pan Tissue & Pathways panel + 46 custom genes = 5,052 genes**, FFPE mouse lung, 5 µm sections, multimodal segmentation (ATP1A1/CD45/E-cadherin, etc.) | 11 sections: uninfected (UN_1, UN_2), WT-infected (WT_1, WT_2), Δgcs1 granuloma day 30 (×3), day 60 (×2), day 90 (×2) | `cell_feature_matrix.h5` (14 to 95 MB), `cells.parquet.gz`, `cell_boundaries`/`nucleus_boundaries.parquet.gz`, `transcripts.parquet.gz` (0.4 to 1.9 GB), `morphology.ome.tif.gz` (4.8 to 10.6 GB) |
| **GSE157627** (Respir Res 2024) | 10x 5′ scRNA-seq, **CD45+ lung cells**, C57BL/6, Pneumocystis weeks 0 to 5 | 6 | `filtered_feature_bc_matrix.tar.gz` (mtx; about 30 to 55 MB) |
| **GSE225246** (Front Immunol 2023) | 10x scRNA-seq, **whole lung** (contains epithelial, endothelial and stromal cells): control, PCP, dexamethasone, dex+PCP (×3 each) | 12 | `<sample>.tar.gz` (mtx; about 85 to 170 MB) |

Totals for GSE328931: matrices 0.6 GB, cells and boundaries about 0.45 GB, transcripts 11.5 GB, morphology images **95 GB**.
The 100 GB `GSE328931_RAW.tar` is not needed because every file can be downloaded per sample.
Cells per section, read from the `.h5` matrices, range from 81k to 323k. The total is about **2.5 million cells**.

**Important format note.** The GEO deposit does **not** include `experiment.xenium` or `gene_panel.json`, so
`spatialdata_io.xenium()` cannot open these folders directly. `scripts/common.py` therefore reads
`cell_feature_matrix.h5` and `cells.parquet` itself, and keeps a `spatialdata_io` path for complete XOA `outs/`
folders (for example the 10x fallback dataset or future lab data).

## Gene-panel coverage (checked against the actual matrix + 10x panel CSV)
Full list: `docs/panel_genes_GSE328931.txt`. Custom add-ons: `docs/custom_addon_genes.txt`. Key genes: `docs/panel_coverage.tsv`.
- **Present:** Trem2, Mmp12, Il17a, Il17f, Il22, Il23r, Rorc, Cd4, Cd8a, Arg1, Nos2, Itgax, Cxcl13, Ccl19, Ccr7,
  Ifng, Cxcl9/10, Gata3, Il13, Foxp3, Clec4n (Dectin-2), Clec4e (Mincle), Cd209a, Il22ra1, Marco, Chil3, Gpnmb,
  Ccr2. The custom add-ons include Il1b, Cd74, H2-Ab1, Dcstamp and sphingolipid genes.
- **Absent:** Spp1, **Clec7a (Dectin-1)**, Siglecf, Il6, Ly6c2, Lyz2, S100a8/9, Sftpc/Sftpb, Scgb1a1, Krt5, Cldn5, Acta2,
  C1qa/b, Nkg7, Trdc, Apoe, Lcn2, Retnla, Ear2. Several of these are very highly expressed genes that the 5K panel excludes.
  This has consequences: neutrophils, AT2 cells, alveolar macrophages and γδ T cells must be called from other markers.
  Dectin-1 biology cannot be measured directly. Every marker list in `scripts/common.py` uses only genes that are on the panel.
- **No fungal probes.** Lesions have to be inferred from host tissue (dense macrophage aggregates). They can be
  checked against H&E or morphology images if needed.

## Aims
1. **Atlas.** QC, cluster and annotate the 11 Xenium sections, and report a fungal lung cell atlas over the infection time course.
2. **PCP label transfer.** Build a combined PCP scRNA-seq reference (GSE157627 + GSE225246). Transfer PAM, Mmp12+
   macrophage and Th17 states onto Xenium cells using genes on the panel (CellTypist, with scANVI as an optional cross-check).
3. **Niches.** Define fungal lesion cores and distance rings. Quantify which PCP-defined states and programs (Th17/IL-17/IL-22,
   IFN-γ, type 2, TLS/iBALT) are enriched at lesion borders compared with distant lung.
4. **Time course and robustness.** Track niche composition at day 30, 60 and 90. Stress-test the results for segmentation
   spillover, label-transfer confidence, leave-one-sample-out stability and panel limits.

## Pipeline
| Step | Script | Output |
|---|---|---|
| 0 | `slurm/00_download_geo.sbatch` (+ `scripts/check_downloads.py`) | per-file download, size/gzip/HDF5 checks |
| 1 | `scripts/01_load_qc.py` | `xenium_qc.h5ad`, QC table, UMAP, lineage dot plot |
| 2 | `scripts/02_reference_scRNA.py` | `pcp_reference.h5ad` (Harmony-integrated, `ref_label`) |
| 3 | `scripts/03_annotate_label_transfer.py` | `xenium_annotated.h5ad` (`lineage`, `cell_type`, `ct_conf`, `prog_*` scores) |
| 4 | `scripts/04_spatial_niches.py` | niches (kNN composition + k-means), lesion cores/rings, squidpy neighbourhood enrichment |
| 5 | `scripts/05_timecourse_stats.py` | per-sample time-course tables, Spearman trend vs day, robustness R1 to R4 |

Run on Cypress:
```bash
module load anaconda            # TODO: confirm module name (`module avail`)
mamba env create -f env/environment.yml   # or conda env create
sbatch --export=ALL,MODE=core,DATA=/lustre/project/<lab>/<you>/pcp_data slurm/00_download_geo.sbatch
# later, for robustness R1 (nucleus-only counts):   MODE=transcripts  (+11.5 GB)
# edit config.yaml paths, then:
for s in 01_load_qc.py 02_reference_scRNA.py 03_annotate_label_transfer.py 04_spatial_niches.py 05_timecourse_stats.py; do
  sbatch --wait --export=ALL,STEP=$s slurm/run_step.sbatch; done
```
The SLURM partition, QOS and module names are placeholders (`TODO`). Check them in the Tulane HPC documentation.

## Expected figures (to make, not results)
1. Study design, panel coverage and QC per section.
2. UMAP and spatial maps of lineages and transferred PCP states in all 11 sections.
3. PCP reference UMAP showing PAM, Mmp12+ and Th17 states, and their transfer confidence on Xenium cells.
4. Lesion cores and distance rings, with program scores (PAM, Th17/IL-17/IL-22, IFN-γ, TLS) plotted against distance from lesion.
5. Neighbourhood-enrichment heat maps and niche composition, for uninfected, WT and granuloma day 30/60/90.
6. Robustness panel: spillover (whole-cell vs nucleus-only), confidence-filtered agreement, leave-one-sample-out, and a list of missing genes.

## Timeline (about 8 to 10 weeks, part-time)
| Week | Work |
|---|---|
| 1 | Set up the environment, download core files (about 2.5 GB), run QC, read the GSE328931 paper |
| 2 to 3 | Build the PCP reference, manually curate labels (PAM, Mmp12+, Th17) |
| 4 | Label transfer, review annotation, optional scANVI |
| 5 to 6 | Lesion detection (tuned on uninfected sections), niches, neighbourhood statistics |
| 7 | Time course and robustness checks; download transcripts for R1; optional Proseg re-segmentation |
| 8 | Figures, GitHub repo and a Quarto/GitHub Pages write-up |
| 9 to 10 | Buffer, discuss with Dr. Kolls, write a CV/fellowship blurb |

## Honest limitations
- **Different fungus and lesion type.** *Cryptococcus* granulomas (made with a Δgcs1 mutant strain) are not Pneumocystis
  alveolar infection. PCP is mainly alveolar and does not form classic granulomas. Results show where PCP-defined programs
  *can* appear in fungal lung lesions; they are not direct evidence about PCP.
- **Mouse and small n.** There are 2 or 3 sections per group (n = 7 granuloma sections). All statistics are per section
  and exploratory. Cells are never treated as replicates.
- **Panel limits.** See above: Dectin-1, Spp1, Il6, Siglecf, Ly6c2, S100a8/9 and Sftpc are not measured. A targeted panel loses
  genes during label transfer, so assignments with low confidence are flagged as `Immune_lowconf`.
- **Reference mismatch.** GSE157627 contains CD45+ cells only. Transferred labels are trusted only on immune lineages.
  Reference labels are made by marker scoring and **must be curated by hand** (TODOs in the scripts).
- **Lesion calls are inferred** from host-cell density, because there are no fungal probes. Validate them against images or the paper's figures.
- The WT-infected timepoint should be confirmed in the paper (marked TODO in `metadata/xenium_samples.tsv`).

## Fallback and extension
- 10x public **Xenium Prime 5K human lung cancer FFPE** dataset (complete `outs/`, 22.6 GB zip):
  `https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/3.0.0/Xenium_Prime_Human_Lung_Cancer_FFPE/Xenium_Prime_Human_Lung_Cancer_FFPE_outs.zip`.
  `common.load_xenium_outs()` can read it.
- The same pipeline can run directly on future Kolls lab Pneumocystis Xenium sections.
