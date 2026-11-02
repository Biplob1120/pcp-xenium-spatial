# Mapping immune niches in *Cryptococcus*-infected lung with Xenium spatial transcriptomics

**Website:** https://biplob1120.github.io/pcp-xenium-spatial/

Md Biplob Hosen · PhD student, Biomedical Sciences · Kolls Lab, Tulane University

A spatial analysis workflow that combines cell-type annotation, neighborhood profiling, and lesion-associated gene programs across infection conditions and granuloma time points. It uses publicly available datasets only.

## Overview
This project develops and evaluates a spatial analysis workflow on a public 10x Xenium Prime 5K dataset of mouse lungs infected with *Cryptococcus neoformans* (GSE328931). Single-cell RNA-seq datasets from *Pneumocystis*-infected lungs (GSE157627, GSE225246) are used as the reference for cell-type annotation, and canonical markers are reviewed alongside them. The longer-term goal is a computational foundation for future spatial studies of *Pneumocystis* infection.

- 2,416,057 cells after QC
- 11 tissue sections
- 5 experimental groups: UN, WT (H99, 15 dpi), and Δ*gcs1* granuloma at G30, G60 and G90
- 5,052 genes on the panel

The full methods, figures, key findings, limitations and references are on the [website](https://biplob1120.github.io/pcp-xenium-spatial/).

## Data
| Dataset | Content | Use |
|---|---|---|
| [GSE328931](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE328931) | Xenium Prime 5K + custom panel, FFPE mouse lung, 11 sections | Spatial query |
| [GSE157627](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157627) | scRNA-seq, CD45+ lung cells, *Pneumocystis* weeks 0–5 | Annotation reference |
| [GSE225246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225246) | scRNA-seq, whole lung, control / *Pneumocystis* / dexamethasone | Annotation reference |

## Repository contents
```
index.html            website
figures/              the eight result figures
code/config.yaml      all analysis parameters (QC, niches, lesion criteria)
code/env/             conda environment (environment.yml)
code/metadata/        sample sheets
code/docs/            panel gene lists and coverage of key genes
code/scripts/         analysis steps (Python, scanpy)
code/slurm/           SLURM job scripts for an HPC cluster
```

## Pipeline
| Step | Script | Purpose |
|---|---|---|
| 0 | `slurm/00_download_geo.sbatch`, `scripts/check_downloads.py` | Download the GEO files and check them |
| 1 | `scripts/01_load_qc.py` | Load the sections, run QC, normalize, cluster |
| 2 | `scripts/02_reference_scRNA.py` | Build the *Pneumocystis* scRNA-seq reference (Harmony) |
| 3 | `scripts/03_annotate_label_transfer.py` | Two-stage annotation (marker lineage + CellTypist), with marker-reviewed corrections |
| 4 | `scripts/04_spatial_niches.py` | Neighborhood clusters, putative lesion regions, distance bands |
| 5 | `scripts/05_timecourse_stats.py` | Section-level summaries and exploratory time trends |
| 6 | `scripts/06_portfolio_figures.py` | Final figures |

## How to run
```bash
conda env create -f code/env/environment.yml
cd code
# set data paths in config.yaml, then run each step as a SLURM job, e.g.
sbatch --export=ALL,STEP=01_load_qc.py slurm/run_step.sbatch
```

## Limitations
The spatial data come from *Cryptococcus*, not *Pneumocystis*. Lesion regions are defined computationally, because there are no fungal probes. The targeted panel limits which genes can be measured. Each group has only 2–3 sections. The results are therefore descriptive; see the website for details.
