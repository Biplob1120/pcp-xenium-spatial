#!/bin/bash
#SBATCH --job-name=pcp_dl
#SBATCH --partition=TODO_PARTITION     # TODO: Cypress partition name (check Tulane HPC docs)
#SBATCH --qos=TODO_QOS                 # TODO: remove if your account does not use QOS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#
# Downloads public GEO data for the project, file by file (no 100 GB tar).
#   MODE=core  (default) : Xenium cell_feature_matrix.h5, cells, boundaries (~1 GB total)
#                          + scRNA-seq references (~1.5 GB)
#   MODE=transcripts      : adds transcripts.parquet.gz (~11.5 GB) - needed for spillover/robustness
#   MODE=all              : adds morphology.ome.tif.gz (~95 GB!) - only if you need images
# Usage:  sbatch --export=ALL,MODE=core,DATA=/path/to/data slurm/00_download_geo.sbatch
# GEO publishes no md5 sums for these files; integrity is checked against the byte sizes in
# GEO's filelist.txt, plus gzip -t / HDF5 open tests in scripts/check_downloads.py.
set -euo pipefail
MODE=${MODE:-core}
DATA=${DATA:-$PWD/data}           # TODO: point to Lustre/scratch, not $HOME (morphology is huge)
mkdir -p "$DATA"/xenium "$DATA"/scrna logs
DATA=$(cd "$DATA" && pwd)
PROJ=${SLURM_SUBMIT_DIR:-$PWD}   # submit from the project root
WGET="wget -c -nv --tries=10 --waitretry=30 --retry-on-http-error=429,500,503 --user-agent=Mozilla/5.0"

# --- Xenium GSE328931 (per-sample files) ---
cd "$DATA"/xenium
$WGET -O filelist.txt https://ftp.ncbi.nlm.nih.gov/geo/series/GSE328nnn/GSE328931/suppl/filelist.txt
SUFFIXES="cell_feature_matrix.h5 cells.parquet.gz cell_boundaries.parquet.gz nucleus_boundaries.parquet.gz"
[[ $MODE == transcripts || $MODE == all ]] && SUFFIXES="$SUFFIXES transcripts.parquet.gz"
[[ $MODE == all ]] && SUFFIXES="$SUFFIXES morphology.ome.tif.gz"
while IFS=$'\t' read -r gsm sample rest; do
  [[ $gsm == gsm ]] && continue
  for s in $SUFFIXES; do
    f=${gsm}_${sample}_${s}
    $WGET -O "$f" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9694nnn/${gsm}/suppl/${f}"
    sleep 3   # be polite to NCBI (it rate-limits bursts)
  done
done < "$PROJ/metadata/xenium_samples.tsv"

# --- scRNA-seq references ---
cd "$DATA"/scrna
while IFS=$'\t' read -r gsm ds sample cond week file; do
  [[ $gsm == gsm ]] && continue
  pre=${gsm:0:7}nnn
  $WGET -O "$file" "https://ftp.ncbi.nlm.nih.gov/geo/samples/${pre}/${gsm}/suppl/${file}"
  mkdir -p "$ds/$sample" && tar -xzf "$file" -C "$ds/$sample"
  sleep 3
done < "$PROJ/metadata/scrna_samples.tsv"

cd "$PROJ"
module load anaconda 2>/dev/null || true      # TODO: exact module name on Cypress
source activate pcp_spatial 2>/dev/null || conda activate pcp_spatial
python scripts/check_downloads.py --data "$DATA"
