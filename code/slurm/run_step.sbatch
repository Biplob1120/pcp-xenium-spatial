#!/bin/bash
#SBATCH --job-name=pcp_step
#SBATCH --partition=TODO_PARTITION     # TODO: Cypress partition (a GPU partition only if using scANVI)
#SBATCH --qos=TODO_QOS                 # TODO: remove if not used
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G                     # 11 Xenium samples x 5k genes fit comfortably; raise for transcripts work
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
# Usage: sbatch --export=ALL,STEP=01_load_qc.py slurm/run_step.sbatch
set -euo pipefail
module load anaconda 2>/dev/null || true      # TODO: exact module name on Cypress
source activate pcp_spatial 2>/dev/null || conda activate pcp_spatial
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p logs
python scripts/${STEP:?set STEP=xx_name.py} --config config.yaml
