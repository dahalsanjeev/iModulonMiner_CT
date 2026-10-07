#!/bin/bash
#SBATCH --job-name=iModulonMiner_optICA
#SBATCH --output=slurm_logs/optICA_%j.out
#SBATCH --error=slurm_logs/optICA_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --account=[ACCOUNT_NAME]
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[EMAIL]

# ---------------------------------------------------------------------------
# Create log directory
# ---------------------------------------------------------------------------
mkdir -p slurm_logs

# ---------------------------------------------------------------------------
# Load modules (uncomment / edit for your cluster)
# ---------------------------------------------------------------------------
module purge
# load your anaconda module within HPC
ml anaconda/2023
conda activate modulome_mac


# set path
export GROOVY_HOME=~/.groovy/grapes/
export PATH=$GROOVY_HOME:$PATH
export MYPATH=[PATH_TO_4_optICA]

# RUN

$MYPATH/4_optICA/run_ica.sh -n 8 -min 30 -max 100 -i 100 -s 5 -v -time 21600  -o $MYPATH/3_quality_control/output_qc/interim/ $MYPATH/3_quality_control/output_qc/processed_data/log_tpm_norm.csv

