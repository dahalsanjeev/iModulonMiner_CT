#!/bin/bash
#SBATCH --job-name=iModulonMiner
#SBATCH --output=slurm_logs/nextflow_%j.out
#SBATCH --error=slurm_logs/nextflow_%j.err
#SBATCH --time=7-00:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=[PARTITION]
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
export MYPATH=[PATH_TO_2_process_data]

## if keeping intermediate files, run this way
nextflow run $MYPATH/2_process_data/main_fixed_no_delete.nf -profile slurm --organism ctherm --metadata $MYPATH/ctherm_nf/ctherm_metadata_fixed.tsv --sequence_dir $MYPATH/ctherm_nf/sequence_files/ --outdir $MYPATH/ctherm_nf/nf_results

## if deleting intermediate files to save space
## Caution: Subsequent runs of iModulonMiner were not tested with this workflow
# nextflow run $MYPATH/2_process_data/main_fixed_w_delete.nf -resume -profile slurm --organism ctherm --metadata $MYPATH/ctherm_nf/ctherm_metadata_fixed.tsv --sequence_dir $MYPATH/ctherm_nf/sequence_files/ --outdir $MYPATH/ctherm_nf/nf_results