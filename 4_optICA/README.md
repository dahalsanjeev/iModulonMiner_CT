# OptICA
This folder contains all scripts required to run ICA at the optimal dimensionality to identify robust components (optICA). Description and benchmarking of optICA is coming soon.

OptICA may take dozens of hours to run using default arguments, depending on the size of your dataset. This can be accelerated by
1. using more processors (i.e. a supercomputer),
1. loosening the tolerance (e.g. `-t 1e-3`), or
1. increasing the dimensionality step size (e.g. `--step-size 20`).

Also, if your dataset has over 500 datasets, we recommend limiting the maximum dimensionality to the number of unique conditions in your dataset.

The `run_ica.sh` script produces three files and a subdirectory:
- `M.csv`: The **M** matrix
- `A.csv`: The **A** matrix
- `dimension_analysis.pdf`: Plot showing the optimal ICA dimensionality
- `ica_runs/`: A subdirectory containing all the **M** and **A** matrices for all dimensions

## Usage
```bash
Usage: run_ica.sh [ARGS] FILE

Arguments
  -i|--iter <n_iter>	       Number of random restarts (default: 100)
  -t|--tolerance <tol>         Tolerance (default: 1e-7)
  -n|--n-cores <n_cores>       Number of cores to use (default: 8)
  -max|--max-dim <max_dim>     Maximum dimensionality for search (default: n_samples)
  -min|--min-dim <min_dim>     Minimum dimensionality for search (default: 20)
  -s|--step-size <step_size>   Dimensionality step size
  -o|--outdir <path>           Output directory for files (default: current directory)
  -l|--logfile                 Name of log file to use if verbose is off (default: ica.log)
  -v|--verbose                 Send output to stdout rather than writing to file
  -h|--help                    Display help information
  -time|--time-out             Timeout for each ICA run in seconds (default: 7200)
```
## Example Usage
```bash
./run_ica.sh -n 16 -min 100 -max 300 -i 96 -v -time 7200 -o ../data/interim/ ../data/processed_data/log_tpm_norm.csv
```
## Alternatives
**The ICA step is the core step for iModulonMiner and we used FastICA[1] in the workflow. It can be substituted with other ICA implementations like InfoMax ICA[2] or Picard ICA[3], as well as sparse ICA[4]**

[1] Hyvarinen, Aapo. "Fast and robust fixed-point algorithms for independent component analysis." IEEE transactions on Neural Networks 10.3 (1999): 626-634.

[2] Lee, Te-Won, Mark Girolami, and Terrence J. Sejnowski. "Independent component analysis using an extended infomax algorithm for mixed subgaussian and supergaussian sources." Neural computation 11.2 (1999): 417-441.

[3] Ablin, Pierre, Jean-François Cardoso, and Alexandre Gramfort. "Faster independent component analysis by preconditioning with Hessian approximations." IEEE Transactions on Signal Processing 66.15 (2018): 4040-4049.

[4] Wang, Zihang, et al. "Sparse Independent Component Analysis with an Application to Cortical Surface fMRI Data in Autism." Journal of the American Statistical Association (2024): 1-13.