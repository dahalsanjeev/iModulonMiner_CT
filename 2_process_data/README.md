# Process Raw Data
Nextflow pipeline to download and process microbial RNA-seq data from NCBI SRA

## Setup
1. Create the environment with all the requirements with the `nextflow_environment.yaml` file:
    1. `conda env create -f nextflow_environment.yml --name nextflow`
2. Install [Docker](https://docs.docker.com/get-docker/)
3. Download your sequence files for the organism of interest:
    1. Download FASTA and GFF3 files for your genome and plasmids (if relevant) from NCBI.
    2. Put these in a folder named `sequence_files`, and make sure that this folder only contains files for one organism.
    3. Rename the genome files to `genome.fasta` and `genome.gff3`.
    4. Rename plasmid files to `plasmid_<name>.fasta` and `plasmid_<name>.gff3`.
4. [Optional] Update the following fields in `conf/user.conf`. This can also be entered in the command line:
    1. `params.organism`: Name of your organism, including strain information if relevant
    2. `params.metadata`: File path for your metadata file
    3. `params.sequence_dir`: Location of FASTA/GFF3 files

## Modifications to the workflow
- The nextflow workflow has been heavily updated to ensure it runs well on HPC. Furthermore, fastq-dl has replaced fastqc in this iteration of iModulonMiner.


## Run Nextflow on HPC (use run_process_data_imodulonminer_cdm.slurm)
Required Arguments:
  -profile              Executor profile name (e.g. local)
  --organism            Name of organism
  --metadata            Path to metadata file
  --sequence_dir        Directory containing *.fasta and *.gff3 files

Optional Arguments:
  --outdir              Directory to place outputs
  --force               Overwrite existing processed data
  --resume              To resume workflow based on the information in its cache

# Example
```
nextflow run main_fixed_no_delete.nf -profile slurm --organism ctherm --metadata ../ctherm_nf/ctherm_metadata_fixed.tsv --sequence_dir ../ctherm_nf/sequence_files/ --outdir ../ctherm_nf/nf_results
```

## Once it's finished running, you may delete the `work` folder in this root directory to save space.


## Common errors

### Exceeding requirements
If you get the error `Process requirement exceed available CPUs` or `Process requirements exceed available memory` when using `-profile local`, then edit `conf/local.config` and change the CPU and memory requirements to ensure these are within your local computer's parameters.

### Missing R1/R2 columns
If you get the error `Cannot invoke method split() on null object`, this means you are missing the R1 and R2 columns from your metadata file.

## Pipeline Alternatives

### Use of other pipelines
Other alignment pipelines can be used for alignment and quantification of RNA-seq data. Example of this can include [NF-Core](https://nf-co.re/rnaseq/). Alignment of eukaryotic organisms is recommended to be done using alternative pipelines. 
