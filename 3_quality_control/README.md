# Step 3: Quality Control

This folder contains two Jupyter notebooks outlining five quality control steps. Between the first and second notebooks, you will need to run curation. For the C. thermocellum, we used data from GOLD and NCBI and GOLD data, available inside the excel file: sra_growth_metadata_ctherm_modified.xlsx. These datasets were modified further in P2 and P3 notebooks.


To run these notebooks, you can either install [Jupyter Notebook](https://jupyter.org/install) or install all requirements listed in the [`environment.yml`](../environment.yml) file on the main page, or use a pre-existing docker container.

## Use of quality control data from alternative alignment pipelines

Alternative pipelines for RNA-seq data processing may provide quality control information in a different format to step 2 of this workflow. Users should consolidate quality control information into the format required by the provided notebooks in order to use these quality control notebooks.
