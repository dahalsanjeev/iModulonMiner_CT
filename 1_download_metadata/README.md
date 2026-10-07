# Download RNAseq metadata

This folder contains scripts to download and format all RNA-seq metadata for an organism from NCBI Sequence Read Archive. To simplify the process, we have created a docker container with all pre-requisite software.

## Example usage

1. Run download_metadata.
```bash
./download_metadata.sh "Clostridium thermocellum" > ctherm_metadata.tsv
```

2. Fix downloaded metadata using fix_metadata.py
```bash
	python fix_metadata.py ctherm_metadata.tsv -o ctherm_metadata_fixed.tsv
```

- ctherm_metadata_fixed.tsv is used for subsequent steps.
