# Phenotypic variation of transcriptomic cell types in mouse motor cortex
![patch-seq coverage](cover.png)
https://www.biorxiv.org/content/10.1101/2020.02.03.929158v1

This repository contains the analysis code and the preprocessed data for the above manuscript. 

------------

## Preprocessed data and meta data

All meta data and preprocessed data are located in the `data/` folder. This includes exonic and intronic gene counts, extracted electrophysiological features, extracted morphological features and z-profiles. There are two separate datasests: the main dataset recorded at room temperature and the follow-up dataset recorded at physiological temperature (files with `phys_temp` in the filenames).

## Downloading the raw data

* Raw morphological reconstructions in SWC format can be downloaded from [ftp://download.brainimagelibrary.org:8811/3a/88/3a88a7687ab66069/](https://bit.ly/2VIvP3Z). Download three folders (`excitatory`, `inhibitory`, and `vip_dendrites_only`) and put them into `data/raw/morph/`. These files are needed to create figures that show morphological reconstructions. Extracted morphological features are provided in this repository (see above).

* Raw electrophysiological traces in NWB format can be downloaded from https://dandiarchive.org/dandiset/000008/. Download all folders (using `dandi` command line tool as described [here](https://github.com/dandi/dandiarchive/issues/385#issuecomment-639803133)] and put them into `data/raw/ephys/`. These files are needed to create figures that illustrate the extraction procedure of electrophysiological features, and to create figures that show electrophysiological traces. Extracted electrophysiological features are provided in this repository (see above). NOTE: Currently the raw data that was prepared for the revision are still missing from the DANDI archive; this will be resolved shortly.

* Raw transcriptomic data in FASTQ format can be downloaded from http://data.nemoarchive.org/biccn/grant/zeng/tolias/. NOTE: the data prepared for revision is currently missing, but should appear shortly.

## The order of running scripts to reproduce our figures

Instructions will follow shortly.
