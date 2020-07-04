# Phenotypic variation within and across transcriptomic cell types in mouse motor cortex
![patch-seq coverage](cover.png)
https://www.biorxiv.org/content/10.1101/2020.02.03.929158v1

This repository contains the analysis code and the preprocessed data for the above manuscript. Instructions will follow shortly.

------------

## Preprocessed data

All preprocessed data are located in the `data/` folder. This includes meta data, exonic and intronic gene counts, extracted electrophysiological features, extracted morphological features and z-profiles. There are two separate datasests: the main dataset recorded at room temperature and the follow-up dataset recorded at physiological temperature (files with `phys_temp` in the filenames).

## Downloading the raw data

* Raw morphological reconstructions in SWC format can be downloaded from [ftp://download.brainimagelibrary.org:8811/3a/88/3a88a7687ab66069/](https://bit.ly/2VIvP3Z). Download three folders (`excitatory`, `inhibitory`, and `vip_dendrites_only`) and put them into `data/raw/morph/`. These files are needed to create figures that show morphological reconstructions.
