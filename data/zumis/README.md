# How to process the raw data from NeMO into expression tables

The raw data processing was performed with the processing pipeline zUMIs [1], a flexible pipeline for processing RNA-sequencing data. The latest version of zUMIs can be downloaded from https://github.com/sdparekh/zUMIs.

Before running zUMIs, make sure that you have downloaded the reference genome and a GTF annotation file,  as well as having indexed the reference genome with STAR.

In this paper, we used the mouse mm10 reference genome, STAR version 2.5.4b and the Gencode M23 annotation release.

zUMIs allows a quick and simple way to run the analysis with a miniconda environment that comes with all required dependencies installed. Note that this way of running the pipeline uses a different version of STAR (2.7.3a) compared to what was used in the paper. Consult the zUMIs wiki if you want to run the pipeline using own dependencies. 

The final item needed to run the pipeline is a yaml-file specifying the run parameters. A template yaml file is provided here containing the settings we used to process the data. If using the template yaml file, please make sure to change the file paths so that they suit your setup. The paths that need to be set are paths to the run folder, STAR index, annotation file, path to ERCC fasta file and the path to the output folder. Settings for computer resource usage is also specified in the yaml-file with default settings of 50 CPU cores and RAM limit of 50 Gb.

## Here are step-by-step instructions to creating the count tables:

Step 1. Download the fastq-files from the NeMO repository at http://data.nemoarchive.org/biccn/grant/zeng/tolias/ to the run folder

Step 2. Untar the files using the following command:

```
for i in `ls | grep 'fastq.tar'`; do tar -xvf $i --strip-components 1; done
```

Step 3. Combine the individual fastq-files into one.

```
Rscript zUMIs/misc/merge_demultiplexed_fastq.R --dir path/to/runfolder
```

Note: The zUMIs program uses one fastq-file for all cells as input.  We have provided our raw data as demultiplexed fastq files, one file per cell on the NeMO repository. Before the program can be run, the fastq-files need to be concatenated into one file. Luckily, zUMIs comes with a function to merge fastq-files. When merging, unique barcodes are created for each cell along with a table relating barcodes to cell names for later demultiplexing.

Step 4. Run zUMIs using dependencies from the miniconda environment.

```
zUMIs.sh -c -y run_params_yaml_file.yaml
```

Note: Depending on your system, this might take a few hours. The output will be created in the output directory specified in the yaml file.

## References

* Parekh S, Ziegenhain C, Vieth B, Enard W, Hellmann I. zUMIs - A fast and flexible pipeline to process RNA sequencing data with UMIs. Gigascience. 2018;7(6). doi:10.1093/gigascience/giy059

