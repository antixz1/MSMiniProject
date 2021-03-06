# MSMiniProject

# Introduction
This is a Linux-based python wrapper for a pipeline that takes in an sra accession for E.coli Illumina Reads, assembles the genome, analyzes contigs to predict CDS and their functionalties, and analyzes the transcriptome expression.
# Installation
## External file requirement
This pipeline encorporates a multi-fasta file named Ecoli.fasta to serve as the subject being queried against in a Blastp run. This file is obtained through prokka, but is not automated in this pipeline, so it must be in your home directory to run correctly.
## Python3 and Blast+
Python3 and Blast+ must be installed on your system to properly run.
The file MSMiniProject.py should be downloaded from this repository to run with python3.

## Sratoolkit
Sra toolkit can be downloaded for ubuntu from the [github page](https://github.com/ncbi/sra-tools) by running
```
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
```
The [installation page](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) has info for installing to other machines, but this wrapper depends on Linux software.

The file should be extracted into your home directory using
```
tar -vxzf sratoolkit.tar.gz
```

## SPAdes
SPAdes can be downloaded for Linux from the [github page](https://github.com/ablab/spades) by running
```
 wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz
```
The [installation page](https://github.com/ablab/spades#sec2) has info for installing to other machines, but this wrapper depends on Linux software

The file should be extracted into your home directory using
```
tar -xzf SPAdes-3.15.4-Linux.tar.gz
```

## GeneMarkS-2
GeneMarkS-2 can be downloaded for Linux [here](http://exon.gatech.edu/GeneMark/license_download.cgi).
The first of the two files downloaded should be extracted into your home directory using
```
tar -xf gms2_linux_64.tar.gz
```

The other file is a key that must be in the home directory and named '.gmhmmp2_key' to function correctly. This can be done by running
```
cp gm_key_64 ~/.gmhmmp2_key
```
Note: The key may be hidden in your home directory because it starts with '.' 
## Tuxedo Software
### TopHat2 (2.2.1) can be installed to the home directory with 
```
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar -xf tophat-2.1.1.Linux_x86_64.tar.gz
```

### Cufflinks (2.2.1) can be installed to the home directory with 
```
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -xf cufflinks-2.2.1.Linux_x86_64.tar.gz
```

### Bowtie2 (2.4.5) can be installed with
```
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip/download
unzip bowtie2-2.4.5-linux-x86_64.zip
```
Note: Bowtie2 files need to be downloaded to the same path as TopHat2 to function correctly. This can be done by creating a Path or moving the unzipped contents to tophat:
```
mv ~/bowtie2-2.4.5-linux-x86_64/* ~/tophat-2.2.1.Linux_x86_64
```

# Utility
Note: Example output/input for each pipeline step can be found in the examples folder within this repository.
## Running
The python code automates the use of the installed software to efficiently analyze genomes. This specifically works with E. coli K-12 data, but can be edited for other types of genetic data. After correct installation of prerequisites, the pipeline can be run with
```
python3 MSMiniProject.py
```
All output files will be placed into a results folder, which is created in your home directory automatically. A .log file called miniproject.log can be found in this folder with useful data from the pipeline.
## Data Acquisition 
The pipeline uses sratoolkit and a predefined accession code to obtain sra data from ncbi and translate it to .fastq format. 

## Genome Assembly
SPAdes then assembles the genome and all irrelevent contigs (<1000bp) are parsed out.

## Gene Prediction
GeneMarkS-2 takes the .fasta file generated with SPAdes and predicts coding sequences within the genome. However, GeneMarkS-2 fails to predict what these CDS are. Along with the pre-required multi-fasta Ecoli file from prokka, as well as the predicted CDS, a blastp is run to align and annotate these CDS.
The predicted functionality is then compared against the RefSeq for E. coli K-12.

## Transcriptomics
Sratoolkit is used to obtain the Rna-Seq data for E. coli K-12, but instead of SPAdes assembly, TopHat maps the .fastq formatted data to a Bowtie2 index, automatically created with the complete annotated genome NC_000913. The output is pushed through Cufflinks to quantify this expression.

