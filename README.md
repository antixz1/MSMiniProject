# MichaelMiniProject

# Introduction
This is a Linux-based python wrapper for a pipeline that takes in an sra accession for E.coli Illumina Reads, assemble the genome, analyze contigs to predict CDS and their functionalties, 
# Installation
## Python3
Python3 must be installed to properly run
## Blast+
Blast+ must be installed to properly run
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
