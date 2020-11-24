#  rapvis: a tool for RNAseq processing and visualization
  
  
##  Dependency
  
  
Required python version:
  
+ python >= 3.6
  
Sevral external software were depended for rapvis:
  
+ trimmomatic
+ hisat2
+ stringtie
+ bwa
+ samtools
+ STAR
+ bedtools
  
##  Mandatory
  
  
+ pandas >= 1.1.2
+ numpy
+ matplotlib
+ seaborn
+ GSEApy
+ rpy2
  
##  Install
  
  
###  Installing from github
  
  
```bash
# Clone remote repository
$ git clone https://github.com/liuwell/rapvis.git
  
# Install required python pacakge
$ cd rapvis
$ pip install -r requirements.txt
  
# Add execution path
# The path of current dir can get by shell command "pwd"
$ echo "export PATH=$PATH:current_dir/rapvis" >> ~/.bashrc
$ source ~/.bashrc
```
  
```bash{code_chunk_offset=0,
# Then you can type -h option to check whether the installation is successful,  
# If the output as follows, it means your installation is successful
$ rapvis_run.py -h
```

```
usage: rapvis_run.py [-h] -i INPUT [-o OUTPUT] [-p THREADS] [-lib LIBRARYPATH]
                     [-m {hisat2,STAR}] [-a {nextera,universal}]
                     [--minlen MINLEN] [--trim5 TRIM5] [--rRNA] [-v]

A tool for RNAseq processing and visualization

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the input data
  -o OUTPUT, --output OUTPUT
                        output directory (default: processed_data)
  -p THREADS, --threads THREADS
                        number of threads (CPUs) to use (default: 5)
  -lib LIBRARYPATH, --libraryPath LIBRARYPATH
                        choose reference species for mapping and annotaion
  -m {hisat2,STAR}, --mapper {hisat2,STAR}
                        choose the mapping program (default: hisat2)
  -a {nextera,universal}, --adapter {nextera,universal}
                        choose illumina adaptor (default: nextera)
  --minlen MINLEN       discard reads shorter than LEN (default: 35)
  --trim5 TRIM5         remove bases from the begining of each read
                        (default:0)
  --rRNA                whether mapping to rRNA(Human)
  -v, --version         show program's version number and exit
```

  
##  Download annotation
  
  
```bash
wget www.test.com/index.gz
```
  
##  Usage
  
  
###  1. Submit the tasks to cluster
  
  
```bash
rapvis_submit.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```
  
###  2. Run in local
  
  
```bash
rapvis_run.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```
  
###  3. Caculated differently expressed genes
  
  
```bash
rapvis_DE.py -i input -p output
```
  
###  4. The Correlation coefficient between samples
  
  
```bash
rapvids_corr.py -i input -o output
```
  