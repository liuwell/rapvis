# rapvis: a tool for RNAseq processing and visualization
  
## Dependency
  
Required python version:
  
+ python >= 3.6
  
Sevral external software were depended for rapvis:
  
+ trimmomatic
+ STAR
+ hisat2
+ stringtie
+ bwa
+ samtools
  
## Mandatory
  
+ pandas >= 1.1.2
+ numpy
+ matplotlib
+ seaborn
+ GSEApy
+ rpy2
  
## Installation
  
### Installing from github
  
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
  
```bash
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

## Build genome index

You can download genome sequence and annotations GTF file from GENCODE. Strongly recommended for mouse and human (files marked with PRI) : <https://www.gencodegenes.org/>.

Other species can download from ENSEMBL, such as Zebrafish,  
genome sequences: <ftp://ftp.ensembl.org/pub/release-101/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz>  
GTF file: <ftp://ftp.ensembl.org/pub/release-101/gtf/danio_rerio/Danio_rerio.GRCz11.101.gtf.gz>

rapvis support **STAR** and **hisat2** for mapping.

### 1. build STAR index
  
```bash
$ rapvis_build.py -mapper STAR -genome GRCh38.primary_assembly.genome.fa.gz -gtf gencode.v35.primary_assembly.annotation.gtf.gz
```

### 2. build hisat2 index

```bash
$ rapvis_build.py -mapper hisat2 -genome GRCh38.primary_assembly.genome.fa.gz -gtf gencode.v35.primary_assembly.annotation.gtf.gz
```

## Usage
  
### 1. Run in local
  
```bash
$ rapvis_run.py -i tests/data1/ -o TestsResult -p 5 -lib STAR_index -m STAR
```
  
### 2. Submit the tasks to cluster
  
```bash
$ rapvis_submit.py -i tests/data1/ -o TestsResult -lib STAR_index -m STAR -p 5 -t 2
```

### 3. Caculated differently expressed genes

rapvis can caculated different expressed genes, based on **R limma**:

```bash
$ rapvis_DE.py -i input_TPM.txt -wt 0:3 -ko 3:6 -p output:
```
  
We can perform gene ontology enrichment analysis by **-go** aption, and the **-s** also needed for determining species:

```bash
$ rapvis_DE.py -i input_TPM.txt -wt 0:3 -ko 3:6 -p output -go -s Human
```

If the input gene matrix not be normalized, we can use **-norm** option to normalize, it based on **limma voom**:

```bash
$ rapvis_DE.py -i input_counts.txt -wt 0:3 -ko 3:6 -p output -norm
```

### 4. The Correlation coefficient between samples
  
We can get the correlation coeffcient heatmap of gene expresstion between samples:

```bash
$ rapvis_corr.py -i input_gene_TPM.txt
```

## Output

Several files included in the output directory:

+ **merge_gene_TPM.txt**  
*the gene expression profiles for all* samples, normalized by TPM
+ **merge_qc_percent.pdf**
*a quality contrl plot for trimmomatic*
+ **merge_mapping_percent.pdf**
*a barplot of the mapping percent in each sample*
+ **merge_gene_TPM_species_type.pdf**
*a stat a detected gene species in each sample, group by gene type*
+ **merge_gene_TPM_species_EI.pdf**
*a stat a detected gene species in each sample, group by expression interval*
+ **merge_gene_TPM_density.pdf**
*a density plot for gene expression distribution*
