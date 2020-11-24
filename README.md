### rapvis: a tool for RNAseq processing and visualization 

***
#### Dependency 
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


#### Mandatory
+ pandas >= 1.1.2
+ numpy
+ matplotlib
+ seaborn
+ GSEApy
+ rpy2

***
#### Install

<<<<<<< HEAD
```bash {cmd=false}
pip install rapvis
```
=======
```
# Clone remote repository
$ git clone https://github.com/liuwell/rapvis.git
>>>>>>> c3008a3502ada0deae1935f9d4d049fb04499fcf

# Install required python pacakge
$ cd rapvis
$ pip install -r requirements.txt

# Add execution path
# The path of current dir can get by "pwd"
$ echo "current_dir/rapvis" >> ~/.bashrc
$ source ~/.bashrc
```
#### Download annotation
```bash
wget www.test.com/index.gz
```
#### Usage
<<<<<<< HEAD
```bash {cmd=true}
/home/liuwei/miniconda2/envs/python36/bin/python3 ./rapvis/rapvis_submit.py -h
=======
#### 1. Submit the tasks to cluster
>>>>>>> c3008a3502ada0deae1935f9d4d049fb04499fcf
```

```bash {cmd=true}
/home/liuwei/miniconda2/envs/python36/bin/python3 ./rapvis/rapvis_run.py -h
```

<<<<<<< HEAD

```bash{cmd=false}
rapvis_submit.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```
=======
#### 2. Run in local
```bash
rapvis_run.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```
#### 3. Caculated differently expressed genes
```bash
rapvis_DE.py -i input -p output
```

#### 4. The Correlation coefficient between samples
```bash
rapvids_corr.py -i inout -o output
```

>>>>>>> c3008a3502ada0deae1935f9d4d049fb04499fcf
