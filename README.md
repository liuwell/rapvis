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

```
# Clone remote repository
$ git clone https://github.com/liuwell/rapvis.git

# Install required python pacakge
$ cd rapvis
$ pip install -r requirements.txt

# Add execution path
# The path of current dir can get by "pwd"
$ echo "current_dir/rapvis" >> ~/.bashrc
$ source ~/.bashrc
```

#### Usage
```
rapvis_submit.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```

