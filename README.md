### rapvis: a tool for RNAseq processing and visualization 

***
#### Dependency 
Sevral external software were depended for rapvis:

+ trimmomatic 
+ hisat2 
+ stringtie
+ bwa
+ samtools
+ STAR
+ bedtools

#### Mandatory
+ python >= 3.6.5
+ pandas >= 1.1.2
+ numpy
+ matplotlib
+ seaborn
+ GSEApy

***
#### Install

```bash {cmd=false}
pip install rapvis
```

#### Usage
```bash {cmd=true}
/home/liuwei/miniconda2/envs/python36/bin/python3 ./rapvis/rapvis_submit.py -h
```

```bash {cmd=true}
/home/liuwei/miniconda2/envs/python36/bin/python3 ./rapvis/rapvis_run.py -h
```


```bash{cmd=false}
rapvis_submit.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```