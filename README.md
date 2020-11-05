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
+ python3.6+
+ pandas >= 1.1.2
+ numpy
+ matplotlib
+ seaborn
+ GSEApy

***
#### Install

```
pip install rapvis
```

#### Usage
```
rapvis_submit.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```

