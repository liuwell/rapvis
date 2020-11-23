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
+ rpy2

***
#### Install

```
# clone remote repository
git clone https://github.com/liuwell/rapvis.git
# install required python pacakge
pip install -r requirements.txt
# add execution path
echo "path_to_rapvis" >> ~/.bashrc
source ~/.bashrc
```

#### Usage
```
rapvis_submit.py -i rawdata/ -o processed -s Human -a universal -p 5 -t 2 --minlen 25 --trim5 3 --merge --rRNA
```

