#!/usr/bin/env python

import glob
import argparse
import re
import os
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns

###############################
def quality(fi):

	'''
	stat trimmomatic QC summary result

	'''
	
	files = sorted(glob.glob("%s/*/*trimmomatic.log" % fi))
	sList = []
	for f in files:
		with open(f) as handle:
			for line in handle:
				if re.match("Input Read Pairs", line):
					name = f.split("/")[-1].split("_")[0]
					line=line.strip().split(' ')
					total, bs, bsRatio, r1, r1Ratio, r2, r2Ratio, drop, dropRatio = line[3], line[6], line[7], line[11],line[12], line[16],line[17], line[19], line[20]
					s = Series([int(bs), int(r1), int(r2), int(drop)], index=["Both Surviving", "R1 Surviving", "R2 Surviving", "Dropped"], name=name)
					sList.append(s)

	df = DataFrame(sList)
	prefix = os.path.join(fi, 'merge_qc')
	bplot(df, prefix)

###############################
def mapping(fi):
	
	'''
	stat hisat2 mapping summary result

	'''
	
	files = sorted(glob.glob("%s/*/*hisat_summary.txt" % fi))
	sList = []
	for f in files:
		name = f.split("/")[-1].split("_")[0]
		f=open(f,'r')
		lines = f.readlines()
		x = [x.strip().split(' ')[-2] for x in lines]
		unmap = int(x[2])
		uniqmap = int(x[3]) + int(x[5])
		multimap = int(x[4])
		 
		s = Series([uniqmap, multimap, unmap], index=["UniqMapped", "MultipleMapped", "UnMapped"], name=name)
		sList.append(s)
		f.close()
	
	df = DataFrame(sList)
	prefix = os.path.join(fi, 'merge_mapping')
	bplot(df, prefix)

###############################
def rRNAratio(fi):
	'''
	stat the ratio of rRNA

	'''
	files = sorted(glob.glob("%s/*/*rRNA.count" % fi))
	d = {}
	sList = []
	for i in files:
		name = i.split("/")[-1].split("_")[0]
		with open(i) as f:
			for line in f:
				line = line.strip().split()
				d[line[0]] = int(line[1])
		s = Series(d, name=name)
		sList.append(s)

	df = DataFrame(sList)
	prefix = os.path.join(fi, 'merge_rRNA')
	bplot(df, prefix)

###############################
def bplot(df, prefix):
	'''
	barplot the summary result from trimmomatic and hisat2

	'''
	width = int(df.shape[0])
	height = 6
	fontsize =20
	if width >= 8 :
		width = math.log(width, 2) * 2 ### adjust the width of barplot

	### rawcount
	df.plot.bar(fontsize=fontsize, width=0.8, stacked=True, figsize=(width, height))
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=20)
	plt.ylabel('Counts', fontsize=20)
	plt.legend(loc=0, borderaxespad=0.5, fontsize = 'large', bbox_to_anchor=(1,1))
	out_box = prefix + "_raw.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()

	### percent
	df_pcts = df.div(df.sum(1).astype(float), axis=0) * 100
	df_pcts.plot(kind='bar', fontsize=fontsize, width=0.8, stacked=True, figsize=(width, height))
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=20)
	plt.ylabel('Percent(%)', fontsize=20)
	plt.legend(loc=0, borderaxespad=0.5, fontsize='large', bbox_to_anchor=(1,1))
	out_box = prefix + "_percent.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='For RNAseq processing')

	parser.add_argument('-i', '--input', required=True, help='the input directory, contain qc and mapping result')
	#parser.add_argument('-p', '--prefix', required=True, help='the prefix of output file')
	
	args = parser.parse_args()
	quality(args.input)
	mapping(args.input)



