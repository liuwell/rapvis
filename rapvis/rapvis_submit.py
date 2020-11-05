#!/usr/bin/env python3

import re
import os
import sys
import time
import glob
import datetime
import argparse
import subprocess
import numpy as np

### local function
from rapvis_quality import *

def current_time():
	'''
	get current time
	'''
	return datetime.datetime.now().strftime('%b-%d-%Y %H:%M:%S')


def GetRunningTasks(name):
	'''
	get the number of running tasks in the server
	'''
	subCounts = 'subCounts.txt'
	os.system("qstat > %s" % subCounts)
	n=0
	with open(subCounts) as sc:
		for line in sc:
			if re.search(name, line):
				n+=1
	os.system("rm %s" % subCounts)

	return n	


def SubmitTask(fi, output, adapter, threads, species, tasks, name, minlen, trim5, queue, rRNA):
	'''
	submit tasks to the server
	'''
	### get the data with fastq format
	files=[]
	fAll = glob.glob("%s/*" % fi)
	for f in fAll:
		if f.endswith('fastq') or f.endswith('fastq.gz') or f.endswith('fq.gz') or f.endswith('fq'):
			files.append(f)
	files=sorted(files)
	f_index=list(np.arange(0, len(files), 2))

	f_num=0
	for i in f_index:
		
		while True:
			### get task number
			n = GetRunningTasks(name)
			### check task number
			if n >= tasks :
				time.sleep(10)
				print("%s, Submitted Tasks: %d, total: %d" % (current_time(), f_num, len(f_index)))
	
			else:
			
				f_num+=1 # samples number
				R1 = files[i]
				R2 = files[i+1]
			
				tmp = "tmp.sh"
				f = open(tmp, 'w')
				f.write("#!/bin/bash\n")
				f.write("#$ -o %s.o\n" % name)
				f.write("#$ -e %s.e\n" % name)

				jobName = name + '.' + str(f_num)
				f.write("#$ -N %s\n" % jobName)
				
				f.write("source ~/.bashrc\n")
				f.write("source ~/.bash_profile\n")

				realpath = sys.path[0]
				f.write("python %s/rapvis_process.py -f1 %s -f2 %s -o %s -a %s -p %d -s %s --minlen %d --trim5 %d\n" %(realpath, R1, R2, output, adapter, threads, species, minlen, trim5))
				
				if rRNA:
					f.write("python %s/rapvis_rRNA.py -f1 %s -f2 %s -o %s -p %d\n" % (realpath, R1, R2, output, threads))
				
				f.close()
			
				subprocess.call("qsub -cwd -q %s %s" % (queue, tmp), shell=True)
				#subprocess.call("qsub -cwd -l node=4 -q %s %s" % (queue, tmp), shell=True)
				#subprocess.call("qsub -cwd -l mem_free=150G -q %s %s" % (queue, tmp), shell=True)
				subprocess.call("rm %s" % tmp, shell=True)
				time.sleep(1)

				break
	


###############################################
import math
from pandas import DataFrame
import pandas as pd
from itertools import islice

def merge_profiles(name, output):
	
	while True:
		
		n = GetRunningTasks(name)
		if n==0:
			print("%s, Merging profiles ... " % current_time())
			files = glob.glob("%s/*/*gene_abund.tab" % (output))
			
			if files  :
				files = sorted(files)
				dict_merge = {}
				for f in files:
					with open(f) as handle:
						for line in islice(handle, 1, None):
							line = line.strip().split("\t")
							k_map = line[1]
							k_RNA = f.split("/")[-2]
							count = float(line[8]) # TPM
							if k_map in dict_merge:
								dict_merge[k_map][k_RNA] = count
					
							else:
								tmp_dic = {}
								tmp_dic[k_RNA] = count
								dict_merge[k_map] = tmp_dic
		
				df = DataFrame(dict_merge).T
				df = df.fillna(value=0) ### fill NA to 0
				df_sum = DataFrame(df.sum(axis=1), columns=['sum'])
				df = df.join(df_sum)
				df = df.sort_values(by="sum", ascending=False) ### sort by sum
				df.drop(['sum'], axis=1, inplace=True)
				merge_out = os.path.join(output, "merge_gene_TPM.txt")
				df.to_csv(merge_out, sep="\t", header=True, index=True, index_label="gene", float_format="%.2f")
				return merge_out
				break
			
			else:
				print("\n### Merge profiles failed, it is not exsit in %s/*/ \n" %(output))
				exit(1)
		else:
			print("%s, Waitiing for task finished, remaining %d tasks" % (current_time(), n))
			time.sleep(10)

###############################################
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
from pandas.api.types import CategoricalDtype
from scipy import stats
import seaborn as sns

def gene_dis(fi, output, species):
	
	print("%s, Caculating gene expression pattern ... " % current_time())

	data = pd.read_table(fi, header=0)
	prefix = os.path.join(output, 'merge_gene_TPM')
	###
	data_melt = data.melt('gene', var_name='sample')
	data_melt = data_melt.query('value>0')
	data_melt.index = data_melt['gene']
	
	### Gene species by gene type
	index_path = "/home/liuwei/genome/hisat2_index"
	gene_type = {}
	with open("%s/%s/gene.type" % (index_path,species)) as f:
		for line in f:
			line=line.strip().split("\t")
			gene_type[line[0]] = line[1]

	gene_type = pd.Series(gene_type, name='gene_type', dtype="string")
	
	type_list = ["protein_coding", "pseudogene", "lincRNA", "antisense"]
	for i in range(0, len(gene_type.index)):
		if gene_type[i] in type_list:
			pass
		elif re.search("pseudogene", gene_type[i]):
			gene_type[i] = 'pseudogene'
			
		else:
			gene_type[i] = 'others'

	# Categories
	data_melt2 = pd.merge(data_melt, gene_type, how='left', sort=False, right_index=True,left_index=True)
	cat_type = CategoricalDtype(categories=data.columns[1:], ordered =True)
	data_melt2['sample'] = data_melt2['sample'].astype(cat_type)

	aspect = int(data.shape[1])
	if aspect >3:
		aspect = np.log(aspect) - 1 
	else :
		aspect = aspect/3

	colors = list(reversed(sns.color_palette()[0:5]))
	hue_order = ["others", "pseudogene", "antisense", "lincRNA", "protein_coding"]
	sns.displot(data_melt2, x="sample", hue="gene_type", hue_order=hue_order, palette=colors, multiple="stack", shrink=.8, height=4, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=15)
	plt.ylabel('Gene species', fontsize=15)

	out_box = prefix + "_species_type.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	
	### Gene species
	sns.displot(data_melt, x="sample", shrink=.8, height=4, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=15)
	plt.ylabel('Gene numbers', fontsize=15)

	out_box = prefix + "_species.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()

	### Gene species by expression interval
	values = pd.cut(data_melt['value'], [0, 1, 5, 10, 50, 100, 1000, 1000000], labels=['0~1', '1~5', '5~10', '10~50', '50~100', '100~1000', '>1000'])
	data_melt = data_melt.copy() ### For SettingWithCopyWarning
	data_melt['ExpressionInterval'] = values
	#data_melt.loc[:,'ExpressionInterval'] = values
	
	sns.displot(data_melt, x="sample", hue="ExpressionInterval", multiple="stack", shrink=.8, height=4, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=15)
	plt.ylabel('Gene species', fontsize=15)

	out_box = prefix + "_species_EI.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	
	### expression density
	#data_melt['log2value'] = np.log2(data_melt['value'])
	#sns.displot(data=data_melt, x="log2value", kind="kde", hue='sample', height=4, aspect=1.4, common_norm=False)
	sns.kdeplot(data=data_melt, x="value", hue='sample', log_scale=True, common_norm=False)
	plt.xlabel('log10(TPM)', fontsize=15)
	out_box = prefix + "_density.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()


# =========================================== #
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='A tool for RNAseq processing and visualization')
	#rapvis
	parser.add_argument('-i', '--input', required=True, help='the input data')
	parser.add_argument('-o', '--output', default = 'processed_data', help = 'output directory (default: processed_data)')
	parser.add_argument('-s', '--species', default='Human', choices=['Human', 'Mouse', 'Rat', 'Rabbit', 'GoldenHamster', 'Zebrafish'], type=str, help='choose reference species for mapping and annotaion (default: Human)')
	parser.add_argument('-a', '--adapter', default='nextera', choices=['nextera', 'universal'], type=str, help='choose illumina adaptor (default: nextera)')
	parser.add_argument('-p', '--threads', default=5, type=int, help='number of threads (CPUs) to use (default: 5)')
	parser.add_argument('-t', '--tasks', default=2, type=int, help='number of submitted tasks (default: 2)')
	parser.add_argument('-n', '--name', default = 'RNAseq', type=str, help = 'project name (default: RNAseq)')
	parser.add_argument('--minlen', default=35, type=int, help='discard reads shorter than minlen (default: 35)')
	parser.add_argument('--trim5', default=0, type=int, help='remove bases from the begining of each read (default:0)')
	parser.add_argument('-q', default='b1.q', choices=['b1.q', 'g1.q'], type=str, help='bind job to queue(s)')
	parser.add_argument('--merge', action='store_true', help='merge gene expression profiles and plot distribution pattern')
	parser.add_argument('--rRNA', action='store_true', help='whether mapping to rRNA')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.2')

	args = parser.parse_args()

	print("\n%s ..... Start RNAseq processing" % (current_time()))
	start_time = time.time()

	SubmitTask(args.input, args.output, args.adapter, args.threads, args.species, args.tasks, args.name, args.minlen, args.trim5, args.q, args.rRNA)
	
	if args.merge:
		fi = merge_profiles(args.name, args.output)
		gene_dis(fi, args.output, args.species)
		
		quality(args.output)
		mapping(args.output)

	if args.rRNA:
		rRNAratio(args.output)

	end_time = time.time()
	run_time = round((end_time - start_time)/60, 5)
	print("\n%s ..... Finished all. Used time: %s m\n" % (current_time(), run_time))


