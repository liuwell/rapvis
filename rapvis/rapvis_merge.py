#!/usr/bin/env python3

import math
from pandas import DataFrame
import pandas as pd
from itertools import islice
import glob
import os

### local function
from rapvis_general import current_time


def merge_profiles(output):
	
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
	
	else:
		print("\n### Merge profiles failed, it is not exsit in %s/*/ \n" %(output))
		exit(1)

###############################################
