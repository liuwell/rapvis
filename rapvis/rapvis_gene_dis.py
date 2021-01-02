#!/usr/bin/env python3
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from scipy import stats
import seaborn as sns
import datetime
import os
import re
import math

###
from rapvis_general import current_time


###
def gene_dis(fi, output, libpath):
	
	print("%s, Caculating gene expression pattern ... " % current_time())

	data = pd.read_table(fi, header=0)
	prefix = os.path.join(output, 'merge_gene_TPM')
	###
	data_melt = data.melt('gene', var_name='sample')
	data_melt = data_melt.query('value>0')
	data_melt.index = data_melt['gene']
	
	### Gene species by gene type
	gene_type = {}
	with open("%s/gene_type.txt" % libpath) as f:
		x = str(data_melt.index[0])
		if x.startswith("ENS"):
			for line in f:
				line=line.strip().split("\t")
				gene_type[line[1]] = line[3]
		else:
			for line in f:
				line=line.strip().split("\t")
				gene_type[line[2]] = line[3]


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

	
	# set width and height
	width = int(data.shape[0])
	height = 6
	fontsize =15
	if width >= 8 :
		width = math.log(width, 2) * 2 ### adjust the width of barplot
	else :
		width = width/1.5

	aspect = width/width
	#aspect = int(data.shape[1])
	#if aspect >3:
		#aspect = np.log(aspect) - 1 
	#if aspect >1:
	#	aspect = np.log(aspect)
	#else :
	#	aspect = aspect

	colors = list(reversed(sns.color_palette()[0:5]))
	hue_order = ["others", "pseudogene", "antisense", "lincRNA", "protein_coding"]
	sns.displot(data_melt2, x="sample", hue="gene_type", hue_order=hue_order, palette=colors, multiple="stack", shrink=.8, height=height, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=fontsize)
	plt.ylabel('Gene species', fontsize=fontsize)

	out_box = prefix + "_species_type.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	
	'''
	### Gene species
	sns.displot(data_melt, x="sample", shrink=.8, height=height, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=fontsize)
	plt.ylabel('Gene numbers', fontsize=fontsize)

	out_box = prefix + "_species.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	'''

	### Gene species by expression interval
	values = pd.cut(data_melt['value'], [0, 1, 5, 10, 50, 100, 1000, 1000000], labels=['0~1', '1~5', '5~10', '10~50', '50~100', '100~1000', '>1000'])
	data_melt = data_melt.copy() ### For SettingWithCopyWarning
	data_melt['ExpressionInterval'] = values
	#data_melt.loc[:,'ExpressionInterval'] = values
	
	sns.displot(data_melt, x="sample", hue="ExpressionInterval", multiple="stack", shrink=.8, height=height, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=fontsize)
	plt.ylabel('Gene species', fontsize=fontsize)

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


