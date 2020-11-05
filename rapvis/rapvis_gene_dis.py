#!/usr/bin/env python3
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
from pandas.api.types import CategoricalDtype
from scipy import stats
import seaborn as sns
import datetime
import os
import re

###
def current_time():
	'''
	get current time
	'''
	return datetime.datetime.now().strftime('%b-%d-%Y %H:%M:%S')

###
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


