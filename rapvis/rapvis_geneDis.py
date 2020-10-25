#!/usr/bin/env python3
# python3.6
# ref link: https://www.jianshu.com/p/91c98585b79b
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os,re
import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from scipy import stats
import seaborn as sns
import argparse

def gene_dis(fi, prefix):

	data = pd.read_table(fi, header=0)
	###
	if data.columns[0] != 'gene':
		data.rename(columns={data.columns[0]:'gene'}, inplace=True)
	
	print(data.columns)
	
	data_melt = data.melt('gene', var_name='sample')
	data_melt = data_melt.query('value>0')
	data_melt.index = data_melt['gene']
	#print(data)
	#data_melt = data.melt(var_name='sample')
	#print(data_melt)
	#data_melt = data_melt.query('value>0')
	#data_melt.index = data_melt['gene']
	### Gene species by gene type
	gene_type = {}
	with open("/sibcb1/wuliganglab1/liuwei/genome/hisat2_index/Ensembl96_GRCh38_GRcm38.gene.type") as f:
		for line in f:
			line=line.strip().split("\t")
			gene_type[line[1]] = line[2]

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
	#cat_type = CategoricalDtype(categories=["protein_coding", "pseudogene", "lincRNA", "antisense", "others"], ordered =True)
	#cat_type = CategoricalDtype(categories=["others", "pseudogene", "antisense", "lincRNA", "protein_coding"], ordered =True)
	#gene_type = gene_type.astype(cat_type)
	data_melt2 = pd.merge(data_melt, gene_type, how='left', sort=False, right_index=True,left_index=True)
	cat_type = CategoricalDtype(categories=data.columns[1:], ordered =True)
	data_melt2['sample'] = data_melt2['sample'].astype(cat_type)

	aspect = np.log(int(data.shape[1])) - 1 
	colors = list(reversed(sns.color_palette()[0:5]))
	hue_order = ["others", "pseudogene", "antisense", "lincRNA", "protein_coding"]
	sns.displot(data_melt2, x="sample", hue="gene_type", hue_order=hue_order, palette=colors, multiple="stack", shrink=.8, height=4, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=15)
	plt.ylabel('Gene species', fontsize=15)

	out_box = prefix + "_gene_species_type.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	
	### Gene species
	sns.displot(data_melt, x="sample", shrink=.8, height=4, aspect=aspect)
	plt.xticks(rotation=90)
	plt.xlabel('Samples', fontsize=15)
	plt.ylabel('Gene numbers', fontsize=15)

	out_box = prefix + "_gene_species.pdf"
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

	out_box = prefix + "_gene_species_EI.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	
	### expression density
	#data_melt['log2value'] = np.log2(data_melt['value'])
	#sns.displot(data=data_melt, x="log2value", kind="kde", hue='sample', height=4, aspect=1.4, common_norm=False)
	sns.kdeplot(data=data_melt, x="value", hue='sample', log_scale=True, common_norm=False)
	plt.xlabel('log10(TPM)', fontsize=15)
	out_box = prefix + "_gene_density.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()

	'''
	# Boxplot of the expression data
	color = {'boxes': 'DarkGreen', 'whiskers': 'DarkOrange', 'medians': 'DarkBlue', 'caps': 'Gray'}
	data.plot(kind='box', color=color, sym='r.', title=prefix)
	plt.xticks(rotation=40)
	plt.xlabel('Samples', fontsize=15)
	plt.ylabel('log2(TPM+1)', fontsize=15)

	out_box = prefix + "_boxplot.pdf"
	plt.savefig(out_box, bbox_inches='tight')
	plt.close()
	
	# Density plot of the expression data
	data.plot(kind='density', title=prefix)
	#data.plot(kind='hist', alpha=0.5, title=prefix, stacked=True)
	#data.plot.kde()
	out_density = prefix + "_density.pdf"
	plt.xlabel('log2(TPM+1)', fontsize=15)
	plt.savefig(out_density)
	plt.close()
	'''

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='For the distribution of gene expression')

	parser.add_argument('-i', '--input', required=True, help='the input data')
	parser.add_argument('-p', '--prefix', required=True, type=str, help='the prefix of output')

	args = parser.parse_args()

	gene_dis(args.input, args.prefix)




