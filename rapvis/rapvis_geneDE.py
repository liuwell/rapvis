#!/usr/bin/env python3
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import os
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import argparse

import gseapy as gp
from gseapy.plot import barplot, dotplot

def DE(fi, wt, ko, prefix, species):

	data = pd.read_table(fi, header=0, index_col=0)
	data = np.log2(data + 1)
	
	#####
	wt = wt.split(':')
	wt1 = int(wt[0])
	wt2 = int(wt[1])
	ko = ko.split(':')
	ko1 = int(ko[0])
	ko2 = int(ko[1])

	# The mean expression of wt samples for each genes
	wt = data.iloc[:, wt1:wt2].mean(axis=1)
	# The mean expression of ko samples for each genes
	ko = data.iloc[:, ko1:ko2].mean(axis=1)
	
	# FoldChange, ko vs wt
	foldchange = ko - wt
	
	# P value
	pvalue = []
	gene_number = len(data.index)
	for i in range(0, gene_number):
		ttest = stats.ttest_ind(data.iloc[i,wt1:wt2], data.iloc[i, ko1:ko2])
		pvalue.append(ttest[1])
	
	### Cut Off
	fold_cutoff = 1
	pvalue_cutoff = 0.05
	
	### vocano plot
	pvalue_arr = np.asarray(pvalue)
	result = pd.DataFrame({'log2MeanTreat': ko, 'log2MeanCtrl': wt, 'log2FoldChange': foldchange, 'pvalue': pvalue_arr,})
	result['-log10(pvalue)'] = -np.log10(result['pvalue'])
	
	result['sig'] = 'normal'
	result.loc[(result.log2FoldChange > fold_cutoff)&(result.pvalue < pvalue_cutoff), 'sig'] = 'up'
	result.loc[(result.log2FoldChange < -fold_cutoff)&(result.pvalue < pvalue_cutoff), 'sig'] = 'down'
	
	sns.scatterplot(x="log2FoldChange", y="-log10(pvalue)", hue='sig', hue_order=('down', 'normal', 'up'), palette=("#377EB8", "grey", "#E41A1C"), data=result)
	plt.xlabel('log2(FoldChange)')
	plt.ylabel('-log10(P value)')
	
	out_volcano = prefix + "_DE_volcano.pdf"
	plt.savefig(out_volcano)
	plt.close()

	#out_result = prefix + "_FoldChange.txt"
	#x = result.loc[result.sig!='normal', :] 
	#x.to_csv(out_result, header=True, sep='\t')

	### Filter
	filtered_ids=[]
	for i in range(0, gene_number):
		### FoldChange and P value
		if(abs(foldchange[i]) >= fold_cutoff) and (pvalue[i] <= pvalue_cutoff):
			### mean value
			if wt[i] > 1 or ko[i]>1:
				filtered_ids.append(i)
	
	out_result = prefix + "_DE_FoldChange.txt"
	x = result.iloc[filtered_ids, :] 
	x.to_csv(out_result, header=True, sep='\t')
	
	
	### Heatmap
	filt1 = data.iloc[filtered_ids, wt1:wt2]
	filt2 = data.iloc[filtered_ids, ko1:ko2]
	filtered = filt1.join(filt2)

	out_heatmap_txt = prefix + "_DE_heatmap.txt"
	filtered.to_csv(out_heatmap_txt, header=True, sep='\t')
	
	ncols = int(data.shape[1])
	#print(ncols)
	sns.clustermap(filtered, cmap='RdBu_r', standard_scale=0, figsize=(ncols*0.5, 10))
	
	out_heatmap = prefix + "_DE_heatmap.pdf"
	plt.savefig(out_heatmap)
	plt.close()

	### GO, enrochr
	gene_up = list(x[x.sig=='up'].index)
	gene_down = list(x[x.sig=='down'].index)
	

	### local model of GO analysis
	# un-regulated
	outdir1 = prefix + '_enrichr_gene_up'
	enr_up = gp.enrichr(gene_list=gene_up,
				 organism=species, description='up-regulated genes',
                 #gene_sets="c5.bp.v7.0.symbols.gmt", # ['KEGG_2016','KEGG_2013']
                 gene_sets='GO_Biological_Process_2018',
                 background=gene_number, # or the number of genes, e.g 20000
                 outdir=outdir1,
                 cutoff=0.05, # only used for testing.
				 format='pdf', figsize=(8, 7), top_term=20, no_plot=False, 
                 verbose=True)
	
	# down-regulated
	outdir2 = prefix + '_enrichr_gene_down'
	enr_down = gp.enrichr(gene_list=gene_down,
				 organism=species, description='down-regulated genes',
                 #gene_sets="c5.bp.v7.0.symbols.gmt", # ['KEGG_2016','KEGG_2013']
                 gene_sets='GO_Biological_Process_2018',
				 background=gene_number, # or the number of genes, e.g 20000
                 outdir=outdir2,
                 cutoff=0.05, # only used for testing.
				 format='pdf', figsize=(8, 7), top_term=20, no_plot=False, 
                 verbose=True)

	# to save your figure, make sure that ``ofname`` is not None
	#barplot(enr_ip.res2d,title='KEGG_2013', ofname='test_GO.pdf')

	### enrichr, Command line
	#gseapy enrichr -i data/gene_list.txt --ds BP2017 -g GO_Biological_Process_2017 -v -o test_enrichr_BP
	### GSEA, Command line
	#gseapy gsea -d data/P53_resampling_data.txt -g KEGG_2016 -c data/P53.cls -o gsea_report -v -t phenotype

	###
	### other tools for GO:
	### https://github.com/GuipengLi/SharePathway
	### https://github.com/jdrudolph/goenrich
	###

	print("\n### Finished the different expressed analysis")
	print("\n### Input file: %s" % fi)
	print("\n### Output txt: %s, %s" % (out_result, out_heatmap_txt))
	print("\n### Output figure: %s, %s" % (out_volcano, out_heatmap))
	print("\n### Output dir: %s, %s\n" % (outdir1, outdir2))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='For caculating differently expressed data, such as miRNA expression, gene expression')

	parser.add_argument('-i', '--input', required=True, help='the input data')
	parser.add_argument('-wt', required=True, type=str, help='the position of wildtype samples, such as 0:3')
	parser.add_argument('-ko', required=True, type=str, help='the position of knowkout samples, such as 3:6')
	parser.add_argument('-p', '--prefix', required=True, type=str, help='the prefix of output')
	parser.add_argument('-s', '--species', default='Human', choices=['Human', 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm'], type=str, help="choose enrichr library for GO analysis, default: Human")

	args = parser.parse_args()

	DE(args.input, args.wt, args.ko, args.prefix, args.species)


