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

### import R package
'''
from rpy2.robjects.packages import importr
# import R's "base" package
base = importr('base')
# import R's "utils" package
utils = importr('utils')
'''

### install R package
'''
# import rpy2's package module
import rpy2.robjects.packages as rpackages

# import R's utility package
utils = rpackages.importr('utils')

# select a mirror for R packages
utils.chooseCRANmirror(ind=1) # select the first mirror in the list

# R package names
packnames = ('ggplot2', 'limma')

# R vector of strings
from rpy2.robjects.vectors import StrVector

# Selectively install what needs to be install.
# We are fancy, just because we can.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))
'''


### get differently expressed gene by R package 'limma'
from rpy2 import robjects
robjects.r('''
			r_limma <- function(fi, fo, fo2, nwt, nko, norm){
			
			library(limma)
			
			data <- read.csv(fi, sep="\t", header=T, row.names=1)
			condition <-factor(c(rep("wt", nwt),rep("ko", nko)), levels=c("wt", "ko"))
			design <- model.matrix(~condition)
			
			data <- data[rowMeans(data)>0.5,]
			if(norm){
				data <- voom(data, design, plot = FALSE)
			}else{
				data <- log2(data+1)
			}
			
			fit <- lmFit(data, design)
			fit <- eBayes(fit)
			 
			DE <- topTable(fit, coef = 2, sort.by = "P", number = Inf)
			write.table(DE, file = fo, sep="\t", quote=F)

			### output normalized data or log2 transformed data
			#data <- round(data)
			#print(data)
			write.table(data, file = fo2, sep="\t", quote=F)
	}
	'''
)

r_func = robjects.r['r_limma']

def DE(fi, wt, ko, prefix, norm, go, species):
	
	wt = wt.split(':')
	wt = int(wt[1]) - int(wt[0])
	
	ko = ko.split(':')
	ko = int(ko[1]) - int(ko[0])
	
	out = prefix + "_DE_all.txt"
	out2 = prefix + "_norm_all.txt"
	r_func(fi, out, out2, wt, ko, norm)
	
	data = pd.read_table(out2, header=0, index_col=0)
	DE = pd.read_table(out, header=0, index_col=0)
	DE.rename(columns={'P.Value':'Pvalue', 'adj.P.Val':'AdjPval'}, inplace=True)

	### Cut Off
	fold_cutoff = 0.6
	pvalue_cutoff = 0.05
	
	DE['-log10(Pvalue)'] = -np.log10(DE['Pvalue'])
	DE['sig'] = 'normal'
	DE.loc[(DE.logFC > fold_cutoff)&(DE.Pvalue < pvalue_cutoff), 'sig'] = 'up'
	DE.loc[(DE.logFC < -fold_cutoff)&(DE.Pvalue < pvalue_cutoff), 'sig'] = 'down'
	
	sns.scatterplot(data=DE, x="logFC", y="-log10(Pvalue)", hue='sig', hue_order=('down', 'normal', 'up'), palette=("#377EB8", "grey", "#E41A1C"), alpha=0.8)
	plt.xlabel('log2(FoldChange)')
	plt.ylabel('-log10(P value)')
	
	out_volcano = prefix + "_DE_volcano.pdf"
	plt.savefig(out_volcano)
	plt.close()

	### Filter
	DE_filters = DE[DE['sig'].isin(['up', 'down'])]
	out_result = prefix + "_DE_FoldChange.txt"
	DE_filters.to_csv(out_result, header=True, sep='\t')
	
	### Heatmap
	filtered_genes = DE_filters.index
	data_filters = data[data.index.isin(filtered_genes)]
	
	out_heatmap_txt = prefix + "_DE_heatmap.txt"
	data_filters.to_csv(out_heatmap_txt, header=True, sep='\t')
	
	ncols = int(data_filters.shape[1])
	sns.clustermap(data_filters, cmap='RdBu_r', standard_scale=0, figsize=(ncols*0.6, 10))
	
	out_heatmap = prefix + "_DE_heatmap.pdf"
	plt.savefig(out_heatmap)
	plt.close()

	
	### GO, enrochr
	if go:
		gene_number = len(DE.index)
		gene_up = list(DE[DE.sig=='up'].index)
		gene_down = list(DE[DE.sig=='down'].index)
		
		### local model of GO analysis
		# un-regulated
		outdir1 = prefix + '_enrichr_gene_up'
		enr_up = gp.enrichr(gene_list=gene_up,
					 organism=species, description='up-regulated genes',
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
	
		### other tools for GO:
		### https://github.com/GuipengLi/SharePathway
		### https://github.com/jdrudolph/goenrich

	print("\n### Finished the different expressed analysis")
	print("\n### Input file: %s" % fi)
	print("\n### Output txt: %s, %s" % (out_result, out_heatmap_txt))
	print("\n### Output figure: %s, %s" % (out_volcano, out_heatmap))
	#print("\n### Output dir: %s, %s\n" % (outdir1, outdir2))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='For caculating differently expressed data, such as miRNA expression, gene expression, based on R package: limma')

	parser.add_argument('-i', '--input', required=True, help='the input data')
	parser.add_argument('-wt', required=True, type=str, help='the position of wildtype samples, such as 0:3')
	parser.add_argument('-ko', required=True, type=str, help='the position of knowkout samples, such as 3:6')
	parser.add_argument('-p', '--prefix', required=True, type=str, help='the prefix of output')
	parser.add_argument('-norm', action='store_true', help='perform the normalization by limma voom')
	parser.add_argument('-go', action='store_true', help='perform the Gene Ontology analysis')

	parser.add_argument('-s', '--species', default='none', type=str, help="choose enrichr library for GO analysis, such as: Human, Mouse, Yeast, Fly, Fish, Worm")
	#parser.add_argument('-s', '--species', default='Human', choices=['Human', 'Mouse', 'Yeast', 'Fly', 'Fish', 'Worm'], type=str, help="choose enrichr library for GO analysis, default: Human")

	args = parser.parse_args()

	DE(args.input, args.wt, args.ko, args.prefix, args.norm, args.go, args.species)

