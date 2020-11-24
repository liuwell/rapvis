#!/usr/bin/env python3

import subprocess
import os
import sys
import argparse
import re
import time
import datetime
import pandas as pd

def current_time():
	return datetime.datetime.now().strftime('%b-%d-%Y %H:%M:%S')


def STAR_index(threads, genome, gtf, sjdbOverhang):
	os.system("STAR --runMode genomeGenerate --runThreadN %d --genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang %d --genomeDir ." %(threads, genome, gtf, sjdbOverhang))


def hisat2_index(threads, genome, gtf):
	os.system("hisat2_extract_splice_sites.py %s > genome.ss" % gtf)
	os.system("hisat2_extract_exons.py %s > genome.exon" % gtf)
	os.system("hisat2-build -p %d %s --ss genome.ss --exon genome.exon genome_tran" % (threads, genome))


def GTF(gtf):
	d_name = {}
	d_id = {}

	sList = []
	with open(gtf) as f:
		for line in f:
			if not line.startswith('#'):
				(chrid, source, genetype, start, end, score, strand, phase, attribute_string) = line.rstrip().split('\t')
				attributes = {}
	
				if genetype=='gene':
					key_value_pair_set = attribute_string.split('; ')
					for key_value_pair in key_value_pair_set[: -1]:
					    key, value = key_value_pair.split(' ')
					    attributes[key] = value[1: -1] ### remove the first str and last str
					### specially for last field
					key, value = key_value_pair_set[-1][: -1].split(' ')
					attributes[key] = value[1: -1]
					
					gene_id = attributes.get('gene_id')
					gene_name = attributes.get('gene_name')
					gene_type = attributes.get('gene_biotype')
	
					#print(gene_id, gene_name, gene_type)
					d_name[gene_name] = gene_type
					d_id[gene_id] = gene_type

					sList.append([gene_id, gene_name, gene_type])
	'''
	df_name =pd.Series(d_name)
	df_name_csv = "gene_name.type.txt"
	df_name.to_csv(df_name_csv, header=False, sep='\t')

	df_id=pd.Series(d_id)
	df_id_csv = "gene_id.type.txt"
	df_id.to_csv(df_id_csv, header=False, sep='\t')
	'''
	df = pd.DataFrame(sList)
	df_csv = "gene_type.txt"
	df.to_csv(df_csv, header=False, index=True, sep='\t')



def main(mapper, threads, genome, gtf, sjdbOverhang):
	if mapper == 'STAR':
		STAR_index(threads, genome, gtf, sjdbOverhang)
	elif mapper == 'hisat2':
		hisat2_index(threads, genome, gtf)
	os.system("cp %s annotation.gtf" % gtf)


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='build genome index')

	parser.add_argument('-mapper', choices=['hisat2', 'STAR'], type=str, help='choose the mapping program')
	parser.add_argument('-genome', type=str, help='reference genome sequences')
	parser.add_argument('-gtf', type=str, help='annotation GTF file')
	parser.add_argument('-threads', default=5, type=int, help='number of threads (CPUs) to use (default: 5)')
	parser.add_argument('-sjdbOverhang', default=149, type=int, help='max(ReadLength)-1, the default value of 149 will work for Illumina 2x150b paired-end reads')
	
	args = parser.parse_args()
	
	print("\n%s ..... Start RNAseq processing" % (current_time()))
	start_time = time.time()

	main(args.mapper, args.threads, args.genome, args.gtf, args.sjdbOverhang)
	GTF(args.gtf)
	
	end_time = time.time()
	run_time = round((end_time - start_time)/60, 5)
	print("\n%s ..... Finished all. Used time: %s m\n" % (current_time(), run_time))






