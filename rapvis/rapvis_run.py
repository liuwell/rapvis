#!/usr/bin/env python3

import glob
import argparse
import numpy as np
import subprocess
import os
import time
import sys

from rapvis_merge import merge_profiles
from rapvis_submit import current_time
from rapvis_quality import *
from rapvis_gene_dis import *
import rapvis_rRNA


def process(fi, output, adapter, threads, species, minlen, trim5, rRNA):
	
	'''
	trim adapter by trimmomatic, mapping to genome by hisat2, transcript asemble by stringtie
	'''

	### get the data with fastq format
	files=[]
	fAll = glob.glob("%s/*" % fi)
	for f in fAll:
		if f.endswith('fastq') or f.endswith('fastq.gz') or f.endswith('fq.gz') or f.endswith('fq'):
			files.append(f)
	files=sorted(files)
	f_index=list(np.arange(0, len(files), 2))

	#f_num=0
	
	# set the path
	realpath = sys.path[0]
	index_path = "/home/liuwei/genome/hisat2_index/"
	
	for i in f_index:

		R1 = files[i]
		R2 = files[i+1]

		file_name = R1.split("/")[-1].split("_")[0]
		outdir = os.path.join(output, file_name)

		### make directory
		if not os.path.exists(outdir):
			try:
				os.makedirs(outdir)
			except Exception as e:
				pass

		prefix = os.path.join(outdir, file_name)
	
		out_R1_p = prefix + "_R1.fq.gz"
		out_R1_u = prefix + "_R1_unpaired.gz"
		out_R2_p = prefix + "_R2.fq.gz"
		out_R2_u = prefix + "_R2_unpaired.gz"
	
		out_log = prefix + "_trimmomatic.log"
	
		print("\nProcessing: %s, %s" % (R1,R2))
	
		### trimmomatic
		subprocess.call("trimmomatic PE -threads %d -phred33 %s %s %s %s %s %s ILLUMINACLIP:%s/../library/adapter/%s:1:30:10:5 SLIDINGWINDOW:4:20 MINLEN:%d HEADCROP:%d 2> %s" % (threads, R1, R2, out_R1_p, out_R1_u, out_R2_p, out_R2_u, realpath, adapter, minlen, trim5, out_log), shell=True)
		
		### Mapping by hisat2
		SummaryFile = prefix + "_hisat_summary.txt"
		HisatOut = prefix + "_hisat_sort.bam"
		subprocess.call("hisat2 -p %d -x %s/%s/genome_tran -1 %s -2 %s -U %s,%s -t --dta --summary-file %s --new-summary|samtools sort -@ %d -m 10G -o %s" % (threads, index_path, species, out_R1_p, out_R2_p, out_R1_u, out_R2_u, SummaryFile, threads, HisatOut), shell=True)
	
		### Asemble by stringtie
		stringtieGTF = prefix + '_stringtie.gtf'
		stringtieGene = prefix + '_gene_abund.tab'
		subprocess.call("stringtie %s -e -G %s/%s/annotation.gtf -p %d -o %s -A %s" % (HisatOut, index_path, species, threads, stringtieGTF, stringtieGene), shell=True)
		
		### for rRNA
		if rRNA:
			rapvis_rRNA.rRNA(R1, R2, output, threads)	

	else:
	
		fi = merge_profiles(args.output)
		gene_dis(fi, output, species)
		quality(output)
		mapping(output)
		
		### for rRNA
		if rRNA:
			rRNAratio(output)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='For RNAseq processing')

	parser.add_argument('-i', '--input', required=True, help='the input data')
	parser.add_argument('-o', '--output', default = 'processed_data', help = 'output directory (default: processed_data)')
	parser.add_argument('-p', '--threads', default=5, type=int, help='number of threads (CPUs) to use (default: 5)')
	parser.add_argument('-s', '--species', default='Human', choices=['Human', 'Mouse', 'Rat', 'Rabbit', 'GoldenHamster', 'Zebrafish'], type=str, help='choose reference species for mapping and annotaion (default: Human)')
	parser.add_argument('-a', '--adapter', default='nextera', choices=['nextera', 'universal'], type=str, help='choose illumina adaptor (default: nextera)')
	parser.add_argument('--minlen', default=35, type=int, help='discard reads shorter than LEN (default: 35)')
	parser.add_argument('--trim5', default=0, type=int, help='remove bases from the begining of each read (default:0)')
	parser.add_argument('--rRNA', action='store_true', help='whether mapping to rRNA')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.2')

	args = parser.parse_args()

	print("\n%s ..... Start RNAseq processing" % (current_time()))
	start_time = time.time()

	process(args.input, args.output, args.adapter, args.threads, args.species, args.minlen, args.trim5, args.rRNA)

	###
	end_time = time.time()
	run_time = round((end_time - start_time)/60, 5)
	print("\n%s ..... Finished all. Used time: %s m\n" % (current_time(), run_time))
