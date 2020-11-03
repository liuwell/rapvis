#!/usr/bin/env python3

import glob
import argparse
import numpy as np
import subprocess
import os

def process(R1, R2, output, adapter, threads, species, minlen, trim5):
	
	'''
	trim adapter by trimmomatic, mapping to genome by hisat2, transcript asemble by stringtie
	'''

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
	subprocess.call("trimmomatic PE -threads %d -phred33 %s %s %s %s %s %s ILLUMINACLIP:/sibcb1/wuliganglab1/liuwei/genome/adaptor/%s:1:30:10:5 SLIDINGWINDOW:4:20 MINLEN:%d HEADCROP:%d 2> %s" % (threads, R1, R2, out_R1_p, out_R1_u, out_R2_p, out_R2_u, adapter, minlen, trim5, out_log), shell=True)
	
	### Mapping by hisat2
	SummaryFile = prefix + "_hisat_summary.txt"
	HisatOut = prefix + "_hisat_sort.bam"
	subprocess.call("hisat2 -p %d -x /sibcb1/wuliganglab1/liuwei/genome/hisat2_index/%s/genome_tran -1 %s -2 %s -U %s,%s -t --dta --summary-file %s --new-summary|samtools sort -@ %d -m 10G -o %s" % (threads, species, out_R1_p, out_R2_p, out_R1_u, out_R2_u, SummaryFile, threads, HisatOut), shell=True)

	### Asemble by stringtie
	stringtieGTF = prefix + '_stringtie.gtf'
	stringtieGene = prefix + '_gene_abund.tab'
	subprocess.call("stringtie %s -e -G /sibcb1/wuliganglab1/liuwei/genome/hisat2_index/%s/annotation.gtf -p %d -o %s -A %s" % (HisatOut, species, threads, stringtieGTF, stringtieGene), shell=True)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='For RNAseq processing')

	parser.add_argument('-f1', required=True, help='the input data, R1')
	parser.add_argument('-f2', required=True, help='the input data, R2')
	parser.add_argument('-o', '--output', default = 'processed_data', help = 'output directory (default: processed_data)')
	parser.add_argument('-p', '--threads', default=5, type=int, help='number of threads (CPUs) to use (default: 5)')
	parser.add_argument('-s', '--species', default='Human', choices=['Human', 'Mouse', 'Rat', 'Rabbit', 'GoldenHamster', 'Zebrafish'], type=str, help='choose reference species for mapping and annotaion (default: Human)')
	parser.add_argument('-a', '--adapter', default='nextera', choices=['nextera', 'universal'], type=str, help='choose illumina adaptor (default: nextera)')
	parser.add_argument('--minlen', default=35, type=int, help='discard reads shorter than LEN (default: 35)')
	parser.add_argument('--trim5', default=0, type=int, help='remove bases from the begining of each read (default:0)')

	args = parser.parse_args()

	process(args.f1, args.f2, args.output, args.adapter, args.threads, args.species, args.minlen, args.trim5)


