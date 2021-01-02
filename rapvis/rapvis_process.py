#!/usr/bin/env python3

import glob
import argparse
import numpy as np
import subprocess
import os
import time
import sys

#from rapvis_merge import merge_profiles, merge_gene_counts
#from rapvis_gene_dis import gene_dis
#from rapvis_quality import rRNAratio
from rapvis_general import current_time
import rapvis_rRNA

def process2(R1, R2, output, adapter, threads, libpath, mapper, minlen, trim5, counts, rRNA):

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

	print("\n%s  Processing: %s, %s" % (current_time(), R1,R2))

	realpath = sys.path[0]
	### trimmomatic
	subprocess.call("trimmomatic PE -threads %d -phred33 %s %s %s %s %s %s ILLUMINACLIP:%s/../library/adapter/%s:1:30:10:5 SLIDINGWINDOW:4:20 MINLEN:%d HEADCROP:%d 2> %s" % (threads, R1, R2, out_R1_p, out_R1_u, out_R2_p, out_R2_u, realpath, adapter, minlen, trim5, out_log), shell=True)
	
	### Mapping by hisat2
	if mapper == 'hisat2':
		SummaryFile = prefix + "_hisat_summary.txt"
		MapOut = prefix + "_hisat_sort.bam"
		subprocess.call("hisat2 -p %d -x %s/genome_tran -1 %s -2 %s -U %s,%s -t --dta --summary-file %s --new-summary|samtools sort -@ %d -m 10G -o %s" % (threads, libpath, out_R1_p, out_R2_p, out_R1_u, out_R2_u, SummaryFile, threads, MapOut), shell=True)
	
	### Mapping by STAR
	elif mapper == 'STAR':
		STARprefix = prefix + "_STAR_"
		subprocess.call("STAR --runThreadN %d --outSAMtype BAM SortedByCoordinate --genomeDir %s --readFilesIn %s %s --readFilesCommand zcat --outFileNamePrefix %s --quantMode GeneCounts --outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 --outFilterMatchNmin 0  --outFilterMismatchNmax 2" % (threads, libpath, out_R1_p, out_R2_p, STARprefix), shell=True)		

		MapOut = prefix + "_STAR_Aligned.sortedByCoord.out.bam" ## sorted bam file
		
	### Asemble by stringtie
	print("%s Asemble ..." % current_time())
	stringtieGTF = prefix + '_stringtie.gtf'
	stringtieGene = prefix + '_gene_abund.tab'
	subprocess.call("stringtie %s -e -G %s/annotation.gtf -p %d -o %s -A %s" % (MapOut, libpath, threads, stringtieGTF, stringtieGene), shell=True)
	
	### Gene counts
	if counts:
		countOut = prefix + '_gene_counts.txt'
		subprocess.call("featureCounts -a %s/annotation.gtf -o %s %s -t exon -g gene_name -T %d -Q 30 -p" % (libpath, countOut, MapOut, threads), shell=True)
	
	### rRNA
	if rRNA:
		rapvis_rRNA.rRNA(R1, R2, output, threads)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A tool for RNAseq processing and visualization')

	parser.add_argument('-R1', required=True, help='the input data R1')
	parser.add_argument('-R2', required=True, help='the input data R2')
	parser.add_argument('-o', '--output', default = 'processed_data', help = 'output directory (default: processed_data)')
	parser.add_argument('-p', '--threads', default=5, type=int, help='number of threads (CPUs) to use (default: 5)')
	#parser.add_argument('-s', '--species', default='Human', choices=['Human', 'Mouse', 'Rat', 'Rabbit', 'GoldenHamster', 'Zebrafish'], type=str, help='choose reference species for mapping and annotaion (default: Human)')
	parser.add_argument('-lib', '--libraryPath', type=str, help='choose reference species for mapping and annotaion')
	parser.add_argument('-m', '--mapper', default='hisat2', choices=['hisat2', 'STAR'], type=str, help='choose the mapping program (default: hisat2)')
	parser.add_argument('-a', '--adapter', default='nextera', type=str, help='choose illumina adaptor (default: nextera), choices {nextera, universal, pAAAAA}')
	parser.add_argument('--minlen', default=35, type=int, help='discard reads shorter than LEN (default: 35)')
	parser.add_argument('--trim5', default=0, type=int, help='remove bases from the begining of each read (default:0)')
	parser.add_argument('--counts', action='store_true', help='Get gene counts')
	parser.add_argument('--rRNA', action='store_true', help='whether mapping to rRNA(Human)')
	parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.2')

	args = parser.parse_args()

	#print("\n%s ..... Start RNAseq processing" % (current_time()))
	#start_time = time.time()

	process2(args.R1, args.R2, args.output, args.adapter, args.threads, args.libraryPath, args.mapper, args.minlen, args.trim5, args.counts, args.rRNA)

	###
	#end_time = time.time()
	#run_time = round((end_time - start_time)/60, 5)
	#print("\n%s ..... Finished all. Used time: %s m\n" % (current_time(), run_time))
