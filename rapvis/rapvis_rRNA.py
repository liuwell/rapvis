#!/usr/bin/env python3

import subprocess
import os
import argparse
import sys
def rRNA(R1, R2, output,  threads):
	'''
	for rRNA mapping by bwa
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

	output = prefix + '_rRNA_sort.bam'
	### use tht build-in index
	realpath = sys.path[0]
	subprocess.call("bwa mem %s/../library/rRNA_index/rRNA_45S_MT %s %s -t %d |samtools sort -@ 5 -m 10G -o %s" % (realpath, R1, R2, threads, output), shell=True)

	output2 = prefix + '_rRNA.bedgraph'
	subprocess.call("genomeCoverageBed -ibam %s -bg > %s" % (output, output2), shell=True)

	output3 = prefix + '_rRNA.count'
	subprocess.call("samtools view %s|awk '{a[$3]+=1}END{for(i in a)print i,a[i]}'|sort -k2nr > %s" %(output, output3), shell=True)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='For rRNA mapping, most in EV RNAseq')

	parser.add_argument('-f1', required=True, help='the input data, R1')
	parser.add_argument('-f2', required=True, help='the input data, R2')
	parser.add_argument('-o', '--output', default = 'processed_data', help = 'output directory (default: processed_data)')
	parser.add_argument('-p', '--threads', default=5, type=int, help='number of threads (CPUs) to use (default: 5)')


	args = parser.parse_args()
	rRNA(args.f1, args.f2, args.output, args.threads)
