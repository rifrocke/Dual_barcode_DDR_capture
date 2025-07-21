#!/usr/bin/env python3

"""
Remove the 5_end and 3_end adapter of all FASTQ files in input folder to FASTQ files in target folder
"""

import sys, os
import glob
import argparse

def parse_args():
	parser = argparse.ArgumentParser(description='Remove the 5_end and 3_end adapter of all FASTQ files in input folder to FASTQ files in target folder in batch mode')
	parser.add_argument("-i", dest = "fastq", type = str, required = True,
                              help = "Input fastq folder" )
	parser.add_argument("-a", dest = "adapter", type = str, required = True,
                              help = "Specify the adapters in 5 and 3 ends" )
	parser.add_argument("-o", dest = "fastq_noadpt", type = str, default = "fastq_no_adapter",
                              help = "Output fastq files without adapters folder (Default: 'fastq_no_adapter' in current location)" )
	args = parser.parse_args()
	return args

def parseFILE(args):
	# create output folder
	os.mkdir(args.fastq_noadpt)
	print("Reading FASTQ files from '%s' folder" %args.fastq)
	print("Writing Normalized FASTQ files to '%s' folder" %args.fastq_noadpt)
	
	# call cutadapt on each input FASTQ file
	for file in sorted(glob.glob(args.fastq + '/*.fastq')):
		input_file = os.path.basename(os.path.normpath(file))
		output_file = input_file.replace(".fastq", "_nodapt.fastq")
		output_location =  os.path.join(args.fastq_noadpt, output_file)
		os.system('cutadapt -g %s -o %s %s' % (args.adapter, output_location, file))

def main():
	args = parse_args()
	parseFILE(args)

main()
