#!/usr/bin/env python3

"""
Convert all FASTQ files in input folder to csv files in target folder
"""

import sys, os
import glob
import argparse

def parse_args():
	parser = argparse.ArgumentParser(description='Convert all FASTQ files in input folder to csv files in target folder in batch mode')
	parser.add_argument("-i", dest = "fastqfolder", type = str, required = True,
                              help = "Input fastq folder" )
	parser.add_argument("-o", dest = "csvfolder", type = str, default = "csv_converted",
                              help = "Output csv folder (Default: 'csv_converted' in current location)" )
	args = parser.parse_args()
	return args

def parseFILE(args):
	# create csv_converted folder
	os.mkdir(args.csvfolder)
	print("Reading FASTQ files from '%s' folder" %args.fastqfolder)
	print("Writing converted csv files to '%s' folder" %args.csvfolder)
	
	# call fastq_preprocess.py on each input FASTQ file
	for fl in sorted(glob.glob(args.fastqfolder + '/*.fastq')):
		input_file = os.path.basename(os.path.normpath(fl))
		output_file = input_file.replace(".fastq", ".csv")
		output_location =  os.path.join(args.csvfolder, output_file)
		os.system('python fastq_preprocess.py %s %s' % (fl, output_location))

def main():
	os.system('python --version')
	args = parse_args()
	parseFILE(args)

main()
