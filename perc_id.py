#!/usr/bin/python

# This script compares an alignment file and calculates percent identity
# Author: Alex Chase
# Date: 2017-07-23
# Usage: python perc_id.py input.aln

from Bio import AlignIO
import sys
import os
import argparse

parser = argparse.ArgumentParser(description="This script will compare aligned sequences.")
parser.add_argument("-i", "--input", help="Input file name", required = True)

args = parser.parse_args()

def main():
	try: 
		input_file = args.input
		align = AlignIO.read(input_file, "fasta")
		print perc_identity(align)
		
		
	except IndexError:
		print("***ERROR USAGE:***\nCorrect Format for Use:\nperc_id.py -i input.aln\n")



def perc_identity(aln):
	i = 0
	for a in range(0,len(aln[0])):
		s = aln[:,a]
		if s == len(s) * s[0]:
			i += 1
	return 100*i/float(len(aln[0]))


if __name__ == '__main__':
	main()