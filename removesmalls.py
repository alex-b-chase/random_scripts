#! /usr/bin/python

# This script creates a new fasta file with sequences greater than a threshold minimum value
# Author: Alex Chase
# Date: 2015-09-22
# Usage: python removesmalls.py input_file.fasta output_file.fasta

from Bio import SeqIO
import sys
import os
import numpy as np
import argparse

__author__ = "AB Chase"

"""
Creates a new fasta file with sequences greater than a threshold minimum value.
Also, removes sequences with gaps comprising more than x% of the sequence length.

Call script:
removesmalls.py input_file.fasta output_file.fasta
"""

parser = argparse.ArgumentParser(description="This script will subset a fasta file by user parameters.")
parser.add_argument("-i", "--input", help="Input file name", required = True)
parser.add_argument("-o", "--output", help="Directs the output to a name of your choice", required = True)
parser.add_argument("-l", "--minlen", help='Minimum base pair length to keep', required = True)
parser.add_argument("-g", "--gcper", help='Percent Minimum GC content for of each contig', required = False, default=0)
parser.add_argument("-G", "--gap", help='Percent of reads with x percent of gaps allowed', required = False, default=0)

args = parser.parse_args()


def main():
	try:
		removesmalls()

	except IndexError:
		print("***ERROR USAGE:***\n\n")
		print("Correct Format for Use:\n")
		print(" %s -i <input.fasta> -o <output.fasta> -l <minlength> \n" % sys.argv[0])

def removesmalls():

	# GET ALL THE USER INPUTS FROM ABOVE

	# get the input / output files
	in_file = args.input
	output_file = args.output

	# get the percent %GC of output
	gc_per = args.gcper
	# get the minimum contig length
	min_len = args.minlen
	# get the percent of acceptable gap characters
	per_thresh = args.gap

	# check to see whether the user parameters are good to go
	while True:

		gc_per = args.gcper
		try:
			val = int(gc_per)
			break
		except ValueError:
			print("GC Percent is not an integer! Please do not include anything but a number 0-100")

		per_thresh = args.gap
		try:
			val = int(per_thresh)
			break
		except ValueError:
			print("Gapped Percent is not an integer! Please do not include anything but a number 0-100")

		min_len = args.minlen
		try:
			val = int(min_len)
			break
		except ValueError:
			print("Minimum contig length is not an integer! Please do not include anything but a number")

	if output_file.lower().endswith('.fa'):
		output_file = args.output
		reject_file = args.output.replace('.fa', '_rejected.fa')

	elif output_file.lower().endswith('.fasta'):
		output_file = args.output
		reject_file = args.output.replace('.fasta', '_rejected.fasta')

	elif output_file.lower().endswith('.fna'):
		output_file = args.output
		reject_file = args.output.replace('.fna', '_rejected.fna')

	elif output_file.lower().endswith('.fas'):
		output_file = args.output
		reject_file = args.output.replace('.fas', '_rejected.fas')

	elif output_file.lower().endswith('.faa'):
		output_file = args.output
		reject_file = args.output.replace('.faa', '_rejected.faa')

	else:
		print "\n***ERROR USAGE:***\n\n"
		print "Not a recognized .fasta file\nMust be in either '.fa' or '.fasta' or '.fna' or '.faa' format."
		print "This script only works for nucleotide or amino acid sequences. No fastq files allowed!"
		sys.exit()


	# Setup an empty list
	good_sequences = []
	bad_sequences = []


	fasta_file_read = open(in_file, "r")
	fasta_file_new = open(output_file, "w")
	fasta_file_reject = open(reject_file, "w")



	for record in SeqIO.parse(fasta_file_read, "fasta") :

		sequence = str(record.seq).upper()

		count = int(len(sequence))

		A_count = sequence.count('A')
		G_count = sequence.count('G')
		T_count = sequence.count('T')
		C_count = sequence.count('C')
		U_count = sequence.count('U')
		T_total = T_count + U_count
		GC_total = G_count + C_count
		N_count = sequence.count('N')

		cg_per = float(GC_total) / float(len(sequence)) * 100

		# remove sequences less than the threshold and throw out aligned sequences with more than 75% gaps
		if (gc_per == 0 and per_thresh == 0):

			if (len(sequence) >= int(min_len)):
				good_sequences.append(record)

			else:
				bad_sequences.append(record)

		elif (gc_per != 0 and per_thresh == 0):

			if (len(sequence) >= int(min_len) and cg_per >= float(gc_per)):
				good_sequences.append(record)

			else:
				bad_sequences.append(record)

		else:

			if (len(sequence) >= int(min_len) and cg_per >= float(gc_per) and (float(sequence.count("-")) / float(len(sequence))) * 100 <= int(per_thresh)):
				good_sequences.append(record)

			else:
				bad_sequences.append(record) 


 
	
	SeqIO.write(good_sequences, fasta_file_new, "fasta")
	SeqIO.write(bad_sequences, fasta_file_reject, "fasta")
	

	fasta_file_new.close()
	fasta_file_read.close()

	print "The %i good sequences are in %s" % (len(good_sequences), output_file)
	print "Please refer to %s for sequences that were discarded\n" % reject_file

if __name__ == '__main__':
	main()