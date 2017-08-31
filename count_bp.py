#!/usr/bin/python

# This script generates a text file with the base pair length of fasta sequence(s)
# Author: Alex Chase
# Date: 2015-09-23
# Usage: python count_bp.py input.fasta

from Bio import SeqIO
import sys
import os
import numpy


"""
This script creates an output text file showing the length of the sequences from an input fasta file.
It will also calculate the percent of GC of each sequnce.
"""


def main():
	try:
		count_bp()

	except IndexError:
		print("***ERROR USAGE:***\n\n")
		print("Correct Format for Use:\n")
		print("count_bp.py input.fasta\n")


def count_bp():

	in_file = sys.argv[1]

	out_file = sys.argv[1]

	if out_file.lower().endswith('.fa'):
		out_file = sys.argv[1].replace('.fa', '_bp-count.txt')

	elif out_file.lower().endswith('.fasta'):
		out_file = sys.argv[1].replace('.fasta', '_bp-count.txt')

	elif out_file.lower().endswith('.fna'):
		out_file = sys.argv[1].replace('.fna', '_bp-count.txt')

	elif out_file.lower().endswith('.faa'):
		out_file = sys.argv[1].replace('.faa', '_bp-count.txt')

	else:
		print "\n***ERROR USAGE:***\n\n"
		print "Not a recognized .fasta file\nMust be in either '.fa' or '.fasta' or '.fna' format."
		print "This script only works for nucleotide sequences."
		sys.exit()


	temp = []

	fasta_file_read = open(in_file, "r")
	text_file_out = open(out_file, "w")


	header = "Sequence_ID\tTotal_A\tTotal_G\tTotal_T\tTotal_C\tTotal_N\tWeird_chars\tTotal_bp\tGC_Percent"

	print >> text_file_out, header

	for record in SeqIO.parse(fasta_file_read, "fasta"):

		sequence = str(record.seq).upper()

		count = int(len(sequence))

		seq_id = str(record.name)

		A_count = sequence.count('A')
		G_count = sequence.count('G')
		T_count = sequence.count('T')
		C_count = sequence.count('C')
		U_count = sequence.count('U')
		T_total = T_count + U_count
		N_count = sequence.count('N')
		weird_char = numpy.float64(count - A_count - G_count - C_count - T_count - U_count - N_count)

		cg_per = numpy.float64(C_count + G_count) / count 

		output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f' % \
		(seq_id, A_count, G_count, T_total, C_count, N_count, weird_char, count, cg_per)

		print >> text_file_out, output_line

		temp.append(output_line)

	return temp

	text_file_out.close()
	fasta_file_read.close()


if __name__ == '__main__':
	main()