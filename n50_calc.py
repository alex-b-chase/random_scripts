#!/usr/bin/python

# This script calculates n50 and other information from genomes
# should be embedded in shell script to loop through many files
# Author: Alex Chase
# Date: 2015-12-01
# Usage: python n50_calc.py input_fasta

from Bio import SeqIO

import sys
import numpy
import os

"""
N50 is metric that is often associated with the length of contigs post-assembly


Given a set of contigs, each with its own length, the N50 is the:

Length N for which 50% of all bases in the sequences are in a sequence of length L < N.  
This can be found mathematically as follows: Take a list L of positive integers. 
Create another list N , which is identical to L, except that every element n in L has been replaced with n copies of itself. 
Then the median of N is the N50 of L. 

For example: If L = {2, 2, 2, 3, 3, 4, 8, 8}, then N consists of six 2s, six 3s, four 4s, and sixteen 8s; the N50 of L is the median of L , which is 6.


To execute the script: 

python n50_calc.py input_fasta

"""


def main():

	try: 
		n50_calc()
		
		
	except IndexError:
		print("***ERROR USAGE:***\n\nCorrect Format for Use:\nn50_calc.py input_fasta_file.fasta\n")


	except IOError:
			print "***ERROR - INCORRECT USAGE***\n"
			print "Do not put '-i' or any other parameter marker.\n"
			print "Script Usage:\nn50_calc.py input_file.fasta\n"


def n50_calc():

	in_file = sys.argv[1]

	out_file = sys.argv[1]

	if out_file.lower().endswith('.fa'):
		out_file = sys.argv[1].replace('.fa', '_n50.txt')

	elif out_file.lower().endswith('.fasta'):
		out_file = sys.argv[1].replace('.fasta', '_n50.txt')

	elif out_file.lower().endswith('.fna'):
		out_file = sys.argv[1].replace('.fna', '_n50.txt')

	else:
		print "\n***ERROR USAGE:***\n\n"
		print "Not a recognized .fasta file\nMust be in either '.fa' or '.fasta' or '.fna' format."
		print "This script only works for nucleotide sequences."
		sys.exit()


	fasta_file_read = open(in_file, "rU")
	text_file_out = open(out_file, "w")


	seqlength = []
	unique = []
	n50 = []

	header = "Sequence_ID\tTotal_Contigs\tGenome_length\tLargest_Contig\tn50\tGC_Percent"

	print >> text_file_out, header

	a = 0
	c = 0
	t = 0
	g = 0
	n = 0
	total_count = 0

	# calculate the largest contig in the file
	for record in SeqIO.parse(in_file, 'fasta'):
		
		bp = len(record.seq)
		seqlength.append(bp)

		# for GC calculation
		sequence = str(record.seq).upper()
		count = int(len(sequence))
		total_count += count

		A_count = sequence.count('A')
		anew = numpy.float64(A_count)
		a += anew

		G_count = sequence.count('G')
		gnew = numpy.float64(G_count)
		g += gnew
		
		C_count = sequence.count('C')
		cnew = numpy.float64(C_count)
		c += cnew

		T_count = sequence.count('T')
		U_count = sequence.count('U')
		T_total = T_count + U_count
		tnew = numpy.float64(T_count)
		t += tnew

		N_count = sequence.count('N')
		nnew = numpy.float64(N_count)
		n += nnew


	seqlength = sorted(seqlength)
	total_contig = len(seqlength)

	largest_contig = seqlength[-1]

	# calculate the %GC for the file
	gcpercent = numpy.float64(c + g) / total_count * 100
	
	# calculate the n50

	# create a unique set of contigs
	for entry in seqlength:
		if not entry in unique:
			unique.append(entry)

	for entry in unique:
		multiplier = seqlength.count(entry) * entry
		for i in range(multiplier):
			n50.append(entry)

	index = len(n50) / 2
	avg = []

	if index % 2==0:
		first = n50[index - 1]
		second = n50[index]
		avg.append(first)
		avg.append(second)
		n50 = numpy.mean(avg)

	else:
		n50 = n50[index - 1]

	seq_id = os.path.splitext(sys.argv[1])[0]

	output_line = '%s\t%d\t%d\t%d\t%d\t%d' % \
	(seq_id, total_contig, total_count, largest_contig, n50, gcpercent)

	print >> text_file_out, output_line

	text_file_out.close()
	fasta_file_read.close()





if __name__ == '__main__':
	main()
