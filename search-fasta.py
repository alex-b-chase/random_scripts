#!/usr/bin/python

# This script searches against a reference file and pulls out those from larger fasta file
# Author: Alex Chase
# Date: 2014-11-22
# Usage: python search-fasta.py -i <input.fasta> -o <output.fasta> -m <mappingfile>

from Bio import SeqIO
import os
import sys
import argparse

"""
This script takes a large fasta file and creates a smaller fasta file with desired sequences.
These desired sequences should be in a text file to create a list <mappingfile>.
The <mappingfile> should be one sequence per line
"""

parser = argparse.ArgumentParser(description="This script will subset a fasta file by user parameters.")
parser.add_argument("-i", "--input", help="Input file name", required = True)
parser.add_argument("-o", "--output", help="Directs the output to a name of your choice", required = True)
parser.add_argument("-m", "--mapping", help="Mapping file with list of sequences", required = True)

args = parser.parse_args()


def main():
	try: 
		taxa_only()

		
	except IndexError:
		print "*** ERROR USAGE: *** \n Correct Format for Use:\nsearch-fasta.py -i <input.fasta> -o <output.fasta> -m <mappingfile>\n\n"
		print "desired_ids.txt should be in the format and must match the exact fasta header in the input file\n"
		print "desired_id1\ndesired_id2\ndesired_idN\n"

	except KeyError:
		print "Fasta header does not match the text file IDs. Check to make sure the fasta header is equal to desired IDs in text file."

        except IOError:
                print "*** ERROR - INCORRECT USAGE ***\n"
                print "Script Usage: \n search-fasta.py -i <input.fasta> -o <output.fasta> -m <mappingfile> \n"

def taxa_only():
	fasta_file = args.input
	wanted_file = args.mapping
	out_file = args.output

	wanted = set()

	with open(wanted_file) as f:
		for line in f:
			line = line.strip('\n')
		#lines = lines.rstrip('\n')
			if line != "":
				wanted.add(line)

	fasta_sequences = SeqIO.parse(open(fasta_file), "fasta")

	count = 0

	with open(out_file, "w") as h:
		for seq in fasta_sequences:

			#if seq.id.split("_1_")[0] in wanted:
			if seq.id in wanted:

				count = count + 1
				SeqIO.write([seq], h, "fasta")



if __name__ == '__main__':
	main()
