#!/usr/bin/python

# This script renames a fasta file with a given reference file (tab-delimited)
# Author: Alex Chase
# Date: 2014-11-23
# Usage: python fasta-rename.py input_fasta_file.fasta text_file_with_ids.txt output_fasta_file_name.fasta

from Bio import SeqIO
import sys
import os

"""
Create a new fasta file with new IDs.
IDs should be referenced in a txt file in the following format:

old name on fasta header    new name for fasta header
11189                       11189_Microbacteriaceae_Curtobacterium

Call script:
fasta.rename.py input_fasta_file.fasta text_file_with_ids.txt output_fasta_file_name.fasta
"""

def main():
	try: 
		rename()
		
		
	except IndexError:
		print("***ERROR USAGE:***\n\nCorrect Format for Use:\nfasta.rename.py input_fasta_file.fasta text_file_with_ids.txt output_fasta_file_name.fasta\n")
		print "Text file should be in the following format (tab delimited):\n"
		print "Old Name\tNew Name for Fasta Header"
		print "11189\t\t11189_curtobacterium\n"	

        except IOError:
                print "***ERROR - INCORRECT USAGE***\n"
                print "Do not put '-i' or any other parameter marker.\n"
                print "Script Usage:\nfasta-rename.py input_fasta_file.fasta text_file_with_ids.txt output_fasta_file_name.fasta\n"

def rename():

	#input fasta file
	in_file = sys.argv[1]

	# text file with desired names in this format:
	# create a dict to hold our old and new names mapping
		# old name on fasta header    new name for fasta header
		# 11189                       11189_Microbacteriaceae_Curtobacterium
	wanted_file = sys.argv[2]

	#name of fasta file for output
	output_file = sys.argv[3]

	
	fasta_file_read = open(in_file, "r")
	fasta_wanted_file = open(wanted_file, "rU")
	fasta_file_new = open(output_file, "w")
 



	data={}
 
	# read in this file, split it into a dict
	for line in fasta_wanted_file.readlines():
		line_split=line.strip('\n').split('\t')
		data[line_split[0]]=line_split[1]

 
	# create an iterable to hold the new data
	new_seq=[]
 
	# iterate over seq, updating the name.  This is going to give us something like:
	#
	# >11189_Microbacteriaceae_Curtobacterium

	count=0

	for record in SeqIO.parse(fasta_file_read,'fasta'):
		new_record_name = data[record.id]
		record.id=new_record_name
		record.name=''
		record.description=''
		new_seq.append(record)
		count=count+1
 
	# write the whole thing out
	SeqIO.write(new_seq, fasta_file_new, 'fasta')

	
	fasta_wanted_file.close()
	fasta_file_new.close()
	fasta_file_read.close()

if __name__ == '__main__':
	main()
