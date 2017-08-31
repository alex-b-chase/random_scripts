#!/usr/bin/python

# This script extracts data from a genbank file and outputs a text file with needed information
# Author: Alex Chase
# Date: 2014-09-23
# Usage: python extract-data-from-genbank.py genbank_file.gb

from Bio import SeqIO
import sys

"""
File prints out the data associated with GenBank formatted files.
Outputs a tab delimited file in the following format:
organism\tname\ttitle_of_paper\tisolation_source\thost\tsequence\tlen_of_seq
"""

def main():
	try:
		extract_data()
		
	except IndexError:
		print "***ERROR - INCORRECT USAGE***\n\nCorrect Usage:\nextract-data-from-genbank.py genbank_file.gb\n\n"
		print "HELP:\nOutputs a tab delimited file in the following format:\n"
		print "organism\tname\ttitle_of_paper\tisolation_source\thost\tsequence\tlen_of_seq"
		print "\nIf you want to create an output file, just add '> output_file.txt' to command line.\n"

        except IOError:
                print "***ERROR - INCORRECT USAGE***\n"
                print "Do not put '-i' or any other parameter marker.\n"
                print "Script Usage:\nextract-data-from-genbank.py genbank_file.gb\n"


def extract_data():
	for r in SeqIO.parse(sys.argv[1], 'genbank'):
	
		title = isolation_source = host = ""
	
		#get organism
		organism = r.annotations['source']
	
		#get title
		if 'references' in r.annotations:
			ref = r.annotations['references'][0]
			title = ref.title
	
		#get source features
		source = r.features[0]
	
		#get isolation_source
		if 'isolation_source' in source.qualifiers:
			isolation_source = source.qualifiers['isolation_source'][0]
		if 'isolation_source' not in source.qualifiers:
			isolation_source = 'no iso source in file'
	
		#get host
		if 'host' in source.qualifiers:
			host = source.qualifiers['host'][0]
		if 'host' not in source.qualifiers:
			host = 'no host in file'
	
		#get sequence and length
		sequence = repr(r.seq)
		length = len(r)
	
		#print output
		sys.stdout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (organism, r.name, title, isolation_source, host, sequence, length))

if __name__ == '__main__':
	main()
