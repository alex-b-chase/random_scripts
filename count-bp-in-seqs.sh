#!/bin/bash

FILES=$1



for f in $FILES

do
	echo "######################################"
	echo "This will only count the number of base pairs for every sequence of ALL fasta files in current directory..."
	echo

	echo "Counting number of base pairs per sequence in $f..."
	echo  

	filename="${f%.*}"


	awk '/^>/ {if (seqlen){print seqlen};print;seqtotal+=seqlen;seqlen=0;seq+=1;next;}
		{seqlen=seqlen+length($0)}
		END{print seqlen;print seq" sequences, total length " seqtotal+seqlen}' "$f" > "$filename-bp-count.txt"


	echo "DONE with $filename!"
	echo ""
	echo "######################################"

done