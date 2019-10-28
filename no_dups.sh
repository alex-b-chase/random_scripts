#!/bin/bash

FILES=$1

for f in $FILES

do

    filename="${f%.*}"

    awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' "$f" | awk '!seen[$1]++' | awk -v OFS="\n" '{print $1,$2}' > "$filename-no-dups.fasta"

done
