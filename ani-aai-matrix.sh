#!/bin/bash

# This script generates an ANI/AAI matrix with clustered PDF output of a set of genomes
# Mostly a wrapper from the Enveomics tools
# Author: Alex Chase
# Date: 2017-07-23
# Usage: ani-aai-matrix.sh -i .fasta -m ani -t 8

####################	DISCLAIMER!!!	####################
# must have:
# NCBI BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and 
# Enveomics (https://github.com/lmrodriguezr/enveomics) 
# installed and in your path for this to work!!!
####################	DISCLAIMER!!!	####################



usage(){
	echo "$(basename "$0") [-h] 
This program to calculate pairwise ANI or AAI values for all genomes in current directory and output matrix
All files in current directory with the same file extension MUST be same format
(i.e. either all nucleotide or all amino acid)
In other words, you cannot have .fa files that represent both nucleotide genomes and translated amino acid genomes
correct use:
	ani-aai-matrix.sh -i .fasta [.fna|.fa] -m [aai|ani] -t [1-8]
where:
	--help or -h  			show this help text
	--input or -i 			input files with extension for genomes (either .fna|.faa|.fasta|.fa )
	--method or -m 			degignate whether you wish to run AAI or ANI calculations
	--thread or -t 			number of threads to use 
					(i.e. how many pairwise comparisons can be run at once) (default=8)"
}

# get user defined options
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit 1
fi

# default thread option
THREAD=8

# get user input for other parameters
while test $# -gt 0; do
		case "$1" in
				-i)
						shift
						if test $# -gt 0; then
								export files=$1
						else
								echo "no input file specified"
								exit 1
						fi
						shift
						;;
				--input*)
						export files=`echo $1 | sed -e 's/^[^=]*=//g'`
						shift
						;;
				-m)
						shift
						if test $# -gt 0; then
								export METHOD=$1
						else
								echo "no method specified"
								exit 1
						fi
						shift
						;;
				--method*)
						export METHOD=`echo $1 | sed -e 's/^[^=]*=//g'`
						shift
						;;
				-t)
						shift
						if test $# -gt 0; then
								export THREAD=$1
						else
								echo "no threads specified"
								exit 1
						fi
						shift
						;;
				--thread*)
						export THREAD=`echo $1 | sed -e 's/^[^=]*=//g'`
						shift
						;;
				*)
						break
						;;
		esac
done

FILES='*'$files

# program the ANI or AAI process!

# get current working directory and make temporary working file, remove older files
REFDIR=$PWD
TEMP=$REFDIR/tempfiles
rm -rf $TEMP
rm -f $REFDIR/total.txt; touch $REFDIR/total.txt
mkdir -p $TEMP
# move everything to new temporary directory
cp $REFDIR/$FILES $TEMP
cd $TEMP

# get an estimate how long this will take...
COUNT=$(ls $TEMP/$FILES | wc -l)

echo -e "\nRunning $METHOD on $COUNT genomes"

# get number of pairwise comparisons needed for matrix
sum=0
for (( i=0;i<=COUNT;i++ ))
do
	(( sum = sum + i ))
done

echo -e "There are $sum pairwise comparisons to make. Sit back and relax!\n\n"

# set thread limit for parallel processes
# set the max background jobs to be run at a given time to not overkill processor
# you can set it to more than 

function max_bg_procs {
	if [[ $# -eq 0 ]] ; then
			echo "Usage: max_bg_procs NUM_PROCS.  Will wait until the number of background (&)"
			echo "  bash processes (as determined by 'jobs -pr') falls below NUM_PROCS"
			return
	fi
	local max_number=$((0 + ${1:-0}))
	while true; do
			local current_number=$(jobs -pr | wc -l)
			if [[ $current_number -lt $max_number ]]; then
					break
			fi
			sleep 1
	done
}


# define functions for both AAI and ANI processes

function aai {
	for f in $FILES
	do
		# to get names of genomes without file extension
		gen1=$( basename $f $files )
		echo "Starting genome ${gen1} ..."

		for h in $FILES
		do
			max_bg_procs $THREAD

			gen2=$( basename $h $files )
			# no need to do pairwise if reciprocal ANI values were aleady calculated
			if grep -Fq "${gen2};${gen1}" ${REFDIR}/total.txt; then
				continue

			else
				aai.rb -1 $f -2 $h -T ${gen1}'@'${gen2}.tab -q &>/dev/null; fi & done 

		wait
		# wait for parellization to finish for genome1 or limit from user defined threads

		# combine recently processed .tab files to master - limit number of files created overall...
		# this will combine all .tab files for pairwise comparisons for one genome at a time (for the loop)
		for i in *.tab
		do
			filename=$( basename $i .tab )
			gen1=$(echo $filename | cut -f1 -d'@')
			gen2=$(echo $filename | cut -f2 -d'@')
			aai=$(cat $i | cut -f1)
			echo ${gen1}';'${gen2}';'${aai} > ${gen1}.${gen2}.final.txt & done
		wait

		cat *.final.txt >> $REFDIR/total.txt
		

		# remove old files
		rm *.tab
		rm *.final.txt
		
	done
}

function ani {
	for f in $FILES
	do
		# to get names of genomes without file extension
		gen1=$( basename $f $files )
		echo "Starting genome ${gen1} ..."

		for h in $FILES
		do
			max_bg_procs $THREAD

			gen2=$( basename $h $files )
			# no need to do pairwise if reciprocal ANI values were aleady calculated
			if grep -Fq "${gen2};${gen1}" ${REFDIR}/total.txt; then
				continue

			else
				ani.rb -1 $f -2 $h -T ${gen1}'@'${gen2}.tab -q &>/dev/null; fi & done 

		wait
		# wait for parellization to finish for genome1 or limit from user defined threads

		# combine recently processed .tab files to master - limit number of files created overall...
		# this will combine all .tab files for pairwise comparisons for one genome at a time (for the loop)
		for i in *.tab
		do
			filename=$( basename $i .tab )
			gen1=$(echo $filename | cut -f1 -d'@')
			gen2=$(echo $filename | cut -f2 -d'@')
			ani=$(cat $i | cut -f1)
			echo ${gen1}';'${gen2}';'${ani} > ${gen1}.${gen2}.final.txt & done
		wait

		cat *.final.txt >> $REFDIR/total.txt
		

		# remove old files
		rm *.tab
		rm *.final.txt
		
	done
}

# make sure user inputs are correct and then initialize functions

# if file extension does not exist, exit
file1=$(echo $FILES | cut -f1 -d' ')

if [ -f $file1 ]; then	
	:
else
	printf "ERROR!\n
	There are no files in current directory with that extension\n"
	exit 1
fi

# if method is not valid, exit
if [[ $METHOD == "ani" ]]
	then
		ani
		NUM=95 
elif [[ $METHOD == "aai" ]]
	then
		aai
		NUM=90
else
	printf "ERROR!\n
	Not a valid selection for method for -m \n"
	exit 1
fi

# replace ';' with tabs to make final output file for matrix construction
tr ';' '\t' < $REFDIR/total.txt > $REFDIR/total2.txt
mv $REFDIR/total2.txt $REFDIR/total.txt


# # combine ANI or AAI results by pairwise interaction and combine results
# for i in *.tab
# do
# 	gen1=$(echo $i | cut -f1 -d'@')
# 	gen2=$(echo $i | cut -f2 -d'@')
# 	ani=$(cat $i | cut -f1); echo -e ${gen1}'\t'${gen2}'\t'${ani} > ${gen1}'@'${gen2}.final.txt & done
# wait

# cat *.final.txt > $REFDIR/total.txt

cd $REFDIR

# generate the matrix

echo "
library(reshape2);
library(gplots);
library(igraph);
dat <- read.table('total.txt', sep = '\\t', header = FALSE, as.is = TRUE);
g <- graph.data.frame(dat, directed = FALSE);
mat <- as.matrix(get.adjacency(g, attr = 'V3', sparse = FALSE));
write.table(
	mat,
	file = '${METHOD}-matrix.txt', quote = F, sep = '\\t',
	col.names = NA, row.names = T);
attach(dat)
breaks = seq(min(V3), max(V3), length.out = 100);
gradient1 = colorpanel( sum( breaks[-1] < ${NUM} ), 'firebrick', 'indianred', 'lightpink');
gradient2 = colorpanel( sum( breaks[-1] >= ${NUM} ), 'palegreen', 'green', 'forestgreen' );
hm.colors = c(gradient1,gradient2);
pdf('${METHOD}.heatmap.pdf', height = 10, width = 15);
heatmap.2(mat,scale='none',breaks=breaks,col=hm.colors,
          trace = 'none', dendrogram = 'row', 
          margin=c(5,10), na.rm = T, na.color = 'white');
dev.off();
" | R --vanilla > summary.log 2>&1

rm -rf $TEMP


