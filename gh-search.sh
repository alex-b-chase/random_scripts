#!/bin/bash

####################	DISCLAIMER!!!	####################
# must have:
# HMMer installed and 
# Pfam reference database
# installed and in your path for this to work!!!
####################	DISCLAIMER!!!	####################



usage(){
	echo "$(basename "$0") [-h] 

This program will search annotated genomes against the Pfam reference database and subset known GH/CB genes

All files in current directory with the same file extension MUST be same format

correct use:
	gh-search.sh -i .fasta [.fna|.fa] -t [1-8] -c [0|1]

where:
	--help or -h  			show this help text
	--input or -i 			input files with extension for genomes (either .fna|.faa|.fasta|.fa )
	--conservative or -c 	[0 or 1] enable filtering of predicted Pfam proteins to reduce false positives [0 = no; 1 = yes] (default=1)
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
FILTER=1

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
						export files=`echo $1 | sed -e 's/[=]*=//g'`
						shift
						;;
				-c)
						shift
						if test $# -gt 0; then
								export FILTER=$1
						else
								echo "no method specified"
								exit 1
						fi
						shift
						;;
				--conservative*)
						export FILTER=`echo $1 | sed -e 's/[=]*=//g'`
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
						export THREAD=`echo $1 | sed -e 's/[=]*=//g'`
						shift
						;;
				*)
						break
						;;
		esac
done

FILES='*'$files
PFAMBASE=/Users/alexchase/software/PfamScan

# get current working directory and make temporary working file, remove older files
REFDIR=$PWD
TEMP=$REFDIR/tempfiles
PFAMOUT=$REFDIR/pfam
rm -rf $TEMP
mkdir -p $TEMP
mkdir -p $PFAMOUT

# move everything to new temporary directory
cp $REFDIR/$FILES $TEMP
cd $TEMP

cat $PFAMBASE/GH_classification.txt | cut -f2 > $TEMP/gh-reference.txt

# get an estimate how long this will take...
COUNT=$(ls $TEMP/$FILES | wc -l)

echo -e "\nRunning $METHOD on $COUNT genomes"
echo -e "Starting annotation with HMMer. 
This step can take a while. 
A rough estimate for Pfam annotation of a genome ~3 bp is ~10 min per genome.\n"

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

function sortoutput {
	
	if [ "$FILTER" -eq "1" ]; then
		checkout=$(echo "best")
	elif [ "$FILTER" -eq "0" ]; then
		checkout=$(echo "pfam")	
	else
		echo "Not a valid option for -c"
		exit 1
	fi

	for i in *${checkout}.dat
	do
		gen1=$(echo $i | rev | cut -f3- -d'.' | rev)
		cp $i $PFAMOUT
		rm -f $gen1.gh-families.txt; touch $gen1.gh-families.txt

		while read target; 
		do 

			cat $i | grep -w "$target" >> $gen1.gh-families.txt

		done < $TEMP/gh-reference.txt

		cat $gen1.gh-families.txt | awk 'BEGIN{OFS="\t"}{$1 = v} 1' v="$gen1" > $gen1.temp.txt
		cat $gen1.gh-families.txt | awk -v var="$gen1" 'BEGIN{OFS="\t"}{print $0, var}' | awk -v OFS="\t" '$1=$1' > $gen1.temp2.txt
		cp $gen1.temp2.txt $gen1.fix.txt
		cat $gen1.temp.txt | awk -v OFS="\t" '$1=$1' | cut -f1,4,5 | sort | uniq -c | sed -e 's/ *//' | tr ' ' '\t' | grep 'PF' > $REFDIR/$gen1.condensed.txt
	done
}


function gh {
	for f in $FILES
	do
		# to get names of genomes without file extension
		gen1=$(echo $f | rev | cut -f 2- -d '.' | rev)
		echo "Starting genome ${gen1} ..."

		max_bg_procs $THREAD

		hmmsearch --cpu 1 --cut_ga --domtblout $gen1.pfam.dat $PFAMBASE/Pfam-A.hmm > /dev/null $f & done

	wait

	if [ "$FILTER" -eq "1" ]; then

		echo "Filtering HMMer output using bit scores..."
		OUTFILTER=$(echo "best-")

		for h in *.pfam.dat
		do

			gen1=${h%.pfam.dat}

			hmmer2filtered_best.py $h $gen1.pfam-best.dat pfam & done
		wait

		echo "Combining dataset - almost done!"
		sortoutput

	elif [ "$FILTER" -eq "0" ]; then
		OUTFILTER=$(echo "")

		echo "Combining dataset - almost done!"
		sortoutput

	else
		echo "Not a valid option for -c"
		exit 1
	fi

	# remove old files
	rm *.temp*
	rm *.gh-families.txt
}

# if file extension does not exist, exit
file1=$(echo $FILES | cut -f1 -d' ')

if [ -f $file1 ]; then	
	gh
else
	printf "ERROR!\n
	There are no files in current directory with that extension\n"
	exit 1
fi

cd $REFDIR


echo "
library(gplots);
library(reshape2);
library(tidyr);
file_list = list.files(pattern = 'condensed');
dataset <- do.call('rbind', lapply(file_list, FUN=function(files){read.table(files, header = F, sep = '\\t')}));
names(dataset) <- c('count', 'GH_Strain', 'target', 'pfam_id');	
dataset[,4] <- gsub('\\\..*','', dataset[,4])
gh_class <- read.table('${PFAMBASE}/GH_classification.txt', sep = '\\t', header = T);
gh_dataset <- merge(dataset, gh_class, by=c('pfam_id'));
gh_dataset <- aggregate(count ~ GH_Strain + GH, gh_dataset, sum);
gh_all_genomes <- spread(gh_dataset, GH, count);
gh_all_genomes[is.na(gh_all_genomes)] <- 0;
allD <-melt(gh_all_genomes, id = 'GH_Strain');
allD <- acast(allD, GH_Strain ~ variable, value.var = 'value');
breaks = seq(min(allD), max(allD), length.out = 100);
gradient1 = colorpanel( sum( breaks[-1] <= max(allD) / 2 ), 'white', 'indianred', 'black');
gradient2 = colorpanel( sum( breaks[-1] > max(allD) / 2), 'black', 'green', 'forestgreen' );
hm.colors = c(gradient1,gradient2);
pdf('${OUTFILTER}GH.analysis_clustering.pdf', height = 10, width = 15);
heatmap.2(allD, col = hm.colors, scale='none',
          trace = 'none', dendrogram = 'row', 
          margin=c(5,10), na.rm = T, na.color = 'white');
dev.off();
final <- cbind(TOTAL = rowSums(gh_all_genomes[, -1]), gh_all_genomes);
write.table(final, file = '${OUTFILTER}total-gh-data.csv', sep = ',', quote=F, row.names=F);
" | R --vanilla > summary.log 2>&1

rm *.condensed.txt
rm -r $TEMP
