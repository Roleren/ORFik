#!/bin/bash

#HT 12/02/19
# script to process and align genomic fastq datasets
#1 STAR INDICES (If not existing, must be done in seperate script for now)
#2 make directories
#3 Trim adaptor
#4 (alternative:) Remove merged contaminants
#4 Remove PhiX
#5 Remove rRNA (from silva)
#6 Remove organism specific ncRNA (from ensembl, zebrafish does not contain tRNA)
#7 Remove organism specific tRNA (from tRNAscan-SE)
#8 Map to organism specific reference

usage(){
cat << EOF
usage: $0 options

script to process and align genomic single end or paired end fastq datasets

OPTIONS:
	Important options:
	-f	path to input folder (only fasta/q files in folder allowed for now)
		Must be file types of: fasta, fa, fastq, fq or gz (single or paired end reads)
	-o	path to output dir
	-p	paired end? (yes, defualt: no)
	-l	minimum length of reads (default: 20)
	-T	max mismatches of reads (default: 3)
	-g	genome dir for all STAR indices
	-s	steps of depletion and alignment wanted:
		(a string: which steps to do? (default: "tr-ge", write "all" to get all: "tr-ph-rR-nc-tR-ge",
		   or tr-co-ge, depending on if you merged contaminants or not.)
			 tr: trim, co: contaminants, ph: phix, rR: rrna, nc: ncrna, tR: trna, ge: genome)
		Write your wanted steps, seperated by "-". Order does not matter.
		To just do trim and alignment to genome write -s "tr-ge"
	-a	adapter sequence for trim (found automaticly if not given), also you can write -a "disable",
		to disable it
	-t	trim front (default 3) How many bases to pre trim reads 5' end,
	        as it frequently represents an untemplated addition during reverse transcription.
	-z	trim tail (default 0) How many bases to pre trim reads on 3' end.
	-A	Alignment type: (default Local, EndToEnd (Local is Local, EndToEnd is force Global))
  -B Allow introns (default yes (1), else no (0))
  -M  Max multimapping (default 10) Set to 1 to get only unique reads. Only applies for genome
      step, not the depletion step.
	Path arguments:
	-S   path to STAR (default: ~/bin/STAR-2.7.0c/source/STAR)
	-P   path to fastp (trimmer) (default: ~/bin/fastp)
	-C	 path to cleaning script, internal

	Less important options:
	-r	resume?: a character (defualt n) (n for new start fresh with file f from point s,
			             (if you want a continue from crash specify the step you want to start
				      from, as defined in -s, start on genome, do -r "ge")
	-m	max cpus allowed (defualt 90)
	-i	include subfolders (defualt n, for no), if you want subfolder do y, for yes.
	-q	a character, Do quality filtering:
	    yes: "default" no: "disable". Uses default fastp QF.
	-k	a character, Keep loaded genomes STAR index:
	    yes (y), no: Remove loaded (n), no shared genome (noShared),
	    default (n)
	-K  Keep contaminant aligned files. Default: "no", alternative: "yes".
	-X  Kept contaminant output file type ("bam", "fastq")
	-u  Keep unaligned reads from genome alignment step.
	      Default "None" (no), alternative: "Fastx""
	-h	this help message

fastp location must be: ~/bin/fastp
STAR location must be: ~/bin/STAR-2.7.0c/source/STAR

NOTE: if STAR is stuck on load, run this line:
STAR --genomeDir /export/valenfs/data/references/Zv10_zebrafish_allsteps --genomeLoad Remove

example usage: RNA_Align_pipeline_folder.sh -f <in.fastq.gz> -o <out_dir>

EOF
}

# Default arguments:
echo "##############################################"
echo $'\nArguments for folder run are the following:'
min_length=20
mismatches=3
gen_dir=""
allSteps="tr-ge"
steps=$allSteps
resume="n"
alignment="Local"
allow_introns=0
adapter="auto"
quality_filtering="disable"
maxCPU=90
multimap=10
subfolders="n"
trim_front=3
trim_tail=0
paired="no"
align_single=""
cleaning=""
keepContam="no"
keepContamType="bam"
keep_unmapped_genome="None"
STAR="~/bin/STAR-2.7.0c/source/STAR"
fastp="~/bin/fastp"
while getopts ":f:o:p:l:T:g:s:a:t:A:B:r:m:k:K:M:S:i:P:I:X:C:q:u:z:h:" opt; do
    case $opt in
    f)
        in_dir=$OPTARG
        echo "-f input folder: $OPTARG"
	      ;;
    o)
        out_dir=$OPTARG
        echo "-o output folder: $OPTARG"
        ;;
    p)
	      paired=$OPTARG
        echo "-p paired end: $OPTARG"
        ;;
    l)
        min_length=$OPTARG
        echo "-l minimum length of reads: $OPTARG"
        ;;
    T)
        mismatches=$OPTARG
        echo "-T max mismatches of reads: $OPTARG"
        ;;
    g)
        gen_dir=$OPTARG
        echo "-g genome dir for all STAR indices: $OPTARG"
        ;;
    s)
        steps=$OPTARG
        echo "-s steps to do: $OPTARG"
        ;;
    a)
        adapter=$OPTARG
        echo "-a adapter sequence: $OPTARG"
        ;;
    q)
	      quality_filtering=$OPTARG
	      echo "-q quality filtering: $quality_filtering"
        ;;
    t)
        trim_front=$OPTARG
        echo "-t trim front number: $OPTARG"
        ;;
    z)
        trim_tail=$OPTARG
        echo "-z trim tail (nt): $OPTARG"
        ;;
    A)
        alignment=$OPTARG
        echo "-A alignment type: $OPTARG"
        ;;
    B)
        allow_introns=$OPTARG
        echo "-B allow_introns: $OPTARG"
        ;;
    r)
      	resume=$OPTARG
      	echo "-r resume (new: n, or step as ge): $OPTARG"
        ;;
    m)
      	maxCPU=$OPTARG
      	echo "-m maxCPU: $OPTARG"
        ;;
    M)
      	multimap=$OPTARG
      	echo "-m max multimap: $OPTARG"
        ;;
    i)
      	subfolders=$OPTARG
      	echo "-i subfolders: $OPTARG"
        ;;
    S)
      	STAR=$OPTARG
      	echo "-S STAR location: $OPTARG"
        ;;
    P)
      	fastp=$OPTARG
      	echo "-P fastp location: $OPTARG"
        ;;
    C)
      	cleaning=$OPTARG
      	echo "-C cleaning location: $OPTARG"
        ;;
    I)
      	align_single=$OPTARG
      	echo "-I align_single location: $OPTARG"
        ;;
    K)
      	keepContam=$OPTARG
      	echo "-K Keep contamination reads: $OPTARG"
        ;;
    k)
      	keepLast=$OPTARG
      	echo "-k Keep Star Index loaded: $OPTARG"
        ;;
    X)
      	keepContamType=$OPTARG
      	echo "-X Contamination reads type: $OPTARG"
        ;;
    u)
      	keep_unmapped_genome=$OPTARG
      	echo "-u Keep unmapped genome reads: $OPTARG"
        ;;
    h)
        usage
        exit
        ;;
    ?)
        echo "Invalid option: -$OPTARG"
        usage
        exit 1
        ;;
    esac
done

echo $'\n'
# steps == "all" is never called from R
if [ "$steps" == "all" ]; then
	steps="tr-ph-rR-nc-tR-ge"
fi
if [ -z "$out_dir" ]; then
	echo "Error, out directory (-o) must be speficied!"
	exit 1
fi

# Normal paired end read run
# 1 check if there exists a grouping, and the are equal number in each grouping
# 2 make a list for each end
# send each end into the paired end function
#TODO: Fix checks, incase files are not in order
function findPairs()
{
	f=($1) && shift
	myArray=($@)

	if [ $(( ${#myArray[@]} % 2 )) == 1 ]; then
		echo "folder must have even number of fasta/q files, for paired end run!"
		exit 1
	fi

  #i=0
	#suffix="_001.fastq.gz"
        # remove suffix to get matched pairs
	#for i in "${!myArray[@]}"; do

		#bn=${myArray[i]} | sed  -e "s/$suffix$//"
		#echo ${myArray[i]} | sed  -e "s/$suffix$//"
        	#echo $(basename ${bn})

    	#done
  echo "#############################################"
  keep="y" # Start with keeping genome loaded
	for ((x=0; x<${#myArray[@]}; x = x + 2));
	do
    echo "Paired end mode for files:"
    echo "Forward: $f/${myArray[x]}"
    echo "Reverse: $f/${myArray[x+1]}"
    echo "Files  $((x + 1)) and $((x + 2)) / $numOfFiles"
		a="$f/${myArray[x]}"
		b="$f/${myArray[x+1]}"
		i=$((x + 2))
		if [[ $i == $numOfFiles ]];then
		  echo "Starting last run:"
			keep=${keepLast} # i.e. last file
		fi

		eval $align_single -o "$out_dir" -f "$a" -F "$b"  -a "$adapter" -q "$quality_filtering" -s "$steps" -r "$current" -l "$min_length" -T $mismatches -g "$gen_dir" -m "$maxCPU" -A "$alignment" -B "$allow_introns" -t "$trim_front" -z "$trim_tail" -k $keep -K $keepContam -u $keep_unmapped_genome -P "$fastp" -S "$STAR"
    echo "-------------------------------------------"
	done
}

# Find pairs with subfolders allowed
function findPairsSub()
{
	myArray=($@)

	if [ $(( ${#myArray[@]} % 2 )) == 1 ]; then
		echo "Folder must have even number of fasta/q files, for paired end run!"
		exit 1
	fi
  echo "#############################################"
	for ((x=0; x<${#myArray[@]}; x = x + 2));
	do
    echo "Paired end subfolder mode for files:"
    echo "Forward: ${myArray[x]}"
    echo "Reverse: ${myArray[x+1]}"
    echo "Files  $((x + 1)) and $((x + 2)) / $numOfFiles"
		a="${myArray[x]}"
		b="${myArray[x+1]}"
		i=$((x + 2))
		if [[ $i == $numOfFiles ]];then
		  echo "Starting last run:"
			keep=${keepLast} # i.e. last file
		fi

		eval $align_single -o "$out_dir" -f "$a" -F "$b"  -a "$adapter" -q "$quality_filtering" -s "$steps" -r "$current" -l "$min_length" -T $mismatches -g "$gen_dir" -m "$maxCPU" -M "$multimap" -A "$alignment" -B "$allow_introns" -t "$trim_front" -z "$trim_tail" -k $keep -K $keepContam -X $keepContamType -u $keep_unmapped_genome -P "$fastp" -S "$STAR"
		echo "-------------------------------------------"

	done
}
# Run per file / pair
# Relative path for non subfolder, full for subfolder
formats='\.\(fasta\|fa\|fastq\|fq\)\(\.gz\)\?$'
if [ $subfolders == "n" ]; then
  	 listOfFiles=$(ls "${in_dir}" | grep -E "${formats}")
  	else
  	 listOfFiles=$(find "${in_dir}" | grep -E "${formats}" | sort)
fi
numOfFiles=$(echo "$listOfFiles" | wc -l)

echo "Total number of files are:"
echo $numOfFiles

# Check if resume, if true, jump to given step
declare -i X
X=0
if [ "$resume" != "n" ]; then
   echo "Resume mode"
   X=$(echo "$steps" | grep -b -o $resume | cut -d: -f1)
fi
length=${#steps}
# For each type in tr-co-ge (do one step at a time)
# This is to keep genome loaded through all samples
# Also easier to continue on crash if done this way
while [ $X -lt $length ]
do
  current=${steps:$X:2}
  echo "Current step:"
  echo $current
  X=$((X + 3))
  i=0
  keep="y"
  if [ $paired == "yes" ]; then
  	if [ $subfolders == "n" ]; then
  		findPairs $in_dir $listOfFiles
  	else
  		findPairsSub $listOfFiles
  	fi
  else
  	for x in $listOfFiles
  	do
  		echo "Single end mode for file: $x"
  		i=$((i + 1))
  		echo "File  $i / $numOfFiles"

  		if [[ $i == $numOfFiles ]];then
  		  echo "Starting last run:"
  			keep=${keepLast}
  		fi
  		if [ $subfolders == "n" ]; then
        x=$in_dir/$x
      fi

  		eval $align_single -o "$out_dir" -f "$x"  -a "$adapter" -q "$quality_filtering" -s "$steps" -r "$current" -l "$min_length" -T $mismatches -g "$gen_dir" -m "$maxCPU" -M "$multimap" -A "$alignment" -B "$allow_introns" -t "$trim_front" -z "$trim_tail" -k $keep -K $keepContam -X $keepContamType -u $keep_unmapped_genome -P "$fastp" -S "$STAR"
  		echo "----------------------------------------------"
  	done
  fi
done
echo "done"
### Cleanup
# Folder cleanup
eval $cleaning $out_dir
# Log command run
echo "./RNA_Align_pipeline_folder.sh $@" > $out_dir/runCommand.log
